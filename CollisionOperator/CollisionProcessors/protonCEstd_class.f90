module protonCEstd_class

  use numPrecision
  use endfConstants
  use universalVariables,            only : REJECTED
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_PROTON
  use particleDungeon_class,         only : particleDungeon

  ! Abstarct interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data Interfaces
  use nuclearDataReg_mod,            only : ndReg_getprotonCE => getprotonCE
  use nuclearDatabase_inter,         only : nuclearDatabase
  use ceProtonDatabase_inter,       only : ceProtonDatabase
  use ceProtonMaterial_class,       only : ceProtonMaterial, ceProtonMaterial_CptrCast
  use ceProtonNuclide_inter,        only : ceProtonNuclide, ceProtonNuclide_CptrCast

  ! Nuclear reactions
  use reactionHandle_inter,          only : reactionHandle
  use uncorrelatedReactionCE_inter,  only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use neutronScatter_class,          only : neutronScatter, neutronScatter_TptrCast
  use fissionCE_class,               only : fissionCE, fissionCE_TptrCast

  ! Cross-Section Packages
  use neutronXsPackages_class,       only : neutronMicroXSs

  ! Scattering procedures
  use scatteringKernels_func, only : asymptoticScatter, targetVelocity_constXS, &
                                     asymptoticInelasticScatter, targetVelocity_DBRCXS, &
                                     relativeEnergy_constXS

  ! Tally interfaces
  use tallyAdmin_class,       only : tallyAdmin

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for CE neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE    -> minimum energy cut-off [MeV] (default = 0.05)
  !!  maxE    -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!             (default = 200.0)
  !!  threshE -> Energy threshold for explicit treatment of target nuclide movement [-].
  !!             Target movement is sampled if neutron energy E < kT * threshE where
  !!             kT is target material temperature in [MeV]. (default = 400.0)
  !!  threshA -> Mass threshold for explicit treatment of target nuclide movement [Mn].
  !!             Target movement is sampled if target mass A < threshA. (default = 1.0)
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type             protonCEstd;
  !!   #minEnergy       <real>;#
  !!   #maxEnergy       <real>;#
  !!   #energyThreshold <real>;#
  !!   #massThreshold   <real>;#
  !!   }
  !!
  type, public, extends(collisionProcessor) :: protonCEstd
    private
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(ceProtonDatabase), pointer, public :: xsData => null()
    class(ceProtonMaterial), pointer, public :: mat    => null()
    class(ceProtonNuclide),  pointer, public :: nuc    => null()

    !! Settings - private
    real(defReal) :: minE
    real(defReal) :: maxE
    real(defReal) :: threshE
    real(defReal) :: threshA

  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: elastic
    procedure :: inelastic
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs

    ! Local procedures
    procedure,private :: scatterFromFixed
    procedure,private :: scatterFromMoving
    procedure,private :: scatterInLAB

  end type protonCEstd

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(protonCEstd), intent(inout) :: self
    class(dictionary), intent(in)     :: dict
    character(100), parameter :: Here = 'init (protonCEstd_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Read settings for protonCEstd
    ! Maximum and minimum energy
    call dict % getOrDefault(self % minE,'minEnergy',0.05_defReal)
    call dict % getOrDefault(self % maxE,'maxEnergy',200.0_defReal)

    ! Thermal scattering kernel thresholds
    call dict % getOrDefault(self % threshE, 'energyThreshold', 400.0_defReal)
    call dict % getOrDefault(self % threshA, 'massThreshold', 1.0_defReal)

    ! Verify settings
    if (self % minE < ZERO) call fatalError(Here,'-ve minEnergy')
    if (self % maxE < ZERO) call fatalError(Here,'-ve maxEnergy')
    if (self % minE >= self % maxE) call fatalError(Here,'minEnergy >= maxEnergy')
    if (self % threshE < 0) call fatalError(Here,' -ve energyThreshold')
    if (self % threshA < 0) call fatalError(Here,' -ve massThreshold')

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, tally, collDat, thisCycle, nextCycle)
    class(protonCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMicroXSs)                :: microXSs
    real(defReal)                        :: r
    character(100),parameter :: Here = 'sampleCollision (protonCEstd_class.f90)'

    ! Verify that particle is CE neutron
    if (p % isMG .or. p % type /= P_PROTON) then
      call fatalError(Here, 'Supports only CE Neutron. Was given MG '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getProtonCE()
    if (.not.associated(self % xsData)) call fatalError(Here, 'There is no active Neutron CE data!')

    ! Verify and load material pointer
    self % mat => ceProtonMaterial_CptrCast(self % xsData % getMaterial(p % matIdx()))
    if (.not.associated(self % mat)) call fatalError(Here, 'Material is not ceProtonMaterial')

    ! Select collision nuclide
    call self % mat % sampleNuclide(p % E, p % pRNG, collDat % nucIdx, collDat % E)

    ! If nuclide was rejected in TMS loop return to tracking
    if (collDat % nucIdx == REJECTED) then
      collDat % MT = noInteraction
      return
    end if

    self % nuc => ceProtonNuclide_CptrCast(self % xsData % getNuclide(collDat % nucIdx))
    if (.not.associated(self % nuc)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

    ! Select Main reaction channel
    call self % nuc % getMicroXSs(microXss, collDat % E, self % mat % kT, p % pRNG)
    r = p % pRNG % get()
    collDat % MT = microXss % invert(r)


  end subroutine sampleCollision

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, p, tally, collDat, thisCycle, nextCycle)
    class(protonCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    character(100),parameter             :: Here = 'implicit (protonCEstd_class.f90)'

    ! Does nothing

  end subroutine implicit

  !!
  !! Process capture reaction
  !!
  subroutine capture(self, p, tally, collDat, thisCycle, nextCycle)
    class(protonCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead =.true.

  end subroutine capture

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, p, tally, collDat, thisCycle, nextCycle)
    class(protonCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead =.true.

  end subroutine fission

  !!
  !! Process elastic scattering
  !!
  !! All CE elastic scattering happens in the CM frame
  !!
  subroutine elastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(protonCEstd), intent(inout)     :: self
    class(particle), intent(inout)         :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon),intent(inout)   :: thisCycle
    class(particleDungeon),intent(inout)   :: nextCycle
    class(uncorrelatedReactionCE), pointer :: reac
    logical(defBool)                       :: isFixed
    character(100),parameter :: Here = 'elastic (protonCEstd_class.f90)'

    ! Get reaction
    reac => uncorrelatedReactionCE_CptrCast( self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if (.not.associated(reac)) call fatalError(Here,'Failed to get elastic neutron scatter')


    ! Scatter particle
    collDat % A =  self % nuc % getMass()

    ! Retrieve kT from either material or nuclide
    collDat % kT = self % nuc % getkT()

    isFixed = (p % E > collDat % kT * self % threshE) .and. (collDat % A > self % threshA)
    isFixed = isFixed .or. (collDat % kT == ZERO)

    ! Apply criterion for Free-Gas vs Fixed Target scattering
    if (.not. reac % inCMFrame()) then
      call self % scatterInLAB(p, collDat, reac)
    elseif (isFixed) then
      call self % scatterFromFixed(p, collDat, reac)
    else
      call self % scatterFromMoving(p, collDat, reac)
    end if


  end subroutine elastic

  !!
  !! Process inelastic scattering
  !!
  subroutine inelastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(protonCEstd), intent(inout)     :: self
    class(particle), intent(inout)         :: p
    type(tallyAdmin), intent(inout)        :: tally
    type(collisionData), intent(inout)     :: collDat
    class(particleDungeon),intent(inout)   :: thisCycle
    class(particleDungeon),intent(inout)   :: nextCycle
    class(uncorrelatedReactionCE), pointer :: reac
    character(100),parameter  :: Here =' inelastic (protonCEstd_class.f90)'

    ! Invert inelastic scattering and get reaction
    collDat % MT = self % nuc % invertInelastic(collDat % E, p % pRNG)
    reac => uncorrelatedReactionCE_CptrCast(self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if (.not.associated(reac)) call fatalError(Here, "Failed to get scattering reaction")

    ! Scatter particle
    if (reac % inCMFrame()) then
      collDat % A =  self % nuc % getMass()
      call self % scatterFromFixed(p, collDat, reac)
    else
      call self % scatterInLAB(p, collDat, reac)
    end if
    ! Apply weigth change
    if (p % E < self % minE) then
      p % w = p % w
      p % isDead = .true.
    else if (p % E > self % maxE) then
      p % E = self % maxE-10e-3
      p % w = p % w * reac % release(collDat % E)
    else
      p % w = p % w * reac % release(collDat % E)
    end if

    

  end subroutine inelastic

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, tally, collDat, thisCycle, nextCycle)
    class(protonCEstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    if (p % E < self % minE) p % isDead = .true.
    if (p % E > self % maxE) p % E = self % maxE-10e-3

  end subroutine cutoffs

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self, p, collDat, reac)
    class(protonCEstd), intent(inout)         :: self
    class(particle), intent(inout)            :: p
    type(collisionData), intent(inout)        :: collDat
    class(uncorrelatedReactionCE), intent(in) :: reac
    real(defReal)                             :: phi    ! Azimuthal scatter angle
    real(defReal)                             :: E_out, mu

    ! Sample scattering angles and post-collision energy
    call reac % sampleOut(mu, phi, E_out, p % E, p % pRNG)

    ! Update neutron state
    p % E = E_out
    call p % rotate(mu, phi)
    collDat % muL = mu

  end subroutine scatterInLAB

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self, p, collDat, reac)
    class(protonCEstd), intent(inout)         :: self
    class(particle), intent(inout)             :: p
    type(collisionData), intent(inout)         :: collDat
    class(uncorrelatedReactionCE), intent(in)  :: reac
    real(defReal)                              :: phi
    real(defReal)                              :: E_out
    real(defReal)                              :: E_outCM, mu
    integer(shortInt)                          :: MT

    ! Read data
    MT = collDat % MT

    ! Sample mu, phi and outgoing energy
    call reac % sampleOut(mu, phi, E_outCM, p % E, p % pRNG)

    ! Save incident energy
    E_out = p % E

    if (MT == N_N_elastic) then
      call asymptoticScatter(E_out, mu, collDat % A)
    else
      call asymptoticInelasticScatter(E_out, mu, E_outCM, collDat % A)
    end if

    ! Update particle state
    call p % rotate(mu, phi)
    p % E = E_out
    collDat % muL = mu

  end subroutine scatterFromFixed

  !!
  !! Subroutine to perform scattering from moving target
  !! Supports only elastic collisions
  !!
  subroutine scatterFromMoving(self, p, collDat, reac)
    class(protonCEstd), intent(inout)          :: self
    class(particle), intent(inout)             :: p
    type(collisionData),intent(inout)          :: collDat
    class(uncorrelatedReactionCE), intent(in)  :: reac
    integer(shortInt)                          :: nucIdx
    real(defReal)                              :: A, kT, mu
    real(defReal),dimension(3)                 :: V_n           ! Neutron velocity (vector)
    real(defReal)                              :: U_n           ! Neutron speed (scalar)
    real(defReal),dimension(3)                 :: dir_pre       ! Pre-collision direction
    real(defReal),dimension(3)                 :: dir_post      ! Post-collicion direction
    real(defReal),dimension(3)                 :: V_t, V_cm     ! Target and CM velocity
    real(defReal)                              :: phi, dummy
    character(100), parameter :: Here = 'ScatterFromMoving (protonCEstd_class.f90)'

    ! Read collision data
    A      = collDat % A
    kT     = collDat % kT
    nucIdx = collDat % nucIdx

    ! Get neutron direction and velocity
    dir_pre = p % dirGlobal()
    V_n     = dir_pre * sqrt(p % E)

    ! Sample target velocity with constant XS
    V_t = targetVelocity_constXS(p % E, dir_pre, A, kT, p % pRNG)

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t *A)/(A+1)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    U_n = norm2(V_n)
    V_n = V_n / U_n

    ! Sample mu and phi in CM frame
    call reac % sampleOut(mu, phi, dummy, p % E, p % pRNG)

    ! Obtain post collision speed
    V_n = rotateVector(V_n, mu, phi) * U_n

    ! Return to LAB frame
    V_n = V_n + V_cm

    ! Calculate new neutron speed and direction
    U_n = norm2(V_n)
    dir_post = V_n / U_n

    ! Update particle state and calculate mu in LAB frame
    p % E = U_n * U_n
    call p % point(dir_post)
    collDat % muL = dot_product(dir_pre, dir_post)

  end subroutine scatterFromMoving


end module protonCEstd_class
