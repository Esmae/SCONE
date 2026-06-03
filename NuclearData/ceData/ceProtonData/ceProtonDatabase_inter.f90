module ceProtonDatabase_inter

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, numToChar
  use RNG_class,         only : RNG
  use particle_class,    only : particle, P_PROTON, printType
  use charMap_class,     only : charMap
  use intMap_class,      only : intMap

  ! Nuclear Data Handles
  use nuclideHandle_inter,   only : nuclideHandle
  use materialHandle_inter,  only : materialHandle
  use reactionHandle_inter,  only : reactionHandle
  use nuclearDatabase_inter, only : nuclearDatabase

  ! Material Menu
  use materialMenu_mod,      only : mm_nMat => nMat

  ! Cache
  use ceCache_mod,           only : materialCache, trackingCache

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public ceProtonDatabase_CptrCast

  !!
  !! An abstract base class for all nulcear databases that support CE Proton
  !!
  !! Its primary goal is to contain CE Proton caching logic so there is
  !! no need to reproduce it in each database implementation.
  !!
  !! It is also used by material and nuclide handles for CE Proton data to order an
  !! update of XSs on the cache
  !!
  !! Interface:
  !!   nuclearDatabase Interface
  !!   energyBounds       -> return maximum and minimum energy
  !!   updateTotalMatXS   -> update Total Material XS on CE Proton Cache
  !!   updateMacroXSs     -> update Macroscopic XSs for a selected material
  !!   updateTotalXS      -> update Total XS for a selected nuclide
  !!   updateMicroXSs     -> update Microscopic XSs for a selected nuclide
  !!
  type, public, abstract, extends(nuclearDatabase) :: ceProtonDatabase

  contains

    ! nuclearDatabase Interface Implementation
    procedure :: getTrackingXS
    procedure :: getTrackMatXS
    procedure :: getTotalMatXS
    procedure :: getMajorantXS

    ! Overrides procedure
    procedure :: getEnergyLoss
    procedure :: getMoliereVolatility

    ! Procedures implemented by a specific CE Proton Database
    procedure(updateTotalMatXS), deferred      :: updateTrackMatXS
    procedure(updateTotalMatXS), deferred      :: updateTotalMatXS
    procedure(updateMacroXSs), deferred        :: updateMacroXSs
    procedure(updateTotalXS), deferred         :: updateTotalNucXS
    procedure(updateMicroXSs), deferred        :: updateMicroXSs
    procedure(energyBounds), deferred          :: energyBounds
    procedure(getMaterialMTxs), deferred       :: getMaterialMTxs
    procedure(getMaterialEnLoss), deferred     :: getMaterialEnLoss
    procedure(getMaterialVolatility), deferred :: getMaterialVolatility

  end type ceProtonDatabase

  abstract interface
    !!
    !! Return energy bounds for data in the database
    !!
    !! eMin and eMax are minimun and maximumum energy such that data
    !! for ALL nuclides if avalible
    !!
    !! Args:
    !!   eMin [out] -> minimum value of energy [MeV]
    !!   eMax [out] -> maximum value of energy [MeV]
    !!
    !! Errors:
    !!   None
    !!
    subroutine energyBounds(self, eMin, eMax)
      import :: ceProtonDatabase, defReal
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(out)           :: eMin
      real(defReal), intent(out)           :: eMax
    end subroutine energyBounds

    !!
    !! Make sure that trackXS of material with matIdx is at energy E = E_track
    !! in ceProtonChache
    !!
    !! The tracking xs correspons to the material total cross section unless TMS
    !! is used. In that case, this is the material temperature majorant xs.
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateTrackMatXS(self, E, matIdx, rand)
      import :: ceProtonDatabase, defReal, shortInt, RNG
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: matIdx
      class(RNG), optional, intent(inout)  :: rand
    end subroutine updateTrackMatXS

    !!
    !! Make sure that totalXS of material with matIdx is at energy E
    !! in ceProtonChache
    !!
    !! ANY CHANGE in ceProtonChache is POSSIBLE
    !!   E.G. All material XSs may be updated to energy E
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateTotalMatXS(self, E, matIdx, rand)
      import :: ceProtonDatabase, defReal, shortInt, RNG
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: matIdx
      class(RNG), optional, intent(inout)  :: rand
    end subroutine updateTotalMatXS

    !!
    !! Make sure that the macroscopic XSs for the material with matIdx are set
    !! to energy E in ceProtonCache
    !!
    !! ANY CHANGE in ceProtonChache is POSSIBLE
    !!   E.G. Extra materials may be set to energy E as well
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMacroXSs(self, E, matIdx, rand)
      import :: ceProtonDatabase, defReal, shortInt, RNG
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: matIdx
      class(RNG), optional, intent(inout)  :: rand
    end subroutine updateMacroXSs

    !!
    !! Make sure that totalXS of nuclide with nucIdx is at energy E
    !! in ceProtonChache
    !!
    !! ANY CHANGE in ceProtonChache is POSSIBLE
    !!   E.G. All nuclid XSs may be updated to energy E
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   nucIdx [in]  -> nuclide index that needs to be updated
    !!   kT [in]      -> thermal energy of material [MeV]
    !!   rand [inout] -> random number generator
    !!
    subroutine updateTotalXS(self, E, nucIdx, kT, rand)
      import :: ceProtonDatabase, defReal, shortInt, RNG
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: nucIdx
      real(defReal), intent(in)            :: kT
      class(RNG), optional, intent(inout)  :: rand
    end subroutine updateTotalXS

    !!
    !! Make sure that the microscopic XSs for the nuclide with nucIdx are set
    !! to energy E in ceProtonCache
    !!
    !! ANY CHANGE in ceProtonChache is POSSIBLE
    !!   E.G. Extra nuclides may be set to energy E as well
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   nucIdx [in]  -> material index that needs to be updated
    !!   kT [in]      -> thermal energy of material [MeV]
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMicroXSs(self, E, nucIdx, kT, rand)
      import :: ceProtonDatabase, defReal, shortInt, RNG
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: nucIdx
      real(defReal), intent(in)            :: kT
      class(RNG), optional, intent(inout)  :: rand
    end subroutine updateMicroXSs

    !!
    !! Retrieves the material cross sections for an MT number
    !!
    !! That is, it sums up the contributions of the constituent nuclides for
    !! that MT number
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!
    function getMaterialMTxs(self, E, matIdx, MT) result(xs)
      import :: ceProtonDatabase, defReal, shortInt
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: matIdx
      integer(shortInt), intent(in)        :: MT
      real(defReal)                        :: xs
    end function getMaterialMTxs

    !!
    !! Calculates the energy loss for a proton at a given energy in a material
    !!
    !! It uses the Bethe Block formula. It sums up the contributions of the constituent
    !! nuclides for that material number
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index
    !!   rand [inout] -> random number generator
    !!
    function getMaterialEnLoss(self, E, matIdx, dist, rand) result(deltaE)
      import :: ceProtonDatabase, defReal, shortInt, RNG
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(in)           :: E
      integer(shortInt), intent(in)       :: matIdx
      real(defReal), intent(in)           :: dist
      class(RNG), intent(inout)           :: rand
      real(defReal)                       :: deltaE
    end function getMaterialEnLoss

    !!
    !! Calculates the Brownian motion volatiliry for a proton at a given energy in a material
    !!
    !! It uses the Moliere elastic scattering small angle model. It sums up the contributions
    !! of the constituent nuclides for that material number
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index
    !!   rand [inout] -> random number generator
    !!
    function getMaterialVolatility(self, E, matIdx, dist, rand) result(sigma)
      import :: ceProtonDatabase, defReal, shortInt, RNG
      class(ceProtonDatabase), intent(in) :: self
      real(defReal), intent(in)           :: E
      integer(shortInt), intent(in)       :: matIdx
      real(defReal), intent(in)           :: dist
      class(RNG), intent(inout)           :: rand
      real(defReal)                       :: sigma
    end function getMaterialVolatility

  end interface
contains

  !!
  !! Return tracking XS requested
  !!
  !! See nuclearDatabase_inter for details!
  !!
  !! Error:
  !!   fatalError if particle is not CE Proton
  !!
  function getTrackingXS(self, p, matIdx, what) result(xs)
    class(ceProtonDatabase), intent(inout) :: self
    class(particle), intent(in)             :: p
    integer(shortInt), intent(in)           :: matIdx
    integer(shortInt), intent(in)           :: what
    real(defReal)                           :: xs
    character(100),parameter :: Here = 'getTrackingXS (ceProtonDatabase_inter.f90)'

    ! Process request
    select case(what)

      case (MATERIAL_XS)
        if (matIdx == VOID_MAT) then
          xs = self % collisionXS
        else
          xs = max(self % getTrackMatXS(p, matIdx), self % collisionXS)
        end if

      case (TRACKING_XS)

        ! READ ONLY - read from previously updated cache
        if (p % E == trackingCache(1) % E) then
          xs = trackingCache(1) % xs
          return
        else
          call fatalError(Here, 'Tracking cache failed to update during tracking')
        end if

      case default
        call fatalError(Here, 'Neither material xs nor tracking xs was asked')

    end select

    ! Update Cache
    trackingCache(1) % E  = p % E
    trackingCache(1) % xs = xs

  end function getTrackingXS

  !!
  !! Return tracking XS for matIdx
  !!
  !! This is the regular material total cross section.
  !!
  !! See nuclearDatabase_inter for details!
  !!
  !! Error:
  !!   fatalError if particle is not CE Proton
  !!
  function getTrackMatXS(self, p, matIdx) result(xs)
    class(ceProtonDatabase), intent(inout) :: self
    class(particle), intent(in)             :: p
    integer(shortInt), intent(in)           :: matIdx
    real(defReal)                           :: xs
    character(100),parameter :: Here = 'getTrackMatXS (ceProtonDatabase_inter.f90)'

    ! Check dynamic type of the particle
    if (p % isMG .or. p % type /= P_PROTON) then
      call fatalError(Here, 'Dynamic type of the partcle is not CE Proton but:'//p % typeToChar())
    end if

    ! Check that matIdx exists
    if (matIdx == VOID_MAT) then
      xs = ZERO
      return
    elseif (matIdx < 1 .or. matIdx > mm_nMat()) then
      print *,'Particle location: ', p % rGlobal()
      call fatalError(Here, 'Particle is in an undefined material with index: '&
              //numToChar(matIdx))
    end if

    ! Check Cache and update if needed
    if (materialCache(matIdx) % E_track /= p % E) call self % updateTrackMatXS(p % E, matIdx, p % pRNG)

    ! Return Cross-Section
    xs = materialCache(matIdx) % trackXS

  end function getTrackMatXS

  !!
  !! Return Total XS for matIdx
  !!
  !! See nuclearDatabase_inter for details!
  !!
  !! Error:
  !!   fatalError if particle is not CE Proton
  !!
  function getTotalMatXS(self, p, matIdx) result(xs)
    class(ceProtonDatabase), intent(inout) :: self
    class(particle), intent(in)             :: p
    integer(shortInt), intent(in)           :: matIdx
    real(defReal)                           :: xs
    character(100),parameter :: Here = 'getTotalMatXS (ceProtonDatabase_inter.f90)'

    ! Check dynamic type of the particle
    if (p % isMG .or. p % type /= P_PROTON) then
      call fatalError(Here, 'Dynamic type of the partcle is not CE Proton but:'//p % typeToChar())
    end if

    ! Check that matIdx exists
    if (matIdx == VOID_MAT) then
      xs = ZERO
      return
    elseif (matIdx < 1 .or. matIdx > mm_nMat()) then
      print *,'Particle location: ', p % rGlobal()
      call fatalError(Here, 'Particle is in an undefined material with index: '&
              //numToChar(matIdx))
    end if

    ! Check Cache and update if needed
    if (materialCache(matIdx) % E_tot /= p % E) call self % updateTotalMatXS(p % E, matIdx, p % pRNG)

    ! Return Cross-Section
    xs = materialCache(matIdx) % xss % total

  end function getTotalMatXS

  !!
  !! Return Majorant XS
  !!
  !! See nuclearDatabase_inter for details
  !!
  !! Error:
  !!   fatalError: this shouldn't be called
  !!
  function getMajorantXS(self, p) result(xs)
    class(ceProtonDatabase), intent(inout) :: self
    class(particle), intent(in)             :: p
    real(defReal)                           :: xs
    character(100),parameter :: Here = 'getMajorantXS (ceProtonDatabase_inter.f90)'

    call fatalError(Here, 'This function should not be called!')
    xs = huge(ONE)

  end function getMajorantXS

  !!
  !! Return Material Energy Loss
  !!
  !! See nuclearDatabase_inter for details
  !!
  !! Error:
  !!   fatalError if particle is not CE Proton, is MG, and the material is undefined
  !!
  function getEnergyLoss(self, p, matIdx, dist) result(deltaE)
    class(ceProtonDatabase), intent(inout) :: self
    class(particle), intent(in)            :: p
    integer(shortInt), intent(in)          :: matIdx
    real(defReal), intent(in)              :: dist
    real(defReal)                          :: deltaE
    character(100), parameter :: Here = 'getEnergyLoss (ceProtonDatabase_inter.f90)'

    ! Check dynamic type of the particle
    if (p % isMG .or. p % type /= P_PROTON) then
      call fatalError(Here, 'Dynamic type of the partcle is not CE Proton but:'//p % typeToChar())
    end if

    ! Check that matIdx exists
    if (matIdx == VOID_MAT) then
      deltaE = ZERO
      return
    elseif (matIdx < 1 .or. matIdx > mm_nMat()) then
      print *,'Particle location: ', p % rGlobal()
      call fatalError(Here, 'Particle is in an undefined material with index: '&
              //numToChar(matIdx))
    end if

    deltaE = self % getMaterialEnLoss(p % E, matIdx, dist, p % pRNG)

  end function getEnergyLoss

  !!
  !! Return Material Energy Loss
  !!
  !! See nuclearDatabase_inter for details
  !!
  !! Error:
  !!   fatalError if particle is not CE Proton, is MG, and the material is undefined
  !!
  function getMoliereVolatility(self, p, matIdx, dist) result(sigma)
    class(ceProtonDatabase), intent(inout) :: self
    class(particle), intent(in)            :: p
    integer(shortInt), intent(in)          :: matIdx
    real(defReal), intent(in)              :: dist
    real(defReal)                          :: sigma
    character(100),parameter :: Here = 'getMoliereVolatility (ceProtonDatabase_inter.f90)'

    ! Check dynamic type of the particle
    if (p % isMG .or. p % type /= P_PROTON) then
      call fatalError(Here, 'Dynamic type of the partcle is not CE Proton but:'//p % typeToChar())
    end if

    ! Check that matIdx exists
    if (matIdx == VOID_MAT) then
      sigma = ZERO
      return
    elseif (matIdx < 1 .or. matIdx > mm_nMat()) then
      print *,'Particle prePath location: ', p % prePath % r
      call fatalError(Here, 'Particle prePath in an undefined material with index: '&
              //numToChar(matIdx))
    end if

    sigma = self % getMaterialVolatility(p % E, matIdx, dist, p % pRNG)

  end function getMoliereVolatility

  !!
  !! Cast nuclearDatabase pointer to ceProtonDatabase pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of ceProtonDatabase class
  !!   Target points to source if source is ceProtonDatabase class
  !!
  pure function ceProtonDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(ceProtonDatabase), pointer           :: ptr

    select type(source)
      class is(ceProtonDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function ceProtonDatabase_CptrCast


end module ceProtonDatabase_inter
