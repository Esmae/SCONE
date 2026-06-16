module aceProtonDatabase_class

  use numPrecision
  use endfConstants
  use universalVariables
  use errors_mod,         only : fatalError
  use genericProcedures,  only : numToChar, removeDuplicatesSorted, binarySearch, sampleNormal
  use dictionary_class,   only : dictionary
  use RNG_class,          only : RNG
  use charMap_class,      only : charMap
  use intMap_class,       only : intMap

  ! Nuclear Data Interfaces
  use nuclearDatabase_inter,        only : nuclearDatabase
  use materialHandle_inter,         only : materialHandle
  use nuclideHandle_inter,          only : nuclideHandle
  use reactionHandle_inter,         only : reactionHandle
  use ceProtonDatabase_inter,       only : ceProtonDatabase, ceProtonDatabase_CptrCast
  use neutronXSPackages_class,      only : neutronMicroXSs
  use ceProtonMaterial_class,       only : ceProtonMaterial

  ! Material Menu
  use materialMenu_mod,             only : materialItem, nuclideInfo, mm_nMat => nMat, &
                                           mm_getMatPtr => getMatPtr, mm_nameMap => nameMap

  ! ACE CE Nuclear Data Objects
  use aceLibrary_mod,               only : new_protonAce, aceLib_load => load, aceLib_kill => kill
  use aceCard_class,                only : aceCard
  use aceProtonNuclide_class,       only : aceProtonNuclide


  ! CE CACHE
  use ceCache_mod,                  only : cache_nuclideCache => nuclideCache, &
                                           cache_materialCache => materialCache, &
                                           cache_init => init

  ! Scattering procedures
  use scatteringKernels_func,  only : relativeEnergy_constXS, dopplerCorrectionFactor

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: aceProtonDatabase_TptrCast
  public :: aceProtonDatabase_CptrCast

  !!
  !! A CE Neutron Database based on ACE file format
  !!
  !! It's possible to use probability tables in the unresolved resonance range if
  !! ures is included in the input file
  !!
  !! Sample input:
  !!   nuclearData {
  !!   handles {
  !!   ce { type aceProtonDatabase; aceLibrary <nuclear data path>;
  !!        #avgDist 3.141;# }
  !!
  !! Public Members:
  !!   nuclides    -> array of aceProtonNuclides with data
  !!   materials   -> array of ceProtonMaterials with data
  !!   eBounds     -> array with bottom (1) and top (2) energy bound
  !!   eGridUnion  -> unionised energy grid
  !!   activeMat   -> array of materials present in the geometry
  !!   nucToZaid   -> map to link nuclide index to zaid index
  !!   hasUrr      -> ures probability tables flag, it's false by default
  !!   hasDBRC     -> DBRC flag, it's false by default
  !!
  !! Interface:
  !!   nuclearData Interface
  !!   ceProtonDatabase Interface
  !!
  type, public, extends(ceProtonDatabase) :: aceProtonDatabase
    type(aceProtonNuclide),dimension(:),pointer  :: nuclides  => null()
    type(ceProtonMaterial),dimension(:),pointer  :: materials => null()
    real(defReal), dimension(:), allocatable     :: eGridUnion
    real(defReal), dimension(2)                  :: eBounds   = ZERO
    integer(shortInt),dimension(:),allocatable   :: activeMat

  contains

    ! nuclearDatabase Procedures
    procedure :: kill
    procedure :: matNamesMap
    procedure :: getMaterial
    procedure :: getNuclide
    procedure :: getReaction
    procedure :: init
    procedure :: activate
    procedure :: initMajorant


    ! ceProtonDatabase Procedures
    procedure :: energyBounds
    procedure :: updateTrackMatXS
    procedure :: updateTotalMatXS
    procedure :: updateMacroXSs
    procedure :: updateTotalNucXS
    procedure :: updateMicroXSs
    procedure :: getMaterialMTxs
    procedure :: getMaterialEnLoss
    procedure :: getMaterialVolatility

  end type aceProtonDatabase


contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(aceProtonDatabase), intent(inout) :: self

    ! Clean
    if (associated(self % nuclides)) then
      call self % nuclides % kill()
      deallocate(self % nuclides)
    end if

    if (associated(self % materials)) then
      call self % materials % kill()
      deallocate(self % materials)
    end if

    self % eBounds = ZERO

    if(allocated(self % activeMat)) deallocate(self % activeMat)

  end subroutine kill

  !!
  !! Return pointer to material names map
  !!
  !! See nuclearData_inter for  more details
  !!
  function matNamesMap(self) result(map)
    class(aceProtonDatabase), intent(in) :: self
    type(charMap), pointer                :: map

    map => mm_nameMap

  end function matNamesMap

  !!
  !! Return pointer to material in a database
  !!
  !! See nuclearData_inter for  more details
  !!
  function getMaterial(self, matIdx) result(mat)
    class(aceProtonDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: matIdx
    class(materialHandle), pointer        :: mat

    ! Check bounds and return
    if( 1 <= matIdx .and. matIdx <= size(self % materials)) then
      mat => self % materials(matIdx)
    else
      mat => null()
    end if

  end function getMaterial

  !!
  !! Return pointer to nuclide in a database
  !!
  !! See nuclearData_inter for  more details
  !!
  function getNuclide(self, nucIdx) result(nuc)
    class(aceProtonDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: nucIdx
    class(nuclideHandle), pointer         :: nuc

    ! Check bounds and return
    if( 1 <= nucIdx .and. nucIdx <= size(self % nuclides)) then
      nuc => self % nuclides(nucIdx)
    else
      nuc => null()
    end if

  end function getNuclide

  !!
  !! Return a pointer to a reaction
  !!
  !! See nuclearData_inter for  more details
  !!
  function getReaction(self, MT, idx) result(reac)
    class(aceProtonDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: MT
    integer(shortInt), intent(in)         :: idx
    class(reactionHandle), pointer        :: reac
    integer(shortInt)                     :: idxMT

    ! Catch case of invalid reaction
    !   MT < 0 -> material reaction
    !   MT = 0 -> does not exist
    !   MT = 1 -> N_total has no reaction object
    if (MT <= 1) then
      reac => null()
      return
    end if

    ! Detect invalid indices
    if (idx < 1 .or. idx > size(self % nuclides)) then
      reac => null()
      return
    end if

    ! Get nuclide reaction
    if (MT == N_N_elastic) then
      reac => self % nuclides(idx) % elasticScatter

    else if (MT == N_fission) then
      reac => self % nuclides(idx) % fission

    else
      ! Find index of MT reaction
      idxMT = self % nuclides(idx) % idxMT % getOrDefault(MT, 0)
      ! See if the MT is present or not
      if (idxMT == 0) then
        reac => null()
      else
        reac => self % nuclides(idx) % MTdata(idxMT) % kinematics
      end if

    end if

  end function getReaction

  !!
  !! Calculates the energy loss per cm
  !!
  !! See nuclearDatabase for more details
  !!
  function getMaterialEnLoss(self, E, matIdx, dist, rand) result(deltaE)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(in)            :: E
    integer(shortInt), intent(in)        :: matIdx
    real(defReal), intent(in)            :: dist
    class(RNG), intent(inout)            :: rand
    real(defReal)                        :: deltaE, deltaE1, deltaE2
    integer(shortInt)                    :: i, nucIdx
    real(defReal)                        :: dens, sigma1, sigma2, sumSigma1, &
                                            sumSigma2, r1, r2

    sumSigma1 = ZERO
    sumSigma2 = ZERO

    associate (mat => self % materials(matIdx))

      ! Construct macro XS for the MT number requested
      do i = 1, size(mat % nuclides)
        dens   = mat % dens(i)
        nucIdx = mat % nuclides(i)

        call self % nuclides(nucIdx) % betheBloch(sigma1, sigma2, E, mat % excitEn)

        sumSigma1 = sumSigma1 + sigma1 * dens
        sumSigma2 = sumSigma2 + sigma2 * dens

      end do

    end associate

    r1 = rand % get()
    r2 = rand % get()
    deltaE1 = sumSigma1 * dist ! Deterministic betheBloch part of energy loss
    deltaE2 = sampleNormal(r1, r2) * sqrt(sumSigma2 * dist) ! Brownian motion part of energy loss 
    ! Energy loss is always positive (i.e. energy doesn't increase)
    ! Energy loss must be symmetric about the deterministic component
    ! There is a max and min cutoff
    deltaE = MIN(MAX(deltaE1+deltaE2,0.0_defReal),2*deltaE1)

  end function getMaterialEnLoss

  !!
  !! Calculates the energy loss per cm
  !!
  !! See nuclearDatabase for more details
  !!
  function getMaterialVolatility(self, E, matIdx, dist, rand) result(sigma)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(in)            :: E
    integer(shortInt), intent(in)        :: matIdx
    real(defReal), intent(in)            :: dist
    class(RNG), intent(inout)            :: rand
    real(defReal)                        :: sigma
    integer(shortInt)                    :: i, nucIdx
    real(defReal)                        :: dens, chiAnum, chiAden, chiC, sumChiANum, &
                                            sumChiADen, sumChiA, sumChiC, F, omega

    sumChiANum = ZERO
    sumChiADen = ZERO
    sumChiC    = ZERO
    F = 0.98_defReal  ! Hardcoded here, could be an input later on

    associate (mat => self % materials(matIdx))

      ! Construct macro XS for the MT number requested
      do i = 1, size(mat % nuclides)
        dens   = mat % dens(i)
        nucIdx = mat % nuclides(i)

        call self % nuclides(nucIdx) % moliereScattering(chiAnum, chiAden, chiC, E)

        sumChiC    = sumChiC + chiC * dist * dens
        sumChiANum = sumChiANum + chiAnum
        sumChiADen = sumChiADen + chiAden

      end do

    end associate

    sumChiA = exp(sumChiANum/sumChiADen)
    omega   = sumChiC / (TWO * sumChiA * (ONE - F))
    sigma   = sqrt( sumChiC * ((ONE + omega) * log(ONE + omega) / omega - ONE) / (ONE + F**2) )

  end function getMaterialVolatility

  !!
  !! Return energy bounds for data in the database
  !!
  !! See ceProtonDatabase for more details
  !!
  subroutine energyBounds(self, eMin, eMax)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(out)            :: eMin
    real(defReal), intent(out)            :: eMax

    eMin = self % eBounds(1)
    eMax = self % eBounds(2)

  end subroutine energyBounds

  !!
  !! Retrieves the material cross sections for an MT number
  !!
  !! See ceProtonDatabase for more details
  !!
  function getMaterialMTxs(self, E, matIdx, MT) result(xs)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    integer(shortInt), intent(in)         :: MT
    real(defReal)                         :: xs
    integer(shortInt)                     :: i, nucIdx
    real(defReal)                         :: dens

    xs = ZERO

    associate (mat => self % materials(matIdx))

      ! Construct macro XS for the MT number requested
      do i = 1, size(mat % nuclides)
        dens   = mat % dens(i)
        nucIdx = mat % nuclides(i)

        xs = xs + self % nuclides(nucIdx) % xsOf(MT, E) * dens

      end do

    end associate

  end function getMaterialMTxs

  !!
  !! Make sure that trackXS of material with matIdx is at energy E = E_track
  !! in ceProtonChache
  !!
  !! See ceProtonDatabase for more details
  !!
  subroutine updateTrackMatXS(self, E, matIdx, rand)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), optional, intent(inout)   :: rand

    associate (matCache => cache_materialCache(matIdx), &
               mat      => self % materials(matIdx))

      ! Set new energy
      matCache % E_track = E
      call self % updateTotalMatXS(E, matIdx, rand)
      matCache % trackXS = matCache % xss % total

    end associate

  end subroutine updateTrackMatXS

  !!
  !! Make sure that totalXS of material with matIdx is at energy E
  !! in ceCache
  !!
  !! See ceProtonDatabase for more details
  !!
  subroutine updateTotalMatXS(self, E, matIdx, rand)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), optional, intent(inout)   :: rand
    integer(shortInt)                     :: i, nucIdx
    real(defReal)                         :: dens

    associate (matCache => cache_materialCache(matIdx), &
               mat      => self % materials(matIdx))

      ! Set new energy and clean current total XS
      matCache % E_tot = E
      matCache % xss % total = ZERO

      ! Construct total macro XS
      do i = 1, size(mat % nuclides)
        dens   = mat % dens(i)
        nucIdx = mat % nuclides(i)

        ! Update if needed
        if (cache_nuclideCache(nucIdx) % E_tot /= E) then
          call self % updateTotalNucXS(E, nucIdx, mat % kT, rand)
        end if

        ! Add microscopic XSs
        matCache % xss % total = matCache % xss % total + &
                                 dens * cache_nuclideCache(nucIdx) % xss % total
      end do

    end associate

  end subroutine updateTotalMatXS

  !!
  !! Make sure that the macroscopic XSs for the material with matIdx are set
  !! to energy E in ceProtonCache
  !!
  !! See ceProtonDatabase for more details
  !!
  subroutine updateMacroXSs(self, E, matIdx, rand)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), optional, intent(inout)   :: rand
    integer(shortInt)                     :: i, nucIdx
    real(defReal)                         :: dens

    associate(mat      => self % materials(matIdx), &
              matCache => cache_materialCache(matIdx))

      ! Clean current xss
      call matCache % xss % clean()

      ! Set new energy
      matCache % E_tot  = E
      matCache % E_tail = E

      ! Construct microscopic XSs
      do i = 1, size(mat % nuclides)
        dens   = mat % dens(i)
        nucIdx = mat % nuclides(i)

        ! Update if needed
        if (cache_nuclideCache(nucIdx) % E_tail /= E .or. cache_nuclideCache(nucIdx) % E_tot /= E) then
          call self % updateMicroXSs(E, nucIdx, mat % kT, rand)
        end if

        ! Add microscopic XSs
        call matCache % xss % add(cache_nuclideCache(nucIdx) % xss, dens)
      end do

    end associate

  end subroutine updateMacroXSs

  !!
  !! Make sure that totalXS of nuclide with nucIdx is at energy E
  !! in ceProtonCache
  !!
  !! See ceProtonDatabase for more details
  !!
  subroutine updateTotalNucXS(self, E, nucIdx, kT, rand)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: nucIdx
    real(defReal), intent(in)             :: kT
    class(RNG), optional, intent(inout)   :: rand

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      ! Check if the nuclide needs ures probability tables or S(a,b) at this energy
      nucCache % E_tot  = E
      call nuc % search(nucCache % idx, nucCache % f, E)
      nucCache % xss % total = nuc % totalXS(nucCache % idx, nucCache % f)

    end associate

  end subroutine updateTotalNucXS

  !!
  !! Make sure that the microscopic XSs for the nuclide with nucIdx are set
  !! to energy E in ceProtonCache
  !!
  !! See ceProtonDatabase for more details
  !!
  subroutine updateMicroXSs(self, E, nucIdx, kT, rand)
    class(aceProtonDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: nucIdx
    real(defReal), intent(in)             :: kT
    class(RNG), optional, intent(inout)   :: rand

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      nucCache % E_tail = E

      ! In case the total XS hasn't been retrieved before (during tracking)
      if (nucCache % E_tot /= E) then
        nucCache % E_tot  = E
        call nuc % search(nucCache % idx, nucCache % f, E)
      end if

      call nuc % microXSs(nucCache % xss, nucCache % idx, nucCache % f)

    end associate

  end subroutine updateMicroXSs

  !!
  !! Initialise Database from dictionary and pointer to self
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine init(self, dict, ptr, silent)
    class(aceProtonDatabase), target, intent(inout) :: self
    class(dictionary), intent(in)                    :: dict
    class(nuclearDatabase), pointer, intent(in)      :: ptr
    logical(defBool), optional, intent(in)           :: silent
    logical(defBool)                                 :: loud
    type(materialItem), pointer                      :: mat
    class(ceProtonDatabase), pointer                 :: ptr_ceDatabase
    type(charMap)                                    :: nucSet
    type(aceCard)                                    :: ACE
    character(pathLen)                               :: aceLibPath
    character(nameLen)                               :: name
    integer(shortInt)                                :: i, j, envFlag, nucIdx, Z
    integer(shortInt)                                :: maxNuc
    logical(defBool)                                 :: isFissileMat
    integer(shortInt),dimension(:),allocatable       :: nucIdxs
    real(defReal)                                    :: temp
    integer(shortInt), parameter :: IN_SET = 1, NOT_PRESENT = 0
    character(100), parameter :: Here = 'init (aceProtonDatabase_class.f90)'

    ! Set build console output flag
    if(present(silent)) then
      loud = .not.silent
    else
      loud = .true.
    end if

    ! Verify pointer
    if (.not.associated(ptr, self)) then
      call fatalError(Here,"Pointer needs to be associated with the self")
    end if

    ! Cast pointer to ceProtonDatabase
    ptr_ceDatabase => ceProtonDatabase_CptrCast(ptr)
    if(.not.associated(ptr_ceDatabase)) call fatalError(Here,"Should not happen. WTF?!")

    ! Create list of all nuclides. Loop over materials
    ! Find maximum number of nuclides: maxNuc
    do i = 1, mm_nMat()
      mat => mm_getMatPtr(i)
      maxNuc = max(maxNuc, size(mat % nuclides))

      ! Add all nuclides in material to the map
      do j = 1, size(mat % nuclides)
        name = trim(mat % nuclides(j) % toChar())
        call nucSet % add(name, mat % nuclides(j) % Z)
      end do
    end do

    ! Check for a minimum average collision distance
    if (dict % isPresent('avgDist')) then
      call dict % get(temp, 'avgDist')

      if (temp <= ZERO) then
        call fatalError(Here, 'Must have a finite, positive minimum average collision distance')
      end if

      self % collisionXS = ONE / temp

    end if

    ! Get path to ACE library
    call dict % get(aceLibPath,'aceLibrary')

    if(aceLibPath == '$SCONE_ACE') then
      ! Get Path from enviromental variable
      call get_environment_variable("SCONE_ACE", aceLibPath, status = envFlag)

      ! Process potential errors
      if(envFlag == -1) then
        call fatalError(Here,'$SCONE_ACE EnVar must have length smaller then: '//numToChar(pathLen))

      else if(envFlag == 1) then
        call fatalError(Here,"EnVar $SCONE_ACE does not exist! Need to point to ACE Library")

      else if(envFlag == 2) then
        call fatalError(Here,"Compiler does not support EnVariables. &
                              &Replace $SCONE_ACE with path in input file!")
      else if(envFlag /= 0) then
        call fatalError(Here,"Impossible value of envFlag:"//numToChar(envFlag))

      end if
    end if

    ! Load library
    call aceLib_load(aceLibPath)

    ! Build nuclide definitions
    allocate(self % nuclides(nucSet % length()))
    i = nucSet % begin()
    nucIdx = 1
    do while (i /= nucSet % end())

      name = nucSet % atKey(i)
      Z = nucSet % atVal(i)

      if (loud) print '(A)', "Building: "// trim(name)// " with index: " //numToChar(nucIdx)

      call new_protonACE(ACE, name)
      call self % nuclides(nucIdx) % init(ACE, nucIdx, Z, ptr_ceDatabase)

      ! Store nucIdx in the dictionary
      call nucSet % atSet(nucIdx, i)
      nucIdx = nucIdx + 1
      i = nucSet % next(i)
    end do

    ! Calculate energy bounds
    self % eBounds(1) = self % nuclides(1) % eGrid(1)
    j = size(self % nuclides(1) % eGrid)
    self % eBounds(2) = self % nuclides(1) % eGrid(j)

    do i = 2, size(self % nuclides)
      self % eBounds(1) = max(self % eBounds(1), self % nuclides(i) % eGrid(1))
      j = size(self % nuclides(i) % eGrid)
      self % eBounds(2) = min(self % eBounds(2), self % nuclides(i) % eGrid(j))
    end do

    ! Build Material definitions
    allocate(self % materials(mm_nMat()))
    allocate(nucIdxs(maxNuc))
    do i = 1, mm_nMat()
      mat => mm_getMatPtr(i)

      ! Load nuclide indices on storage space
      ! Find if material is fissile and if stochastic
      ! mixing temperature bounds are respected
      isFissileMat = .false.
      ! Loop over nuclides
      do j = 1, size(mat % nuclides)
        name = trim(mat % nuclides(j) % toChar())

        ! Find nuclide definition to see if fissile
        ! Also used for checking stochastic mixing bounds
        nucIdxs(j) = nucSet % get(name)
        isFissileMat = isFissileMat .or. self % nuclides(nucIdxs(j)) % isFissile()

      end do

      ! Load data into material
      call self % materials(i) % set( name     = mat % name,     &
                                      matIdx   = i,              &
                                      database = ptr_ceDatabase, &
                                      temp     = mat % T,        &
                                      fissile  = isFissileMat,   &
                                      excitation = mat % excitEn )
      call self % materials(i) % setComposition( mat % dens, nucIdxs(1:size(mat % nuclides)))

    end do

    !! Clean up
    call aceLib_kill()

  end subroutine init

  !!
  !! Activate this nuclearDatabase
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine activate(self, activeMat, silent)
    class(aceProtonDatabase), intent(inout)     :: self
    integer(shortInt), dimension(:), intent(in) :: activeMat
    logical(defBool), optional, intent(in)      :: silent

    ! Load active materials
    if(allocated(self % activeMat)) deallocate(self % activeMat)
    self % activeMat = activeMat

    ! Configure Cache
    call cache_init(size(self % materials), size(self % nuclides))

  end subroutine activate

  subroutine initMajorant(self, loud, maxTemp, scaleDensity)
    class(aceProtonDatabase), intent(inout) :: self
    logical(defBool), intent(in), optional   :: loud
    real(defReal), intent(in), optional      :: maxTemp
    real(defReal), intent(in), optional      :: scaleDensity

    ! Dummy routine just to satisfy the interface of ceProtonDatabase. Majorant is not used for protons.
  end subroutine initMajorant

  !!
  !! Cast nuclearDatabase pointer to aceProtonDatabase type pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of aceProtonDatabase type
  !!   Target points to source if source is aceProtonDatabase type
  !!
  pure function aceProtonDatabase_TptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    type(aceProtonDatabase), pointer           :: ptr

    select type(source)
      type is(aceProtonDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function aceProtonDatabase_TptrCast

  !!
  !! Cast nuclearDatabase pointer to aceProtonDatabase class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of aceProtonDatabase class
  !!   Target points to source if source is aceProtonDatabase class
  !!
  pure function aceProtonDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(aceProtonDatabase), pointer          :: ptr

    select type(source)
      class is(aceProtonDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function aceProtonDatabase_CptrCast


end module aceProtonDatabase_class
