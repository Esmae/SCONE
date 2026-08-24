module protonEnergyClerk_class

  use numPrecision
  use tallyCodes
  use universalVariables
  use genericProcedures,          only : fatalError
  use display_func,               only : statusMsg
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill

  ! Nuclear Data interface
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Tally Filters
  use tallyFilter_inter,          only : tallyFilter
  use tallyFilterFactory_func,    only : new_tallyFilter

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Responses
  use tallyResponseSlot_class,    only : tallyResponseSlot

  implicit none
  private

  !!
  !! Collision and path estimator of proton energy deposition
  !!
  !! Private Members:
  !!   filter   -> Space to store tally Filter
  !!   map      -> Space to store tally Map
  !!
  !! Interface
  !!   tallyClerk Interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myProtonEnergyClerk {
  !!   type protonEnergyClerk;
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !! }
  !!
  type, public, extends(tallyClerk) :: protonEnergyClerk
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable                  :: filter
    class(tallyMap), allocatable                     :: map
    type(tallyResponseSlot),dimension(:),allocatable :: response


  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportOutColl
    procedure  :: reportPath

    ! Output procedures
    procedure  :: display
    procedure  :: print

  end type protonEnergyClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(protonEnergyClerk), intent(inout)     :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name

    ! Assign name
    call self % setName(name)

    ! Load filetr
    if( dict % isPresent('filter')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Load map
    if( dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if


  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(protonEnergyClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate filter
    if (allocated(self % filter)) then
      deallocate(self % filter)
    end if

    ! Kill and deallocate map
    if (allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if



  end subroutine kill

  !!
  !! Process incoming track length report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportPath(self,p, L, xsData,mem)
    class(protonEnergyClerk), intent(inout)  :: self
    class(particle), intent(in)              :: p
    real(defReal), intent(in)                :: L
    class(nuclearDatabase), intent(inout)    :: xsData
    type(scoreMemory), intent(inout)         :: mem
    type(particleState)                      :: state
    integer(shortInt)                        :: binIdx
    integer(longInt)                         :: addr
    real(defReal)                            :: scoreVal
    character(100), parameter :: Here = 'reportPath (protonEnergyClerk_class.f90)'

    ! Get current particle state
    state = p

    ! Check if within filter
    if (allocated(self % filter)) then
      if (self % filter % isFail(state)) return
    end if

    ! Find bin index
    if (allocated(self % map)) then
      binIdx = self % map % map(state)
    else
      binIdx = 1
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate bin address
    addr = self % getMemAddress() + (binIdx - 1)


    scoreVal = - state % E + p % prePath % E
    call mem % score(scoreVal, addr)



  end subroutine reportPath

  !!
  !! Returns array of codes that represent different reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(protonEnergyClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable    :: validCodes

    validCodes = [outColl_CODE, path_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(protonEnergyClerk), intent(in) :: self
    integer(shortInt)                    :: S

    if(allocated(self % map)) S = self % map % bins(0)

  end function getSize

  !!
  !! Process outgoing collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportOutColl(self, p, MT, muL, xsData, mem)
    class(protonEnergyClerk), intent(inout)  :: self
    class(particle), intent(in)              :: p
    integer(shortInt), intent(in)            :: MT
    real(defReal), intent(in)                :: muL
    class(nuclearDatabase), intent(inout)    :: xsData
    type(scoreMemory), intent(inout)         :: mem
    type(particleState)                      :: state
    integer(shortInt)                        :: binIdx
    integer(longInt)                         :: addr
    real(defReal)                            :: scoreVal
    character(100), parameter :: Here = 'reportOutColl (protonEnergyClerk_class.f90)'


    ! Get current particle state
    state = p

    ! Check if within filter
    if (allocated(self % filter)) then
      if (self % filter % isFail(state)) return
    end if

    ! Find bin index
    if (allocated(self % map)) then
      binIdx = self % map % map(state)
    else
      binIdx = 1
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate bin address
    addr = self % getMemAddress() + (binIdx - 1)

    
    scoreVal = - state % E + p % preCollision % E
    call mem % score(scoreVal, addr)

  end subroutine reportOutColl

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(protonEnergyClerk), intent(in)  :: self
    type(scoreMemory), intent(in)         :: mem

    call statusMsg('protonEnergyClerk does not support display yet')

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(protonEnergyClerk), intent(in)          :: self
    class(outputFile), intent(inout)              :: outFile
    type(scoreMemory), intent(in)                 :: mem
    real(defReal)                                 :: val, std
    integer(shortInt)                             :: i
    integer(shortInt),dimension(:),allocatable    :: resArrayShape
    character(nameLen)                            :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If proton energy clerk has map print map information
    if (allocated(self % map)) then
      call self % map % print(outFile)
    end if

    ! Write results.
    ! Get shape of result array
    resArrayShape = [self % map % binArrayShape()]

    ! Start array
    name ='Res'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i = 1, product(resArrayShape)
      call mem % getResult(val, std, self % getMemAddress() - 1 + i)
      call outFile % addResult(val,std)

    end do

    call outFile % endArray()
    call outFile % endBlock()

  end subroutine print

end module protonEnergyClerk_class
