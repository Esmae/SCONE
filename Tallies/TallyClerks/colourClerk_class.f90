module colourClerk_class

  use numPrecision
  use tallyCodes
  use universalVariables
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

  ! Basic tally modules
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  implicit none
  private

  !!
  !! Colour Clerk
  !! Estimate the rate at which particles change colour
  !! Particles change colour when they enter certain specified regions of space in 1D
  !! Scores (Number of switches from colour 1 -> colour 2)/Number of particles of colour 1
  !! (Is a cumulative score - this is currently turned off)
  !! Scores for a user-specified number of cycles 
  !!
  !! Sample dictionary input:
  !!
  !!  clerkName {
  !!      type colourClerk;
  !!      cycles 500;
  !!      axis x;
  !!      bins (bin1 bin2);
  !!      bin1 {min -250; max -50}
  !!      bin2 {min 50; max 250}
  !!  }
  !!
  type, public, extends(tallyClerk) :: colourClerk
    private
    real(defReal),dimension(:,:,:),allocatable, public :: value            !! cycle-wise value of tally
    !! Value(i,j,k)=Number of colour switches i->j during cycle k
    !! Value(i,i,k)=Number of particles of colour i during cycle k
    integer(shortInt)                              :: numBins = 0            !! Number of bins
    integer(shortInt)                              :: maxCycles = 0    !! Number of tally cycles
    integer(shortInt)                              :: currentCycle = 0 !! track current cycle
    real(defReal), dimension(:,:), allocatable     :: colourBins        !! bin boundaries
    integer(shortInt)                              :: axis             !! Axis of bin 


  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportCycleEnd
    procedure  :: getBinIdx
    procedure  :: reportDungeonUpdate

    ! Output procedures
    procedure  :: display
    procedure  :: print

    ! Deconstructor
    procedure  :: kill
  end type colourClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(colourClerk), intent(inout) :: self
    class(dictionary), intent(in)             :: dict
    character(nameLen), intent(in)            :: name
    character(nameLen),dimension(:),allocatable :: binNames
    integer(shortInt)                 :: i
    character(nameLen)                 :: str
    character(100), parameter     :: Here = 'init (colourClerk_class.f90)'

    if(.not.dict % isPresent('axis')) call fatalError(Here,"Keyword 'axis' must be present")

    ! Find axis of tally
    call dict % get(str,'axis')
    select case(str)
      case('x')
        self % axis = X_axis

      case('y')
        self % axis = Y_axis

      case('z')
        self % axis = Z_axis

      case default
        call fatalError(Here,'Unrecognised axis: '//trim(str)//' must be x, y or z')
    end select


    ! Assign name
    call self % setName(name)

    ! Read bin names
    call dict % get(binNames,'bins')
    self % numBins = size(binNames)

    ! Load bins
    allocate(self % colourBins(self % numBins,2))

    ! Load bin boundaries
    do i=1, self % numBins
       associate( binName => dict % getDictPtr(binNames(i)))
       call binName % get(self % colourBins(i,1), 'min')
       call binName % get(self % colourBins(i,2), 'max')
       end associate
    end do
    
    ! Read number of cycles for which to track colour changes
    call dict % get(self % maxCycles, 'cycles')


    ! Allocate space for storing colour numbers and colour changes
    allocate(self % value(self % numBins, self % numBins, self % maxCycles))
    self % value = ZERO

  end subroutine init

  !!
  !! Returns array of codes that represent different reports
  !!
  function validReports(self) result(validCodes)
    class(colourClerk),intent(in)      :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [cycleEnd_Code, dungeonUpdate_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(colourClerk), intent(in) :: self
    integer(shortInt)                      :: S

    S = self % maxCycles * self % numBins * self % numBins

  end function getSize




  !!
  !! Returns bin index of particle
  !! BIdx is zero if particle is not in a defined bin
  !!
  subroutine getBinIdx(self, state, bIdx) 
    class(colourClerk), intent(in) :: self
    class(particleState), intent(in)    :: state
    integer(shortInt)              ::  i
    integer(shortInt), intent(inout) :: bIdx
    logical                        :: found
    real(defReal)                  :: pos

    character(100), parameter     :: Here = 'init (colourClerk_class.f90)'
   
    found = .false. 
    ! Loop through bins 
    do i=1, self % numBins
       ! Get position of particle, along axis of tally
       pos = state % r(self % axis)
       if ((self % colourBins(i,1) .lt. pos) .and. (self % colourBins(i,2) .gt. pos)) then
           if (found) then
               call fatalError(Here,'Particle lies in more than one bin')
           else
               found = .true.
               bIdx = i
           end if
       end if
    end do
    if (.NOT. found) then
      bIdx = 0
    end if

    end subroutine getBinIdx

  !! 
  !! Update particle colour flags in Dungeon
  !! Called after reportCycleEnd in eigenphysicsPackage
  !! Also called before cycles start to initialise particle colours
  !!
  subroutine reportDungeonUpdate(self,end)
     class(colourClerk), intent(inout) :: self
     class(particleDungeon), intent(inout) :: end
     integer(shortInt)               :: i, pc, bIdx
    
      ! Loop through particles in Dungeon, and update colour flat
      do i = 1,end % popSize()
        associate( state => end % get(i) )
          ! Each particle state has an associated colour
          ! Default is 0
          ! Find the bin the particle ends the cycle in
          call self % getBinIdx(state, bIdx)
          ! Find current colour of particle
          pc = state % colour
          
         if ((bIdx .ne. 0) .and. (pc .ne. bIdx)) then
              ! Particle has changed colour
              end % prisoners(i) % colour = bIdx
          end if
        end associate
      end do

  end subroutine reportDungeonUpdate


  !!
  !! Process cycle end
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(colourClerk), intent(inout) :: self
    class(particleDungeon), intent(in)        :: end
    type(scoreMemory), intent(inout)          :: mem
    integer(shortInt)                         :: i, j, cc, idx, bIdx, pc
    real(defReal)                             :: totWgt, one_log2, totWgt2
    
    if (self % currentCycle < self % maxCycles) then

      self % currentCycle = self % currentCycle + 1
      cc = self % currentCycle

      ! Copy last cycles data over to this cycle
      ! Because scoring is cumulative
      !if (self % currentCycle .ne. 1) then
      !    self % value(:,:,cc) = self % value(:,:,cc-1)
      !end if

      ! Loop through population scoring:
      ! 1) Particles according to the colour they started the cycle with
      ! 2) Particles that have changed colour.
      ! Particle colours are updated as appropriate
      do i = 1,end % popSize()
        associate( state => end % get(i) )
          ! Each particle state has an associated colour
          ! Default is 0
          ! Find the bin the particle ends the cycle in
          call self % getBinIdx(state, bIdx)
          ! Find current colour of particle
          pc = state % colour
          ! Tally current colour of particle
          if (pc .ne. 0) then
              self % value(pc,pc,cc) = self % value(pc,pc,cc) + state % wgt
          end if
          ! Tally particle switching
          if ((bIdx .ne. 0) .and. (pc .ne. bIdx)) then
              ! Particle has changed colour
              ! Option 1
              ! Particle has changed from white (0) -> new colour
              ! Particle colour change isn't tallied
              ! Particle colour is still updated (changed in updateDungeonReport)
              ! Option 2
              ! Particle has changed from a non-white colour -> new colour
              ! Particle colour change is tallied
              ! Particle colour is updated (changed in updateDungeonReport)
              ! Note: Particles can start 'white", i.e. not in a defined bin
              ! But particles cannot become white. They keep their current colour until they
              ! move into another defined bin
              if (pc .ne. 0) then 
                  ! particle was colour pc, is now colour bIdx
                  self % value(pc, bIdx, cc) = self % value(pc, bIdx, cc) + state % wgt
              end if
          end if
        end associate
      end do
    end if

  end subroutine reportCycleEnd

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(colourClerk), intent(in) :: self
    type(scoreMemory), intent(in)    :: mem

    print *, 'colourClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  subroutine print(self, outFile, mem)
    class(colourClerk), intent(in) :: self
    class(outputFile), intent(inout) :: outFile
    type(scoreMemory), intent(in)    :: mem
    integer(shortInt)                :: i, j, k
    character(nameLen)               :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print name
    name = 'colourTally'


    call outFile % startArray(name, [self % numBins, self % numBins, self % maxCycles])
    do i=1,self % maxCycles
      do j=1, self % numBins
        do k=1, self % numBins
          call outFile % addValue(self % value(j,k,i))
        end do
      end do
    end do
    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  elemental subroutine kill(self)
    class(colourClerk), intent(inout) :: self

    if(allocated(self % value)) deallocate(self % value)
    if(allocated(self % colourBins)) deallocate(self % colourBins)
    self % numBins = 0
    self % currentCycle = 0
    self % maxCycles = 0

  end subroutine kill

end module colourClerk_class
