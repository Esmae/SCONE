module law67Pdf_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError, numToChar, binarySearch, isSorted
  use aceCard_class,      only : aceCard
  use RNG_class,          only : RNG

  use tabularEnergy_class, only : tabularEnergy

  ! Correlated PDFs

  implicit none
  private

  ! Define public interface
  type, public :: law67Pdf
    private
    type(tabularEnergy), dimension(:), allocatable :: ePdfs
    real(defReal), dimension(:), allocatable       :: muGrid
    integer(shortInt)                              :: interpolation

  contains
    ! Public Interface
    generic   :: init => init_fromACE
    procedure :: sample
    procedure :: probabilityOf
    procedure :: bounds
    procedure :: kill
    procedure :: getInterF

    ! Private procedures
    procedure, private :: init_fromACE
  end type law67Pdf

  contains

  !
  ! Helper function to get interpolation factor for a given x and bin
  ! Shamelessly copied from tabularPdf_class.f90
  !
  elemental function getInterF(self, x, bin) result(f)
    class(law67Pdf), intent(in)   :: self
    real(defReal), intent(in)     :: x
    integer(shortInt), intent(in) :: bin
    real(defReal)                 :: f

    f = (x - self % muGrid(bin)) / (self % muGrid(bin + 1) - self % muGrid(bin))

  end function getInterF

  !!
  !! Returns probability that neutron was emmited at E_out knowing mu
  !! Note difference with LAW 61, where mu and E_out are sampled together
  !!
  function probabilityOf(self, mu, E_out) result(prob)
    class(law67Pdf), intent(in) :: self
    real(defReal), intent(in)   :: mu
    real(defReal), intent(in)   :: E_out
    real(defReal)               :: prob
    real(defReal)               :: f
    integer(shortInt)           :: bin

    bin = binarySearch(self % muGrid, mu)

    ! Propability of a given angle at the energy
    f = self % getInterF(mu, bin)
    prob = ( (ONE-f) * self % ePdfs(bin) % probabilityOf(E_out) + &
                          f * self % ePdfs(bin+1) % probabilityOf(E_out))

  end function probabilityOf

  ! 
  ! For LAW 67, bounds return the bounds of the energy distribution for a given mu
  ! Note the difference with LAW 61, where it is independent of mu
  !
  subroutine bounds(self, mu, E_min, E_max)
    class(law67Pdf), intent(in)  :: self
    real(defReal), intent(in)    :: mu
    real(defReal), intent(out)   :: E_min, E_max
    integer(shortInt)            :: bin

    ! Select the bin corresponding to mu
    bin = binarySearch(self % muGrid, mu)

    ! Get bounds of the energy pdf corresponding to the selected bin
    call self % ePdfs(bin) % bounds(E_min, E_max)

  end subroutine bounds

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(law67Pdf), intent(inout) :: self

    if(allocated(self % ePdfs)) then
      call self % ePdfs % kill()
      deallocate(self % ePdfs)
    end if

    if (allocated(self % muGrid)) deallocate(self % muGrid)

  end subroutine kill


  ! 
  ! Sample E_out given mu and random number generator
  ! Note the difference with LAW 61, where mu AND E_out are sampled together
  !
  subroutine sample(self, mu, E_out, rand)
    class(law67Pdf), intent(in) :: self
    real(defReal), intent(in)   :: mu
    real(defReal), intent(out)  :: E_out
    class(RNG), intent(inout)   :: rand
    real(defReal)               :: eps, r
    integer(shortInt)           :: bin

    ! Select the bin corresponding to mu
    bin = binarySearch(self % muGrid, mu)

    ! Sample Energy
    eps = self % getInterF(mu, bin)
    r = rand % get()

    if (r < eps) then
      E_out = self % ePdfs(bin+1) % sample(rand)
    else
      E_out = self % ePdfs(bin) % sample(rand)
    end if

  end subroutine sample

  !!
  !! Initialise law67Pdf from ACE card
  !! Note that for LAW 67, the angular distribution in not read here
  !!
  subroutine init_fromACE(self, ACE)
    class(law67Pdf), intent(inout)                  :: self
    type(aceCard), intent(inout)                    :: ACE
    integer(shortInt)                               :: Nmu, i
    integer(shortInt), dimension(:), allocatable    :: LMU
    character(100), parameter                       :: Here = 'init_fromACE (law67Pdf_class.f90)'

    ! Read interpolation and number of mu bins
    self % interpolation = ACE % readInt()
    if (self % interpolation /= 2) call fatalError(Here, 'Only linear interpolation is supported for LAW 67')
    
    ! Need the mu grid, cause it's not sampled in there
    Nmu = ACE % readInt()
    allocate(self % muGrid(Nmu))
    self % muGrid = ACE % readRealArray(Nmu)

    ! Check bounds and monotonicity of mu grid
    if (any(abs(self % muGrid) > ONE)) call fatalError(Here, 'mu grid must be in [-1,1]')

    ! Check non decreasing mu grid
    if (.not. isSorted(self % muGrid)) call fatalError(Here, 'mu grid must be strictly increasing')

    ! Initialise energy pdfs, one for each mu
    allocate(self % ePdfs(Nmu))
    LMU = ACE % readIntArray(Nmu)

    do i=1, Nmu
      call ACE % setToEnergyLaw(LMU(i))
      call self % ePdfs(i) % init(ACE)
    end do

  end subroutine init_fromACE

end module law67Pdf_class