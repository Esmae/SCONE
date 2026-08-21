module endfLaw67_class

  use numPrecision
  use endfConstants
  use genericProcedures,        only : fatalError, binarySearch, searchError, interpolate, isSorted
  use RNG_class,                only : RNG
  use aceCard_class,            only : aceCard
  use correlatedLawENDF_inter,  only : correlatedLawENDF
  use law67Pdf_class,           only : law67Pdf
  use muEndfPdfSlot_class,      only : muEndfPdfSlot
  use angleLawENDF_inter,       only : angleLawENDF
  use angleLawENDFfactory_func, only : new_angleLawENDF

  implicit none
  private

  !!
  !! Constructors
  !!
  interface endfLaw67
    module procedure new_endfLaw67_fromACE
  end interface

  type, public, extends(correlatedLawENDF) :: endfLaw67
    private
    real(defReal), dimension(:), allocatable :: eGrid
    type(law67Pdf),dimension(:), allocatable :: pdfs
    class(angleLawENDF), allocatable         :: muLaw
 
  contains
    ! Interface implementation
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Private procedures
    procedure :: init_fromACE

    end type endfLaw67

    contains

    !!
    !! Samples mu and E_out givent incident energy E_in and random nummber generator
    !! 
    !!
    subroutine sample(self, mu, E_out, E_in, rand)
        class(endfLaw67), intent(in) :: self
        real(defReal), intent(out)   :: mu
        real(defReal), intent(out)   :: E_out
        real(defReal), intent(in)    :: E_in
        class(RNG), intent(inout)    :: rand
        integer(shortInt)            :: idx
        real(defReal)                :: E_min_low, E_max_low
        real(defReal)                :: E_min_up, E_max_up
        real(defReal)                :: E_min, E_max
        real(defReal)                :: factor
        real(defReal)                :: r, eps
        character(100),parameter     :: Here='sample (endfLaw67_class.f90)'

        ! Sample mu from the angular distribution first
        mu = self % muLaw % sample(E_in, rand)

        ! Find Interval index
        idx = binarySearch(self % eGrid, E_in)
        call searchError(idx,Here)

        ! Calculate threshold and sample random number
        eps = (E_in - self % eGrid(idx)) / (self % eGrid(idx+1) - self % eGrid(idx))
        r = rand % get()

        ! NOTE: Unlike in the equivalent tabular distribution for mu bounds of the
        !       energy distribution change between bins. We need to interpolate them
        !       as the result.

        ! Obtain bounds of the lower and upper bin, for the given mu
        call self % pdfs(idx)   % bounds(mu, E_min_low, E_max_low)
        call self % pdfs(idx+1) % bounds(mu, E_min_up,  E_max_up )

        ! Calculate interpolated bounds
        E_min = E_min_low * (ONE - eps) + eps * E_min_up
        E_max = E_max_low * (ONE - eps) + eps * E_max_up

        ! Calculate interpolation between bounds of the distribution from which
        ! outgoing energy was sampled
        if(r < eps) then
          call self % pdfs(idx+1) % sample(mu, E_out, rand)
          factor = (E_out- E_min_up)/(E_max_up - E_min_up)
        else
          call self % pdfs(idx) % sample(mu, E_out, rand)
          factor = (E_out- E_min_low)/(E_max_low - E_min_low)
        end if

        ! Interpolate outgoing energy
        E_out = E_min *(ONE - factor) + factor * E_max

    end subroutine sample


    !!
    !! Returns probability that neutron was emitted at mu & E_out given incident energy E_in
    !! Not entirely if the double interpolation is fine (in muLaw and here)
    !!
    function probabilityOf(self, mu, E_out, E_in) result(prob)
        class(endfLaw67), intent(in) :: self
        real(defReal), intent(in)    :: mu
        real(defReal), intent(in)    :: E_out
        real(defReal), intent(in)    :: E_in
        real(defReal)                :: prob
        integer(shortInt)            :: idx
        real(defReal)                :: prob_1, prob_0, E_1, E_0
        character(100),parameter     :: Here='probabilityOf (endfLaw67_class.f90)'

        ! First get probability of mu from the angular distribution
        prob = self % muLaw % probabilityOf(mu, E_in)

        ! Find interval index
        idx = binarySearch(self % eGrid,E_in)
        call searchError(idx,Here)

        ! Obtain probabilities & energies at boundaries of the interval
        prob_0 = prob * self % pdfs(idx)   % probabilityOf(mu, E_out)
        prob_1 = prob * self % pdfs(idx+1) % probabilityOf(mu, E_out)

        E_0 = self % eGrid(idx)
        E_1 = self % eGrid(idx+1)

        ! Interpolate
        prob = interpolate(E_0, E_1, prob_0, prob_1, E_in)

    end function probabilityOf

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
        class(endfLaw67), intent(inout) :: self

        if(allocated(self % eGrid)) deallocate(self % eGrid)

        if(allocated(self % pdfs)) then
            deallocate(self % pdfs)
        end if

        call self % muLaw % kill()

    end subroutine kill

    !!
    !! Initialise endfLaw67 from ACE card
    !! Need root and offset to go back to energy distiribution in ACE card
    !! after reading angular distribution
    !!
    subroutine init_fromACE(self, ACE, MT, root, offset)
      class(endfLaw67), intent(inout)             :: self
      type(aceCard), intent(inout)                :: ACE
      integer(shortInt), intent(in)               :: MT
      integer(shortInt), intent(in)               :: root
      integer(shortInt), intent(in)               :: offset
      integer(shortInt)                           :: NP, NR, numE, i
      integer(shortInt),dimension(:),allocatable  :: L
      character(100), parameter                   :: Here='init_fromACE (endfLaw67_class.f90)'

      ! Read angular data
      call new_angleLawENDF(self % muLaw, ACE, MT)
      ! Go back to the energy distribution
      call ACE % setRelativeTo(root, offset)

      ! Read number of points without moving read head
      NR = ACE % readInt()
      if (NR /= 0) call fatalError(Here, 'Multiple interpolation regions not supported yet')
      
      numE = ACE % readInt()
      self % eGrid = ACE % readRealArray(numE)

      ! Verify energy grid
      if (any(self % eGrid <= ZERO)) call fatalError(Here,'-ve values in energy grid!')
      if (.not.isSorted(self % eGrid)) call fatalError(Here,'energy grid is not sorted-increasing')

      ! Read locators of outgoing mu-E PDFs and allocate space
      L = ACE % readIntArray(numE)
      allocate(self % pdfs(numE))
    
      do i=1, numE
        call ACE % setToEnergyLaw(L(i))
        call self % pdfs(i) % init(ACE)
      end do

    end subroutine init_fromACE

    !!
    !! Build new instance of endfLaw67 from ACE
    !!
    function new_endfLaw67_fromACE(ACE, MT, root, offset) result(new)
      class(aceCard), intent(inout) :: ACE
      integer(shortInt), intent(in) :: MT
      integer(shortInt), intent(in) :: root
      integer(shortInt), intent(in) :: offset
      type(endfLaw67)               :: new

      call new % init_fromACE(ACE, MT, root, offset)

    end function new_endfLaw67_fromACE

end module endfLaw67_class