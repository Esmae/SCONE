module wrightFisher_func
  !! This module contains functions required to simulate spherical Brownian motion
  !! i.e. Brownian motion where the angle of movement undergoes BM on a sphere
  !!
  !! The algorithms can be found in
  !!
  !! "A note on the exact simulation of spherical Brownian motion"
  !!  Mijatović, Aleksandar, Veno Mramor, and Gerónimo Uribe Bravo. 
  !!  Statistics & Probability Letters 165 (2020): 108836.
  !! 
  !! 
  !! "Exact simulation of the Wright–Fisher diffusion"
  !! Jenkins, Paul A., and Dario Spano. (2017): 1478-1509.
  !!

  use numPrecision
  ! use genericProcedures, only : fatalError, sampleNormalBoxMuller
  use RNG_class,         only : RNG
  use genericProcedures, only : fatalError

  implicit none
  ! private

  ! Public Interface
  ! public :: sampleWrightFisher


  !!
  !! Samples wright fisher diffusion
  !!
  !! interface sampleWrightFisher
  !!   module procedure getSampleWrightFisher
  !! end interface sampleWrightFisher

contains
  

  !! 
  !! TODO: Delete after testing
  !!
  !! Copied from genericProcedures
  !!
  !!
  !! Standard normal distribution sampling given two random numbers
  !! It uses the Box-Muller transformation
  !!
  !! Args:
  !!   r1, r2 [in] -> uniform random numbers
  !!
  !! Returns:
  !!   sample -> sample from a standard Gaussian distribution
  !!
  function sampleNormalBoxMuller(r1, r2) result(sample)
    real(defReal), intent(in) :: r1
    real(defReal), intent(in) :: r2
    real(defReal)             :: sample

    sample = sqrt(-TWO * log(r1)) * cos(TWO_PI * r2)

  end function sampleNormalBoxMuller


  !!   
  !! Function required to sample from the ancestral process of Kingmans coalescent with mutation
  !! Called by func2
  !!               
  function func1(nm, numBlocks, thetaSBM) result(s)
      integer(shortInt), intent(in)   :: nm, numBlocks
      real(defReal), intent(in)       :: thetaSBM
      real(defReal)                   :: s
      character(100), parameter :: Here = 'func1 (wrightFisher_func.f90)'
      s = log(thetaSBM+2*nm-1)+log_gamma(thetaSBM+numBlocks+nm-1)
      s = s-log_gamma(thetaSBM+numBlocks)-log_gamma(numBlocks+1.0_defReal)-log_gamma(nm-numBlocks+1.0_defReal)
  end function func1

  !!
  !! Function required to sample from the ancestral process of Kingmans coalescent with mutation
  !! Called by numberOfBlocks and func3
  !! 
  function func2(nm, numBlocks, sigEsq, thetaSBM) result(s)
      integer(shortInt), intent(in)     :: nm, numBlocks
      real(defReal), intent(in)         :: sigEsq, thetaSBM
      real(defReal)                     :: s, tmp
      character(100), parameter :: Here = 'func2 (wrightFisher_func.f90)'

      if (nm>0) then
          tmp = func1(nm, numBlocks, thetaSBM)
          s = exp(tmp-nm*(nm+thetaSBM-1)*sigEsq*0.5)
      else
          s = 1
      end if
  end function func2

  !!
  !! Function required to sample from the ancestral process of Kingmans coalescent with mutation
  !! Called by numberOfBlocks
  !!
  function func3(numBlocks, sigEsq, thetaSBM) result(n)
      integer(shortInt), intent(in)   :: numBlocks
      real(defReal), intent(in)       :: sigEsq, thetaSBM
      integer(shortInt)               :: n
      real(defReal)                   :: b_curr, b_next
      character(100), parameter :: Here = 'func3 (wrightFisher_func.f90)'
      
      n = 0
      b_curr = func2(numBlocks, numBlocks, sigEsq, thetaSBM)
      b_next = func2(numBlocks+1, numBlocks, sigEsq, thetaSBM)

      do while (b_next >= b_curr)
          n = n + 1
          b_curr = b_next
          b_next = func2(n+numBlocks+1, numBlocks, sigEsq, thetaSBM)
      end do
  end function func3


   !!
   !! Simulates from the ancestral process of Kingmans coelescent with mutation
   !! Algorithm 2 from Jenkin's et al
   !!
   function numberOfBlocks(sigEsq, rand) result(numBlocks)
    real(defReal), intent(in)         :: sigEsq ! t
    class(RNG), intent(inout)         :: rand
    integer                           :: numBlocks
    integer                           :: i
    logical                           :: proceed
    real(defReal)                     :: thetaSBM, sigma, r1, r2
    real(defReal)                     :: smin, smax, increment
    integer, allocatable              :: k(:), k_tmp(:)
    character(100), parameter         :: Here = 'numberOfBlocks (wrightFisher_func.f90)'

    numBlocks = 0
    thetaSBM = 1.0

    if (sigEsq < 0.07) then
        sigma = sqrt(2.0 / (3.0 * sigEsq))
        r1 = rand % get()
        r2 = rand % get()
        numBlocks = nint(2/sigEsq + sigma * sampleNormalBoxMuller(r1,r2))  
    else
        ! allocate array
        allocate(k(1))
        k = 0
        proceed = .true.
        r1 = rand % get()
        ! r1 shouldn't be zero
        r1 = max(r1, tiny(r1)) 
        smin = 0.0
        smax = 0.0
        increment = 0

        do while (proceed)
            if (size(k)<(numBlocks+1)) then
                ! Need to increase size of array
                ! Double k in size
                allocate(k_tmp(2*numBlocks))
                k_tmp(1:size(k)) = k
                k_tmp(size(k)+1:) = 0
                call move_alloc(k_tmp, k)
            end if
            k(numBlocks+1) = ceiling(func3(numBlocks,sigEsq,thetaSBM)/2.0_defReal)
            do i = 0, k(numBlocks+1)-1
                increment = func2(numBlocks+2*i,numBlocks,sigEsq,thetaSBM) - func2(numBlocks+2*i+1,numBlocks,sigEsq,thetaSBM)
                smin = smin + increment
                smax = smax + increment
            end do
      
            increment = func2(numBlocks+2*k(numBlocks+1),numBlocks,sigEsq,thetaSBM)
            smin = smin + increment - func2(numBlocks+2*k(numBlocks+1)+1,numBlocks,sigEsq,thetaSBM)
            smax = smax + increment

            do while (smin < r1 .and. r1 < smax)
                do i = 0, numBlocks
                    k(i+1) = k(i+1) + 1
                    increment = func2(i+2*k(i+1),i,sigEsq,thetaSBM)
                    smax = smin + increment
                    smin = smin + increment - func2(i+2*k(i+1) + 1,i,sigEsq,thetaSBM)
                end do
            end do

            if (smin > r1) then
                proceed = .false.
            else
                numBlocks = numBlocks + 1
            end if

        end do
    end if

  end function numberOfBlocks

  !!
  !! Samples from beta distribution beta(a1,a2)
  !! Uses gamma distribution to do so
  !!
  function getBeta(a1, a2, rand) result(b)
      integer(shortInt), intent(in)   :: a1, a2
      class(RNG), intent(inout)       :: rand
      real(defReal)                   :: ga1, ga2, b

      ga1 = getGamma(a1, rand)
      ga2 = getGamma(a2, rand)
      b = ga1/(ga1+ga2)
  end function getBeta

  !!
  !! Only for sampling gamma(a,1)
  !! Samples from gamma distribution
  !!
  function getGamma(a, rand) result(g)
      integer(shortInt), intent(in)   :: a
      class(RNG), intent(inout)       :: rand
      real(defReal)                   :: g, c, d, rn1, r1, r2, b, tmp

      d = a - 1.0/3.0
      c = 1/(3*sqrt(d))
      do 
          r1 = rand % get()
          r2 = rand % get()
          rn1 = sampleNormalBoxMuller(r1,r2)
          b = 1+c*rn1
          if (b .lt. 0.0) cycle
          b = b**3
          r1 = rand % get()
          tmp = 1-0.0331*rn1**4
          if (r1 < tmp) exit
          tmp = 0.5*rn1**2+d*(1-b*log(b))
          if (log(r1)<tmp) exit
      end do
      g = d*tmp
  end function getGamma

  !! 
  !! Takes in sigEsq (std from moliere's theory for small-angle scattering)
  !! and a random number generator
  !!
  function wrightFisher(sigEsq,rand) result(wF)
      real(defReal), intent(in)         :: sigEsq
      class(RNG), intent(inout)         :: rand
      real(defReal)                     :: wF, r1, r2, rn1
      integer(shortInt)                 :: numBlocks
      character(100), parameter         :: Here = 'wrightFisher (wrightFisher_func.f90)'


      if (sigEsq > 1E-9) then
          numBlocks = numberOfBlocks(sigEsq, rand)
          wF = getBeta(1,1+numBlocks, rand) ! Samples from beta distribution
      else ! sigEsq too small for standard algorithm
          wF = sigEsq/2
          r1 = rand % get()
          r2 = rand % get()
          rn1 = sampleNormalBoxMuller(r1,r2)
          wF = abs(rn1*sqrt(sigEsq * wF * (1 - wF)))
      end if

  end function wrightFisher

  !!
  !! Get new direction from current direction
  !! Given radial component of wright fisher diffusion wF
  !! dir: (x,y,z) direction vector, should be normalised
  !! Uses algorithm 1 in Mijatović et al.
  !!
  subroutine getSBMDirection(dirIn, dirOut, wF, rand)
      real(defReal), dimension(3), intent(in)    :: dirIn
      real(defReal), dimension(3), intent(inout) :: dirOut
      real(defReal), intent(in)                  :: wF
      class(RNG), intent(inout)                  :: rand
      real(defReal), dimension(3)                :: arr1, arr2
      real(defReal)                              :: theta, rxy2, denom
      character(100), parameter         :: Here = 'getSBMDirection (wrightFisher_func.f90)'

      ! Sample random angle, for angular component of wright fisher diffusion
      theta = rand % get() * TWO_PI

      ! Compute arr1 (called u in algo 1)
      ! arr1 = ((0,0,1)-dirIn)/|(0,0,1)-dirIn|
      denom = sqrt(dirIn(1)**2+dirIn(2)**2+(dirIn(3)-1)**2)
      if (denom > 1E-10) then
          !! original particle direction isn't too close to z axis
          arr1 = ([0,0,1]-dirIn)/denom
      else
          !! original particle direction is close to z axis. 
          if (dirIn(1) .ne. 0) then
              rxy2 = dirIn(2)/dirIn(1) 
              if (rxy2 < 1E-10) then
                  ! x << y
                  arr1 = [SIGN(1.0_defReal,-1.0*dirIn(1)),0.0_defReal,0.0_defReal]
              else if (rxy2 > 1E+10) then
                  ! y << x
                  arr1 = [0.0_defReal,SIGN(1.0_defReal,-1*dirIn(2)),0.0_defReal]
              else
                  arr1 = [SIGN(sqrt(1/(1+rxy2)), -1*dirIn(1)),SIGN(sqrt(rxy2/(rxy2+1)), -1*dirIn(2)),0.0_defReal]
              end if
          else 
              if (dirIn(2) .ne. 0) then
                  arr1 = [0.0_defReal,SIGN(1.0_defReal,-1*dirIn(2)),0.0_defReal]
              else
                  ! DirIn is (0,0,1)
                  ! Take as the default case the result in the limit 1-z->0, where x=y
                  arr1 = [SIGN(1/sqrt(2.0_defReal), -1*dirIn(1)),SIGN(1/sqrt(2.0_defReal), -1*dirIn(2)),0.0_defReal]
              end if
          end if
      end if
          


      arr2(1) = 2*sqrt(wF*(1-wF))*cos(theta)
      arr2(2) = 2*sqrt(wF*(1-wF))*sin(theta)
      arr2(3) = 1-2*wF

      ! Performing the operation I-2*arr1*arr1^T on arr2
      ! dirOut = arr2 - 2*MATMUL(MATMUL(reshape(arr1,[3,1]),reshape(arr1,[1,3])),reshape(arr2,[3,1]))

      ! Check the above works and then remove TODO
      dirOut(1) = (1-2*arr1(1)**2)*arr2(1)-2*arr1(1)*arr1(2)*arr2(2)-2*arr1(1)*arr1(3)*arr2(3)
      dirOut(2) = (1-2*arr1(2)**2)*arr2(2)-2*arr1(2)*arr1(1)*arr2(1)-2*arr1(2)*arr1(3)*arr2(3)
      dirOut(3) = (1-2*arr1(3)**2)*arr2(3)-2*arr1(3)*arr1(1)*arr2(1)-2*arr1(3)*arr1(2)*arr2(2)

      !dirOut should still be normalised to 1
      if (abs(norm2(dirOut)-1) .gt. 0.0001) then
          call fatalError(here,'New direction vector is not a unit vector')
      end if
      dirOut = dirOut/norm2(dirOut)


  end subroutine getSBMDirection
  



end module wrightFisher_func
