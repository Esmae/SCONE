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
  !! "Exact simulation of the Wright–Fisher diffusion"
  !! Jenkins, Paul A., and Dario Spano. (2017): 1478-1509.
  !!

  use numPrecision
  use genericProcedures, only : fatalError, sampleNormalBoxMuller
  use RNG_class,         only : RNG

  implicit none
  private

  ! Public Interface
  public :: sampleWrightFisher


  !!
  !! Samples wright fisher diffusion
  !!
  interface sampleWrightFisher
    module procedure getSampleWrightFisher
  end interface sampleWrightFisher

contains


  !!
  !!   
  !!  
  !!               
  !!  
  !! 

  function log_a(k, m, theta) result(x)
      integer(shortInt), intent(in)   :: k, m
      real(defReal), intent(in)       :: theta
      real(defReal)                   :: x
      character(100), parameter :: Here = 'log_a (wrightFisher_func.f90)'

      x = log(theta+2*k-1)+log_gamma(theta+m+k-1)-log_gamma(theta+m)-log_gamma(m+1)-log_gamma(k-m+1)
  end function log_a

  !!
  !!
  !!
  !!
  !!
  !!

  function b(k, m, t, theta) result(x)
      integer(shortInt), intent(in)     :: k, m
      real(defReal), intent(in)         :: t, theta
      real(refReal)                     :: x, tmp
      character(100), parameter :: Here = 'b (wrightFisher_func.f90)'

      if (k>0) then
          tmp = log_a(k, m, theta)
          x = exp(xtmp-k*(k+theta-1)*t*0.5)
      else
          x = 1
  end function b


  function c(m, t, theta) result(n)
      integer(shortInt), intent(in)   :: m
      real(refReal), intent(in)       :: t, theta
      integer(shortInt)               :: n
      real(refReal)                   :: b_curr, b_next
      character(100), parameter :: Here = 'c (wrightFisher_func.f90)'
      
      n = 0
      b_curr = b(m, m, t, theta)
      b_next = b(m+1, m, t, theta)

      do while (b_next >= b_curr)
          n = n + 1
          b_curr = b_next
          b_next = b(i+m+1, m, t, theta)
      end do
  end function c



   !!!!!!!!!!!

   function number_of_blocks(t, rand) result(m)
    real(refReal), intent(in)         :: t
    class(RNG), intent(inout)         :: rand
    integer                           :: i
    logical                           :: proceed
    real(refReal)                     :: theta, mu, sigma, r1, r2
    real(refReal)                     :: u, smin, smax, increment
    integer, allocatable              :: k(:), k_tmp(:)
    character(100), parameter         :: Here = 'number_of_blocks (wrightFisher_func.f90)'

    m = 0
    theta = 1.0

    if (t < 0.07) then
        mu    = 2.0 / t
        sigma = sqrt(2.0 / (3.0 * t))
        r1 = rand % get()
        r2 = rand % get()
        m = nint(mu + sigma * sampleNormalBoxMuller(r1,r2))  
    else
        ! allocate array
        allocate(k(1))
        k = 0
        proceed = .true.
        u = rand % get()
        ! u shouldn't be zero
        u = max(u, tiny(u)) 
        smin = 0.0
        smax = 0.0
        increment = 0

        do while (proceed)
            if (size(k)<(m+1)) then
                ! Need to increase size of array
                ! Double k in size
                allocate(k_tmp(2*m))
                k_tmp(1:size(k)) = k
                k_tmp(size(k)+1:) = 0
                call move_alloc(k_tmp, k)
            end if
            k(m+1) = ceiling(c(m,t,theta)/2)
            do i = 0, k(m+1)-1
                increment = b(m+2*i,m,t,theta) - b(m+2*i+1,m,t,theta)
                smin = smin + increment
                smax = smax + increment
            end do
      
            increment = b(m+2*k(m+1),m,t,theta)
            smin = smin + increment - b(m+2*k(m+1)+1,m,t,theta)
            smax = smax + increment

            do while (smin < u .and. u < smax)
                do i = 0, m
                k(i+1) = k(i+1) + 1
                increment = b(i+2*k(i+1),i,t,theta)
                smax = smin + increment
                smin = smin + increment - b(i+2*k(i+1) + 1,i,t,theta)
                end do
            end do

            if (smin > u) then
                proceed = .false.
            else
                m = m + 1
            end if

        end do

  end function number_of_blocks

  function get_beta(a, b, rand) result(x)
      real(defReal), intent(in)       :: a, b
      class(RNG), intent(intout)      :: rand
      real(defReal)                   :: ga, gb

      ga = get_gamma(a, rand)
      gb = get_gamma(b, rand)
      x = ga/(ga+gb)
  end function get_beta


  !! Only for sampling gamma(a,1)
  function get_gamma(a, rand) result(x)
      real(defReal), intent(in)       :: a
      class(RNG), intent(inout)       :: rand
      real(defReal)                   :: x, c, d, y, r1, r2, v, tmp

      d = a - 1/3
      c = 1/(3*sqrt(d))
      do 
          r1 = rand % get()
          r2 = rand % get()
          y = sampleNormalBoxMuller(r1,r2)
          v = 1+c*y
          if v<0: cycle

          v = v**3
          r1 = rand % get()
          tmp = 1-0.0331*y**4
          if (r1 < tmp) exit
          tmp = 0.5*y**2+d*(1-v*log(v))
          if (log(r1)<tmp) exit
      end do
      x = d*v
  end function get_gamma


  function wright_fisher(t,rand) result(x)
      real(defReal), intent(in)         :: t, r1, r2
      class(RNG), intent(inout)         :: rand
      real(refReal)                     :: x, y
      character(100), parameter         :: Here = 'wright_fisher (wrightFisher_func.f90)'

      if (t > 1E-9) then
          num_blocks = number_of_blocks(t, rand)
          x = get_beta(1,1+m, rand)
      else 
          x = r/2
          r1 = rand % get()
          r2 = rand % get()
          y = sampleNormalBoxMuller(r1,r2)
          x = abs(y*sqrt(t * x * (1 - x)))

  end function wright_fisher
















end module wrightFisher_func
