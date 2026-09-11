!!
!! Transport operator for surface tracking
!!
module transportOperatorSBM_class
  use numPrecision
  use universalVariables

  use errors_mod,                 only : fatalError
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use wrightFisher_func,          only : wrightFisher, getSBMDirection

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,             only : geometry, distCache

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with surface tracking
  !!
  !! Sample Input Dictionary:
  !!   trans { type transportOperatorST; cache 0;}
  !!
  type, public, extends(transportOperator) :: transportOperatorSBM
    logical(defBool)  :: cache = .true.
    real(defReal)     :: dx = 0.05 !TODO: Sort default, and init
    real(defReal)     :: minE = 0.05
  contains
    procedure :: transit => sphericalBrownianMotion
    ! Override procedure
    procedure :: init
  end type transportOperatorSBM

contains


  !!
  !! Performs spherical brownian motion by increment dx, or until surface is hit
  !!
  subroutine sphericalBrownianMotion(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorSBM), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    integer(shortInt)                         :: event
    real(defReal)                             :: sigmaT, dist, sigma, rad_wr, eLoss, sigmaTrack, distST, invSigmaTrack
    real(defReal), dimension(3)               :: dirOut
    type(distCache)                           :: cache
    real(defReal), parameter                  :: tol  = 1.0E-12
    character(100), parameter :: Here = 'sphericalBrownianMotion (transportOperatorST_class.f90)'

    STLoop: do
      
     ! ST DIST CALCULATION
      sigmaTrack = self % xsData % getTrackingXS(p, p % matIdx(), MATERIAL_XS)

      ! Obtain the local cross-section, depending on the material
      ! This branch is called in the case of voids with no imposed XS
      if (sigmaTrack < tol) then

        distST = INFINITY
        invSigmaTrack = INFINITY
        sigmaT = ZERO

      else

        invSigmaTrack = ONE / sigmaTrack
        distST = -log( p % pRNG % get()) * invSigmaTrack

        ! Obtain the local cross-section
        sigmaT = self % xsData % getTrackMatXS(p, p % matIdx())

        ! Should never happen! Catches NaN distances
        if (distST /= distST) call fatalError(Here, "Distance is NaN")

      end if
     ! END OF ST DIST CALCULATION



      ! Save state before movement
      call p % savePrePath()
      dist = self % dx
      ! Move to the next stop.
      if (self % cache) then
        call self % geom % move_withCache(p % coords, dist, event, cache)

      else
        call self % geom % move(p % coords, dist, event)

      end if
      ! dist is now the actual distance travelled
      ! Particle has been moved, but direction of particle has yet to change
      ! Need to sample SBM to get angle change
      sigma  = self % xsData % getMoliereVolatility(p, p % prePath % matIdx, dist) 
      rad_wr = wrightFisher(sigma**2,p % pRNG) 
      !Get new direction
      call getSBMDirection(p % dirGlobal(), dirOut, rad_wr, p % pRNG)
      call p % point(dirOut)
      !Do energy update given travel distance dist
      eLoss = self % xsData % getEnergyLoss(p, p % prePath % matIdx, dist)
      ! TODO: Remove hard coded minimum energy 0.05MeV
      if (eLoss .gt. (p % E - self % minE)) then
          p % E = 0
          p % isDead = .true.
          p % fate = ABS_FATE
          call tally % reportPath(p, dist)
      else 
          p % E = p % E - eLoss
          ! Send tally report for a path moved
          call tally % reportPath(p, dist)

          select case(p % matIdx())
              
            ! Kill particle if it has leaked
            case(OUTSIDE_FILL)
                p % isDead = .true.
                p % fate = LEAK_FATE

            ! Give error if the particle somehow ended in an undefined material
            case(UNDEF_MAT)
                print *, "Particle location: ", p % rGlobal()
                call fatalError(Here, "Particle is in undefined material")
                
            ! Give error if the particle is in a region with overlapping cells
            case(OVERLAP_MAT)
                print *, "Particle location: ", p % rGlobal()
                call fatalError(Here, "Particle is in overlapping cells")
                
            case default
                ! All is well

            end select
      end if

      if (p % isDead) exit STLoop

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real, report collision if virtual
      if (event == COLL_EV) then
        if ((p % pRNG % get() < sigmaT*invSigmaTrack) .and. (p % pRNG % get() <1-exp(-sigmaT*dist))) then
          exit STLoop
        else
          call tally % reportInColl(p, .true.)
        end if
      end if

    end do STLoop

    call tally % reportTrans(p)

  end subroutine sphericalBrownianMotion

  !!
  !! Initialise ST operator from a dictionary
  !!
  !! See transportOperator_inter for details
  !!
  subroutine init(self, dict)
    class(transportOperatorSBM), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Initialise superclass
    call init_super(self, dict)

    if (dict % isPresent('cache')) then
      call dict % get(self % cache, 'cache')
    end if
    if (dict % isPresent('minEnergy')) then
      call dict % get(self % minE, 'minEnergy')
    end if
    if (dict % isPresent('dx')) then
      call dict % get(self % dx, 'dx')
    end if

  end subroutine init

end module transportOperatorSBM_class
