! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!-----------------------------------------------------------------------!
!                                                                       !
! Spherical-wave module                                                 !
!                                                                       !
!  Peter Haynes, 22-Sep-2004                                            !
!  Modified by Mark Robinson, Dec-2007                                  !
!  Modified by Quintin Hill, July 2008 (minor)                          !
!  Gradient routines by Nicholas Hine, August 2010                      !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
! This module provides public routines for calculating the following    !
!  quantities involving the spherical-wave basis set:                   !
!                                                                       !
! * The value of a given basis function at a particular point in space  !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
! The use of this module should be acknowledged in all code and         !
!  resulting publications by referencing:                               !
!                                                                       !
! "Localised spherical-wave basis set for O(N) total-energy             !
!  pseudopotential calculations", P.D.Haynes and M.C.Payne,             !
!  Comput. Phys. Commun. 102 (1997) 17-27.                              !
!                                                                       !
!-----------------------------------------------------------------------!


module spherical_wave

  use constants, only: DP

  implicit none
  private

  ! Public variables

  logical, public                    :: sw_initialised = .false.
  real(kind=DP), public, allocatable :: sw_bessel_zeros(:,:)

  ! Public subroutines

  public :: sw_init
  public :: sw_eval_harmonic
  public :: sw_legendre_fast
  public :: sw_bessel_fast,sw_deriv_bessel_fast
  public :: sw_bessel_accurate
  public :: sw_bessel_zeros_init !qoh
  public :: sw_exit !qoh
  public :: sw_real_sph_harm !qoh
  public :: sw_grad_real_sph_harm !ndmh
  public :: sw_real_sph_harm_unit !qoh
  public :: sw_recp_generate !ars
  public :: sw_recp_generate_in_tb !ars
 ! public :: sw_ngwfs_swcoeff !ars

  ! Private parameters

  real(kind=DP), parameter :: eps = epsilon(1.0_dp)
  real(kind=DP), parameter :: tny = 1.0e-100_dp

  ! Private variables and workspace

  integer :: lmax_init = -1
  integer, save :: num_bessel_points            ! num of pts for Bessel spline
  integer :: num_legendre_points                ! num of pts for Legendre spline
  real(kind=DP), save :: bessel_inv_spacing     ! inv spacing for Bessel spline
  real(kind=DP), save :: legendre_inv_spacing   ! inv spacing for Legendre spl
  real(kind=DP), save :: trig_inv_spacing       ! inv spacing for trig spline
  real(kind=DP), allocatable :: bessel_spline(:,:,:) ! spline fit for Bessel fn
  real(kind=DP), allocatable :: legendre_spline(:,:,:) ! spline fit for Legendre fn
  real(kind=DP), allocatable :: trig_spline(:,:,:)   ! spline fit for trig fns
  real(kind=DP), allocatable :: harmonic_norm(:)     ! normalisation for sph harm

contains

  !-----------------------------------------------------------------------!
  ! sw_init                                                               !
  !                                                                       !
  ! This subroutine must be called to initialise the module before any    !
  !  other routines are called.                                           !
  !                                                                       !
  !  subroutine sw_init(lmax,nmax,rmaxsq)                                 !
  !    integer, intent(in)  :: lmax  - maximum angular momentum component !
  !    integer, intent(in)  :: nmax  - maximum radial quantum number      !
  !    real(kind=DP), intent(in) :: rmaxsq - overides default max xsq     !
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_init(lmax,nmax,rmaxsq)

    use comms, only: comms_abort
    use constants, only: stdout, DP
    use timer, only: timer_clock
    use utils, only: utils_dealloc_check !qoh

    implicit none

    ! Arguments
    integer, intent(in)            :: lmax   ! maximum angular momentum component
    integer, intent(in)            :: nmax   ! maximum radial quantum number
    real(kind=DP), optional, intent(in) :: rmaxsq ! overides default max xsq

    ! Local variables
    integer       :: ierr               ! error flag
    integer       :: nmax_init = -1
    real(kind=DP) :: rmaxsq_init = 0.0_DP
    real(kind=DP) :: xmaxsq    ! Maximum squared arguement
    logical       :: largerxmaxsq ! Is new xmaxsq larger

    call timer_clock('sw_init',1)

    ! Check arguments
    if (lmax < 0) then
       write(stdout,'(a)') 'Error in sw_init: lmax < 0'
       call comms_abort
    end if
    if (nmax < 1) then
       write(stdout,'(a)') 'Error in sw_init: nmax < 1'
       call comms_abort
    end if

    ! qoh: Check for larger rmaxsq
    if (present(rmaxsq)) then
       largerxmaxsq = (rmaxsq > rmaxsq_init )
    else
       largerxmaxsq = .false.
    end if

    ! if basis size larger than that of a previous initialisation deallocate
    if (sw_initialised .and. &
       (lmax > lmax_init .or. nmax > nmax_init .or. largerxmaxsq)) then
       deallocate(sw_bessel_zeros,stat=ierr)
       call utils_dealloc_check('sw_init','sw_bessel_zeros',ierr)
       deallocate(bessel_spline,stat=ierr)
       call utils_dealloc_check('sw_init','bessel_spline',ierr)
       deallocate(legendre_spline,stat=ierr)
       call utils_dealloc_check('sw_init','legendre_spline',ierr)
       deallocate(trig_spline,stat=ierr)
       call utils_dealloc_check('sw_init','trig_spline',ierr)
       deallocate(harmonic_norm,stat=ierr)
       call utils_dealloc_check('sw_init','harmonic_norm',ierr)
       sw_initialised = .false.
    end if

    ! intialise
    if (.not. sw_initialised) then
       ! Calculate spherical Bessel function zeros
       call sw_bessel_zeros_init(nmax,max(lmax,1))

       ! Set up interpolation scheme for spherical Bessel functions
       if (present(rmaxsq)) then
          call sw_bessel_init(nmax,lmax,xmaxsq,rmaxsq)
       else
          call sw_bessel_init(nmax,lmax,xmaxsq)
       endif

       ! Set up interpolation scheme for Legendre functions
       call sw_legendre_init(lmax)

       ! Set up interpolation scheme for trigonometric functions
       call sw_trig_init(lmax)

       ! Set up normalisation factors for spherical harmonics
       call sw_harmonic_norm_init(lmax)

       ! Successful initialisation
       sw_initialised = .true.
       lmax_init = lmax
       nmax_init = nmax
       rmaxsq_init = xmaxsq
    end if

    call timer_clock('sw_init',2)

  end subroutine sw_init


  !-----------------------------------------------------------------------!
  ! sw_eval_harmonic                                                      !
  !                                                                       !
  ! Evaluates a spherical harmonic for a given l and m                    !
  !                                                                       !
  !  subroutine sw_eval_recip(l,m,pos,value)                              !
  !    integer, intent(in)     :: l         <- angular momentum           !
  !    integer, intent(in)     :: m         <- azimuthal angular momentum !
  !    type(POINT), intent(in) :: pos       <- position for evaluation    !
  !-----------------------------------------------------------------------!

  real(kind=DP) function sw_eval_harmonic(l,m,pos)

    use comms, only: comms_abort
    use constants, only: stdout, DP
    use geometry, only: point

    implicit none

    ! Arguments
    integer, intent(in)     :: l,m       ! quantum numbers
    type(POINT), intent(in) :: pos       ! position where values required

    ! Local variables
    real(kind=DP) :: z2,xy,xyz2
    real(kind=DP) :: cos_phi,sin_phi,cos_theta2
    integer  :: absm,lm
    real(kind=DP) :: legendre,trig
    complex(kind=DP) :: e_imphi

    ! Check arguments
    if (l > lmax_init) then
        write(stdout,'(a)') 'Error in sw_eval_harmonic: l too high'
        call comms_abort
    end if

    ! Find position info
    xy = pos%X**2 + pos%Y**2
    z2 = pos%Z**2
    xyz2 = xy + z2
    xy = sqrt(xy)

    ! Calculate phi and theta
    if (xy > eps) then
       cos_phi = pos%X / xy
       sin_phi = pos%Y / xy
    else
       cos_phi = 1.0_DP
       sin_phi = 0.0_DP
    end if
    if (xyz2 > eps) then
       cos_theta2 = z2 / xyz2
    else
       cos_theta2 = 1.0_DP
    end if

    ! Calculate associated legendre polynomial
    absm = abs(m)
    call sw_legendre_fast(l,absm,cos_theta2,legendre)
    if(pos%Z < 0) legendre = legendre * (-1)**(l+absm)

    ! Calculate trig functions
    ! should use spline fits?
!    if (m < 0) then
!       call sw_sin_fast(absm,sin_phi2,trig)
!       if (pos%Y < 0) trig = trig * (-1)**m
!    else if (m > 0) then
!       call sw_cos_fast(absm,cos_phi2,trig)
!       if (pos%X < 0) trig = trig * (-1)**m
!    else
!       trig = 1.0_DP
!    end if

    e_imphi = cmplx(cos_phi, sin_phi, DP) **m
    if (m < 0) then
       trig = -aimag(e_imphi)
    else if (m > 0) then
       trig = real(e_imphi)
    else
       trig = 1.0_DP
    endif

    ! Calculate spherical harmonic
    lm = absm + l*(l+1)/2 +1

    sw_eval_harmonic = harmonic_norm(lm) * legendre * trig

  end function sw_eval_harmonic


  !-----------------------------------------------------------------------!
  ! sw_bessel_accurate                                                    !
  !                                                                       !
  ! This subroutine returns the values of the spherical Bessel function   !
  ! of the first kind and its derivative for a given order l and argument !
  ! x.                                                                    !
  !                                                                       !
  !  subroutine sw_bessel_accurate(l,x,j,jp)                              !
  !    integer, intent(in)        :: l  - the order (l >= 0)              !
  !    real(kind=DP), intent(in)  :: x  - the argument (x >= 0)           !
  !    real(kind=DP), intent(out) :: j  - sBf of first kind               !
  !    real(kind=DP), intent(out) :: jp - derivative of sBf of first kind !
  !                                                                       !
  ! On exit:  j contains the spherical Bessel function of the first kind  !
  !             of order l and argument x                                 !
  !           jp contains the first derivative of the spherical Bessel    !
  !             function of the first kind of order l and argument x      !
  !                                                                       !
  ! This subroutine uses Steed's method adapted for spherical Bessel      !
  ! functions, which solves four simultaneous equations for the functions !
  ! and their derivatives. The starting point is the evaluation of a      !
  ! continued fraction by the modified Lentz method. For more details,    !
  ! see Numerical Recipes in Fortran77 sec. 6.7 and references therein.   !
  !-----------------------------------------------------------------------!

  subroutine sw_bessel_accurate(l,x,j,jp)

    use comms, only: comms_abort
    use constants, only: stdout, DP, pi

    implicit none

    ! Arguments

    integer, intent(in)   :: l    ! The order
    real(kind=DP), intent(in)  :: x    ! The argument
    real(kind=DP), intent(out) :: j    ! Spherical Bessel function of the first kind
    real(kind=DP), intent(out) :: jp   ! Derivative of the first kind

    ! Local variables

    integer, parameter :: max_iter = 10000 ! The maximum number of iterations
                                           ! attempted to calculate CF1
    integer :: k                           ! Loop variable
    integer :: isign                       ! Sign
    real(kind=DP), parameter :: halfpi = 0.5_dp * pi
    real(kind=DP) :: xinv                       ! 1/x
    real(kind=DP) :: factor                     ! Temporary multiplicative factor
    real(kind=DP) :: cf1,b,c,d,delta            ! Variables used in calculation of
                                           ! continued fraction (CF1)
    real(kind=DP) :: rj,rjp,rj0,rjp0,rjnew      ! Variables used in downward
                                           ! recurrence of first kind
    real(kind=DP) :: f,w,gamma,j0!/,jp0           ! Variables used in solution of
                                           ! simultaneous equations

    ! Check arguments
    if (l < 0) then
       write(stdout,'(a)') 'Error in spherical_bessel: l < 0'
       call comms_abort
    end if
    if (x < 0.0_dp) then
       write(stdout,'(a)') 'Error in spherical_bessel: x < 0'
       call comms_abort
    end if

    ! Special case when x=0 (to avoid division by zero)
    if (x < eps) then
       j = 0.0_dp
       jp = 0.0_dp
       if (l == 0) j = 1.0_dp
       if (l == 1) jp = 1.0_dp / 3.0_dp
       return
    end if
    xinv = 1.0_dp / x

    ! Evaluate continued fraction CF1 by Lentz's method
    !
    !                       l       1          1
    ! CF1 = j'(x) / j (x) = - - ---------- ---------- ...
    !        l       l      x   (2l+3)/x - (2l+5)/x -

    delta = 0.0_dp
    ! qoh: Always initialise cf1 and isign, otherwise they could be used below
    ! qoh: uninitialised
    cf1 = max(l * xinv, tny)
    isign = 1
    if (x < real(max_iter + l,dp)) then ! Do not attempt for very large x
       c = cf1
       d = 0.0_dp
       b = (2*l+1) * xinv
       do k=1,max_iter
          b = b + 2.0_dp * xinv
          d = b - d
          if (abs(d) < tny) d = tny
          c = b - 1.0_dp / c
          if (abs(c) < tny) c = tny
          d = 1.0_dp / d
          delta = c * d
          cf1 = delta * cf1
          if (d < 0.0_dp) isign = -isign
          if (abs(delta - 1.0_dp) < eps) exit
       end do
    end if

    ! If CF1 does not converge (or not attempted), use asymptotic expansions:
    !
    ! j (z) ~ cos[x - (l+1) pi/2] / x
    !  l

    if (abs(delta - 1.0_dp) >= eps) then

       j = cos(x - halfpi * (l+1)) * xinv
       jp = -sin(x - halfpi * (l+1)) * xinv
       write(stdout,'(a)') 'Warning in sw_bessel_accurate: &
            &asymptotic expansion used'

    else ! Otherwise use Steed's method

       ! Recur j and j' downwards to order 0
       !
       !           l+1
       ! j   (z) = --- j (z) + j'(z)
       !  l-1       z   l       l
       !
       !           l-1
       ! j'  (z) = --- j   (z) - j (z)
       !  l-1       z   l-1       l

       rj = isign * tny
       rjp = cf1 * rj

       rj0 = rj
       rjp0 = rjp

       factor = (l+1) * xinv
       do k=l,1,-1
          rjnew = factor * rj + rjp
          factor = factor - xinv
          rjp = (factor - xinv) * rjnew - rj
          rj = rjnew
       end do

       ! Continued fraction CF1 for l=0:

       if (rj == 0.0_dp) rj = eps
       f = rjp / rj

       ! Continued fraction CF2 is trivially evaluated for l=0:
       !
       !                                                   1
       ! CF2 = [j'(x) + i y'(x)] / [j (x) + i y (z)] = i - -
       !         0         0         0         0           x

       ! Wronskian gives fourth simultaneous equation:
       !                                  -2
       ! W = j (x) y'(x) - y (x) j'(x) = x
       !      l     l       l     l

       w = xinv * xinv

       ! Solve four equations linking j (x), j'(x), y (x) and y'(x)
       !                               0      0      0         0

       gamma = -xinv - f
       j0 = sign(sqrt(w / (1.0_dp + gamma * gamma)),rj)
       !jp0 = f * j0 !qoh: unused

       ! Rescale j and j'

       factor = j0 / rj
       j = rj0 * factor
       jp = rjp0 * factor

    end if

  end subroutine sw_bessel_accurate


  !-----------------------------------------------------------------------!
  ! sw_bessel_init                                                        !
  !                                                                       !
  ! This subroutine fits cubic splines to the spherical Bessel function   !
  ! of the first kind and its derivative, for use when the square of the  !
  ! argument is passed.                                                   !
  !                                                                       !
  !  subroutine sw_bessel_init(nmax,lmax,rmaxsq)                          !
  !    integer, intent(in)   :: nmax   - maximum radial quantum number    !
  !    integer, intent(in)   :: lmax   - maximum angular momentum         !
  !    real(kind=DP), intent(in)  :: rmaxsq - maximum distance from centre!
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_bessel_init(nmax,lmax_in,xmaxsq,rmaxsq)

    use utils, only: utils_alloc_check !qoh

    implicit none

    ! Arguments

    integer, intent(in)            :: nmax    ! maximum radial quantum number
    integer, intent(in)            :: lmax_in ! maximum angular momentum
    real(kind=DP), intent(out)     :: xmaxsq  ! maximum argument squared
    real(kind=DP), optional, intent(in) :: rmaxsq ! maximum xsq to be included

    ! Local variables

    real(kind=DP), parameter :: tol = 1.0e-9_dp   ! minimum required accuracy
    real(kind=DP), parameter :: third = 1.0_dp / 3.0_dp
    integer :: ierr                          ! error flag
    integer :: l                             ! angular momentum
    integer :: ipoint                        ! point loop counter
    integer :: lmax                          ! maximum angular momentum
    real(kind=DP) :: xmax                         ! maximum argument
    real(kind=DP) :: x,xsq                        ! (squared) argument
    real(kind=DP) :: j,jp                         ! sBf and first derivative
    real(kind=DP) :: spacing                      ! spacing

    ! qoh: We usually need l=-1 spherical bessel when generating spherical waves
    ! qoh: so we need at least l=1 spherical bessels.
    lmax = max(lmax_in,1)

    ! Set up uniform spacing for spline
    xmax = sw_bessel_zeros(nmax,lmax)
    xmaxsq = xmax*xmax
    if (present(rmaxsq)) xmaxsq = max (xmaxsq,rmaxsq)
    bessel_inv_spacing = 1.0_dp / tol**third
    num_bessel_points = xmaxsq * bessel_inv_spacing + 2
    spacing = 1.0_dp / bessel_inv_spacing

    ! Allocate module workspace
    allocate(bessel_spline(num_bessel_points,2,0:2*lmax+1),stat=ierr)
    call utils_alloc_check('sw_bessel_init','bessel_spline',ierr) !qoh

    ! Perform cubic spline fitting
    do l=0,lmax
       do ipoint=1,num_bessel_points
          xsq = (ipoint - 1) * spacing
          x = sqrt(xsq)
          call sw_bessel_accurate(l,x,j,jp)
          bessel_spline(ipoint,1,2*l) = j
          bessel_spline(ipoint,1,2*l+1) = jp
       end do
    end do
    call sw_spline_fit(num_bessel_points,2*(lmax+1),bessel_spline)

  end subroutine sw_bessel_init


  !-----------------------------------------------------------------------!
  ! sw_harmonic_norm_init                                                 !
  !                                                                       !
  ! This subroutine evaluates normalisation factors for the spherical     !
  ! harmonics.                                                            !
  !                                                                       !
  !  subroutine sw_harmonic_norm_init(lmax)                               !
  !    integer, intent(in)   :: lmax  - maximum angular momentum          !
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_harmonic_norm_init(lmax)

    use constants, only: DP, PI
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments

    integer, intent(in) :: lmax  ! maximum angular momentum

    ! Local variables

    integer :: ierr           ! Error flag
    integer :: lmmax          ! Number of functions
    integer :: l              ! Angular momentum loop counter
    integer :: m              ! Azimuthal quantum number loop counter
    integer :: lm             ! Function loop counter
    real(kind=DP) :: factor   ! Quotient of factorials
    real(kind=DP) :: prefactor     ! (2l+1)/(2 pi)

    ! Number of functions to normalise
    lmmax = ((lmax+1)*(lmax+2))/2

    ! Allocate module workspace
    allocate(harmonic_norm(lmmax),stat=ierr)
    call utils_alloc_check('sw_harmonic_norm_init','harmonic_norm',ierr) !qoh

    ! Calculate normalisation factors
    lm = 1
    do l=0,lmax
       prefactor = (2*l+1) / (2.0_dp * pi)
       harmonic_norm(lm) = sqrt(0.5_dp * prefactor)
       lm = lm + 1
       factor = 1.0_DP
       do m=1,l
          factor = factor * real((l+m) * (l-m+1),kind=DP)
          harmonic_norm(lm) = sqrt(prefactor / factor)
          lm = lm + 1
       end do
    end do

  end subroutine sw_harmonic_norm_init


  !-----------------------------------------------------------------------!
  ! sw_bessel_fast                                                        !
  !                                                                       !
  ! This subroutine evaluates cubic spline fits to the spherical Bessel   !
  ! function of the first kind, for use when the square of the argument   !
  ! is passed.                                                            !
  !                                                                       !
  !  real(kind=DP) function sw_bessel_fast(li,xsq)                        !
  !    integer, intent(in)   :: li   - angular momentum                   !
  !    real(kind=DP), intent(in)  :: xsq  - squared argument              !
  !                                                                       !
  !-----------------------------------------------------------------------!

  real(kind=DP) function sw_bessel_fast(li,xsq)

    use comms, only: comms_abort
    use constants, only: DP, stdout

    implicit none

    ! Arguments

    integer, intent(in)  :: li   ! Angular momentum (order)
    real(kind=DP), intent(in) :: xsq  ! Squared argument

    ! Local variables

    integer :: l      ! abs(angular momentum)
    integer :: ipos   ! Largest spline point less than argument
    integer :: twol   ! Twice l
    real(kind=DP) :: pos   ! Argument in units of spline spacing
    real(kind=DP) :: a,b   ! Interpolation factors

    l = abs(li)
    twol = 2 * l
    pos = xsq * bessel_inv_spacing
    ipos = int(pos)
    b = pos - real(ipos,dp)
    a = 1.0_dp - b

    ! check bessel_spline is large enough for this point
    if (ipos+2 > num_bessel_points) then
       write(stdout,'(a,i10,a,i10)') 'Error in sw_bessel_fast: &
            &num_bessel_points too small: ',ipos,' > ',num_bessel_points
       call comms_abort
    endif

    sw_bessel_fast = a * bessel_spline(ipos+1,1,twol) + &
         b * bessel_spline(ipos+2,1,twol) + &
         (a*a*a-a) * bessel_spline(ipos+1,2,twol) + &
         (b*b*b-b) * bessel_spline(ipos+2,2,twol)
    if (li < 0) sw_bessel_fast = sw_bessel_fast * (-1)**l

  end function sw_bessel_fast


  !-----------------------------------------------------------------------!
  ! sw_deriv_bessel_fast                                                  !
  !                                                                       !
  ! This subroutine evaluates cubic spline fits to the first derivative   !
  ! of the spherical Bessel function of the first kind, for use when the  !
  ! square of the argument is passed.                                     !
  !                                                                       !
  !  real(kind=DP) function sw_deriv_bessel_fast(li,xsq)                  !
  !    integer, intent(in)   :: li   - angular momentum                   !
  !    real(kind=DP), intent(in)  :: xsq  - squared argument              !
  !                                                                       !
  !-----------------------------------------------------------------------!

  real(kind=DP) function sw_deriv_bessel_fast(li,xsq)

    implicit none

    ! Arguments

    integer, intent(in)  :: li   ! Angular momentum (order)
    real(kind=DP), intent(in) :: xsq  ! Squared argument

    ! Local variables

    integer :: l      ! abs(angular momentum)
    integer :: ipos   ! Largest spline point less than argument
    integer :: twolp1 ! Twice l plus one
    real(kind=DP) :: pos   ! Argument in units of spline spacing
    real(kind=DP) :: a,b   ! Interpolation factors

    l = abs(li)
    twolp1 = 2 * l + 1
    pos = xsq * bessel_inv_spacing
    ipos = int(pos)
    b = pos - real(ipos,dp)
    a = 1.0_dp - b
    sw_deriv_bessel_fast = a * bessel_spline(ipos+1,1,twolp1) + &
         b * bessel_spline(ipos+2,1,twolp1) + &
         (a*a*a-a) * bessel_spline(ipos+1,2,twolp1) + &
         (b*b*b-b) * bessel_spline(ipos+2,2,twolp1)
    if (li < 0) sw_deriv_bessel_fast = sw_deriv_bessel_fast * (-1)**l

  end function sw_deriv_bessel_fast


  !-----------------------------------------------------------------------!
  ! sw_legendre_fast                                                      !
  !                                                                       !
  ! This subroutine evaluates cubic spline fits to the Legendre           !
  ! function, for use when the square of the argument is passed.          !
  !                                                                       !
  ! Fixed by Jacek Dziedzic in 2011.09 so that it doesn't become imprecise!
  ! for corner cases.                                                     !
  !                                                                       !
  !  subroutine sw_legendre_fast(l,m,xsq,plm)                             !
  !    integer, intent(in)        :: l,m    - quantum numbers needed      !
  !    real(kind=DP), intent(in)  :: xsq    - squared argument            !
  !    real(kind=DP), intent(out) :: plm    - results                     !
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_legendre_fast(l,m,xsq,plm)

    implicit none

    ! Arguments

    integer, intent(in)   :: l,m    ! Quantum numbers
    real(kind=DP), intent(in)  :: xsq    ! Squared argument
    real(kind=DP), intent(out) :: plm    ! Results

    ! Local variables

    integer :: lm          ! Function loop counter
    integer :: ipos        ! Largest spline point less than argument
    real(kind=DP) :: pos        ! Argument in units of spline spacing
    real(kind=DP) :: a,b,f1,f2  ! Interpolation factors

    ! jd: First take care of corner cases
    ! : 1.0 (+/- rounding error)
    if(xsq >= 1.0_DP-eps) then
       if(m==0) then
          plm = 1.0_DP
       else
          plm = 0.0_DP
       end if
       return
    end if
    ! : very close to 1.0, where the spline interpolation rapidly becomes
    !   inaccurate (for certain {l, m}'s) -- use the accurate approach.
    if(xsq >= 0.999_DP) then
       call sw_legendre_accurate(l,m,sqrt(xsq),plm)
       return
    end if
    ! : very close to 0.0, where the spline interpolation rapidly becomes
    !   inaccurate (for remaining {l, m}'s) -- use the accurate approach.
    if(xsq <= 0.001_DP) then
       call sw_legendre_accurate(l,m,sqrt(xsq),plm)
       return
    end if

    ! jd: Now deal with the usual cases

    ! Work out number of functions required
    lm = m + l*(l+1)/2 +1

    ! Work out position in spline grid and related interpolation factors
    pos = xsq * legendre_inv_spacing
    ipos = int(pos)
    b = pos - real(ipos,dp)

    a = 1.0_dp - b
    f1 = a*a*a - a
    f2 = b*b*b - b

    ! Evaluate splines
    plm = a * legendre_spline(ipos+1,1,lm) + &
          b * legendre_spline(ipos+2,1,lm) + &
          f1 * legendre_spline(ipos+1,2,lm) + &
          f2 * legendre_spline(ipos+2,2,lm)

  end subroutine sw_legendre_fast


  !-----------------------------------------------------------------------!
  ! sw_cos_fast                                                           !
  !                                                                       !
  ! This subroutine evaluates cubic spline fits to the cos trigonometric  !
  ! function, for use when the square of the argument is passed.          !
  !                                                                       !
  !  subroutine sw_cos_fast(m,xsq,cosmphi)                                !
  !    integer, intent(in)        :: m       - quantum numbers            !
  !    real(kind=DP), intent(in)  :: xsq     - squared argument           !
  !    real(kind=DP), intent(out) :: cosmphi - results                    !
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_cos_fast(m,xsq,cosmphi)

    implicit none

    ! Arguments

    integer, intent(in)   :: m           ! Max azimuthal quantum number
    real(kind=DP), intent(in)  :: xsq         ! Squared argument
    real(kind=DP), intent(out) :: cosmphi     ! Results

    ! Local variables

    integer :: twomm1      ! Twice m minus one
    integer :: ipos        ! Largest spline point less than argument
    real(kind=DP) :: pos        ! Argument in units of spline spacing
    real(kind=DP) :: a,b,f1,f2  ! Interpolation factors

    ! Work out position in spline grid and related interpolation factors
    pos = xsq * trig_inv_spacing
    ipos = int(pos)
    b = pos - real(ipos,dp)
    a = 1.0_dp - b
    f1 = a*a*a - a
    f2 = b*b*b - b

    ! Evaluate splines
    twomm1 = 2*m-1
    cosmphi = a * trig_spline(ipos+1,1,twomm1) + &
              b * trig_spline(ipos+2,1,twomm1) + &
              f1 * trig_spline(ipos+1,2,twomm1) + &
              f2 * trig_spline(ipos+2,2,twomm1)

  end subroutine sw_cos_fast


  !-----------------------------------------------------------------------!
  ! sw_sin_fast                                                           !
  !                                                                       !
  ! This subroutine evaluates cubic spline fits to the sin trigonometric  !
  ! function, for use when the square of the argument is passed.          !
  !                                                                       !
  !  subroutine sw_sin_fast(m,xsq,sinmphi)                                !
  !    integer, intent(in)   :: m       - max azimuthal quantum number    !
  !    real(kind=DP), intent(in)  :: xsq     - squared argument           !
  !    real(kind=DP), intent(out) :: sinmphi - results                    !
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_sin_fast(m,xsq,sinmphi)

    implicit none

    ! Arguments

    integer, intent(in)   :: m           ! Max azimuthal quantum number
    real(kind=DP), intent(in)  :: xsq         ! Squared argument
    real(kind=DP), intent(out) :: sinmphi     ! Results

    ! Local variables

    integer :: twom        ! Twice m
    integer :: ipos        ! Largest spline point less than argument
    real(kind=DP) :: pos        ! Argument in units of spline spacing
    real(kind=DP) :: a,b,f1,f2  ! Interpolation factors

    ! Work out position in spline grid and related interpolation factors
    pos = xsq * trig_inv_spacing
    ipos = int(pos)
    b = pos - real(ipos,dp)
    a = 1.0_dp - b
    f1 = a*a*a - a
    f2 = b*b*b - b

    ! Evaluate splines
    twom = 2*m
    sinmphi = a * trig_spline(ipos+1,1,twom) + &
              b * trig_spline(ipos+2,1,twom) + &
              f1 * trig_spline(ipos+1,2,twom) + &
              f2 * trig_spline(ipos+2,2,twom)

  end subroutine sw_sin_fast


  !-----------------------------------------------------------------------!
  ! sw_legendre_accurate                                                  !
  !                                                                       !
  ! This subroutine returns the value of the associated Legendre          !
  ! function P_l^m(x).                                                    !
  !                                                                       !
  !  subroutine sw_legendre_accurate(l,m,x,plm)                           !
  !    integer, intent(in)   :: l      - the orbital angular momentum     !
  !    integer, intent(in)   :: m      - the azimuthal quantum number     !
  !    real(kind=DP), intent(in)  :: x      - argument                    !
  !    real(kind=DP), intent(out) :: plm    - the legendre function       !
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_legendre_accurate(l,m,x,plm)

    use comms, only: comms_abort
    use constants, only: DP, stdout

    implicit none

    ! Arguments

    integer, intent(in) :: l        ! angular momentum
    integer, intent(in) :: m        ! azimuthal quantum number
    real(kind=DP), intent(in) :: x       ! argument
    real(kind=DP), intent(out) :: plm    ! legendre function

    ! Local variables

    integer  :: i                   ! loop counter
    real(kind=DP) :: pmm                 ! legendre function P_m^m(x)
    real(kind=DP) :: rho                 ! sqrt(1-x^2)
    real(kind=DP) :: factor              ! factorial
    real(kind=DP) :: pmm1                ! temporary variable for recurrence

    ! Check arguments
    if (l < 0) then
       write(stdout,'(a)') 'Error in sw_legendre_accurate: l < 0'
       call comms_abort
    end if
    if (m < 0) then
       write(stdout,'(a)') 'Error in sw_legendre_accurate: m < 0'
       call comms_abort
    end if
    if (abs(x) > 1.0_dp) then
       write(stdout,'(a)') 'Error in sw_legendre_accurate: |x| > 1'
       call comms_abort
    end if

    !                                         m
    ! Calculate associated Legendre function P (x)
    !                                         m
    pmm = 1.0_dp
    if (m > 0) then
       rho = sqrt((1.0_dp - x)*(1.0_dp + x))
       factor = 1.0_dp
       do i=1,m
          pmm = pmm * factor * rho
          factor = factor + 2.0_dp
       end do
    end if

    !                                         m
    ! Calculate associated Legendre function P (x)
    !                                         l
    if (l == m) then
       plm = pmm
    else
       pmm1 = x * (2*m+1) * pmm
       if (l == m+1) then
          plm = pmm1
       else
          do i=m+2,l
             plm = (x*(2*i-1)*pmm1 - (i+m-1)*pmm) / real(i-m,dp)
             pmm = pmm1
             pmm1 = plm
          end do
       end if
    end if

  end subroutine sw_legendre_accurate


  !-----------------------------------------------------------------------!
  ! sw_legendre_init                                                      !
  !                                                                       !
  ! This subroutine fits cubic splines to the associated Legendre         !
  ! functions,  for use when the square of the argument is passed.        !
  !                                                                       !
  !  subroutine sw_legendre_init(lmax)                                    !
  !    integer, intent(in)   :: lmax  - maximum angular momentum          !
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_legendre_init(lmax)

    use utils, only: utils_alloc_check !qoh

    implicit none

    ! Arguments

    integer, intent(in) :: lmax  ! maximum angular momentum

    ! Local variables

    real(kind=DP), parameter :: tol = 1.0e-9_dp   ! minimum required accuracy
    real(kind=DP), parameter :: third = 1.0_dp / 3.0_dp
    integer :: ierr                          ! error flag
    integer :: l                             ! angular momentum
    integer :: m                             ! azimuthal quantum number
    integer :: ipoint                        ! point loop counter
    integer :: lmmax                         ! number of distinct fns
    integer :: lm                            ! function loop counter
    real(kind=DP) :: x,xsq                        ! (squared) argument
    real(kind=DP) :: plm                          ! Legendre function
    real(kind=DP) :: spacing                      ! spacing

    ! Set up uniform spacing for spline
    legendre_inv_spacing = nint(max(lmax,1)**2 / tol**third)
    num_legendre_points = legendre_inv_spacing + 2
    spacing = 1.0_dp / legendre_inv_spacing

    ! Number of funtions to fit
    lmmax = ((lmax+1) * (lmax+2))/2

    ! Allocate module workspace
    allocate(legendre_spline(num_legendre_points,2,lmmax),stat=ierr)
    call utils_alloc_check('sw_legendre_init','legendre_spline',ierr) !qoh

    ! Perform cubic spline fitting
    lm = 1
    do l=0,lmax
       do m=0,l
          do ipoint=1,num_legendre_points-1
             xsq = (ipoint - 1) * spacing
             x = sqrt(xsq)
             call sw_legendre_accurate(l,m,x,plm)
             legendre_spline(ipoint,1,lm) = plm
          end do
          legendre_spline(num_legendre_points,1,lm) = plm
          lm = lm + 1
       end do
    end do
    call sw_spline_fit(num_legendre_points,lmmax,legendre_spline)

  end subroutine sw_legendre_init


  !-----------------------------------------------------------------------!
  ! sw_trig_init                                                          !
  !                                                                       !
  ! This subroutine fits cubic splines to the trigonometric functions,    !
  ! for use when the square of the argument is passed.                    !
  !                                                                       !
  !  subroutine sw_trig_init(lmax)                                        !
  !    integer, intent(in)   :: lmax  - maximum angular momentum          !
  !                                                                       !
  !-----------------------------------------------------------------------!

  subroutine sw_trig_init(lmax)

    use utils, only: utils_alloc_check

    implicit none

    ! Arguments

    integer, intent(in) :: lmax  ! maximum angular momentum

    ! Local variables

    real(kind=DP), parameter :: tol = 1.0e-9_dp   ! minimum required accuracy
    real(kind=DP), parameter :: third = 1.0_dp / 3.0_dp
    integer :: ierr                          ! error flag
    integer :: m                             ! azimuthal quantum number
    integer :: ipoint                        ! point loop counter
    integer :: num_trig_points                    ! num of pts for trig spline
    real(kind=DP) :: x,xsq                        ! (squared) argument
    real(kind=DP) :: phic,phis                    ! azimuthal angles
    real(kind=DP) :: spacing                      ! spacing

    ! Set up uniform spacing for spline
    trig_inv_spacing = nint(max(lmax,1)**2 / tol**third)
    num_trig_points = trig_inv_spacing + 2
    spacing = 1.0_dp / trig_inv_spacing

    ! Allocate module workspace
    allocate(trig_spline(num_trig_points,2,2*max(lmax,1)),stat=ierr)
    call utils_alloc_check('sw_trig_init','trig_spline',ierr)

    ! jd: Take care of the points untouched by the loop that follows
    trig_spline = 0.0_DP

    ! Calculate function values for spline
    do ipoint=1,num_legendre_points-1
       xsq = (ipoint - 1) * spacing
       x = sqrt(xsq)
       phic = acos(x)
       phis = asin(x)
       do m=1,max(lmax,1)
          trig_spline(ipoint,1,2*m-1) = cos(m*phic)
          trig_spline(ipoint,1,2*m) = sin(m*phis)
       end do
    end do

    ! Perform cubic spline fitting
    call sw_spline_fit(num_trig_points,2*max(lmax,1),trig_spline)

  end subroutine sw_trig_init


  !-----------------------------------------------------------------------!
  ! sw_bessel_zeros_init                                                  !
  !                                                                       !
  ! This subroutine returns the values of the first n zeros of the        !
  ! spherical Bessel functions of the first kind up to order lmax.        !
  !                                                                       !
  !  subroutine sw_bessel_zeros_init(n,lmax)                              !
  !    integer, intent(in)  :: n     - the number of zeros required       !
  !    integer, intent(in)  :: lmax  - the maximum order required         !
  !                                                                       !
  ! On exit:  The results are returned in the double precision module     !
  !           array sw_bessel_zeros(n,0:lmax)                             !
  !                                                                       !
  ! This subroutine uses the Newton-Raphson method to find the zeros of   !
  ! the spherical Bessel functions, using known results to bracket the    !
  ! zeros.                                                                !
  !-----------------------------------------------------------------------!

  subroutine sw_bessel_zeros_init(n,lmax)

    use comms, only: comms_abort
    use constants, only: DP, PI, stdout
    use utils, only: utils_alloc_check, utils_dealloc_check ! qoh

    implicit none !qoh

    ! Arguments

    integer, intent(in)  :: n     ! The number of zeros to be found
    integer, intent(in)  :: lmax  ! The maximum order required

    ! Local variables

    integer, parameter :: max_iter = 100     ! The maximum number of
                                             ! Newton-Raphson iterations
    integer :: k                             ! Loop variable for zero
    integer :: l                             ! Loop variable for order
    integer :: iter                          ! Iteration counter for
                                             ! Newton-Raphson
    integer :: ierr                          ! Error flag
    logical :: toohigh,toolow                ! Flags for root finding
    real(kind=DP) :: z,zold,lower,upper      ! Variables for root finding
    real(kind=DP) :: j,jp                    ! Spherical Bessel values
    real(kind=DP), allocatable :: zeros(:)   ! Temporary array for zeros

    ! Allocate required space for zeros
    if (allocated(sw_bessel_zeros)) then
       if (n > size(sw_bessel_zeros,1) .or. &
            lmax > size(sw_bessel_zeros,2)-1) then
          deallocate(sw_bessel_zeros,stat=ierr)
          call utils_dealloc_check('sw_bessel_zeros_init','sw_bessel_zeros',&
               &ierr)    !qoh
       end if
    end if
    if (.not. allocated(sw_bessel_zeros)) then
       allocate(sw_bessel_zeros(n,0:lmax),stat=ierr)
       call utils_alloc_check('sw_bessel_zeros_init','sw_bessel_zeros',ierr)
       !qoh
    end if

    ! The zeros of Bessel functions interlace according to:
    !
    ! z    < z      < z      < z
    !  l,n    l+1,n    l,n+1    l+1,n+1
    !
    ! Therefore we solve for more than n zeros of the lower order functions
    ! in order to provide bounds for the higher orders - allocate local
    ! array for this:
    allocate(zeros(n+lmax),stat=ierr)
    call utils_alloc_check('sw_bessel_zeros_init','zeros',ierr) !qoh

    ! Zeros of j (x) are trivial: nth zero is z   = n pi
    !           0                              0,n
    !
    do k=1,n+lmax
       zeros(k) = k * pi
    end do
    sw_bessel_zeros(1:n,0) = zeros(1:n)

    ! Set status flag to 0 (success) initially
    ierr = 0

    ! Loop over higher orders (angular momenta) required
    do l=1,lmax

       ! Loop over number of zeros (including extras to provide bounds for
       ! higher orders) required
       do k=1,n+lmax-l

          ! Take bounds from previous set
          lower = zeros(k)
          upper = zeros(k+1)

          ! Initial guess is half-way between bounds
          z = 0.5_dp * (lower + upper)

          ! Set flags
          toohigh = .false.
          toolow = .false.

          ! Start of Newton-Raphson iteration
          do iter=1,max_iter

             ! Get values of j and j' at this z
             call sw_bessel_accurate(l,z,j,jp)

             ! Update position
             if (abs(jp) < tny) then
                write(stdout,'(a)') 'Error in sw_bessel_zeros_init: &
                     &j'' too small'
                call comms_abort
             end if
             zold = z
             z = zold - j / jp

             ! Check bounds
             if (z < lower) then
                if (toolow) then
                   write(stdout,'(a)') 'Error in sw_bessel_zeros_init: &
                        &hit lower bound twice'
                   call comms_abort
                else
                   z = lower
                   toolow = .true.
                end if
             end if
             if (z > upper) then
                if (toohigh) then
                   write(stdout,'(a)') 'Error in sw_bessel_zeros_init: &
                        &hit upper bound twice'
                   call comms_abort
                else
                   z = upper
                   toohigh = .true.
                end if
             end if

             ! Check convergence
             if (abs(j) < eps .and. abs(z-zold) < eps*z) exit

          end do ! Newton-Raphson

          zeros(k) = z

       end do ! Zeros

       sw_bessel_zeros(1:n,l) = zeros(1:n)

    end do ! Angular momentum

    !qoh: Deallocate the internal zeros array

    deallocate(zeros, stat=ierr)
    call utils_dealloc_check('sw_bessel_zeros_init','zeros',ierr)

  end subroutine sw_bessel_zeros_init


  subroutine sw_spline_fit(n,m,y)

    use comms, only: comms_abort
    use constants, only: DP, stdout
    use utils, only: utils_alloc_check, utils_dealloc_check !qoh

    implicit none

    ! Arguments

    integer, intent(in) :: n            ! Length of data set to fit
    integer, intent(in) :: m            ! Number of data sets to fit
    real(kind=DP), intent(inout) :: y(n,2,m) ! The data to fit/second derivative

    ! Local variables

    integer :: i                             ! Point loop counter
    integer :: j                             ! Data set loop counter
    integer :: ierr                          ! Error flag
    real(kind=DP), allocatable :: diag_work(:)    ! Diagonal elements
    real(kind=DP), allocatable :: subdiag_work(:) ! Sub-diagonal elements

    ! Allocate temporary workspace
    allocate(diag_work(n-2),stat=ierr)
    call utils_alloc_check('sw_spline_fit','diag_work',ierr) !qoh

    allocate(subdiag_work(n-3),stat=ierr)
    call utils_alloc_check('sw_spline_fit','subdiag_work',ierr) !qoh

    ! Set up linear equation set to solve
    diag_work = 4.0_dp
    subdiag_work = 1.0_dp
    do j=1,m
       do i=3,n-1
          y(i,2,j) = y(i+1,1,j) - 2.0_dp*y(i,1,j) + y(i-1,1,j)
       end do
       ! Linear interpolate between first two points to avoid problems with
       ! square root divergence at origin
       y(1,2,j) = 0.0_dp ; y(2,2,j) = 0.0_dp ; y(n,2,j) = 0.0_dp
    end do

    ! Solve set of linear equations
    call dptsv(n-3,m,diag_work,subdiag_work,y(3,2,1),2*n,ierr)
    if (ierr /= 0) then
       write(stdout,'(a,i6)') &
            'Error in sw_spline_fit: dptsv failed with code ',ierr
       call comms_abort
    end if

    ! Deallocate temporary workspace
    deallocate(subdiag_work,stat=ierr)
    call utils_dealloc_check('sw_spline_fit','subdiag_work',ierr) !qoh

    deallocate(diag_work,stat=ierr)
    call utils_dealloc_check('sw_spline_fit','diag_work',ierr) !qoh

  end subroutine sw_spline_fit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sw_real_sph_harm(xpt, ypt, zpt, rad, lval, em)

    !=====================================================================!
    ! This function calculates the real spherical harmonic at a point.    !
    ! Originally based on code from function ngwf_at_point written by     !
    ! Chris-Kriton Skylaris.  Now rewritten and optimised for speed.      !
    !---------------------------------------------------------------------!
    ! Written by Quintin Hill in September 2009.                          !
    ! Fixed by Jacek Dziedzic in 2011.09 -- l=4,m=-3 was flipped with     !
    ! l=4,m=3 and l=4,m=-1 was flipped with l=4,m=1.                      !
    !=====================================================================!

    use comms, only: comms_abort
    use constants, only: stdout
    use geometry, only: point

    implicit none

    real(kind=DP) :: sw_real_sph_harm !Value of RSH at point
    real(kind=DP), intent(in) :: xpt, ypt, zpt  ! Point coordinates
    real(kind=DP), intent(in) :: rad ! distance of point from origin
    integer, intent(in) :: lval ! l
    integer, intent(in) :: em  ! m

    real(kind=DP) :: cos_theta, rsq

    sw_real_sph_harm = 0.0_DP

    select case (lval)
    case (0)
       ! qoh: Z_{00} function
       ! z_00 = 0.5_DP/sqrt(pi)
       sw_real_sph_harm=0.282094791773878_DP

    case(1)
       ! qoh: Z_{1m} function
       ! normz_1m = sqrt(0.75_DP)/sqrtpi
       if (rad.gt.(1.d-10) ) then
          select case (em)
          case (-1)
             sw_real_sph_harm=0.48860251190292_DP*ypt/rad
          case(0)
             sw_real_sph_harm=0.48860251190292_DP*zpt/rad
          case(1)
             sw_real_sph_harm=0.48860251190292_DP*xpt/rad
          end select
       else if (em == 0) then
          ! qoh: is this the correct behaviour?
          sw_real_sph_harm = 0.48860251190292_DP
       end if
    case(2)
       ! qoh: Z_{2m} function
       if (rad.gt.(1.d-10) ) then
          select case (em)
          case(-2)
             sw_real_sph_harm=1.0925484305920790_DP*xpt*ypt/(rad*rad)
          case(-1)
             sw_real_sph_harm=1.0925484305920790_DP*zpt*ypt/(rad*rad)
          case(0)
             sw_real_sph_harm=0.31539156525251999_DP*&
                  (3.0_DP*(zpt*zpt/(rad*rad))-1.0_DP)
          case(1)
             sw_real_sph_harm=1.0925484305920790_DP*xpt*zpt/(rad*rad)
          case(2)
             sw_real_sph_harm=0.54627421529603948_DP*&
                  (xpt*xpt-ypt*ypt)/(rad*rad)
          end select
       else if (em == 0) then
          ! qoh: is this the correct behaviour?
          sw_real_sph_harm = 0.63078313050503998_DP
       end if

    case(3)
       ! qoh: Z_{3m} function
       if (rad.gt.(1.d-10) ) then
          select case (em)
          case(-3)
             !norm =sqrt(35.0_DP/(pi*32.0_DP) )
             sw_real_sph_harm = 0.59004358992664352_DP*&
                  (3.0_DP*xpt*xpt - ypt*ypt)*ypt/(rad*rad*rad)
          case (-2)
             !norm =sqrt(105.0_DP/(pi*4.0_DP) )
             sw_real_sph_harm = 2.8906114426405538_DP*&
                  xpt*ypt*zpt/(rad*rad*rad)
          case (-1)
             !norm =sqrt(21.0_DP/(pi*32.0_DP) )
             sw_real_sph_harm = 0.45704579946446572_DP*&
                  (5.0_DP*zpt*zpt/(rad*rad) - 1.0_DP)*ypt/rad
          case(0)
             !norm =sqrt(7.0_DP/(pi*16.0_DP) )
             cos_theta= zpt/rad
             sw_real_sph_harm = 0.37317633259011540_DP*&
                  (5.0_DP*cos_theta*cos_theta-3.0_DP)*cos_theta
          case(1)
             !norm =sqrt(21.0_DP/(pi*32.0_DP) )
             sw_real_sph_harm = 0.45704579946446572_DP*&
                  (5.0_DP*zpt*zpt/(rad*rad) - 1.0_DP)*xpt/rad
          case(2)
             !norm =sqrt(105.0_DP/(pi*16.0_DP) )
             sw_real_sph_harm = 1.4453057213202769_DP*&
                  (xpt*xpt-ypt*ypt)*zpt/(rad*rad*rad)
          case(3)
             !norm =sqrt(35.0_DP/(pi*32.0_DP) )
             sw_real_sph_harm = 0.59004358992664352_DP*&
                  (xpt*xpt - 3.0_DP*ypt*ypt)*xpt/(rad*rad*rad)
          end select
       else if (em == 0) then
          ! qoh: is this the correct behaviour?
          sw_real_sph_harm = 0.74635266518023080_DP
       end if

    case(4)
       ! qoh: Z_{4m} function
       if (rad.gt.(1.d-10) ) then
          rsq = rad*rad
          select case (em)
          case(-4)
             ! norm = 0.75_DP*sqrt(35.0_DP/PI)
             sw_real_sph_harm = 2.5033429417967046_DP * &
                  xpt*ypt*(xpt*xpt - ypt*ypt)/(rsq*rsq)
          case(-3) ! jd: Fixed
             ! norm = 0.75_DP*sqrt(35.0_DP/(2.0_DP*PI))
             sw_real_sph_harm = 1.7701307697799304_DP * &
                  (3.0_DP *xpt*xpt - ypt*ypt) *ypt*zpt/(rsq*rsq)
          case(-2)
             ! norm = 0.75_DP*sqrt(5.0_DP/PI)
             sw_real_sph_harm =  0.94617469575756008_DP * &
                  xpt*ypt*(7.0_DP*zpt*zpt - rsq)/(rsq*rsq)
          case(-1) ! jd: Fixed
             ! norm = 0.75_DP*sqrt(5.0_DP/(2.0_DP*PI))
             sw_real_sph_harm = 0.66904654355728921_DP * &
                  ypt*zpt*(7.0_DP*zpt*zpt - 3.0_DP*rsq)/(rsq*rsq)
          case(0)
             ! norm = 0.1875_DP*sqrt(1.0_DP/PI)
             sw_real_sph_harm = 0.10578554691520431_DP * &
                  (zpt*zpt*(35.0_DP*zpt*zpt-30.0_DP*rsq)/(rsq*rsq) + 3.0_DP)
          case(1) ! jd: Fixed
             ! norm = 0.75_DP*sqrt(5.0_DP/(2.0_DP*PI))
             sw_real_sph_harm = 0.66904654355728921_DP * &
                  xpt*zpt*(7.0_DP*zpt*zpt - 3.0_DP*rsq)/(rsq*rsq)
          case(2)
             ! norm = 0.375_DP*sqrt(5.0_DP/PI)
             sw_real_sph_harm = 0.47308734787878004_DP * &
                  (xpt*xpt - ypt * ypt)*(7.0_DP*zpt*zpt - rsq)/(rsq*rsq)
          case(3) ! jd: Fixed
             ! norm = 0.75_DP*sqrt(35.0_DP/(2.0_DP*PI))
             sw_real_sph_harm = 1.7701307697799304_DP * &
                  (xpt*xpt - 3.0_DP*ypt*ypt) *xpt*zpt/(rsq*rsq)
          case(4)
             ! norm = 0.1875_DP*sqrt(35.0_DP/PI)
             sw_real_sph_harm = 0.62583573544917614_DP * &
                  ((xpt*xpt - 6.0_DP*ypt*ypt) *xpt*xpt + ypt**4)/(rsq*rsq)
          end select
       else if (em == 0) then
          ! qoh: is this the correct behaviour?
          sw_real_sph_harm = 0.84628437532163447_DP
       end if

    case default
       !qoh: Try to use slow function if possible
       if (sw_initialised) then
          sw_real_sph_harm = sw_eval_harmonic(lval,em,point(xpt,ypt,zpt))
       else
          write (stdout,'(a,i4,a)') 'Error in sw_real_sph_harm: &
               &Angular momentum l=',lval,' too high.'
          call comms_abort
       end if
    end select

  end function sw_real_sph_harm
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine sw_grad_real_sph_harm(grad, xpt, ypt, zpt, rad, lval, em)

    !=====================================================================!
    ! This function calculates the gradients in real space of a real      !
    ! spherical harmonic at a point.                                      !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine on 29/07/2010, building on the routine     !
    ! sw_real_sph_harm, which was written by Quintin Hill in September    !
    ! 2009.                                                               !
    ! This was all done with Mathematica, so it's possible I could have   !
    ! mis-transcribed a formula or two for the higher angular momenta.    !
    ! jd 2011.09: Fixed the bug in the l=4, m={-3,-1,1,3} case, propagated!
    ! from sw_real_sph_harm.                                              !
    !=====================================================================!

    use comms, only: comms_abort
    use constants, only: stdout

    implicit none

    real(kind=DP) :: grad(3)  !Value of RSH at point
    real(kind=DP), intent(in) :: xpt, ypt, zpt  ! Point coordinates
    real(kind=DP), intent(in) :: rad ! distance of point from origin
    integer, intent(in) :: lval ! l
    integer, intent(in) :: em  ! m

    real(kind=DP) :: rsq
    
    grad(:) = -9999999999999.9999_DP

    select case (lval)
    case (0)
       ! ndmh: dZ_{00}/d{x,y,z}
       ! z_00 = 0.5_DP/sqrt(pi)
       !sw_real_sph_harm=0.282094791773878_DP
       grad(:) = 0.0_DP

    case(1)
       ! ndmh: dZ_{1m}/d{x,y,z}
       ! normz_1m = sqrt(0.75_DP)/sqrtpi
       if (rad.gt.(1.d-10) ) then
          select case (em)
          case (-1)
             !sw_real_sph_harm=0.48860251190292_DP*ypt/rad
             grad(1) = 0.48860251190292_DP*(-xpt*ypt/rad**3)
             grad(2) = 0.48860251190292_DP*(rad**2 - ypt**2)/rad**3
             grad(3) = 0.48860251190292_DP*(-ypt*zpt/rad**3)
          case(0)
             !sw_real_sph_harm=0.48860251190292_DP*zpt/rad
             grad(1) = 0.48860251190292_DP*(-xpt*zpt/rad**3)
             grad(2) = 0.48860251190292_DP*(-ypt*zpt/rad**3)
             grad(3) = 0.48860251190292_DP*(rad**2 - zpt**2)/rad**3
          case(1)
             !sw_real_sph_harm=0.48860251190292_DP*xpt/rad
             grad(1) = 0.48860251190292_DP*(rad**2 - xpt**2)/rad**3
             grad(2) = 0.48860251190292_DP*(-xpt*ypt/rad**3)
             grad(3) = 0.48860251190292_DP*(-xpt*zpt/rad**3)
          end select
       else
          ! ndmh: NB NB NB: In the limit r->0 for L=1, one component of 
          ! ndmh: the gradient wrt {x,y,z} diverges for each M value.
          ! ndmh: This routine therefore returns sqrt(3/4pi) for the
          ! ndmh: relevant component, and one must take 
          ! ndmh: lim(f(r)/r) * sqrt(3/4pi) to get the limit of the product
          ! ndmh: d (f(r) S_LM(rhat)) / d x_alpha as r->0.
          select case (em)
          case (-1)
             grad(1) = 0.0_DP
             grad(2) = 0.48860251190292_DP
             grad(3) = 0.0_DP
          case(0)
             grad(1) = 0.0_DP
             grad(2) = 0.0_DP
             grad(3) = 0.48860251190292_DP
          case(1)
             grad(1) = 0.48860251190292_DP
             grad(2) = 0.0_DP
             grad(3) = 0.0_DP
          end select
       end if
    case(2)
       ! ndmh: dZ_{2m}/d{x,y,z}
       if (rad.gt.(1.d-10) ) then
          select case (em)
          case(-2)
             !sw_real_sph_harm=1.0925484305920790_DP*xpt*ypt/(rad*rad)
             grad(1) = 1.0925484305920790_DP*((rad**2 - 2.0_DP*xpt**2)*ypt)/rad**4
             grad(2) = 1.0925484305920790_DP*(xpt*(rad**2 - 2.0_DP*ypt**2))/rad**4
             grad(3) = 1.0925484305920790_DP*(-2.0_DP*xpt*ypt*zpt)/rad**4
          case(-1)
             !sw_real_sph_harm=1.0925484305920790_DP*zpt*ypt/(rad*rad)
             grad(1) = 1.0925484305920790_DP*(-2.0_DP*xpt*ypt*zpt)/rad**4
             grad(2) = 1.0925484305920790_DP*((rad**2 - 2.0_DP*ypt**2)*zpt)/rad**4
             grad(3) = 1.0925484305920790_DP*(ypt*(rad**2 - 2.0_DP*zpt**2))/rad**4
          case(0)
             !sw_real_sph_harm=0.31539156525251999_DP*&
             !     (3.0_DP*(zpt*zpt/(rad*rad))-1.0_DP)
             grad(1) = 0.31539156525251999_DP*(-6.0_DP*xpt*zpt**2)/rad**4
             grad(2) = 0.31539156525251999_DP*(-6.0_DP*ypt*zpt**2)/rad**4
             grad(3) = 0.31539156525251999_DP*(6.0_DP*rad**2*zpt - 6.0_DP*zpt**3)/rad**4
          case(1)
             !sw_real_sph_harm=1.0925484305920790_DP*xpt*zpt/(rad*rad)
             grad(1) = 1.0925484305920790_DP*((rad**2 - 2.0_DP*xpt**2)*zpt)/rad**4
             grad(2) = 1.0925484305920790_DP*(-2.0_DP*xpt*ypt*zpt)/rad**4
             grad(3) = 1.0925484305920790_DP*(xpt*(rad**2 - 2.0_DP*zpt**2))/rad**4
          case(2)
             !sw_real_sph_harm=0.54627421529603948_DP*&
             !     (xpt*xpt-ypt*ypt)/(rad*rad)
             grad(1) = 0.54627421529603948_DP*(2.0_DP*xpt*(rad**2 - xpt**2 + ypt**2))/rad**4
             grad(2) = 0.54627421529603948_DP*(-2.0_DP*ypt*(rad**2 + xpt**2 - ypt**2))/rad**4
             grad(3) = 0.54627421529603948_DP*(-2.0_DP*(xpt**2 - ypt**2)*zpt)/rad**4
          end select
       else
          ! ndmh: is this the correct behaviour?
          grad(:) = 0.0_DP
       end if

    case(3)
       ! ndmh: dZ_{3m}/d{x,y,z}
       if (rad.gt.(1.d-10) ) then
          select case (em)
          case(-3)
             !norm =sqrt(35.0_DP/(pi*32.0_DP) )
             !sw_real_sph_harm = 0.59004358992664352_DP*&
             !     (3.0_DP*xpt*xpt - ypt*ypt)*ypt/(rad*rad*rad)
             grad(1) = 0.59004358992664352_DP*(3.0_DP*xpt*ypt* &
                  (2.0_DP*rad**2-3*xpt**2+ypt**2))/rad**5
             grad(2) = 0.59004358992664352_DP*3.0_DP* &
                  (-3.0_DP*xpt**2*ypt**2+ypt**4+rad**2*(xpt**2-ypt**2))/rad**5
             grad(3) = 0.59004358992664352_DP*(3.0_DP*ypt*zpt* &
                  (-3.0_DP*xpt**2+ypt**2))/rad**5
          case (-2)
             !norm =sqrt(105.0_DP/(pi*4.0_DP) )
             !sw_real_sph_harm = 2.8906114426405538_DP*&
             !     xpt*ypt*zpt/(rad*rad*rad)
             grad(1) = 2.8906114426405538_DP * ypt*zpt* &
                  (rad**2-3.0_DP*xpt**2)/rad**5
             grad(2) = 2.8906114426405538_DP * xpt*zpt* &
                  (rad**2-3.0_DP*ypt**2)/rad**5
             grad(3) = 2.8906114426405538_DP * xpt*ypt* &
                  (rad**2-3.0_DP*zpt**2)/rad**5
          case (-1)
             !norm =sqrt(21.0_DP/(pi*32.0_DP) )
             !sw_real_sph_harm = 0.45704579946446572_DP*&
             !     (5.0_DP*zpt*zpt/(rad*rad) - 1.0_DP)*ypt/rad
             !sw_real_sph_harm = 0.45704579946446572_DP*&
             !     (4.0_DP*zpt*zpt - xpt*xpt - ypt*ypt)*ypt/rad**3
             grad(1) = 0.45704579946446572_DP * xpt*ypt* &
                  (rad**2-15.0_DP*zpt**2)/rad**5
             grad(2) = 0.45704579946446572_DP * (-rad**4 &
                  - 15.0_DP*ypt**2*zpt**2 + rad**2*(ypt**2 + 5.0_DP*zpt**2)) &
                  /rad**5
             grad(3) = 0.45704579946446572_DP * ypt*zpt* &
                  (11.0_DP*rad**2 - 15.0_DP*zpt**2)/rad**5
          case(0)
             !norm =sqrt(7.0_DP/(pi*16.0_DP) )
             !cos_theta= zpt/rad
             !sw_real_sph_harm = 0.37317633259011540_DP*&
             !     (5.0_DP*cos_theta*cos_theta-3.0_DP)*cos_theta
             !sw_real_sph_harm = 0.37317633259011540_DP*&
             !     (5.0_DP*zpt**3/rad**3-3.0_DP*zpt/rad)
             !sw_real_sph_harm = 0.37317633259011540_DP*&
             !     (2.0_DP*zpt*zpt-3.0_DP*(xpt*xpt+ypt*ypt))*zpt/rad**3
             grad(1) = 0.37317633259011540_DP * 3.0_DP*xpt*zpt* &
                  (rad**2-5.0_DP*zpt**2)/rad**5
             grad(2) = 0.37317633259011540_DP * 3.0_DP*ypt*zpt* &
                  (rad**2-5.0_DP*zpt**2)/rad**5
             grad(3) = 0.37317633259011540_DP * (-3.0_DP*(rad**4 -  &
                  6.0_DP*zpt**2*rad**2 + 5.0_DP*zpt**4))/rad**5
          case(1)
             !norm =sqrt(21.0_DP/(pi*32.0_DP) )
             !sw_real_sph_harm = 0.45704579946446572_DP*&
             !     (5.0_DP*zpt*zpt/(rad*rad) - 1.0_DP)*xpt/rad
             !sw_real_sph_harm = 0.45704579946446572_DP*&
             !     (4.0_DP*zpt*zpt - xpt*xpt - ypt*ypt)*xpt/rad**3
             grad(1) = 0.45704579946446572_DP * (-rad**4 -  &
                  15.0_DP*xpt**2*zpt**2 + rad**2*(xpt**2+5.0_DP*zpt**2))/rad**5
             grad(2) = 0.45704579946446572_DP *xpt*ypt* &
                  (rad**2-15.0_DP*zpt**2)/rad**5
             grad(3) = 0.45704579946446572_DP *xpt*zpt* &
                  (11.0_DP*rad**2-15.0_DP*zpt**2)/rad**5
          case(2)
             !norm =sqrt(105.0_DP/(pi*16.0_DP) )
             !sw_real_sph_harm = 1.4453057213202769_DP*&
             !     (xpt*xpt-ypt*ypt)*zpt/(rad*rad*rad)
             grad(1) = 1.4453057213202769_DP * xpt * (2.0_DP*rad**2 &
                  - 3.0_DP*xpt**2 + 3.0_DP*ypt**2)*zpt / rad**5
             grad(2) = 1.4453057213202769_DP * ypt * (-2.0_DP*rad**2 &
                  - 3.0_DP*xpt**2 + 3.0_DP*ypt**2)*zpt / rad**5
             grad(3) = 1.4453057213202769_DP * (xpt**2-ypt**2)*(rad**2 &
                  - 3.0_DP*zpt**2) / rad**5
          case(3)
             !norm =sqrt(35.0_DP/(pi*32.0_DP) )
             !sw_real_sph_harm = 0.59004358992664352_DP*&
             !     (xpt*xpt - 3.0_DP*ypt*ypt)*xpt/(rad*rad*rad)
             grad(1) = 0.59004358992664352_DP * (-3.0_DP*xpt**4 &
                  + 9.0_DP * xpt**2*ypt**2 + 3.0_DP*rad**2*(xpt**2-ypt**2)) &
                  / rad**5
             grad(2) = 0.59004358992664352_DP * (-3.0_DP*xpt*ypt * &
                  (2.0_DP*rad**2 + xpt**2 - 3.0_DP*ypt**2)) / rad**5
             grad(3) = 0.59004358992664352_DP * (-3.0_DP*xpt*zpt* &
                  (xpt**2-3.0_DP*ypt**2)) / rad**5
          end select
       else
          ! ndmh: is this the correct behaviour?
          grad(:) = 0.0_DP
       end if

    case(4)
       ! qoh: Z_{4m} function
       if (rad.gt.(1.d-10) ) then
          rsq = rad*rad
          select case (em)
          case(-4)
             ! norm = 0.75_DP*sqrt(35.0_DP/PI)
             !sw_real_sph_harm = 2.5033429417967046_DP * &
             !     xpt*ypt*(xpt*xpt - ypt*ypt)/(rsq*rsq)
             grad(1) = 2.5033429417967046_DP*(ypt*(-4.0_DP*xpt**4 + 4.0_DP*xpt**2*ypt**2 + rad**2*(3.0_DP*xpt**2 - ypt**2)))/rad**6
             grad(2) = 2.5033429417967046_DP*(xpt*(-4.0_DP*xpt**2*ypt**2 + 4.0_DP*ypt**4 + rad**2*(xpt**2 - 3.0_DP*ypt**2)))/rad**6
             grad(3) = 2.5033429417967046_DP*(-4.0_DP*xpt*ypt*(xpt**2 - ypt**2)*zpt)/rad**6
          case(-3)
             ! jd: Fixed
             ! norm = 0.75_DP*sqrt(35.0_DP/(2.0_DP*PI))
             !sw_real_sph_harm = 1.7701307697799304_DP * &
             !     (3.0_DP *xpt*xpt - ypt*ypt) *ypt*zpt/(rsq*rsq)
             grad(1) = 1.7701307697799304_DP*(xpt*ypt*(6.0_DP*rad**2 - &
                  12.0_DP*xpt**2 + 4.0_DP*ypt**2)*zpt)/rad**6
             grad(2) = 1.7701307697799304_DP*((-12.0_DP*xpt**2*ypt**2 + &
                  4.0_DP*ypt**4 + rad**2*(3.0_DP*xpt**2 - 3.0_DP*ypt**2))*zpt)/rad**6
             grad(3) = 1.7701307697799304_DP*(ypt*(3.0_DP*xpt**2 - &
                  1.0_DP*ypt**2)*(rad**2 - 4.0_DP*zpt**2))/rad**6
          case(-2)
             ! norm = 0.75_DP*sqrt(5.0_DP/PI)
             !sw_real_sph_harm =  0.94617469575756008_DP * &
             !     xpt*ypt*(7.0_DP*zpt*zpt - rsq)/(rsq*rsq)
             grad(1) = 0.94617469575756008_DP*(ypt*(xpt**2*(4.0_DP*xpt**2 + &
                  4.0_DP*ypt**2 - 24.0_DP*zpt**2) + rad**2*(-3.0_DP*xpt**2 - 1.0_DP*ypt**2 + 6.0_DP*zpt**2)))/rad**6
             grad(2) = 0.94617469575756008_DP*(xpt*(ypt**2*(4.0_DP*xpt**2 + &
                  4.0_DP*ypt**2 - 24.0_DP*zpt**2) + rad**2*(-1.0_DP*xpt**2 - 3.0_DP*ypt**2 + 6.0_DP*zpt**2)))/rad**6
             grad(3) = 0.94617469575756008_DP*(xpt*ypt*zpt*(12.0_DP*rad**2 + &
                  4.0_DP*xpt**2 + 4.0_DP*ypt**2 - 24.0_DP*zpt**2))/rad**6
          case(-1)
             ! jd: Fixed
             ! norm = 0.75_DP*sqrt(5.0_DP/(2.0_DP*PI))
             !sw_real_sph_harm = 0.66904654355728921_DP * &
             !     ypt*zpt*(7.0_DP*zpt*zpt - 3.0_DP*rsq)/(rsq*rsq)
             grad(1) = 0.66904654355728921_DP*(xpt*ypt*zpt*(6.0_DP*rad**2 - 28.0_DP*zpt**2))/rad**6
             grad(2) = 0.66904654355728921_DP*(-3.0_DP*rad**4*zpt - &
                  28.0_DP*ypt**2*zpt**3 + rad**2*(6.0_DP*ypt**2*zpt + 7.0_DP*zpt**3))/rad**6
             grad(3) = 0.66904654355728921_DP*(ypt*(-3.0_DP*rad**4 + &
                  27.0_DP*rad**2*zpt**2 - 28.0_DP*zpt**4))/rad**6
          case(0)
             ! norm = 0.1875_DP*sqrt(1.0_DP/PI)
             !sw_real_sph_harm = 0.10578554691520431_DP * &
             !     (zpt*zpt*(35.0_DP*zpt*zpt-30.0_DP*rsq)/(rsq*rsq) + 3.0_DP)
             grad(1) = 0.10578554691520431_DP*(60.0_DP*rad**2*xpt*zpt**2 - 140.0_DP*xpt*zpt**4)/rad**6
             grad(2) = 0.10578554691520431_DP*(60.0_DP*rad**2*ypt*zpt**2 - 140.0_DP*ypt*zpt**4)/rad**6
             grad(3) = 0.10578554691520431_DP*(-60.0_DP*rad**4*zpt + 200.0_DP*rad**2*zpt**3 - 140.0_DP*zpt**5)/rad**6
          case(1)
             ! jd: Fixed
             ! norm = 0.75_DP*sqrt(5.0_DP/(2.0_DP*PI))
             !sw_real_sph_harm = 0.66904654355728921_DP * &
             !     xpt*zpt*(7.0_DP*zpt*zpt - 3.0_DP*rsq)/(rsq*rsq)
             grad(1) = 0.66904654355728921_DP*(-3*rad**4*zpt - 28*xpt**2*zpt**3 + rad**2*(6*xpt**2*zpt + 7*zpt**3))/rad**6
             grad(2) = 0.66904654355728921_DP*(2*xpt*ypt*zpt*(3*rad**2 - 14*zpt**2))/rad**6
             grad(3) = 0.66904654355728921_DP*(xpt*(-3*rad**4 + 27*rad**2*zpt**2 - 28*zpt**4))/rad**6
          case(2)
             ! norm = 0.375_DP*sqrt(5.0_DP/PI)
             !sw_real_sph_harm = 0.47308734787878004_DP * &
             !     (xpt*xpt - ypt * ypt)*(7.0_DP*zpt*zpt - rsq)/(rsq*rsq)
             grad(1) = 0.47308734787878004_DP*(xpt*(4.0_DP*xpt**4 - &
                  4.0_DP*ypt**4 - 24.0_DP*xpt**2*zpt**2 + 24.0_DP*ypt**2*zpt**2 + rad**2*(-4.0_DP*xpt**2 + 12.0_DP*zpt**2)))/rad**6
             grad(2) = 0.47308734787878004_DP*(ypt*(4.0_DP*xpt**4 - &
                  4.0_DP*ypt**4 - 24.0_DP*xpt**2*zpt**2 + 24.0_DP*ypt**2*zpt**2 + rad**2*(4.0_DP*ypt**2 - 12.0_DP*zpt**2)))/rad**6
             grad(3) = 0.47308734787878004_DP*(zpt*(4.0_DP*xpt**4 - &
                  4.0_DP*ypt**4 + rad**2*(12.0_DP*xpt**2 - 12.0_DP*ypt**2) - 24.0_DP*xpt**2*zpt**2 + 24.0_DP*ypt**2*zpt**2))/rad**6
          case(3)
             ! jd: Fixed
             ! norm = 0.75_DP*sqrt(35.0_DP/(2.0_DP*PI))
             !sw_real_sph_harm = 1.7701307697799304_DP * &
             !     (xpt*xpt - 3.0_DP*ypt*ypt) *xpt*zpt/(rsq*rsq)
             grad(1) = 1.7701307697799304_DP*((-4.0_DP*xpt**4 + &
                  12.0_DP*xpt**2*ypt**2 + rad**2*(3.0_DP*xpt**2 - 3.0_DP*ypt**2))*zpt)/rad**6
             grad(2) = 1.7701307697799304_DP*(xpt*ypt*(-6.0_DP*rad**2 - 4.0_DP*xpt**2 + 12.0_DP*ypt**2)*zpt)/rad**6
             grad(3) = 1.7701307697799304_DP*(xpt*(xpt**2 - 3.0_DP*ypt**2)*(rad**2 - 4.0_DP*zpt**2))/rad**6
          case(4)
             ! norm = 0.1875_DP*sqrt(35.0_DP/PI)
             !sw_real_sph_harm = 0.62583573544917614_DP * &
             !     ((xpt*xpt - 6.0_DP*ypt*ypt) *xpt*xpt + ypt**4)/(rsq*rsq)
             grad(1) = 0.62583573544917614_DP*(-4.0_DP*xpt**5 + &
                  24.0_DP*xpt**3*ypt**2 - 4.0_DP*xpt*ypt**4 + rad**2*(4.0_DP*xpt**3 - 12.0_DP*xpt*ypt**2))/rad**6
             grad(2) = 0.62583573544917614_DP*(-16.0_DP*xpt**4*ypt + &
                  4.0_DP*ypt**3*zpt**2 + xpt**2*(16.0_DP*ypt**3 - 12.0_DP*ypt*zpt**2))/rad**6
             grad(3) = 0.62583573544917614_DP*(-4.0_DP*(xpt**4 - &
                  6.0_DP*xpt**2*ypt**2 + ypt**4)*zpt)/rad**6
          end select
       else
          ! ndmh: is this the correct behaviour?
          grad(:) = 0.0_DP
       end if

    case default
       ! ndmh: abort if L>4
       write (stdout,'(a,i4,a)') 'Error in sw_grad_real_sph_harm: &
             &Angular momentum l=',lval,' too high.'
       call comms_abort
    end select

  end subroutine sw_grad_real_sph_harm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sw_real_sph_harm_unit(xpt, ypt, zpt, lval, em)

    !=====================================================================!
    ! This function calculates the real spherical harmonic at a point,    !
    ! where (x,y,z) is a unit vector. Based on sw_real_sph_harm.          !
    ! This routine avoids division and reduces branching.                 !
    !---------------------------------------------------------------------!
    ! Written by Quintin Hill on 21/04/2010.                              !
    !=====================================================================!

    use comms, only: comms_abort
    use constants, only: stdout
    use geometry, only: point

    implicit none

    real(kind=DP) :: sw_real_sph_harm_unit ! Value of RSH at point
    real(kind=DP), intent(in) :: xpt, ypt, zpt  ! Point coordinates
    integer, intent(in) :: lval ! l
    integer, intent(in) :: em  ! m

    sw_real_sph_harm_unit = 0.0_DP

    select case (lval)
    case (0)
       ! qoh: Z_{00} function
       ! z_00 = 0.5_DP/sqrt(pi)
       sw_real_sph_harm_unit=0.282094791773878_DP

    case(1)
       ! qoh: Z_{1m} function
       ! normz_1m = sqrt(0.75_DP)/sqrtpi
       select case (em)
       case (-1)
          sw_real_sph_harm_unit=0.48860251190292_DP*ypt
       case(0)
          sw_real_sph_harm_unit=0.48860251190292_DP*zpt
       case(1)
          sw_real_sph_harm_unit=0.48860251190292_DP*xpt
       end select

    case(2)
       ! qoh: Z_{2m} function
       select case (em)
       case(-2)
          sw_real_sph_harm_unit=1.0925484305920790_DP*xpt*ypt
       case(-1)
          sw_real_sph_harm_unit=1.0925484305920790_DP*zpt*ypt
       case(0)
          sw_real_sph_harm_unit=0.31539156525251999_DP*(3.0_DP*zpt*zpt - 1.0_DP)
       case(1)
          sw_real_sph_harm_unit=1.0925484305920790_DP*xpt*zpt
       case(2)
          sw_real_sph_harm_unit=0.54627421529603948_DP*(xpt*xpt-ypt*ypt)
       end select

    case(3)
       ! qoh: Z_{3m} function
       select case (em)
       case(-3)
          !norm =sqrt(35.0_DP/(pi*32.0_DP) )
          sw_real_sph_harm_unit = 0.59004358992664352_DP*&
               (3.0_DP*xpt*xpt - ypt*ypt)*ypt
       case (-2)
          !norm =sqrt(105.0_DP/(pi*4.0_DP) )
          sw_real_sph_harm_unit = 2.8906114426405538_DP*&
               xpt*ypt*zpt
       case (-1)
          !norm =sqrt(21.0_DP/(pi*32.0_DP) )
          sw_real_sph_harm_unit = 0.45704579946446572_DP*&
               (5.0_DP*zpt*zpt - 1.0_DP)*ypt
       case(0)
          !norm =sqrt(7.0_DP/(pi*16.0_DP) )
          sw_real_sph_harm_unit = 0.37317633259011540_DP*&
               (5.0_DP*zpt*zpt - 3.0_DP)*zpt
       case(1)
          !norm =sqrt(21.0_DP/(pi*32.0_DP) )
          sw_real_sph_harm_unit = 0.45704579946446572_DP*&
               (5.0_DP*zpt*zpt - 1.0_DP)*xpt
       case(2)
          !norm =sqrt(105.0_DP/(pi*16.0_DP) )
          sw_real_sph_harm_unit = 1.4453057213202769_DP*&
               (xpt*xpt-ypt*ypt)*zpt
       case(3)
          !norm =sqrt(35.0_DP/(pi*32.0_DP) )
          sw_real_sph_harm_unit = 0.59004358992664352_DP*&
               (xpt*xpt - 3.0_DP*ypt*ypt)*xpt
       end select

    case(4)
       ! qoh: Z_{4m} function
       select case (em)
       case(-4)
          ! norm = 0.75_DP*sqrt(35.0_DP/PI)
          sw_real_sph_harm_unit = 2.5033429417967046_DP * &
               xpt*ypt*(xpt*xpt - ypt*ypt)
       case(-3) ! jd: Fixed
          ! norm = 0.75_DP*sqrt(35.0_DP/(2.0_DP*PI))
          sw_real_sph_harm_unit = 1.7701307697799304_DP * &
               (3.0_DP*xpt*xpt - ypt*ypt) *ypt*zpt
       case(-2)
          ! norm = 0.75_DP*sqrt(5.0_DP/PI)
          sw_real_sph_harm_unit =  0.94617469575756008_DP * &
               xpt*ypt*(7.0_DP*zpt*zpt - 1.0_DP)
       case(-1) ! jd: Fixed
          ! norm = 0.75_DP*sqrt(5.0_DP/(2.0_DP*PI))
          sw_real_sph_harm_unit = 0.66904654355728921_DP * &
               ypt*zpt*(7.0_DP*zpt*zpt - 3.0_DP)
       case(0)
          ! norm = 0.1875_DP*sqrt(1.0_DP/PI)
          sw_real_sph_harm_unit = 0.10578554691520431_DP * &
               (zpt*zpt*(35.0_DP*zpt*zpt-30.0_DP) + 3.0_DP)
       case(1) ! jd: Fixed
          ! norm = 0.75_DP*sqrt(5.0_DP/(2.0_DP*PI))
          sw_real_sph_harm_unit = 0.66904654355728921_DP * &
               xpt*zpt*(7.0_DP*zpt*zpt - 3.0_DP)
       case(2)
          ! norm = 0.375_DP*sqrt(5.0_DP/PI)
          sw_real_sph_harm_unit = 0.47308734787878004_DP * &
               (xpt*xpt - ypt * ypt)*(7.0_DP*zpt*zpt - 1.0_DP)
       case(3) ! jd: Fixed
          ! norm = 0.75_DP*sqrt(35.0_DP/(2.0_DP*PI))
          sw_real_sph_harm_unit = 1.7701307697799304_DP * &
               (xpt*xpt - 3.0_DP*ypt*ypt) *xpt*zpt
       case(4)
          ! norm = 0.1875_DP*sqrt(35.0_DP/PI)
          sw_real_sph_harm_unit = 0.62583573544917614_DP * &
               ((xpt*xpt - 6.0_DP*ypt*ypt) *xpt*xpt + ypt**4)
       end select

    case default
       !qoh: Try to use slow function if possible
       if (sw_initialised) then
          sw_real_sph_harm_unit = sw_eval_harmonic(lval,em,point(xpt,ypt,zpt))
       else
          write (stdout,'(2a)')&
               'Angular momentum too high in sw_real_sph_harm_unit',&
               'Program execution stops.'
          call comms_abort
       end if
    end select

  end function sw_real_sph_harm_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sw_exit

    !==================================================================!
    ! This subroutine deallocates the module arrays.                   !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in July 2008                             !
    !==================================================================!

    use utils, only: utils_dealloc_check

    implicit none
    integer :: ierr

    if (sw_initialised) then
       deallocate(sw_bessel_zeros,stat=ierr)
       call utils_dealloc_check('sw_exit','sw_bessel_zeros',ierr)
       deallocate(bessel_spline,stat=ierr)
       call utils_dealloc_check('sw_exit','bessel_spline',ierr)
       deallocate(legendre_spline,stat=ierr)
       call utils_dealloc_check('sw_exit','legendre_spline',ierr)
       deallocate(trig_spline,stat=ierr)
       call utils_dealloc_check('sw_exit','trig_spline',ierr)
       deallocate(harmonic_norm,stat=ierr)
       call utils_dealloc_check('sw_exit','harmonic_norm',ierr)
       sw_initialised = .false.
    else if (allocated(sw_bessel_zeros)) then
       deallocate(sw_bessel_zeros,stat=ierr)
       call utils_dealloc_check('sw_exit','sw_bessel_zeros',ierr)
    end if

  end subroutine sw_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sw_recp_generate(sw_out,sw_work,lval,mval,qnl,radius,disp)

    !====================================================================!
    ! This subroutine generates spherical waves in reciprocal space and  !
    ! performs an inverse Fourier transform in an FFT box to put them in !
    ! real space                                                         !
    !--------------------------------------------------------------------!
    ! Originally written by Mark Robinson, November 2007, as an internal !
    ! subroutine of properties_ngwfs_char                                !
    ! Moved by Alvaro Ruiz Serrano to spherical waves module in          !
    ! February 2009.                                                     !
    !====================================================================!

    use constants, only: cmplx_i, PI
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, operator(.DOT.)
    use simulation_cell, only: pub_fftbox
    use timer, only: timer_clock

    ! arguments
    real(kind=DP), intent(out)      :: sw_out(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    complex(kind=DP), intent(inout) :: sw_work(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)
    integer,       intent(in)       :: lval,mval
    real(kind=DP), intent(in)       :: qnl,radius
    type(POINT),   intent(in)       :: disp

    ! variables
    integer          :: nx,ny,nz
    type(POINT)      :: gvec
    real(kind=DP)    :: rad2,qnl2
    real(kind=DP)    :: gmag, gmagsq ! G vector magnitude and its square
    real(kind=DP)    :: sph_bess     ! Spherical bessel function at point
    real(kind=DP)    :: sph_harm     ! Spherical harmonic at point
    complex(kind=DP) :: fac1,fac2    ! Common factors

    call timer_clock('sw_recp_generate',1)

    ! calculate factor common to all points
    rad2 = radius**2
    qnl2 = qnl**2
    sph_bess = sw_bessel_fast(lval-1,rad2*qnl2)
    fac1 = 4*PI * cmplx_i**lval * qnl*rad2 * sph_bess
    fac2 = fac1 * radius * sph_bess

    ! calculate a SW in the reciprocal FFT box
    do nz=1,pub_fftbox%total_pt3
       do ny=1,pub_fftbox%total_ld2
          do nx=1,pub_fftbox%total_ld1
             ! create gvec
             gvec%X = pub_fftbox%recip_grid(1,nx,ny,nz)
             gvec%Y = pub_fftbox%recip_grid(2,nx,ny,nz)
             gvec%Z = pub_fftbox%recip_grid(3,nx,ny,nz)

             gmag   = pub_fftbox%recip_grid(4,nx,ny,nz)
             gmagsq = 2 * pub_fftbox%recip_grid(5,nx,ny,nz)

             ! calculate spherical harmonic
             sph_harm = sw_real_sph_harm(gvec%x,gvec%y,gvec%z,gmag,&
                 lval,mval)

             if(abs(gmagsq-qnl2) < epsilon(1.0_DP)) then
                ! calculate spherical wave
                sw_work(nx,ny,nz) = fac2 / (gmag + qnl) * sph_harm
             else
                ! calculate spherical bessel function
                sph_bess = sw_bessel_fast(lval,gmagsq*rad2)
                ! calculate spherical wave
                sw_work(nx,ny,nz) = fac1 / (gmagsq - qnl2) * sph_harm * sph_bess
             endif

             ! shift atomic center
             sw_work(nx,ny,nz) = sw_work(nx,ny,nz)* &
                  exp(-cmplx_i*(gvec.dot.disp))
          enddo
       enddo
    enddo

    ! FFT to real space and place into output array
    call fourier_apply_box('Coarse','Backwards',sw_work)

    do nz=1,pub_fftbox%total_pt3
       do ny=1,pub_fftbox%total_ld2
          do nx=1,pub_fftbox%total_ld1
             sw_out(nx,ny,nz) = real(sw_work(nx,ny,nz),kind=DP)
          end do
       end do
    end do

    call timer_clock('sw_recp_generate',2)

  end subroutine sw_recp_generate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sw_recp_generate_in_tb(sw_out,sw_work,lval,mval,qnl,aval,disp, &
       tb_n1,tb_n2,tb_n3)

    !====================================================================!
    ! This subroutine generates spherical waves in reciprocal space and  !
    ! performs an inverse Fourier transform in a tightbox to put them in !
    ! real space.                                                        !
    !--------------------------------------------------------------------!
    ! Originally written by Mark Robinson, November 2007, as an internal !
    ! subroutine of properties_ngwfs_char                                !
    ! Moved by Alvaro Ruiz Serrano to spherical waves module in          !
    ! February 2009.                                                     !
    ! Converted to use tightboxes by Quintin Hill and Alvaro Ruiz        !
    ! Serrano in March 2009.                                             !
    !====================================================================!

    use constants, only: cmplx_i, PI
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, operator(.DOT.)
    use simulation_cell, only: pub_tb_recip_grid
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer, intent(in) :: tb_n1,tb_n2,tb_n3
    real(kind=DP), intent(out)    :: sw_out(tb_n1,tb_n2,tb_n3)
    complex(kind=DP), intent(out) :: sw_work(tb_n1,tb_n2,tb_n3)
    integer,       intent(in) :: lval
    integer,       intent(in) :: mval
    real(kind=DP), intent(in) :: qnl  ! j_l(qr)
    real(kind=DP), intent(in) :: aval ! Cutoff radius of spherical wave
    type(POINT),   intent(in) :: disp ! Tightbox origin relative to origin of sw

    ! Internal variables
    integer          :: tbgp1, tbgp2, tbgp3 ! Tightbox grid points
    type(POINT)      :: gvec         !G vector
    real(kind=DP)    :: asq          ! a squared (square of cutoff radius)
    real(kind=DP)    :: qsq          ! qnl squared
    real(kind=DP)    :: gmag, gmagsq ! G vector magnitude and its square
    real(kind=DP)    :: sph_bess     ! Spherical Bessel function at point
    real(kind=DP)    :: sph_harm     ! Spherical harmonic at point
    complex(kind=DP) :: fac1,fac2    ! Common factors

    call timer_clock('sw_recp_generate_in_tb',1)

    ! calculate factor common to all points
    qsq = qnl*qnl
    asq = aval*aval
    sph_bess = sw_bessel_fast(lval-1,qsq*asq)
    fac1 = 4.0_DP*PI * qnl *asq * sph_bess * cmplx_i**lval
    fac2 = fac1 * aval * sph_bess

    ! calculate a SW in the reciprocal FFT box
    do tbgp3=1,tb_n3
       do tbgp2=1,tb_n2
          do tbgp1=1,tb_n1
             ! create gvec
             gvec%X = pub_tb_recip_grid(1,tbgp1,tbgp2,tbgp3)
             gvec%Y = pub_tb_recip_grid(2,tbgp1,tbgp2,tbgp3)
             gvec%Z = pub_tb_recip_grid(3,tbgp1,tbgp2,tbgp3)

             ! calculate gmag
             gmag =  pub_tb_recip_grid(4,tbgp1,tbgp2,tbgp3)
             gmagsq =  2.0_DP * pub_tb_recip_grid(5,tbgp1,tbgp2,tbgp3)

             ! calculate spherical harmonic
             sph_harm = sw_real_sph_harm(gvec%x,gvec%y,gvec%z,gmag,&
                 lval,mval)
             if(abs(gmagsq-qsq) < epsilon(1.0_DP)) then
                ! calculate spherical wave
                sw_work(tbgp1,tbgp2,tbgp3) = fac2 / (gmag + qnl) * sph_harm
             else
                ! calculate spherical bessel function
                sph_bess = sw_bessel_fast(lval,gmagsq*asq)
                ! calculate spherical wave
                sw_work(tbgp1,tbgp2,tbgp3) = sph_harm * sph_bess * fac1 / &
                     (gmagsq - qsq)
             endif

             ! shift atomic center
             sw_work(tbgp1,tbgp2,tbgp3) = sw_work(tbgp1,tbgp2,tbgp3)&
                  *exp(-cmplx_i*(gvec.DOT.disp))
          enddo
       enddo
    enddo

    call fourier_apply_box('Coarse','Backwards',sw_work,tb=.true.)

    sw_out = real(sw_work,kind=DP)

    call timer_clock('sw_recp_generate_in_tb',2)

  end subroutine sw_recp_generate_in_tb

end module spherical_wave

