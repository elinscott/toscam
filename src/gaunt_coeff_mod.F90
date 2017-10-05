! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                    Gaunt coefficients module                   !
!                                                                !
! This module calculates Gaunt coefficients and Real Gaunt       !
! coefficients in an indexed form for efficient storage.         !
!----------------------------------------------------------------!
! Written by Nicholas Hine in 2010.                              !
! Cleaned up by Nicholas Hine in April 2011.                     !
!================================================================!

module gaunt_coeff

  use constants, only: DP

  implicit none

  private

  ! Type to store evaluated Gaunt coefficients efficiently
  type GAUNT_COEFFS
     integer :: lmax
     integer :: ngaunt
     integer :: ngradcoeff
     integer,allocatable :: jsel(:,:)
     integer,allocatable :: ksel(:,:)
     integer,allocatable :: gsel(:,:,:)
     real(kind=DP), allocatable :: coeff(:)
     real(kind=DP), allocatable :: grad_coeff(:,:)
  end type GAUNT_COEFFS

  public :: GAUNT_COEFFS

  type(GAUNT_COEFFS) :: gaunt

  public :: gaunt_init
  public :: gaunt_exit
  public :: realgaunt
  public :: gaunt_grad_integrals

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function realgaunt(lup,mup,li,mi,lj,mj)

    implicit none

    ! Arguments
    integer, intent(in) :: lup,mup
    integer, intent(in) :: li,mi
    integer, intent(in) :: lj,mj

    ! Local Variables
    integer :: isel, jsel, ksel

    isel = lup*lup + lup + 1 + mup
    jsel = gaunt%jsel(li*li+li+1+mi,lj*lj+lj+1+mj)
    ksel = gaunt%ksel(isel,jsel)
    realgaunt = gaunt%coeff(ksel)

  end function realgaunt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine gaunt_grad_integrals(angint,li,mi,lj,mj)

    use constants, only: SQRT_PI

    implicit none

    ! Arguments
    real(kind=DP) :: angint(3,3)
    integer, intent(in) :: li,mi
    integer, intent(in) :: lj,mj

    ! Local Variables
    integer :: i,j

    i = li*li + li + mi + 1
    j = lj*lj + lj + mj + 1

    angint(1,1) = realgaunt(1,1,li,mi,lj,mj) * SQRT_PI * 2.0_DP / sqrt(3.0_DP)
    angint(2,1) = realgaunt(1,-1,li,mi,lj,mj) * SQRT_PI * 2.0_DP / sqrt(3.0_DP)
    angint(3,1) = realgaunt(1,0,li,mi,lj,mj) * SQRT_PI * 2.0_DP / sqrt(3.0_DP)
    angint(1,2) = internal_get_grad_int(4,i,j)
    angint(2,2) = internal_get_grad_int(5,i,j)
    angint(3,2) = internal_get_grad_int(6,i,j)
    angint(1,3) = internal_get_grad_int(7,i,j)
    angint(2,3) = internal_get_grad_int(8,i,j)
    angint(3,3) = 0.0_DP

contains

    real(kind=DP) function internal_get_grad_int(w,i,j)

       integer :: w,i,j
       
       internal_get_grad_int = 0.0_DP
    
    end function internal_get_grad_int

  end subroutine gaunt_grad_integrals


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine gaunt_init(lmax)

    !==================================================================!
    ! This subroutine initialises arrays for storing the real Gaunt    !
    ! coefficients up to a maximum angular momentum lmax.              !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 26/05/10.           !
    !==================================================================!

    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: lmax

    ! Local Variables
    integer :: ierr
    integer :: isel,jsel,k
    real(kind=DP),allocatable :: tmp(:)

    ! Set highest angular momentum over all species
    gaunt%lmax = lmax

    ! Allocate storage for selection rules
    allocate(gaunt%ksel((2*gaunt%lmax-1)**2, &
         gaunt%lmax**2*(gaunt%lmax**2+1)/2),stat=ierr)
    call utils_alloc_check('gaunt_init','gaunt%ksel',ierr)

    allocate(gaunt%jsel((gaunt%lmax)**2,(gaunt%lmax)**2),stat=ierr)
    call utils_alloc_check('gaunt_init','gaunt%jsel',ierr)

    ! Make list of indices to li,mi components
    k = 0
    do isel=1,(gaunt%lmax)**2
       do jsel=1,isel
          k = k + 1
          gaunt%jsel(isel,jsel) = k
          gaunt%jsel(jsel,isel) = k
       end do
    end do

    ! Allocate temporary array for Gaunt coefficients
    allocate(tmp((2*gaunt%lmax-1)**2*gaunt%lmax**4),stat=ierr)
    call utils_alloc_check('gaunt_init','tmp',ierr)

    ! Calculate the coefficients
    call gaunt_calculate(gaunt%lmax,gaunt%ngaunt,gaunt%ksel,tmp)

    ! Allocate array for storage of Gaunt coefficients
    allocate(gaunt%coeff(0:gaunt%ngaunt),stat=ierr)
    call utils_alloc_check('gaunt_init','gaunt%coeff',ierr)

    ! Copy coefficients from temporary array
    gaunt%coeff(0) = 0.0_DP
    gaunt%coeff(1:gaunt%ngaunt) = tmp(1:gaunt%ngaunt)
    ! Deallocate temporary storage
    deallocate(tmp,stat=ierr)
    call utils_dealloc_check('gaunt_init','tmp',ierr)

  end subroutine gaunt_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine gaunt_exit

    !==================================================================!
    ! This subroutine deallocates arrays for storing the real Gaunt    !
    ! coefficients.                                                    !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 26/05/10.           !
    !==================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Local Variables
    integer :: ierr

    ! Deallocate Gaunt coefficient arrays
    deallocate(gaunt%coeff,stat=ierr)
    call utils_dealloc_check('gaunt_exit','gaunt%coeff',ierr)
    deallocate(gaunt%jsel,stat=ierr)
    call utils_dealloc_check('gaunt_exit','gaunt%jsel',ierr)
    deallocate(gaunt%ksel,stat=ierr)
    call utils_dealloc_check('gaunt_exit','gaunt%ksel',ierr)

  end subroutine gaunt_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function fact(n)

    !==================================================================!
    ! This function calculates the factorial of an integer n.          !
    !------------------------------------------------------------------!
    ! This function was written by Nicholas Hine on 26/05/10.          !
    !==================================================================!

    implicit none

    ! Arguments
    integer,intent(in) :: n

    ! Local Variables
    integer :: i

    fact = 1.0_DP
    do i=2,n
       fact = fact*real(i,kind=DP)
    end do

  end function fact


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function perm(n,k)

    !==================================================================!
    ! This function calculates the number of permutations of k objects !
    ! in a set of n.                                                   !
    !------------------------------------------------------------------!
    ! This function was written by Nicholas Hine on 26/05/10.          !
    !==================================================================!

    implicit none

    ! Arguments
    integer,intent(in) :: n,k

    ! Local Variables
    integer :: i

    if ((n>=0).and.(n>=k)) then
       perm = 1.0_DP
       do i=n-k+1,n
          perm = perm*real(i,kind=DP)
       end do
    else
       perm = 0.0_DP
    end if

  end function perm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function rgaunt(ll,mm,l1,m1,l2,m2)

    !==================================================================!
    ! This function calculates the Gaunt coefficient for a given set   !
    ! of values of l1,m1,ll,mm,l2,m2, where G(l1,m1,ll,mm,l2,m2) =     !
    ! \int \sqrt{4\pi} Y*(l1,m1) Y*(ll,mm) Y(l2,m2) d\hat{r}           !
    !------------------------------------------------------------------!
    ! This function was written by Nicholas Hine on 26/05/10, with     !
    ! inspiration provided by the corresponding ABINIT routine.        !
    !==================================================================!

    implicit none

    ! Arguments
    integer,intent(in) :: ll,mm
    integer,intent(in) :: l1,m1
    integer,intent(in) :: l2,m2

    ! Local Variables
    integer :: i1,i2
    integer :: ls2
    integer :: lsl
    integer :: ls1
    integer :: lsum
    integer :: j1,j2
    integer :: n1,n2
    real(kind=DP) :: sgn,sum
    real(kind=DP) :: fac2,fac1
    logical :: nonzero

    ! Initialisations
    rgaunt = 0.0_DP
    nonzero = .true.

    ! Check if this coefficient is nonzero by applying selection rules
    if (abs(m1) > l1) nonzero = .false.
    if (abs(mm) > ll) nonzero = .false.
    if (abs(m2) > l2) nonzero = .false.
    if ((-m1 - mm + m2) /= 0) nonzero = .false.
    lsum = l1 + ll + l2
    if (mod(lsum,2) /= 0) nonzero = .false.
    ls2 = lsum - 2*l2
    lsl = lsum - 2*ll
    ls1 = lsum - 2*l1
    if ((ls2 < 0) .or. (lsl < 0) .or. (ls1 < 0)) nonzero = .false.

    ! Coefficient is nonzero, so evaluate it
    if (nonzero) then

       ! G = -1^(ls2/2)*sqrt((2l1+1)*(2ll+1)*(2l2+1))*lsl!*ls1!/(lsum!*(lsl/2)!*(ls1/2)!)
       fac1 = fact(l2 + m2) * fact(l2 - m2)
       fac2 = (2*l1 + 1)*(2*ll + 1)*(2*l2 + 1)
       if (mm >= 0) then
          fac1 = fac1*perm(ll+mm,2*mm)
       else
          fac1 = fac1/perm(ll-mm,-2*mm)
       end if
       if (m1 >= 0) then
          fac1 = fac1/perm(l1+m1,2*m1)
       else
          fac1 = fac1*perm(l1-m1,-2*m1)
       end if
       rgaunt = (-1)**(ls2/2) * sqrt(fac2) &
            * fact(lsl) * fact(ls1) / fact(lsum+1) &
            * fact(lsum/2) / (fact(ls2/2) &
            * fact(lsl/2) * fact(ls1/2)) * sqrt(fac1)

       ! Set up shorthand variables
       i1 = l2 - ll - m1;    i2 = l2 - l1 + mm
       n1 = l1 + m1;         n2 = ll - mm
       j1 = -min(0, i1, i2); j2 = min(ls2, n1, n2)

       if (j1>0) then
          sgn = (-1.0_DP)**j1
       else
          sgn = 1.0_DP
       end if

       ! Compute \sum_j1^j2 (-1)^j1 *
       ! P(n1,j1)*P(n2,j1)*P(l1+ll-l2,j1)/(j1!*(i1+j1)!*(i2+j1)!
       sum = 0.0_DP
       do while(j1 <= j2)
          sum = sum + sgn * perm(n1,j1)/fact(j1) &
               * perm(n2,j1)/fact(i1 + j1) &
               * perm(ls2,j1)/fact(i2 + j1)
          j1 = j1 + 1
          sgn = -sgn
       end do
       rgaunt = rgaunt * sum

    end if

  end function rgaunt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine gaunt_calculate(lmax,n_gaunt,gaunt_select,real_gaunt)

    !===================================================================!
    ! This function calculates a set of Real Gaunt coefficients up to   !
    ! a given lmax value, for all the possible combinations of l1,m1,   !
    ! ll,mm,l2 and m2 allowed, and indexes them in an efficient manner. !
    !-------------------------------------------------------------------!
    ! Arguments:                                                        !
    !  lmax     (input)      : maximum angular momentum value to use:   !
    !      args up to (2*lmax-1,m),(lmax,m),(lmax,m) will be calculated !
    !  n_gaunt (output)      : number of nonzero coefficients found     !
    !  gaunt_select (output) : selection rules indicating nonzero coeffs!
    !  real_gaunt (output)   : array of nonzero real Gaunt coefficients !
    !-------------------------------------------------------------------!
    ! This function was written by Nicholas Hine on 26/05/10, with      !
    ! inspiration provided by the corresponding ABINIT routine.         !
    !===================================================================!

    use constants, only: DP, PI
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: lmax
    integer,intent(out) :: n_gaunt
    integer,intent(out) :: gaunt_select((2*lmax-1)**2,lmax**2*(lmax**2+1)/2)
    real(kind=DP),intent(out) :: real_gaunt((2*lmax-1)**2*(lmax)**4)

    ! Local Variables
    integer :: idx1,idx2,idxup
    integer :: l1,l2,m1,m2,lup,mup
    integer :: m1p,mupp,m2p
    integer :: sgn1,sgnup
    real(kind=DP) :: ccprod
    real(kind=DP) :: real_gaunt_tmp
    integer :: ierr
    logical :: nonzero

    ! Type to store a set of arrays of different lengths
    type c3_arr
       complex(kind=DP), allocatable :: val(:,:)
    end type c3_arr
    type(c3_arr), allocatable :: cc(:)

    ! Pre-compute coeffs relating Real Spherical Harmonics to normal Spherical
    ! Harmonics: S_lm = cc(l,m',m)*Ylm'
    ! cc(l,m',m) = { \delta_m,m'  for  m = 0                                 (1)
    !              { i/sqrt(2)(\delta_m,m' - (-1)^m \delta_-m,m'  for m < 0  (2)
    !              { 1/sqrt(2)(\delta_-m,m' + (-1)^m \delta_m,m'  for m > 0  (3)
    allocate(cc(0:4*(lmax-1)),stat=ierr)
    call utils_alloc_check('gaunt_calculate','cc',ierr)
    do lup=0,4*(lmax-1)
       allocate(cc(lup)%val(-lup:lup,-lup:lup),stat=ierr)
       call utils_alloc_check('gaunt_calculate','cc%val',ierr)
       cc(lup)%val(:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)
       cc(lup)%val(0,0) = cmplx(1.0_DP,0.0_DP,kind=DP)    ! (1)
       do mup=-lup,-1                                     ! (2)
          cc(lup)%val(-mup,mup) = cmplx(0.0_DP,-(-1.0_DP)**mup/sqrt(2.0_DP),kind=DP)
          cc(lup)%val(mup,mup) = cmplx(0.0_DP,1.0_DP/sqrt(2.0_DP),kind=DP)
       end do
       do mup=1,lup                                       ! (3)
          cc(lup)%val(mup,mup) = cmplx((-1.0_DP)**mup/sqrt(2.0_DP),0.0_DP,kind=DP)
          cc(lup)%val(-mup,mup) = cmplx(1.0_DP/sqrt(2.0_DP),0.0_DP,kind=DP)
       end do
    end do

    ! Initialisation of count of nonzero coeffs
    n_gaunt = 0

    ! Loop over lup,mup
    do lup=0,lmax-1
       do mup=-lup,lup

          ! Find index in lup,mup arrays
          idxup = lup**2 + lup + mup + 1

          ! Find sgnup
          sgnup = 1
          if ((lup > 0).and.(lup < lmax).and.(mup < 0)) sgnup = -1

          ! Loop over l1,m1
          do l1=0,lmax-1
             do m1=-l1,l1

                ! Find index in l1,m1 arrays
                idx1 = l1**2 + l1 + m1 + 1

                ! Find sgn1
                sgn1 = 1
                if ((l1 > 0).and.(l1 < lmax).and.(m1 < 0)) sgn1 = -1

                ! Need l1,m1 <= lup,mup only
                if (idx1>idxup) cycle
                gaunt_select(:,idxup*(idxup-1)/2+idx1) = 0

                ! Loop over l2,m2 (skipping l2's which are zero by sel. rules)
                do l2=abs(l1-lup),l1+lup,2
                   do m2=-l2,l2

                      ! Find index in l2,m2 arrays
                      idx2 = l2**2 + l2 + m2 + 1

                      ! Apply Real Gaunt coeff selection rules
                      nonzero = .false.
                      if ((l2<=l1+lup).and.&
                           (((m1==mup).and.((m2==0).or.(m2==2*abs(mup)))).or. &
                           ((m1==-mup).and.(m2==-abs(m1)-abs(mup))).or. &
                           ((abs(m1)/=(abs(mup)).and. &
                           ((m2==sgn1*sgnup*(abs(m1)+abs(mup))).or. &
                           (m2==sgn1*sgnup*abs(abs(m1)-abs(mup))) &
                           ))))) nonzero = .true.
                      if (.not.nonzero) cycle

                      ! Compute 'Real' Gaunt coefficient from cc coeffs and
                      ! Gaunt coefficients for normal Spherical Harmonics
                      ! obtained from the rgaunt function

                      ! Loop over m1',mup',m2'
                      real_gaunt_tmp = 0.0_DP
                      do mupp=-lup,lup
                         do m1p=-l1,l1
                            do m2p=-l2,l2

                               ! Apply selection rule on m1',mup',m2'
                               if (m2p /= -mupp-m1p) cycle

                               ! cc(l1,m1',m1)*cc(lup,mup',mup)*cc(l2,m2',m2)
                               ccprod = real(cc(l1)%val(m1p,m1) * &
                                    cc(lup)%val(mupp,mup) * &
                                    cc(l2)%val(m2p,m2),kind=DP)

                               ! Add to total, taking into account signs
                               if ((abs(ccprod) >= 1d-12)) real_gaunt_tmp = &
                                    real_gaunt_tmp + ccprod * (-1)**mupp * &
                                    rgaunt(l2,m2p,l1,m1p,lup,-mupp)

                            end do  ! m2p
                         end do  ! m1p
                      end do  ! mupp

                      ! Store Gaunt coeff in indexed array and keep track of
                      ! number of nonzero coeffs
                      if (abs(real_gaunt_tmp) >= 1d-12) then
                         n_gaunt = n_gaunt + 1
                         real_gaunt(n_gaunt) = real_gaunt_tmp/sqrt(4.0_DP*PI)
                         gaunt_select(idx2,idxup*(idxup-1)/2+idx1) = n_gaunt
                      end if

                   end do  ! m2
                end do  ! l2
             end do  ! m1
          end do  ! l1
       end do  ! mup
    end do  ! lup

    ! Deallocate memory
    do lup=0,4*(lmax-1)
       deallocate(cc(lup)%val,stat=ierr)
       call utils_dealloc_check('gaunt_calculate','cc%val',ierr)
    end do
    deallocate(cc,stat=ierr)
    call utils_dealloc_check('gaunt_calculate','cc',ierr)

  end subroutine gaunt_calculate

end module gaunt_coeff
