! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!      Projector-Augmented Wave Shape Function module            !
!                                                                !
! This module contains routines relating to the shape functions  !
! used in the construction of the compensation charge nhat(r)    !
! within the PAW formalism.                                      !
!----------------------------------------------------------------!
! This module was created by Nicholas Hine in June 2010.         !
!================================================================!

module paw_shape

  use constants, only: DP, stdout

  implicit none

  private

  ! Type describing compensation density shape functions
  type PAW_SHAPE_INFO

     ! ndmh: number of points in grid
     integer :: npts

     ! ndmh: The type of g_L shape function to use
     integer :: shape_type

     ! ndmh: parameters for g_L
     real(kind=DP) :: rshape
     real(kind=DP), pointer :: norm(:)

     ! ndmh: parameters for shape_type:3 Bessel function shapes
     real(kind=DP), pointer :: shape_alpha(:,:)
     real(kind=DP), pointer :: shape_q(:,:)

     ! ndmh: parameters for exponential shape functions
     real(kind=DP), pointer :: shape_lambda(:)
     real(kind=DP), pointer :: shape_sigma(:)

     real(kind=DP), pointer :: shape_rad(:,:)
     real(kind=DP), pointer :: grad_shape_rad(:,:)

  end type PAW_SHAPE_INFO

  public PAW_SHAPE_INFO

  ! Public subroutines
  public :: paw_shape_init
  public :: paw_shape_exit
  public :: paw_shape_gLSLM_real
  public :: paw_shape_gLSLM_recip
  public :: paw_shape_grad_gLSLM_real
  public :: paw_shape_grad_gLSLM_recip
  public :: paw_shape_calculate

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_init(shape,rad,rab,npt,lmax)

    !=====================================================================!
    ! This subroutine initialises the arrays storing the shape functions  !
    ! for each atom evaluated on a radial grid.                           !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !  shape (inout) : PAW_SHAPE_INFO type describing this shape function !
    !  rad (in) : radial grid for the shape function                      !
    !  rab (in) : radial grid "spacings" for the shape function           !
    !  npt (in) : number of points in radial grid                         !
    !  lmax (in) : highest angular momentum requested for shape functions !
    !---------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.              !
    !=====================================================================!

    use services, only: services_radial_integral_rmax
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(PAW_SHAPE_INFO), intent(inout) :: shape
    integer, intent(in) :: lmax
    integer, intent(in) :: npt
    real(kind=DP), intent(in) :: rad(npt)
    real(kind=DP), intent(in) :: rab(npt)

    ! Local Variables
    integer :: ierr
    integer :: l
    real(kind=DP), allocatable :: work(:,:)
    real(kind=DP) :: test_norm
    character(80) :: error

    ! Check arguments
    if ((shape%shape_type<1).or.(shape%shape_type>3)) then
       write(error,'(a,i2,a)') 'Error in paw_shape_init: shape_type = ', &
            shape%shape_type,' is not currently supported'
       call utils_abort(error)
    end if
    shape%npts = npt

    ! Allocate storage for shape function
    allocate(shape%shape_alpha(2,0:lmax),stat=ierr)
    call utils_alloc_check('paw_shape_init','shape%shape_alpha',ierr)
    allocate(shape%shape_q(2,0:lmax),stat=ierr)
    call utils_alloc_check('paw_shape_init','shape%shape_q',ierr)
    allocate(shape%shape_lambda(0:lmax),stat=ierr)
    call utils_alloc_check('paw_shape_init','shape%shape_lambda',ierr)
    allocate(shape%shape_sigma(0:lmax),stat=ierr)
    call utils_alloc_check('paw_shape_init','shape%shape_sigma',ierr)
    allocate(shape%shape_rad(npt,0:lmax),stat=ierr)
    call utils_alloc_check('paw_shape_init','shape%shape_rad',ierr)
    allocate(shape%grad_shape_rad(npt,0:lmax),stat=ierr)
    call utils_alloc_check('paw_shape_init','shape%grad_shape_rad',ierr)
    allocate(shape%norm(0:lmax),stat=ierr)
    call utils_alloc_check('paw_shape_init','shape%norm',ierr)
    allocate(work(npt,2),stat=ierr)
    call utils_alloc_check('paw_shape_init','work',ierr)

    shape%shape_alpha(:,:) = 0.0_DP
    shape%shape_q(:,:) = 0.0_DP
    shape%shape_lambda(:) = 0.0_DP
    shape%shape_sigma(:) = 0.0_DP

    do l=0,lmax

       if (shape%shape_type==3) then

          ! Evaluate parameters of Bessel shape function
          call paw_shape_bessel_params(shape%shape_alpha(:,l), &
               shape%shape_q(:,l),l,shape%rshape)

       else if (shape%shape_type==1) then

          ! Evaluate parameters of exponential shape function
          call paw_shape_expon_params(shape%shape_lambda(l), &
               shape%shape_sigma(l),l,shape%rshape)

       end if

       ! Evaluate values of shape function up to rc
       call paw_shape_calculate(shape%shape_rad(:,l), &
            shape%grad_shape_rad(:,l),rad,npt,l,shape%shape_type, &
            shape%rshape,shape%shape_alpha(:,l),shape%shape_q(:,l), &
            shape%shape_lambda(l),shape%shape_sigma(l))

       ! Find normalisation of shape function
       work(:,1) = shape%shape_rad(:,l) * rad**(l+2)
       test_norm = services_radial_integral_rmax(npt,rab,rad, &
               shape%rshape,work(:,1),work(:,2))
       shape%norm(l) = test_norm

       if ((shape%shape_type==1).or.(shape%shape_type==2)) then
          shape%shape_rad(:,l) = shape%shape_rad(:,l) / shape%norm(l)
          shape%grad_shape_rad(:,l) = shape%grad_shape_rad(:,l) / shape%norm(l)

          ! Find normalisation of shape function again
          work(:,1) = shape%shape_rad(:,l) * rad**(l+2)
          test_norm = services_radial_integral_rmax(npt,rab,rad, &
               shape%rshape,work(:,1),work(:,2))

       end if

       ! Check normalisation integrates sufficiently accurately to 1 up to rc
       if ((test_norm > 1.0001_DP).or.(test_norm<0.9999_DP)) then
          write(error,'(a,i1,a)') 'Error in paw_shape_init: Shape function &
               &for l=',l,' is not properly normalised'
          call utils_abort(error)
       end if

    end do

    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_shape_init','work',ierr)

  end subroutine paw_shape_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_exit(shape)

    !====================================================================!
    ! This subroutine deallocates arrays relating to the shape functions !
    ! in the PAW_SHAPE_INFO type.                                        !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  shape (inout) : PAW_SHAPE_INFO type describing this shape function!
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.             !
    !====================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PAW_SHAPE_INFO), intent(inout) :: shape

    ! Local Variables
    integer :: ierr

    ! Deallocate storage for shape function
    deallocate(shape%norm,stat=ierr)
    call utils_dealloc_check('paw_shape_exit','shape%norm',ierr)
    deallocate(shape%grad_shape_rad,stat=ierr)
    call utils_dealloc_check('paw_shape_exit','shape%grad_shape_rad',ierr)
    deallocate(shape%shape_rad,stat=ierr)
    call utils_dealloc_check('paw_shape_exit','shape%shape_rad',ierr)
    deallocate(shape%shape_q,stat=ierr)
    call utils_dealloc_check('paw_shape_init','shape%shape_q',ierr)
    deallocate(shape%shape_alpha,stat=ierr)
    call utils_dealloc_check('paw_shape_init','shape%shape_alpha',ierr)
    deallocate(shape%shape_sigma,stat=ierr)
    call utils_dealloc_check('paw_shape_init','shape%shape_sigma',ierr)
    deallocate(shape%shape_lambda,stat=ierr)
    call utils_dealloc_check('paw_shape_init','shape%shape_lambda',ierr)

  end subroutine paw_shape_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_gLSLM_real(gLSLM,shape,lup,mup, &
      cell_start1,cell_start2,cell_start3,grid,box_n1,box_n2,box_n3, &
      atom_origin)

    !=====================================================================!
    ! This subroutine calculates the function g_L(r)*S_LM(\hat{r})        !
    ! which is required to form the compensation charge \hat{n}.          !
    ! It does so on the simulation cell fine grid in a box centered on    !
    ! the atom.                                                           !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !  gLSLM (out) : real-space box for shape function for a given L,M    !
    !  shape (inout) : PAW_SHAPE_INFO type describing this shape function !
    !  lup (in) : angular momentum L                                      !
    !  mup (in) : azimuthal angular momentum M                            !
    !  cell_start1 (in) : grid point index at which to start box in 1-dir !
    !  cell_start2 (in) : grid point index at which to start box in 2-dir !
    !  cell_start3 (in) : grid point index at which to start box in 3-dir !
    !  grid (in) : GRID_INFO type describing grid for charge density      !
    !  box_n1 (in) : size of augmentation box in grid points in 1-dir     !
    !  box_n2 (in) : size of augmentation box in grid points in 2-dir     !
    !  box_n3 (in) : size of augmentation box in grid points in 3-dir     !
    !  atom_origin (in) : POINT type describing location of atom          !
    !---------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.              !
    !=====================================================================!

    use basis, only: basis_box_origin_to_atom
    use cell_grid, only: GRID_INFO
    use constants, only: PI, cmplx_i
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_real_sph_harm, sw_bessel_accurate, sw_init

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    real(kind=DP),intent(out) :: gLSLM(box_n1,box_n2,box_n3)
    type(PAW_SHAPE_INFO), intent(in) :: shape
    integer,intent(in) :: lup
    integer,intent(in) :: mup
    integer,intent(in) :: cell_start1
    integer,intent(in) :: cell_start2
    integer,intent(in) :: cell_start3
    type(GRID_INFO),intent(in) :: grid
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: i1,i2,i3
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    type(POINT) :: r_cell
    type(POINT) :: r_sphere
    real(kind=DP) :: rmag,q(1:2)
    real(kind=DP) :: sinc
    real(kind=DP) :: gval
    real(kind=DP) :: slmval
    real(kind=DP) :: j(1:2),jp(1:2)

    if (lup>4) call sw_init(lup+1,1)

    call basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
         cell_start1,cell_start2,cell_start3,grid%da1,grid%da2,grid%da3)

    do i3=1,box_n3
       do i2=1,box_n2
          do i1=1,box_n1

             r_cell = box_origin + (i1-1)*grid%da1 &
                                 + (i2-1)*grid%da2 &
                                 + (i3-1)*grid%da3

             r_sphere = r_cell - atom_origin
             rmag = geometry_magnitude(r_sphere)

             if (rmag>shape%rshape) then
                gLSLM(i1,i2,i3) = 0.0_DP
                cycle
             end if

             if (shape%shape_type==3) then
                q(1:2) = shape%shape_q(1:2,lup)
                call sw_bessel_accurate(lup,rmag*q(1),j(1),jp(1))
                call sw_bessel_accurate(lup,rmag*q(2),j(2),jp(2))
                gval = shape%shape_alpha(1,lup)*j(1) + &
                     shape%shape_alpha(2,lup)*j(2)
             else if (shape%shape_type==2) then
                if (abs(PI*rmag/shape%rshape)<1e-10_DP) then
                   sinc = 1.0_DP
                else
                   sinc = sin(PI*rmag/shape%rshape)/(PI*rmag/shape%rshape)
                end if
                gval = sinc**2 * rmag**lup / shape%norm(lup)
             else if (shape%shape_type==1) then

             end if

             slmval = sw_real_sph_harm(r_sphere%x, &
                  r_sphere%y,r_sphere%z,rmag,lup,mup)

             gLSLM(i1,i2,i3) = slmval * gval

          end do
       end do
    end do

  end subroutine paw_shape_gLSLM_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_gLSLM_recip(gLSLM_recip,shape,lup,mup, &
      cell_start1,cell_start2,cell_start3,grid,box_n1,box_n2,box_n3, &
      atom_origin)

    !=====================================================================!
    ! This subroutine calculates the function g_L(r)*S_LM(\hat{r})        !
    ! which is required to form the compensation charge \hat{n}.          !
    ! It does so on the simulation cell fine grid in a box centered on    !
    ! the atom.                                                           !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !  gLSLM_recip (out) : recip-space box for shape func for this L,M    !
    !  shape (inout) : PAW_SHAPE_INFO type describing this shape function !
    !  lup (in) : angular momentum L                                      !
    !  mup (in) : azimuthal angular momentum M                            !
    !  cell_start1 (in) : grid point index at which to start box in 1-dir !
    !  cell_start2 (in) : grid point index at which to start box in 2-dir !
    !  cell_start3 (in) : grid point index at which to start box in 3-dir !
    !  grid (in) : GRID_INFO type describing grid for charge density      !
    !  box_n1 (in) : size of augmentation box in grid points in 1-dir     !
    !  box_n2 (in) : size of augmentation box in grid points in 2-dir     !
    !  box_n3 (in) : size of augmentation box in grid points in 3-dir     !
    !  atom_origin (in) : POINT type describing location of atom          !
    !---------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.              !
    !=====================================================================!

    use basis, only: basis_box_origin_to_atom
    use cell_grid, only: GRID_INFO
    use constants, only: PI, cmplx_i
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_real_sph_harm, sw_bessel_accurate, sw_init

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    complex(kind=DP),intent(out) :: gLSLM_recip(box_n1,box_n2,box_n3)
    type(PAW_SHAPE_INFO), intent(in) :: shape
    integer,intent(in) :: lup
    integer,intent(in) :: mup
    integer,intent(in) :: cell_start1
    integer,intent(in) :: cell_start2
    integer,intent(in) :: cell_start3
    type(GRID_INFO),intent(in) :: grid
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: i1,i2,i3
    integer :: k1,k2,k3
    integer :: ii
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    real(kind=DP) :: rmag,q(1:2)
    real(kind=DP) :: gval
    real(kind=DP) :: slmval
    real(kind=DP) :: j(1:2),jp(1:2),jq,jpq
    type(POINT) :: g_vector3,g_vector2,g_vector1
    type(POINT) :: pairvec3,pairvec2,pairvec1
    real(kind=DP) :: g_length
    complex(kind=DP) :: phase_fac

    if (lup>4) call sw_init(lup+1,1)

    call basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
         cell_start1,cell_start2,cell_start3,grid%da1,grid%da2,grid%da3)

    ! initialisations
    gLSLM_recip = 0.0_DP
    pairvec3=( real(grid%n3,DP)/real(box_n3,DP) ) * pub_cell%b3
    pairvec2=( real(grid%n2,DP)/real(box_n2,DP) ) * pub_cell%b2
    pairvec1=( real(grid%n1,DP)/real(box_n1,DP) ) * pub_cell%b1
    phase_fac = cmplx(0.0_DP,-1.0_DP)**lup
    rmag = shape%rshape

    if (shape%shape_type==3) then

       q(1:2) = shape%shape_q(1:2,lup)
       if (lup>0) then
          call sw_bessel_accurate(lup-1,rmag*q(1),j(1),jp(1))
          call sw_bessel_accurate(lup-1,rmag*q(2),j(2),jp(2))
          ! store j(lup,q_i*r_c)*q_i in j_i
          j(1:2) = j(1:2) * q(1:2)
       else
          j(1) = cos(rmag*q(1)) / rmag
          j(2) = cos(rmag*q(2)) / rmag
       end if
    end if

    ! loop over reciprocal-space grid points of augmentation box
    do i3=1,box_n3
       k3 = i3 - 1
       if (k3>box_n3/2) k3 = k3 - box_n3
       g_vector3 = real(k3,kind=DP)*pairvec3

       do i2=1,box_n2
          k2 = i2 - 1
          if (k2>box_n2/2) k2 = k2 - box_n2
          g_vector2 = real(k2,kind=DP)*pairvec2 + g_vector3

          do i1=1,box_n1
             k1 = i1 - 1
             if (k1>box_n1/2) k1 = k1 - box_n1
             g_vector1 = real(k1,kind=DP)*pairvec1 + g_vector2

             g_length = geometry_magnitude(g_vector1)

             ! Calculate sw transform of shape function at this G
             if (shape%shape_type==3) then
                call sw_bessel_accurate(lup,rmag*g_length,jq,jpq)

                gval = 0.0_DP
                do ii=1,2
                   if (abs(q(ii)-g_length)>epsilon(1.0_DP)) then
                      gval = gval + shape%shape_alpha(ii,lup) * &
                           shape%rshape**2 * jq * j(ii) / &
                           (g_length**2-q(ii)**2)
                   else
                      gval = gval + shape%shape_alpha(ii,lup) * &
                           shape%rshape**3 * j(ii)**2 / q(ii) / &
                           (g_length + q(ii))
                   end if
                end do
             else if (shape%shape_type==2) then
                gval = 0.0_DP
             else if (shape%shape_type==1) then
                gval = 0.0_DP
             end if

             ! multiply value of shapefunc at each g-point by the appropriate
             ! phase factor (n.b. real part and imaginary part stored as
             ! separate consecutive elements in the x-direction of the array),
             ! and the appropriate real spherical harmonic factor.
             slmval = 4.0_DP*PI * sw_real_sph_harm(g_vector1%x, &
                  g_vector1%y, g_vector1%z, g_length, lup, mup)

             gLSLM_recip(i1,i2,i3) = gLSLM_recip(i1,i2,i3) + &
                  gval * slmval * phase_fac * &
                  exp(-cmplx_i*(g_vector1.dot.box_to_atom))

          end do
       end do
    end do

  end subroutine paw_shape_gLSLM_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_grad_gLSLM_real(grad_gLSLM,shape,lup,mup, &
      cell_start1,cell_start2,cell_start3,grid,box_n1,box_n2,box_n3, &
      atom_origin)

    !=====================================================================!
    ! This subroutine calculates the gradient of the function             !
    ! g_L(r)*S_LM(\hat{r}) in each of the three cartesian directions.     !
    ! This is required to form the forces on each atom resulting from     !
    ! the interaction of the compensation charge \hat{n} with the         !
    ! effective potential. It does so on the simulation cell fine grid    !
    ! in a box centered on the atom.                                      !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !  grad_gLSLM (out) : gradients of shape function for a given L,M     !
    !  shape (inout) : PAW_SHAPE_INFO type describing this shape function !
    !  lup (in) : angular momentum L                                      !
    !  mup (in) : azimuthal angular momentum M                            !
    !  cell_start1 (in) : grid point index at which to start box in 1-dir !
    !  cell_start2 (in) : grid point index at which to start box in 2-dir !
    !  cell_start3 (in) : grid point index at which to start box in 3-dir !
    !  grid (in) : GRID_INFO type describing grid for charge density      !
    !  box_n1 (in) : size of augmentation box in grid points in 1-dir     !
    !  box_n2 (in) : size of augmentation box in grid points in 2-dir     !
    !  box_n3 (in) : size of augmentation box in grid points in 3-dir     !
    !  atom_origin (in) : POINT type describing location of atom          !
    !---------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.              !
    !=====================================================================!

    use basis, only: basis_box_origin_to_atom
    use cell_grid, only: GRID_INFO
    use constants, only: PI
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_real_sph_harm, sw_grad_real_sph_harm, &
         sw_bessel_accurate, sw_init

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    real(kind=DP),intent(out) :: grad_gLSLM(box_n1,box_n2,box_n3,3)
    type(PAW_SHAPE_INFO), intent(in) :: shape
    integer,intent(in) :: lup
    integer,intent(in) :: mup
    integer,intent(in) :: cell_start1
    integer,intent(in) :: cell_start2
    integer,intent(in) :: cell_start3
    type(GRID_INFO), intent(in) :: grid
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: i1,i2,i3
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    type(POINT) :: r_cell
    type(POINT) :: r_sphere
    type(POINT) :: r_sphere_unit
    real(kind=DP) :: sinc,gradsinc
    real(kind=DP) :: rmag
    real(kind=DP) :: grad_gval
    real(kind=DP) :: gval
    real(kind=DP) :: slmval
    real(kind=DP) :: grad_SLM_R(3)
    real(kind=DP) :: j1,j2,jp1,jp2

    if ((lup+1)>4) call sw_init(lup+1,1)

    call basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
         cell_start1,cell_start2,cell_start3,grid%da1,grid%da2,grid%da3)

    do i3=1,box_n3
       do i2=1,box_n2
          do i1=1,box_n1

             r_cell = box_origin + (i1-1)*grid%da1 &
                                 + (i2-1)*grid%da2 &
                                 + (i3-1)*grid%da3

             r_sphere = r_cell - atom_origin
             rmag = geometry_magnitude(r_sphere)

             if (rmag>shape%rshape) then
                grad_gLSLM(i1,i2,i3,:) = 0.0_DP
                cycle
             end if

             if (rmag>1d-10) then
                r_sphere_unit = (1.0_DP/rmag)*r_sphere
             else
                r_sphere_unit = 0.0_DP*r_sphere
             end if

             if (shape%shape_type==3) then
                ! Spherical Bessel Shape Function
                if (rmag>1d-10) then
                   call sw_bessel_accurate(lup,rmag*shape%shape_q(1,lup),j1,jp1)
                   call sw_bessel_accurate(lup,rmag*shape%shape_q(2,lup),j2,jp2)
                   gval = shape%shape_alpha(1,lup)*j1 + &
                        shape%shape_alpha(2,lup)*j2
                   grad_gval = &
                        shape%shape_q(1,lup)*shape%shape_alpha(1,lup)*jp1 + &
                        shape%shape_q(2,lup)*shape%shape_alpha(2,lup)*jp2
                else
                   if (lup==1) then ! special value for L=1: lim_{r->0} g(r)/r
                      gval = (shape%shape_alpha(1,lup) * shape%shape_q(1,lup) +&
                           shape%shape_alpha(2,lup) * shape%shape_q(2,lup)) &
                           / 3.0_DP
                      grad_gval = 0.0_DP
                   else
                      gval = 0.0_DP
                      grad_gval = 0.0_DP
                   end if
                end if
             else if (shape%shape_type==2) then
                ! Sinc-Squared Shape Function
                if (rmag>1d-10) then
                   sinc = sin(PI*rmag/shape%rshape)/(PI*rmag/shape%rshape)
                   gval = sinc**2 * rmag**lup / shape%norm(lup)
                   gradsinc = (cos(PI*rmag/shape%rshape)-sinc)/rmag
                   grad_gval = lup * rmag**(lup-1) * sinc**2 + &
                         2.0_DP * rmag**lup * sinc * gradsinc
                else
                   if (lup<=1) then ! special value for L=1: lim_{r->0} g(r)/r
                      gval = 1.0_DP / shape%norm(lup)
                      grad_gval = 0.0_DP
                   else
                      gval = 0.0_DP
                      grad_gval = 0.0_DP
                   end if
                end if
             else if (shape%shape_type==1) then
                ! Negative-Exponential Shape Function
                ! TODO
             end if

             ! Get spherical harmonic at this \hat{r}
             slmval = sw_real_sph_harm(r_sphere%x,r_sphere%y,r_sphere%z, &
                  rmag,lup,mup)

             ! Get gradients of spherical harmonic at this \hat{r}
             call sw_grad_real_sph_harm(grad_SLM_R(:),r_sphere%x,r_sphere%y, &
                  r_sphere%z,rmag,lup,mup)

             ! S_LM(\hat{r}) * d/dR(g(|r-R|) part
             grad_gLSLM(i1,i2,i3,1) = grad_gval * slmval * r_sphere_unit%x
             grad_gLSLM(i1,i2,i3,2) = grad_gval * slmval * r_sphere_unit%y
             grad_gLSLM(i1,i2,i3,3) = grad_gval * slmval * r_sphere_unit%z

             ! g(|r-R|) * d S_LM(\hat{r}) / dR(g(|r-R|) part
             grad_gLSLM(i1,i2,i3,:) = grad_gLSLM(i1,i2,i3,:) + &
                  grad_SLM_R(:) * gval

          end do
       end do
    end do

  end subroutine paw_shape_grad_gLSLM_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_grad_gLSLM_recip(grad_gLSLM,shape,lup,mup, &
      cell_start1,cell_start2,cell_start3,grid,box_n1,box_n2,box_n3, &
      atom_origin)

    !=====================================================================!
    ! This subroutine calculates the gradient of the function             !
    ! g_L(r)*S_LM(\hat{r}) in each of the three cartesian directions.     !
    ! This is required to form the forces on each atom resulting from     !
    ! the interaction of the compensation charge \hat{n} with the         !
    ! effective potential. It does so on the simulation cell fine grid    !
    ! in a box centered on the atom.                                      !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !  gLSLM_recip (out) : recip-space box for shape func for this L,M    !
    !  shape (inout) : PAW_SHAPE_INFO type describing this shape function !
    !  lup (in) : angular momentum L                                      !
    !  mup (in) : azimuthal angular momentum M                            !
    !  cell_start1 (in) : grid point index at which to start box in 1-dir !
    !  cell_start2 (in) : grid point index at which to start box in 2-dir !
    !  cell_start3 (in) : grid point index at which to start box in 3-dir !
    !  grid (in) : GRID_INFO type describing grid for charge density      !
    !  box_n1 (in) : size of augmentation box in grid points in 1-dir     !
    !  box_n2 (in) : size of augmentation box in grid points in 2-dir     !
    !  box_n3 (in) : size of augmentation box in grid points in 3-dir     !
    !  atom_origin (in) : POINT type describing location of atom          !
    !---------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.              !
    !=====================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: cmplx_i
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: pub_cell

    use fourier, only: fourier_apply_box

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    complex(kind=DP),intent(out) :: grad_gLSLM(box_n1,box_n2,box_n3,3)
    type(PAW_SHAPE_INFO), intent(in) :: shape
    integer,intent(in) :: lup
    integer,intent(in) :: mup
    integer,intent(in) :: cell_start1
    integer,intent(in) :: cell_start2
    integer,intent(in) :: cell_start3
    type(GRID_INFO), intent(in) :: grid
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: i1,i2,i3
    integer :: k1,k2,k3
    type(POINT) :: g_vector3,g_vector2,g_vector1
    type(POINT) :: pairvec3,pairvec2,pairvec1

    ! Get shape function in reciprocal space into third grad array
    call paw_shape_gLSLM_recip(grad_gLSLM(:,:,:,3),shape,lup,mup, &
         cell_start1,cell_start2,cell_start3,grid,box_n1,box_n2,box_n3, &
         atom_origin)

    ! initialisations
    pairvec3=( real(grid%n3,DP)/real(box_n3,DP) ) * pub_cell%b3
    pairvec2=( real(grid%n2,DP)/real(box_n2,DP) ) * pub_cell%b2
    pairvec1=( real(grid%n1,DP)/real(box_n1,DP) ) * pub_cell%b1

    ! loop over reciprocal-space grid points of augmentation box
    do i3=1,box_n3
       k3 = i3 - 1
       if (k3>box_n3/2) k3 = k3 - box_n3
       g_vector3 = real(k3,kind=DP)*pairvec3

       do i2=1,box_n2
          k2 = i2 - 1
          if (k2>box_n2/2) k2 = k2 - box_n2
          g_vector2 = real(k2,kind=DP)*pairvec2 + g_vector3

          do i1=1,box_n1
             k1 = i1 - 1
             if (k1>box_n1/2) k1 = k1 - box_n1
             g_vector1 = real(k1,kind=DP)*pairvec1 + g_vector2

             grad_gLSLM(i1,i2,i3,1) = &
                  grad_gLSLM(i1,i2,i3,3)*g_vector1%x * cmplx_i
             grad_gLSLM(i1,i2,i3,2) = &
                  grad_gLSLM(i1,i2,i3,3)*g_vector1%y * cmplx_i
             grad_gLSLM(i1,i2,i3,3) = &
                  grad_gLSLM(i1,i2,i3,3)*g_vector1%z * cmplx_i

          end do
       end do
    end do

  end subroutine paw_shape_grad_gLSLM_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_calculate(shape_rad,grad_shape_rad,rad,npt,lup, &
       shape_type,rcut,alpha,q,lambda,sigma)

    !==================================================================!
    ! This subroutine evaluates the shape function on a radial grid.   !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  shape_rad (out) : shape function on radial grid                 !
    !  grad_shape_rad (out) : gradient of shape function on radial grid!
    !  rad (in) : radial grid points                                   !
    !  npt (in) : number of radial grid points                         !
    !  lup (in) : angular momentum L                                   !
    !  shape_type (in) : integer describing type of shape function     !
    !  rcut (in) : shape function radius (usually same as PAW rcut)    !
    !  alpha (in) : parameter for Bessel shape functions               !
    !  q (in) : parameter for Bessel shape functions                   !
    !  lambda (in) : parameter for exponential shape functions         !
    !  sigma (in) : parameter for exponential shape functions          !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    !==================================================================!

    use constants, only: PI
    use spherical_wave, only: sw_bessel_accurate

    implicit none

    ! Arguments
    integer, intent(in) :: npt
    real(kind=DP), intent(in) :: rad(npt)
    real(kind=DP), intent(out) :: shape_rad(npt)
    real(kind=DP), intent(out) :: grad_shape_rad(npt)
    integer, intent(in) :: lup
    integer, intent(in) :: shape_type
    real(kind=DP), intent(in) :: alpha(2), q(2)
    real(kind=DP), intent(in) :: sigma, lambda
    real(kind=DP), intent(in) :: rcut

    ! Local Variables
    integer :: i
    real(kind=DP) :: j1,j2,jp1,jp2
    real(kind=DP) :: sinc, gradsinc

    shape_rad(:) = 0.0_DP
    grad_shape_rad(:) = 0.0_DP

    select case (shape_type)

    ! Exponential shape
    case (1)

      do i=1,npt

         shape_rad(i) = exp(-lambda*rad(i)/sigma) * rad(i)**lup
         grad_shape_rad(i) = -lambda * exp(-lambda*rad(i)/sigma) &
              / sigma * rad(i)**lup
         if (lup > 0) grad_shape_rad(i) = grad_shape_rad(i) + &
              real(lup,kind=DP) * exp(-lambda*rad(i)/sigma) * rad(i)**(lup-1)
         if (rad(i)>rcut) then
            shape_rad(i) = 0.0_DP
            grad_shape_rad(i) = 0.0_DP
         end if

      end do

    ! Psinc function shapes
    case (2)

      do i=1,npt-1

         if (rad(i)>0.0_DP) then
            sinc = (sin(PI*rad(i)/rcut)/(PI*rad(i)/rcut))
            gradsinc = (cos(PI*rad(i)/rcut)-sinc)/rad(i)

            shape_rad(i) = rad(i)**lup * sinc**2
            if (rad(i+1) > rcut) shape_rad(i) = 0.0_DP

            grad_shape_rad(i) = lup * rad(i)**(lup-1) * sinc**2 + &
                 2.0_DP * rad(i)**lup * sinc * gradsinc
            if (rad(i+1) > rcut) grad_shape_rad(i) = 0.0_DP
         else
            ! ddor: Accommodate zero l and r
            if (lup==0) then
               shape_rad(i) = 1.0_DP
               grad_shape_rad(i) = 0.0_DP
            else
               shape_rad(i) = 0.0_DP
               grad_shape_rad(i) = 0.0_DP
            endif
         end if

      end do

    ! Bessel function shapes
    case (3)

      do i=1,npt-1
         call sw_bessel_accurate(lup,rad(i)*q(1),j1,jp1)
         call sw_bessel_accurate(lup,rad(i)*q(2),j2,jp2)
         shape_rad(i) = alpha(1)*j1 + alpha(2)*j2
         if (rad(i+1) > rcut) shape_rad(i) = 0.0_DP
         grad_shape_rad(i) = alpha(1)*q(1)*jp1 + alpha(2)*q(2)*jp2
         if (rad(i+1) > rcut) grad_shape_rad(i) = 0.0_DP
      end do

    end select

    shape_rad(npt) = 0.0_DP
    grad_shape_rad(npt) = 0.0_DP

  end subroutine paw_shape_calculate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_expon_params(lambda_l,sigma_l,lup,rcut)

    !==================================================================!
    ! This subroutine calculates the parameters of a shape function    !
    ! expressed as a negative exponential.                             !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  lambda (in) : parameter for exponential shape functions         !
    !  sigma (in) : parameter for exponential shape functions          !
    !  lup (in) : angular momentum L                                   !
    !  rcut (in) : shape function radius (usually same as PAW rcut)    !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 19/07/10.           !
    !==================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP),intent(out) :: lambda_l
    real(kind=DP),intent(out) :: sigma_l
    integer,intent(in) :: lup
    real(kind=DP),intent(in) :: rcut

    if (lup==1) lambda_l = 0.0_DP
    if (rcut>0.0_DP) sigma_l = 0.0_DP

    call utils_abort('Error in paw_shape_init: shape_type = 1 &
         &is not currently supported')

  end subroutine paw_shape_expon_params


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_bessel_params(alpha_l,q_l,lup,rcut)

    !==================================================================!
    ! This subroutine calculates the parameters of a shape function    !
    ! expressed as the sum of two Bessel functions.                    !
    ! Ensures that:   1) g_L(r=rc) = 0                                 !
    !                 2) d/dr g_L'(r=rc) = 0                           !
    !                 3) \int_0^rc r^L+2 g_L(r) dr = 1                 !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  alpha (in) : parameter for Bessel shape functions               !
    !  q (in) : parameter for Bessel shape functions                   !
    !  lup (in) : angular momentum L                                   !
    !  rcut (in) : shape function radius (usually same as PAW rcut)    !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    ! Mimics behaviour of ABINIT routine shapebes                      !
    !==================================================================!

    use spherical_wave, only: sw_bessel_accurate

    implicit none

    ! Arguments
    real(kind=DP),intent(out) :: alpha_l(2)
    real(kind=DP),intent(out) :: q_l(2)
    integer,intent(in) :: lup
    real(kind=DP),intent(in) :: rcut

    ! Local Variables
    integer :: i
    real(kind=DP) :: D
    real(kind=DP) :: jl,jlp
    real(kind=DP) :: amat(2,2),b(2) ! For solving simultaneous eq

    ! Solve for first two solutions q1 and q2 of j_l(qi*rcut) = 0
    ! Condition (1) is then ensured
    call paw_shape_solve_bessel(q_l,1.0_DP,0.0_DP,lup,2)
    q_l(1:2) = q_l(1:2)/rcut

    ! Evaluate: q1*j_l'(q1*rcut) = d/dr (j_l(q1r))_r=rc
    !           q2*j_l'(q2*rcut) = d/dr (j_l(q2r))_r=rc
    !           j_(l+1)(q1*rcut)*rcut**(l+2)/q1 = \int_0^rc j_l(q1r) r^l+2 dr
    !           j_(l+1)(q2*rcut)*rcut**(l+2)/q2 = \int_0^rc j_l(q2r) r^l+2 dr
    do i=1,2
       call sw_bessel_accurate(lup,q_l(i)*rcut,jl,jlp)
       amat(1,i) = jlp*q_l(i)
       call sw_bessel_accurate(lup+1,q_l(i)*rcut,jl,jlp)
       amat(2,i) = jl*rcut**(lup+2)/q_l(i)
    end do

    ! Solution vector for conditions (2) and (3)
    b(1) = 0.0_DP
    b(2) = 1.0_DP

    ! Solve simultaneous equations A.(alpha1,alpha2) = (0,1)
    ! for alpha1 and alpha2 so that conditions (2) and (3) are then satisfied
    D = amat(1,1)*amat(2,2) - amat(1,2)*amat(2,1)
    alpha_l(1) = (amat(2,2)*b(1) - amat(1,2)*b(2))/D
    alpha_l(2) = (amat(1,1)*b(2) - amat(2,1)*b(1))/D

  end subroutine paw_shape_bessel_params


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_shape_solve_bessel(zeros,alpha,beta,lup,nzeros)

    !==================================================================!
    ! This subroutine finds the lowermost solutions x to the equation  !
    !    f(x) = 0    where     f(x) = alpha*j_l(x)+beta*x*j_l'(x)      !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  zeros (out) : zeros of function                                 !
    !  alpha (in) : coefficient of bessel function                     !
    !  beta (in) : coefficient of bessel function gradient times x     !
    !  lup (in) : angular momentum l                                   !
    !  nzeros (in) : number of zeros requested                         !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    ! Mimics behaviour of the ABINIT routine solvbes                   !
    !==================================================================!

    use spherical_wave, only: sw_bessel_accurate

    implicit none

    ! Arguments
    integer, intent(in) :: lup
    integer, intent(in) :: nzeros
    real(kind=DP), intent(in) :: alpha,beta
    real(kind=DP), intent(out) :: zeros(nzeros)

    ! Local Variables
    integer :: izero
    real(kind=DP) :: bdx
    real(kind=DP) :: jl,jlp
    real(kind=DP) :: f1,f2
    real(kind=DP) :: x,xp
    real(kind=DP),parameter :: dx = 0.1_dp
    real(kind=DP),parameter :: tol = 1e-12_DP

    ! Starting value of x
    x = dx

    ! Loop until enough solutions found
    do izero=1,nzeros

       ! f1 = alpha*j_l(x) + beta*x*j_l'(x)
       call sw_bessel_accurate(lup,x,jl,jlp)
       f1 = alpha * jl + beta * x * jlp

       ! f2 = alpha*j_l(x+dx) + beta*x*j_l'(x+dx)
       x = x + dx
       call sw_bessel_accurate(lup,x,jl,jlp)
       f2 = alpha * jl + beta * x * jlp

       ! Keep evaluating f2 = alpha*j_l(x+dx)+beta*x*j_l'(x+dx) and increasing
       ! dx until interval containing the solution is found
       f2 = f1
       do while (f1*f2 >= 0.0_DP)
          x = x + dx
          call sw_bessel_accurate(lup,x,jl,jlp)
          f2 = alpha * jl + beta * x * jlp
       end do

       ! Bisection search until solution f2=0 is found to within tol
       bdx = dx
       xp = x
       do while (bdx > tol)
          bdx = 0.5_DP * bdx

          ! Check which side of the answer this evaluation lies and update x'
          ! accordingly
          if (f1*f2 < 0.0_DP) then
             xp = xp - bdx
          else
             xp = xp + bdx
          end if

          ! Update f2 with new value
          call sw_bessel_accurate(lup,xp,jl,jlp)
          f2 = alpha * jl + beta * xp * jlp

       end do

       ! Store this solution
       zeros(izero) = xp

    end do

  end subroutine paw_shape_solve_bessel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module paw_shape
