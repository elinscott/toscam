! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                                                                             !
! This module contains data structures and routines to perform FFTs within    !
! the ONETEP code                                                             !
!                                                                             !
! Original version by Peter Haynes, 7/7/03                                    !
! Parallel version by Peter Haynes, 21/6/04                                   !
!                                                                             !
! Modified by Nicholas Hine, July 2008 to improve cache performance, reduce   !
! memory usage, and remove unnecessary memory-copying. Also removed           !
! overloaded interface to fourier_apply to make calls explicit                !
!=============================================================================!

module fourier

  use constants, only: DP, LONG

!CW
#ifdef FFTW3GPU
 use fortran_cuda
#endif
#define GPU_ALL____
#define GPU_ALL2___
!ENDCW

  implicit none

  private

  public :: fourier_size
  public :: fourier_init_cell
  public :: fourier_exit_cell
  public :: fourier_apply_cell_forward
  public :: fourier_apply_cell_backward
  public :: fourier_filter_cell
  public :: fourier_interpolate_cell

  public :: fourier_init
  public :: fourier_apply_box
  public :: fourier_apply_box_pair
  public :: fourier_exit
  public :: fourier_filter
  public :: fourier_interpolate
  public :: fourier_interpolate_product

  ! *********************************
  ! ***  P r i v a t e   d a t a  ***
  ! *********************************

#ifdef DEBUG
  interface fourier_write
     module procedure fourier_write_real_2
     module procedure fourier_write_real_3
     module procedure fourier_write_complex_2
     module procedure fourier_write_complex_3
  end interface
#endif

  type serial_fft3d_info
     integer :: n1,n2,n3  ! The dimensions of the FFT grid
     integer :: ld1,ld2   ! The dimensions of the arrays holding the grid
     integer :: iwork_len ! Length of integer workspace
     integer :: dwork_len ! Length of real workspace
     integer :: zwork_len ! Length of complex workspace
     integer, allocatable, dimension(:) :: iwork ! integer workspace
     real(kind=DP), allocatable, dimension(:) :: dwork ! real workspace
     complex(kind=DP), allocatable, dimension(:) :: zwork ! complex workspace
#ifdef FFTW
     integer(kind=LONG) :: forward_plan, backward_plan
#endif
#ifdef FFTW3
     integer(kind=LONG) :: forward_plan, backward_plan
     integer(kind=LONG) :: n1_f_plan, n1_b_plan
     integer(kind=LONG) :: n2_f_plan, n2_b_plan
     integer(kind=LONG) :: n3_f_plan, n3_b_plan
     integer(kind=LONG) :: n3_f_plan_1, n3_b_plan_1
     integer(kind=LONG) :: n3_f_plan_2, n3_b_plan_2
     integer(kind=LONG) :: n3_f_plan_3, n3_b_plan_3
#endif
!CW
#ifdef FFTW3GPU
     integer(kind=LONG) :: forward_plan, backward_plan
#endif
!END CW
#ifdef ACML
     integer :: tbl_size
     complex(kind=DP), allocatable, dimension(:) :: tbl
#endif
  end type serial_fft3d_info

  type parallel_fft3d_info
     integer :: n1,n2,n3     ! The dimensions of the whole FFT grid
     integer :: ld1,ld2,ld3  ! The dimensions of the arrays holding the grids
     integer :: iwork_len    ! Length of integer workspaces
     integer :: dwork_len    ! Length of real workspaces
     integer :: zwork_len    ! Length of complex workspaces
     integer :: num12slabs   ! Number of 12-slabs for real space
     integer :: max12slabs   ! Maximum number of 12-slabs
     integer :: num23slabs   ! Number of 23-slabs for real space
     integer :: max23slabs   ! Maximum number of 23-slabs
     integer, allocatable, dimension(:) :: idx12slab ! index of 12-slabs
     integer, allocatable, dimension(:) :: idx23slab ! index of 23-slabs
     complex(kind=DP), allocatable, dimension(:,:) :: buf12slab ! buffer for 12-slab
     complex(kind=DP), allocatable, dimension(:,:) :: buf23slab ! buffer for 23-slab
     complex(kind=DP), allocatable, dimension(:,:,:) :: fsendbuf ! forward send buf
     complex(kind=DP), allocatable, dimension(:,:,:) :: frecvbuf ! forward send buf
     complex(kind=DP), allocatable, dimension(:,:,:) :: bsendbuf ! backw'd recv buf
     complex(kind=DP), allocatable, dimension(:,:,:) :: brecvbuf ! forward recv buf
     integer, allocatable, dimension(:,:) :: iwork ! integer workspaces
     real(kind=DP), allocatable, dimension(:,:) :: dwork ! real workspaces
     complex(kind=DP), allocatable, dimension(:,:) :: zwork ! complex workspaces
#ifdef FFTW
     integer(kind=LONG) :: forward_plan(2), backward_plan(2)
     real(kind=DP), allocatable, dimension(:,:) :: rfftwbuf ! buffer for rFFTw
#endif
#ifdef FFTW3
     integer(kind=LONG) :: forward_plan(2), backward_plan(2)
#endif
!CW
#ifdef FFTW3GPU
     integer(kind=LONG) :: forward_plan(2), backward_plan(2)
#endif
!END CW
#ifdef ACML
     integer :: tbl_size
     real(kind=DP), allocatable, dimension(:) :: rbuf ! buffer for zd/dz fft
     complex(kind=DP), allocatable, dimension(:,:) :: tbl
#endif
  end type parallel_fft3d_info

  ! qoh: Information arrays for tightboxes
  type(serial_fft3d_info) :: internal_tb_coarse_info
  type(serial_fft3d_info) :: internal_tb_fine_info
  ! qoh: Information arrays for FFT boxes
  type(serial_fft3d_info) :: internal_box_coarse_info
  type(serial_fft3d_info) :: internal_box_dbl_info
  ! ndmh: Information arrays for augmentation box on fine grid
  type(serial_fft3d_info) :: internal_aug_fine_info
  ! qoh: Information arrays for cell
  integer, parameter :: MAX_CELL_INFOS = 6
  type(parallel_fft3d_info) :: internal_cell_info(MAX_CELL_INFOS)
  integer :: cell_info_index = 0

  ! Workspace for interpolation and filtering
  complex(kind=DP), allocatable, dimension(:,:,:) :: coarse_work
  complex(kind=DP), allocatable, dimension(:,:,:) :: fine_work

#ifdef FFTW
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1
  integer, parameter :: FFTW_REAL_TO_COMPLEX = -1
  integer, parameter :: FFTW_COMPLEX_TO_REAL = 1
  integer, parameter :: FFTW_ESTIMATE = 0
  integer, parameter :: FFTW_MEASURE = 1
  integer, parameter :: FFTW_OUT_OF_PLACE = 0
  integer, parameter :: FFTW_IN_PLACE = 8
  integer, parameter :: FFTW_USE_WISDOM = 16
#endif

#ifdef FFTW3
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1
  integer, parameter :: FFTW_MEASURE = 0
  integer, parameter :: FFTW_ESTIMATE = 64
#endif

!CW
#ifdef FFTW3GPU
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1
  integer,parameter  :: batch=1
  integer, parameter :: FFTW_MEASURE = 0
  integer, parameter :: FFTW_ESTIMATE = 64
#endif
!END CW

#ifdef ACML
  integer, parameter :: ACML_MODE_PLAN = 100
  integer, parameter :: ACML_MODE_FORWARD = 1
  integer, parameter :: ACML_MODE_BACKWARD = -1
#endif

#ifdef VENDOR
#ifdef ALPHA
  integer, parameter :: dxml_structure_size = 256
  integer, external :: zfft_init_3d,zfft_apply_3d,zfft_exit_3d
  integer, external :: zfft_init_2d,zfft_apply_2d,zfft_exit_2d
  integer, external :: dfft_init,dfft_apply,dfft_exit
#endif
#endif

contains

  ! ***************************************
  ! ***  P u b l i c   r o u t i n e s  ***
  ! ***************************************

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_size(n,flag)

    implicit none

    ! Arguments
    integer, intent(inout) :: n           ! The FFT grid size
    integer, intent(in), optional :: flag ! Flag specifying whether result
    !                                       must be even or odd

    ! The following (FFT library dependent) parameters determine which
    ! factors are allowed:
    !   num_factors : the number of factors available
    !   factor(:)   : the factors themselves (in ascending order)
    !   limit(:)    : the limits on the power of the factors (-ve -> unlimited)

#ifdef FFTW
    integer, parameter :: num_factors = 6
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,1,1 /)
#endif
#ifdef FFTW3
    integer, parameter :: num_factors = 6
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,1,1 /)
#endif
!CW
#ifdef FFTW3GPU
    integer, parameter :: num_factors = 6
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,1,1 /)
#endif
!END CW
#ifdef ACML
    integer, parameter :: num_factors = 6
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,1,1 /)
#endif
#ifdef VENDOR
#ifdef ALPHA
    integer, parameter :: num_factors = 5
    integer, parameter :: factor(num_factors) = (/ 2,3,5,7,11 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,1,1 /)
#endif
#ifdef SUN
    integer, parameter :: num_factors = 7
    integer, parameter :: factor(num_factors) = (/ 2,3,4,5,7,11,13 /)
    integer, parameter :: limit(num_factors) = (/ -1,-1,-1,-1,-1,-1,-1 /)
#endif
#endif

    ! Local variables
    integer :: ifac
    integer :: fac
    integer :: power
    integer :: ntest
    integer :: stride

    ! Set starting point
    if (present(flag)) then
       if (mod(flag,2) == 0) then ! force n to be even
          n = n+mod(n,2)
       else                       ! force n to be odd
          n = n-mod(n,2)+1
       end if
    end if

    ! Set stride
    if (present(flag)) then
       stride = 2
    else
       stride = 1
    end if

    ! Check whether n is an allowed product of factors
    do
       ntest = n
       do ifac=num_factors,1,-1 ! take out largest factors first
          fac = factor(ifac)
          power = limit(ifac)
          do
             if (mod(ntest,fac) /= 0 .or. power == 0) exit
             ntest = ntest / fac
             power = power - 1
          end do
       end do
       if (ntest == 1) exit
       n = n + stride
    end do

  end subroutine fourier_size

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_init

    use rundat, only: pub_tightbox_fft_coarse, pub_tightbox_fft_fine, &
         pub_aug
    use simulation_cell, only: pub_fftbox, &
         pub_maxtight_pts1, pub_maxtight_pts2, pub_maxtight_pts3, &
         pub_aug_box_n1, pub_aug_box_n2, pub_aug_box_n3
    use utils, only : utils_alloc_check

    implicit none

    ! Local variables

    integer :: ierr ! Status flag

    ! qoh: tightbox on coarse grid
    if (pub_tightbox_fft_coarse) &
         call internal_serial_init(pub_maxtight_pts1, pub_maxtight_pts2,&
         pub_maxtight_pts3, pub_maxtight_pts1,pub_maxtight_pts2, &
         internal_tb_coarse_info)

    ! qoh: tightbox on double grid
    if (pub_tightbox_fft_fine) &
         call internal_serial_init(2*pub_maxtight_pts1, 2*pub_maxtight_pts2,&
         2*pub_maxtight_pts3, 2*pub_maxtight_pts1, 2*pub_maxtight_pts2, &
         internal_tb_fine_info)

    ! ndmh: augmentation box on fine grid
    if (pub_aug) &
         call internal_serial_init(pub_aug_box_n1, pub_aug_box_n2, &
         pub_aug_box_n3, pub_aug_box_n1, pub_aug_box_n2, internal_aug_fine_info)

    ! Coarse-Grid FFTbox
    call internal_serial_init(pub_fftbox%total_pt1,pub_fftbox%total_pt2, &
         pub_fftbox%total_pt3,pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         internal_box_coarse_info)

    ! Double-Grid FFTbox
    call internal_serial_init(pub_fftbox%total_pt1_dbl, &
         pub_fftbox%total_pt2_dbl,pub_fftbox%total_pt3_dbl, &
         pub_fftbox%total_ld1_dbl,pub_fftbox%total_ld2_dbl, &
         internal_box_dbl_info)

    ! Workspaces for FFT box filtering and interpolation
    ! Coarse Grid
    allocate(coarse_work(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('fourier_init','coarse_work',ierr)
    ! Double Grid
    allocate(fine_work(pub_fftbox%total_ld1_dbl,pub_fftbox%total_ld2_dbl, &
         pub_fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('fourier_init','fine_work',ierr)

  end subroutine fourier_init

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_init_cell(grid)

    use cell_grid, only: GRID_INFO
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid

    if (grid%fft_index/=0) then
       call utils_abort('Error in fourier_init: grid already initialised')
    end if

    cell_info_index = cell_info_index + 1

    if (cell_info_index>MAX_CELL_INFOS) then
       call utils_abort('Error in fourier_init: too many grids allocated')
    end if

    grid%fft_index = cell_info_index
    call internal_parallel_init(grid,internal_cell_info(grid%fft_index))

  end subroutine fourier_init_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_apply_cell_forward(rspc,gspc,grid)
    !=========================================================================!
    ! Performs a Fast Fourier Transform from real to reciprocal space.        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rspc (input):  Slab-distributed array that contains the real space    !
    !                  data to be transformed.                                !
    !   gspc (output): Slab-distributed array that will contain the result    !
    !                  of the transform in reciprocal space.                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 7/7/03                                         !
    ! Parallel version by Peter Haynes, 21/6/04                               !
    ! Improved for const-correctness and documented, Jacek Dziedzic, 05/2010  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only : comms_free
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: rspc(grid%ld1,grid%ld2,grid%max_slabs12)
    complex(kind=DP), intent(out) :: gspc(grid%ld3,grid%ld2,grid%max_slabs23)

    ! -----------------------------------------------------------------------

    call timer_clock('fourier_apply_cell_forward',1)

    ! Make a license heartbeat call, if ACCELRYS defined
    call make_license_heartbeat

    ! Do the hard work
    call internal_fft3d_cell_forward(internal_cell_info(grid%fft_index),rspc, &
         gspc)

    call timer_clock('fourier_apply_cell_forward',2, &
         calculate_flops(internal_cell_info(grid%fft_index)))

    ! Free up comms resources
    call comms_free

  end subroutine fourier_apply_cell_forward

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_apply_cell_backward(rspc,gspc,grid)
    !=========================================================================!
    ! Performs a Fast Fourier Transform from reciprocal to real space.        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rspc (output): Slab-distributed array that contains the reciprocal    !
    !                  space data to be transformed.                          !
    !   gspc (input):  Slab-distributed array that will contain the result    !
    !                  of the transform in real space.                        !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 7/7/03                                         !
    ! Parallel version by Peter Haynes, 21/6/04                               !
    ! Improved for const-correctness and documented, Jacek Dziedzic, 05/2010  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only : comms_free
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out) :: rspc(grid%ld1,grid%ld2,grid%max_slabs12)
    complex(kind=DP), intent(in) :: gspc(grid%ld3,grid%ld2,grid%max_slabs23)

    ! -----------------------------------------------------------------------

    call timer_clock('fourier_apply_cell_backward',1)

    ! Make a license heartbeat call, if ACCELRYS defined
    call make_license_heartbeat

    ! Do the hard work
    call internal_fft3d_cell_backward(internal_cell_info(grid%fft_index), &
         rspc,gspc)

    call timer_clock('fourier_apply_cell_backward',2, &
         calculate_flops(internal_cell_info(grid%fft_index)))

    ! Free up comms resources
    call comms_free

  end subroutine fourier_apply_cell_backward

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_interpolate_cell(rspc1,rspc2,grid1,grid2,apply_nyquist)
    !=========================================================================!
    ! Performs a Fast Fourier Transform Interpolation from a specified grid   !
    ! to a larger grid.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rspc1 (input):  Slab-distributed array that contains the real space   !
    !                   data on grid 1 to be interpolated to grid 2.          !
    !   rspc2 (output): Slab-distributed array that will contain the result   !
    !                   of the filtering in real space on grid 2.             !
    !   grid1 (input):  GRID_INFO describing grid 1 (the smaller grid).       !
    !   grid2 (input):  GRID_INFO describing grid 2 (the larger grid).        !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 29/06/10                                      !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only : comms_free, comms_recv, comms_send, pub_my_node_id
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid1
    type(GRID_INFO), intent(in) :: grid2
    real(kind=DP), intent(in) :: rspc1(grid1%ld1,grid1%ld2,grid1%max_slabs12)
    real(kind=DP), intent(out) :: rspc2(grid2%ld1,grid2%ld2,grid2%max_slabs12)
    logical, intent(in), optional :: apply_nyquist

    ! Local Variables
    integer :: ierr
    integer :: islab23,jslab23,kslab23
    integer :: n3,n2,m3,m2,e3,e2,l3,l2,k3,k2
    integer :: j1
    integer :: sendnode, recvnode
    complex(kind=DP), allocatable :: gspc1(:,:,:)
    complex(kind=DP), allocatable :: gspc2(:,:,:)
    complex(kind=DP), allocatable :: gbuf1(:,:)
    logical :: loc_apply_nyquist


    ! -----------------------------------------------------------------------

    call timer_clock('fourier_interpolate_cell',1)

    ! Optional arguments
    loc_apply_nyquist = .false.
    if (present(apply_nyquist)) loc_apply_nyquist = apply_nyquist

    ! Allocate workspaces
    allocate(gspc1(grid1%ld3,grid1%ld2,grid1%max_slabs23),stat=ierr)
    call utils_alloc_check('fourier_interpolate_cell','gspc1',ierr)
    allocate(gspc2(grid2%ld3,grid2%ld2,grid2%max_slabs23),stat=ierr)
    call utils_alloc_check('fourier_interpolate_cell','gspc2',ierr)
    allocate(gbuf1(grid1%ld3,grid1%ld2),stat=ierr)
    call utils_alloc_check('fourier_interpolate_cell','gbuf1',ierr)

    ! Make a license heartbeat call, if ACCELRYS defined
    call make_license_heartbeat

    ! Transform from real to reciprocal space on grid 1
    call internal_fft3d_cell_forward(internal_cell_info(grid1%fft_index), &
         rspc1,gspc1)

    ! Nyquist filter (grid is always going to be even)
    if (grid1%num_slabs23 > 0) then
       gspc1(grid1%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       gspc1(:,grid1%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if
    ! Nyquist filter for last slab
    if (pub_my_node_id==grid1%node_slab23(grid1%n1/2+1)) &
         gspc1(:,:,grid1%num_slabs23) = (0.0_DP,0.0_DP)

    ! Scale grid1 recip data
    gspc1 = gspc1 * (real(grid2%n1,kind=DP) * real(grid2%n2,kind=DP) &
         * real(grid2%n3,kind=DP)) / (real(grid1%n1,kind=DP) &
         * real(grid1%n2,kind=DP) * real(grid1%n3,kind=DP))

    ! Find limits in reciprocal i3,i2 directions on grids 1 and 2
    n3 = grid1%n3;         n2 = grid1%n2
    e3 = grid2%n3;         e2 = grid2%n2
    m3 = (n3+1)/2;         m2 = (n2+1)/2
    l3 = n3/2+1+mod(n3,2); l2 = n2/2+1+mod(n2,2)
    k3 = e3-n3+l3;         k2 = e2-n2+l2
    gspc2(:,:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)

    ! Loop over all slabs of grid2 (the larger grid)
    do islab23=1,grid2%max_slabs23

       ! Find if there are any slabs needing to be sent to other nodes
       ! on this step
       do jslab23=1,grid1%num_slabs23
          ! Find the global index of this slab, the node holding this slab
          ! in grid2, and the local index of this slab on that node
          j1 = jslab23 + grid1%first_slab23(pub_my_node_id) - 1
          sendnode = grid2%node_slab23(j1)
          kslab23 = j1 - grid2%first_slab23(sendnode) + 1
          ! If the local index on that node is the same as the step we
          ! are currently on, then send the slab
          if (kslab23==islab23) then
             if (sendnode/=pub_my_node_id) then
                call comms_send(sendnode,gspc1(:,:,jslab23),tag=j1)
             end if
          end if
       end do

       ! Find index of this slab in global array
       j1 = islab23 + grid2%first_slab23(pub_my_node_id) - 1
       ! As long as there is a slab left on this node
       if (islab23<=grid2%num_slabs23) then
          ! And as long as we are still within the slabs of grid1
          if (j1<=(grid1%n1/2+1)) then
             ! Find which node holds the slab required on this node on this step
             recvnode = grid1%node_slab23(j1)

             ! Receive slab, or copy if local
             if (recvnode/=pub_my_node_id) then
                call comms_recv(recvnode,gbuf1(:,:),tag=j1)
             else
                jslab23 = j1 - grid1%first_slab23(pub_my_node_id) + 1
                gbuf1(:,:) = gspc1(:,:,jslab23)
             end if

             gspc2(1:m3,1:m2,islab23)   = gbuf1(1:m3,1:m2)
             gspc2(k3:e3,1:m2,islab23)  = gbuf1(l3:n3,1:m2)
             gspc2(1:m3,k2:e2,islab23)  = gbuf1(1:m3,l2:n2)
             gspc2(k3:e3,k2:e2,islab23) = gbuf1(l3:n3,l2:n2)
          end if
       end if

    end do

    ! Transform from reciprocal to real space on grid 2
    call internal_fft3d_cell_backward(internal_cell_info(grid2%fft_index), &
         rspc2,gspc2)

    call timer_clock('fourier_interpolate_cell',2, &
         calculate_flops(internal_cell_info(grid1%fft_index)) &
         + calculate_flops(internal_cell_info(grid2%fft_index)))

    ! Free up comms resources
    call comms_free

    ! Dellocate workspaces
    deallocate(gbuf1,stat=ierr)
    call utils_dealloc_check('fourier_interpolate_cell','gbuf1',ierr)
    deallocate(gspc2,stat=ierr)
    call utils_dealloc_check('fourier_interpolate_cell','gspc2',ierr)
    deallocate(gspc1,stat=ierr)
    call utils_dealloc_check('fourier_interpolate_cell','gspc1',ierr)

  end subroutine fourier_interpolate_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_filter_cell(rspc1,rspc2,grid1,grid2,apply_nyquist)
    !=========================================================================!
    ! Performs a Fast Fourier Transform Filtering from a specified grid to a  !
    ! smaller grid.                                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rspc1 (input):  Slab-distributed array that contains the real space   !
    !                   data on grid 1 to be filtered to grid 2.              !
    !   rspc2 (output): Slab-distributed array that will contain the result   !
    !                   of the filtering in real space on grid 2.             !
    !   grid1 (input):  GRID_INFO describing grid 1 (the larger grid).        !
    !   grid2 (input):  GRID_INFO describing grid 2 (the smaller grid).       !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 29/06/10                                      !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only : comms_free, comms_recv, comms_send, pub_my_node_id
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid1
    type(GRID_INFO), intent(in) :: grid2
    real(kind=DP), intent(in) :: rspc1(grid1%ld1,grid1%ld2,grid1%max_slabs12)
    real(kind=DP), intent(out) :: rspc2(grid2%ld1,grid2%ld2,grid2%max_slabs12)
    logical, intent(in), optional :: apply_nyquist

    ! Local Variables
    integer :: ierr
    integer :: islab23,jslab23,kslab23
    integer :: n3,n2,m3,m2,e3,e2,l3,l2,k3,k2
    integer :: j1
    integer :: sendnode, recvnode
    complex(kind=DP), allocatable :: gspc1(:,:,:)
    complex(kind=DP), allocatable :: gspc2(:,:,:)
    complex(kind=DP), allocatable :: gbuf1(:,:)
    logical :: loc_apply_nyquist

    ! -----------------------------------------------------------------------

    call timer_clock('fourier_filter_cell',1)

    ! Optional arguments
    loc_apply_nyquist = .false.
    if (present(apply_nyquist)) loc_apply_nyquist = apply_nyquist

    ! Allocate workspaces
    allocate(gspc1(grid1%ld3,grid1%ld2,grid1%max_slabs23),stat=ierr)
    call utils_alloc_check('fourier_filter_cell','gspc1',ierr)
    allocate(gspc2(grid2%ld3,grid2%ld2,grid2%max_slabs23),stat=ierr)
    call utils_alloc_check('fourier_filter_cell','gspc2',ierr)
    allocate(gbuf1(grid1%ld3,grid1%ld2),stat=ierr)
    call utils_alloc_check('fourier_filter_cell','gbuf1',ierr)

    ! Make a license heartbeat call, if ACCELRYS defined
    call make_license_heartbeat

    ! Transform from real to reciprocal space on grid 1 (the larger grid)
    call internal_fft3d_cell_forward(internal_cell_info(grid1%fft_index), &
         rspc1,gspc1)

    ! Scale grid1 recip data
    gspc1 = gspc1 * (real(grid2%n1,kind=DP) * real(grid2%n2,kind=DP) &
         * real(grid2%n3,kind=DP)) / (real(grid1%n1,kind=DP) * &
         real(grid1%n2,kind=DP) * real(grid1%n3,kind=DP))

    ! Find limits in reciprocal i3,i2 directions on grids 1 and 2
    n3 = grid1%n3;         n2 = grid1%n2
    e3 = grid2%n3;         e2 = grid2%n2
    m3 = (e3+1)/2;         m2 = (e2+1)/2
    k3 = e3/2+1+mod(e3,2); k2 = e2/2+1+mod(e2,2)
    l3 = n3-e3+k3;         l2 = n2-e2+k2
    gspc2(:,:,:) = cmplx(0.0_DP,0.0_DP,kind=DP)

    ! Loop over all slabs of grid2 (the smaller grid)
    do islab23=1,grid2%max_slabs23

       ! Find if there are any slabs needing to be sent to other nodes
       ! on this step
       do jslab23=1,grid1%num_slabs23
          ! Find the global index of this slab, the node holding this slab
          ! in grid2, and the local index of this slab on that node
          j1 = jslab23 + grid1%first_slab23(pub_my_node_id) - 1
          if (j1>(grid2%n1/2+1)) cycle
          sendnode = grid2%node_slab23(j1)
          kslab23 = j1 - grid2%first_slab23(sendnode) + 1
          ! If the local index on that node is the same as the step we
          ! are currently on, then send the slab
          if (kslab23==islab23) then
             if (sendnode/=pub_my_node_id) then
                call comms_send(sendnode,gspc1(1,1,jslab23),grid1%ld3*grid1%ld2,tag=j1)
             end if
          end if
       end do

       ! Find which node holds the slab required for this node on this step
       j1 = islab23 + grid2%first_slab23(pub_my_node_id) - 1
       if (islab23<=grid2%num_slabs23) then
          recvnode = grid1%node_slab23(j1)

          jslab23 = j1 - grid1%first_slab23(recvnode) + 1
          if (recvnode/=pub_my_node_id) then
             call comms_recv(recvnode,gbuf1(1,1),grid1%ld3*grid1%ld2,tag=j1)
          else
             gbuf1(:,:) = gspc1(:,:,jslab23)
          end if

          gspc2(1:m3,1:m2,islab23)   = gbuf1(1:m3,1:m2)
          gspc2(k3:e3,1:m2,islab23)  = gbuf1(l3:n3,1:m2)
          gspc2(1:m3,k2:e2,islab23)  = gbuf1(1:m3,l2:n2)
          gspc2(k3:e3,k2:e2,islab23) = gbuf1(l3:n3,l2:n2)

       end if

    end do

    ! Nyquist filter (grid is always going to be even)
    if (grid2%num_slabs23 > 0) then
       gspc2(grid2%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       gspc2(:,grid2%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if
    ! Nyquist filter for last slab
    if (pub_my_node_id==grid2%node_slab23(grid2%n1/2+1)) &
         gspc2(:,:,grid2%num_slabs23) = (0.0_DP,0.0_DP)

    ! Transform from reciprocal to real space on grid 2
    call internal_fft3d_cell_backward(internal_cell_info(grid2%fft_index), &
         rspc2,gspc2)

    call timer_clock('fourier_filter_cell',2, &
         calculate_flops(internal_cell_info(grid1%fft_index)) &
         + calculate_flops(internal_cell_info(grid2%fft_index)))

    ! Free up comms resources
    call comms_free

    ! Dellocate workspaces
    deallocate(gbuf1,stat=ierr)
    call utils_dealloc_check('fourier_filter_cell','gbuf1',ierr)
    deallocate(gspc2,stat=ierr)
    call utils_dealloc_check('fourier_filter_cell','gspc2',ierr)
    deallocate(gspc1,stat=ierr)
    call utils_dealloc_check('fourier_filter_cell','gspc1',ierr)

  end subroutine fourier_filter_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_apply_box(grid,dir,rspc,gspc,scale_in,tb,aug,padbox)

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    use timer, only : timer_clock

    implicit none

    character, intent(in) :: grid
    character, intent(in) :: dir
    complex(kind=DP), dimension(:,:,:), intent(inout) :: rspc
    complex(kind=DP), dimension(:,:,:), intent(inout), optional :: gspc
    logical, optional, intent(in) :: scale_in

    ! Switches for specific types of transform (tightbox, augbox, padded box)
    logical, optional, intent(in) :: tb
    logical, optional, intent(in) :: aug
    logical, optional, intent(in) :: padbox

    ! Local variables
    integer :: npts ! Total number of points in FFT grid
    real(kind=DP) :: flops
    logical :: scale
    logical :: tightbox ! Whether we are doing an FFT on a tightbox
    logical :: augbox ! Whether we are doing an FFT on the augmentation box

    ! Check arguments

    ! Optional supression of scaling for backwards transform
    if (present(scale_in)) then
       scale = scale_in
    else
       scale = .true.
    end if

    ! qoh: Optionally do FFT on tightbox instead of FFT box
    tightbox = .false.
    if (present(tb)) tightbox = tb

    ! ndmh: Optionally do FFT on augmentation box
    augbox = .false.
    if (present(aug)) augbox = aug

    if (grid /= 'C' .and. grid /= 'c' .and. grid /= 'F' .and. grid /= 'f') then
       if (pub_on_root) write(stdout,'(3a)') 'Error in fourier_apply_box: &
            &unknown grid flag "',grid,'"'
       call comms_abort
    end if

    if (dir /= 'F' .and. dir /= 'f' .and. dir /= 'B' .and. dir /= 'b') then
       if (pub_on_root) write(stdout,'(3a)') 'Error in fourier_apply_box: &
            &unknown direction flag "',dir,'"'
       call comms_abort
    end if

    if (present(gspc)) then

       if (size(rspc,1) /= size(gspc,1)) then
          if (pub_on_root) write(stdout,'(a)') 'Error in fourier_apply_box: &
               &array size mismatch in dimension 1'
          call comms_abort
       end if

       if (size(rspc,2) /= size(gspc,2)) then
          if (pub_on_root) write(stdout,'(a)') 'Error in fourier_apply_box: &
               &array size mismatch in dimension 2'
          call comms_abort
       end if

       if (size(rspc,3) /= size(gspc,3)) then
          if (pub_on_root) write(stdout,'(a)') 'Error in fourier_apply_box: &
               &array size mismatch in dimension 3'
          call comms_abort
       end if

    end if

    ! Start timer
    call timer_clock('fourier_apply_box',1)

    ! Copy input data into output array for in-place transform

    if (present(gspc)) then
       if (dir == 'F' .or. dir == 'f') then
          gspc = rspc
       else
          rspc = gspc
       end if
    end if

    ! Select appropriate routine
    if (tightbox) then
       if (grid == 'C' .or. grid == 'c') then
          if (present(gspc) .and. (dir == 'F' .or. dir == 'f')) then
             call internal_fft3d_box(internal_tb_coarse_info,dir,gspc,scale)
          else
             call internal_fft3d_box(internal_tb_coarse_info,dir,rspc,scale)
          end if
          npts = internal_tb_coarse_info%n1 * internal_tb_coarse_info%n2 * &
               internal_tb_coarse_info%n3
       else
          if (present(gspc) .and. (dir == 'F' .or. dir == 'f')) then
             call internal_fft3d_box(internal_tb_fine_info,dir,gspc,scale)
          else;
             call internal_fft3d_box(internal_tb_fine_info,dir,rspc,scale)
          end if
          npts = internal_tb_fine_info%n1 * internal_tb_fine_info%n2 * &
               internal_tb_fine_info%n3
       end if
    else if (augbox) then
       if (grid /= 'F' .and. grid /= 'f') then
          if (pub_on_root) write(stdout,'(3a)') 'Error in fourier_apply_box: &
               &can only use grid flag "F" with aug box transform'
          call comms_abort
          npts = 0 ! supress warning
       else
          if (present(gspc) .and. (dir == 'F' .or. dir == 'f')) then
             call internal_fft3d_box(internal_aug_fine_info,dir,gspc,scale)
          else;
             call internal_fft3d_box(internal_aug_fine_info,dir,rspc,scale)
          end if
          npts = internal_aug_fine_info%n1 * internal_aug_fine_info%n2 * &
               internal_aug_fine_info%n3
       end if
    else
       ! qoh: FFT box
       if (grid == 'C' .or. grid == 'c') then
          if (present(gspc) .and. (dir == 'F' .or. dir == 'f')) then
             call internal_fft3d_box(internal_box_coarse_info,dir,gspc,scale, &
                  padbox)
          else
             call internal_fft3d_box(internal_box_coarse_info,dir,rspc,scale, &
                  padbox)
          end if
          npts = internal_box_coarse_info%n1 * internal_box_coarse_info%n2 * &
               internal_box_coarse_info%n3
       else
          if (present(gspc) .and. (dir == 'F' .or. dir == 'f')) then
             call internal_fft3d_box(internal_box_dbl_info,dir,gspc,scale, &
                  padbox)
          else
             call internal_fft3d_box(internal_box_dbl_info,dir,rspc,scale, &
                  padbox)
          end if
          npts = internal_box_dbl_info%n1 * internal_box_dbl_info%n2 * &
               internal_box_dbl_info%n3
       end if
    end if


    ! Stop timer
    flops = 5.0_DP * npts * log(real(npts,kind=DP)) / log(2.0_DP)
    call timer_clock('fourier_apply_box',2,flops)

  end subroutine fourier_apply_box

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_apply_box_pair(grid,dir,rspc1,rspc2,gspc)

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    use timer, only : timer_clock

    implicit none

    character, intent(in) :: grid
    character, intent(in) :: dir
    real(kind=DP), dimension(:,:,:), intent(inout) :: rspc1, rspc2
    complex(kind=DP), dimension(:,:,:), intent(inout) :: gspc

    ! Local variables

    integer :: npts ! Total number of points in FFT grid
    real(kind=DP) :: flops

    ! Check arguments

    if (grid /= 'C' .and. grid /= 'c' .and. grid /= 'F' .and. grid /= 'f') then
       if (pub_on_root) write(stdout,'(3a)') 'Error in fourier_apply_box_pair: &
            &unknown grid flag "',grid,'"'
       call comms_abort
    end if

    if (dir /= 'F' .and. dir /= 'f' .and. dir /= 'B' .and. dir /= 'b') then
       if (pub_on_root) write(stdout,'(3a)') 'Error in fourier_apply_box_pair: &
            &unknown direction flag "',dir,'"'
       call comms_abort
    end if

    if (size(rspc1,1) /= size(gspc,1) .or. size(rspc2,1) /= size(gspc,1)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_apply_box_pair: &
            &array size mismatch in dimension 1'
       call comms_abort
    end if

    if (size(rspc1,2) /= size(gspc,2) .or. size(rspc2,2) /= size(gspc,2)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_apply_box_pair: &
            &array size mismatch in dimension 2'
       call comms_abort
    end if

    if (size(rspc1,3) /= size(gspc,3) .or. size(rspc2,3) /= size(gspc,3)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_apply_box_pair: &
            &array size mismatch in dimension 3'
       call comms_abort
    end if

    ! Start timer
    call timer_clock('fourier_apply_box_pair',1)

    ! For forward routine, copy data into gspc array
    if (dir == 'F' .or. dir == 'f') gspc = cmplx(rspc1,rspc2,kind=DP)

    ! Select appropriate routine

    if (grid == 'C' .or. grid == 'c') then
       call internal_fft3d_box(internal_box_coarse_info,dir,gspc,.true.)
       npts = internal_box_coarse_info%n1 * internal_box_coarse_info%n2 * &
            internal_box_coarse_info%n3
    else
       call internal_fft3d_box(internal_box_dbl_info,dir,gspc,.true.)
       npts = internal_box_dbl_info%n1 * internal_box_dbl_info%n2 * &
            internal_box_dbl_info%n3
    end if

    ! For backward routine, copy data out of gspc array
    if (dir == 'B' .or. dir == 'b') then
       rspc1 = real(gspc,kind=DP)
       rspc2 = aimag(gspc)
    end if

    ! Stop timer
    flops = 5.0_DP * npts * log(real(npts,kind=DP)) / log(2.0_DP)
    call timer_clock('fourier_apply_box_pair',2,flops)

  end subroutine fourier_apply_box_pair

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_exit

    use rundat, only: pub_tightbox_fft_coarse, pub_tightbox_fft_fine, pub_aug
    use utils, only : utils_dealloc_check

    implicit none

    ! Local variables
    integer :: ierr ! Status flag

    ! Workspace for filtering and interpolation
    deallocate(fine_work,stat=ierr)
    call utils_dealloc_check('fourier_exit','fine_work',ierr)
    deallocate(coarse_work,stat=ierr)
    call utils_dealloc_check('fourier_exit','coarse_work',ierr)

    ! Double FFTbox grid
    call internal_serial_exit(internal_box_dbl_info)

    ! Coarse FFTbox grid
    call internal_serial_exit(internal_box_coarse_info)

    ! ndmh: augmentation box on fine grid
    if (pub_aug) &
         call internal_serial_exit(internal_aug_fine_info)

    ! qoh: tightbox on fine grid
    if (pub_tightbox_fft_fine) &
         call internal_serial_exit(internal_tb_fine_info)

    ! qoh: tightbox on coarse grid
    if (pub_tightbox_fft_coarse) &
         call internal_serial_exit(internal_tb_coarse_info)

  end subroutine fourier_exit

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_exit_cell

    implicit none

    ! Local variables
    integer :: igrid ! whole-cell grid counter

    ! Whole-simulation-cell grids
    do igrid=cell_info_index,1,-1
       call internal_parallel_exit(internal_cell_info(igrid))
    end do
    cell_info_index = 0

#ifdef FFTW3
    call dfftw_cleanup
#endif

  end subroutine fourier_exit_cell

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_interpolate(in1,in2,out1,out2)

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out1(:,:,:)
    real(kind=DP), intent(out) :: out2(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3
    real(kind=DP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef DEBUG
    ! Check array sizes
    if (size(in1,1) < n1 .or. size(in1,2) < n2 .or. size(in1,3) < n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array in1'
       call comms_abort
    end if
    if (size(in2,1) < n1 .or. size(in2,2) < n2 .or. size(in2,3) < n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array in2'
       call comms_abort
    end if
    if (size(out1,1) < 2*n1 .or. size(out1,2) < 2*n2 .or. &
         size(out1,3) < 2*n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array out1'
       call comms_abort
    end if
    if (size(out2,1) < 2*n1 .or. size(out2,2) < 2*n2 .or. &
         size(out2,3) < 2*n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array out2'
       call comms_abort
    end if
#endif

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scale = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)
    coarse_work(1:n1,1:n2,1:n3) = scale * cmplx(in1(1:n1,1:n2,1:n3), &
         in2(1:n1,1:n2,1:n3),kind=DP)

    ! Fourier transform to reciprocal space on coarse grid

    call fourier_apply_box('Coarse','Forward',coarse_work,scale_in=.false.)

    ! Copy Fourier components to fine reciprocal space grid

    fine_work(1:m1,1:m2,1:m3) = coarse_work(1:m1,1:m2,1:m3)
    fine_work(k1:2*n1,1:m2,1:m3) = coarse_work(l1:n1,1:m2,1:m3)
    fine_work(1:m1,k2:2*n2,1:m3) = coarse_work(1:m1,l2:n2,1:m3)
    fine_work(k1:2*n1,k2:2*n2,1:m3) = coarse_work(l1:n1,l2:n2,1:m3)
    fine_work(1:m1,1:m2,k3:2*n3) = coarse_work(1:m1,1:m2,l3:n3)
    fine_work(k1:2*n1,1:m2,k3:2*n3) = coarse_work(l1:n1,1:m2,l3:n3)
    fine_work(1:m1,k2:2*n2,k3:2*n3) = coarse_work(1:m1,l2:n2,l3:n3)
    fine_work(k1:2*n1,k2:2*n2,k3:2*n3) = coarse_work(l1:n1,l2:n2,l3:n3)

    ! Zero all other components

    fine_work(:,:,l3:j3) = 0.0_DP
    fine_work(:,l2:j2,1:l3-1) = 0.0_DP
    fine_work(:,l2:j2,j3+1:) = 0.0_DP
    fine_work(l1:j1,1:l2-1,1:l3-1) = 0.0_DP
    fine_work(l1:j1,j2+1:,1:l3-1) = 0.0_DP
    fine_work(l1:j1,1:l2-1,j3+1:) = 0.0_DP
    fine_work(l1:j1,j2+1:,j3+1:) = 0.0_DP

    ! Fourier transform to real space on fine grid

    call fourier_apply_box('Fine','Backward',fine_work,scale_in=.false., &
         padbox=.true.)

    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
    do j3 = 1,2*n3
       do j2 = 1,2*n2
          out1(1:2*n1,j2,j3) = real(fine_work(1:2*n1,j2,j3),DP)
          out2(1:2*n1,j2,j3) = aimag(fine_work(1:2*n1,j2,j3))
       end do
    end do

  end subroutine fourier_interpolate

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_interpolate_product(in1,in2,out)

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3,j1,j2,j3
    real(kind=DP) :: scale


    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef DEBUG
    ! Check array sizes
    if (size(in1,1) < n1 .or. size(in1,2) < n2 .or. size(in1,3) < n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array in1'
       call comms_abort
    end if
    if (size(in2,1) < n1 .or. size(in2,2) < n2 .or. size(in2,3) < n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array in2'
       call comms_abort
    end if
    if (size(out,1) < 2*n1 .or. size(out,2) < 2*n2 .or. &
         size(out,3) < 2*n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array out'
       call comms_abort
    end if
#endif

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1+l1 ; k2 = n2+l2 ; k3 = n3+l3
    j1 = k1-2+mod(n1,2) ; j2 = k2-2+mod(n2,2) ; j3 = k3-2+mod(n3,2)

    ! Pack data into real and imaginary parts of complex array
    ! ndmh: and normalise in advance as it uses fewer operations
    scale = 8.0_DP / (internal_box_dbl_info%ld1*internal_box_dbl_info%ld2* &
         internal_box_dbl_info%n3)
    coarse_work(1:n1,1:n2,1:n3) = scale * cmplx(in1(1:n1,1:n2,1:n3), &
         in2(1:n1,1:n2,1:n3),kind=DP)

    ! Fourier transform to reciprocal space on coarse grid

    call fourier_apply_box('Coarse','Forward',coarse_work,scale_in=.false.)

    ! Copy Fourier components to fine reciprocal space grid

    fine_work(1:m1,1:m2,1:m3) = coarse_work(1:m1,1:m2,1:m3)
    fine_work(k1:2*n1,1:m2,1:m3) = coarse_work(l1:n1,1:m2,1:m3)
    fine_work(1:m1,k2:2*n2,1:m3) = coarse_work(1:m1,l2:n2,1:m3)
    fine_work(k1:2*n1,k2:2*n2,1:m3) = coarse_work(l1:n1,l2:n2,1:m3)
    fine_work(1:m1,1:m2,k3:2*n3) = coarse_work(1:m1,1:m2,l3:n3)
    fine_work(k1:2*n1,1:m2,k3:2*n3) = coarse_work(l1:n1,1:m2,l3:n3)
    fine_work(1:m1,k2:2*n2,k3:2*n3) = coarse_work(1:m1,l2:n2,l3:n3)
    fine_work(k1:2*n1,k2:2*n2,k3:2*n3) = coarse_work(l1:n1,l2:n2,l3:n3)

    ! Zero all other components

    fine_work(:,:,l3:j3) = 0.0_DP
    fine_work(:,l2:j2,1:l3-1) = 0.0_DP
    fine_work(:,l2:j2,j3+1:) = 0.0_DP
    fine_work(l1:j1,1:l2-1,1:l3-1) = 0.0_DP
    fine_work(l1:j1,j2+1:,1:l3-1) = 0.0_DP
    fine_work(l1:j1,1:l2-1,j3+1:) = 0.0_DP
    fine_work(l1:j1,j2+1:,j3+1:) = 0.0_DP

    ! Fourier transform to real space on fine grid
    call fourier_apply_box('Fine','Backward',fine_work,scale_in=.false., &
         padbox=.true.)

    ! ndmh: unpack data from real and imaginary parts of complex array
    ! ndmh: and multiply them together to find product
    do j3 = 1,2*n3
       do j2 = 1,2*n2
          out(1:2*n1,j2,j3) = real(fine_work(1:2*n1,j2,j3),DP) * &
               aimag(fine_work(1:2*n1,j2,j3))
       end do
    end do

!CW
    call resample
    contains
    subroutine resample
    do j3 = 1,2*n3
      do j2 = 1,2*n2
          out(1:2*n1,j2,j3) = real(fine_work(1:2*n1,j2,j3),DP) * aimag(fine_work(1:2*n1,j2,j3))
      end do
    end do
    end subroutine
!END CW

  end subroutine fourier_interpolate_product

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_filter(in1,in2,out1,out2)

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: in1(:,:,:)
    real(kind=DP), intent(in) :: in2(:,:,:)
    real(kind=DP), intent(out) :: out1(:,:,:)
    real(kind=DP), intent(out) :: out2(:,:,:)

    ! Local variables
    integer :: n1,n2,n3
    integer :: m1,m2,m3,l1,l2,l3,k1,k2,k3
    real(kind=DP) :: scale

    ! Coarse FFT box dimensions
    n1 = internal_box_coarse_info%n1
    n2 = internal_box_coarse_info%n2
    n3 = internal_box_coarse_info%n3

#ifdef DEBUG
    ! Check array sizes
    if (size(in1,1) < 2*n1 .or. size(in1,2) < 2*n2 .or. &
         size(in1,3) < 2*n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array in1'
       call comms_abort
    end if
    if (size(in2,1) < 2*n1 .or. size(in2,2) < 2*n2 .or. &
         size(in2,3) < 2*n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array in2'
       call comms_abort
    end if
    if (size(out1,1) < n1 .or. size(out1,2) < n2 .or. size(out1,3) < n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array out1'
       call comms_abort
    end if
    if (size(out2,1) < n1 .or. size(out2,2) < n2 .or. size(out2,3) < n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in fourier_interpolate: &
            &invalid dimensions for array out2'
       call comms_abort
    end if
#endif

    m1 = (n1+1)/2 ; m2 = (n2+1)/2 ; m3 = (n3+1)/2
    l1 = n1/2+2 ; l2 = n2/2+2 ; l3 = n3/2+2
    k1 = n1 + l1 ; k2 = n2 + l2 ; k3 = n3 + l3

    ! Pack data into real and imaginary parts of complex array

    fine_work(1:2*n1,1:2*n2,1:2*n3) = cmplx(in1(1:2*n1,1:2*n2,1:2*n3), &
         in2(1:2*n1,1:2*n2,1:2*n3),dp)

    ! Fourier transform to reciprocal space on fine grid

    call fourier_apply_box('Fine','Forward',fine_work,scale_in=.false., &
         padbox=.true.)

    ! Copy Fourier components to coarse reciprocal space grid

    coarse_work(1:m1,1:m2,1:m3) = fine_work(1:m1,1:m2,1:m3)
    coarse_work(l1:n1,1:m2,1:m3) = fine_work(k1:2*n1,1:m2,1:m3)
    coarse_work(1:m1,l2:n2,1:m3) = fine_work(1:m1,k2:2*n2,1:m3)
    coarse_work(l1:n1,l2:n2,1:m3) = fine_work(k1:2*n1,k2:2*n2,1:m3)
    coarse_work(1:m1,1:m2,l3:n3) = fine_work(1:m1,1:m2,k3:2*n3)
    coarse_work(l1:n1,1:m2,l3:n3) = fine_work(k1:2*n1,1:m2,k3:2*n3)
    coarse_work(1:m1,l2:n2,l3:n3) = fine_work(1:m1,k2:2*n2,k3:2*n3)
    coarse_work(l1:n1,l2:n2,l3:n3) = fine_work(k1:2*n1,k2:2*n2,k3:2*n3)

    ! Fourier transform to real space on coarse grid

    call fourier_apply_box('Coarse','Backward',coarse_work,scale_in=.false.)

    ! Normalise
    scale = 0.125_DP / (internal_box_coarse_info%ld1*&
         internal_box_coarse_info%ld2*internal_box_coarse_info%n3)
    coarse_work(1:n1,1:n2,1:n3) = scale * coarse_work(1:n1,1:n2,1:n3)

    ! Unpack data from real and imaginary parts of complex array
    ! ndmh: a loop saves on cache misses when repeating for aimag
    do k3 = 1,n3
       do k2 = 1,n2
          do k1 = 1,n1
             out1(k1,k2,k3) = real(coarse_work(k1,k2,k3),dp)
             out2(k1,k2,k3) = aimag(coarse_work(k1,k2,k3))
          end do
       end do
    end do

    !out1(1:n1,1:n2,1:n3) = real(coarse_work(1:n1,1:n2,1:n3),dp)
    !out2(1:n1,1:n2,1:n3) = aimag(coarse_work(1:n1,1:n2,1:n3))


  end subroutine fourier_filter

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ! *****************************************
  ! ***  P r i v a t e   r o u t i n e s  ***
  ! *****************************************

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_serial_init(n1,n2,n3,ld1,ld2,fft3d_info)

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    use utils, only : utils_alloc_check

    implicit none

    integer, intent(in) :: n1,n2,n3 ! Dimensions of FFT grid
    integer, intent(in) :: ld1,ld2  ! Array dimensions containing grid
    type(serial_fft3d_info), intent(out) :: fft3d_info
!CW
    integer :: i
!END CW

    ! Local variables

    integer :: ierr                 ! Error flag
#ifdef SUN
    complex(kind=DP) :: zdum
#endif

    ! Simple checks on validity of input arguments

    if (n1 < 1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_serial_init &
            &(fourier_mod.F90): n1 < 1'
       call comms_abort
    end if

    if (n2 < 1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_serial_init &
            &(fourier_mod.F90): n2 < 1'
       call comms_abort
    end if

    if (n3 < 1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_serial_init &
            &(fourier_mod.F90): n3 < 1'
       call comms_abort
    end if

    if (ld1 < n1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_serial_init &
            &(fourier_mod.F90): ld1 < n1'
       call comms_abort
    end if

    if (ld2 < n2) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_serial_init &
            &(fourier_mod.F90): ld2 < n2'
       call comms_abort
    end if

    ! Platform-dependent checking

#ifdef FFTW
    if (ld1 /= n1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_serial_init &
            &(fourier_mod.F90): ld1 must equal n1 for FFTw'
       call comms_abort
    end if
    if (ld2 /= n2) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_serial_init &
            &(fourier_mod.F90): ld2 must equal n2 for FFTw'
       call comms_abort
    end if
#endif

    ! Simple initialisation

    ! Store basic information in fft3d_info structure

    fft3d_info%n1 = n1
    fft3d_info%n2 = n2
    fft3d_info%n3 = n3
    fft3d_info%ld1 = ld1
    fft3d_info%ld2 = ld2
    fft3d_info%iwork_len = 0
    fft3d_info%dwork_len = 0
    fft3d_info%zwork_len = 0

    ! Platform-dependent workspace allocation

#ifdef FFTW3
    fft3d_info%iwork_len = 12
    fft3d_info%zwork_len = ld1*ld2*n3
#endif

!CW
#ifdef FFTW3GPU
    fft3d_info%iwork_len = 6
    fft3d_info%zwork_len = ld1*ld2*n3
#endif
!END CW

#ifdef ACML
    fft3d_info%zwork_len = ld1*ld2*n3
    fft3d_info%tbl_size = ld1*ld2*n3 + 4*(ld1+ld2+n3) + 400
    allocate(fft3d_info%tbl(fft3d_info%tbl_size),stat=ierr)
    call utils_alloc_check('internal_serial_init (fourier_mod.F90)', &
         'fft3d_info%tbl',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    fft3d_info%iwork_len = dxml_structure_size
#endif
#ifdef SUN
    fft3d_info%iwork_len = 3*128
    fft3d_info%dwork_len = 2*(n1+n2+n3)
    fft3d_info%zwork_len = max(n1,n2,n3) + 16*n3
#endif
#endif

    if (fft3d_info%iwork_len > 0) then
       allocate(fft3d_info%iwork(fft3d_info%iwork_len),stat=ierr)
       call utils_alloc_check('internal_serial_init (fourier_mod.F90)', &
            'fft3d_info%iwork',ierr)
    end if

    if (fft3d_info%dwork_len > 0) then
       allocate(fft3d_info%dwork(fft3d_info%dwork_len),stat=ierr)
       call utils_alloc_check('internal_serial_init (fourier_mod.F90)', &
            'fft3d_info%dwork',ierr)
    end if

    if (fft3d_info%zwork_len > 0) then
       allocate(fft3d_info%zwork(fft3d_info%zwork_len),stat=ierr)
       call utils_alloc_check('internal_serial_init (fourier_mod.F90)', &
            'fft3d_info%zwork',ierr)
    end if

    ! Platform-dependent initialisation

#ifdef FFTW
    call fftw3d_f77_create_plan(fft3d_info%forward_plan,n1,n2,n3, &
         FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
    call fftw3d_f77_create_plan(fft3d_info%backward_plan,n1,n2,n3, &
         FFTW_BACKWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
#endif

#ifdef FFTW3
    fft3d_info%iwork(1) = n1
    fft3d_info%iwork(2) = n2
    fft3d_info%iwork(3) = n3
    fft3d_info%iwork(4) = ld1
    fft3d_info%iwork(5) = ld2
    fft3d_info%iwork(6) = n3
    fft3d_info%iwork(7) = n2*n3
    fft3d_info%iwork(8) = n3*n1
    fft3d_info%iwork(9) = n1*n2
    fft3d_info%iwork(10) = n2*n3
    fft3d_info%iwork(11) = n3*n1
    fft3d_info%iwork(12) = n1*n2
    ! 3D FFTs
    call dfftw_plan_many_dft(fft3d_info%forward_plan,3,fft3d_info%iwork(1),1, &
         fft3d_info%zwork,fft3d_info%iwork(4),1,0,fft3d_info%zwork, &
         fft3d_info%iwork(4),1,0,FFTW_FORWARD,FFTW_MEASURE)
    call dfftw_plan_many_dft(fft3d_info%backward_plan,3,fft3d_info%iwork(1),1, &
         fft3d_info%zwork,fft3d_info%iwork(4),1,0,fft3d_info%zwork, &
         fft3d_info%iwork(4),1,0,FFTW_BACKWARD,FFTW_MEASURE)

    ! 1D FFTs

    ! n2*n3 transforms along '1'-dir, stride 1, dist n1
    call dfftw_plan_many_dft(fft3d_info%n1_f_plan, 1, & ! plan, rank
         fft3d_info%iwork(1), &   ! n
         fft3d_info%iwork(7), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(4), &   ! inembed
         1, &                     ! istride
         fft3d_info%iwork(4), &   ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(4), &   ! onembed
         1, &                     ! ostride
         fft3d_info%iwork(4), &   ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n1_b_plan, 1, & ! plan, rank
         fft3d_info%iwork(1), &   ! n
         fft3d_info%iwork(7), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(4), &   ! inembed
         1, &                     ! istride
         fft3d_info%iwork(4), &   ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(4), &   ! onembed
         1, &                     ! ostride
         fft3d_info%iwork(4), &   ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! n3 transforms along '2'-dir, stride n1, dist n1*n2
    call dfftw_plan_many_dft(fft3d_info%n2_f_plan, 1, & ! plan, rank
         fft3d_info%iwork(2), &   ! n
         fft3d_info%iwork(3), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(5), &   ! inembed
         fft3d_info%iwork(4), &   ! istride
         fft3d_info%iwork(12), &  ! dist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(5), &   ! onembed
         fft3d_info%iwork(4), &   ! ostride
         fft3d_info%iwork(12), &  ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n2_b_plan, 1, & ! plan, rank
         fft3d_info%iwork(2), &   ! n
         fft3d_info%iwork(3), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(5), &   ! inembed
         fft3d_info%iwork(4), &   ! istride
         fft3d_info%iwork(12), &  ! dist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(5), &   ! onembed
         fft3d_info%iwork(4), &   ! ostride
         fft3d_info%iwork(12), &  ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! n1*n2 transforms along '3'-dir, stride n1*n2, dist 1
    call dfftw_plan_many_dft(fft3d_info%n3_f_plan, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(9), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n3_b_plan, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(9), &   ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! 1 transform along '3'-dir, stride n1*n2, dist 1
    call dfftw_plan_many_dft(fft3d_info%n3_f_plan_1, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         1, &                     ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n3_b_plan_1, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         1, &                     ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! n1/4+1 transforms along '3'-dir, stride n1*n2, dist 1
    call dfftw_plan_many_dft(fft3d_info%n3_f_plan_2, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(1)/4+1, & ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n3_b_plan_2, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(1)/4+1, & ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

    ! n1/4 transforms along '3'-dir, stride n1*n2, dist 1
    call dfftw_plan_many_dft(fft3d_info%n3_f_plan_3, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(1)/4, & ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_FORWARD, &          ! sign
         FFTW_MEASURE)            ! flags
    call dfftw_plan_many_dft(fft3d_info%n3_b_plan_3, 1, & ! plan, rank
         fft3d_info%iwork(3), &   ! n
         fft3d_info%iwork(1)/4, & ! howmany
         fft3d_info%zwork, &      ! in
         fft3d_info%iwork(6), &   ! inembed
         fft3d_info%iwork(9), &   ! istride
         1, &                     ! idist
         fft3d_info%zwork, &      ! out
         fft3d_info%iwork(6), &   ! onembed
         fft3d_info%iwork(9), &   ! ostride
         1, &                     ! odist
         FFTW_BACKWARD, &         ! sign
         FFTW_MEASURE)            ! flags

#endif

!CW
#if defined(FFTW3GPU) 
    fft3d_info%iwork(1) = n1
    fft3d_info%iwork(2) = n2
    fft3d_info%iwork(3) = n3
    fft3d_info%iwork(4) = ld1
    fft3d_info%iwork(5) = ld2
    fft3d_info%iwork(6) = n3
    call cufftplan3d_Z2Z(i,n1,n2,n3,BATCH)
    fft3d_info%forward_plan=i
    call cufftplan3d_Z2Z(i,n1,n2,n3,BATCH)
    fft3d_info%backward_plan=i
#endif
!END CW

#ifdef ACML
    call zfft3dy(ACML_MODE_PLAN,1.0_DP,.true.,n1,n2,n3, &
         fft3d_info%zwork(1),1,ld1,ld1*ld2,fft3d_info%zwork(1),1,ld1,ld1*ld2, &
         fft3d_info%tbl(1),fft3d_info%tbl_size,ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_serial_init &
            &(fourier_mod.F90): zfft3dy failed with code ',ierr
       call comms_abort
    end if
#endif

#ifdef VENDOR
#ifdef ALPHA
    ierr = zfft_init_3d(n1,n2,n3,fft3d_info%iwork(1),.true.)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_serial_init &
            &(fourier_mod.F90): zfft_init_3d failed with code ',ierr
       call comms_abort
    end if
#endif
#ifdef SUN
    call zfftz3(0,n1,n2,n3,1.0_DP,zdum,ld1,ld2,zdum,ld1,ld2, &
         fft3d_info%dwork,fft3d_info%iwork,fft3d_info%zwork, &
         2*fft3d_info%zwork_len,ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_serial_init &
            &(fourier_mod.F90): zfftz3 failed with code ',ierr
       call comms_abort
    end if
#endif
#endif

  end subroutine internal_serial_init

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_parallel_init(grid,fft3d_info)

    use cell_grid, only: GRID_INFO
    use comms, only: comms_abort, pub_total_num_nodes, pub_on_root
    use constants, only: DP, stdout
    use utils, only : utils_alloc_check,utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    type(parallel_fft3d_info), intent(out) :: fft3d_info

    ! Local variables
    integer :: n1,n2,n3    ! Dimensions of FFT grid
    integer :: ld1,ld2,ld3 ! Array dimensions containing grids
    integer :: ierr        ! Error flag
#ifdef SUN
    real(kind=DP) :: ddum
    complex(kind=DP) :: zdum
#endif
!CW
    integer :: i
!END CW

    n1 = grid%n1
    n2 = grid%n2
    n3 = grid%n3
    ld1 = grid%ld1
    ld2 = grid%ld2
    ld3 = grid%ld3

    ! Simple checks on validity of input arguments
    if (n1 < 1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): n1 < 1'
       call comms_abort
    end if

    if (n2 < 1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): n2 < 1'
       call comms_abort
    end if

    if (n3 < 1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): n3 < 1'
       call comms_abort
    end if

    if (ld1 < n1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld1 < n1'
       call comms_abort
    end if

    if (ld2 < n2) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld2 < n2'
       call comms_abort
    end if

    if (ld3 < n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld3 < n3'
       call comms_abort
    end if

    ! Check parallel strategy module has distributed fine cell FFTs

    if (.not. grid%distributed) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): no parallel strategy'
       call comms_abort
    end if

    ! Platform-dependent checking

#ifdef FFTW
    if (ld1 /= n1+2) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld1 must equal n1+2 for FFTw'
       call comms_abort
    end if
    if (ld2 /= n2) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld2 must equal n2 for FFTw'
       call comms_abort
    end if
    if (ld3 /= n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): ld3 must equal n3 for FFTw'
       call comms_abort
    end if
#endif

    ! Simple initialisation

    ! Store basic information in fft3d_info structure

    fft3d_info%n1 = n1
    fft3d_info%n2 = n2
    fft3d_info%n3 = n3
    fft3d_info%ld1 = ld1
    fft3d_info%ld2 = ld2
    fft3d_info%ld3 = ld3
    fft3d_info%iwork_len = 0
    fft3d_info%dwork_len = 0
    fft3d_info%zwork_len = 0

    ! Store parallel strategy information in fft3d_info structure

    fft3d_info%num12slabs = grid%num_my_slabs12
    fft3d_info%max12slabs = grid%max_slabs12

    allocate(fft3d_info%idx12slab(0:pub_total_num_nodes),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%idx12slab',ierr)

    fft3d_info%idx12slab(0:pub_total_num_nodes-1) = &
         grid%first_slab12(0:pub_total_num_nodes-1)
    fft3d_info%idx12slab(pub_total_num_nodes) = n3 + 1

    fft3d_info%num23slabs = grid%num_slabs23
    fft3d_info%max23slabs = grid%max_slabs23

    allocate(fft3d_info%idx23slab(0:pub_total_num_nodes),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%idx23slab',ierr)

    fft3d_info%idx23slab(0:pub_total_num_nodes-1) = &
         grid%first_slab23(0:pub_total_num_nodes-1)
    fft3d_info%idx23slab(pub_total_num_nodes) = ld1/2 + 1

    ! General workspace allocation

    allocate(fft3d_info%buf12slab(fft3d_info%ld1/2,fft3d_info%n2),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%buf12slab',ierr)

    allocate(fft3d_info%buf23slab(fft3d_info%ld3,fft3d_info%n2),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%buf23slab',ierr)

    allocate(fft3d_info%fsendbuf(fft3d_info%max12slabs,fft3d_info%n2, &
         fft3d_info%ld1/2),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%fsendbuf',ierr)

    allocate(fft3d_info%bsendbuf(fft3d_info%max23slabs,fft3d_info%n2, &
         fft3d_info%n3),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%bsendbuf',ierr)

    allocate(fft3d_info%frecvbuf(fft3d_info%max12slabs,fft3d_info%n2, &
         fft3d_info%num23slabs),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%frecvbuf',ierr)

    allocate(fft3d_info%brecvbuf(fft3d_info%max23slabs,fft3d_info%n2, &
         fft3d_info%num12slabs),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%brecvbuf',ierr)

    ! Platform-dependent workspace allocation

#ifdef FFTW
    allocate(fft3d_info%rfftwbuf(fft3d_info%ld1,fft3d_info%n2),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%rfftwbuf',ierr)
#endif

#ifdef FFTW3
    fft3d_info%iwork_len = 7
    fft3d_info%dwork_len = ld1*n2
    fft3d_info%zwork_len = max(n2*ld1/2,ld2*ld3)
#endif

!CW
#ifdef FFTW3GPU
    fft3d_info%iwork_len = 7
    fft3d_info%dwork_len = ld1*n2
    fft3d_info%zwork_len = max(n2*ld1/2,ld2*ld3)
#endif
!END CW

#ifdef ACML
    fft3d_info%zwork_len = max(n2*ld1/2,ld2*ld3)
    fft3d_info%dwork_len = max(n2*ld1/2,ld2*ld3)
    fft3d_info%tbl_size = ld1*ld2*n3 + 4*(ld1+ld2+n3) + 300
    allocate(fft3d_info%tbl(fft3d_info%tbl_size,3),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%tbl',ierr)
    allocate(fft3d_info%rbuf(ld1),stat=ierr)
    call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
         'fft3d_info%rbuf',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    fft3d_info%iwork_len = dxml_structure_size
#endif
#ifdef SUN
    fft3d_info%iwork_len = 2*128
    fft3d_info%dwork_len = 2*max(n1,n2+n3)
    fft3d_info%zwork_len = 2*max(n1,n2,n3)
#endif
#endif

    if (fft3d_info%iwork_len > 0) then
       allocate(fft3d_info%iwork(fft3d_info%iwork_len,3),stat=ierr)
       call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
            'fft3d_info%iwork',ierr)
    end if

    if (fft3d_info%dwork_len > 0) then
       allocate(fft3d_info%dwork(fft3d_info%dwork_len,3),stat=ierr)
       call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
            'fft3d_info%dwork',ierr)
    end if

    if (fft3d_info%zwork_len > 0) then
       allocate(fft3d_info%zwork(fft3d_info%zwork_len,3),stat=ierr)
       call utils_alloc_check('internal_parallel_init (fourier_mod.F90)', &
            'fft3d_info%zwork',ierr)
    end if

    ! Platform-dependent initialisation

#ifdef FFTW
    call rfftw_f77_create_plan(fft3d_info%forward_plan(1),n1, &
         FFTW_REAL_TO_COMPLEX,FFTW_MEASURE+FFTW_USE_WISDOM)
    call fftw2d_f77_create_plan(fft3d_info%forward_plan(2),n3,n2, &
         FFTW_FORWARD,FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
    call fftw2d_f77_create_plan(fft3d_info%backward_plan(2),n3,n2, &
         FFTW_BACKWARD,FFTW_MEASURE+FFTW_USE_WISDOM)
    call rfftw_f77_create_plan(fft3d_info%backward_plan(1),n1, &
         FFTW_COMPLEX_TO_REAL,FFTW_MEASURE+FFTW_USE_WISDOM)
#endif

#ifdef FFTW3
    fft3d_info%iwork(1,1) = n1
    fft3d_info%iwork(2,1) = n3
    fft3d_info%iwork(3,1) = n2
    fft3d_info%iwork(4,1) = ld1
    fft3d_info%iwork(5,1) = ld3
    fft3d_info%iwork(6,1) = ld2
    fft3d_info%iwork(7,1) = ld1/2
    call dfftw_plan_many_dft_r2c(fft3d_info%forward_plan(1),1, &
         fft3d_info%iwork(1,1),n2,fft3d_info%dwork,fft3d_info%iwork(4,1),1, &
         ld1,fft3d_info%zwork,fft3d_info%iwork(7,1),1,ld1/2,FFTW_MEASURE)
    call dfftw_plan_many_dft(fft3d_info%forward_plan(2),2, &
         fft3d_info%iwork(2,1),1,fft3d_info%zwork,fft3d_info%iwork(5,1),1,0, &
         fft3d_info%zwork,fft3d_info%iwork(5,1),1,0,FFTW_FORWARD,FFTW_MEASURE)
    call dfftw_plan_many_dft(fft3d_info%backward_plan(2),2, &
         fft3d_info%iwork(2,1),1,fft3d_info%zwork(1,1),fft3d_info%iwork(5,1), &
         1,0,fft3d_info%zwork(1,2),fft3d_info%iwork(5,1),1,0,FFTW_BACKWARD, &
         FFTW_MEASURE)
    call dfftw_plan_many_dft_c2r(fft3d_info%backward_plan(1),1, &
         fft3d_info%iwork(1,1),n2,fft3d_info%dwork,fft3d_info%iwork(7,1),1, &
         ld1/2,fft3d_info%dwork,fft3d_info%iwork(4,1),1,ld1,FFTW_MEASURE)
#endif

!CW
#ifdef FFTW3GPU
   fft3d_info%iwork(1,1) = n1
   fft3d_info%iwork(2,1) = n3
   fft3d_info%iwork(3,1) = n2
   fft3d_info%iwork(4,1) = ld1
   fft3d_info%iwork(5,1) = ld3
   fft3d_info%iwork(6,1) = ld2
   fft3d_info%iwork(7,1) = ld1/2
!GPU
#ifdef GPU_ALL
  call cufftplan1d_D2Z(i,  n1, BATCH)
  fft3d_info%forward_plan(1)=i
#endif
   call cufftplan2d_Z2Z(i, n3,n2,BATCH)
   fft3d_info%forward_plan(2)=i
   call cufftplan2d_Z2Z(i, n3,n2,BATCH)
   fft3d_info%backward_plan(2)=i
#ifdef GPU_ALL2
  call cufftplan1d_Z2D(i, n1,  BATCH)
  fft3d_info%backward_plan(1)=i
#endif
!CPU
#ifndef GPU_ALL
   call dfftw_plan_many_dft_r2c(fft3d_info%forward_plan(1),1, &
         fft3d_info%iwork(1,1),n2,fft3d_info%dwork,fft3d_info%iwork(4,1),1, &
         ld1,fft3d_info%zwork,fft3d_info%iwork(7,1),1,ld1/2,FFTW_MEASURE)
#endif
!  call dfftw_plan_many_dft(fft3d_info%forward_plan(2),2, &
!       fft3d_info%iwork(2,1),1,fft3d_info%zwork,fft3d_info%iwork(5,1),1,0, &
!       fft3d_info%zwork,fft3d_info%iwork(5,1),1,0,FFTW_FORWARD,FFTW_MEASURE)
!  call dfftw_plan_many_dft(fft3d_info%backward_plan(2),2, &
!       fft3d_info%iwork(2,1),1,fft3d_info%zwork(1,1),fft3d_info%iwork(5,1), &
!       1,0,fft3d_info%zwork(1,2),fft3d_info%iwork(5,1),1,0,FFTW_BACKWARD, &
!       FFTW_MEASURE)
#ifndef GPU_ALL2
   call dfftw_plan_many_dft_c2r(fft3d_info%backward_plan(1),1, &
         fft3d_info%iwork(1,1),n2,fft3d_info%dwork,fft3d_info%iwork(7,1),1, &
         ld1/2,fft3d_info%dwork,fft3d_info%iwork(4,1),1,ld1,FFTW_MEASURE)
#endif
#endif
!END CW

#ifdef ACML
    call dzfft(ACML_MODE_PLAN,n1,fft3d_info%dwork(1,1),fft3d_info%tbl(1,1),ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): dzfft failed with code ',ierr
       call comms_abort
    end if
    call zfft2dx(ACML_MODE_PLAN,1.0_DP,.true.,.true.,n3,n2,fft3d_info%zwork(1,1), &
         1,ld3,fft3d_info%zwork(1,1),1,ld3,fft3d_info%tbl(1,2),ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): dzfft failed with code ',ierr
       call comms_abort
    end if
    call zdfft(ACML_MODE_PLAN,n1,fft3d_info%dwork(1,1),fft3d_info%tbl(1,3),ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): zdfft failed with code ',ierr
       call comms_abort
    end if
#endif

#ifdef VENDOR
#ifdef ALPHA
    ierr = dfft_init(n1,fft3d_info%iwork(1,1),.true.)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): dfft_init failed with code ',ierr
       call comms_abort
    end if
    ierr = zfft_init_2d(n3,n2,fft3d_info%iwork(1,2),.true.)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): zfft_init_2d failed with code ',ierr
       call comms_abort
    end if
#endif
#ifdef SUN
    call dfftzm(0,n1,n2,1.0_DP,ddum,ld1,zdum,ld1/2,fft3d_info%dwork(:,1), &
         fft3d_info%iwork(:,1),fft3d_info%zwork(:,1),2*fft3d_info%zwork_len, &
         ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): dfftzm failed with code ',ierr
       call comms_abort
    end if
    call zfftz2(0,n3,n2,1.0_DP,zdum,ld3,zdum,ld3,fft3d_info%dwork(:,2), &
         fft3d_info%iwork(:,2),fft3d_info%zwork(:,2),2*fft3d_info%zwork_len, &
         ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): zfftz2 failed with code ',ierr
       call comms_abort
    end if
    call zfftdm(0,n1,n2,1.0_DP,zdum,ld1/2,ddum,ld1,fft3d_info%dwork(:,3), &
         fft3d_info%iwork(:,3),fft3d_info%zwork(:,3),2*fft3d_info%zwork_len, &
         ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_init &
            &(fourier_mod.F90): zfftdm failed with code ',ierr
       call comms_abort
    end if
#endif
#endif

  end subroutine internal_parallel_init

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_fft3d_cell_forward(fft3d_info,rspc,gspc)

#ifdef VENDOR
    use comms, only: pub_on_root, comms_abort, pub_my_node_id, &
         pub_total_num_nodes, comms_barrier
    use constants, only: DP, stdout
#else
    use comms, only:  pub_my_node_id, pub_total_num_nodes, comms_barrier, &
         comms_send, comms_recv, comms_abort, pub_on_root
    use constants, only: DP, stdout
#endif

    implicit none

    type(parallel_fft3d_info), intent(inout) :: fft3d_info
    real(kind=DP), intent(in) :: &
         rspc(fft3d_info%ld1,fft3d_info%ld2,fft3d_info%max12slabs)
    complex(kind=DP), intent(out) :: &
         gspc(fft3d_info%ld3,fft3d_info%ld2,fft3d_info%max23slabs)

    ! Local variables
#ifdef ACML
    integer :: ierr                ! Error flag
    real(kind=DP) :: scale         ! Scale Factor
#endif
#ifdef VENDOR
    integer :: ierr                ! Error flag
#endif
    integer :: i1,i2,i3            ! Grid loop counters
    integer :: islab12,islab23     ! Slab loop counters
    integer :: i1s,i1e             ! Start/end of 1-rods to send
    integer :: inode               ! Node loop counter
    integer :: send_id,recv_id     ! Send/receive node

    ! Zero output array
    gspc = (0.0_DP,0.0_DP)

    !--------------------------------------------------------------------!
    !                                                                    !
    ! real   : rspc(ld1,ld2,num12slabs)            =  n(r1,r2,r3)        !
    !                       |                                            !
    !                       | 1D FFT along 1: r1 <-> g1                  !
    !                       V                                            !
    ! complex: buf12slab(ld1/2,n2)                 =  n(g1,r2)           !
    !                       |                         per 12-slab        !
    !                       | transpose                                  !
    !                       V                                            !
    ! complex: fsendbuf(max12slabs,n2,ld1/2)       =  n(r3,r2,g1)        !
    !                       |                                            !
    !                       | communicate                                !
    !                       V                                            !
    ! complex: frecvbuf(max12slabs,n2,num23slabs)  =  n(r3,r2,g1)        !
    !                       |                         per node           !
    !                       | copy                                       !
    !                       V                                            !
    ! complex: gspc(ld3,ld2,num23slabs)            =  n(r3,r2,g1)        !
    !                       |                                            !
    !                       | 2D FFT in (2,3): r2,r3 <-> g2,g3           !
    !                       V                                            !
    ! complex: gspc(ld3,ld2,num23slabs)            =  n(g3,g2,g1)        !
    !                                                                    !
    !--------------------------------------------------------------------!

    ! Loop over 12-slabs on this node

    do islab12=1,fft3d_info%num12slabs

       ! 1D FFT of 1-rods in this 12-slab on this node

#ifdef FFTW
       call rfftw_f77(fft3d_info%forward_plan(1),fft3d_info%n2, &
            rspc(1,1,islab12),1,fft3d_info%ld1,fft3d_info%rfftwbuf, &
            1,fft3d_info%ld1)
       do i2=1,fft3d_info%n2
          fft3d_info%buf12slab(1,i2) = &
               cmplx(fft3d_info%rfftwbuf(1,i2),0.0_DP,kind=DP)
          do i1=2,fft3d_info%ld1/2-1
             fft3d_info%buf12slab(i1,i2) = &
                  cmplx(fft3d_info%rfftwbuf(i1,i2), &
                  fft3d_info%rfftwbuf(fft3d_info%ld1-i1,i2),kind=DP)
          end do
          fft3d_info%buf12slab(fft3d_info%ld1/2,i2) = &
               cmplx(fft3d_info%rfftwbuf(fft3d_info%ld1/2,i2),0.0_DP,kind=DP)
       end do
#endif

#ifdef FFTW3
       call dfftw_execute_dft_r2c(fft3d_info%forward_plan(1), &
            rspc(1,1,islab12),fft3d_info%buf12slab)
#endif

!CW
#ifdef FFTW3GPU
#ifndef GPU_ALL
      call dfftw_execute_dft_r2c(fft3d_info%forward_plan(1),rspc(1,1,islab12),fft3d_info%buf12slab)
#else
      call cufftexecd2z(int(fft3d_info%forward_plan(1)), rspc(1,1,islab12), fft3d_info%buf12slab, FFTW_FORWARD,fft3d_info%n2,BATCH)
#endif
#endif
!END CW

#ifdef ACML
       do i2=1,fft3d_info%n2
          fft3d_info%rbuf(:) = rspc(:,i2,islab12)
          call dzfft(ACML_MODE_FORWARD,fft3d_info%n1,fft3d_info%rbuf(1), &
               fft3d_info%tbl(1,1),ierr)
          fft3d_info%buf12slab(1,i2) = &
               cmplx(fft3d_info%rbuf(1),0.0_DP,kind=DP)
          do i1=2,fft3d_info%ld1/2-1
             fft3d_info%buf12slab(i1,i2) = &
                  cmplx(fft3d_info%rbuf(i1), &
                  fft3d_info%rbuf(fft3d_info%ld1-i1),kind=DP)
          end do
          fft3d_info%buf12slab(fft3d_info%ld1/2,i2) = &
               cmplx(fft3d_info%rbuf(fft3d_info%ld1/2),0.0_DP,kind=DP)
       end do
       if (ierr /= 0) then
          if (pub_on_root.or..true.) write(stdout,'(a,i6)') &
               'Error in internal_fft3d_cell_forward (fourier_mod.F90): &
               &dzfft failed with code ',ierr
          call comms_abort
       end if
       scale = sqrt(real(fft3d_info%n1,kind=DP))
       fft3d_info%buf12slab = fft3d_info%buf12slab * scale
#endif

#ifdef VENDOR
#ifdef ALPHA
       do i2=1,fft3d_info%n2
          ierr = dfft_apply('R','C','F',rspc(1,i2,islab12), &
               fft3d_info%buf12slab(1,i2),fft3d_info%iwork(1,1),1)
          if (ierr /= 0) then
             if (pub_on_root) write(stdout,'(a,i6)') &
                  'Error in internal_fft3d_cell_forward (fourier_mod.F90): &
                  &dfft_apply failed with code ',ierr
             call comms_abort
          end if
       end do
#endif
#ifdef SUN
       call dfftzm(-1,fft3d_info%n1,fft3d_info%n2,1.0_DP,rspc(:,:,islab12), &
            fft3d_info%ld1,fft3d_info%buf12slab,fft3d_info%ld1/2, &
            fft3d_info%dwork(:,1),fft3d_info%iwork(:,1), &
            fft3d_info%zwork(:,1),2*fft3d_info%zwork_len,ierr)
       if (ierr /= 0) then
          if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_fft3d_cell&
               &_forward (fourier_mod.F90): dfftzm failed with code ',ierr
          call comms_abort
       end if
#endif
#endif

       ! Transpose 1-rods in this 12-slab into 3-rods in 23-slab
       do i2=1,fft3d_info%n2
          do i1=1,fft3d_info%ld1/2
             fft3d_info%fsendbuf(islab12,i2,i1) = fft3d_info%buf12slab(i1,i2)
          end do
       end do

    end do

    ! Local communication phase
    i3 = fft3d_info%idx12slab(pub_my_node_id)
    do islab12=1,fft3d_info%num12slabs
       do i2=1,fft3d_info%n2
          i1 = fft3d_info%idx23slab(pub_my_node_id)
          do islab23=1,fft3d_info%num23slabs
             gspc(i3,i2,islab23) = fft3d_info%fsendbuf(islab12,i2,i1)
             i1 = i1 + 1
          end do
       end do
       i3 = i3 + 1
    end do

    ! Global communication phase
    call comms_barrier
    do inode=1,pub_total_num_nodes-1
       send_id = modulo(pub_my_node_id + inode,pub_total_num_nodes)
       recv_id = modulo(pub_my_node_id - inode,pub_total_num_nodes)

       ! Send packet to node send_id
       i1s = fft3d_info%idx23slab(send_id)
       i1e = fft3d_info%idx23slab(send_id+1)-1
       call comms_send(send_id,fft3d_info%fsendbuf(:,:,i1s:i1e))

       ! Receive packet from node recv_id
       call comms_recv(recv_id,fft3d_info%frecvbuf)

       ! Copy received data into gspc array
       do islab23=1,fft3d_info%num23slabs
          do i2=1,fft3d_info%n2
             i3 = fft3d_info%idx12slab(recv_id)
             do islab12=1,fft3d_info%idx12slab(recv_id+1)-i3
                gspc(i3,i2,islab23) = fft3d_info%frecvbuf(islab12,i2,islab23)
                i3 = i3 + 1
             end do
          end do
       end do

    end do

    ! Loop over 23-slabs on this node

    do islab23=1,fft3d_info%num23slabs

       ! 2D FFT of this 23-slab on this node

#ifdef FFTW
       call fftwnd_f77_one(fft3d_info%forward_plan(2),gspc(1,1,islab23), &
            gspc(1,1,islab23))
#endif

#ifdef FFTW3
       call dfftw_execute_dft(fft3d_info%forward_plan(2),gspc(1,1,islab23), &
            gspc(1,1,islab23))
#endif

!CW
#ifdef FFTW3GPU
      !call dfftw_execute_dft(fft3d_info%forward_plan(2),gspc(1,1,islab23), gspc(1,1,islab23))
       call cufftexecz2z(int(fft3d_info%forward_plan(2)), gspc(1,1,islab23), gspc(1,1,islab23),FFTW_FORWARD,fft3d_info%n3*fft3d_info%n2,BATCH)
#endif
!END CW

#ifdef ACML
       call zfft2dx(ACML_MODE_BACKWARD,1.0_DP,.true.,.true.,fft3d_info%n3,fft3d_info%n2, &
            gspc(1,1,islab23),1,fft3d_info%ld3,gspc(1,1,islab23),1,fft3d_info%ld3, &
            fft3d_info%tbl(1,2),ierr)
       if (ierr /= 0) then
          if (pub_on_root.or..true.) write(stdout,'(a,i6)') &
               'Error in internal_fft3d_cell_forward (fourier_mod.F90): &
               &zfft2dx failed with code ',ierr
          call comms_abort
       end if
#endif

#ifdef VENDOR
#ifdef ALPHA
       ierr = zfft_apply_2d('C','C','F',gspc(1,1,islab23), &
            gspc(1,1,islab23),fft3d_info%ld3,fft3d_info%iwork(1,2),1,1)
       if (ierr /= 0) then
          if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_fft3d_cell&
               &_forward (fourier_mod.F90): zfft_apply_2d failed with code ', &
               ierr
          call comms_abort
       end if
#endif
#ifdef SUN
       call zfftz2(-1,fft3d_info%n3,fft3d_info%n2,1.0_DP,gspc(:,:,islab23), &
            fft3d_info%ld3,gspc(:,:,islab23),fft3d_info%ld3, &
            fft3d_info%dwork(:,2),fft3d_info%iwork(:,2), &
            fft3d_info%zwork(:,2),2*fft3d_info%zwork_len,ierr)
       if (ierr /= 0) then
          if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_fft3d_cell&
               &_forward (fourier_mod.F90): zfftz2 failed with code ',ierr
          call comms_abort
      end if
#endif
#endif

    end do

    ! Synchronise all nodes
    call comms_barrier

  end subroutine internal_fft3d_cell_forward

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_fft3d_cell_backward(fft3d_info,rspc,gspc)

#ifdef VENDOR
    use comms, only: pub_on_root, comms_abort, pub_my_node_id, &
         pub_total_num_nodes, comms_barrier
    use constants, only: DP, stdout
#else
    use comms, only:  pub_my_node_id, pub_total_num_nodes, comms_barrier, &
         comms_send, comms_recv, comms_abort, pub_on_root
    use constants, only: stdout
#endif

    implicit none

    type(parallel_fft3d_info), intent(inout) :: fft3d_info
    real(kind=DP), intent(out) :: &
         rspc(fft3d_info%ld1,fft3d_info%ld2,fft3d_info%max12slabs)
    complex(kind=DP), intent(in) :: &
         gspc(fft3d_info%ld3,fft3d_info%ld2,fft3d_info%max23slabs)

    ! Local variables
#ifdef ACML
    integer :: ierr                ! Error flag
#endif
#ifdef VENDOR
    integer :: ierr                ! Error flag
#endif
    integer :: i1,i2,i3            ! Grid loop counters
    integer :: islab12,islab23     ! Slab loop counters
    integer :: i3s,i3e             ! Start/end of 3-rods to send
    integer :: inode               ! Node loop counter
    integer :: send_id,recv_id     ! Send/receive node
    real(kind=DP) :: scale

    ! Zero output array
    rspc = 0.0_DP

    !--------------------------------------------------------------------!
    !                                                                    !
    ! complex: gspc(ld3,ld2,num23slabs)            =  n(g3,g2,g1)        !
    !                       |                                            !
    !                       | 2D FFT in (2,3): g2,g3 <-> r2,r3           !
    !                       V                                            !
    ! complex: buf23slab(ld3,n2)                   =  n(r3,r2)           !
    !                       |                         per 23-slab        !
    !                       | transpose                                  !
    !                       V                                            !
    ! complex: bsendbuf(max23slabs,n2,n3)          =  n(g1,r2,r3)        !
    !                       |                                            !
    !                       | communicate                                !
    !                       V                                            !
    ! complex: brecvbuf(max23slabs,n2,num12slabs)  =  n(g1,r2,r3)        !
    !                       |                         per node           !
    !                       | copy                                       !
    !                       V                                            !
    ! real   : rspc(ld1,ld2,num12slabs)            =  n(g1,r2,r3)        !
    !                       |                                            !
    !                       | 1D FFT along 1: g1 <-> r1                  !
    !                       V                                            !
    ! real   : rspc(ld1,ld2,num12slabs)            =  n(r1,r2,r3)        !
    !                                                                    !
    !--------------------------------------------------------------------!

    ! Loop over 23-slabs on this node

    do islab23=1,fft3d_info%num23slabs

       ! 2D FFT of this 23-slab on this node

#ifdef FFTW
       call fftwnd_f77_one(fft3d_info%backward_plan(2),gspc(1,1,islab23), &
            fft3d_info%buf23slab)
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       fft3d_info%buf23slab = fft3d_info%buf23slab * scale
#endif

#ifdef FFTW3
       call dfftw_execute_dft(fft3d_info%backward_plan(2), &
            gspc(1,1,islab23),fft3d_info%buf23slab)
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       fft3d_info%buf23slab = fft3d_info%buf23slab * scale
#endif

!CW
#ifdef FFTW3GPU
        !call dfftw_execute_dft(fft3d_info%backward_plan(2), gspc(1,1,islab23),fft3d_info%buf23slab)
       call cufftexecz2z(int(fft3d_info%backward_plan(2)), gspc(1,1,islab23),fft3d_info%buf23slab,FFTW_BACKWARD,fft3d_info%n3*fft3d_info%n2,BATCH)
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       fft3d_info%buf23slab = fft3d_info%buf23slab * scale
#endif
!END CW

#ifdef ACML
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       call zfft2dx(ACML_MODE_FORWARD,scale,.true.,.false.,fft3d_info%n3,fft3d_info%n2, &
            gspc(1,1,islab23),1,fft3d_info%ld3,fft3d_info%buf23slab(1,1),1,fft3d_info%ld3, &
            fft3d_info%tbl(1,2),ierr)
       if (ierr /= 0) then
          if (pub_on_root.or..true.) write(stdout,'(a,i6)') &
               'Error in internal_fft3d_cell_backward (fourier_mod.F90): &
               &zfft2dx failed with code ',ierr
          call comms_abort
       end if
#endif

#ifdef VENDOR
#ifdef ALPHA
       ierr = zfft_apply_2d('C','C','B',gspc(:,:,islab23), &
            fft3d_info%buf23slab,fft3d_info%ld3,fft3d_info%iwork(:,2),1,1)
       if (ierr /= 0) then
          if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_fft3d_cell&
               &_backward (fourier_mod.F90): zfft_apply_2d failed with code ', &
               ierr
          call comms_abort
       end if
#endif
#ifdef SUN
       scale = 1.0_DP / (fft3d_info%n2*fft3d_info%n3)
       call zfftz2(1,fft3d_info%n3,fft3d_info%n2,scale,gspc(:,:,islab23), &
            fft3d_info%ld3,fft3d_info%buf23slab,fft3d_info%ld3, &
            fft3d_info%dwork(:,2),fft3d_info%iwork(:,2), &
            fft3d_info%zwork(:,2),2*fft3d_info%zwork_len,ierr)
       if (ierr /= 0) then
          if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_fft3d_cell&
               &_backward (fourier_mod.F90): zfftz2 failed with code ',ierr
          call comms_abort
       end if
#endif
#endif

       ! Transpose 3-rods in this 23-slab into 1-rods in 12-slab
       do i2=1,fft3d_info%n2
          do i3=1,fft3d_info%n3
             fft3d_info%bsendbuf(islab23,i2,i3) = fft3d_info%buf23slab(i3,i2)
          end do
       end do

    end do

    ! Local communication phase
    i3 = fft3d_info%idx12slab(pub_my_node_id)
    do islab12=1,fft3d_info%num12slabs
       do i2=1,fft3d_info%n2
          i1 = fft3d_info%idx23slab(pub_my_node_id)
          do islab23=1,fft3d_info%num23slabs
             rspc(2*i1-1,i2,islab12) = &
                  real(fft3d_info%bsendbuf(islab23,i2,i3),kind=DP)
             rspc(2*i1,i2,islab12) = &
                  aimag(fft3d_info%bsendbuf(islab23,i2,i3))
             i1 = i1 + 1
          end do
       end do
       i3 = i3 + 1
    end do

    ! Global communication phase
    call comms_barrier
    do inode=1,pub_total_num_nodes-1
       send_id = modulo(pub_my_node_id + inode,pub_total_num_nodes)
       recv_id = modulo(pub_my_node_id - inode,pub_total_num_nodes)

       ! Send packet to node send_id
       i3s = fft3d_info%idx12slab(send_id)
       i3e = fft3d_info%idx12slab(send_id+1)-1
       call comms_send(send_id,fft3d_info%bsendbuf(:,:,i3s:i3e))

       ! Receive packet from node recv_id
       call comms_recv(recv_id,fft3d_info%brecvbuf)

       ! Copy received data into rspc array
       do islab12=1,fft3d_info%num12slabs
          do i2=1,fft3d_info%n2
             i1 = fft3d_info%idx23slab(recv_id)
             do islab23=1,fft3d_info%idx23slab(recv_id+1)-i1
                rspc(2*i1-1,i2,islab12) = &
                     real(fft3d_info%brecvbuf(islab23,i2,islab12),kind=DP)
                rspc(2*i1,i2,islab12) = &
                     aimag(fft3d_info%brecvbuf(islab23,i2,islab12))
                i1 = i1 + 1
             end do
          end do
       end do

    end do

    ! Loop over 12-slabs on this node

    do islab12=1,fft3d_info%num12slabs

       ! 1D FFT of 1-rods in this 12-slab on this node

#ifdef FFTW
       do i2=1,fft3d_info%n2
          fft3d_info%rfftwbuf(1,i2) = rspc(1,i2,islab12)
          do i1=2,fft3d_info%ld1/2-1
             fft3d_info%rfftwbuf(i1,i2) = rspc(2*i1-1,i2,islab12)
             fft3d_info%rfftwbuf(fft3d_info%ld1-i1,i2) = &
                  rspc(2*i1,i2,islab12)
          end do
          fft3d_info%rfftwbuf(fft3d_info%ld1/2,i2) = &
               rspc(fft3d_info%ld1-1,i2,islab12)
       end do
       call rfftw_f77(fft3d_info%backward_plan(1),fft3d_info%n2, &
            fft3d_info%rfftwbuf,1,fft3d_info%ld1,rspc(1,1,islab12), &
            1,fft3d_info%ld1)
       scale = 1.0_DP / fft3d_info%n1
       rspc(:,:,islab12) = rspc(:,:,islab12) * scale
#endif

#ifdef FFTW3
       call dfftw_execute_dft_c2r(fft3d_info%backward_plan(1), &
            rspc(1,1,islab12),rspc(1,1,islab12))
       scale = 1.0_DP / fft3d_info%n1
       rspc(:,:,islab12) = rspc(:,:,islab12) * scale
#endif

!CW
#ifdef FFTW3GPU
#ifdef GPU_ALL2
       call cufftexecz2d(int(fft3d_info%backward_plan(1)),rspc(1,1,islab12),rspc(1,1,islab12),FFTW_BACKWARD,fft3d_info%n1,BATCH)
#else
       call dfftw_execute_dft_c2r(fft3d_info%backward_plan(1), rspc(1,1,islab12),rspc(1,1,islab12))
#endif
       scale = 1.0_DP / fft3d_info%n1
       rspc(:,:,islab12) = rspc(:,:,islab12) * scale
#endif
!END CW

#ifdef ACML
       do i2=1,fft3d_info%n2
          fft3d_info%rbuf(1) = rspc(1,i2,islab12)
          do i1=2,fft3d_info%ld1/2-1
             fft3d_info%rbuf(i1) = rspc(2*i1-1,i2,islab12)
             fft3d_info%rbuf(fft3d_info%ld1-i1) = &
                  -rspc(2*i1,i2,islab12)
          end do
          fft3d_info%rbuf(fft3d_info%ld1/2) = &
               rspc(fft3d_info%ld1-1,i2,islab12)
          call zdfft(ACML_MODE_FORWARD,fft3d_info%n1,fft3d_info%rbuf(1), &
               fft3d_info%tbl(1,3),ierr)
          scale = sqrt(1.0_DP / fft3d_info%n1)
          rspc(:,i2,islab12) = fft3d_info%rbuf(:) * scale
       end do
       if (ierr /= 0) then
          if (pub_on_root.or..true.) write(stdout,'(a,i6)') &
               'Error in internal_fft3d_cell_backward (fourier_mod.F90): &
               &zdfft failed with code ',ierr
          call comms_abort
       end if
#endif

#ifdef VENDOR
#ifdef ALPHA
       do i2=1,fft3d_info%n2
          ierr = dfft_apply('C','R','B',rspc(1,i2,islab12), &
               rspc(1,i2,islab12),fft3d_info%iwork(1,1),1)
          if (ierr /= 0) then
             if (pub_on_root) write(stdout,'(a,i6)') &
                  'Error in internal_fft3d_cell_backward (fourier_mod.F90): &
                  &dfft_apply failed with code ',ierr
             call comms_abort
          end if
       end do
#endif
#ifdef SUN
       scale = 1.0_DP / fft3d_info%n1
       call zfftdm(1,fft3d_info%n1,fft3d_info%n2,scale,rspc(:,:,islab12), &
            fft3d_info%ld1/2,rspc(:,:,islab12),fft3d_info%ld1, &
            fft3d_info%dwork(:,3),fft3d_info%iwork(:,3), &
            fft3d_info%zwork(:,3),2*fft3d_info%zwork_len,ierr)
       if (ierr /= 0) then
          if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_fft3d_cell&
               &_backward (fourier_mod.F90): zfftdm failed with code ',ierr
          call comms_abort
       end if
#endif
#endif

    end do

    ! Synchronise all nodes
    call comms_barrier

  end subroutine internal_fft3d_cell_backward

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_fft3d_box(fft3d_info,dir,data,apply_scale,padbox)

    use comms, only : pub_on_root, comms_abort
    use constants, only: stdout

    implicit none

    type(serial_fft3d_info), intent(inout) :: fft3d_info
    character, intent(in) :: dir
    complex(kind=DP), intent(inout) :: &
         data(fft3d_info%ld1,fft3d_info%ld2,fft3d_info%n3)
    logical, intent(in) :: apply_scale
    logical, intent(in), optional :: padbox

    ! Local variables
#ifdef ACML
    integer :: ierr
#endif
#ifdef VENDOR
    integer :: ierr             ! Error flag
#endif
    real(kind=DP) :: scale
#ifdef FFTW
    !   complex(kind=DP) :: zdum
    complex(kind=DP) :: zdum(1) !cks, 13/10/2004: Hack to compile with NAG f95
#endif

    integer :: i1,i2,j1
    logical :: do_padbox

    do_padbox = .false.
    if (present(padbox)) do_padbox = padbox

    ! Platform-dependent call

#ifdef FFTW
    if (dir == 'F' .or. dir == 'f') then
       call fftwnd_f77_one(fft3d_info%forward_plan,data,zdum)
    else
       call fftwnd_f77_one(fft3d_info%backward_plan,data,zdum)
       scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
       if (apply_scale) then
          data = data * scale
       end if
    end if
#endif

#ifdef FFTW3
    if (dir == 'F' .or. dir == 'f') then
       if (.not.do_padbox) then
          call dfftw_execute_dft(fft3d_info%forward_plan,data,data)
       else

          ! ndmh: padded version: assumes G-vectors where:
          ! ndmh:      n1/4+2 <= i1 <= 3*n1/4+1
          ! ndmh:  or  n2/4+2 <= i2 <= 3*n2/4+1
          ! ndmh:  or  n3/4+2 <= i3 <= 3*n3/4+1
          ! ndmh: are not required, so avoids as many of the 1D transforms
          ! ndmh: as possible

          ! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
          call dfftw_execute_dft(fft3d_info%n1_f_plan,data(1,1,1),data(1,1,1))
          ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
          do i1=1,fft3d_info%n1/4+1
             call dfftw_execute_dft(fft3d_info%n2_f_plan,data(i1,1,1),data(i1,1,1))
          end do
          ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
          do i1=3*fft3d_info%n1/4+2,fft3d_info%n1
             call dfftw_execute_dft(fft3d_info%n2_f_plan,data(i1,1,1),data(i1,1,1))
          end do
          ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
          ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
          j1 = 3*fft3d_info%n1/4+2
          do i2=1,fft3d_info%n2/4+1
              call dfftw_execute_dft(fft3d_info%n3_f_plan_2,data(1,i2,1),data(1,i2,1))
              call dfftw_execute_dft(fft3d_info%n3_f_plan_3,data(j1,i2,1),data(j1,i2,1))
          end do
          ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
          ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
          do i2=3*fft3d_info%n2/4+2,fft3d_info%n2
              call dfftw_execute_dft(fft3d_info%n3_f_plan_2,data(1,i2,1),data(1,i2,1))
              call dfftw_execute_dft(fft3d_info%n3_f_plan_3,data(j1,i2,1),data(j1,i2,1))
          end do

#if 0
          ! ndmh: full multi-1D-FFT version, equivalent to full 3D transform
          call dfftw_execute_dft(fft3d_info%n3_f_plan,data(1,1,1),data(1,1,1))
          do i1=1,fft3d_info%n1
             call dfftw_execute_dft(fft3d_info%n2_f_plan,data(i1,1,1),data(i1,1,1))
          end do
          call dfftw_execute_dft(fft3d_info%n1_f_plan,data(1,1,1),data(1,1,1))
#endif

       end if
    else
       if (.not.do_padbox) then
          ! ndmh: normal full 3D FFT version
          call dfftw_execute_dft(fft3d_info%backward_plan,data,data)
       else

          ! ndmh: padded version: assumes G-vectors where:
          ! ndmh:      n1/4+2 <= i1 <= 3*n1/4+1
          ! ndmh:  or  n2/4+2 <= i2 <= 3*n2/4+1
          ! ndmh:  or  n3/4+2 <= i3 <= 3*n3/4+1
          ! ndmh: are all zero, so avoids as many of the 1D transforms
          ! ndmh: as possible

          ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
          ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1 = 1
          j1 = 3*fft3d_info%n1/4+2
          do i2=1,fft3d_info%n2/4+1
              call dfftw_execute_dft(fft3d_info%n3_b_plan_2,data(1,i2,1),data(1,i2,1))
              call dfftw_execute_dft(fft3d_info%n3_b_plan_3,data(j1,i2,1),data(j1,i2,1))
          end do
          ! ndmh: do n1/4+1 transforms along '3' dir for each i2 where
          ! ndmh: 1 <= i2 <= n2/4+1 or 3*n2/4+2 <= i2 <= n2, starting at i1=3*n1/4+2
          do i2=3*fft3d_info%n2/4+2,fft3d_info%n2
              call dfftw_execute_dft(fft3d_info%n3_b_plan_2,data(1,i2,1),data(1,i2,1))
              call dfftw_execute_dft(fft3d_info%n3_b_plan_3,data(j1,i2,1),data(j1,i2,1))
          end do
          ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
          do i1=1,fft3d_info%n1/4+1
             call dfftw_execute_dft(fft3d_info%n2_b_plan,data(i1,1,1),data(i1,1,1))
          end do
          ! ndmh: do n3 transforms along the '2' dir, starting at i1 = 1
          do i1=3*fft3d_info%n1/4+2,fft3d_info%n1
             call dfftw_execute_dft(fft3d_info%n2_b_plan,data(i1,1,1),data(i1,1,1))
          end do
          ! ndmh: do n2*n3 transforms along the '1' dir, starting at i1 = 1,i2 = 1
          call dfftw_execute_dft(fft3d_info%n1_b_plan,data(1,1,1),data(1,1,1))

#if 0
          ! ndmh: full multi-1D-FFT version, equivalent to full 3D transform
          call dfftw_execute_dft(fft3d_info%n3_b_plan,data(1,1,1),data(1,1,1))
          do i1=1,fft3d_info%n1
             call dfftw_execute_dft(fft3d_info%n2_b_plan,data(i1,1,1),data(i1,1,1))
          end do
          call dfftw_execute_dft(fft3d_info%n1_b_plan,data(1,1,1),data(1,1,1))
#endif
       end if
       scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
       if (apply_scale) then
          data = data * scale
       end if
    end if
#endif

!CW
#if defined(FFTW3GPU) 
    if (dir == 'F' .or. dir == 'f') then
       call cufftExecZ2Z(int(fft3d_info%forward_plan),data,data,FFTW_FORWARD,fft3d_info%n1*fft3d_info%n2*fft3d_info%n3,batch)
    else
       call cufftExecZ2Z(int(fft3d_info%backward_plan),data,data,FFTW_BACKWARD,fft3d_info%n1*fft3d_info%n2*fft3d_info%n3,batch)
       scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
       if (apply_scale) then
          data = data * scale
       end if
    end if
#endif
!END CW

#ifdef ACML
    ! ndmh: Note ACML's definition of forward and backward seem to be reversed
    ! ndmh: with respect to everyone elses' definitions, hence the flags are
    ! ndmh: passed in seemingly the wrong way round.
    if (dir == 'F' .or. dir == 'f') then
       call zfft3dy(ACML_MODE_BACKWARD,1.0_DP,.true.,fft3d_info%n1, &
            fft3d_info%n2,fft3d_info%n3,data(1,1,1),1,fft3d_info%ld1, &
            fft3d_info%ld1*fft3d_info%ld2,data(1,1,1),1,fft3d_info%ld1, &
            fft3d_info%ld1*fft3d_info%ld2,fft3d_info%tbl(1), &
            fft3d_info%tbl_size,ierr)
    else
       if (apply_scale) then
          scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
       else
          scale = 1.0_DP
       end if
       call zfft3dy(ACML_MODE_FORWARD,scale,.true.,fft3d_info%n1, &
            fft3d_info%n2,fft3d_info%n3,data(1,1,1),1,fft3d_info%ld1, &
            fft3d_info%ld1*fft3d_info%ld2,data(1,1,1),1,fft3d_info%ld1, &
            fft3d_info%ld1*fft3d_info%ld2,fft3d_info%tbl(1), &
            fft3d_info%tbl_size,ierr)
    end if
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_serial_init &
            &(fourier_mod.F90): zfft3dy failed with code ',ierr
       call comms_abort
    end if
#endif

#ifdef VENDOR
#ifdef ALPHA
    if (dir == 'F' .or. dir == 'f') then
       ierr = zfft_apply_3d('C','C','F',data,data,fft3d_info%ld1, &
            fft3d_info%ld2,fft3d_info%iwork,1,1,1)
    else
       ierr = zfft_apply_3d('C','C','B',data,data,fft3d_info%ld1, &
            fft3d_info%ld2,fft3d_info%iwork,1,1,1)
    end if
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_fft3d_box &
            &(fourier_mod.F90): zfft_apply_3d failed with code ',ierr
       call comms_abort
    end if
#endif
#ifdef SUN
    if (dir == 'F' .or. dir == 'f') then
       call zfftz3(-1,fft3d_info%n1,fft3d_info%n2,fft3d_info%n3,1.0_DP, &
            data,fft3d_info%ld1,fft3d_info%ld2,data,fft3d_info%ld1, &
            fft3d_info%ld2,fft3d_info%dwork,fft3d_info%iwork, &
            fft3d_info%zwork,2*fft3d_info%zwork_len,ierr)
    else
       if (apply_scale) then
          scale = 1.0_DP / (fft3d_info%n1*fft3d_info%n2*fft3d_info%n3)
       else
          scale = 1.0_DP
       end if
       call zfftz3(1,fft3d_info%n1,fft3d_info%n2,fft3d_info%n3,scale, &
            data,fft3d_info%ld1,fft3d_info%ld2,data,fft3d_info%ld1, &
            fft3d_info%ld2,fft3d_info%dwork,fft3d_info%iwork, &
            fft3d_info%zwork,2*fft3d_info%zwork_len,ierr)
    end if
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_fft3d_box &
            &(fourier_mod.F90): zfftz3 failed with code ',ierr
       call comms_abort
    end if
#endif
#endif

  end subroutine internal_fft3d_box

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_serial_exit(fft3d_info)

    use comms, only : pub_on_root, comms_abort
    use utils, only : utils_dealloc_check

    implicit none

    type(serial_fft3d_info), intent(inout) :: fft3d_info

    ! Local variables

    integer :: ierr             ! Error flag

    ! Platform-dependent finalisation

#ifdef FFTW
    call fftwnd_f77_destroy_plan(fft3d_info%forward_plan)
    call fftwnd_f77_destroy_plan(fft3d_info%backward_plan)
#endif

#ifdef FFTW3
    call dfftw_destroy_plan(fft3d_info%forward_plan)
    call dfftw_destroy_plan(fft3d_info%backward_plan)
#endif

!CW
#if defined(FFTW3GPU)
    call cufftdestroy(int(fft3d_info%forward_plan))
    call cufftdestroy(int(fft3d_info%backward_plan))
#endif
!END CW

#ifdef ACML
    deallocate(fft3d_info%tbl,stat=ierr)
    call utils_dealloc_check('internal_serial_exit (fourier_mod.F90)', &
         'fft3d_info%tbl',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    ierr = zfft_exit_3d(fft3d_info%iwork)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_serial_exit &
            &(fourier_mod.F90): zfft_exit_3d failed with code ',ierr
       call comms_abort
    end if
#endif
#endif

    ! Workspace deallocation

    if (allocated(fft3d_info%zwork)) then
       deallocate(fft3d_info%zwork,stat=ierr)
       call utils_dealloc_check('internal_serial_exit (fourier_mod.F90)', &
            'fft3d_info%zwork',ierr)
    end if
    if (allocated(fft3d_info%dwork)) then
       deallocate(fft3d_info%dwork,stat=ierr)
       call utils_dealloc_check('internal_serial_exit (fourier_mod.F90)', &
            'fft3d_info%dwork',ierr)
    end if
    if (allocated(fft3d_info%iwork)) then
       deallocate(fft3d_info%iwork,stat=ierr)
       call utils_dealloc_check('internal_serial_exit (fourier_mod.F90)', &
            'fft3d_info%iwork',ierr)
    end if

  end subroutine internal_serial_exit

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine internal_parallel_exit(fft3d_info)

    use comms, only : pub_on_root, comms_abort
    use utils, only : utils_dealloc_check

    implicit none

    type(parallel_fft3d_info), intent(inout) :: fft3d_info

    ! Local variables

    integer :: ierr             ! Error flag

    ! Platform-dependent finalisation

#ifdef FFTW
    call rfftw_f77_destroy_plan(fft3d_info%backward_plan(1))
    call fftwnd_f77_destroy_plan(fft3d_info%backward_plan(2))
    call fftwnd_f77_destroy_plan(fft3d_info%forward_plan(2))
    call rfftw_f77_destroy_plan(fft3d_info%forward_plan(1))

    deallocate(fft3d_info%rfftwbuf,stat=ierr)
    call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
         'fft3d_info%rfftwbuf',ierr)
#endif

#ifdef FFTW3
    call dfftw_destroy_plan(fft3d_info%backward_plan(1))
    call dfftw_destroy_plan(fft3d_info%backward_plan(2))
    call dfftw_destroy_plan(fft3d_info%forward_plan(2))
    call dfftw_destroy_plan(fft3d_info%forward_plan(1))
#endif

!CW
#ifdef FFTW3GPU
!GPU
#ifdef GPU_ALL
  call cufftdestroy(int(fft3d_info%backward_plan(1)))
#endif
  call cufftdestroy(int(fft3d_info%backward_plan(2)))
  call cufftdestroy(int(fft3d_info%forward_plan(2)))
#ifdef GPU_ALL2
 call cufftdestroy(int(fft3d_info%forward_plan(1)))
#endif
!CPU
#ifndef GPU_ALL
  call dfftw_destroy_plan(fft3d_info%backward_plan(1))
#endif
 !call dfftw_destroy_plan(fft3d_info%backward_plan(2))
 !call dfftw_destroy_plan(fft3d_info%forward_plan(2))
#ifndef GPU_ALL2
  call dfftw_destroy_plan(fft3d_info%forward_plan(1))
#endif
#endif
!END CW

#ifdef ACML
    deallocate(fft3d_info%rbuf,stat=ierr)
    call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
         'fft3d_info%rbuf',ierr)
    deallocate(fft3d_info%tbl,stat=ierr)
    call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
         'fft3d_info%tbl',ierr)
#endif

#ifdef VENDOR
#ifdef ALPHA
    ierr = zfft_exit_2d(fft3d_info%iwork(1,2))
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_exit &
            &(fourier_mod.F90): zfft_exit_2d failed with code ',ierr
       call comms_abort
    end if
    ierr = dfft_exit(fft3d_info%iwork(1,1))
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i6)') 'Error in internal_parallel_exit &
            &(fourier_mod.F90): dfft_exit failed with code ',ierr
       call comms_abort
    end if
#endif
#endif

    ! Workspace deallocation

    if (allocated(fft3d_info%zwork)) then
       deallocate(fft3d_info%zwork,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%zwork',ierr)
    end if

    if (allocated(fft3d_info%dwork)) then
       deallocate(fft3d_info%dwork,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%dwork',ierr)
    end if

    if (allocated(fft3d_info%iwork)) then
       deallocate(fft3d_info%iwork,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%iwork',ierr)
    end if

    ! General workspace deallocation

    if (allocated(fft3d_info%brecvbuf)) then
       deallocate(fft3d_info%brecvbuf,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%brecvbuf',ierr)
    end if

    if (allocated(fft3d_info%frecvbuf)) then
       deallocate(fft3d_info%frecvbuf,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%frecvbuf',ierr)
    end if

    if (allocated(fft3d_info%bsendbuf)) then
       deallocate(fft3d_info%bsendbuf,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%bsendbuf',ierr)
    end if

    if (allocated(fft3d_info%fsendbuf)) then
       deallocate(fft3d_info%fsendbuf,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%fsendbuf',ierr)
    end if

    if (allocated(fft3d_info%buf23slab)) then
       deallocate(fft3d_info%buf23slab,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%buf23slab',ierr)
    end if

    if (allocated(fft3d_info%buf12slab)) then
       deallocate(fft3d_info%buf12slab,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%buf12slab',ierr)
    end if

    if (allocated(fft3d_info%idx23slab)) then
       deallocate(fft3d_info%idx23slab,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%idx23slab',ierr)
    end if

    if (allocated(fft3d_info%idx12slab)) then
       deallocate(fft3d_info%idx12slab,stat=ierr)
       call utils_dealloc_check('internal_parallel_exit (fourier_mod.F90)', &
            'fft3d_info%idx12slab',ierr)
    end if

  end subroutine internal_parallel_exit

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================!
  ! The following routines are provided for debugging purposes !
  !============================================================!
#ifdef DEBUG
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_write_real_2(baseunit,darray)

    use comms, only: pub_my_node_id

    implicit none

    integer, intent(in) :: baseunit
    real(kind=DP), intent(in) :: darray(:,:)

    integer :: i1,i2

    write(baseunit+pub_my_node_id,'(i6)') 2
    write(baseunit+pub_my_node_id,'(2i6)') size(darray,1),size(darray,2)
    do i2=1,size(darray,2)
       do i1=1,size(darray,1),2
          write(baseunit+pub_my_node_id,'(2e24.16)') darray(i1,i2), &
               darray(i1+1,i2)
       end do
    end do

  end subroutine fourier_write_real_2

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_write_real_3(baseunit,darray)

    use comms, only: pub_my_node_id

    implicit none

    integer, intent(in) :: baseunit
    real(kind=DP), intent(in) :: darray(:,:,:)

    integer :: i1,i2,i3

    write(baseunit+pub_my_node_id,'(i6)') 3
    write(baseunit+pub_my_node_id,'(2i6)') size(darray,1),size(darray,2), &
         size(darray,3)
    do i3=1,size(darray,3)
       do i2=1,size(darray,2)
          do i1=1,size(darray,1),2
             write(baseunit+pub_my_node_id,'(2e24.16)') darray(i1,i2,i3), &
                  darray(i1+1,i2,i3)
          end do
       end do
    end do

  end subroutine fourier_write_real_3

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_write_complex_2(baseunit,zarray)

    use comms, only: pub_my_node_id

    implicit none

    integer, intent(in) :: baseunit
    complex(kind=DP), intent(in) :: zarray(:,:)

    integer :: i1,i2

    write(baseunit+pub_my_node_id,'(i6)') 2
    write(baseunit+pub_my_node_id,'(2i6)') 2*size(zarray,1),size(zarray,2)
    do i2=1,size(zarray,2)
       do i1=1,size(zarray,1)
          write(baseunit+pub_my_node_id,'(2e24.16)') real(zarray(i1,i2), &
               kind=DP),aimag(zarray(i1,i2))
       end do
    end do

  end subroutine fourier_write_complex_2

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_write_complex_3(baseunit,zarray)

    use comms, only: pub_my_node_id

    implicit none

    integer, intent(in) :: baseunit
    complex(kind=DP), intent(in) :: zarray(:,:,:)

    integer :: i1,i2,i3

    write(baseunit+pub_my_node_id,'(i6)') 3
    write(baseunit+pub_my_node_id,'(2i6)') 2*size(zarray,1),size(zarray,2), &
         size(zarray,3)
    do i3=1,size(zarray,3)
       do i2=1,size(zarray,2)
          do i1=1,size(zarray,1)
             write(baseunit+pub_my_node_id,'(2e24.16)') real(zarray(i1,i2,i3), &
                  kind=DP),aimag(zarray(i1,i2,i3))
          end do
       end do
    end do

  end subroutine fourier_write_complex_3
#endif

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  real(kind=DP) function calculate_flops(info)
    ! jd: Arguments
    type(parallel_fft3d_info), intent(in) :: info

    ! ------------------------------------------------------------------------
    calculate_flops = (2.5_DP * info%num12slabs * info%n2 * info%n1 * &
         log(real(info%n1,kind=DP)) + 5.0_DP * info%num23slabs * info%n2 * &
         info%n3 * log(real(info%n2*info%n3,kind=DP))) / log(2.0_DP)

  end function calculate_flops

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine make_license_heartbeat()

#ifdef ACCELRYS
    use license, only: license_sendHeartbeat,LIC_SUCCESS
    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    integer :: ierr

    ! Make a license heartbeat call
    call license_sendHeartbeat(ierr)
    if (ierr /= LIC_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i6)') &
            'Error in fourier_apply_cell_{forward|backward}: &
            &license_sendHeartbeat failed with code ',ierr
       call comms_abort
    end if
#endif

  end subroutine make_license_heartbeat
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module fourier
