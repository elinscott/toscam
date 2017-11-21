! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris and Arash A. Mostofi
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module kinetic

  implicit none

  private

  public :: kinetic_apply_on_box
  public :: kinetic_apply_to_func_batch
  public :: kinetic_apply_to_funcs_on_grid
  public :: kinetic_grad_on_box
  public :: kinetic_grad_to_func_batch

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_apply_to_func_batch(kinetic_fftbox_batch, &  ! output
       funcs_on_grid, fbasis, batch_size, local_start, local_end) ! input

    !======================================================================!
    ! This subroutine applies the kinetic energy operator to a batch       !
    ! of functions (eg NGWFs) in fftboxes.                                 !
    !======================================================================!
    ! Arguments:                                                           !
    ! kinetic_fftbox_batch (output): Batch of fftboxes with kinetic-NGWFs  !
    ! funcs_on_grid (input) : functions on this node in ppd representation !
    ! fbasis (input)        : Function basis type for this function set    !
    ! batch_size (input)    : Number of fftboxes in the batch              !
    ! local_start(input)    : First function of the batch                  !
    ! local_end  (input)    : Last function of the batch                   !
    !======================================================================!
    ! Written by Chris-Kriton Skylaris on 20/11/2003.                      !
    ! Modified to use new basis routines by Nicholas Hine, May 2008        !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type.  !
    ! Modified to not use workspace_mod by Nicholas Hine, November 2009.   !
    !======================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_copy_function_to_box
    use constants, only: DP
    use fourier, only: fourier_apply_box, fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    integer, intent(in) :: batch_size, local_start, local_end
    real(kind=DP), intent(in)  :: &
         funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(out) :: &
         kinetic_fftbox_batch(pub_fftbox%total_ld1, pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3, batch_size)

    ! Local Variables
    integer :: row1, row2
    integer :: iii
    integer :: ierr
    integer :: row_start1, row_start2, row_start3
    integer :: batch_count
    real(kind=DP), allocatable :: dwork_box1(:,:,:), dwork_box2(:,:,:)
    complex(kind=DP), allocatable :: zwork_box(:,:,:)


    ! ndmh: allocate workspace
    allocate(dwork_box1(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_apply_to_func_batch','dwork_box1',ierr)
    allocate(dwork_box2(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_apply_to_func_batch','dwork_box2',ierr)
    allocate(zwork_box(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_apply_to_func_batch','zwork_box',ierr)

    ! cks: determine where 'row' tightbox begins wrt fftbox
    call basis_ket_start_wrt_fftbox(row_start1, row_start2, row_start3, &
         pub_fftbox%total_pt1, pub_fftbox%total_pt2, pub_fftbox%total_pt3)

    batch_count =0
    do iii =local_start, local_end, 2
       row1 =iii ; row2 =iii +1

       ! ndmh: copy first function to FFTbox
       call basis_copy_function_to_box(dwork_box1, pub_fftbox%total_ld1, &
            pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            row_start1, row_start2, row_start3,&
            fbasis%tight_boxes(row1), funcs_on_grid, &
            fbasis%spheres(row1))

       if (row2 <= fbasis%node_num) then

          ! ndmh: copy second function to FFTbox
          call basis_copy_function_to_box(dwork_box2, pub_fftbox%total_ld1, &
               pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
               row_start1, row_start2, row_start3,&
               fbasis%tight_boxes(row2), funcs_on_grid, &
               fbasis%spheres(row2))

       else
          dwork_box2 = 0.0_DP
       end if

       ! FFT 'row' functions to reciprocal space:
       ! aam: use new fourier routines
       call fourier_apply_box_pair('Coarse', 'Forward', dwork_box1, &
            dwork_box2, zwork_box)

       ! Apply kinetic energy operator
       call kinetic_apply_on_box(zwork_box)

       ! FFT back to real space:
       call fourier_apply_box('Coarse', 'Backward', zwork_box)

       ! cks: put functions in batch of fftboxes after the kinetic
       ! cks: operator has been applied on them
       batch_count = batch_count + 1
       kinetic_fftbox_batch(:,:,:,batch_count) = real(zwork_box,kind=DP)
       if (row2 <= local_end) then
          batch_count = batch_count + 1
          kinetic_fftbox_batch(:,:,:,batch_count) = aimag(zwork_box)
       end if

    end do

    ! ndmh: deallocate workspace
    deallocate(zwork_box,stat=ierr)
    call utils_dealloc_check('kinetic_apply_to_func_batch','zwork_box',ierr)
    deallocate(dwork_box2,stat=ierr)
    call utils_dealloc_check('kinetic_apply_to_func_batch','dwork_box2',ierr)
    deallocate(dwork_box1,stat=ierr)
    call utils_dealloc_check('kinetic_apply_to_func_batch','dwork_box1',ierr)

  end subroutine kinetic_apply_to_func_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_apply_to_funcs_on_grid(kin_funcs_on_grid, &  ! output
       funcs_on_grid, fbasis)                                     ! input

    !======================================================================!
    ! This subroutine applies the kinetic energy operator to a batch       !
    ! of functions (eg NGWFs) in fftboxes.                                 !
    !======================================================================!
    ! Arguments:                                                           !
    ! kin_funcs_on_grid (output): kinetic energy operator applied to funcs !
    ! funcs_on_grid (input) : functions on this node in ppd representation !
    ! fbasis (input)        : Function basis type for this function set    !
    ! batch_size (input)    : Number of fftboxes in the batch              !
    ! local_start(input)    : First function of the batch                  !
    ! local_end  (input)    : Last function of the batch                   !
    !======================================================================!
    ! Written by Nicholas Hine, March 2011.                                !
    !======================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_copy_function_to_box, &
         basis_extract_function_from_box
    use constants, only: DP
    use fourier, only: fourier_apply_box, fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(out)  :: &
         kin_funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in)  :: &
         funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    integer :: row1, row2
    integer :: iii
    integer :: ierr
    integer :: row_start1, row_start2, row_start3
    real(kind=DP), allocatable :: dwork_box1(:,:,:), dwork_box2(:,:,:)
    complex(kind=DP), allocatable :: zwork_box(:,:,:)


    ! ndmh: allocate workspace
    allocate(dwork_box1(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_apply_to_funcs_on_grid','dwork_box1',ierr)
    allocate(dwork_box2(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_apply_to_funcs_on_grid','dwork_box2',ierr)
    allocate(zwork_box(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_apply_to_funcs_on_grid','zwork_box',ierr)

    ! cks: determine where 'row' tightbox begins wrt fftbox
    call basis_ket_start_wrt_fftbox(row_start1, row_start2, row_start3, &
         pub_fftbox%total_pt1, pub_fftbox%total_pt2, pub_fftbox%total_pt3)

    do iii=1,fbasis%node_num, 2
       row1 = iii
       row2 = iii + 1

       ! ndmh: copy first function to FFTbox
       call basis_copy_function_to_box(dwork_box1, pub_fftbox%total_ld1, &
            pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            row_start1, row_start2, row_start3,&
            fbasis%tight_boxes(row1), funcs_on_grid, &
            fbasis%spheres(row1))

       if (row2 <= fbasis%node_num) then

          ! ndmh: copy second function to FFTbox
          call basis_copy_function_to_box(dwork_box2, pub_fftbox%total_ld1, &
               pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
               row_start1, row_start2, row_start3,&
               fbasis%tight_boxes(row2), funcs_on_grid, &
               fbasis%spheres(row2))

       else
          dwork_box2 = 0.0_DP
       end if

       ! FFT 'row' functions to reciprocal space:
       ! aam: use new fourier routines
       call fourier_apply_box_pair('Coarse', 'Forward', dwork_box1, &
            dwork_box2, zwork_box)

       ! Apply kinetic energy operator
       call kinetic_apply_on_box(zwork_box)

       ! FFT back to real space:
       call fourier_apply_box_pair('Coarse', 'Backward', dwork_box1, &
            dwork_box2, zwork_box)

       ! ndmh: extract first function from FFTbox
       call basis_extract_function_from_box(kin_funcs_on_grid, &
            pub_fftbox%total_ld1, pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            dwork_box1, fbasis%spheres(row1), fbasis%tight_boxes(row1), &
            row_start1, row_start2, row_start3, fbasis%spheres(row1)%offset)

       if (row2 <= fbasis%node_num) then

          ! ndmh: extract second function from FFTbox
          call basis_extract_function_from_box(kin_funcs_on_grid, &
               pub_fftbox%total_ld1, pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
               dwork_box2, fbasis%spheres(row2), fbasis%tight_boxes(row2), &
               row_start1, row_start2, row_start3, fbasis%spheres(row2)%offset)

       else
          dwork_box2 = 0.0_DP
       end if

    end do

    ! ndmh: deallocate workspace
    deallocate(zwork_box,stat=ierr)
    call utils_dealloc_check('kinetic_apply_to_funcs_on_grid','zwork_box',ierr)
    deallocate(dwork_box2,stat=ierr)
    call utils_dealloc_check('kinetic_apply_to_funcs_on_grid','dwork_box2',ierr)
    deallocate(dwork_box1,stat=ierr)
    call utils_dealloc_check('kinetic_apply_to_funcs_on_grid','dwork_box1',ierr)

  end subroutine kinetic_apply_to_funcs_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_apply_on_box(data_complex)

    !=========================================================!
    ! This subroutine applies the kinetic energy operator in  !
    ! reciprocal space to a function in an FFTbox.            !
    !---------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris in 2000.    !
    ! Rewritten by Peter D. Haynes                            !
    !=========================================================!

    use constants, only: DP
    use simulation_cell, only: pub_fftbox

    implicit none

    ! Argument
    complex(kind=DP), intent(inout) :: data_complex(:,:,:)

    ! Local variable
    integer :: i1,i2,i3

    ! apply the kinetic energy operator
    do i3=1,pub_fftbox%total_pt3
       do i2=1,pub_fftbox%total_pt2
          do i1=1,pub_fftbox%total_pt1
             data_complex(i1,i2,i3) = data_complex(i1,i2,i3) * &
                  pub_fftbox%recip_grid(5,i1,i2,i3)
          end do
       end do
    end do

  end subroutine kinetic_apply_on_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_grad_to_func_batch(grad_fftbox_batch, &      ! output
       funcs_on_grid, fbasis, batch_size, local_start, local_end) ! input

    !=========================================================================!
    ! This subroutine places a batch of functions (eg NGWFs) in FFTboxes and  !
    ! applies the grad operator to them, by transforming to reciprocal space, !
    ! multiplying by iG and transforming back.                                !
    !=========================================================================!
    ! Arguments:                                                              !
    ! grad_fftbox_batch (output): Batch of fftboxes                           !
    ! funcs_on_grid (input) : functions on this node in ppd representation    !
    ! fbasis (input)        : Function basis type this function set           !
    ! batch_size (input)    : Number of fftboxes in the batch                 !
    ! local_start(input)    : First function of the batch                     !
    ! local_end  (input)    : Last function of the batch                      !
    !=========================================================================!
    ! Written by Peter Haynes on 16/7/2005.                                   !
    ! Based on kinetic_apply_to_ngwf_batch by Chris-Kriton Skylaris on        !
    !    20/11/2003.                                                          !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type.     !
    !=========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_copy_function_to_box
    use constants, only: DP
    use fourier, only: fourier_apply_box, fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    integer, intent(in) :: batch_size, local_start, local_end
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(in) :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(out) :: &
         grad_fftbox_batch(pub_fftbox%total_ld1, pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3, batch_size, 3)

    ! <<< local variables >>>
    integer :: row1, row2
    integer :: iii
    integer :: row_start1, row_start2, row_start3
    integer :: batch_count
    integer :: ierr
    integer :: dim
    real(kind=DP), allocatable :: dwork_box1(:,:,:), dwork_box2(:,:,:)
    complex(kind=DP), dimension(:,:,:), allocatable :: zwork_box
    complex(kind=DP), dimension(:,:,:,:), allocatable :: zwork_grad_box

    ! ndmh: allocate workspace
    allocate(dwork_box1(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_grad_to_func_batch','dwork_box1',ierr)
    allocate(dwork_box2(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_grad_to_func_batch','dwork_box2',ierr)
    allocate(zwork_box(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('kinetic_grad_to_func_batch','zwork_box',ierr)
    allocate(zwork_grad_box(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3,3),stat=ierr)
    call utils_alloc_check('kinetic_grad_to_func_batch','zwork_grad_box',ierr)

    zwork_box = (0.0_DP,0.0_DP)

    ! Position of 'row' function in fftbox
    ! cks: determine where 'row' tightbox begins wrt fftbox
    call basis_ket_start_wrt_fftbox(row_start1, row_start2, row_start3, &
        pub_fftbox%total_pt1, pub_fftbox%total_pt2, pub_fftbox%total_pt3)

    batch_count = 0
    do iii=local_start, local_end, 2
       row1=iii ; row2=iii+1

       ! ndmh: use new basis_copy routine
       call basis_copy_function_to_box(dwork_box1, pub_fftbox%total_ld1, &
            pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            row_start1, row_start2, row_start3, &
            fbasis%tight_boxes(row1), funcs_on_grid, &
            fbasis%spheres(row1))

       if (row2 <= fbasis%node_num) then
          ! ndmh: use new basis_copy routine
          call basis_copy_function_to_box(dwork_box2, pub_fftbox%total_ld1, &
               pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
               row_start1, row_start2, row_start3, &
               fbasis%tight_boxes(row2), funcs_on_grid, &
               fbasis%spheres(row2))
       else
          dwork_box2 = 0.0_DP
       end if

       ! FFT 'row' functions to reciprocal space:
       ! aam: use new fourier routines
       call fourier_apply_box_pair('Coarse', 'Forward', dwork_box1, &
            dwork_box2, zwork_box)

       call kinetic_grad_on_box(zwork_box,zwork_grad_box)

       do dim=1,3
          ! FFT back to real space:
          call fourier_apply_box('Coarse','Backward',zwork_grad_box(:,:,:,dim))
       end do

       ! cks: put functions in batch of fftboxes after the grad
       ! cks: operator has been applied on them
       batch_count = batch_count + 1
       grad_fftbox_batch(:,:,:,batch_count,:) = real(zwork_grad_box,kind=DP)
       if (row2 <= local_end) then
          batch_count = batch_count + 1
          grad_fftbox_batch(:,:,:,batch_count,:) = aimag(zwork_grad_box)
       end if

    end do

    ! ndmh: deallocate workspace
    deallocate(zwork_grad_box,stat=ierr)
    call utils_dealloc_check('kinetic_grad_to_func_batch','zwork_grad_box',ierr)
    deallocate(zwork_box,stat=ierr)
    call utils_dealloc_check('kinetic_grad_to_func_batch','zwork_box',ierr)
    deallocate(dwork_box2,stat=ierr)
    call utils_dealloc_check('kinetic_grad_to_func_batch','dwork_box2',ierr)
    deallocate(dwork_box1,stat=ierr)
    call utils_dealloc_check('kinetic_grad_to_func_batch','dwork_box1',ierr)

  end subroutine kinetic_grad_to_func_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine kinetic_grad_on_box(data_complex_in,data_complex_out)

    !=========================================================!
    ! This subroutine applies the grad operator in            !
    ! reciprocal space to a function in an FFTbox.            !
    !---------------------------------------------------------!
    ! Written by Peter Haynes on 16/7/2005.                   !
    ! Based on kinetic_apply_on_box by Chris-Kriton Skylaris  !
    ! in 2000.                                                !
    !=========================================================!

    use constants, only: DP
    use simulation_cell, only: pub_fftbox

    implicit none

    ! Argument
    complex(kind=DP), intent(in) :: data_complex_in(:,:,:)
    complex(kind=DP), intent(out) :: data_complex_out(:,:,:,:)

    ! Local variable
    integer :: i1,i2,i3,dim

    ! apply the grad operator
    do i3=1,pub_fftbox%total_pt3
       do i2=1,pub_fftbox%total_pt2
          do i1=1,pub_fftbox%total_pt1
             do dim=1,3
                data_complex_out(i1,i2,i3,dim) = data_complex_in(i1,i2,i3) * &
                     cmplx(0.0_DP,pub_fftbox%recip_grid(dim,i1,i2,i3),kind=DP)
             end do
          end do
       end do
    end do

  end subroutine kinetic_grad_on_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module kinetic

