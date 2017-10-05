! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!   with subsequent optimisations by Nicholas D.M. Hine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module potential

  implicit none

  private

  public :: potential_apply_to_ngwf_batch
  public :: potential_apply_to_ppd_funcs
  public :: potential_sawtooth_efield
  public :: potential_add_efield_ion_energy
  public :: potential_input_to_workspace

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine potential_sawtooth_efield(locpot_on_grid,grid)

    !==========================================================================!
    ! This subroutine adds to the local potential a "sawtooth" contribution    !
    ! due to a uniform external electric field. The electric field is given    !
    ! in the input file in terms of its cartesian coordinates in atomic units. !
    ! 1.0 volt/angstrom = 0.01944655 Eh/(e*a0).                                !
    !==========================================================================!
    ! Arguments:                                                               !
    ! locpot_on_grid (input-output): local potential in simulation cell file   !
    !                             to which the electric field contribution     !
    !                             is to be added.                              !
    !==========================================================================!
    ! Written by Chris-Kriton Skylaris on 9/5/2004 for the ONETEP              !
    ! linear-scaling DFT program.                                              !
    ! Rewritten by Peter Haynes 1/7/2004 for fourier parallelisation.          !
    ! Modified 22/01/2011 by Nicholas Hine to use cell_grid_real_pt routine to !
    ! save on storage.                                                         !
    !==========================================================================!

    use cell_grid,         only: GRID_INFO, cell_grid_real_pt
    use constants,         only: DP
    use rundat,            only: pub_constant_efield
    use simulation_cell,   only: pub_cell

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    real(kind=DP), intent(inout) :: locpot_on_grid(grid%ld1,grid%ld2, &
         grid%max_slabs12)


    ! Local variables
    integer :: i1,i2,islab12      ! Loop counters
    real(kind=DP) :: rpt(3)

    ! Loop over real-space grid on this node and add contribution to locpot
    ! at each point
    do islab12=1,grid%num_my_slabs12
       do i2=1,grid%n2
          do i1=1,grid%n1
             call cell_grid_real_pt(rpt,i1,i2,islab12,grid)
             locpot_on_grid(i1,i2,islab12) = locpot_on_grid(i1,i2,islab12) + &
                  sum(pub_constant_efield(1:3)*rpt(1:3))
          end do
       end do
    end do

  end subroutine potential_sawtooth_efield


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine potential_add_efield_ion_energy(ewald_energy)

    !==========================================================================!
    ! This subroutine adds to the ion-ion energy the contribution from the     !
    ! interaction of a uniform external electric field with the ion charges.   !
    ! The electric field is given in the input file in terms of its cartesian  !
    ! coordinates in atomic units. 1.0 volt/angstrom = 0.01944655 Eh/(e*a0).   !
    !==========================================================================!
    ! Arguments:                                                               !
    ! ewald_energy (inout) : ion-ion energy to which to add efield-ion energy. !
    !==========================================================================!
    ! Written by Nicholas Hine on 10/06/2011.                                  !
    !==========================================================================!

    use comms,             only: comms_reduce, pub_my_node_id
    use constants,         only: DP
    use parallel_strategy, only: pub_num_atoms_on_node, pub_elements_on_node
    use rundat,            only: pub_constant_efield

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: ewald_energy

    ! Local variables
    integer :: iat
    real(kind=DP) :: efield_energy

    ! Loop over atoms on this node and calculate energy of interaction
    ! for each
    efield_energy = 0.0_DP
    do iat=1,pub_num_atoms_on_node(pub_my_node_id)
       efield_energy = efield_energy - pub_elements_on_node(iat)%ion_charge &
            *(pub_constant_efield(1)*pub_elements_on_node(iat)%centre%x &
            + pub_constant_efield(2)*pub_elements_on_node(iat)%centre%y &
            + pub_constant_efield(3)*pub_elements_on_node(iat)%centre%z)
    end do
    
    call comms_reduce('SUM',efield_energy)

    ewald_energy = ewald_energy + efield_energy

  end subroutine potential_add_efield_ion_energy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine potential_apply_to_ngwf_batch(potential_fftbox_batch, &  ! output
       ngwfs_on_grid, ngwf_basis, potential_dbl, &                    ! input
       batch_size, local_start, local_end, max_current_size)          ! input

    !==========================================================================!
    ! This subroutine applies the local potential operator to a batch          !
    ! of NGWFs in fftboxes.                                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! potential_fftbox_batch (output): Batch off fftboxes with potential-NGWFs !
    ! ngwfs_on_grid (input) : NGWFs on this node in ppd representation         !
    ! ngwf_basis (input)    : Function basis type for NGWFs                    !
    ! potential_dbl (input) : Local potential on dbl grid in simulation cell   !
    ! batch_size (input)    : Number of fftboxes in the batch                  !
    ! local_start(input)    : First NGWF of the batch                          !
    ! local_end  (input)    : Last NGWF of the batch                           !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/11/2003 for the ONETEP program.   !
    ! Modified by Chris-Kriton Skylaris on 10/7/2004 so that                   !
    ! it works with the data-parallel local potential on the fine grid.        !
    ! Modified by Chris-Kriton Skylaris on 16/11/2004 to fix bug that could    !
    ! cause trouble when pub_cell%node_num=1.                                  !
    ! Modified by Nicholas Hine in 2008 to replace two-stage copy of each row  !
    ! function with single call to basis_copy_function_to_box.                 !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type and   !
    ! to deposit result of fourier_filter straight into potential_fftbox array.!
    !==========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, &
         basis_copy_function_to_box, basis_location_func_wrt_cell
    use cell_grid, only: cell_grid_extract_box, pub_dbl_grid
    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use fourier, only: fourier_filter, fourier_interpolate
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell, pub_fftbox
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size, local_start, local_end
    integer, intent(in) :: max_current_size ! max size of batch over all nodes
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(out) :: potential_fftbox_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3, batch_size)
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in) :: potential_dbl(pub_dbl_grid%ld1, &
         pub_dbl_grid%ld2, pub_dbl_grid%max_group_slabs12)

    ! Local Variables
    real(kind=DP), dimension(:,:,:), allocatable :: row1_box, row2_box
    real(kind=DP), dimension(:,:,:), allocatable :: potential_box_dbl
    real(kind=DP), dimension(:,:,:), allocatable :: row1_box_dbl
    real(kind=DP), dimension(:,:,:), allocatable :: row2_box_dbl
    real(kind=DP), dimension(:,:,:), allocatable :: buffer_dbl

    integer :: ierr
    integer :: row1, row2
    integer :: n1, n2, n3, ld1, ld2
    integer :: n1_dbl, n2_dbl, n3_dbl, ld1_dbl, ld2_dbl
    integer :: iii
    integer :: row_start1, row_start2, row_start3
    integer :: batch_count
    integer :: row1_cell_start1, row1_cell_start2, row1_cell_start3
    integer :: row2_cell_start1, row2_cell_start2, row2_cell_start3
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: prev_start1, prev_start2, prev_start3

    logical :: i_need_potential ! whether pub_my_node_id needs pot'l in fftbox

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering potential_apply_to_ngwf_batch'
#endif

    ! Start timer
    call timer_clock('potential_apply_to_ngwf_batch', 1)

    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2
    n1_dbl = pub_fftbox%total_pt1_dbl
    n2_dbl = pub_fftbox%total_pt2_dbl
    n3_dbl = pub_fftbox%total_pt3_dbl
    ld1_dbl = pub_fftbox%total_ld1_dbl
    ld2_dbl = pub_fftbox%total_ld2_dbl

    ! ndmh: initialisations to prevent compiler warnings
    prev_start1 = -1111111
    prev_start2 = -2222222
    prev_start3 = -3333333
    row2_cell_start1 = -1111111
    row2_cell_start2 = -2222222
    row2_cell_start3 = -3333333

    ! Allocate workspace arrays
    allocate(row1_box(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch','row1_box',ierr)
    allocate(row2_box(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch','row2_box',ierr)
    allocate(potential_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch',&
         'potential_box_dbl',ierr)
    allocate(row1_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch','row1_box_dbl',&
         ierr)
    allocate(row2_box_dbl(ld1_dbl, ld2_dbl, 2*n3), stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch','row2_box_dbl',&
         ierr)
    allocate(buffer_dbl(ld1_dbl, ld2_dbl, pub_dbl_grid%max_group_slabs12), &
         stat=ierr)
    call utils_alloc_check('potential_apply_to_ngwf_batch','buffer_dbl',&
         ierr)

    ! cks: initialisation
    buffer_dbl = 0.0_DP
    row1_box_dbl = 0.0_DP
    row2_box_dbl = 0.0_DP
    potential_box_dbl = 0.0_DP

    ! cks: determine where 'row' tightbox begins wrt fftbox
    call basis_ket_start_wrt_fftbox(row_start1, row_start2, row_start3, &
         n1, n2, n3)

    batch_count = 0
    do iii=local_start, local_start+max_current_size-1, 2
       row1 = iii ; row2 = iii + 1

       if (row1 <= local_end) then

          ! ndmh: use new basis function copying straight to fftbox from ppds
          call basis_copy_function_to_box(row1_box, ld1, ld2, n3, &
               row_start1, row_start2, row_start3,&
               ngwf_basis%tight_boxes(row1), ngwfs_on_grid, &
               ngwf_basis%spheres(row1))

          ! Find position of 'row1' functions wrt simulation cell:
          call basis_location_func_wrt_cell(row1_cell_start1, &
               row1_cell_start2,row1_cell_start3,ngwf_basis%tight_boxes(row1))
       else
          row1_box = 0.0_DP

          ! ndmh: prevent unassigned use of row1_cell_start's in debug mode
          row1_cell_start1 = -1234
          row1_cell_start2 = -1234
          row1_cell_start3 = -1234
       end if

       if (row2 <= local_end) then

          ! ndmh: use new basis function copying straight to fftbox from ppds
          call basis_copy_function_to_box(row2_box, ld1, ld2, n3, &
               row_start1, row_start2, row_start3,&
               ngwf_basis%tight_boxes(row2), ngwfs_on_grid, &
               ngwf_basis%spheres(row2))

          ! Find position of 'row2' functions wrt simulation cell:
          call basis_location_func_wrt_cell(row2_cell_start1, &
               row2_cell_start2,row2_cell_start3,ngwf_basis%tight_boxes(row2))

       else
          row2_box = 0.0_DP

          ! ndmh: prevent unassigned use of row2_cell_start's in debug mode
          row2_cell_start1 = -1234
          row2_cell_start2 = -1234
          row2_cell_start3 = -1234
       end if


       if (row1 <= local_end) then
          ! cks: interpolate row1 and row2 ngwfs
          call fourier_interpolate(row1_box, row2_box, row1_box_dbl, &
               row2_box_dbl)
       else
          row1_box_dbl = 0.0_DP
          row2_box_dbl = 0.0_DP
       end if


       ! cks:----- Get potential for FIRST NGWF if needed ------------
       i_need_potential = ((row1_cell_start1 /= prev_start1 .or. &
            row1_cell_start2 /= prev_start2 .or. &
            row1_cell_start3 /= prev_start3) .and. &
            row1 <= local_end)

       ! cks: Put potential into FFTbox for first NGWF
       ! ndmh: routine moved to basis_mod
       fftbox_start1_dbl = 2*(row1_cell_start1 - row_start1) + 1
       fftbox_start2_dbl = 2*(row1_cell_start2 - row_start2) + 1
       fftbox_start3_dbl = 2*(row1_cell_start3 - row_start3) + 1

       call cell_grid_extract_box(potential_box_dbl,&
            buffer_dbl, potential_dbl, pub_dbl_grid, n1_dbl, n2_dbl, n3_dbl, &
            ld1_dbl, ld2_dbl, fftbox_start1_dbl, fftbox_start2_dbl, &
            fftbox_start3_dbl, i_need_potential, .true.)

       if (i_need_potential) then
          prev_start1 = row1_cell_start1
          prev_start2 = row1_cell_start2
          prev_start3 = row1_cell_start3
       end if
       ! cks:-- END Get potential for FIRST NGWF if needed ------------

       ! Apply potential to 'row1' function:
       if (row1 <= local_end) row1_box_dbl = potential_box_dbl * row1_box_dbl


       ! cks:----- Get potential for SECOND NGWF if needed -----------
       i_need_potential = ((row2_cell_start1 /= prev_start1 .or. &
            row2_cell_start2 /= prev_start2 .or. &
            row2_cell_start3 /= prev_start3) .and. &
            row2 <= local_end)

       ! cks: Put potential into FFTbox for second NGWF
       fftbox_start1_dbl = 2*(row2_cell_start1 - row_start1) + 1
       fftbox_start2_dbl = 2*(row2_cell_start2 - row_start2) + 1
       fftbox_start3_dbl = 2*(row2_cell_start3 - row_start3) + 1
       call cell_grid_extract_box(potential_box_dbl, &
            buffer_dbl, potential_dbl, pub_dbl_grid, n1_dbl, n2_dbl, n3_dbl, &
            ld1_dbl, ld2_dbl, fftbox_start1_dbl, fftbox_start2_dbl, &
            fftbox_start3_dbl, i_need_potential, .true.)

       if (i_need_potential) then
          prev_start1 = row2_cell_start1
          prev_start2 = row2_cell_start2
          prev_start3 = row2_cell_start3
       endif
       ! cks:-- END Get potential for SECOND NGWF if needed -----------

       ! Apply potential to 'row2' function:
       if (row2 <= local_end) row2_box_dbl = potential_box_dbl * row2_box_dbl

       ! ndmh: deposit one result of interpolation in potential_fftbox_batch
       ! ndmh: and discard other result
       if ((row1 <= local_end).and.(row2 > local_end)) then

          batch_count = batch_count + 1

          ! Filter result to standard grid
          call fourier_filter(row1_box_dbl, row2_box_dbl, &
               potential_fftbox_batch(:,:,:,batch_count), row2_box)

       ! ndmh: deposit both results of interpolation in potential_fftbox_batch
       else if ((row1 <= local_end).and.(row2 <= local_end)) then

          batch_count = batch_count + 2

          ! Filter result to standard grid
          call fourier_filter(row1_box_dbl, row2_box_dbl, &
               potential_fftbox_batch(:,:,:,batch_count-1), &
               potential_fftbox_batch(:,:,:,batch_count))

       end if

    end do

    deallocate(buffer_dbl, stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch','buffer_dbl',&
         ierr)
    deallocate(row2_box_dbl, stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch','row2_box_dbl',&
         ierr)
    deallocate(row1_box_dbl, stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch','row1_box_dbl',&
         ierr)
    deallocate(potential_box_dbl, stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch',&
         'potential_box_dbl',ierr)
    deallocate(row2_box, stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch','row2_box',ierr)
    deallocate(row1_box, stat=ierr)
    call utils_dealloc_check('potential_apply_to_ngwf_batch','row1_box',ierr)

    call timer_clock('potential_apply_to_ngwf_batch', 2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving potential_apply_to_ngwf_batch'
#endif

  end subroutine potential_apply_to_ngwf_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine potential_apply_to_ppd_funcs(pot_funcs_on_grid, & ! in-output
       fbasis, potential_std)                                  ! input

    !==========================================================================!
    ! This subroutine calculates applies the potential to a set of functions   !
    ! in PPD storage.                                                          !
    !==========================================================================!
    !  Arguments:                                                              !
    !    pot_funcs_on_grid (inout) : functions in PPD format to multiply by    !
    !                                potential                                 !
    !    fbasis (input)      : The function basis type for the functions       !
    !    potential_dbl (in)  : local potential on the "dbl" sim cell grid,     !
    !                          which is actually coarse as dbl_grid_scale=1    !
    !==========================================================================!
    ! Written by Nicholas Hine on 28/02/2011.                                  !
    !==========================================================================!

    use basis, only: basis_location_func_wrt_cell, &
         basis_multiply_function_by_box
    use cell_grid, only: cell_grid_extract_box, pub_std_grid
    use comms, only: comms_barrier, pub_on_root, pub_my_node_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(inout) :: pot_funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in) :: potential_std(pub_std_grid%ld1, &
         pub_std_grid%ld2, pub_std_grid%max_group_slabs12)

    ! Local Variables
    integer :: ierr
    integer :: local_func
    integer :: cell_start1,cell_start2,cell_start3
    integer :: box_n1,box_n2,box_n3
    real(kind=DP), allocatable :: potential_box(:,:,:)
    real(kind=DP), allocatable :: potential_buffer(:,:,:)
    logical :: i_need_box

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &potential_apply_to_ppd_funcs'
#endif

    ! Start timer
    call timer_clock('potential_apply_to_ppd_funcs',1)

    box_n1 = fbasis%maxtight_pts1
    box_n2 = fbasis%maxtight_pts2
    box_n3 = fbasis%maxtight_pts3

    ! Allocate storage for tightbox of potential and buffer
    allocate(potential_box(box_n1,box_n2,box_n3),stat=ierr)
    call utils_alloc_check('potential_apply_to_ppd_funcs','potential_box',ierr)
    allocate(potential_buffer(box_n1,box_n2,pub_std_grid%max_group_slabs12), &
         stat=ierr)
    call utils_alloc_check('potential_apply_to_ppd_funcs','potential_buffer', &
         ierr)

    ! Loop over ket functions up to max on any node
    do local_func=1,fbasis%max_on_node

       if (local_func <= fbasis%num_on_node(pub_my_node_id)) then
          ! Find position of ket function wrt simulation cell:
          call basis_location_func_wrt_cell(cell_start1, &
               cell_start2,cell_start3,fbasis%tight_boxes(local_func))
          i_need_box = .true.
       else
          cell_start1 = -1234
          cell_start2 = -1234
          cell_start3 = -1234
          i_need_box = .false.
       end if

       ! Extract tightbox of locpot data from whole-cell grid
       call cell_grid_extract_box(potential_box,&
            potential_buffer, potential_std, pub_std_grid, &
            box_n1, box_n2, box_n3, box_n1, box_n2, &
            cell_start1, cell_start2, cell_start3, i_need_box, .true.)

       ! Multiply the ket function by the locpot box directly in PPDs
       if (local_func <= fbasis%num_on_node(pub_my_node_id)) then
          call basis_multiply_function_by_box(pot_funcs_on_grid, &
               box_n1, box_n2, box_n3, potential_box, &
               fbasis%spheres(local_func), fbasis%tight_boxes(local_func), &
               1, 1, 1, fbasis%spheres(local_func)%offset)
       end if

    end do

    ! Deallocate workspace
    deallocate(potential_buffer,stat=ierr)
    call utils_dealloc_check('potential_apply_to_ppd_funcs', &
         'potential_buffer',ierr)
    deallocate(potential_box,stat=ierr)
    call utils_dealloc_check('potential_apply_to_ppd_funcs', &
         'potential_box',ierr)

    ! pdh: re-sync nodes
    call comms_barrier

    ! Stop timer for total
    call timer_clock('potential_apply_to_ppd_funcs',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &potential_apply_to_ppd_funcs'
#endif

  end subroutine potential_apply_to_ppd_funcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine potential_input_to_workspace(potential_work,potential_in, &
       work_grid,in_grid)

    !==========================================================================!
    ! This subroutine transfers a potential from the input grid to a version   !
    ! on the 'workspace' grid, where all the slabs of the local node's comms   !
    ! group are duplicated on each node of the group. The grids may be of      !
    ! different scales, in which case the result is filtered from the input    !
    ! down to the workspace grid.                                              !
    !==========================================================================!
    !  Arguments:                                                              !
    !    potential_work (out) : Workspace array (sized for one group's slabs)  !
    !    potential_in   (in)  : Input array (sized for one node's slabs)       !
    !    work_grid      (in)  : GRID_INFO defining workspace grid              !
    !    in_grid        (in)  : GRID_INFO defining input grid                  !
    !==========================================================================!
    ! Written by Nicholas Hine on 16/03/2011.                                  !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root, comms_allgather, pub_comms_group_size, &
         pub_group_comm, pub_first_node_in_group, pub_my_node_id, &
         pub_my_rank_in_group
    use constants, only: DP
    use fourier, only: fourier_filter_cell
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: work_grid
    type(GRID_INFO), intent(in) :: in_grid
    real(kind=DP), intent(out) :: potential_work(work_grid%ld1,work_grid%ld2, &
         work_grid%max_group_slabs12)
    real(kind=DP), intent(in) :: potential_in(in_grid%ld1,in_grid%ld2, &
         in_grid%max_slabs12)

    ! Local Variables
    integer :: i3_start, i3_finish
    logical :: in_is_work

    in_is_work = .false.
    if ((work_grid%n1==in_grid%n1).and.(work_grid%n2==in_grid%n2).and. &
         (work_grid%n3==in_grid%n3)) in_is_work = .true.

    ! ndmh: get data from potential_fine to potential_dbl, depending on scales
    i3_start = work_grid%my_first_slab12_in_group
    i3_finish = work_grid%my_last_slab12_in_group
    if (in_is_work) then
       ! ndmh: copy potential from fine grid to double grid and allgather the
       ! ndmh: other nodes' data within the group if required
       if (pub_comms_group_size>1) then
          call comms_allgather(potential_work,potential_in, &
               length_src=work_grid%group_lengths_slabs12(pub_my_rank_in_group), &
               lengths_dest=work_grid%group_lengths_slabs12, &
               displs_dest=work_grid%group_displs_slabs12,comm=pub_group_comm)
       else
          potential_work(:,:,i3_start:i3_finish) = &
               potential_in(:,:,1:work_grid%num_my_slabs12)
       end if
    else
       ! ndmh: filter potential from fine grid to double grid and allgather the
       ! ndmh: other nodes' data within the group if required
       i3_finish = i3_start + work_grid%max_slabs12 - 1
       call fourier_filter_cell(potential_in(:,:,:), &
            potential_work(:,:,i3_start:i3_finish),in_grid,work_grid, &
            apply_nyquist=.true.)
       if (pub_comms_group_size>1) then
          call comms_allgather(potential_work,potential_work(:,:,i3_start:i3_finish), &
               length_src=work_grid%group_lengths_slabs12(pub_my_rank_in_group), &
               lengths_dest=work_grid%group_lengths_slabs12, &
               displs_dest=work_grid%group_displs_slabs12,comm=pub_group_comm)
       end if
    end if

  end subroutine potential_input_to_workspace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module potential
