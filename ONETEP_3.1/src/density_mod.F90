! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Arash A. Mostofi and Nicholas D.M. Hine
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module density

  use constants, only: DP

  implicit none

  private

  type RADIAL_DENSITY_TYPE

    ! ndmh: flag to show whether each initial guess density exists
    logical :: present

    ! ndmh: number of points of radial grid
    integer :: npts

    ! ndmh: radial grid positions
    real(kind=DP), allocatable, dimension(:) :: rad

    ! ndmh: density on radial grid
    real(kind=DP), allocatable, dimension(:,:) :: den

    ! ndmh: augmentation charge on radial grid (if present)
    real(kind=DP), allocatable, dimension(:,:) :: aug_den

  end type RADIAL_DENSITY_TYPE

  type(RADIAL_DENSITY_TYPE), allocatable :: radial_densities(:)

  public :: density_on_grid
  public :: density_initial_guess_recip
  public :: density_initial_guess_real
  public :: density_check_possemidefinite
  public :: density_render_possemidefinite
  public :: density_on_grid_renorm
  public :: density_plot_slice
  public :: density_radial_init
  public :: density_radial_store
  public :: density_radial_exit

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_on_grid(density_fine, &                         ! output
       grid, denskern, overlap, ngwfs_on_grid, ngwf_basis)           ! input

    !==========================================================================!
    ! This subroutine calculates the charge density of a fine grid in the      !
    ! simulation cell. Wrapper for calculating the density on the fine grid    !
    ! (which may be any scale), by calling the appropriate routine to generate !
    ! the density on the coarse or double grids and then interpolating up (if  !
    ! necessary) to the grid scale of the fine grid.                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! density_fine (output) : The total charge density of the system on the    !
    ! fine grid of the simulation cell                                         !
    ! denskern      (input) : Density kernel sparse matrix                     !
    ! overlap       (input) : Overlap matrix                                   !
    ! ngwfs_on_grid (input) : The NGWFs on the real-space ppd representation   !
    !  that belong to this node                                                !
    ! ngwf_basis            : The function basis for the NGWFs                 !
    !--------------------------------------------------------------------------!
    ! Current form by Nicholas Hine 02/28/11. Old code moved to new routine    !
    ! density_on_dbl_grid (see below).                                         !
    !==========================================================================!

    use cell_grid, only: GRID_INFO, pub_dbl_grid, pub_std_grid
    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_dbl_is_std
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins,1)
    type(SPAM3), intent(in) :: overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out) :: density_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    integer :: ierr
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_work

    ! ndmh: if the output is to be on the standard grid, use that as workspace
    if (((grid%n1==pub_std_grid%n1).and.(grid%n2==pub_std_grid%n2).and. &
         (grid%n3==pub_std_grid%n3)).or.(pub_dbl_is_std)) then

       ! ndmh: allocate coarse-grid whole cell workspace array
       allocate(density_work(pub_std_grid%ld1,pub_std_grid%ld2,&
            pub_std_grid%max_group_slabs12,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('density_on_grid','density_work',ierr)

       ! ndmh: calculate the density on the coarse grid
       call density_on_coarse_grid(density_work, &              ! output
            denskern, overlap, ngwfs_on_grid, ngwf_basis)       ! input

       ! ndmh: copy or upscale to output grid as required
       call density_workspace_to_output(density_work,density_fine, &
            pub_std_grid,grid)

       deallocate(density_work,stat=ierr)
       call utils_dealloc_check('density_on_grid','density_work',ierr)

    else ! ndmh: otherwise, default is to use dbl grid as workspace

       ! ndmh: allocate double-grid whole cell workspace array
       allocate(density_work(pub_dbl_grid%ld1,pub_dbl_grid%ld2,&
            pub_dbl_grid%max_group_slabs12,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('density_on_grid','density_work',ierr)

       ! ndmh: calculate the density on the double grid
       call density_on_dbl_grid(density_work, &                 ! output
            denskern, overlap, ngwfs_on_grid, ngwf_basis)       ! input

       ! ndmh: copy or upscale to 'fine' grid as required
       call density_workspace_to_output(density_work,density_fine, &
            pub_dbl_grid,grid)

       deallocate(density_work,stat=ierr)
       call utils_dealloc_check('density_on_grid','density_work',ierr)

    end if

  end subroutine density_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_workspace_to_output(density_work,density_out, &
       work_grid,out_grid)

    !==========================================================================!
    ! This subroutine transfers a density from the workspace grid to a version !
    ! on the output grid. On the workspace grid all the slabs of the local     !
    ! node's comms group are duplicated on each node of the group. The grids   !
    ! may be of different scales, in which case the result is interpolated     !
    ! from the workspace up to the output grid.                                !
    !==========================================================================!
    !  Arguments:                                                              !
    !    density_work (in)  : Workspace array (sized for one group's slabs)    !
    !    density_out  (out) : Output array (sized for one node's slabs)        !
    !    work_grid    (in)  : GRID_INFO defining workspace grid                !
    !    out_grid     (in)  : GRID_INFO defining output grid                   !
    !==========================================================================!
    ! Written by Nicholas Hine on 16/03/2011.                                  !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root, comms_reduce, pub_comms_group_size, &
         pub_group_comm, pub_first_node_in_group, pub_my_node_id
    use constants, only: DP
    use fourier, only: fourier_interpolate_cell
    use simulation_cell, only: pub_cell

    ! Arguments
    type(GRID_INFO), intent(in) :: work_grid
    type(GRID_INFO), intent(in) :: out_grid
    real(kind=DP), intent(inout) :: density_work(work_grid%ld1,work_grid%ld2, &
         work_grid%max_group_slabs12,pub_cell%num_spins)
    real(kind=DP), intent(out) :: density_out(out_grid%ld1,out_grid%ld2, &
         out_grid%max_slabs12,pub_cell%num_spins)

    ! Local Variables
    integer :: is
    integer :: node
    integer :: i3_start, i3_finish
    logical :: work_is_out

    work_is_out = .false.
    if ((work_grid%n1==out_grid%n1).and.(work_grid%n2==out_grid%n2).and. &
         (work_grid%n3==out_grid%n3)) work_is_out = .true.

    ! ndmh: transfer data from density_work to density_out
    do is=1,pub_cell%num_spins
       if (pub_comms_group_size>1) then
          ! ndmh: sum density data over all nodes in this group
          do node=0,pub_comms_group_size-1
             i3_start = work_grid%first_slab12(node + &
                  pub_first_node_in_group) - &
                  work_grid%first_slab12(pub_first_node_in_group) + 1
             i3_finish = work_grid%last_slab12(node + &
                  pub_first_node_in_group) - &
                  work_grid%first_slab12(pub_first_node_in_group) + 1
             ! ndmh: if workspace grid is output grid, sum result
             ! ndmh: directly in the density_out array, otherwise
             ! ndmh: sum it in-place and transfer
             if (work_is_out) then
                call comms_reduce('SUM',density_out(:,:, &
                     1:(i3_finish-i3_start+1),is), &
                     comm=pub_group_comm,root=node, &
                     d_array_src=density_work(:,:,i3_start:i3_finish,is))
             else
                call comms_reduce('SUM',density_work(:,:, &
                     i3_start:i3_finish,is),comm=pub_group_comm,root=node)
             end if
          end do
       else ! ndmh: no summation required - copy straight to fine grid
          if (work_is_out) then
             i3_start = work_grid%my_first_slab12_in_group
             i3_finish = work_grid%my_last_slab12_in_group
             density_out(:,:,1:work_grid%num_my_slabs12,is) = &
                  density_work(:,:,i3_start:i3_finish,is)
          end if
       end if
       ! ndmh: now interpolate from workspace grid to output grid if required
       if (.not.work_is_out) then
          i3_start = work_grid%my_first_slab12_in_group
          i3_finish = i3_start + work_grid%max_slabs12 - 1
          call fourier_interpolate_cell(density_work(:,:,i3_start:i3_finish,is), &
               density_out(:,:,:,is),work_grid,out_grid,apply_nyquist=.true.)
       end if
       ! ndmh: zero any padding in density_out
       density_out(:,:,(out_grid%num_my_slabs12+1):out_grid%max_slabs12,is) &
            = 0.0_DP

    end do

  end subroutine density_workspace_to_output


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_on_coarse_grid(density_std, &                     ! output
       denskern, overlap, ngwfs_on_grid, ngwf_basis)                   ! input

    !==========================================================================!
    ! This subroutine calculates the charge density on the coarse grid in the  !
    ! simulation cell.                                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! density  (output) : The total charge density of the system on the        !
    ! coarse grid of the simulation cell                                       !
    ! denskern      (input) : Density kernel sparse matrix                     !
    ! overlap       (input) : Overlap matrix                                   !
    ! ngwfs_on_grid (input) : The NGWFs on the real-space ppd representation   !
    !  that belong to this node                                                !
    ! ngwf_basis            : The function basis for the NGWFs                 !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/02/2011.                                  !
    !==========================================================================!

    use basis, only: basis_location_func_wrt_cell, basis_copy_function_to_box
    use cell_grid, only: cell_grid_deposit_box, pub_std_grid
    use comms, only: comms_barrier, pub_on_root, pub_my_node_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS, function_basis_sum_ppd_funcs
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(out) :: density_std(pub_std_grid%ld1, &
         pub_std_grid%ld2, pub_std_grid%max_group_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    integer :: ierr
    integer :: is
    integer :: local_col
    integer :: col_cell_start1,col_cell_start2,col_cell_start3
    integer :: box_n1,box_n2,box_n3
    real(kind=DP), allocatable :: row_sum_on_grid(:,:)
    real(kind=DP), allocatable :: density_box(:,:,:)
    real(kind=DP), allocatable :: density_buffer(:,:,:)
    logical :: i_have_box

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering density_on_coarse_grid'
#endif

    call timer_clock('density_on_coarse_grid', 1)

    box_n1 = pub_maxtight_pts1
    box_n2 = pub_maxtight_pts2
    box_n3 = pub_maxtight_pts3

    ! Allocate temporary storage
    allocate(row_sum_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('density_on_coarse_grid','row_sum_on_grid',ierr)
    allocate(density_box(box_n1,box_n2,box_n3),stat=ierr)
    call utils_alloc_check('density_on_coarse_grid','density_box',ierr)
    allocate(density_buffer(box_n1,box_n2,pub_std_grid%max_group_slabs12), &
         stat=ierr)
    call utils_alloc_check('density_on_coarse_grid','density_buffer',ierr)

    row_sum_on_grid = 0.0_DP
    density_box = 0.0_DP
    density_std = 0.0_DP

    ! Calculate \sum_b K^ab \phi_b(r)
    call function_basis_sum_ppd_funcs(row_sum_on_grid,ngwf_basis,denskern, &
         1,pub_cell%num_spins,overlap,ngwfs_on_grid,ngwf_basis)

    do is=1,pub_cell%num_spins

       ! Multiply by \phi_a(r)
       row_sum_on_grid(:,is) = row_sum_on_grid(:,is)*ngwfs_on_grid(:)

       ! Loop over col functions up to max on any node
       do local_col=1,ngwf_basis%max_on_node

          ! If there are any cols left on this node...
          if (local_col <= ngwf_basis%num_on_node(pub_my_node_id)) then

             ! Copy density for this column to box
             call basis_copy_function_to_box(density_box,box_n1,box_n2,box_n3, &
                  1,1,1,ngwf_basis%tight_boxes(local_col), &
                  row_sum_on_grid(:,is),ngwf_basis%spheres(local_col))

             ! Find position of col function wrt simulation cell:
             call basis_location_func_wrt_cell(col_cell_start1, &
                  col_cell_start2,col_cell_start3, &
                  ngwf_basis%tight_boxes(local_col))

             i_have_box = .true.
          else
             col_cell_start1 = -1234
             col_cell_start2 = -1234
             col_cell_start3 = -1234
             i_have_box = .false.
          end if

          ! Deposit box of density to whole-cell grid
          call cell_grid_deposit_box(density_std(:,:,:,is),&
               density_box, density_buffer, pub_std_grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               col_cell_start1, col_cell_start2, col_cell_start3, &
               i_have_box, .true.)

       end do

    end do

    ! Deallocate workspace
    deallocate(density_buffer,stat=ierr)
    call utils_dealloc_check('density_on_coarse_grid','density_buffer',ierr)
    deallocate(density_box,stat=ierr)
    call utils_dealloc_check('density_on_coarse_grid','density_box',ierr)
    deallocate(row_sum_on_grid,stat=ierr)
    call utils_dealloc_check('density_on_coarse_grid','row_sum_on_grid',ierr)

    ! Re-sync nodes
    call comms_barrier

    call timer_clock('density_on_coarse_grid', 2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving density_on_coarse_grid'
#endif

  end subroutine density_on_coarse_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_on_dbl_grid(density_dbl, &                        ! output
       denskern, overlap, ngwfs_on_grid, ngwf_basis)                   ! input

    !==========================================================================!
    ! This subroutine calculates the charge density on the double grid in the  !
    ! simulation cell.                                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! density_dbl  (output) : The total charge density of the system on the    !
    ! double grid of the simulation cell                                       !
    ! denskern      (input) : Density kernel sparse matrix                     !
    ! overlap       (input) : Overlap matrix                                   !
    ! ngwfs_on_grid (input) : The NGWFs on the real-space ppd representation   !
    !  that belong to this node                                                !
    ! ngwf_basis            : The function basis for the NGWFs                 !
    !--------------------------------------------------------------------------!
    ! Key internal variables:                                                  !
    !   batch_size: This is set equal to density_batch_size which comes from   !
    !     the rundat module. The value of batch size determines how large is   !
    !     the batch of accumulated fftboxes. Increasing this number decreases  !
    !     the communication per processor but increases the allocated memory   !
    !     per processor.                                                       !
    !--------------------------------------------------------------------------!
    ! Comment on parallel efficiency:                                          !
    ! This code for computing the charge density in parallel has good parallel !
    ! scalability because it adheres to two important conditions:              !
    ! 1) Interpolation, multiplication and deposition to the fine-grid         !
    !    charge density in the simulation cell happens only once per           !
    !    column of the density kernel.                                         !
    ! 2) Only NGWFs in ppd representation are communicated between the         !
    !    processors in one-to-one non-blocking fashion. This is crucial        !
    !    as communication of NGWFs in fftboxes rather than ppd representation  !
    !    is too costly and destroys parallel scalability.                      !
    ! The above goals are achieved by calculating the charge density in two    !
    ! stages: First density_batch_row_sums (involves one-to-one communication) !
    ! is called to accumulate sums of rows for a batch of columns of the       !
    ! density kernel. Then the subroutine density_batch_interp_deposit         !
    ! (involves no communication) is called to interpolate the accumulated     !
    ! sums and their columns and deposit them in the fine grid of the          !
    ! simulation cell to build the charge density.                             !
    !--------------------------------------------------------------------------!
    ! Originaly written by Chris-Kriton Skylaris in January 2001,              !
    ! to use a "pair-box".                                                     !
    ! Improved by Oswaldo Dieguez in 2001 so that it "cshifts" the "pair-box"  !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box" and    !
    ! interpolates only once per row of the density kernel.                    !
    ! Rewritten and parallelised by Chris-Kriton Skylaris in November 2003.    !
    ! Improved parallel version written by Chris-Kriton Skylaris on 15/1/2004. !
    ! Modified by Peter Haynes, July 2006 to use parallel SPAM 2.              !
    ! Spin polarised by Peter Haynes, July 2006                                !
    ! Modified in July 2009 by Nicholas Hine to use function_basis type and    !
    ! function sum routines.                                                   !
    ! Created in current form by Nicholas Hine on 28/02/2011 - Code all comes  !
    ! from previous version of density_on_fine_grid.                           !
    !==========================================================================!

    use cell_grid, only: pub_dbl_grid
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, stdout
    use fourier, only: fourier_interpolate_cell
    use function_basis, only: FUNC_BASIS, function_basis_sum_fftbox_batch
    use rundat, only: density_batch_size
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins,1)
    type(SPAM3), intent(in) :: overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(out) :: density_dbl(pub_dbl_grid%ld1, &
         pub_dbl_grid%ld2, pub_dbl_grid%max_group_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    integer :: batch_size
    integer :: max_current_size
    integer :: batch_count
    integer :: n_batches
    integer :: local_start, local_end
    integer :: ierr    ! pdh: error flag
    integer :: idx_len
    integer, allocatable :: overlap_idx(:)            ! Index for overlap matrix
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: row_fftbox_sum_batch
    real(kind=DP), allocatable, dimension(:) :: row_on_grid_buffer
    real(kind=DP), allocatable, dimension(:,:,:) :: fftbox_buffer
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: row_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:) :: buffer_dbl

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering density_on_dbl_grid'
#endif

    call timer_clock('density_on_dbl_grid',1)

    ! Obtain index for overlap matrix
    idx_len = sparse_index_length(overlap)
    allocate(overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','overlap_idx',ierr)
    call sparse_generate_index(overlap_idx,overlap)

    ! Set batch size according to input file
    batch_size = density_batch_size

    ! Allocate workspace
    allocate(row_fftbox_sum_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3, &
         pub_cell%num_spins,batch_size),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    allocate(row_on_grid_buffer(ngwf_basis%func_on_grid_buffer_size),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','row_on_grid_buffer',ierr)
    allocate(fftbox_buffer(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','fftbox_buffer',ierr)
    allocate(row_fftbox_dbl(pub_fftbox%total_ld1_dbl, &
         pub_fftbox%total_ld2_dbl,pub_fftbox%total_pt3_dbl, &
         pub_cell%num_spins,2),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','row_fftbox_dbl',ierr)
    allocate(buffer_dbl(pub_fftbox%total_ld1_dbl, pub_fftbox%total_ld2_dbl, &
         pub_dbl_grid%max_group_slabs12, pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('density_on_dbl_grid','buffer_dbl',ierr)

    ! cks: zero before accumulation
    density_dbl = 0.0_DP

    ! pdh: zero workspace
    row_fftbox_dbl = 0.0_DP

    ! cks: number of row-steps per row-block
    n_batches = ngwf_basis%max_on_node / batch_size
    if (mod(ngwf_basis%max_on_node, batch_size) > 0) n_batches = n_batches + 1

    ! cks: loop over batches of NGWFs belonging to pub_my_node_id
    local_start = 1
    do batch_count=1,n_batches

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ', &
         batch_count,' of ',n_batches,' in density_on_fine_grid'
#endif

       ! cks: limits of my current batch in pub_my_node_id-local NGWF counting
       ! cks: scheme
       local_end = min(local_start+batch_size-1,ngwf_basis%node_num)

       ! cks: maximum size of current batch over all nodes
       max_current_size = local_end - local_start + 1
       call comms_reduce('MAX', max_current_size)

       ! cks: zero before accumulation
       row_fftbox_sum_batch = 0.0_DP

       call timer_clock('density_batch_row_sums',1)

       ! ndmh: density-kernel block contributions to fftboxes of current batch
       call function_basis_sum_fftbox_batch(row_fftbox_sum_batch, &
            ngwfs_on_grid, ngwf_basis, ngwf_basis, batch_size, local_start, &
            local_end, overlap_idx, idx_len, denskern, 1, 1,&
            row_on_grid_buffer, 1.0_DP)

       call timer_clock('density_batch_row_sums',2)
       call timer_clock('density_batch_interp_deposit',1)

       ! cks: interpolation of fftboxes of current batch, multiplication
       ! cks: and accumulation to charge-density 12-slabs of each node
       call density_batch_interp_deposit(density_dbl, &            ! inout
            row_fftbox_sum_batch, ngwfs_on_grid, ngwf_basis, &     ! input
            local_start, local_end,batch_size, max_current_size, & ! input
            fftbox_buffer, row_fftbox_dbl, buffer_dbl)             ! work

       call timer_clock('density_batch_interp_deposit',2)

       ! Update first NGWF of next batch
       local_start = local_start + batch_size

    end do

    ! Deallocate workspace
    deallocate(buffer_dbl,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','buffer_dbl',ierr)
    deallocate(row_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','row_fftbox_dbl',ierr)
    deallocate(fftbox_buffer,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','fftbox_buffer',ierr)
    deallocate(row_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','row_on_grid_buffer',ierr)
    deallocate(row_fftbox_sum_batch,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','row_fftbox_sum_batch',ierr)
    deallocate(overlap_idx,stat=ierr)
    call utils_dealloc_check('density_on_dbl_grid','overlap_idx',ierr)

    call timer_clock('density_on_dbl_grid',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving density_on_dbl_grid'
#endif

  end subroutine density_on_dbl_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_batch_interp_deposit(density_dbl, &          ! input-output
       row_fftbox_sum_batch, ngwfs_on_grid, ngwf_basis, &         ! input
       local_start, local_end, batch_size, max_current_size, &    ! input
       col_fftbox, row_fftbox_dbl, buffer_dbl) ! workspace

    !==========================================================================!
    ! This subroutines interpolates each of the sums-in-fftboxes of the        !
    ! batch and its corresponding col NGWF. Then it deposits them in the       !
    ! correct position in the fine grid charge density in the simulation cell. !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !  density_dbl (input/output): Charge density in double-grid real space    !
    !    representation.                                                       !
    !  row_fftbox_sum_batch (input): Batch of fftboxes. Each fftbox            !
    !   contains the quantity phi_b*K^{ba} for the function phi_a.             !
    !  ngwfs_on_grid (input): Current NGWFs for pub_my_node_id in ppd          !
    !    representation.                                                       !
    !  local_start (input): Number of the first NGWF of the current batch      !
    !    in the counting scheme of all NGWFs of pub_my_node_id.                !
    !  local_end (input): Number of the last NGWF of the current batch         !
    !    in the counting scheme of all NGWFs of pub_my_node_id.                !
    !  batch_size (input): Number of fftboxes (phi_a NGWFs) in each batch      !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 15/1/2004 for the ONETEP code.       !
    ! Modified by Chris-Kriton Skylaris on 15/6/2004 so that it works          !
    ! with the data-parallel charge density.                                   !
    ! Modified by Nicholas Hine on 21/07/2009 to use function basis type.      !
    ! Modified by Nicholas Hine on 04/09/2009 to combine fourier_interpolate   !
    ! and multiplication - removes need for col_fftbox_dbl.                    !
    ! Modified by Nicholas Hine in November 2009 to accumulate the density for !
    ! each atom before depositing it to cell, for a parallel-efficiency gain.  !
    !--------------------------------------------------------------------------!

    use basis, only: basis_copy_function_to_box, &
         basis_ket_start_wrt_fftbox, basis_location_func_wrt_cell
    use cell_grid, only: cell_grid_deposit_box, pub_dbl_grid
    use comms, only: comms_abort, comms_barrier, pub_my_node_id, &
         pub_on_root, pub_rank_comm
    use constants, only: DP, stdout
    use fourier, only: fourier_interpolate_product
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell, pub_fftbox
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer, intent(in) :: local_start
    integer, intent(in) :: local_end
    integer, intent(in) :: batch_size
    integer, intent(in) :: max_current_size
    real(kind=DP), intent(inout) :: density_dbl(pub_dbl_grid%ld1, &
         pub_dbl_grid%ld2, pub_dbl_grid%max_group_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(in) :: row_fftbox_sum_batch(&
         pub_fftbox%total_ld1, pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
         pub_cell%num_spins, batch_size)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(out) :: col_fftbox(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)
    real(kind=DP), intent(out) :: row_fftbox_dbl(pub_fftbox%total_ld1_dbl, &
         pub_fftbox%total_ld2_dbl, pub_fftbox%total_pt3_dbl, &
         pub_cell%num_spins, 2)
    real(kind=DP), intent(out) :: buffer_dbl(pub_fftbox%total_ld1_dbl, &
         pub_fftbox%total_ld2_dbl, pub_dbl_grid%max_group_slabs12)

    ! cks: <<local variables>>
    integer :: is
    integer :: col_start1, col_start2, col_start3
    integer :: col_cell_start1, col_cell_start2, col_cell_start3
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: col, local_col
    integer :: atom_of_col
    integer :: first_on_col_atom, last_on_col_atom
    integer :: batch_count
    logical :: i_have_box

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering density_batch_interp_deposit'
#endif

    ! cks: the col function will always be at the centre of the fftbox
    ! cks: or left where it is in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(col_start1,col_start2,col_start3, &
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    ! cks: loop over the members of the largest current batch (from all nodes)
    ! cks: and interpolate each pair of fftboxes of each node, multiply them
    ! cks: together and deposit them accordingly to the 12-slabs of all nodes
    batch_count =0
    do local_col=local_start,local_start+max_current_size-1
       batch_count = batch_count + 1

       call timer_clock('density_fftbox_interpolate_multiply', 1)

       ! cks: interpolate and multiply my col NGWF if pub_my_node_id
       ! cks: posseses an fftbox for this local_col value
       if ( local_col .le. local_end) then

          ! ndmh: find information about this column NGWF
          col = local_col + ngwf_basis%first_on_node(pub_my_node_id) - 1
          atom_of_col = ngwf_basis%atom_of_func(col)
          first_on_col_atom = ngwf_basis%first_on_atom(atom_of_col)
          last_on_col_atom = first_on_col_atom + &
               ngwf_basis%num_on_atom(atom_of_col) - 1

          ! ndmh: copy column NGWF to column fftbox
          call basis_copy_function_to_box(col_fftbox, &
               pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3, &
               col_start1, col_start2, col_start3, &
               ngwf_basis%tight_boxes(local_col), ngwfs_on_grid, &
               ngwf_basis%spheres(local_col))

          ! ndmh: interpolate col and sum of rows together with c2c ffts
          ! ndmh: then multiply them together to find contribution to the
          ! ndmh: charge density
          do is=1,pub_cell%num_spins
             call fourier_interpolate_product( &
                  row_fftbox_sum_batch(:,:,:,is,batch_count), col_fftbox, &
                  row_fftbox_dbl(:,:,:,is,1))
          end do

          ! cks: find where tightbox of local_col is located in simulation cell
          call basis_location_func_wrt_cell( &
               col_cell_start1, col_cell_start2, col_cell_start3, &
               ngwf_basis%tight_boxes(local_col))

          ! ndmh: if this is a new atom or the start of a new batch, reset the
          ! ndmh: accumulated density for this atom to the density for this NGWF
          if ((col == first_on_col_atom).or.(local_col==local_start)) then
             row_fftbox_dbl(:,:,:,:,2) = row_fftbox_dbl(:,:,:,:,1)
          else ! ndmh: add the density for this NGWF to the sum for this atom
             row_fftbox_dbl(:,:,:,:,2) = row_fftbox_dbl(:,:,:,:,2) + &
                  row_fftbox_dbl(:,:,:,:,1)
          end if

          ! ndmh: if this is the last NGWF on this atom or the end of the batch,
          ! ndmh: we will need to deposit the density to the whole-cell array
          if ((col == last_on_col_atom).or.(local_col==local_end)) then
             i_have_box = .true.
          else
             i_have_box = .false.
          end if

       else

          ! ndmh: prevent passing of unassigned variables in debug mode
          col_cell_start1 = -1234
          col_cell_start2 = -1234
          col_cell_start3 = -1234
          i_have_box =.false.

       end if

       ! ndmh: synchronise nodes in this rank so that load-balancing time is
       ! ndmh: reported in density_fftbox_interpolate_multiply rather than as
       ! ndmh: the wait for the alltoall in cell_grid_deposit_box.
       call comms_barrier(pub_rank_comm)
       call timer_clock('density_fftbox_interpolate_multiply', 2)
       call timer_clock('density_fftbox_deposit_to_cell', 1)

       ! cks: get into the depositing fftboxes to 12-slabs of all nodes
       ! cks: regardless of if pub_my_node_id has an fftbox to deposit
       ! ndmh: box_to_cell routine moved to basis_mod
       fftbox_start1_dbl = 2*(col_cell_start1 - col_start1) + 1
       fftbox_start2_dbl = 2*(col_cell_start2 - col_start2) + 1
       fftbox_start3_dbl = 2*(col_cell_start3 - col_start3) + 1
       do is=1,pub_cell%num_spins
          call cell_grid_deposit_box(density_dbl(:,:,:,is), &
               row_fftbox_dbl(:,:,:,is,2), buffer_dbl, pub_dbl_grid, &
               pub_fftbox%total_pt1_dbl, pub_fftbox%total_pt2_dbl, &
               pub_fftbox%total_pt3_dbl, pub_fftbox%total_ld1_dbl, &
               pub_fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_have_box, .true.)
       end do

       call timer_clock('density_fftbox_deposit_to_cell', 2)

    end do

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving density_batch_interp_deposit'
#endif

  end subroutine density_batch_interp_deposit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_initial_guess_recip(density, & ! output
       elements, struct_fac, grid)

    !===============================================================!
    ! This subroutine initialises an approximate (guess) charge     !
    ! density as a superposition of Gaussian functions on atoms.    !
    !---------------------------------------------------------------!
    ! Written by Arash A. Mostofi on 19/2/2004.                     !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 to work        !
    ! with data-parallel structure factor.                          !
    ! Rewritten by Peter D. Haynes on 30/6/2004 to use new Fourier  !
    ! parallelisation.                                              !
    ! Modified by Nicholas Hine on 19/09/2008 for rearranged        !
    ! structure factor array to improve cache performance           !
    ! Moved to density_mod, and dependence on p_species array       !
    ! removed so that it only uses elements array, by Nicholas Hine !
    ! on 04/11/09.                                                  !
    ! Adjusted to allow execution where number of 23-slabs is less  !
    ! than number of nodes by Nicholas Hine, December 2009.         !
    !===============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_abort, pub_my_node_id, pub_total_num_nodes
    use constants, only: UP, DN, DP, stdout
    use fourier, only: fourier_apply_cell_backward
    use integrals, only: integrals_trace_on_grid
    use ion, only: ELEMENT
    use rundat, only: pub_spin
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)   :: grid
    type(ELEMENT), intent(in)     :: elements(pub_cell%nat)
    complex(kind=DP), intent(in)  :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    real(kind=DP), intent(out) :: density(grid%ld1,grid%ld2, &
         grid%max_slabs12, pub_cell%num_spins)

    ! Local variables
    integer :: ierr                            ! Error flag
    real(kind=DP), parameter :: rad_default = 1.8_DP ! Default core radius
    integer :: i3,i2,i1,islab23,islab12
    integer :: atom,species
    real(kind=DP) :: ion_charge
    real(kind=DP) :: g(3),gsq,dens_value,maxrad,intdens
    real(kind=DP) :: up_fac,dn_fac
    real(kind=DP),allocatable :: expon(:),pref(:)
    complex(kind=DP), allocatable :: dens_recip(:,:,:)

    ! Start timer
    call timer_clock('density_initial_guess_recip',1)

    ! jd: Takes care of padding between n1 and ld1 etc.
    density = 0.0_DP

    ! Allocate complex workspace for reciprocal space
    allocate(dens_recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('density_initial_guess_recip','dens_recip',ierr)
    allocate(expon(pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('density_initial_guess_recip','expon',ierr)
    allocate(pref(pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('density_initial_guess_recip','pref',ierr)

    ! ndmh: Set up gaussian exponents and prefactors
    ! ndmh: Loop over atomic species
    do species=1,pub_cell%num_pspecies

       ! Find first example of this atom in elements array
       do atom=1,pub_cell%nat+1
          if (atom>pub_cell%nat) then
             write(stdout,'(a)') 'Error in density_initial_guess_recip: no &
                  &atom of species ',species,' found in elements array'
             call comms_abort
          end if
          if (elements(atom)%pspecies_number==species) exit
       end do

       ! Find ion charge from elements array
       ion_charge = elements(atom)%ion_charge

       ! If projector radii are defined for this species, find maximum radius
       ! amongst all shells, else use default
       maxrad = elements(atom)%max_core_radius
       if (maxrad < 0.01_DP) maxrad = rad_default

       ! If radius is zero, set to default value
       if (maxrad == 0.0_DP) then
          maxrad = rad_default
       end if

       ! Exponent of gaussian
       expon(species) = maxrad * maxrad * 0.25_DP
       expon(species) = expon(species) * 1.75_DP ! fudge factor

       ! Prefactor, normalised to ion charge of species
       pref(species) = ion_charge / grid%weight

    end do ! loop over species

    ! Loop over reciprocal space grid on this node
    do islab23=1,grid%num_slabs23          ! along b1

       dens_recip(:,:,islab23) = (0.0_DP,0.0_DP)

       do i2=1,grid%n2                     ! along b2
          do i3=1,grid%n3                  ! along b3

             ! ndmh: loop over species
             do species=1,pub_cell%num_pspecies

                call cell_grid_recip_pt(g,islab23 + &
                     grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)

                gsq = sum(g(:)**2)

                dens_value = pref(species) * exp(-expon(species)*gsq)

                dens_recip(i3,i2,islab23) = dens_recip(i3,i2,islab23) + &
                     struct_fac(species,i3,i2,islab23) * dens_value

             end do ! loop over species

          end do  ! loop along b3
       end do     ! loop along b2
    end do        ! loop along b1


    ! G=0 element must be real (find first 23-slab)
    if (pub_my_node_id==grid%node_slab23(1)) then
       if (aimag(dens_recip(1,1,1)) /= 0.0_DP) then
          write(stdout,'(a/a)') 'WARNING in density_initial_guess_recip:', &
               '  density not real - setting imaginary part of G=0 term to zero'
          dens_recip(1,1,1) = cmplx(real(dens_recip(1,1,1),kind=DP),0.0_DP, &
               kind=DP)
       end if
    end if

    ! Nyquist filter (fine grid is always going to be even)
    if (grid%num_slabs23>0) then
       dens_recip(grid%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       dens_recip(:,grid%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if

    ! Find last 23-slab and apply Nyquist filter
    if (pub_my_node_id==grid%node_slab23(grid%n1/2+1)) &
         dens_recip(:,:,grid%num_slabs23) = (0.0_DP,0.0_DP)

    ! FFT the local ionic potential from reciprocal to real space
    call fourier_apply_cell_backward(density(:,:,:,1),dens_recip,grid)

    ! Deallocate workspace
    deallocate(pref,stat=ierr)
    call utils_dealloc_check('density_initial_guess_recip','pref',ierr)
    deallocate(expon,stat=ierr)
    call utils_dealloc_check('density_initial_guess_recip','expon',ierr)
    deallocate(dens_recip,stat=ierr)
    call utils_dealloc_check('density_initial_guess_recip','dens_recip',ierr)

    ! Density may have negative values
    call density_render_possemidefinite(density(:,:,:,1),grid)

    intdens = integrals_trace_on_grid(density(:,:,:,1),grid)

    ! Renormalise density
    call density_on_grid_renorm(density(:,:,:,1), &
         grid,intdens,elements)

    ! If spin-polarised, separate into up and down densities
    if (pub_cell%num_spins == 2) then
       intdens = integrals_trace_on_grid(density(:,:,:,1),grid)
       up_fac = 0.5_DP * (intdens + pub_spin) / intdens
       dn_fac = 0.5_DP * (intdens - pub_spin) / intdens
       do islab12=1,grid%num_my_slabs12   ! along a3
          do i2=1,grid%n2                 ! along a2
             do i1=1,grid%n1              ! along a1
                dens_value = density(i1,i2,islab12,1)
                density(i1,i2,islab12,UP) = up_fac * dens_value
                density(i1,i2,islab12,DN) = dn_fac * dens_value
             end do
             density(grid%n1+1:,i2,islab12,DN) = 0.0_DP
          end do
          density(:,grid%n2+1:,islab12,DN) = 0.0_DP
       end do
       density(:,:,grid%num_my_slabs12+1:,DN) = 0.0_DP
    end if

    ! Stop timer
    call timer_clock('density_initial_guess_recip',2)

  end subroutine density_initial_guess_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_initial_guess_real(density_on_grid, & ! output
       elements, grid, add_aug_den)

    !===============================================================!
    ! This subroutine initialises an approximate (guess) charge     !
    ! density as a superposition of Gaussian functions on atoms.    !
    !---------------------------------------------------------------!
    ! Written by Nicholas Hine in May 2011.                         !
    !===============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: comms_abort, pub_my_node_id, pub_total_num_nodes, &
         pub_on_root
    use constants, only: UP, DN, DP, PI, stdout, ANGSTROM
    use integrals, only: integrals_trace_on_grid
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_elements_on_node, pub_max_atoms_on_node, &
         pub_first_atom_on_node, pub_num_atoms_on_node
    use rundat, only: pub_spin, pub_aug, pub_devel_code
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)   :: grid
    type(ELEMENT), intent(in)     :: elements(pub_cell%nat)
    real(kind=DP), intent(out) :: density_on_grid(grid%ld1,grid%ld2, &
         grid%max_slabs12, pub_cell%num_spins)
    logical, intent(in) :: add_aug_den

    ! Local variables
    integer :: ierr
    integer :: iat, loc_iat
    integer :: isp
    integer :: is
    integer :: box_n1, box_n2, box_n3
    integer :: box_start1, box_start2, box_start3
    logical :: i_have_box
    real(kind=DP) :: spin_fac
    real(kind=DP) :: intdens,scale
    real(kind=DP),allocatable :: density_box(:,:,:,:), buffer(:,:,:)

    ! Start timer
    call timer_clock('density_initial_guess_real',1)

    ! jd: Takes care of padding between n1 and ld1 etc.
    density_on_grid = 0.0_DP

    if (pub_cell%num_spins==1) spin_fac = 1.0_DP
    if (pub_cell%num_spins==2) spin_fac = 0.5_DP

    ! ndmh: Loop over atomic species to check densities have been initialised
    do isp=1,pub_cell%num_species

       if (pub_on_root) then

          ! ndmh: check if we already have a density for this element
          if (.not.radial_densities(isp)%present) then

             write(stdout,'(a,i3)') 'Error in density_initial_guess_real: &
                  &atom density for species',isp
             write(stdout,'(a)') 'has not been initialised. All species must &
                  &have atom density calculated to'
             write(stdout,'(a)') 'use real-space density calculation.'
             call comms_abort

          end if

       end if

       ! ndmh: share this radial density with other nodes
       call density_radial_bcast(radial_densities(isp))

    end do ! loop over species

    ! Pick a box size: use tightbox size scaled to this grid
    scale = real(grid%n1,kind=DP)/real(pub_cell%total_pt1,kind=DP)
    box_n1 = min(int(pub_maxtight_pts1*scale + 1),grid%n1)
    scale = real(grid%n2,kind=DP)/real(pub_cell%total_pt2,kind=DP)
    box_n2 = min(int(pub_maxtight_pts2*scale + 1),grid%n2)
    scale = real(grid%n3,kind=DP)/real(pub_cell%total_pt3,kind=DP)
    box_n3 = min(int(pub_maxtight_pts3*scale + 1),grid%n3)

    ! Allocate workspace
    allocate(density_box(box_n1,box_n2,box_n3,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('density_initial_guess_real','density_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('density_initial_guess_real','buffer',ierr)

    ! Loop over atoms on this node
    do loc_iat=1,pub_max_atoms_on_node

       if (loc_iat<=pub_num_atoms_on_node(pub_my_node_id)) then
          iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
          isp = pub_elements_on_node(loc_iat)%species_number

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               pub_elements_on_node(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid)

          ! Reset box to zero
          density_box(:,:,:,:) = 0.0_DP
          do is=1,pub_cell%num_spins

             ! Add in augmentation density if required
             if (pub_aug.and.add_aug_den) radial_densities(isp)%den(:,is) = &
                  radial_densities(isp)%den(:,is) &
                  + radial_densities(isp)%aug_den(:,is)

             ! Transfer the radial density to the 3D box
             call internal_density_to_box(density_box(:,:,:,is), &
                  radial_densities(isp)%den(:,is),radial_densities(isp)%rad, &
                  pub_elements_on_node(loc_iat)%radius,radial_densities(isp)%npts, &
                  box_n1,box_n2,box_n3,box_start1,box_start2,box_start3,grid, &
                  pub_elements_on_node(loc_iat)%centre)

             ! Subtract off augmentation density if required
             if (pub_aug.and.add_aug_den) radial_densities(isp)%den(:,is) = &
                  radial_densities(isp)%den(:,is) &
                  - radial_densities(isp)%aug_den(:,is)

             ! Scale by 1/2 if spin-polarised (atom solver density is
             ! always non-polarised)
             density_box = density_box * spin_fac

          end do  ! is

          i_have_box = .true.
       else
          ! Nothing to deposit on this node
          i_have_box = .false.
       end if

       ! Deposit this box to the simulation cell if present, or just wait for
       ! data from other nodes if no box
       do is=1,pub_cell%num_spins
          call cell_grid_deposit_box(density_on_grid(:,:,:,is), &
               density_box(:,:,:,is), buffer, grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_have_box, .false.)
       end do
    end do

    ! Deallocate
    deallocate(density_box,stat=ierr)
    call utils_dealloc_check('density_initial_guess_real','density_box',ierr)
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('density_initial_guess_real','buffer',ierr)

    ! Density may have negative values
    call density_render_possemidefinite(density_on_grid(:,:,:,1),grid)

    intdens = integrals_trace_on_grid(density_on_grid(:,:,:,1),grid)

    ! Renormalise density
    if (add_aug_den) call density_on_grid_renorm(density_on_grid(:,:,:,1), &
         grid,intdens,elements)

    ! ndmh: show density if required
    if (index(pub_devel_code,'WRITE_DENS_GUESS')>0) then
       if (add_aug_den) then
          call visual_scalarfield(density_on_grid,grid, &
               'Guess density (in e/ang^3) for:', '_guessdensity_aug', &
               elements, ANGSTROM**3)
       else
          call visual_scalarfield(density_on_grid,grid, &
               'Guess density (in e/ang^3) for:', '_guessdensity', &
               elements, ANGSTROM**3)
       end if
    end if
    ! Stop timer
    call timer_clock('density_initial_guess_real',2)

contains

    subroutine internal_density_to_box(box,den,r,rcut,npts, &
         box_n1,box_n2,box_n3,cell_start1,cell_start2,cell_start3,grid, &
         atom_origin)

      use basis, only: basis_box_origin_to_atom
      use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
           OPERATOR(*), geometry_magnitude
      use services, only: services_locate_interp,services_linear_interpolation

      ! Arguments
      integer,intent(in) :: box_n1, box_n2, box_n3
      real(kind=DP),intent(inout) :: box(box_n1,box_n2,box_n3)
      integer,intent(in) :: npts
      real(kind=DP), intent(in) :: den(npts)
      real(kind=DP), intent(in) :: r(npts), rcut
      integer,intent(in) :: cell_start1
      integer,intent(in) :: cell_start2
      integer,intent(in) :: cell_start3
      type(GRID_INFO),intent(in) :: grid
      type(POINT),intent(in) :: atom_origin

      ! Local Variables
      integer :: i1,i2,i3
      integer :: ipt
      type(POINT) :: box_origin
      type(POINT) :: box_to_atom
      type(POINT) :: r_cell, r_sphere
      real(kind=DP) :: rmag

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

               if (rmag>rcut) then
                  box(i1,i2,i3) = 0.0_DP
               else
                  ipt = services_locate_interp(rmag,r,npts)
                  box(i1,i2,i3) = services_linear_interpolation(rmag,den(ipt), &
                       den(ipt+1),r(ipt),r(ipt+1))
               end if

            end do
         end do
      end do

    end subroutine internal_density_to_box

  end subroutine density_initial_guess_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_check_possemidefinite(density_fine,grid)

    !================================================================!
    ! This subroutine is checks if a given density on the fine grid  !
    ! has negative values and stops if does.                         !
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                      !
    ! Modified by Quintin Hill to use pub_cell on 17/10/2008.        !
    !================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_abort
    use constants, only: DP, stdout
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1,&
         grid%ld2,grid%n3)

    ! cks: internal definitions
    integer :: r1, r2, r3


    do r3=1,grid%n3
       do r2=1,grid%n2
          do r1=1,grid%n1

             if ( density_fine(r1,r2,r3).lt.(0.0_DP)) then
                write(stdout,*)'density_fine(r1,r2,r3)=',density_fine(r1,r2,r3)
                write(stdout,*)'r1,r2,r3=',r1,r2,r3
                write(stdout,*)&
                     'NEGATIVE DENSITY value: density_check_possemidefinite'
                write(stdout,*)'ONETEP execution STOPS'
                call comms_abort
             endif

          enddo
       enddo
    enddo


  end subroutine density_check_possemidefinite


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_render_possemidefinite(density,grid)

    !=======================================================================!
    ! This subroutine renders the density of the fine grid positive         !
    ! semidefinite by zeroing out all negative values                       !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in spring 2001.                      !
    ! Modified by Chris-Kriton Skylaris on 18/7/2004 so that it works with  !
    ! the parallel version of ONETEP.                                       !
    !=======================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: density(grid%ld1,grid%ld2,grid%max_slabs12)

    ! cks: internal declarations
    integer :: row1, row2, islab12

    do islab12=1,grid%num_my_slabs12
       do row2=1,grid%n2
          do row1=1,grid%n1
             density(row1,row2,islab12) = &
                  max(density(row1,row2,islab12),0.0_DP)
          end do
       end do
    end do

  end subroutine density_render_possemidefinite


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_on_grid_renorm(density,grid,integrated_density,elements)

    !======================================================================!
    ! This subroutine renormalises the density on the fine grid to the     !
    ! correct number of electrons.                                         !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001.                            !
    ! Modified on 1/9/2004 by Chris-Kriton Skylaris to add external charge !
    ! Modified on 17/10/2008 by Quintin Hill to use pub_cell.              !
    !======================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP
    use ion, only: element
    use rundat, only: pub_charge
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: density(grid%ld1,grid%ld2, &
         grid%max_slabs12)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=DP), intent(in) :: integrated_density

    ! cks: internal declarations
    integer :: row, ne
    real(kind=DP) :: ne_real
    real(kind=DP) :: factor

    ! cks: find the total number of electrons in the simulation cell
    ! ndmh: add as reals
    ne_real = 0
    do row=1,pub_cell%nat
       ne_real = ne_real + elements(row)%ion_charge
    end do

    ! cks, 1/9/2004: addition of external charge
    ! ndmh: take int part of result (for fractional charges)
    ne = int(ne_real - pub_charge)

    ! cks: renormalise density
    factor = real(ne,kind=DP) / integrated_density
    density = factor * density

  end subroutine density_on_grid_renorm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_plot_slice(coord,xvec,yvec,ovec,nx,ny,data,n1,n2,n3, &
       fname,title,xlabel,ylabel)

    !=================================================================!
    ! This routine plots a slice of a 3D array in MTV format          !
    !-----------------------------------------------------------------!
    ! Written by Peter D. Haynes in summer 2004.                      !
    !=================================================================!


    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, stdout, PI
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    character, intent(in) :: coord         ! (C)artesian or (F)ractional
    real(kind=DP), intent(in) :: xvec(:)   ! Vector giving x-axis
    real(kind=DP), intent(in) :: yvec(:)   ! Vector giving y-axis
    real(kind=DP), intent(in) :: ovec(:)   ! Vector giving origin
    integer, intent(in) :: nx              ! Number of points along x
    integer, intent(in) :: ny              ! Number of points along y
    real(kind=DP), intent(in) :: data(:,:,:)  ! The data to plot
    integer, intent(in) :: n1,n2,n3        ! Size of data
    character(len=*), intent(in) :: fname  ! Filename
    character(len=*), optional, intent(in) :: title
    character(len=*), optional, intent(in) :: xlabel
    character(len=*), optional, intent(in) :: ylabel

    ! Local variables
    integer, parameter :: iounit = 20   ! Unit for output file
    integer, parameter :: finesse = 4   ! Accuracy of spline interpolation
    integer :: ierr                     ! Error flag
    integer :: ix,iy                    ! Plot point counters
    integer :: i1,i2,i3                 ! Interpolation counters
    integer :: j1,j2,j3                 ! Interpolation counters
    integer :: i,j,k                    ! Loop counters
    integer :: n(3)                     ! Copy of n1,n2,n3
    integer :: m1,m2,m3                 ! Number of points to interpolate from
    real(kind=DP), parameter :: cutoff = 2.0_DP     ! Interpolation cutoff
    real(kind=DP) :: x,h,delta                      ! Grid spacing/position
    real(kind=DP) :: normfac                        ! PSinc normalisation
    real(kind=DP) :: xaxrel(3),yaxrel(3),orgrel(3)  ! Plot frame definition
    real(kind=DP) :: ptrel(3)                       ! Plot point
    real(kind=DP) :: psinc23,psinc3                 ! Interpolation variables
    real(kind=DP), allocatable :: psinc(:,:)        ! Value of PSinc for spline
    real(kind=DP), allocatable :: d2psinc(:,:)      ! Value of 2nd deriv
    real(kind=DP), allocatable :: psincpt(:,:)      ! Value of PSinc at point
    real(kind=DP), allocatable :: line(:)           ! Line of data for output

    ! Check arguments
    if (coord /= 'C' .and. coord /= 'c' .and. coord /= 'F' .and. &
         coord /= 'f') then
       if (pub_on_root) write(stdout,'(2a)') 'Error in density_plot_slice: &
            &invalid coordinate specifier ',coord
       call comms_abort
    end if
    if (size(xvec) < 3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in density_plot_slice: &
            &x-axis requires a 3-component vector'
       call comms_abort
    end if
    if (size(yvec) < 3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in density_plot_slice: &
            &y-axis requires a 3-component vector'
       call comms_abort
    end if
    if (size(ovec) < 3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in density_plot_slice: &
            &origin requires a 3-component vector'
       call comms_abort
    end if
    if (nx < 1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in density_plot_slice: &
            &number of x-points < 1'
       call comms_abort
    end if
    if (ny < 1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in density_plot_slice: &
            &number of y-points < 1'
       call comms_abort
    end if
    if (size(data,1) < n1) then
       if (pub_on_root) write(stdout,'(a)') 'Error in density_plot_slice: &
            &data array mismatch in index 1'
       call comms_abort
    end if
    if (size(data,2) < n2) then
       if (pub_on_root) write(stdout,'(a)') 'Error in density_plot_slice: &
            &data array mismatch in index 2'
       call comms_abort
    end if
    if (size(data,3) < n3) then
       if (pub_on_root) write(stdout,'(a)') 'Error in density_plot_slice: &
            &data array mismatch in index 3'
       call comms_abort
    end if

    if (pub_on_root) write(stdout,'(2a)') '>>>>>>>>>>>> Writing slice to file: ', &
         trim(fname)

    ! Allocate workspace
    allocate(line(nx),stat=ierr)
    call utils_alloc_check('density_plot_slice','line',ierr)
    allocate(psinc(finesse*max(n1,n2,n3),3),stat=ierr)
    call utils_alloc_check('density_plot_slice','psinc',ierr)
    allocate(d2psinc(finesse*max(n1,n2,n3),3),stat=ierr)
    call utils_alloc_check('density_plot_slice','d2psinc',ierr)
    allocate(psincpt(max(n1,n2,n3),3),stat=ierr)
    call utils_alloc_check('density_plot_slice','psincpt',ierr)

    ! Convert to fractional coordinates if necessary
    if (coord == 'C' .or. coord == 'c') then
       xaxrel(1) = (xvec(1) * pub_cell%b1%x + xvec(2) * pub_cell%b1%y + &
            xvec(3) * pub_cell%b1%z) * 0.5_DP / pi
       xaxrel(2) = (xvec(1) * pub_cell%b2%x + xvec(2) * pub_cell%b2%y + &
            xvec(3) * pub_cell%b2%z) * 0.5_DP / pi
       xaxrel(3) = (xvec(1) * pub_cell%b3%x + xvec(2) * pub_cell%b3%y + &
            xvec(3) * pub_cell%b3%z) * 0.5_DP / pi
       yaxrel(1) = (yvec(1) * pub_cell%b1%x + yvec(2) * pub_cell%b1%y + &
            yvec(3) * pub_cell%b1%z) * 0.5_DP / pi
       yaxrel(2) = (yvec(1) * pub_cell%b2%x + yvec(2) * pub_cell%b2%y + &
            yvec(3) * pub_cell%b2%z) * 0.5_DP / pi
       yaxrel(3) = (yvec(1) * pub_cell%b3%x + yvec(2) * pub_cell%b3%y + &
            yvec(3) * pub_cell%b3%z) * 0.5_DP / pi
       orgrel(1) = (ovec(1) * pub_cell%b1%x + ovec(2) * pub_cell%b1%y + &
            ovec(3) * pub_cell%b1%z) * 0.5_DP / pi
       orgrel(2) = (ovec(1) * pub_cell%b2%x + ovec(2) * pub_cell%b2%y + &
            ovec(3) * pub_cell%b2%z) * 0.5_DP / pi
       orgrel(3) = (ovec(1) * pub_cell%b3%x + ovec(2) * pub_cell%b3%y + &
            ovec(3) * pub_cell%b3%z) * 0.5_DP / pi
    else
       xaxrel = xvec(1:3)
       yaxrel = yvec(1:3)
       orgrel = ovec(1:3)
    end if

    ! Write header to file
    if (pub_on_root) then
       open(unit=iounit,file=trim(fname))
       write(iounit,'(a)') '$DATA=CONTOUR'
       write(iounit,'(a,i4)') '% XMIN=0  XMAX=1  NX=',nx
       write(iounit,'(a,i4)') '% YMIN=0  YMAX=1  NY=',ny
       if (present(title)) write(iounit,'(3a)') '% TOPLABEL="',trim(title),'"'
       write(iounit,'(a)') '% EQUALSCALE'
       if (present(xlabel)) write(iounit,'(3a)') '% XLABEL="',trim(xlabel),'"'
       if (present(ylabel)) write(iounit,'(3a)') '% YLABEL="',trim(ylabel),'"'
       write(iounit,'(a)') '% CONTSTYLE=2'
       write(iounit,'(a)') '% NSTEPS=10'
    end if

    ! Set up spline-fitted PSinc functions
    n(1) = n1 ; n(2) = n2 ; n(3) = n3
    do k=1,3
       h = 1.0_DP / real(finesse*n(k),kind=DP)
       do i=1,finesse*n(k)
          x = (i-1) * h
          psinc(i,k) = 0.5_DP
          do j=1,(n(k)-1)/2   ! omit Nyquist for even grid to keep psinc real
             psinc(i,k) = psinc(i,k) + cos(2.0_DP*pi*j*x)
          end do
       end do
       normfac = 1.0_DP / (0.5_DP + real((n(k)-1)/2,kind=DP))
       psinc(1:finesse*n(k),k) = psinc(1:finesse*n(k),k) * normfac
       call internal_spline_uni_fit(finesse*n(k),h,psinc(:,k),d2psinc(:,k))
    end do

    ! Calculate set of grid points to use in interpolation
    m1 = min(int(2.0_DP * cutoff / pub_cell%d1),(n1-1)/2)
    m2 = min(int(2.0_DP * cutoff / pub_cell%d2),(n2-1)/2)
    m3 = min(int(2.0_DP * cutoff / pub_cell%d3),(n3-1)/2)

    ! Loop over plotting points
    do iy=1,ny
       do ix=1,nx

          ! Get position of point in fractional coordinates
          ptrel = orgrel + (ix-1)*xaxrel/nx + (iy-1)*yaxrel/ny

          ! Interpolate PSinc's for this point
          h = 1.0_DP / real(finesse*n1,kind=DP)
          delta = 1.0_DP / real(n1,kind=DP)
          do j1=-m1,m1
             i1 = modulo(int(ptrel(1)*n1) + j1,n1) + 1
             x = modulo((i1-1) * delta - ptrel(1),1.0_DP)
             psincpt(i1,1) = internal_spline_uni_eval(finesse*n1,h,psinc(:,1), &
                  d2psinc(:,1),x)
          end do
          h = 1.0_DP / real(finesse*n2,kind=DP)
          delta = 1.0_DP / real(n2,kind=DP)
          do j2=-m2,m2
             i2 = modulo(int(ptrel(2)*n2) + j2,n2) + 1
             x = modulo((i2-1) * delta - ptrel(2),1.0_DP)
             psincpt(i2,2) = internal_spline_uni_eval(finesse*n2,h,psinc(:,2), &
                  d2psinc(:,2),x)
          end do
          h = 1.0_DP / real(finesse*n3,kind=DP)
          delta = 1.0_DP / real(n3,kind=DP)
          do j3=-m3,m3
             i3 = modulo(int(ptrel(3)*n3) + j3,n3) + 1
             x = modulo((i3-1) * delta - ptrel(3),1.0_DP)
             psincpt(i3,3) = internal_spline_uni_eval(finesse*n3,h,psinc(:,3), &
                  d2psinc(:,3),x)
          end do

          ! Evaluate function at this point
          line(ix) = 0.0_DP
          do j3=-m3,m3
             i3 = modulo(int(ptrel(3)*n3) + j3,n3) + 1
             psinc3 = psincpt(i3,3)
             do j2=-m2,m2
                i2 = modulo(int(ptrel(2)*n2) + j2,n2) + 1
                psinc23 = psincpt(i2,2) * psinc3
                do j1=-m1,m1
                   i1 = modulo(int(ptrel(1)*n1) + j1,n1) + 1
                   line(ix) = line(ix) + data(i1,i2,i3) * psincpt(i1,1) * &
                        psinc23
                end do
             end do
          end do

       end do

       ! Write this data to file
       if (pub_on_root) write(iounit,'(8(f9.4,1x))') line

    end do

    ! Write footer to file
    if (pub_on_root) then
       write(iounit,'(a)') '$ END'
       close(iounit)
    end if

    ! Deallocate workspace
    deallocate(psincpt,stat=ierr)
    call utils_dealloc_check('density_plot_slice','psincpt',ierr)
    deallocate(d2psinc,stat=ierr)
    call utils_dealloc_check('density_plot_slice','d2psinc',ierr)
    deallocate(psinc,stat=ierr)
    call utils_dealloc_check('density_plot_slice','psinc',ierr)
    deallocate(line,stat=ierr)
    call utils_dealloc_check('density_plot_slice','line',ierr)

  contains

    !==========================================================================
    ! Fit a natural cubic spline to data tabulated on a uniform grid
    !==========================================================================

    subroutine internal_spline_uni_fit(n,h,y,ypp)

      implicit none
      integer, intent(in) :: n
      real(kind=DP), intent(in) :: h,y(n)
      real(kind=DP), intent(out) :: ypp(n)
      real(kind=DP), allocatable :: d(:),e(:)
      integer :: i,info

      allocate(d(n),stat=info)
      call utils_alloc_check('internal_spline_uni_fit (density_mod.F90)',&
           'd',info)
      allocate(e(n-1),stat=info)
      call utils_alloc_check('internal_spline_uni_fit (density_mod.F90)',&
           'e',info)

      d(:) = 4.0_DP*h*h
      e(:) = h*h
      ypp(1) = y(2)-2.0_DP*y(1)+y(n)
      do i=2,n-1
         ypp(i) = y(i+1)-2.0_DP*y(i)+y(i-1)
      end do
      ypp(n) = y(1)-2.0_DP*y(n)+y(n-1)
      ypp(:) = 6.0_DP * ypp(:)
      call dptsv(n,1,d,e,ypp,n,info)
      if (info /= 0) then
         if (pub_on_root) write(stdout,'(a/a,i6)') &
              'Error in internal_spline_uni_fit (density_mod.F90): ', &
              'dptsv failed with code ',info
         call comms_abort
      end if
      deallocate(e,stat=info)
      call utils_dealloc_check('internal_spline_uni_fit (density_mod.F90)',&
           'e',info)
      deallocate(d,stat=info)
      call utils_dealloc_check('internal_spline_uni_fit (density_mod.F90)',&
           'd',info)

    end subroutine internal_spline_uni_fit


    !==========================================================================
    ! Evaluate a natural cubic spline on a uniform grid assumed to start at x=0
    !==========================================================================

    real(kind=DP) function internal_spline_uni_eval(n,h,y,ypp,x)

      implicit none
      integer, intent(in) :: n
      real(kind=DP), intent(in) :: h,x,y(n),ypp(n)
      integer :: ix
      real(kind=DP), parameter :: sixth = 1.0d0 / 6.0d0
      real(kind=DP) :: a,b,pos

      pos = x / h
      ix = int(pos)
      b = pos - real(ix,kind=DP)
      a = 1.0_DP - b
      internal_spline_uni_eval = a*y(ix+1)+b*y(ix+2) + &
           ((a*a*a-a)*ypp(ix+1)+(b*b*b-b)*ypp(ix+2))*h*h*sixth

    end function internal_spline_uni_eval

  end subroutine density_plot_slice


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_radial_init

    !=================================================================!
    ! This subroutine initialises storage for radial densities        !
    !-----------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine in May 2011.                      !
    !=================================================================!

    use comms, only: pub_on_root, comms_abort
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check

    ! Local Variables
    integer :: ierr

    ! Allocate storage for all species
    allocate(radial_densities(pub_cell%num_species),stat=ierr)
    call utils_alloc_check('density_radial_init','radial_densities',ierr)
    radial_densities(:)%present = .false.

  end subroutine density_radial_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_radial_store(isp,npts,rad,den,aug_den)

    !=================================================================!
    ! This subroutine stores a radial density for a given species, to !
    ! be used later in construction of the initial guess density      !
    !-----------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine in May 2011.                      !
    !=================================================================!

    use constants, only: PI
    use rundat, only: pub_aug
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_abort

    ! Arguments
    integer, intent(in) :: isp
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: rad(npts)
    real(kind=DP), intent(in) :: den(npts)
    real(kind=DP), intent(in) :: aug_den(npts)

    ! Local Variables
    integer :: ierr
    integer :: is
    character(len=10) :: isp_str

    ! Check sanity of arguments
    if ((isp<1).or.(isp>pub_cell%num_species)) then
       write(isp_str,'(i9)') isp
       call utils_abort('Error in density_radial_store: invalid species &
            &number supplied: '//adjustl(trim(isp_str)))
    end if
    if ((npts<0).or.(npts>50000)) then
       write(isp_str,'(i9)') isp
       call utils_abort('Error in density_radial_store: invalid radial grid &
            &size for species '//adjustl(trim(isp_str)))
    end if

    ! Allocate storage for this species
    allocate(radial_densities(isp)%rad(npts),stat=ierr)
    call utils_alloc_check('density_radial_store', &
         'radial_densities(isp)%rad',ierr)
    allocate(radial_densities(isp)%den(npts,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('density_radial_store', &
         'radial_densities(isp)%den',ierr)
    if (pub_aug) then
       allocate(radial_densities(isp)%aug_den(npts,pub_cell%num_spins), &
            stat=ierr)
       call utils_alloc_check('density_radial_store', &
            'radial_densities(isp)%aug_den',ierr)
    end if


    ! Copy in arguments
    radial_densities(isp)%npts = npts
    radial_densities(isp)%rad(1:npts) = rad(1:npts)
    do is=1,pub_cell%num_spins
       radial_densities(isp)%den(1:npts,is) = den(1:npts) / 4.0_DP / PI
       if (pub_aug) radial_densities(isp)%aug_den(1:npts,is) = &
            aug_den(1:npts) / 4.0_DP / PI
    end do

    ! Set flag to indicate this density exists
    radial_densities(isp)%present = .true.

  end subroutine density_radial_store


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_radial_bcast(rad_den)

    !=================================================================!
    ! This subroutine stores a radial density for a given species, to !
    ! be used later in construction of the initial guess density      !
    !-----------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine in May 2011.                      !
    !=================================================================!

    use comms, only: comms_bcast, pub_root_node_id
    use constants, only: PI
    use rundat, only: pub_aug
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_abort

    ! Arguments
    type(RADIAL_DENSITY_TYPE), intent(inout) :: rad_den

    ! Local Variables
    integer :: ierr

    ! Broadcast size of arrays
    call comms_bcast(pub_root_node_id,rad_den%npts)

    ! Allocate storage for this species
    if (.not.rad_den%present) then
       allocate(rad_den%rad(rad_den%npts),stat=ierr)
       call utils_alloc_check('density_radial_bcast', &
            'radial_densities(isp)%rad',ierr)
       allocate(rad_den%den(rad_den%npts,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('density_radial_bcast', &
            'radial_densities(isp)%den',ierr)
       if (pub_aug) then
          allocate(rad_den%aug_den(rad_den%npts,pub_cell%num_spins),stat=ierr)
          call utils_alloc_check('density_radial_bcast', &
               'radial_densities(isp)%aug_den',ierr)
       end if
    end if

    ! Broadcast density and radial grid
    call comms_bcast(pub_root_node_id,rad_den%rad)
    call comms_bcast(pub_root_node_id,rad_den%den)
    if (pub_aug) call comms_bcast(pub_root_node_id,rad_den%aug_den)

    ! Set flag to indicate this density exists
    rad_den%present = .true.

  end subroutine density_radial_bcast


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine density_radial_exit

    !=================================================================!
    ! This subroutine deallocates storage for radial densities        !
    !-----------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine in May 2011.                      !
    !=================================================================!

    use rundat, only: pub_aug
    use simulation_cell, only: pub_cell
    use utils, only: utils_dealloc_check

    ! Local Variables
    integer :: isp
    integer :: ierr

    do isp=1,pub_cell%num_species
       if (radial_densities(isp)%present) then
          if (pub_aug) then
             deallocate(radial_densities(isp)%aug_den,stat=ierr)
             call utils_dealloc_check('density_radial_exit', &
                  'radial_densities(isp)%aug_den',ierr)
          end if
          deallocate(radial_densities(isp)%den,stat=ierr)
          call utils_dealloc_check('density_radial_exit', &
               'radial_densities(isp)%den',ierr)
          deallocate(radial_densities(isp)%rad,stat=ierr)
          call utils_dealloc_check('density_radial_exit', &
               'radial_densities(isp)%rad',ierr)
       end if
    end do

    deallocate(radial_densities,stat=ierr)
    call utils_dealloc_check('density_radial_exit', &
         'radial_densities',ierr)

  end subroutine density_radial_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module density
