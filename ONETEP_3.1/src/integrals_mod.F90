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

module integrals

  use constants, only: DP

  implicit none

  private

  public :: integrals_locpot
  public :: integrals_kinetic
  public :: integrals_grad
  public :: integrals_pos
  public :: integrals_brappd_ketfftbox
  public :: integrals_brappd_ketppd
  public :: integrals_trace_on_grid
  public :: integrals_product_on_grid

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine integrals_locpot(locpot,  &                        ! input-output
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, grid, potential_fine)
    !==========================================================================!
    ! This subroutine calculates local potential integrals between two         !
    ! sets of functions.                                                       !
    ! Result is locpot_\alpha\beta = < bra_\alpha | Vloc | ket_\beta >         !
    !==========================================================================!
    !  Arguments:                                                              !
    !    locpot     (inout)  : sparse local potential matrix elements          !
    !    bras_on_grid  (in)  : bra functions on grid in PPD format             !
    !    bra_basis (input)   : The basis type for the bras                     !
    !    kets_on_grid  (in)  : ket functions on grid in PPD format             !
    !    ket_basis (input)   : The basis type for the kets                     !
    !    potential_fine (in) : local potential on fine sim cell grid           !
    !==========================================================================!
    ! Originaly written by Chris-Kriton Skylaris in January 2001,              !
    ! to use a "pair-box".                                                     !
    ! Improved by Oswaldo Dieguez in 2001 so that it "cshifts" the "pair-box". !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".       !
    ! Rewritten by Chris-Kriton Skylaris on 23/11/2003 so that it runs on      !
    ! parallel computers.                                                      !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005                   !
    ! Modified to use planned sums system by Nicholas Hine, May 2008           !
    ! Symmetrisation options added by Nicholas Hine, July 2008.                !
    ! Modified to remove integer buffers by Nicholas Hine, June 2009           !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009                !
    ! Modified to take different function sets for bras and kets by            !
    ! Nicholas Hine on 11/11/2009.                                             !
    ! Modified to skip symmetrisation of the locpot matrix for the case of     !
    ! cross-Hamiltonians (bras/=kets) by Nicholas Hine, Oct 2010.              !
    ! Modified for possibility of performing integrals on coarse grid by       !
    ! Nicholas Hine in March 2011.                                             !
    !==========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_location_func_wrt_cell
    use cell_grid, only: GRID_INFO, pub_std_grid, pub_dbl_grid
    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use potential, only: potential_apply_to_ppd_funcs, &
         potential_input_to_workspace
    use rundat, only: pub_dbl_is_std
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: locpot
    type(FUNC_BASIS), intent(in) :: bra_basis
    real(kind=DP), intent(in) :: bras_on_grid(bra_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: ket_basis
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: kets_on_grid(ket_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in) :: potential_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12)

    ! Local Variables
    integer :: ierr
    real(kind=DP), allocatable :: potential_work(:,:,:)
    real(kind=DP), allocatable :: locpot_kets_on_grid(:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering integrals_locpot'
#endif

    ! ndmh: if the potential is on the standard grid, use that as workspace
    if (((grid%n1==pub_std_grid%n1).and.(grid%n2==pub_std_grid%n2).and. &
         (grid%n3==pub_std_grid%n3)).or.(pub_dbl_is_std)) then

       ! Start timer
       call timer_clock('integrals_locpot_coarse_grid',1)

       ! Allocate coarse-grid workspace
       allocate(potential_work(pub_std_grid%ld1,pub_std_grid%ld2,&
            pub_std_grid%max_group_slabs12),stat=ierr)
       call utils_alloc_check('integrals_locpot','potential_work',ierr)
       allocate(locpot_kets_on_grid(ket_basis%n_ppds*pub_cell%n_pts),stat=ierr)
       call utils_alloc_check('integrals_locpot','locpot_kets_on_grid',ierr)

       ! Get workspace potential from input potential
       call potential_input_to_workspace(potential_work,potential_fine, &
            pub_std_grid,grid)

       ! Copy in the ket functions
       locpot_kets_on_grid = kets_on_grid

       ! Apply the potential to the kets
       call potential_apply_to_ppd_funcs(locpot_kets_on_grid,ket_basis,&
            potential_work)

       ! Now do the overlap matrix of the bras with the locpot times the kets
       call integrals_brappd_ketppd(locpot,bras_on_grid,bra_basis, &
            locpot_kets_on_grid,ket_basis)

       ! Deallocate temporary storage
       deallocate(locpot_kets_on_grid,stat=ierr)
       call utils_dealloc_check('integrals_locpot','locpot_kets_on_grid',ierr)
       deallocate(potential_work,stat=ierr)
       call utils_dealloc_check('integrals_locpot','potential_work',ierr)

       ! Stop timer
       call timer_clock('integrals_locpot_coarse_grid',2)

    else

       ! Allocate double-grid workspace
       allocate(potential_work(pub_dbl_grid%ld1,pub_dbl_grid%ld2,&
            pub_dbl_grid%max_group_slabs12),stat=ierr)
       call utils_alloc_check('integrals_locpot','potential_work',ierr)

       ! Get workspace potential from input potential
       call potential_input_to_workspace(potential_work,potential_fine, &
            pub_dbl_grid,grid)

       ! Calculate the locpot integrals on the double grid, with interpolation
       call integrals_locpot_dbl_grid(locpot,bras_on_grid,bra_basis, &
            kets_on_grid,ket_basis,potential_work)

       ! Deallocate temporary storage
       deallocate(potential_work,stat=ierr)
       call utils_dealloc_check('integrals_locpot','potential_work',ierr)

    end if

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving integrals_locpot'
#endif

  end subroutine integrals_locpot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine integrals_locpot_dbl_grid(locpot,  &                ! input-output
       bras_on_grid, bra_basis, kets_on_grid, ket_basis, potential_dbl)

    !==========================================================================!
    ! This subroutine calculates local potential integrals between two         !
    ! sets of functions.                                                       !
    ! Result is locpot_\alpha\beta = < bra_\alpha | Vloc | ket_\beta >         !
    !==========================================================================!
    !  Arguments:                                                              !
    !    locpot     (inout)  : sparse local potential matrix elements          !
    !    bras_on_grid  (in)  : bra functions on grid in PPD format             !
    !    bra_basis (input)   : The basis type for the bras                     !
    !    kets_on_grid  (in)  : ket functions on grid in PPD format             !
    !    ket_basis (input)   : The basis type for the kets                     !
    !    potential_fine (in) : local potential on fine sim cell grid           !
    !==========================================================================!
    ! Originaly written by Chris-Kriton Skylaris in January 2001,              !
    ! to use a "pair-box".                                                     !
    ! Improved by Oswaldo Dieguez in 2001 so that it "cshifts" the "pair-box". !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".       !
    ! Rewritten by Chris-Kriton Skylaris on 23/11/2003 so that it runs on      !
    ! parallel computers.                                                      !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005                   !
    ! Modified to use planned sums system by Nicholas Hine, May 2008           !
    ! Symmetrisation options added by Nicholas Hine, July 2008.                !
    ! Modified to remove integer buffers by Nicholas Hine, June 2009           !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009                !
    ! Modified to take different function sets for bras and kets by            !
    ! Nicholas Hine on 11/11/2009.                                             !
    ! Modified to skip symmetrisation of the locpot matrix for the case of     !
    ! cross-Hamiltonians (bras/=kets) by Nicholas Hine, Oct 2010.              !
    ! Split off from integrals_locpot for clarity of above by Nicholas Hine in !
    ! March 2011.                                                              !
    !==========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_location_func_wrt_cell
    use cell_grid, only: pub_dbl_grid
    use comms, only: comms_abort, comms_barrier, comms_reduce, pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use potential, only: potential_apply_to_ngwf_batch
    use rundat, only: locpot_int_batch_size, pub_locpot_scheme
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_generate_index, sparse_index_length, &
         sparse_scale, sparse_create, sparse_transpose, sparse_axpy, &
         sparse_destroy, sparse_expand, pattern_lower, pattern_alternate
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: locpot
    type(FUNC_BASIS), intent(in) :: bra_basis
    real(kind=DP), intent(in) :: bras_on_grid(bra_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: ket_basis
    real(kind=DP), intent(in) :: kets_on_grid(ket_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in) :: potential_dbl(pub_dbl_grid%ld1, &
         pub_dbl_grid%ld2, pub_dbl_grid%max_group_slabs12)

    ! Local Variables
    real(kind=DP), allocatable, dimension(:,:,:,:) :: potential_fftbox_batch
    integer :: batch_size
    integer :: batch_count
    integer :: batch_index
    integer :: n_batches
    integer :: local_col, local_start, local_end
    integer :: ket_start(1:3), ket_cell_start(1:3)
    integer :: max_current_size ! maximum batch size ovar all nodes
    integer :: idx_len  ! pdh: sparse matrix index length
    integer :: ierr ! pdh: error flag
    integer, allocatable, dimension(:) :: locpot_idx ! pdh: sparse index
    integer, allocatable :: ket_box_start(:,:)
    real(kind=DP), allocatable, dimension(:) :: bra_on_grid_buffer
    real(kind=DP), allocatable, dimension(:) :: ket_on_bragrid_buffer
    type(SPAM3) :: locpot_transpose  ! ndmh: for symmetrising locpot
    character(len=80) :: locpot_scheme

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering integrals_locpot_dbl_grid'
#endif

    ! Start timer
    call timer_clock('integrals_locpot_dbl_grid',1)

    ! Set batch size
    batch_size = locpot_int_batch_size

    ! pdh: obtain index for locpot
    idx_len = sparse_index_length(locpot)
    allocate(locpot_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','locpot_idx',ierr)
    call sparse_generate_index(locpot_idx,locpot)

    ! Allocate workspace
    allocate(potential_fftbox_batch(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3,batch_size),stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','potential_fftbox_batch',ierr)
    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','ket_box_start',ierr)
    allocate(bra_on_grid_buffer(bra_basis%func_on_grid_buffer_size),stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','bra_on_grid_buffer',ierr)
    allocate(ket_on_bragrid_buffer(bra_basis%max_n_ppds_sphere*pub_cell%n_pts),&
         stat=ierr)
    call utils_alloc_check('integrals_locpot_dbl_grid','ket_on_bragrid_buffer',ierr)

    ! cks: number of row-steps per row-block
    n_batches = ket_basis%max_on_node / batch_size
    if (mod(ket_basis%max_on_node,batch_size) > 0) n_batches = n_batches + 1

    ! cks: ket functions are by definition placed in the centre of the fftbox
    ! cks: and stay there, or they are left where they are in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(ket_start(1),ket_start(2),ket_start(3), &
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    local_start = 1
    do batch_count=1,n_batches

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop',batch_count, &
         ' of',n_batches,' in integrals_locpot_dbl_grid'
#endif

       local_end = min(local_start+batch_size-1,ket_basis%node_num)

       ! cks: maximum size of current batch over all nodes
       max_current_size = local_end - local_start + 1
       call comms_reduce('MAX',max_current_size)

       ! cks: apply locpot operator on to kets
       call potential_apply_to_ngwf_batch(potential_fftbox_batch, &  ! output
            kets_on_grid, ket_basis, potential_dbl, batch_size, &    ! input
            local_start, local_end, max_current_size)                ! input

       ! ndmh: find fftbox start positions for each of the kets in this batch
       batch_index = 1
       do local_col=local_start,local_end

          ! Find position of ket function wrt simulation cell
          call basis_location_func_wrt_cell(ket_cell_start(1), &
               ket_cell_start(2), ket_cell_start(3), &
               ket_basis%tight_boxes(local_col))

          ket_box_start(1,batch_index) = ket_cell_start(1) - ket_start(1) + 1
          ket_box_start(2,batch_index) = ket_cell_start(2) - ket_start(2) + 1
          ket_box_start(3,batch_index) = ket_cell_start(3) - ket_start(3) + 1

          batch_index = batch_index + 1

       end do

       if (bra_basis%name/=ket_basis%name) then
          locpot_scheme = 'FULL'
       else
          locpot_scheme = pub_locpot_scheme
       end if

       ! ndmh: calculate locpot integrals in SPAM3 format
       call integrals_brappd_ketfftbox(locpot,                       & ! inout
            bras_on_grid, bra_basis,                                 & ! input
            potential_fftbox_batch, ket_box_start, batch_size,       & ! input
            local_start, local_end, idx_len, locpot_idx,             & ! input
            pub_locpot_scheme,                                       & ! input
            bra_on_grid_buffer, ket_on_bragrid_buffer)                 ! work

       local_start = local_start + batch_size

    end do

    ! Deallocate workspace
    deallocate(ket_on_bragrid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','ket_on_bragrid_buffer',ierr)
    deallocate(bra_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','bra_on_grid_buffer',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','ket_box_start',ierr)
    deallocate(potential_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','potential_fftbox_batch',ierr)
    deallocate(locpot_idx,stat=ierr)
    call utils_dealloc_check('integrals_locpot_dbl_grid','locpot_idx',ierr)

    ! ndmh: expand to full matrix
    select case (locpot_scheme)
    case ('LOWER','lower')
          call sparse_expand(locpot,PATTERN_LOWER)
    case ('ALTERNATE','alternate')
          call sparse_expand(locpot,PATTERN_ALTERNATE)
    case ('FULL','full')
          ! ndmh: only symmetrise if bra_basis==ket_basis
          if (bra_basis%name==ket_basis%name) then
             call sparse_scale(locpot,0.5_DP)
             call sparse_create(locpot_transpose, locpot)
             call sparse_transpose(locpot_transpose, locpot)
             call sparse_axpy(locpot, locpot_transpose, 1.0d0)
             call sparse_destroy(locpot_transpose)
          end if
    case ('ASYM','asym')
       if (pub_on_root) write(stdout,'(a)') &
            'WARNING: Skipping symmetrisation of locpot in integrals_locpot - '
       if (pub_on_root) write(stdout,'(a)') &
            'Resulting matrix may not be Hermitian.'
    case default
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in integrals_locpot: calculation pattern"',&
            trim(pub_locpot_scheme), '" not recognised'
       call comms_abort
    end select

    ! pdh: re-sync nodes
    call comms_barrier

    ! Stop timer
    call timer_clock('integrals_locpot_dbl_grid',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &integrals_locpot_dbl_grid'
#endif

  end subroutine integrals_locpot_dbl_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine integrals_kinetic(kinet, bras_on_grid, bra_basis, kets_on_grid, &
       ket_basis)

    !========================================================================!
    ! This subroutine calculates kinetic integrals between two function sets.!
    ! Result is kinet_\alpha\beta = < bra_\alpha | \nabla^2 | ket_\beta >    !
    !========================================================================!
    !  Arguments:                                                            !
    !    kinet      (inout) : sparse kinetic energy matrix elements          !
    !    bras_on_grid  (in) : bra functions on grid in PPD format            !
    !    bra_basis (input)  : The basis type for the bras                    !
    !    kets_on_grid  (in) : ket functions on grid in PPD format            !
    !    ket_basis (input)  : The basis type for the kets                    !
    !========================================================================!
    ! Originaly written by Chris-Kriton Skylaris in January 2001,            !
    ! to use a "pair-box".                                                   !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".     !
    ! Rewritten by Chris-Kriton Skylaris on 20/11/2003 so that it runs on    !
    ! parallel computers.                                                    !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005                 !
    ! Further modification to use parallel SPAM 2, July 2006                 !
    ! Modified to remove integer buffers by Nicholas Hine, June 2009         !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009              !
    ! Modified to take different function sets for bras and kets by          !
    ! Nicholas Hine on 11/11/2009.                                           !
    ! Modified to calculate full pattern of the kinetic matrix for the case  !
    ! of valence-conduction Hamiltonians by Laura Ratcliff, Oct 2010.        !
    !========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_location_func_wrt_cell
    use comms, only: comms_barrier, pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use kinetic, only: kinetic_apply_to_func_batch
    use rundat, only: kinetic_int_batch_size
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_generate_index, sparse_index_length, &
         sparse_expand, PATTERN_ALTERNATE
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: kinet
    type(FUNC_BASIS), intent(in) :: bra_basis
    real(kind=DP), intent(in) :: bras_on_grid(bra_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: ket_basis
    real(kind=DP), intent(in) :: kets_on_grid(ket_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    real(kind=DP), allocatable, dimension(:,:,:,:) :: kinetic_fftbox_batch
    integer :: batch_size
    integer :: batch_count
    integer :: batch_index
    integer :: n_batches
    integer :: local_col, local_start, local_end
    integer :: ket_start(1:3), ket_cell_start(1:3)
    integer :: ierr                            ! pdh: error flag
    integer :: idx_len                         ! pdh: length of sparse index
    integer, allocatable :: kinet_idx(:)       ! pdh: sparse index
    integer, allocatable :: ket_box_start(:,:)
    real(kind=DP), allocatable :: bra_on_grid_buffer(:)
    real(kind=DP), allocatable :: ket_on_bragrid_buffer(:)
    character(20) :: kinetic_scheme

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering integrals_kinetic'
#endif

    ! Start timer
    call timer_clock('integrals_kinetic',1)

    ! Obtain index for kinet
    idx_len = sparse_index_length(kinet)
    allocate(kinet_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_kinetic','kinet_idx',ierr)
    call sparse_generate_index(kinet_idx,kinet)

    ! Define shorthand variables
    batch_size = kinetic_int_batch_size

    ! Allocate workspace
    allocate(kinetic_fftbox_batch(pub_fftbox%total_ld1, pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_kinetic','kinetic_fftbox_batch',ierr)
    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_kinetic','ket_box_start',ierr)
    allocate(bra_on_grid_buffer(bra_basis%func_on_grid_buffer_size),stat=ierr)
    call utils_alloc_check('integrals_kinetic','bra_on_grid_buffer',ierr)
    allocate(ket_on_bragrid_buffer(bra_basis%max_n_ppds_sphere*pub_cell%n_pts),&
         stat=ierr)
    call utils_alloc_check('integrals_kinetic','ket_on_bragrid_buffer',ierr)

    ! cks: number of row-steps per row-block
    n_batches = ket_basis%max_on_node / batch_size
    if (mod(ket_basis%max_on_node,batch_size) > 0) n_batches = n_batches + 1

    ! cks: ket functions are by definition placed in the centre of the fftbox
    ! cks: and stay there, or they are left where they are in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(ket_start(1),ket_start(2),ket_start(3), &
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    local_start = 1
    do batch_count=1,n_batches

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ',&
         batch_count, ' of ',n_batches,' in integrals_kinetic'
#endif

       local_end = min(local_start+batch_size-1,ket_basis%node_num)

       ! cks: apply kinetic operator on to kets
       call kinetic_apply_to_func_batch(kinetic_fftbox_batch, &          ! out
            kets_on_grid, ket_basis, batch_size, local_start, local_end) ! in

       ! ndmh: find fftbox start positions for each of the kets in this batch
       batch_index = 1
       do local_col=local_start,local_end

          ! Find position of ket function wrt simulation cell
          call basis_location_func_wrt_cell(ket_cell_start(1), &
               ket_cell_start(2), ket_cell_start(3), &
               ket_basis%tight_boxes(local_col))

          ket_box_start(1,batch_index) = ket_cell_start(1) - ket_start(1) + 1
          ket_box_start(2,batch_index) = ket_cell_start(2) - ket_start(2) + 1
          ket_box_start(3,batch_index) = ket_cell_start(3) - ket_start(3) + 1

          batch_index = batch_index + 1

       end do

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a,2(a))') 'DEBUG:  ', &
            'after kinetic_apply_to_func_batch'
#endif

       ! lr408: Force calculation of the full integral if this is for a cross matrix
       if (bra_basis%name/=ket_basis%name) then
          kinetic_scheme = 'FULL'
       else
          kinetic_scheme = 'ALTERNATE'
       end if

       ! ndmh: calculate kinetic energy integrals in SPAM3 format
       call integrals_brappd_ketfftbox(kinet,                          & ! inout
            bras_on_grid, bra_basis, kinetic_fftbox_batch,             & ! input
            ket_box_start, batch_size, local_start, local_end,         & ! input
            idx_len, kinet_idx, kinetic_scheme,                        & ! input
            bra_on_grid_buffer, ket_on_bragrid_buffer)                   ! work

       local_start = local_start + batch_size

    end do

    ! ndmh: expand to full matrix
    if (kinetic_scheme=='ALTERNATE') then
       call sparse_expand(kinet,PATTERN_ALTERNATE)
    end if

    ! pdh: deallocate workspace
    deallocate(ket_on_bragrid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','ket_on_bragrid_buffer',ierr)
    deallocate(bra_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','bra_on_grid_buffer',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','ket_box_start',ierr)
    deallocate(kinetic_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','kinetic_fftbox_batch',ierr)
    deallocate(kinet_idx,stat=ierr)
    call utils_dealloc_check('integrals_kinetic','kinet_idx',ierr)

    ! pdh: sync nodes
    call comms_barrier

    ! Stop timer
    call timer_clock('integrals_kinetic',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving integrals_kinetic'
#endif

  end subroutine integrals_kinetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine integrals_grad(grad, bras_on_grid, bra_basis, kets_on_grid, &
       ket_basis)

    !=======================================================================!
    ! This subroutine calculates grad integrals between two function sets.  !
    ! Result is grad(i)_\alpha\beta = < bra_\alpha | \nabla_i | ket_\beta > !
    !=======================================================================!
    !  Arguments:                                                           !
    !    grad(3)    (inout) : sparse grad matrix elements                   !
    !    bras_on_grid  (in) : bra functions on grid in PPD format           !
    !    bra_basis (input)  : The basis type for the bras                   !
    !    kets_on_grid  (in) : ket functions on grid in PPD format           !
    !    ket_basis (input)  : The basis type for the kets                   !
    !=======================================================================!
    ! Originally written by Chris-Kriton Skylaris in January 2001,          !
    ! to use a "pair-box".                                                  !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".    !
    ! Rewritten by Chris-Kriton Skylaris on 20/11/2003 so that it runs on   !
    ! parallel computers.                                                   !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005                !
    ! Further modification to use parallel SPAM 2, July 2006                !
    ! Modified to remove integer buffers by Nicholas Hine, June 2009.       !
    ! Adapted for SPAM3 by Nicholas Hine, July 2009.                        !
    ! Modified to take different function sets for bras and kets by         !
    ! Nicholas Hine on 11/11/2009.                                          !
    !=======================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_location_func_wrt_cell
    use comms, only: comms_barrier, pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use kinetic, only: kinetic_grad_to_func_batch
    use rundat, only: kinetic_int_batch_size
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_generate_index, sparse_index_length, &
         sparse_expand, pattern_alternate
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: grad(3)
    type(FUNC_BASIS), intent(in) :: bra_basis
    real(kind=DP), intent(in) :: bras_on_grid(bra_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: ket_basis
    real(kind=DP), intent(in) :: kets_on_grid(ket_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: grad_fftbox_batch
    integer :: batch_size
    integer :: batch_count
    integer :: batch_index
    integer :: n_batches
    integer :: local_col, local_start, local_end
    integer :: ket_start(1:3), ket_cell_start(1:3)
    integer :: ierr         ! pdh: error flag
    integer :: dim          ! pdh: cartesian direction
    integer :: idx_len      ! pdh: length of sparse index
    integer, allocatable :: ket_box_start(:,:)
    integer, allocatable, dimension(:) :: grad_idx   ! pdh: sparse index
    real(kind=DP), allocatable, dimension(:) :: bra_on_grid_buffer
    real(kind=DP), allocatable, dimension(:) :: ket_on_bragrid_buffer
    character(20) :: pattern

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering integrals_grad'
#endif

    ! Start timer
    call timer_clock('integrals_grad',1)

    ! Obtain index for grad
    idx_len = sparse_index_length(grad(1))
    allocate(grad_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_grad','grad_idx',ierr)
    call sparse_generate_index(grad_idx,grad(1))

    ! Define shorthand variables
    batch_size = kinetic_int_batch_size / 3 + 1

    ! Allocate workspace
    allocate(grad_fftbox_batch(pub_fftbox%total_ld1, pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3, batch_size, 3), stat=ierr)
    call utils_alloc_check('integrals_grad','grad_fftbox_batch',ierr)
    allocate(bra_on_grid_buffer(bra_basis%func_on_grid_buffer_size),stat=ierr)
    call utils_alloc_check('integrals_grad','bra_on_grid_buffer',ierr)
    allocate(ket_on_bragrid_buffer(bra_basis%max_n_ppds_sphere*pub_cell%n_pts),&
         stat=ierr)
    call utils_alloc_check('integrals_grad','ket_on_bragrid_buffer',ierr)
    allocate(ket_box_start(3,batch_size),stat=ierr)
    call utils_alloc_check('integrals_grad','ket_box_start',ierr)

    ! cks: number of row-steps per row-block
    n_batches = ket_basis%max_on_node / batch_size
    if (mod(ket_basis%max_on_node,batch_size) > 0) n_batches = n_batches + 1

    ! cks: ket functions are by definition placed in the centre of the fftbox
    ! cks: and stay there, or they are left where they are in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(ket_start(1),ket_start(2),ket_start(3), &
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    local_start = 1
    do batch_count=1,n_batches

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ', &
         batch_count,' of ',n_batches,' in integrals_grad'
#endif

       local_end = min(local_start+batch_size-1,ket_basis%node_num)

       ! cks: apply grad operator on to kets
       call kinetic_grad_to_func_batch(grad_fftbox_batch, &               ! out
            kets_on_grid, ket_basis, batch_size, local_start, local_end)  ! in

       ! ndmh: find fftbox start positions for each of the kets in this batch
       batch_index = 1

       do local_col=local_start,local_end

          ! Find position of ket function wrt simulation cell
          call basis_location_func_wrt_cell(ket_cell_start(1), &
               ket_cell_start(2), ket_cell_start(3), &
               ket_basis%tight_boxes(local_col))

          ket_box_start(1,batch_index) = ket_cell_start(1) - ket_start(1) + 1
          ket_box_start(2,batch_index) = ket_cell_start(2) - ket_start(2) + 1
          ket_box_start(3,batch_index) = ket_cell_start(3) - ket_start(3) + 1

          batch_index = batch_index + 1

       end do

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: after &
            &kinetic_grad_to_func_batch'
#endif

       do dim=1,3

          ! ddor: FULL integral is calculated for non-square matrices
          ! lr408: changed to different basis names
          if (bra_basis%name .eq. ket_basis%name) then
             pattern='ALTERNATE'
          else
             pattern='FULL'
          end if

          ! ndmh: calculate kinetic energy integrals in SPAM3 format
          call integrals_brappd_ketfftbox(grad(dim),&              ! inout
               bras_on_grid, bra_basis, &                          ! in
               grad_fftbox_batch(:,:,:,:,dim), &                   ! in
               ket_box_start, batch_size, local_start, local_end,& ! in
               idx_len, grad_idx, pattern, &                       ! in
               bra_on_grid_buffer, ket_on_bragrid_buffer)          ! work

       end do

       local_start = local_start + batch_size

    end do

    ! pdh: expand to full matrix (which is antisymmetric)
    ! ddor: non-square matrix needs to be calculated fully, otherwise expand
    if (pattern=='ALTERNATE') then
       do dim=1,3
          call sparse_expand(grad(dim),PATTERN_ALTERNATE,.false.)
       end do
    endif

    ! pdh: deallocate workspace
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_grad','ket_box_start',ierr)
    deallocate(ket_on_bragrid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_grad','ket_on_bragrid_buffer',ierr)
    deallocate(bra_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_grad','bra_on_grid_buffer',ierr)
    deallocate(grad_fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_grad','grad_fftbox_batch',ierr)
    deallocate(grad_idx,stat=ierr)
    call utils_dealloc_check('integrals_grad','grad_idx',ierr)

    ! pdh: sync nodes
    call comms_barrier

    ! Stop timer
    call timer_clock('integrals_grad',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving integrals_grad'
#endif

  end subroutine integrals_grad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !=======================================================================!
    ! This subroutine calculates matrix elements of powers of the position  !
    ! operator (within the FFTbox) between two function sets.               !
    ! Gives rmat(i)_\alpha\beta = < bra_\alpha | r^FFT_i^pow | ket_\beta >  !
    ! where r^FFT_i = r_i - R_FFT, with R_FFT being the origin of the beta  !
    ! function's FFTbox.                                                    !
    ! This function does not account for PAW augmentation charges and the   !
    ! result must be post-processed to include them.                        !
    !=======================================================================!
    !  Arguments:                                                           !
    !    rmat(3)    (inout) : sparse position operator matrix elements      !
    !    bras_on_grid  (in) : bra functions on grid in PPD format           !
    !    bra_basis  (input) : The basis type for the bras                   !
    !    kets_on_grid  (in) : ket functions on grid in PPD format           !
    !    ket_basis  (input) : The basis type for the kets                   !
    !    order      (input) : Power of the position operator r_i            !
    !    axis   (opt input) : The axis to calculate, if only doing one.     !
    !=======================================================================!
    ! Originally written by Mark Robinson as internal_matrix_elements in    !
    ! the routine polarisation_calculate.                                   !
    ! Re-worked by Nicholas Hine in 2008 to use SPAM3 and new version of    !
    ! integrals_brappd_ketfftbox                                            !
    ! Moved to integrals_mod by Nicholas Hine on 11/11/2009.                !
    ! Changed from pattern ALTERNATE to pattern FULL, as with the former it !
    ! was impossible to correctly account for the origin shift (12/04/2010) !
    !=======================================================================!

  subroutine integrals_pos(rmat, bras_on_grid, bra_basis, kets_on_grid, &
       ket_basis,order,axis)

    use basis, only: basis_copy_function_to_box, &
         basis_ket_start_wrt_fftbox, basis_location_func_wrt_cell
    use cell_grid, only: pub_dbl_grid
    use comms, only : comms_reduce
    use constants, only : DP
    use fourier, only: fourier_interpolate, fourier_filter
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*), operator(+) !ddor
    use rundat, only: locpot_int_batch_size
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_alloc_check, utils_dealloc_check
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index !ddor

    ! Arguments
    type(SPAM3), intent(inout)  :: rmat(3)
    type(FUNC_BASIS), intent(in) :: bra_basis
    real(kind=DP), intent(in) :: bras_on_grid(bra_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: ket_basis
    real(kind=DP), intent(in) :: kets_on_grid(ket_basis%n_ppds*pub_cell%n_pts)
    ! ddor: Power of position operator for which matrix elements are calculated
    integer, intent(in) :: order
    ! ddor: Specifies if we only need one direction
    integer, optional, intent(in)        :: axis

    ! Local Variables
    type(POINT) :: r1,r2,r3
    type(POINT) :: a1,a2,a3
    integer :: i1,i2,i3
    integer :: n3,ld1,ld2,ld1_dbl,ld2_dbl
    integer :: max_current_size
    integer :: local_start,local_end
    integer :: batch,n_batches,batch_size,batch_count
    integer :: batch_index
    integer :: local_col
    integer :: ket_start(1:3), ket_cell_start(1:3)
    integer, allocatable :: ket_box_start(:,:)
    integer :: ngwf,f1,f2
    integer :: idx_len
    integer :: xyz,ierr
    ! ddor: Integers used to select axis if only one is needed
    integer :: xyz_index,xyz_count

    ! workspace
    real(kind=DP), allocatable, dimension(:,:,:,:) :: fftbox_batch
    integer,       allocatable, dimension(:)       :: idx
    real(kind=DP), allocatable, dimension(:)       :: bra_on_grid_buffer
    real(kind=DP), allocatable, dimension(:)       :: ket_on_bragrid_buffer

    ! position operator, FFTbox and tightbox
    real(kind=DP), dimension(:,:,:,:), allocatable :: r_op
    real(kind=DP), dimension(:,:,:), allocatable :: fftbox1,fftbox2
    real(kind=DP), dimension(:,:,:), allocatable :: fftbox1_dbl,fftbox2_dbl

    ! local copies of fftbox info
    n3 = pub_fftbox%total_pt3
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2
    ld1_dbl = pub_fftbox%total_ld1_dbl
    ld2_dbl = pub_fftbox%total_ld2_dbl

    ! collect batch size and determine number of batches
    batch_size = locpot_int_batch_size
    n_batches = ket_basis%max_on_node / batch_size
    if (mod(ket_basis%max_on_node,batch_size) > 0) n_batches = n_batches + 1

    ! obtain index for rmat (all 3 are the same)
    idx_len = sparse_index_length(rmat(1))
    allocate(idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_pos','idx',ierr)
    call sparse_generate_index(idx,rmat(1))

    ! allocate workspace
    allocate(fftbox_batch(ld1,ld2,n3,batch_size),stat=ierr)
    call utils_alloc_check('integrals_pos','fftbox_batch',ierr)
    allocate(bra_on_grid_buffer(bra_basis%func_on_grid_buffer_size), &
         stat=ierr)
    call utils_alloc_check('integrals_pos','bra_on_grid_buffer',ierr)
    allocate(ket_on_bragrid_buffer( &
         bra_basis%max_n_ppds_sphere*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('integrals_pos','ket_on_bragrid_buffer',ierr)

    ! allocate position operator and FFTboxes
    allocate(ket_box_start(3, batch_size), stat=ierr)
    call utils_alloc_check('integrals_pos','ket_box_start',ierr)
    allocate(fftbox1(ld1,ld2,n3), stat=ierr)
    call utils_alloc_check('integrals_pos','fftbox1',ierr)
    allocate(fftbox2(ld1,ld2,n3), stat=ierr)
    call utils_alloc_check('integrals_pos','fftbox2', ierr)
    allocate(fftbox1_dbl(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('integrals_pos','fftbox1_dbl',ierr)
    allocate(fftbox2_dbl(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('integrals_pos','fftbox2_dbl',ierr)
    allocate(r_op(3,ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('integrals_pos','r_op',ierr)

    ! calculate vectors between fine grid points
    a1 = (1.0_DP / pub_dbl_grid%n1) * pub_cell%a1
    a2 = (1.0_DP / pub_dbl_grid%n2) * pub_cell%a2
    a3 = (1.0_DP / pub_dbl_grid%n3) * pub_cell%a3

    ! ddor: construct position operator to power `order' in FFTbox
    r_op = 0.0_DP
    do i3=1,2*n3
       r3 = real(i3-1,kind=DP) * a3
       do i2=1,ld2_dbl
          r2 = r3 + real(i2-1,kind=DP) * a2
          do i1=1,ld1_dbl
             r1 = r2 + real(i1-1,kind=DP) * a1
             r_op(1,i1,i2,i3) = (r1%X)**order
             r_op(2,i1,i2,i3) = (r1%Y)**order
             r_op(3,i1,i2,i3) = (r1%Z)**order
          enddo
       enddo
    enddo

    ! cks: ket functions are by definition placed in the centre of the fftbox
    ! cks: and stay there, or they are left where they are in the simulation
    ! cks: cell (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(ket_start(1),ket_start(2),ket_start(3), &
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    ! ddor: Do one direction only if axis is present
    if (present(axis)) then
       xyz_count = 1
    else
       xyz_count = 3
    endif

    ! loop over x,y,z (or just axis)
    do xyz_index=1,xyz_count

       ! ddor: xyz is the label on the axis direction
       if (present(axis)) then
          xyz = axis
       else
          xyz = xyz_index
       endif

       ! loop over batches of NGWFs on this node to construct
       ! the matrix elements R_ab = < phi_a | r | phi_b >
       local_start = 1
       do batch=1,n_batches

          ! find last NGWF in batch
          local_end = min(local_start+batch_size-1,ket_basis%node_num)
          max_current_size = local_end - local_start
          call comms_reduce('MAX',max_current_size)

          ! place NGWFs in FFTboxes
          batch_count = 1
          do ngwf=local_start,local_start+max_current_size,2
             f1 = ngwf
             f2 = ngwf +1

             ! copy function 1 into first fftbox if required
             if (f1 <= local_end) then
                call basis_copy_function_to_box(fftbox1,ld1,ld2,n3, &
                     ket_start(1),ket_start(2),ket_start(3), &
                     ket_basis%tight_boxes(f1),kets_on_grid, &
                     ket_basis%spheres(f1))
             else
                fftbox1 = 0.0_DP
             endif

             ! copy function 2 into second fftbox if required
             if (f2 <= local_end) then
                call basis_copy_function_to_box(fftbox2,ld1,ld2,n3, &
                     ket_start(1),ket_start(2),ket_start(3), &
                     ket_basis%tight_boxes(f2),kets_on_grid, &
                     ket_basis%spheres(f2))
             else
                fftbox2 = 0.0_DP
             endif

             if (f1 <= local_end) then
                ! interpolate fftboxes to fine grid
                call fourier_interpolate(fftbox1,fftbox2,fftbox1_dbl, &
                     fftbox2_dbl)

                ! apply FFTbox position operator and increment batch_count
                fftbox1_dbl = r_op(xyz,:,:,:) * fftbox1_dbl
                if (f2 <= local_end) fftbox2_dbl = r_op(xyz,:,:,:) &
                     * fftbox2_dbl

                ! filter fftboxes to standard grid
                call fourier_filter(fftbox1_dbl,fftbox2_dbl,fftbox1,fftbox2)
             endif

             ! copy functions into batch
             fftbox_batch(:,:,:,batch_count) = fftbox1
             batch_count = batch_count +1
             if (f2 <= local_end) then
                fftbox_batch(:,:,:,batch_count) = fftbox2
                batch_count = batch_count +1
             endif
          enddo

          ! ndmh: find fftbox start positions for each of the kets in this
          ! ndmh: batch
          batch_index = 1
          do local_col=local_start,local_end

             ! Find position of ket function wrt simulation cell
             call basis_location_func_wrt_cell(ket_cell_start(1), &
                  ket_cell_start(2), ket_cell_start(3), &
                  ket_basis%tight_boxes(local_col))

             ket_box_start(1,batch_index) = ket_cell_start(1) - &
                  ket_start(1) + 1
             ket_box_start(2,batch_index) = ket_cell_start(2) - &
                  ket_start(2) + 1
             ket_box_start(3,batch_index) = ket_cell_start(3) - &
                  ket_start(3) + 1

             batch_index = batch_index + 1

          end do

          ! ndmh: calculate position operator matrix elements in SPAM3 format
          call integrals_brappd_ketfftbox(rmat(xyz),                    & ! inout
               bras_on_grid, bra_basis,                                 & ! input
               fftbox_batch, ket_box_start, batch_size,                 & ! input
               local_start, local_end, idx_len, idx, 'FULL',            & ! input
               bra_on_grid_buffer, ket_on_bragrid_buffer)                 ! work

          ! increment batch start
          local_start = local_start + batch_size

       enddo  ! end batch loop

    enddo  ! end loop over xyz

    ! deallocate FFTbox and tightbox
    deallocate(r_op,stat =ierr)
    call utils_dealloc_check('integrals_pos','r_op',ierr)
    deallocate(fftbox2_dbl,stat =ierr)
    call utils_dealloc_check('integrals_pos','fftbox2_dbl',ierr)
    deallocate(fftbox1_dbl,stat =ierr)
    call utils_dealloc_check('integrals_pos','fftbox1_dbl',ierr)
    deallocate(fftbox2,stat =ierr)
    call utils_dealloc_check('integrals_pos','fftbox2',ierr)
    deallocate(fftbox1,stat =ierr)
    call utils_dealloc_check('integrals_pos','fftbox1',ierr)
    deallocate(ket_box_start,stat=ierr)
    call utils_dealloc_check('integrals_pos','ket_box_start',ierr)

    ! deallocate workspace
    deallocate(ket_on_bragrid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_pos','ket_on_bragrid_buffer',ierr)
    deallocate(bra_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_pos','bra_on_grid_buffer',ierr)
    deallocate(fftbox_batch,stat=ierr)
    call utils_dealloc_check('integrals_pos','fftbox_batch',ierr)

    ! deallocate index
    deallocate(idx,stat=ierr)
    call utils_dealloc_check('integrals_pos','idx',ierr)

  end subroutine integrals_pos


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine integrals_brappd_ketfftbox(bracket,                  &    ! inout
       bras_on_grid, bra_basis,                                   &    ! input
       ket_fftbox_batch, ket_box_start, batch_size,               &    ! input
       local_start, local_end, idx_len, bracket_idx, pattern,     &    ! input
       bra_on_grid_buffer, ket_on_bragrid_buffer)                      ! work

    !=========================================================================!
    ! This subroutine computes a "bracket" <fa|fb> in sparse matrix (SPAM3)   !
    ! storage form. In general, the "bras" fa, which are represented on the   !
    ! coarse grid in ppd-indexed form, are a different set of functions from  !
    ! the "kets" fb, which are in fftboxes. The matrix is divided into blocks !
    ! calculated by different nodes. In general the ket functions belong      !
    ! to the local node while the bra functions are received from other nodes.!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! bracket (input/output) : SPAM3 structure for matrix elements            !
    ! bracket_idx (input)    : index for above matrix                         !
    ! bras_on_grid (input)   : Bra functions in ppd representation            !
    ! ket_fftbox_batch (input) : Ket functions in fftbox representation       !
    ! local_start (input)    :                                                !
    ! local_end (input)      :                                                !
    ! batch_size (input)     :                                                !
    ! bra_basis (input)      : The basis type for the bra functions           !
    ! scheme (input)         : Which terms of the matrix to calculate         !
    !                          scheme = 'FULL' calculates all terms           !
    !-------------------------------------------------------------------------!
    ! Integral subroutines were originally written by Chris-Kriton Skylaris   !
    ! in 2000.                                                                !
    ! Arash Mostofi rewrote them so that they use a "triple box" and          !
    ! complex-to-complex FFTs.                                                !
    ! This subroutine was written by Chris-Kriton Skylaris on 18/9/2003       !
    ! and is capable of running on parallel computers with an                 !
    ! arbitary number of processors.                                          !
    ! Modifications and speed-ups by Peter Haynes 16/03/2005.                 !
    ! Modifications and speed-ups by Chris-Kriton Skylaris 08/09/2005.        !
    ! Further modification to use parallel SPAM 2, July 2006                  !
    ! Modifications to symmetrise the locpot matrix consistently between      !
    ! upper and lower triangles by Nicholas Hine, 28/04/2008                  !
    ! Modified to use planned exection by Nicholas Hine, 14/05/2008.          !
    ! Modified for SPAM3, function basis type and request-based comms by      !
    ! Nicholas Hine in June 2009.                                             !
    !=========================================================================!

    use simulation_cell, only: pub_cell, pub_fftbox
    use basis, only: basis_extract_function_from_box, basis_find_function_wrt_box
    use comms, only: comms_barrier, comms_free, comms_reduce, pub_my_node_id, &
         pub_total_num_nodes
    use function_basis, only: FUNC_BASIS, function_basis_init_requests, &
         function_basis_request, function_basis_await_requests, &
         function_basis_respond_to_reqs, function_basis_recv, &
         function_basis_batch_row_plan, pub_buffer_sphere
    use sparse, only: SPAM3, sparse_put_element, sparse_node_of_elem, &
         sparse_first_elem_on_node
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in)        :: batch_size, local_start, local_end
    type(SPAM3), intent(inout) :: bracket
    type(FUNC_BASIS), intent(in) :: bra_basis
    integer, intent(in)        :: idx_len
    integer, intent(in)        :: bracket_idx(idx_len)
    real(kind=DP), intent(in)  :: bras_on_grid(bra_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in)  :: ket_fftbox_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3, batch_size)
    integer, intent(in)        :: ket_box_start(3, batch_size)
    character(*), intent(in)   :: pattern
    real(kind=DP), intent(out) :: bra_on_grid_buffer( &
         bra_basis%func_on_grid_buffer_size)
    real(kind=DP), intent(out) :: ket_on_bragrid_buffer( &
         bra_basis%max_n_ppds_sphere*pub_cell%n_pts)

    ! internal declarations
    integer :: col, recv_row
    integer :: local_col, local_row
    integer :: recv_node
    integer :: bra_start
    integer :: bra_start1, bra_start2, bra_start3
    integer :: batch_count
    integer :: dot_npts
    logical :: remote_row                ! pdh: flag to avoid pointers
    integer :: ipt                       ! pdh: point loop counter
    integer :: n_bra_ppds
    real(kind=DP) :: bracket_el          ! pdh: matrix element
    integer :: ierr                      ! ndmh: error flag

    ! ndmh: variables for planned execution
    integer, parameter :: lookahead = 1  ! how many plan_steps in advance to send
    integer :: req_row
    integer :: plan_steps                ! number of steps in plan on this node
    integer :: plan_step                 ! counter for current step in plan
    integer :: prev_req_row              ! Last row function requested
    integer :: prev_recv_row             ! Last row function received
    integer :: req_node
    logical, allocatable :: reqs_out(:)
    logical, allocatable :: reqs_in(:)
    integer, allocatable :: plan(:,:)    ! Plan for this node
    integer, allocatable :: func_requests(:)
    integer, allocatable :: request_handles(:)


    ! cks: synchronize PEs
    call comms_barrier

    ! ndmh: allocate workspace for plan
    plan_steps = bra_basis%num*batch_size
    allocate(plan(2,plan_steps),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketfftbox','plan',ierr)
    allocate(func_requests(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketfftbox','func_requests',ierr)
    allocate(request_handles(0:pub_total_num_nodes+2),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketfftbox','request_handles',ierr)
    allocate(reqs_in(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketfftbox','reqs_in',ierr)
    allocate(reqs_out(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketfftbox','reqs_out',ierr)

    ! ndmh: create a plan from the matrix index
    call function_basis_batch_row_plan(plan_steps,plan,idx_len,bracket_idx, &
         bracket,local_start,local_end,pattern,reqs_in)

    ! ndmh: initializations
    prev_recv_row = -1
    prev_req_row = -1
    local_row = -1
    remote_row = .false.

    ! ndmh: initialisation of send request receive operations
    call function_basis_init_requests(func_requests,request_handles,reqs_in,reqs_out)

    do plan_step=2-lookahead,plan_steps

       ! ndmh: check if we need to request an function to be sent to this node
       if (plan_step+lookahead<=plan_steps) then
          req_row = plan(1,plan_step+lookahead)
          req_node = sparse_node_of_elem(req_row,bracket,'R')
          ! ndmh: if this bra is not local to this node and we do not already
          ! ndmh: have it, then send a request for it to req_node
          if ((req_node /= pub_my_node_id) .and. (req_row /= prev_req_row)) then
             call function_basis_request(req_node,req_row)
             prev_req_row = req_row
          end if
       end if

       ! ndmh: respond to any send requests made of this node
       call function_basis_respond_to_reqs(func_requests,request_handles, &
            bra_basis,bras_on_grid)

       ! ndmh: cycle if still pre-start
       if (plan_step<1) cycle

       ! ndmh: find the row and col of the matrix element to calculate
       recv_row = plan(1,plan_step)
       col = plan(2,plan_step)

       ! ndmh: find the node the row function is stored on
       recv_node = sparse_node_of_elem(recv_row,bracket,'R')

       ! ndmh: receive or copy the ppd list and offset into the buffer sphere
       if (prev_recv_row /= recv_row) then
          if (recv_node==pub_my_node_id) then
             local_row = recv_row - &
                  sparse_first_elem_on_node(pub_my_node_id,bracket,'R') + 1
             n_bra_ppds = bra_basis%spheres(local_row)%n_ppds_sphere
             pub_buffer_sphere%n_ppds_sphere = n_bra_ppds
             pub_buffer_sphere%ppd_list(:,1:n_bra_ppds) = &
                  bra_basis%spheres(local_row)%ppd_list(:,1:n_bra_ppds)
             pub_buffer_sphere%offset = bra_basis%spheres(local_row)%offset
             remote_row = .false.
          else
             call function_basis_recv(recv_node,recv_row,func_requests, &
                  request_handles,pub_buffer_sphere,bra_on_grid_buffer, &
                  bra_basis,bras_on_grid)
             remote_row = .true.
          end if
          prev_recv_row = recv_row
       end if

       ! Find position in batch of col and index on this node
       local_col = col - &
            sparse_first_elem_on_node(pub_my_node_id,bracket,'C') + 1
       batch_count = local_col - local_start + 1

       ! Find position of bra tightbox start wrt ket fftbox
       call basis_find_function_wrt_box(bra_start1, bra_start2, bra_start3, &
            ket_box_start(1,batch_count), ket_box_start(2,batch_count), &
            ket_box_start(3,batch_count), bra_basis%all_tbs(recv_row))

       ! cks: extract ppds belonging to bra function from ket fftbox
       call basis_extract_function_from_box(ket_on_bragrid_buffer(:), &
            pub_fftbox%total_ld1,pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            ket_fftbox_batch(:,:,:,batch_count), pub_buffer_sphere, &
            bra_basis%all_tbs(recv_row), bra_start1, bra_start2, bra_start3, 1)

       ! cks: ddot ppds - these bra and ket representations have the same ppds
       dot_npts = pub_buffer_sphere%n_ppds_sphere * pub_cell%n_pts
       bracket_el = 0.0_DP
       bra_start = pub_buffer_sphere%offset
       if (remote_row) then
          do ipt=0,dot_npts-1
             bracket_el = bracket_el + &
                  bra_on_grid_buffer(bra_start+ipt) * &
                  ket_on_bragrid_buffer(1+ipt)
          end do
       else
          do ipt=0,dot_npts-1
             bracket_el = bracket_el + &
                  bras_on_grid(bra_start+ipt) * &
                  ket_on_bragrid_buffer(1+ipt)
          end do
       end if

       ! cks: scale with grid point weight
       bracket_el = pub_cell%weight * bracket_el

       ! ndmh: deposit in matrix
       call sparse_put_element(bracket_el,bracket,recv_row,col)

    end do

    ! ndmh: send signal indicating completion to all nodes
    do req_node=0,pub_total_num_nodes-1
       if (reqs_in(req_node)) then
          call function_basis_request(req_node,-2000)
       end if
    end do

    call function_basis_await_requests(func_requests,request_handles, &
         bra_basis,bras_on_grid)

    ! ndmh: synchronise again so that non-blocking sends are completed and all
    ! ndmh: handles restored before deallocation
    call comms_free
    call comms_barrier

    ! ndmh: deallocate local workspace
    deallocate(reqs_out,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketfftbox','reqs_out',ierr)
    deallocate(reqs_in,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketfftbox','reqs_in',ierr)
    deallocate(request_handles,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketfftbox','request_handles',ierr)
    deallocate(func_requests,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketfftbox','func_requests',ierr)
    deallocate(plan,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketfftbox','plan',ierr)


  end subroutine integrals_brappd_ketfftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine integrals_brappd_ketppd(bracket,  &            ! inout
       bras_on_grid, bra_basis, kets_on_grid, ket_basis)     ! input

    !=========================================================================!
    ! This subroutine computes a "bracket" <fa|fb> in sparse matrix (SPAM3)   !
    ! storage form. In general, the "bras" fa can be a different set of       !
    ! functions from the "kets" fb, but both are represented on the same grid !
    ! in ppd-indexed form. The matrix is divided into blocks calculated by    !
    ! different nodes. In general the ket functions belong to the local node  !
    ! while the bra functions are received from other nodes.                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! bracket (input/output) : SPAM3 structure for matrix elements            !
    ! bras_on_grid (input)   : Bra functions in ppd representation            !
    ! kets_on_grid (input)   : Ket functions in ppd representation            !
    ! bra_basis (input)      : The basis type for the bra functions           !
    ! ket_basis (input)      : The basis type for the ket functions           !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2009 based on integrals_brackets,      !
    ! originally written by Chris-Kriton Skylaris in 2000 and modified by     !
    ! Oswaldo Dieguez, Arash Mostofi, Peter Haynes and Nicholas Hine in       !
    ! 2003 to 2009.                                                           !
    !=========================================================================!

    use simulation_cell, only: pub_cell
    use comms, only: comms_barrier, comms_free, comms_reduce, pub_my_node_id, &
         pub_total_num_nodes
    use function_basis, only: FUNC_BASIS, function_basis_init_requests, &
         function_basis_request, function_basis_await_requests, &
         function_basis_respond_to_reqs, function_basis_recv, &
         function_basis_batch_row_plan, pub_buffer_sphere
    use sparse, only: SPAM3, sparse_put_element, sparse_node_of_elem, &
         sparse_first_elem_on_node, sparse_num_elems_on_node, &
         sparse_node_num_element, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: bracket
    type(FUNC_BASIS), intent(in) :: bra_basis
    type(FUNC_BASIS), intent(in) :: ket_basis
    real(kind=DP), intent(in) :: bras_on_grid(bra_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in) :: kets_on_grid(ket_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    integer :: col, recv_row
    integer :: local_col, local_row
    integer :: recv_node
    integer :: n_ppds
    logical :: remote_row                ! pdh: flag to avoid pointers
    real(kind=DP) :: bracket_el          ! pdh: matrix element
    integer :: ierr                      ! ndmh: error flag
    integer :: idx_len
    integer,allocatable :: bracket_idx(:)
    real(kind=DP),allocatable :: bra_on_grid_buffer(:)
    ! ndmh: variables for planned execution
    integer, parameter :: lookahead = 1  ! how many plan_steps in advance to send
    integer :: req_row
    integer :: plan_steps                ! number of steps in plan on this node
    integer :: plan_step                 ! counter for current step in plan
    integer :: prev_req_row              ! Last row function requested
    integer :: prev_recv_row             ! Last row function received
    integer :: req_node
    logical, allocatable :: reqs_out(:)
    logical, allocatable :: reqs_in(:)
    integer, allocatable :: plan(:,:)    ! Plan for this node
    integer, allocatable :: func_requests(:)
    integer, allocatable :: request_handles(:)

    call timer_clock('integrals_brappd_ketppd',1)

    ! cks: synchronize PEs
    call comms_barrier

    ! ndmh: allocate workspace
    plan_steps = sparse_node_num_element(bracket)
    allocate(plan(2,plan_steps),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketppd','plan',ierr)
    idx_len = sparse_index_length(bracket)
    allocate(bracket_idx(idx_len),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketppd','bracket_idx',ierr)
    allocate(func_requests(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketppd','func_requests',ierr)
    allocate(request_handles(0:pub_total_num_nodes+2),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketppd','request_handles',ierr)
    allocate(reqs_in(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketppd','reqs_in',ierr)
    allocate(reqs_out(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('integrals_brappd_ketppd','reqs_out',ierr)
    allocate(bra_on_grid_buffer(bra_basis%max_n_ppds_sphere*pub_cell%n_pts), &
         stat=ierr)
    call utils_alloc_check('integrals_brappd_ketppd','bra_on_grid_buffer',ierr)

    ! ndmh: get the matrix index
    call sparse_generate_index(bracket_idx,bracket)

    ! ndmh: create a plan from the matrix index
    local_col = sparse_num_elems_on_node(pub_my_node_id,bracket,'C')
    call function_basis_batch_row_plan(plan_steps,plan,idx_len,bracket_idx, &
         bracket,1,local_col,'FULL',reqs_in)

    ! ndmh: initializations
    prev_recv_row = -1
    prev_req_row = -1
    local_row = -1
    remote_row = .false.

    ! ndmh: initialisation of send request receive operations
    call function_basis_init_requests(func_requests,request_handles,reqs_in, &
         reqs_out)

    do plan_step=1,plan_steps

       ! ndmh: check if we need to request a function to be sent to this node
       if (plan_step+lookahead<=plan_steps) then
          req_row = plan(1,plan_step+lookahead)
          req_node = sparse_node_of_elem(req_row,bracket,'R')
          ! ndmh: if this bra is not local to this node and we do not already
          ! ndmh: have it, then send a request for it to req_node
          if ((req_node /= pub_my_node_id) .and. (req_row /= prev_req_row)) then
             call function_basis_request(req_node,req_row)
             prev_req_row = req_row
          end if
       end if

       ! ndmh: respond to any send requests made of this node
       call function_basis_respond_to_reqs(func_requests,request_handles, &
            bra_basis,bras_on_grid)

       ! ndmh: find the row and col of the matrix element to calculate
       recv_row = plan(1,plan_step)
       col = plan(2,plan_step)

       ! ndmh: find the node the row function is stored on
       recv_node = sparse_node_of_elem(recv_row,bracket,'R')

       ! ndmh: receive or copy the ppd list and offset into the buffer sphere
       if (prev_recv_row /= recv_row) then
          if (recv_node==pub_my_node_id) then
             local_row = recv_row - &
                  sparse_first_elem_on_node(pub_my_node_id,bracket,'R') + 1
             n_ppds = bra_basis%spheres(local_row)%n_ppds_sphere
             pub_buffer_sphere%n_ppds_sphere = &
                  bra_basis%spheres(local_row)%n_ppds_sphere
             pub_buffer_sphere%ppd_list(:,1:n_ppds) = &
                  bra_basis%spheres(local_row)%ppd_list(:,1:n_ppds)
             pub_buffer_sphere%offset = bra_basis%spheres(local_row)%offset
             remote_row = .false.
          else
             call function_basis_recv(recv_node,recv_row,func_requests, &
                  request_handles,pub_buffer_sphere,bra_on_grid_buffer, &
                  bra_basis,bras_on_grid)
             remote_row = .true.
          end if
          prev_recv_row = recv_row
       end if

       ! Find local index of col on this node
       local_col = col - &
            sparse_first_elem_on_node(pub_my_node_id,bracket,'C') + 1

       ! ndmh: calculate bracket matrix element
       if (remote_row) then
          bracket_el = &
               internal_bra_dot_ket_ppds(bra_on_grid_buffer, &
               kets_on_grid, pub_buffer_sphere%n_ppds_sphere, ket_basis%n_ppds, &
               pub_buffer_sphere, ket_basis%spheres(local_col))
       else
          bracket_el = &
               internal_bra_dot_ket_ppds(bras_on_grid, &
               kets_on_grid, bra_basis%n_ppds, ket_basis%n_ppds, &
               bra_basis%spheres(local_row), ket_basis%spheres(local_col))
       end if

       ! cks: scale with grid point weight
       bracket_el = pub_cell%weight * bracket_el

       ! ndmh: deposit in matrix
       call sparse_put_element(bracket_el,bracket,recv_row,col)

    end do

    ! ndmh: send signal indicating completion to all nodes
    do req_node=0,pub_total_num_nodes-1
       if (reqs_in(req_node)) then
          call function_basis_request(req_node,-2000)
       end if
    end do

    ! ndmh: wait until completion messages have been received from all
    ! ndmh: other nodes with which communication is required.
    call function_basis_await_requests(func_requests,request_handles, &
         bra_basis,bras_on_grid)

    ! ndmh: synchronise again so that non-blocking sends are completed and all
    ! ndmh: handles restored before deallocation
    call comms_free
    call comms_barrier

    ! ndmh: deallocate workspace
    deallocate(bra_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketppd','bra_on_grid_buffer', &
         ierr)
    deallocate(reqs_out,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketppd','reqs_out',ierr)
    deallocate(reqs_in,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketppd','reqs_in',ierr)
    deallocate(request_handles,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketppd','request_handles',ierr)
    deallocate(func_requests,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketppd','func_requests',ierr)
    deallocate(bracket_idx,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketppd','bracket_idx',ierr)
    deallocate(plan,stat=ierr)
    call utils_dealloc_check('integrals_brappd_ketppd','plan',ierr)

    call timer_clock('integrals_brappd_ketppd',2)

contains

    real(kind=DP) function internal_bra_dot_ket_ppds(bra_on_grid, ket_on_grid, &
         n_bra_ppds, n_ket_ppds, bra_sphere, ket_sphere)

      !==================================================================!
      ! This function calculates a <fa|fb> matrix element as the dot     !
      ! between the points in the ppds in common between the fa function !
      ! (the "bra") and the fb function (the "ket").                     !
      !------------------------------------------------------------------!
      ! Arguments:                                                       !
      ! bra_on_grid (input): The bra functions in ppd representation     !
      ! ket_on_grid (input): The ket functions in ppd representation     !
      ! n_bra_ppds (input) : Number of ppds of bra_on_grid               !
      ! n_ket_ppds (input) : Number of ppds of ket_on_grid               !
      ! bra_sphere (input) : Sphere for bra function                     !
      ! ket_sphere (input) : Sphere for ket function                     !
      !------------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 5/12/2003                    !
      ! Updated by Nicholas Hine in December 2009 to take advantage of   !
      ! the fact that ppds are listed in ascending order, so there is no !
      ! need for a full double loop over all ppds of the bras and kets.  !
      !==================================================================!

      use basis, only: SPHERE
      use constants, only: DP
      use simulation_cell, only: pub_cell

      implicit none

      ! Arguments
      integer, intent(in) :: n_bra_ppds, n_ket_ppds
      real(kind=DP), intent(in) :: bra_on_grid(n_bra_ppds*pub_cell%n_pts)
      real(kind=DP), intent(in) :: ket_on_grid(n_ket_ppds*pub_cell%n_pts)
      type(SPHERE), intent(in)  :: bra_sphere
      type(SPHERE), intent(in)  :: ket_sphere

      ! Local Variables
      integer :: ibra, iket              ! index in bra and ket ppd lists
      integer :: bra_ppd, ket_ppd        ! ppd index of bra and ket ppds
      integer :: bra_start, ket_start    ! index of first point of ppd
      integer :: i                       ! point loop counter
      real(kind=DP) :: ppd_sum

      ppd_sum = 0.0_DP

      ! Start at the beginning of the ket ppd list
      iket = 1
      ket_ppd = ket_sphere%ppd_list(1,iket)

      ! cks: accumulate bra-ket integral
      do ibra=1,bra_sphere%n_ppds_sphere
         ! ndmh: find ppd number and start position for ibra
         bra_ppd = bra_sphere%ppd_list(1,ibra)
         bra_start = bra_sphere%offset + (ibra-1)*pub_cell%n_pts
         do
            ! ndmh: keep moving on while ket_ppd is less than bra_ppd
            if (ket_ppd < bra_ppd) then
               iket = iket + 1
               if (iket > ket_sphere%n_ppds_sphere) exit
               ket_ppd = ket_sphere%ppd_list(1,iket)
               cycle
            end if

            ! ndmh: ppd numbers match, so add product of these ppds
            if (ket_ppd == bra_ppd) then
               ket_start = ket_sphere%offset + (iket-1)*pub_cell%n_pts
               do i=0,pub_cell%n_pts-1
                  ppd_sum = ppd_sum + bra_on_grid(bra_start+i) * &
                       ket_on_grid(ket_start+i)
               end do
            end if

            ! ndmh: move on to next bra ppd as soon as we have found a match
            ! ndmh: or moved beyond bra_ppd
            if (ket_ppd >= bra_ppd) exit

         end do
      end do

      internal_bra_dot_ket_ppds = ppd_sum

    end function internal_bra_dot_ket_ppds

  end subroutine integrals_brappd_ketppd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function integrals_trace_on_grid(integrand,grid)

    !======================================================================!
    ! Returns the integral of a function on the passed grid.               !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000. Modified by Peter D. Haynes!
    ! on 1/7/2004 for Fourier parallelisation.                             !
    ! Modified to use pub_cell by Quintin Hill on 15/10/2008.              !
    !======================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: integrand(grid%ld1, &
         grid%ld2, grid%max_slabs12)

    ! Local variables
    integer :: i1,i2,islab12
    real(kind=DP) :: integral

    ! Sum over real-space on this node
    integral = 0.0_DP
    do islab12=1,grid%num_my_slabs12
       do i2=1,grid%n2
          do i1=1,grid%n1
             integral = integral + integrand(i1,i2,islab12)
          end do
       end do
    end do

    ! Sum over all nodes
    call comms_reduce('SUM',integral)

    ! Multiply by weight per grid point
    integrals_trace_on_grid = integral * grid%weight

  end function integrals_trace_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function integrals_product_on_grid(grid, x1, x2, m1, m2, m3)
    !======================================================================!
    ! Returns the integral of a product of two functions, x1 and x2 on the !
    ! passed grid. If x2 is omitted, integrates x1 only. Optional arguments!
    ! m1, m2 and m3 can be used to narrow down the integration to a subset !
    ! of the grid.                                                         !
    !----------------------------------------------------------------------!
    ! Arguments:                                                           !
    !   grid (input): The grid on which the integration takes place.       !
    !   x1, x2 (input): The quantities to integrate. x2 is optional.       !
    !   m1, m2, m3 (input, optional): Values to override grid%n1, grid%n2  !
    !                                 and grid%num_my_slabs12 with.        !
    ! Return value: Integral of x1(r)*x2(r) on the grid.                   !
    !----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in June 2010.                              !
    !======================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
    use utils, only: utils_assert

    ! jd: Arguments
    type(GRID_INFO), intent(in)         :: grid
    real(kind=DP), intent(in)           :: x1(grid%ld1, grid%ld2, &
         grid%max_slabs12)
    real(kind=DP), intent(in), optional :: x2(grid%ld1, grid%ld2, &
         grid%max_slabs12)
    integer, intent(in), optional :: m1, m2, m3

    ! jd: Internal variables
    integer :: i1, i2, i3
    real(kind=DP) :: integral
    integer :: p1, p2, p3

    !------------------------------------------------------------------------

    ! jd: Default upper bounds of integration
    p1 = grid%n1
    p2 = grid%n2
    p3 = grid%num_my_slabs12

    ! jd: Override them, if the user desires
    if(present(m1)) p1 = m1
    if(present(m2)) p2 = m2
    if(present(m3)) p3 = m3

    call utils_assert(p1>0 .and. p2>0 .and. p3>=0,'Bad integration bounds in &
         &integrals_product_on_grid')

    integral = 0.0_DP
    ! jd: The following duplicates code in order not to check the
    !     if(present) condition repeatedly in a tight loop
    if(present(x2)) then
       do i3=1, p3
          do i2=1, p2
             do i1=1, p1
                integral = integral + x1(i1,i2,i3) * x2(i1,i2,i3)
             end do
          end do
       end do
    else
       do i3=1, p3
          do i2=1, p2
             do i1=1, p1
                integral = integral + x1(i1,i2,i3)
             end do
          end do
       end do
    end if

    call comms_reduce('SUM',integral)

    ! jd: NB, the weight remains unchanged even if integration bounds overridden
    integrals_product_on_grid = integral * grid%weight

  end function integrals_product_on_grid
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module integrals

