! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by                              !
!                                                                             !
!   Nicholas D.M. Hine                                                        !
!                                                                             !
!   in December 2010.                                                         !
!                                                                             !
!   Based largely on previous code by                                         !
!                                                                             !
!   Chris-Kriton Skylaris, Nicholas Hine, Arash A. Mostofi, Peter D. Haynes   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cell_grid

  use constants, only : DP, stdout, VERBOSE
  use geometry, only: POINT

  implicit none

  private

  ! Public subroutines
  public :: cell_grid_init
  public :: cell_grid_exit
  public :: cell_grid_distribute
  public :: cell_grid_deposit_box
  public :: cell_grid_extract_box
  public :: cell_grid_real_pt
  public :: cell_grid_recip_pt
  public :: cell_grid_box_start_wrt_atom

  ! Data structure defining a 3D grid, parallelised over '12' slabs in real
  ! space and '23' slabs in reciprocal space
  type GRID_INFO
     ! Sizes of the grid and leading dimensions of the arrays
     integer :: n1
     integer :: ld1
     integer :: n2
     integer :: ld2
     integer :: n3
     integer :: ld3
     ! Integration weight per grid point
     real(kind=DP) :: weight
     ! ndmh: lattice vectors divided by number of grid points along
     ! nmdh: that lattice vector
     type(POINT) :: da1, da2, da3
     ! ndmh: reciprocal lattice vectors for the cell this grid spans
     type(POINT) :: b1, b2, b3
     ! Information about parallel FFT:
     !   the parallel FFT is done in two stages - in real-space the data is
     !   arranged in slabs dual to the '3' axis. After a 1D FFT along the
     !   '1' axis it is then transposed to slabs dual to the '1' axis. A
     !   2D FFT in the '23' plane leaves the reciprocal space data with the
     !   fastest running index along the '3' axis.
     ! In real-space the '12' slabs are indexed by their position along the
     !   '3' axis: from 1 to grid%n3
     ! Maximum number of '12' real-space slabs held by any node
     integer :: max_slabs12
     ! Maximum number of '12' real-space slabs held by any group of nodes
     integer :: max_group_slabs12
     ! Number of '12' real-space slabs held by this node
     integer :: num_my_slabs12
     ! Number of '12' real-space slabs held by this group of nodes
     integer :: num_group_slabs12
     ! Index of the first '12' real-space slab of this node in the group slabs
     integer :: my_first_slab12_in_group
     ! Index of the last '12' real-space slab of this node in the group slabs
     integer :: my_last_slab12_in_group
     ! Index of first '12' real-space slab on each node
     integer, pointer, dimension(:) :: first_slab12
     ! Index of last '12' real-space slab on each node
     integer, pointer, dimension(:) :: last_slab12
     ! List of nodes on which all the '12' slabs are held
     integer, pointer, dimension(:) :: node_slab12
     ! List of displacements of the data for each node in the group within array
     integer, pointer, dimension(:) :: group_displs_slabs12
     ! List of lengths of data of 12 slabs for each node in the group
     integer, pointer, dimension(:) :: group_lengths_slabs12

     ! In reciprocal-space the '23' slabs are indexed by their position along the
     !   '1' axis: from 1 to n1
     ! Maximum number of '23' real-space slabs held by any node
     integer :: max_slabs23
     ! Number of '23' real-space slabs held by this node
     integer :: num_slabs23
     ! Index of first '23' real-space slab on each node
     integer, pointer, dimension(:) :: first_slab23
     ! Index of last '23' real-space slab on each node
     integer, pointer, dimension(:) :: last_slab23
     ! List of nodes on which all the '23' slabs are held
     integer, pointer, dimension(:) :: node_slab23
     ! Reciprocal-space Cartesian coordinates of grid-points on this node:
     !   first index is Cartesian component 1-3, magnitude in 4
     !   second index is position along '3' axis
     !   third index is position along '2' axis
     !   fourth index is '23' slab index on this node
     real(kind=DP), pointer, dimension(:,:,:) :: coulomb_recip
     ! Flag to show that FFT distribution has occurred
     logical :: distributed
     ! Index of the entry for this grid in the array of parallel_fft3d_info
     ! structures in fourier_mod (needed to perform whole-cell FFTs)
     integer :: fft_index
  end type GRID_INFO

  public :: GRID_INFO

  type(GRID_INFO), public :: pub_std_grid
  type(GRID_INFO), public :: pub_dbl_grid
  type(GRID_INFO), public :: pub_fine_grid

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine cell_grid_exit(grid)

    !=========================================================================!
    ! This subroutine frees up any memory used by a whole-cell grid structure.!
    !-------------------------------------------------------------------------!
    ! Arguments:							      !
    !   grid: GRID_INFO type defining the grid to be deallocated.             !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 22/06/10.                                     !
    ! Based on code by Peter Haynes and Arash Mostofi (2003-2004).            !
    !=========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid

    ! Local Variables
    integer :: ierr         ! Error flag

    ! ndmh: Deallocate grid%group_lengths_slabs12 on this node
    if (associated(grid%group_lengths_slabs12)) then
       deallocate(grid%group_lengths_slabs12,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%group_lengths_slabs12',ierr)
    end if

    ! ndmh: Deallocate grid%group_displs_slabs12 on this node
    if (associated(grid%group_displs_slabs12)) then
       deallocate(grid%group_displs_slabs12,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%group_displs_slabs12',ierr)
    end if

    ! ndmh: Deallocate grid%node_slab12 on this node
    if (associated(grid%node_slab12)) then
       deallocate(grid%node_slab12,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%node_slab12',ierr)
    end if

    ! ndmh: Deallocate grid%node_slab23 on this node
    if (associated(grid%node_slab23)) then
       deallocate(grid%node_slab23,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%node_slab23',ierr)
    end if

    ! ndmh: Deallocate grid%first_slab12 on this node
    if (associated(grid%first_slab12)) then
       deallocate(grid%first_slab12,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%first_slab12',ierr)
    end if

    ! ndmh: Deallocate grid%last_slab12 on this node
    if (associated(grid%last_slab12)) then
       deallocate(grid%last_slab12,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%last_slab12',ierr)
    end if

    ! ndmh: Deallocate grid%first_slab23 on this node
    if (associated(grid%first_slab23)) then
       deallocate(grid%first_slab23,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%first_slab23',ierr)
    end if

    ! ndmh: Deallocate grid%last_slab23 on this node
    if (associated(grid%last_slab23)) then
       deallocate(grid%last_slab23,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%last_slab23',ierr)
    end if

    ! ndmh: Deallocate grid%coulomb_recip on this node
    if (associated(grid%coulomb_recip)) then
       deallocate(grid%coulomb_recip,stat=ierr)
       call utils_dealloc_check('cell_grid_exit', &
            'grid%coulomb_recip',ierr)
    end if

  end subroutine cell_grid_exit


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  subroutine cell_grid_init(grid,cell,n1,n2,n3,ld1,ld2,ld3)

    !=========================================================================!
    ! This subroutine frees up any memory used by a whole-cell grid structure.!
    !-------------------------------------------------------------------------!
    ! Arguments:							      !
    !   grid: GRID_INFO type defining the grid to be allocated.               !
    !   cell: CELL_INFO type defining the cell to be used.                    !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 22/06/10.                                     !
    ! Based on code by Peter Haynes and Arash Mostofi (2003-2004).            !
    !=========================================================================!

    use comms, only: pub_total_num_nodes, pub_comms_group_size
    use geometry, only: OPERATOR(*), OPERATOR(.dot.), operator(.cross.)
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid
    type(CELL_INFO), intent(in) :: cell
    integer, intent(in) :: n1,n2,n3,ld1,ld2,ld3

    ! Local Variables
    integer :: ierr         ! Error flag

    ! Set sizes
    grid%n1 = n1
    grid%n2 = n2
    grid%n3 = n3
    grid%ld1 = ld1
    grid%ld2 = ld2
    grid%ld3 = ld3

    ! ndmh: initialise lattice vectors divided by number of grid points
    ! nmdh:  along that lattice vector
    grid%da1 = (1.0_DP/real(n1,kind=DP))*cell%a1
    grid%da2 = (1.0_DP/real(n2,kind=DP))*cell%a2
    grid%da3 = (1.0_DP/real(n3,kind=DP))*cell%a3

    ! ndmh: copy in reciprocal lattice vectors from the cell structure
    grid%b1 = cell%b1
    grid%b2 = cell%b2
    grid%b3 = cell%b3

    ! Set integration weight
    grid%weight = (grid%da1.cross.grid%da2).dot.grid%da3

    ! Zero integer describing entry for this grid in fourier_mod
    grid%fft_index = 0

    ! Allocate internal arrays in grid
    allocate(grid%node_slab12(n3),stat=ierr)
    call utils_alloc_check('cell_grid_init','grid%node_slab12',ierr)
    allocate(grid%node_slab23(n1),stat=ierr)
    call utils_alloc_check('cell_grid_init','grid%node_slab23',ierr)
    allocate(grid%first_slab12(0:pub_total_num_nodes),stat=ierr)
    call utils_alloc_check('cell_grid_init','grid%first_slab12',ierr)
    allocate(grid%last_slab12(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('cell_grid_init','grid%last_slab12',ierr)
    allocate(grid%group_displs_slabs12(0:pub_comms_group_size-1),stat=ierr)
    call utils_alloc_check('cell_grid_init','grid%group_displs_slabs12',ierr)
    allocate(grid%group_lengths_slabs12(0:pub_comms_group_size-1),stat=ierr)
    call utils_alloc_check('cell_grid_init','grid%group_lengths_slabs12',ierr)
    allocate(grid%first_slab23(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('cell_grid_init','grid%first_slab23',ierr)
    allocate(grid%last_slab23(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('cell_grid_init','grid%last_slab23',ierr)

  end subroutine cell_grid_init

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine cell_grid_distribute(grid,cell,scale,grid_name,allocate_grids, &
       is_fine)

    !=========================================================================!
    ! This subroutine distributes the data across the nodes for the parallel  !
    ! FFT routine, and sets up appropriate indices.                           !
    !-------------------------------------------------------------------------!
    ! Arguments:							      !
    !   grid : GRID_INFO structure holding information about the grid.        !
    !   cell : CELL_INFO structure holding information about the cell.        !
    !   scale : Factor by which this grid should be finer than standard grid  !
    !   grid_name : Identifying string for this grid (eg 'Fine', 'Standard')  !
    !   allocate_grids (input) : whether to fill the coulomb_recip array      !
    !   is_fine (input, optional) : .true. iff this is the fine grid          !
    !                               or dbl grid and they are commensurate)    !
    !-------------------------------------------------------------------------!
    ! Modules used:						              !
    !   utils           : for memory allocation/deallocation checks           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:					              !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:						      !
    !   grid%n1 > 0 etc.                                                      !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 8/6/04                                         !
    ! Based on an original version written by Peter Haynes, 17/7/03	      !
    ! No longer needs elements array as argument (Nicholas Hine 13/03/08)     !
    ! Name changed and moved to cell_grid mod by Nicholas Hine, 13/12/10.     !
    !=========================================================================!

    use comms, only: comms_abort, comms_reduce, comms_bcast, pub_on_root, &
         pub_my_node_id, pub_total_num_nodes, pub_first_node_in_group, &
         pub_group_comm, pub_comms_group_size, pub_my_rank_in_group
    use constants, only: DP, NORMAL
    use rundat, only: pub_output_detail, pub_multigrid_in_use, &
         pub_is_multigrid_nlevels, pub_is_discretization_order
    use simulation_cell, only: CELL_INFO
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_abort, utils_assert

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid    ! The grid to distribute
    type(CELL_INFO), intent(in)    :: cell    ! The cell for this grid
    logical, intent(in) :: allocate_grids     ! Whether we need to allocate
                                              ! coulomb_recip
    real(kind=DP), intent(in) :: scale
    character(len=*), intent(in) :: grid_name

    logical, intent(in), optional :: is_fine

    ! Local variables
    integer :: n1, n2, n3                     ! Local copy of fine grid sizes
    integer :: ld1, ld2, ld3                  ! Local copy of dimensions
    integer :: n1half, n2half, n3half         ! Num grid points in half-spaces
    integer :: ierr                           ! Error flag
    integer :: n12slabs, n23slabs             ! Checksums on slabs
    integer :: node                           ! Node counter
    integer :: slab                           ! Slab counter
    integer :: i1, i2, i3                     ! Fine grid counters
    integer :: i3start                        ! Trick to avoid origin
    integer :: j2, j3                         ! Auxiliary fine grid counters
    integer :: min_slabs                      ! ndmh: min slabs per node
    integer :: mean_slabs                     ! ndmh: mean slabs per node
    real(kind=DP) :: b1(3),b2(3),b3(3)        ! Local copies of recip lattice
    real(kind=DP) :: g1(3),g12(3),g(3)        ! Temporary recip-space vectors
    real(kind=DP) :: gsq                      ! Squared magnitude of G-vector

    ! Check arguments
    ! Grids are made even if scale>1 (ie not for standard grid, but for others)
    n1 = cell%total_pt1 * scale
    if (scale>1.0_DP) n1 = n1 + mod(n1,2)
    ld1 = n1 + 2
    n1half = n1/2+1
    if (n1 < 1) call utils_abort('Error in cell_grid_distribute: grid%n1 <= 0')
    n2 = cell%total_pt2 * scale
    if (scale>1.0_DP) n2 = n2 + mod(n2,2)
    ld2 = n2
    n2half = n2/2
    if (n2 < 1) call utils_abort('Error in cell_grid_distribute: grid%n2 <= 0')
    n3 = cell%total_pt3 * scale
    if (scale>1.0_DP) n3 = n3 + mod(n3,2)
    ld3 = n3
    n3half = n3/2
    if (n3 < 1) call utils_abort('Error in cell_grid_distribute: grid%n3 <= 0')

    ! Allocate module arrays
    call cell_grid_init(grid,cell,n1,n2,n3,ld1,ld2,ld3)

    !------------!
    ! '12' slabs !
    !------------!

    ! In real-space the '12' slabs are indexed by their position along the
    !   '3' axis: from 1 to grid%n3.

    ! As far as possible, share slabs out equally across all nodes
    grid%num_my_slabs12 = n3 / pub_total_num_nodes
    if (filling_order(pub_my_node_id,pub_total_num_nodes) < &
         mod(n3,pub_total_num_nodes)) &
         grid%num_my_slabs12 = grid%num_my_slabs12 + 1

    ! jd: If the parallel multigrid solver is in use, adjust the distribution
    !     to be compatible with the multigrid
    if(present(is_fine)) then
       if(pub_multigrid_in_use .and. is_fine) then
          call internal_adjust_for_multigrid(pub_is_multigrid_nlevels, &
               pub_is_discretization_order)
       end if
    end if

    ! Check that this allocation has worked
    n12slabs = grid%num_my_slabs12
    call comms_reduce('SUM',n12slabs)
    call utils_assert(n12slabs==n3, 'Error in cell_grid_distribute: ''12'' &
         &slab checksum failed.',n12slabs,n3)

    ! Calculate maximum number of slabs on any node
    grid%max_slabs12 = grid%num_my_slabs12
    call comms_reduce('MAX',grid%max_slabs12)

    ! Calculate number of slabs in this comms group, ensuring there are
    ! at least grid%max_slabs12 slabs beyond the start of the slabs of the
    ! last node in the group (so that arrays do not overrun).
    grid%num_group_slabs12 = grid%num_my_slabs12
    if ((pub_my_rank_in_group==pub_comms_group_size-1) .and. &
         (grid%num_my_slabs12<grid%max_slabs12)) then
       grid%num_group_slabs12 = grid%num_group_slabs12 + 1
    end if
    call comms_reduce('SUM',grid%num_group_slabs12,comm=pub_group_comm)

    ! Calculate maximum number of slabs in any comms group
    grid%max_group_slabs12 = grid%num_group_slabs12
    call comms_reduce('MAX',grid%max_group_slabs12)

    ! Allocate '12' slabs
    slab = 1

    ! Loop over nodes
    do node=0,pub_total_num_nodes-1

       ! Get number of '12' slabs on this node
       n12slabs = grid%num_my_slabs12
       call comms_bcast(node,n12slabs)

       ! Allocate '12' slabs to node
       grid%first_slab12(node) = slab
       grid%last_slab12(node) = slab + n12slabs - 1

       ! Mark these slabs as belonging to node
       grid%node_slab12(slab:slab+n12slabs-1) = node

       slab = slab + n12slabs

    end do
    grid%first_slab12(pub_total_num_nodes) = slab

    ! ndmh: set up lists of lengths and displacements of the 12-slabs of
    ! ndmh: all the nodes in this node's comms group
    do node=0,pub_comms_group_size-1
       grid%group_displs_slabs12(node) = &
            grid%first_slab12(pub_first_node_in_group+node) &
            - grid%first_slab12(pub_first_node_in_group)
       grid%group_lengths_slabs12(node) = &
            grid%first_slab12(pub_first_node_in_group+node+1) &
            - grid%first_slab12(pub_first_node_in_group+node)
    end do
    grid%my_first_slab12_in_group = &
         grid%group_displs_slabs12(pub_my_rank_in_group) + 1
    grid%my_last_slab12_in_group = &
         grid%group_displs_slabs12(pub_my_rank_in_group) &
         + grid%group_lengths_slabs12(pub_my_rank_in_group)
    grid%group_displs_slabs12 = grid%group_displs_slabs12*grid%ld1*grid%ld2
    grid%group_lengths_slabs12 = grid%group_lengths_slabs12*grid%ld1*grid%ld2

    !------------!
    ! '23' slabs !
    !------------!

    ! In reciprocal-space the '23' slabs are indexed by their position along the
    !   '1' axis: from 1 to cell%total_pt1_fine/2+1 (because of the
    !   inversion symmetry of the FFT of a real function)

    ! As far as possible, share slabs out equally across all nodes
    ! ndmh: using filling_order function:
    ! ndmh: eg       node:  0   1   2   3   4   5   6   7
    ! ndmh: filling order:  8   4   2   5   1   6   3   7

    grid%num_slabs23 = n1half / pub_total_num_nodes
    if (filling_order(pub_my_node_id,pub_total_num_nodes) < &
         mod(n1half,pub_total_num_nodes)) &
         grid%num_slabs23 = grid%num_slabs23 + 1

    ! Check that this allocation has worked
    n23slabs = grid%num_slabs23
    call comms_reduce('SUM',n23slabs)
    if (n23slabs /= n1half) then
       call utils_abort('Error in cell_grid_distribute: ''23'' &
            &slab checksum failed.')
    end if

    ! Calculate maximum number of slabs on any node
    grid%max_slabs23 = grid%num_slabs23
    call comms_reduce('MAX',grid%max_slabs23)

    ! Allocate '23' slabs
    slab = 1

    ! Loop over nodes
    do node=0,pub_total_num_nodes-1

       ! Get number of '23' slabs on this node
       n23slabs = grid%num_slabs23
       call comms_bcast(node,n23slabs)

       ! Allocate '23' slabs to node
       grid%first_slab23(node) = slab
       grid%last_slab23(node) = slab + n23slabs - 1

       ! Mark these slabs as belonging to node
       grid%node_slab23(slab:slab+n23slabs-1) = node

       slab = slab + n23slabs

    end do

    ! Warning if there are any nodes with no '23' slabs
    if (n1half / pub_total_num_nodes < 1) then
       if (pub_on_root.and.pub_output_detail>=VERBOSE) then
          write(stdout,*) ' '
          write(stdout,'(a)') 'WARNING in cell_grid_distribute:'
          write(stdout,'(a)') 'Average of less than one 23-slab per node, so &
               &workload cannot be'
          write(stdout,'(a)') 'equally balanced between nodes.'
       end if
    end if

    if (allocate_grids) then

       !-----------------!
       ! Cartesian grids !
       !-----------------!

       ! Reciprocal-space coefficients of the Coulomb interaction:
       !  grid%coulomb_recip(:,:,:)
       !   first index is Cartesian component 1-3, magnitude in 4, inverse
       !       squared magnitude in 5
       !   second index is position along '3' axis
       !   third index is position along '2' axis
       !   fourth index is '23' slab index on this node

       ! Allocate module arrays
       allocate(grid%coulomb_recip(ld3,ld2,grid%max_slabs23),stat=ierr)
       call utils_alloc_check('cell_grid_distribute', &
            'grid%coulomb_recip',ierr)

       ! Make local copies of reciprocal lattice vectors
       b1(1) = cell%b1%x ; b1(2) = cell%b1%y ; b1(3) = cell%b1%z
       b2(1) = cell%b2%x ; b2(2) = cell%b2%y ; b2(3) = cell%b2%z
       b3(1) = cell%b3%x ; b3(2) = cell%b3%y ; b3(3) = cell%b3%z

       ! Reciprocal-space grid
       grid%coulomb_recip = 0.0_DP
       i3start = 1
       if (grid%first_slab23(pub_my_node_id)==1) i3start = 2     ! Avoid G=0
       i1 = grid%first_slab23(pub_my_node_id)
       do slab=1,grid%num_slabs23                ! Loop over '23' slabs on this node
          g1 = (i1-1) * b1
          do i2=1,n2                            ! Loop along '2' in this slab
             if (i2 > n2half+1) then
                j2 = i2 - n2
             else
                j2 = i2
             end if
             g12 = g1 + (j2-1) * b2
             do i3=i3start,n3                   ! Loop along '3' in this slab
                if (i3 > n3half+1) then
                   j3 = i3 - n3
                else
                   j3 = i3
                end if
                g = g12 + (j3-1) * b3
                gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
                grid%coulomb_recip(i3,i2,slab) = 1.0_DP / gsq
             end do
             i3start = 1
          end do
          i1 = i1 + 1
       end do

    end if ! allocate grids

    ! All done...
    grid%distributed = .true.

    ! Temporary diagnostics
    if (pub_on_root.and.pub_output_detail>=NORMAL) then
       write(stdout,'(/a)') &
            '********************** Fourier parallelisation information &
            &*********************'
       write(stdout,'(7x,2a,3i6)') grid_name,' (whole simulation cell) &
            &dimensions: ', n1,n2,n3
    end if
    if (pub_output_detail>=VERBOSE) then
       if (pub_on_root) write(stdout,'(/15x,a)') &
            '  Real-space grid (''12'' slabs):'
       do node=0,pub_total_num_nodes-1
          i1 = grid%first_slab12(node)
          i2 = grid%last_slab12(node)
          n12slabs = grid%num_my_slabs12
          call comms_bcast(node,n12slabs)
          if (pub_on_root) then
             if (n12slabs>0) then
                write(stdout,'(15x,a,4(i4,a))') '    Node ',node, &
                     ': ',n12slabs,' (from ',i1,' to ',i2,')'
             else
                write(stdout,'(15x,a,2(i4,a))') '    Node ',node, &
                     ': ',n12slabs
             end if
          end if
       end do
       if (pub_on_root) write(stdout,'(/15x,a)') &
            '  Reciprocal-space grid (''23'' slabs):'
       do node=0,pub_total_num_nodes-1
          i1 = grid%first_slab23(node)
          i2 = grid%last_slab23(node)
          n23slabs = grid%num_slabs23
          call comms_bcast(node,n23slabs)
          if (pub_on_root) then
             if (n23slabs>0) then
                write(stdout,'(15x,a,4(i4,a))') '    Node ',node, &
                     ': ',n23slabs,' (from ',i1,' to ',i2,')'
             else
                write(stdout,'(15x,a,2(i4,a))') '    Node ',node, &
                     ': ',n23slabs
             end if
          end if
       end do
    else if (pub_output_detail==NORMAL) then
       min_slabs = grid%num_my_slabs12
       call comms_reduce('MIN',min_slabs)
       mean_slabs = int(real(n3,kind=DP) / real(pub_total_num_nodes,kind=DP))
       if (pub_on_root) then
          write(stdout,'(8x,a,3i6)') &
               '  Real-space (''12'') slabs/node (min max mean): ', &
               min_slabs, grid%max_slabs12, mean_slabs
       end if
       min_slabs = grid%num_slabs23
       call comms_reduce('MIN',min_slabs)
       mean_slabs = int(real(n1half,kind=DP) / real(pub_total_num_nodes,kind=DP))
       if (pub_on_root) then
          write(stdout,'(7x,a,3i6)') &
               '  Recip-space (''23'') slabs/node (min max mean): ', &
               min_slabs, grid%max_slabs23, mean_slabs
       end if
    end if
    if (pub_on_root.and.pub_output_detail>=NORMAL) write(stdout,'(a)') &
         '******************************************************************&
          &**************'

contains

    !--------------------------------------------------------------------------

      ! ndmh: returns the order in which node m should be allocated a slab,
      ! ndmh: given N total nodes. Attempts to distribute load efficiently.
      integer function filling_order(m,N)

        ! Arguments
        integer, intent(in) :: m
        integer, intent(in) :: N

        ! Local Variables
        integer :: p, q, s, t

        ! Special case - root node is last to fill
        if (m==0) then
           filling_order = N - 1
           return
        end if

        ! Find highest power of two less than or equal to N
        p = int(log(real(N,DP))/log(2.0_DP))
        if (2**p==N) then
           p = p - 1
        end if

        ! Now count up in filling order until m is found, jumping by
        ! progressively lower powers of two through the order
        filling_order = 0
        do q = p,0,-1
           t = 2**q
           do s=1,N
              if (t>N-1) exit
              ! Found the right node, so return current filling_order
              if (t==m) return
              filling_order = filling_order + 1
              t = t + 2**(q+1)
           end do
        end do

      end function filling_order

    !--------------------------------------------------------------------------

      ! jd: Adjusts the slab distribution for use with multigrid
      subroutine internal_adjust_for_multigrid(nlevels,fd_order)

        ! Arguments
        integer, intent(in) :: nlevels
        integer, intent(in) :: fd_order

        ! Local variables
        integer :: nslabs, nslabs_pq, nslabs_now, nslabs_left
        integer :: multigrid_granularity
        integer :: min_slabs, nfull_pairs, nnodes_to_burden
        integer :: deficit_on_last
        logical :: i_am_last
        integer :: local_fd_order

        ! ---------------------------------------------------

        i_am_last = (pub_my_node_id == pub_total_num_nodes-1)

        ! Determine the multigrid granularity
        multigrid_granularity = 2**(nlevels+1)

        ! Determine nslabs
        nslabs = grid%num_my_slabs12
        call comms_reduce('SUM',nslabs)

        ! Round down to nearest multiple of granularity
        nslabs_pq = ((nslabs-1)/multigrid_granularity)*multigrid_granularity+1

        ! Mind that fd_order == 2 translates into 12th order FDs for 2nd derivs.
        if(fd_order == 2) then
           local_fd_order = 12
        else
           local_fd_order = fd_order
        end if

        ! Determine the smallest number of slabs possible for this nproc
        ! Each proc must have at least (local_fd_order/2)+1, rounded up to the
        ! nearest even.
        min_slabs = ((local_fd_order/2+1)+mod(local_fd_order/2+1,2)) * &
             pub_total_num_nodes
        call utils_assert(nslabs_pq >= min_slabs, &
             'Too few slabs per processor to run a multigrid calculation at &
             &this is_discretization_order. Decrease the number of processors.',min_slabs)

        ! If running on only processor, no point in doing anything else
        if(pub_total_num_nodes == 1) return

        ! If running on more processors, distribute the slabs

        ! Start by giving every proc the smallest number of slabs possible
        ! This is an even number for every proc
        grid%num_my_slabs12 = min_slabs/pub_total_num_nodes

        ! See how many more need to be distributed
        nslabs_now = grid%num_my_slabs12
        call comms_reduce('SUM',nslabs_now)
        nslabs_left = nslabs - nslabs_now

        ! If this is odd, only the last one can get the tail
        if(mod(nslabs_left,2) /= 0) then
           nslabs_left = nslabs_left - 1
           if(i_am_last) grid%num_my_slabs12 = grid%num_my_slabs12 + 1
        end if

        ! See how many pairs of slabs can be given to all nodes
        nfull_pairs = nslabs_left/(2*pub_total_num_nodes)

        ! Give everyone these pairs
        grid%num_my_slabs12 = grid%num_my_slabs12 + nfull_pairs * 2

        ! See how many more need to be distributed
        nslabs_left = nslabs_left - nfull_pairs*2*pub_total_num_nodes

        ! Should definitely be fewer than 2 for everyone
        ! and it should be an even number
        call utils_assert(nslabs_left < 2*pub_total_num_nodes .and. &
             mod(nslabs_left,2)==0, &
             'Logic error (1) in internal_adjust_for_multigrid',nslabs_left)

        ! Add them, in pairs, to the first processors (they are likely to have
        ! an easier job with slowly varying potentials)
        nnodes_to_burden = nslabs_left/2

        if(pub_my_node_id < nnodes_to_burden) &
           grid%num_my_slabs12 = grid%num_my_slabs12 + 2

        ! Last processor might need more slabs to account for margin between
        ! nslabs and nslabs_pq
        deficit_on_last = 0
        if(i_am_last) then
           deficit_on_last = min_slabs/pub_total_num_nodes ! must have this many
           deficit_on_last = deficit_on_last-grid%num_my_slabs12 ! have already
           deficit_on_last = deficit_on_last+(nslabs-nslabs_pq)! don't count
        end if
        call comms_bcast(pub_total_num_nodes-1,deficit_on_last)

        if(deficit_on_last <= 0) return

        ! If the deficit is odd, we cannot easily take it from other nodes
        ! (because we must work in pairs). In that case, *increase* the deficit
        ! to an even number. That often happens if the deficit is 1 -- there
        ! are no processors from which we can steal 1 slab -- we must steal 2.
        if(mod(deficit_on_last,2)/=0) then
           deficit_on_last = deficit_on_last+1
        end if

        ! Take pairs of slabs from first processors until the deficit is
        ! reduced to zero
        nnodes_to_burden = deficit_on_last/2
        if(pub_my_node_id < nnodes_to_burden) &
           grid%num_my_slabs12 = grid%num_my_slabs12 - 2

        ! And give them to the last one
        if(i_am_last) grid%num_my_slabs12 = grid%num_my_slabs12+deficit_on_last


      end subroutine internal_adjust_for_multigrid

    !--------------------------------------------------------------------------

  end subroutine cell_grid_distribute


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  subroutine cell_grid_deposit_box(cell_fine, &                        ! output
       box_fine, buffer_fine, grid, nf1, nf2, nf3, ldf1, ldf2,  &      ! input
       cell_start1, cell_start2, cell_start3, i_have_box, &            ! input
       group_comms)

    !========================================================================!
    ! This (parallel) subroutine deposits into the 12-slabs that belong to   !
    ! pub_my_node_id the parts of the current fine grid box from             !
    ! pub_my_node_id (if it has the box on it) or the data from the          !
    ! box on any other processor that contribute to the 12-slabs             !
    ! of pub_my_node_id.                                                     !
    ! At the same time it sends to all other processors the parts of         !
    ! the fine grid box (if is on pub_my_node_id) which contribute to their  !
    ! 12-slabs.                                                              !
    !------------------------------------------------------------------------!
    ! Key array: buffer_fine is a buffer which has the dimensions of         !
    !   the fine grid box in directions 1 and 2 and the dimensions           !
    !   of pub_max_slabs12 in direction 3. All incoming communications       !
    !   of slab contributions from other nodes are stored into buffer_fine   !
    !   before they are deposited to the charge 12-slabs belonging to        !
    !   pub_my_node_id.                                                      !
    !------------------------------------------------------------------------!
    ! This subroutine was originally written by Chris-Kriton Skylaris        !
    ! in spring 2001 for the ONES program and was originally called          !
    ! "density_put_pair_box_to_cell". It was renamed to                      !
    ! "density_put_super_box_to_cell" by Arash A. Mostofi in September 2002. !
    ! It was rewritten by Chris-Kriton Skylaris on 15/6/2004 so that         !
    ! it works with the data-parallel charge density and is part of          !
    ! the ONETEP program.                                                    !
    ! Spin polarised by Peter Haynes, July 2006                              !
    ! It was moved to basis_mod and renamed to basis_deposit_box_to_cell     !
    ! by Nicholas Hine in February 2008                                      !
    ! Removed references to tightbox and simplified argument list, Nicholas  !
    ! Hine, October 2009.                                                    !
    ! Modified to use comms_alltoall by Nicholas Hine, for improved parallel !
    ! efficiency, and tidied up, November 2009.                              !
    ! Modified for GRID_INFO type by Nicholas Hine, June 2010.               !
    ! Major changes for group comms by Nicholas Hine, December 2010.         !
    ! Moved to new cell_grid_mod and renamed by Nicholas Hine, January 2011. !
    ! Removed spin-polarisation of arrays - now handled as a loop in the     !
    ! calling routine - for simplicity and portability, by Nicholas Hine     !
    ! in February 2011.                                                      !
    !========================================================================!

    use comms, only: comms_recv, comms_send, comms_free, pub_total_num_nodes, &
         pub_my_node_id, pub_my_rank_in_group, pub_comms_group_size
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    ! sizes of cell grid
    type(GRID_INFO),intent(in)   :: grid
    ! sizes of box_fine:
    integer, intent(in)          :: nf1, nf2, nf3, ldf1, ldf2
    ! box start wrt simulation cell:
    integer, intent(in)          :: cell_start1, cell_start2, cell_start3
    ! data in fine-grid box to deposit
    real(kind=DP), intent(in)    :: box_fine(ldf1, ldf2, nf3)
    real(kind=DP), intent(inout) :: buffer_fine(ldf1, ldf2, &
         grid%max_group_slabs12)
    real(kind=DP), intent(inout) :: cell_fine(grid%ld1, grid%ld2, &
         grid%max_group_slabs12)
    ! whether pub_my_node_id has a box to deposit or not:
    logical, intent(in)          :: i_have_box
    ! whether doing group-based comms or not
    logical, intent(in)          :: group_comms

    ! Local Variables
#ifdef ACCELRYS
    integer :: i1,i2           ! point counters in 1,2 directions
#endif
    logical :: middle1, middle2, middle3  ! whether the box extends out of the sim cell in any dimension
    integer :: start1_fine, start2_fine, start3_fine ! start of first block of box wrt sim cell
    integer :: end1_fine, end2_fine, end3_fine       ! end of second block of box wrt sim cell
    integer :: mid_start1, mid_start2, mid_start3    ! start of second block of box wrt sim cell
    integer :: mid_end1, mid_end2, mid_end3          ! end of first block of box wrt sim cell
    integer :: start3_node          ! node which poseseses first slab of first block along direction 3
    integer :: mid_end3_node        ! node which poseseses last slab of first block along direction 3
    integer :: mid_start3_node      ! node which poseseses first slab of second block along direction 3
    integer :: end3_node            ! node which poseseses last slab of second block along direction 3
    integer :: ierr              ! error flag
    integer :: box_first3        ! start of first block of box for current node in direction 3
    integer :: box_second3       ! start of second block of box for current node in direction 3
    integer :: send_node       ! node counter
    integer :: recv_node       ! node counter
    integer :: start3_sent     ! start of first block in dir-3 in resulting buffer_fine just sent
    integer :: mid_end3_sent   ! end of first block in dir-3 in resulting buffer_fine just sent
    integer :: mid_start3_sent ! start of second block in dir-3 in resulting buffer_fine just sent
    integer :: end3_sent       ! end of second block in dir-3 in resulting buffer_fine just sent
    integer :: box_stop3 ! end of first or second block of box for current node in direction 3
    integer :: block_size         ! size counter in dir-3
    integer :: first_block_size   ! size counter in dir-3
    integer :: second_block_size  ! size counter in dir-3
    integer :: first_bf1    ! size of first block of box in direction 1
    integer :: first_bf2    ! size of first block of box in direction 2
    integer :: second_bf1   ! start of second block of box in direction 1
    integer :: second_bf2   ! start of second block of box in direction 2
    integer :: my_first3_start  ! start of first block of box of pub_my_node_id for its 12-slabs
    integer :: my_first3_end    ! end of first block of box of pub_my_node_id for its 12-slabs
    integer :: my_second3_start ! start of second block of box of pub_my_node_id for its 12-slabs
    integer :: my_second3_end   ! end of second block of box of pub_my_node_id for its 12-slabs
    integer :: my_box_second3 ! start of second block of (whole of) box of pub_my_node_id
    integer :: superstep    ! parallel superstep counter
    integer :: node_stride
    integer, allocatable, dimension(:,:) :: others_limits  ! limits of the eight possible blocks from node x
    integer, allocatable, dimension(:,:) :: my_limits ! all limits for the buffer_fines of all nodes
    ! that will result from the box of pub_my_node_id

    call timer_clock('cell_grid_deposit_box', 1)
#ifdef ITC_TRACE
    call VTBEGIN(vt_basis_deposit_box_to_cell, vt_err)
#endif

    ! ndmh: check for nonsensical input
    if (((cell_start1<=-grid%n1).or.(cell_start1>grid%n1).or. &
         (cell_start2<=-grid%n2).or.(cell_start2>grid%n2).or. &
         (cell_start3<=-grid%n3).or.(cell_start3>grid%n3)).and.i_have_box) then
       write(stdout,'(a,i6,a)') 'On node',pub_my_node_id,':'
       write(stdout,'(3(a,i5))') 'cell_start1: ', cell_start1, &
            ', grid%n1: ',grid%n1, ', nf1: ', nf1
       write(stdout,'(3(a,i5))') 'cell_start2: ', cell_start2, &
            ', grid%n2: ',grid%n2, ', nf2: ', nf2
       write(stdout,'(3(a,i5))') 'cell_start3: ', cell_start3, &
            ', grid%n3: ',grid%n3, ', nf3: ', nf3
       call utils_abort('Error in cell_grid_deposit_box: invalid cell start &
            &position for deposited box.')
    end if

    allocate(my_limits(12,0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('cell_grid_deposit_box','my_limits',ierr)
    allocate(others_limits(12,0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('cell_grid_deposit_box','others_limits',ierr)

    call cell_grid_box_limits(my_limits,others_limits, &
         start3_node,mid_end3_node,mid_start3_node,end3_node,my_box_second3, &
         i_have_box,grid,cell_start1,cell_start2,cell_start3,nf1,nf2,nf3, &
         group_comms)

     ! qoh: Initialisations to prevent compiler warnings
    my_first3_start = -1; my_second3_start = -1
    my_first3_end = -2; my_second3_end = -2

    node_stride = 1
    if (group_comms) node_stride = pub_comms_group_size
    do superstep=0,pub_total_num_nodes-1,node_stride

       ! cks: --- FIND NODES FOR CURRENT SUPERSTEP ----------------------
       ! pdh: simplify
       send_node = modulo(pub_my_node_id+superstep,pub_total_num_nodes)
       recv_node = modulo(pub_my_node_id-superstep,pub_total_num_nodes)
       ! cks: --- END FIND NODES FOR CURRENT SUPERSTEP ------------------


       ! cks: >>>>>>>> MPI SEND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       ! cks: Unique tags for possible future use in the sends here.
       ! cks: tag_base =(pub_total_num_nodes * send_node + pub_my_node_id -1)*3 +1
       ! cks: first_tag =tag_base, second_tag =tag_base +1, third_tag= tag_base+2

       start3_sent   = my_limits(9, send_node)
       mid_end3_sent = my_limits(10,send_node)

       ! cks: send first block
       block_size = mid_end3_sent - start3_sent + 1
       if ( block_size .gt. 0) then

          if (send_node .eq. start3_node) then
             box_first3 = 1
          else
             box_first3 = grid%first_slab12(send_node - &
                  modulo(send_node,node_stride)) - &
                  grid%first_slab12(start3_node - &
                  modulo(start3_node,node_stride)) - &
                  my_limits(9, start3_node) + 2
          endif
          box_stop3 = box_first3 + block_size - 1

          if (pub_my_node_id .ne. send_node) then
             call comms_send(send_node, &
                  box_fine(:,:,box_first3:box_stop3))
           else
             ! cks: store where current block starts and ends in direction-3 of box
             ! cks: of pub_my_node_id
             my_first3_start = box_first3
             my_first3_end   = box_stop3
          endif

       endif

       mid_start3_sent = my_limits(11, send_node)
       end3_sent       = my_limits(12, send_node)

       ! cks: send second block
       block_size = end3_sent - mid_start3_sent + 1
       if ( block_size .gt. 0) then

          box_second3 = my_box_second3
          if (send_node .ne. mid_start3_node) then
             box_second3 = box_second3 + &
                  grid%first_slab12(send_node - &
                  modulo(send_node,node_stride)) - &
                  grid%first_slab12(mid_start3_node - &
                  modulo(mid_start3_node,node_stride))
          endif
          box_stop3 = box_second3 + block_size - 1


          if (pub_my_node_id .ne. send_node) then
             call comms_send(send_node, &
                  box_fine(:,:,box_second3:box_stop3))
          else
             ! cks: store where current block starts and ends in direction-3 of box
             ! cks: of pub_my_node_id
             my_second3_start = box_second3
             my_second3_end   = box_stop3

          endif


       endif

       ! cks: >>>>>>>> END MPI SEND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


       ! cks: <<<<<<< MPI RECEIVE - AND IMMEDIATELY DEPOSIT <<<<<<<<<<<<<<<<<<<<

       ! cks: Unique tags for possible future use in the recvs here.
       ! cks: tag_base =(pub_total_num_nodes * pub_my_node_id + recv_node -1)*3 +1
       ! cks: first_tag =tag_base, second_tag =tag_base +1, third_tag= tag_base+2

       ! cks: deposit only stuff from other nodes or my stuff if I have any
       if ( (pub_my_node_id .ne. recv_node) .or. i_have_box) then

          start1_fine     = others_limits(1,recv_node)
          mid_end1        = others_limits(2,recv_node)
          mid_start1      = others_limits(3,recv_node)
          end1_fine       = others_limits(4,recv_node)

          middle1 =.false.
          if ( end1_fine >= mid_start1 ) middle1 = .true.

          start2_fine     = others_limits(5,recv_node)
          mid_end2        = others_limits(6,recv_node)
          mid_start2      = others_limits(7,recv_node)
          end2_fine       = others_limits(8,recv_node)

          middle2 =.false.
          if ( end2_fine >= mid_start2 ) middle2 = .true.

          start3_fine     = others_limits(9,recv_node)
          mid_end3        = others_limits(10,recv_node)
          mid_start3      = others_limits(11,recv_node)
          end3_fine       = others_limits(12,recv_node)

          middle3 =.false.
          if (end3_fine >= mid_start3) middle3 =.true.

          ! cks: receive FIRST block
          first_block_size = mid_end3 - start3_fine + 1
          if ( first_block_size .gt. 0) then

             if (pub_my_node_id .ne. recv_node) then
                call comms_recv(recv_node, &
                     buffer_fine(:,:,start3_fine:mid_end3))
             else
                ! cks: put block from box_fine to my buffer_fine
                buffer_fine(:,:,start3_fine:mid_end3) = &
                     box_fine(:,:,my_first3_start:my_first3_end)
             end if

          end if

          ! cks: receive SECOND block
          second_block_size = end3_fine - mid_start3 + 1
          if (second_block_size .gt. 0) then

             if (pub_my_node_id .ne. recv_node) then
                call comms_recv(recv_node, &
                     buffer_fine(:,:,mid_start3:end3_fine))
             else
                ! cks: put block from box_fine to my buffer_fine
                buffer_fine(:,:,mid_start3:end3_fine) = &
                     box_fine(:,:,my_second3_start:my_second3_end)
             end if

          end if

          first_bf1 = mid_end1 - start1_fine + 1
          first_bf2 = mid_end2 - start2_fine + 1

          second_bf1 = mid_end1 - start1_fine + 2
          second_bf2 = mid_end2 - start2_fine + 2


          ! cks: depositing to 12-slabs of pub_my_node_id starts here
          ! cks: There are two blocks per 12-slab-node-block possible per every direction
          ! cks: so there are eight deposition steps possible in total.
          cell_fine(start1_fine:mid_end1,start2_fine:mid_end2, &
               start3_fine:mid_end3) = &
               cell_fine(start1_fine:mid_end1,start2_fine:mid_end2, &
               start3_fine:mid_end3) + &
               buffer_fine(1:first_bf1,1:first_bf2,start3_fine:mid_end3)

          if (middle1) then
             cell_fine(mid_start1:end1_fine,start2_fine:mid_end2, &
                  start3_fine:mid_end3) = &
                  cell_fine(mid_start1:end1_fine,start2_fine:mid_end2, &
                  start3_fine:mid_end3) + &
                  buffer_fine(second_bf1:nf1,1:first_bf2, &
                  start3_fine:mid_end3)
          end if

          if (middle2) then
             cell_fine(start1_fine:mid_end1,mid_start2:end2_fine, &
                  start3_fine:mid_end3) = &
                  cell_fine(start1_fine:mid_end1,mid_start2:end2_fine, &
                  start3_fine:mid_end3) + &
                  buffer_fine(1:first_bf1,second_bf2:nf2, &
                  start3_fine:mid_end3)
          end if

          if (middle3) then
             cell_fine(start1_fine:mid_end1,start2_fine:mid_end2, &
                  mid_start3:end3_fine) = &
                  cell_fine(start1_fine:mid_end1,start2_fine:mid_end2, &
                  mid_start3:end3_fine) + &
                  buffer_fine(1:first_bf1,1:first_bf2, &
                  mid_start3:end3_fine)
          end if

          if (middle1 .and. middle2) then
             cell_fine(mid_start1:end1_fine,mid_start2:end2_fine, &
                  start3_fine:mid_end3) = &
                  cell_fine(mid_start1:end1_fine,mid_start2:end2_fine, &
                  start3_fine:mid_end3) + &
                  buffer_fine(second_bf1:nf1,second_bf2:nf2, &
                  start3_fine:mid_end3)
          end if

          if (middle1 .and. middle3) then
             cell_fine(mid_start1:end1_fine,start2_fine:mid_end2, &
                  mid_start3:end3_fine) = &
                  cell_fine(mid_start1:end1_fine,start2_fine:mid_end2, &
                  mid_start3:end3_fine) + &
                  buffer_fine(second_bf1:nf1,1:first_bf2, &
                  mid_start3:end3_fine)
          end if

          if (middle2 .and. middle3) then
             cell_fine(start1_fine:mid_end1,mid_start2:end2_fine, &
                  mid_start3:end3_fine) = &
                  cell_fine(start1_fine:mid_end1,mid_start2:end2_fine, &
                  mid_start3:end3_fine) + &
                  buffer_fine(1:first_bf1,second_bf2:nf2, &
                  mid_start3:end3_fine)
          end if

          if (middle1 .and. middle2 .and. middle3) then
#ifdef ACCELRYS
             do i2 = mid_start2,end2_fine
                do i1 = mid_start1,end1_fine
                   cell_fine(i1,i2,mid_start3:end3_fine) = &
                        cell_fine(i1,i2,mid_start3:end3_fine) + &
                        buffer_fine(second_bf1+i1-1,second_bf2+i2-1, &
                        mid_start3:end3_fine)
                end do
             end do
#else
             cell_fine(mid_start1:end1_fine,mid_start2:end2_fine, &
                  mid_start3:end3_fine) = &
                  cell_fine(mid_start1:end1_fine,mid_start2:end2_fine, &
                  mid_start3:end3_fine) + &
                  buffer_fine(second_bf1:nf1,second_bf2:nf2, &
                  mid_start3:end3_fine)
#endif
          end if


       end if
       ! cks: <<<<<<< END MPI RECEIVE - AND IMMEDIATELY DEPOSIT <<<<<<<<<<<<<<<

    end do

    ! ndmh: wait for sends to finish before leaving
    call comms_free

    deallocate(others_limits,stat=ierr)
    call utils_dealloc_check('cell_grid_deposit_box','others_limits',ierr)
    deallocate(my_limits,stat=ierr)
    call utils_dealloc_check('cell_grid_deposit_box','my_limits',ierr)

#ifdef ITC_TRACE
    call VTEND(vt_basis_deposit_box_to_cell, vt_err)
#endif
    call timer_clock('cell_grid_deposit_box', 2)

  end subroutine cell_grid_deposit_box

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine cell_grid_extract_box(box_fine, &                        ! output
       buffer_fine, cell_fine, grid, nf1, nf2, nf3, ldf1, ldf2, &     ! input
       cell_start1, cell_start2, cell_start3, i_need_box, group_comms)! input

    !======================================================================!
    ! This (parallel) subroutine deposits into the 12-slabs that belong to !
    ! pub_my_node_id the parts of the current fine grid box from           !
    ! pub_my_node_id (if the box is to be on it) and the current           !
    ! boxes of all other processors that contribute to the 12-slabs        !
    ! of pub_my_node_id.                                                   !
    ! At the same time it sends to all other processors the parts of       !
    ! the fine grid box (if it is on pub_my_node_id) which                 !
    ! contribute to their 12-slabs.                                        !
    !----------------------------------------------------------------------!
    ! Key array: buffer_fine is a buffer which has the dimensions of       !
    !   the fine grid box in directions 1 and 2 and the dimensions         !
    !   of pub_max_slabs12 in direction 3.                                 !
    !   All slab contributions to cell data from the box from a slab       !
    !   of one node are packed together into buffer_fine before being      !
    !   sent to the node where the box belongs.                            !
    !----------------------------------------------------------------------!
    ! This subroutine was originally written by Chris-Kriton Skylaris      !
    ! in spring 2001 for the ONES program and was originally called        !
    ! "potential_from_cell_to_pair_box". It was modified to work with      !
    ! the triple box and renamed to "basis_extract_fftbox_from_cell"       !
    ! by Arash A. Mostofi in September 2002.                               !
    ! It was rewritten by Chris-Kriton Skylaris on 10/7/2004 so that       !
    ! it works with the data-parallel local potential on the fine grid     !
    ! and is part of the ONETEP program.                                   !
    ! It was moved to basis_mod and renamed to basis_extract_box_from_cell !
    ! by Nicholas Hine in February 2008                                    !
    ! Removed references to tightbox and simplified argument list,         !
    ! Nicholas Hine, October 2009.                                         !
    ! Modified to use GRID_INFO type by Nicholas Hine in June 2010.        !
    ! Major changes for group comms by Nicholas Hine, December 2010.       !
    ! Moved to new cell_grid_mod and renamed by Nicholas Hine, Jan 2011    !
    !======================================================================!

    use comms, only: comms_send, comms_recv, &
         comms_waitany, comms_test, &
         pub_comms_group_size, pub_total_num_nodes, &
         pub_my_node_id, pub_null_handle
    use constants, only: DP
    use utils, only: utils_abort, utils_alloc_check,utils_dealloc_check
    use timer, only: timer_clock
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    integer, intent(in)          :: nf1  ! size of box_fine in dir-1
    integer, intent(in)          :: nf2  ! size of box_fine in dir-2
    integer, intent(in)          :: nf3  ! size of box_fine in dir-3
    integer, intent(in)          :: ldf1 ! full size of box_fine in dir-1
    integer, intent(in)          :: ldf2 ! full size of box_fine in dir-2
    integer, intent(in)          :: cell_start1
    integer, intent(in)          :: cell_start2
    integer, intent(in)          :: cell_start3
    real(kind=DP), intent(inout) :: buffer_fine(ldf1, ldf2, grid%max_group_slabs12)
    real(kind=DP), intent(in)    :: cell_fine(grid%ld1, grid%ld2, &
         grid%max_group_slabs12)
    real(kind=DP), intent(inout) :: box_fine(ldf1, ldf2, nf3)
    logical, intent(in)          :: i_need_box
    ! whether doing group-based comms or single-node comms
    logical, intent(in)          :: group_comms

    ! Internal declarations:
    logical :: middle1, middle2, middle3 ! whether the box extends out of the sim cell in any dimension
    integer :: start1_fine, start2_fine, start3_fine ! start of first block of box wrt sim cell
    integer :: end1_fine, end2_fine, end3_fine       ! end of second block of box wrt sim cell
    integer :: mid_start1, mid_start2, mid_start3    ! start of second block of box wrt sim cell
    integer :: mid_end1, mid_end2, mid_end3          ! end of first block of box wrt sim cell
    integer :: start3_node     ! node which poseseses first slab of first block along direction 3
    integer :: mid_end3_node   ! node which poseseses last slab of first block along direction 3
    integer :: mid_start3_node ! node which poseseses first slab of second block along direction 3
    integer :: end3_node       ! node which poseseses last slab of second block along direction 3
    integer :: send_node ! node counter
    integer :: recv_node ! node counter
    integer :: ierr      ! error flag
    integer :: superstep ! parallel superstep counter
    integer :: first_bf1  ! size of first block of box in direction 1
    integer :: first_bf2  ! size of first block of box in direction 2
    integer :: second_bf1 ! start of second block of box in direction 1
    integer :: second_bf2 ! start of second block of box in direction 2
    integer :: block_size ! size counter in dir-3
    integer :: my_first3_start  ! start of first block of box of pub_my_node_id for its 12-slabs
    integer :: my_first3_end    ! end of first block of box of pub_my_node_id for its 12-slabs
    integer :: my_second3_start ! start of second block of box of pub_my_node_id for its 12-slabs
    integer :: my_second3_end   ! end of second block of box of pub_my_node_id for its 12-slabs
    integer :: my_box_second3 ! start of second block of (whole of) box of pub_my_node_id

    integer, allocatable, dimension(:,:) :: my_limits ! all limits of my box on all nodes (slabs)
    integer, allocatable, dimension(:,:) :: others_limits ! limits of other nodes' boxes on my slab

    ! ndmh: for new comms system
    integer :: send_handles(2)
    logical :: first_send_finished
    logical :: second_send_finished
    integer :: node_stride

    call timer_clock('cell_grid_extract_box', 1)
#ifdef ITC_TRACE
    call VTBEGIN(vt_basis_extract_box_from_cell, vt_err)
#endif

    ! ndmh: check for nonsensical input
    if (((cell_start1<=-grid%n1).or.(cell_start1>grid%n1).or. &
         (cell_start2<=-grid%n2).or.(cell_start2>grid%n2).or. &
         (cell_start3<=-grid%n3).or.(cell_start3>grid%n3)).and.i_need_box) then
       write(stdout,'(a,i6,a)') 'On node',pub_my_node_id,':'
       write(stdout,'(3(a,i5))') 'cell_start1: ', cell_start1, &
            ', grid%n1: ',grid%n1, ', nf1: ', nf1
       write(stdout,'(3(a,i5))') 'cell_start2: ', cell_start2, &
            ', grid%n2: ',grid%n2, ', nf2: ', nf2
       write(stdout,'(3(a,i5))') 'cell_start3: ', cell_start3, &
            ', grid%n3: ',grid%n3, ', nf3: ', nf3
       call utils_abort('Error in cell_grid_deposit_box: invalid cell start &
            &position for deposited box.')
    end if

     allocate(my_limits(12,0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('cell_grid_extract_box','my_limits',ierr)
    allocate(others_limits(12,0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('cell_grid_extract_box','others_limits',ierr)

    call cell_grid_box_limits(my_limits,others_limits, &
         start3_node,mid_end3_node,mid_start3_node,end3_node,my_box_second3, &
         i_need_box,grid,cell_start1,cell_start2,cell_start3,nf1,nf2,nf3, &
         group_comms)

    ! qoh: Initialisations to prevent compiler warnings
    my_first3_start = -1; my_second3_start = -1
    my_first3_end = -2; my_second3_end = -2

#ifdef DEBUG
    ! ndmh: initialise any padding around the box to zero to prevent warnings
    box_fine(nf1+1:ldf1,:,:) = 0.0_DP
    box_fine(:,nf2+1:ldf2,:) = 0.0_DP
#endif

    !---------------------------------------------------------------------------
    ! ##### cks: Fill up boxes with slab-12 potential
    !---------------------------------------------------------------------------
    send_handles = 0
    node_stride = 1
    if (group_comms) node_stride = pub_comms_group_size
    do superstep=0,pub_total_num_nodes-1,node_stride

       ! cks: --- FIND NODES FOR CURRENT SUPERSTEP ----------------------
       send_node = modulo(pub_my_node_id+superstep,pub_total_num_nodes)
       recv_node = modulo(pub_my_node_id-superstep,pub_total_num_nodes)
       ! cks: --- END FIND NODES FOR CURRENT SUPERSTEP ------------------



       ! cks: ---- limits for filling current buffer_fine before send--------
       start1_fine     = others_limits(1, send_node)
       mid_end1        = others_limits(2, send_node)
       mid_start1      = others_limits(3, send_node)
       end1_fine       = others_limits(4, send_node)
       middle1 = (end1_fine >= mid_start1)

       start2_fine     = others_limits(5, send_node)
       mid_end2        = others_limits(6, send_node)
       mid_start2      = others_limits(7, send_node)
       end2_fine       = others_limits(8, send_node)
       middle2 = (end2_fine >= mid_start2)

       start3_fine     = others_limits(9, send_node)
       mid_end3        = others_limits(10, send_node)
       mid_start3      = others_limits(11, send_node)
       end3_fine       = others_limits(12, send_node)
       middle3 = (end3_fine >= mid_start3)


       first_bf1  = mid_end1-start1_fine+1
       first_bf2  = mid_end2-start2_fine+1

       second_bf1 = mid_end1-start1_fine+2
       second_bf2 = mid_end2-start2_fine+2
       ! cks: ---- END limits for filling current buffer_fine before send------


       ! cks: ------- fill up FIRST AND SECOND parts of buffer-----------------

       buffer_fine(1:first_bf1,1:first_bf2,start3_fine:mid_end3) = &
            cell_fine(start1_fine:mid_end1,start2_fine:mid_end2, &
            start3_fine:mid_end3)

       if (middle1) &
            buffer_fine(second_bf1:nf1,1:first_bf2,start3_fine:mid_end3) = &
            cell_fine(mid_start1:end1_fine,start2_fine:mid_end2, &
            start3_fine:mid_end3)

       if (middle2) &
            buffer_fine(1:first_bf1,second_bf2:nf2,start3_fine:mid_end3) = &
            cell_fine(start1_fine:mid_end1,mid_start2:end2_fine, &
            start3_fine:mid_end3)

       if (middle1 .and. middle2) &
            buffer_fine(second_bf1:nf1,second_bf2:nf2,start3_fine:mid_end3) = &
            cell_fine(mid_start1:end1_fine,mid_start2:end2_fine, &
            start3_fine:mid_end3)


       if (middle3) then

          buffer_fine(1:first_bf1,1:first_bf2,mid_start3:end3_fine) = &
               cell_fine(start1_fine:mid_end1,start2_fine:mid_end2, &
               mid_start3:end3_fine)

          if (middle1) &
               buffer_fine(second_bf1:nf1,1:first_bf2,mid_start3:end3_fine) = &
               cell_fine(mid_start1:end1_fine,start2_fine:mid_end2, &
               mid_start3:end3_fine)

          if (middle2) &
               buffer_fine(1:first_bf1,second_bf2:nf2,mid_start3:end3_fine) = &
               cell_fine(start1_fine:mid_end1,mid_start2:end2_fine, &
               mid_start3:end3_fine)

          if (middle1 .and. middle2) &
               buffer_fine(second_bf1:nf1,second_bf2:nf2,mid_start3:end3_fine)=&
               cell_fine(mid_start1:end1_fine,mid_start2:end2_fine, &
               mid_start3:end3_fine)

       end if

       ! cks: --- END fill up FIRST AND SECOND parts of buffer------------------



       ! cks: >>>>>>>> MPI SEND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       ! cks: FIRST BLOCK
       block_size = mid_end3-start3_fine+1
       first_send_finished = .true.
       send_handles(1) = pub_null_handle
       if ((block_size > 0).and.(pub_my_node_id /= send_node)) then
          ! pdh: array slices passed to comms must be overtly contiguous
          ! pdh: to avoid copy in/copy out
          call comms_send(send_node, buffer_fine(:,:,start3_fine:mid_end3), &
               return_handle=send_handles(1),add_to_stack=.false.)
          first_send_finished = .false.
       end if


       ! cks: SECOND BLOCK
       block_size = end3_fine-mid_start3+1
       second_send_finished = .true.
       send_handles(2) = pub_null_handle
       if ((block_size > 0).and.(pub_my_node_id /= send_node)) then
          ! pdh: array slices passed to comms must be overtly contiguous
          ! pdh: to avoid copy in/copy out
          call comms_send(send_node, buffer_fine(:,:,mid_start3:end3_fine), &
               return_handle=send_handles(2),add_to_stack=.false.)
          second_send_finished = .false.
       end if

       ! cks: >>>>>>>> END MPI SEND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>





       ! cks: <<<<<<< MPI RECEIVE - DIRECTLY IN BOX <<<<<<<<<<<<<<<<<<<<


       ! cks: FIRST BLOCK
       block_size = my_limits(10,recv_node) - my_limits(9,recv_node) + 1
       if (block_size > 0) then

          if (recv_node == start3_node) then
             my_first3_start = 1
          else
             my_first3_start = grid%first_slab12(recv_node - &
                  modulo(recv_node,node_stride)) - &
                  grid%first_slab12(start3_node - &
                  modulo(recv_node,node_stride)) - &
                  my_limits(9,start3_node) + 2
          end if
          my_first3_end = my_first3_start + block_size - 1


          if ( pub_my_node_id /= recv_node) then
             ! pdh: array slices passed to comms must be overtly contiguous
             ! pdh: to avoid copy in/copy out
             call comms_recv(recv_node, &
                  box_fine(:,:,my_first3_start:my_first3_end))
          else
             box_fine(:,:,my_first3_start:my_first3_end) = &
                  buffer_fine(:,:,start3_fine:mid_end3)
          end if

       end if


       ! cks: SECOND_BLOCK
       block_size = my_limits(12,recv_node) - my_limits(11,recv_node) + 1
       if (block_size > 0) then

          my_second3_start = my_box_second3
          if (recv_node /= mid_start3_node) &
               my_second3_start = my_second3_start + &
               grid%first_slab12(recv_node - &
               modulo(recv_node,node_stride)) - 1
          my_second3_end = my_second3_start + block_size - 1


          if ( pub_my_node_id /= recv_node) then
             ! pdh: array slices passed to comms must be overtly contiguous
             ! pdh: to avoid copy in/copy out
             call comms_recv(recv_node, &
                  box_fine(:,:,my_second3_start:my_second3_end))
          else
             box_fine(:,:,my_second3_start:my_second3_end) = &
                  buffer_fine(:,:,mid_start3:end3_fine)
          end if

       end if

       ! cks: <<< END MPI RECEIVE - DIRECTLY IN BOX <<<<<<<<<<<<<<<<<<<<



       ! cks: All non-blocking sends should be completed before changing
       ! cks: the contents of buffer_fine
       ! ndmh: now implemented with comms_waitany and manual handle management
       do
          if (.not.first_send_finished) &
               call comms_test(first_send_finished,send_handles(1))
          if (.not.second_send_finished) &
               call comms_test(second_send_finished,send_handles(2))
          if (first_send_finished.and.second_send_finished) exit
          call comms_waitany(2,send_handles)
       end do

    end do
    !---------------------------------------------------------------------------
    ! ##### cks: END Fill up boxes with slab-12 potential
    !---------------------------------------------------------------------------

#ifdef ITC_TRACE
    call VTEND(vt_basis_extract_box_from_cell, vt_err)
#endif
    call timer_clock('cell_grid_extract_box', 2)

    deallocate(others_limits,stat=ierr)
    call utils_dealloc_check('cell_grid_extract_box','others_limits',ierr)
    deallocate(my_limits,stat=ierr)
    call utils_dealloc_check('cell_grid_extract_box','my_limits',ierr)

  end subroutine cell_grid_extract_box


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  subroutine cell_grid_box_limits(my_limits,others_limits, &
       start3_node,mid_end3_node,mid_start3_node,end3_node,my_box_second3, &
       i_have_box,grid,cell_start1,cell_start2,cell_start3,nf1,nf2,nf3, &
       group_comms)

    !======================================================================!
    ! This subroutine determines the limits of a set of boxes within a     !
    ! distributed whole-cell grid. Each node may have a box, and the grid  !
    ! point limits are determined for the communication by all nodes to    !
    ! all nodes, or optionally just within a group of nodes.               !
    !----------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine in December 2010, based !
    ! on older code from basis_mod written by Chris-Kriton Skylaris for    !
    ! the earlier cell extract/deposit code, and subsequently modified by  !
    ! Nicholas Hine.                                                       !
    !======================================================================!

    use comms, only: comms_alltoall, pub_rank_comm, pub_world_comm, &
         pub_comms_group_size, pub_total_num_nodes, &
         pub_my_rank_in_group

    ! Arguments
    integer, intent(out) :: my_limits(12,0:pub_total_num_nodes-1)
    integer, intent(out) :: others_limits(12,0:pub_total_num_nodes-1)
    integer, intent(out) :: start3_node     ! node which poseseses first slab of first block along direction 3
    integer, intent(out) :: mid_end3_node   ! node which poseseses last slab of first block along direction 3
    integer, intent(out) :: mid_start3_node ! node which poseseses first slab of second block along direction 3
    integer, intent(out) :: end3_node       ! node which poseseses last slab of second block along direction 3
    integer, intent(out) :: my_box_second3  ! start of second block of (whole of) box of pub_my_node_id
    logical, intent(in) :: i_have_box
    integer, intent(in) :: cell_start1,cell_start2,cell_start3
    integer, intent(in) :: nf1, nf2, nf3
    type(GRID_INFO), intent(in) :: grid
    logical, intent(in) :: group_comms

    ! Local Variables
    logical :: middle1, middle2, middle3  ! whether the box extends out of the sim cell in any dimension
    integer :: start1_fine, start2_fine, start3_fine ! start of first block of box wrt sim cell
    integer :: end1_fine, end2_fine, end3_fine       ! end of second block of box wrt sim cell
    integer :: mid_start1, mid_start2, mid_start3    ! start of second block of box wrt sim cell
    integer :: mid_end1, mid_end2, mid_end3          ! end of first block of box wrt sim cell
    integer :: first_node_in_group,last_node_in_group
    integer :: the_node                  ! node counter
    integer :: node_start, node_stride   ! node loop start and stride

    !----------------------------------------------------------------------
    ! ***** cks: Limits of box of pub_my_node_id wrt simulation cell
    !----------------------------------------------------------------------
    if (i_have_box) then

       ! cks: determine limits of box wrt sim cell along direction-1
       call cell_grid_box_wrt_cell_1d( &
            start1_fine, end1_fine, middle1, &
            cell_start1, grid%n1, nf1)

       mid_end1 = end1_fine
       mid_start1 = end1_fine + 1
       if (middle1) then
          mid_end1 = grid%n1
          mid_start1 = 1
       endif

       ! cks: determine limits of box wrt sim cell along direction-2
       call cell_grid_box_wrt_cell_1d( &
            start2_fine, end2_fine, middle2, &
            cell_start2, grid%n2, nf2)

       mid_end2 = end2_fine
       mid_start2 = end2_fine + 1
       if (middle2) then
          mid_end2 = grid%n2
          mid_start2 = 1
       endif

       ! cks: determine limits of box wrt sim cell along direction-3
       call cell_grid_box_wrt_cell_1d( &
            start3_fine, end3_fine, middle3, &
            cell_start3, grid%n3, nf3)

       mid_end3 = end3_fine
       mid_start3 = end3_fine + 1
       if (middle3) then
          mid_end3 = grid%n3
          mid_start3 = 1
       endif

       ! cks: find nodes (according to 12-slabs) to which the limits of
       ! cks: the box in direction 3 lie
       start3_node = grid%node_slab12(start3_fine)
       if (group_comms) start3_node = start3_node - &
            modulo(start3_node,pub_comms_group_size) + pub_my_rank_in_group
       mid_end3_node = grid%node_slab12(mid_end3)
       if (group_comms) mid_end3_node = mid_end3_node - &
            modulo(mid_end3_node,pub_comms_group_size) + pub_my_rank_in_group
       if (middle3) then
          mid_start3_node = grid%node_slab12(mid_start3)
          if (group_comms) mid_start3_node = mid_start3_node &
               - modulo(mid_start3_node,pub_comms_group_size) &
               + pub_my_rank_in_group
          end3_node = grid%node_slab12(end3_fine)
          if (group_comms) end3_node = end3_node &
               - modulo(end3_node,pub_comms_group_size) + pub_my_rank_in_group
       else
          ! cks: case where there is no second bit along direction-3
          mid_start3_node = pub_total_num_nodes
          end3_node       = -1
       endif

       my_box_second3 = mid_end3 - start3_fine + 2

    else
       !qoh: Initialisations to avoid compiler warnings
       start1_fine = grid%n1 + 1
       start2_fine = grid%n2 + 1
       start3_fine = grid%n3 + 1
       mid_start1  = grid%n1 + 1
       mid_start2  = grid%n2 + 1
       mid_start3  = grid%n3 + 1
       end1_fine = -1; end2_fine = -1; end3_fine = -1
       mid_end1 = -1; mid_end2 = -1; mid_end3 = -1
       start3_node = pub_total_num_nodes
       end3_node = -1
       mid_end3_node = -1
       mid_start3_node = pub_total_num_nodes
       my_box_second3 = -grid%n3
       middle3 = .false.
       middle2 = .false.
       middle1 = .false.
    endif
    !----------------------------------------------------------------------
    ! ***** cks: END Limits of box of pub_my_node_id wrt simulation cell
    !----------------------------------------------------------------------


    !-----------------------------------------------------------------------------------------
    ! ##### cks: Slab-12 limits in direction-3 of box of pub_my_node_id for all nodes
    !-----------------------------------------------------------------------------------------
    my_limits = 0
    node_start = 0
    if (group_comms) node_start = pub_my_rank_in_group
    node_stride = 1
    if (group_comms) node_stride = pub_comms_group_size
    do the_node=node_start,pub_total_num_nodes-1,node_stride

       if (i_have_box) then

          my_limits(1,the_node) = start1_fine
          my_limits(2,the_node) = mid_end1
          my_limits(3,the_node) = mid_start1
          my_limits(4,the_node) = end1_fine

          my_limits(5,the_node) = start2_fine
          my_limits(6,the_node) = mid_end2
          my_limits(7,the_node) = mid_start2
          my_limits(8,the_node) = end2_fine

          if (group_comms) then
             first_node_in_group = the_node - &
                  modulo(the_node,pub_comms_group_size)
             last_node_in_group = first_node_in_group + pub_comms_group_size - 1
          else
             first_node_in_group = the_node
             last_node_in_group = the_node
          end if

          ! cks: START3
          if (the_node == start3_node) then
             ! cks: exact slab where it starts
             my_limits(9,the_node) = start3_fine - &
                  grid%first_slab12(first_node_in_group) + 1
          elseif (the_node > start3_node) then
             ! cks: start at first slab of node group
             my_limits(9,the_node) = 1
          else
             ! cks: basically it does not start
             my_limits(9,the_node) = grid%n3 + 1
          endif


          ! cks: MID_END3
          if (the_node == mid_end3_node) then
             ! cks: exact slab of first end
             my_limits(10,the_node) = mid_end3 - &
                  grid%first_slab12(first_node_in_group) + 1
          elseif (the_node > mid_end3_node) then
             ! cks: basically it does not end
             my_limits(10,the_node) = -1
          else
             ! cks: end at last slab of node group
             my_limits(10,the_node) = grid%last_slab12(last_node_in_group) - &
                  grid%first_slab12(first_node_in_group) + 1
          endif


          ! cks: MID_START3
          if (the_node == mid_start3_node) then
             ! cks: exact slab where it starts
             my_limits(11,the_node) = mid_start3 - &
                  grid%first_slab12(first_node_in_group) + 1
          elseif (the_node > mid_start3_node) then
             ! cks: start at first slab of node group
             my_limits(11,the_node) = 1

          else
             ! cks: basically it does not start
             my_limits(11,the_node) = grid%n3 + 1
          endif


          ! cks: END3
          if (the_node == end3_node) then
             ! cks: exact slab of first end
             my_limits(12,the_node) = end3_fine - &
                  grid%first_slab12(first_node_in_group) + 1
          elseif (the_node > end3_node) then
             ! cks: basically it does not end
             my_limits(12,the_node) = -1
          else
             ! cks: end at last slab of node group
             my_limits(12, the_node) = grid%last_slab12(last_node_in_group) &
                  - grid%first_slab12(first_node_in_group) + 1
          endif


       else

          ! cks: I don't have a box to send
          my_limits(1, the_node) = grid%n1 + 1
          my_limits(2, the_node) = -1
          my_limits(3, the_node) = grid%n1 + 1
          my_limits(4, the_node) = -1

          my_limits(5, the_node) = grid%n2 + 1
          my_limits(6, the_node) = -1
          my_limits(7, the_node) = grid%n2 + 1
          my_limits(8, the_node) = -1

          my_limits(9, the_node) = grid%n3 + 1
          my_limits(10,the_node) = -1
          my_limits(11,the_node) = grid%n3 + 1
          my_limits(12,the_node) = -1

       endif

    enddo
    !-----------------------------------------------------------------------------------------
    ! ##### cks: END Slab-12 limits in direction-3 of box of pub_my_node_id for all nodes
    !-----------------------------------------------------------------------------------------

    ! ndmh: use alltoall to communicate limits
    others_limits = my_limits
    if (group_comms) then
       call comms_alltoall(others_limits,length=12*pub_comms_group_size, &
            comm=pub_rank_comm)
    else
       call comms_alltoall(others_limits)
    end if

contains

    subroutine cell_grid_box_wrt_cell_1d(start, finish, middle, &
         cell_start, cell_pts, box_pts)

      !=====================================================================!
      ! Returns the location of the start and end points of a box along one !
      ! dimension of a simulation cell grid, such that the origin of the    !
      ! box is always in the simulation cell and not in one of its periodic !
      ! images. Also returned is the logical element 'middle': value is     !
      ! 'true' if the box is split over the end of the simulation cell,     !
      ! otherwise false.                                                    !
      !---------------------------------------------------------------------!
      !                                                                     !
      !---------------------------------------------------------------------!
      ! Originally written by Chris-Kriton Skylaris on 16/10/2001 for the   !
      ! ONES program as "basis_pairbox_wrt_cell_1d_fine".                   !
      ! Adapted to "triple box" form by Arash A. Mostofi in September 2002  !
      ! Adapted from basis_superbox_wrt_cell_1d_fine by Nick Hine in March  !
      ! 2008 to allow basis_extract_box and basis_deposit_box to extract    !
      ! and deposit to points on the grid that do not coincide with         !
      ! points of the coarse grid                                           !
      !=====================================================================!

      implicit none

      ! Arguments
      integer, intent(out) :: start
      integer, intent(out) :: finish
      logical, intent(out) :: middle
      integer, intent(in)  :: cell_start
      integer, intent(in)  :: cell_pts
      integer, intent(in)  :: box_pts

      middle = .false.
      start = cell_start
      finish = start + box_pts - 1
      if (start<1) then
         start = start + cell_pts
         middle = .true.
      elseif (finish>cell_pts) then
         finish = finish - cell_pts
         middle = .true.
      endif

    end subroutine cell_grid_box_wrt_cell_1d

  end subroutine cell_grid_box_limits


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cell_grid_real_pt(rpt,i1,i2,i3,grid)

    !======================================================================!
    ! This subroutine returns the real-space position vector rvec of the   !
    ! (i,j,k)-th point of a given grid. Grids are distributed over nodes,  !
    ! so the k-index refers to the segment of the array local to this node !
    ! in real space. This subroutine replaces the array grid%real_grid,    !
    ! which was rather large in memory and rarely used.                    !
    !----------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine in January 2011.        !
    !======================================================================!

    use comms, only: pub_my_node_id

    ! Arguments
    real(kind=DP), intent(out) :: rpt(3)
    integer, intent(in) :: i1,i2,i3
    type(GRID_INFO), intent(in) :: grid

    ! Local Variables
    integer :: i3t

    i3t = i3 + grid%first_slab12(pub_my_node_id) - 1
    rpt(1) = (i3t-1)*grid%da3%x + (i2-1)*grid%da2%x + (i1-1)*grid%da1%x
    rpt(2) = (i3t-1)*grid%da3%y + (i2-1)*grid%da2%y + (i1-1)*grid%da1%y
    rpt(3) = (i3t-1)*grid%da3%z + (i2-1)*grid%da2%z + (i1-1)*grid%da1%z

  end subroutine cell_grid_real_pt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cell_grid_recip_pt(recpt,i1,i2,i3,grid)

    !======================================================================!
    ! This subroutine returns the reciprocal-space vector recpt of the     !
    ! (i,j,k)-th point of a given grid. Grids are distributed over nodes,  !
    ! so the k-index refers to the segment of the array local to this node !
    ! in real space. This subroutine replaces the array grid%recip_grid,   !
    ! which was rather large in memory.                                    !
    !----------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine in January 2011.        !
    !======================================================================!

    use utils, only: utils_abort

    ! Arguments
    real(kind=DP), intent(out) :: recpt(3)
    integer, intent(in) :: i1,i2,i3
    type(GRID_INFO), intent(in) :: grid

    ! Local Variables
    integer :: j2,j3

    ! ndmh: Sanity Check
    if ((i1<1).or.(i1>grid%n1/2+1)) then
       call utils_abort('Error in cell_grid_recip_pt: invalid i1 value')
    end if
    if ((i2<1).or.(i2>grid%n2)) then
       call utils_abort('Error in cell_grid_recip_pt: invalid i2 value')
    end if
    if ((i3<1).or.(i3>grid%n3)) then
       call utils_abort('Error in cell_grid_recip_pt: invalid i3 value')
    end if

    j2 = i2
    if (i2 > grid%n2/2+1) j2 = j2 - grid%n2
    j3 = i3
    if (i3 > grid%n3/2+1) j3 = j3 - grid%n3

    recpt(1) = (i1-1)*grid%b1%x + (j2-1)*grid%b2%x + (j3-1)*grid%b3%x
    recpt(2) = (i1-1)*grid%b1%y + (j2-1)*grid%b2%y + (j3-1)*grid%b3%y
    recpt(3) = (i1-1)*grid%b1%z + (j2-1)*grid%b2%z + (j3-1)*grid%b3%z

  end subroutine cell_grid_recip_pt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cell_grid_box_start_wrt_atom(cell_start1, cell_start2, &
       cell_start3,atom_origin,box_n1,box_n2,box_n3,grid)

    use comms, only: pub_my_node_id
    use constants, only: PI
    use geometry, only: POINT, OPERATOR(+), OPERATOR(-), OPERATOR(*), &
         OPERATOR(.DOT.)
    use simulation_cell, only: pub_cell
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(out) :: cell_start1
    integer, intent(out) :: cell_start2
    integer, intent(out) :: cell_start3
    type(POINT), intent(in) :: atom_origin
    integer, intent(in) :: box_n1
    integer, intent(in) :: box_n2
    integer, intent(in) :: box_n3
    type(GRID_INFO), intent(in) :: grid

    ! Local Variables
    type(POINT) :: box_origin
    type(POINT) :: half_box_span
    real(kind=DP) :: bs1,bs2,bs3
    real(kind=DP),parameter :: inv_two_pi = 0.5_DP / PI

    if ((box_n1>grid%n1).or.(box_n1<1).or. &
         (box_n2>grid%n2).or.(box_n2<1).or. &
         (box_n3>grid%n3).or.(box_n3<1)) then
       write(stdout,'(a,i6,a)') 'On node',pub_my_node_id,':'
       write(stdout,'(a,i5,a,i5)') 'box_n1: ', box_n1, &
            ', grid%n1: ',grid%n1
       write(stdout,'(a,i5,a,i5)') 'box_n2: ', box_n2, &
            ', grid%n2: ',grid%n2
       write(stdout,'(a,i5,a,i5)') 'box_n3: ', box_n3, &
            ', grid%n3: ',grid%n3
       call utils_abort('Error in cell_grid_box_start_wrt_atom: invalid &
            &box size')
    end if

    ! Find vector spanning half the box diagonal
    half_box_span = 0.5_DP*real(box_n1,kind=DP)*grid%da1 + &
                    0.5_DP*real(box_n2,kind=DP)*grid%da2 + &
                    0.5_DP*real(box_n3,kind=DP)*grid%da3

    ! Find vector to origin of box
    box_origin = atom_origin - half_box_span

    ! Calculate origin of box in terms of grid points
    bs1 = (box_origin.DOT.pub_cell%b1)*inv_two_pi*real(grid%n1,kind=DP)
    bs2 = (box_origin.DOT.pub_cell%b2)*inv_two_pi*real(grid%n2,kind=DP)
    bs3 = (box_origin.DOT.pub_cell%b3)*inv_two_pi*real(grid%n3,kind=DP)

    ! Grid point indices start at 1 rather than 0
    cell_start1 = floor(bs1) + 1
    cell_start2 = floor(bs2) + 1
    cell_start3 = floor(bs3) + 1

    ! Loop back into cell if origin is outside cell
    if (cell_start1<1) cell_start1 = cell_start1 + grid%n1
    if (cell_start2<1) cell_start2 = cell_start2 + grid%n2
    if (cell_start3<1) cell_start3 = cell_start3 + grid%n3
    if (cell_start1>grid%n1) cell_start1 = cell_start1 - grid%n1
    if (cell_start2>grid%n2) cell_start2 = cell_start2 - grid%n2
    if (cell_start3>grid%n3) cell_start3 = cell_start3 - grid%n3

  end subroutine cell_grid_box_start_wrt_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module cell_grid
