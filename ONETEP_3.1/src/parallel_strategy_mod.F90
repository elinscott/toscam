! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                     Parallel strategy module                   !
!----------------------------------------------------------------!
! Written by Peter Haynes, 10/7/03                               !
!================================================================!

module parallel_strategy


  use constants, only : DP, stdout, VERBOSE
  use geometry, only: POINT
  use ion, only : element

  implicit none

  private

  ! Public subroutines

  public :: parallel_strategy_distr_atoms
  public :: parallel_strategy_distr_funcs
  public :: parallel_strategy_check_atoms
  public :: parallel_strategy_list_overlaps
  public :: parallel_strategy_exit

  ! Public variables

  ! Information about an overlap list created by the module
  integer, public :: pub_max_overlaps                   ! Max # overlaps
  integer, public, allocatable :: pub_num_overlaps(:)   ! Num overlaps
  integer, public, allocatable :: pub_overlap_list(:,:) ! Overlap list

  ! Information about distribution of atoms over nodes
  integer, public :: pub_max_atoms_on_node  ! Max # atoms on any node

  ! ndmh: public arrays giving information on the distribution of atoms
  integer, public, allocatable :: pub_num_atoms_on_node(:)  ! Num atoms on node
  integer, public, allocatable :: pub_first_atom_on_node(:) ! First atom on node
  integer, public, allocatable :: pub_node_of_atom(:)       ! Node of each atom

  ! ddor: public arrays giving information on the distribution of Hubbard atoms
  integer, public, allocatable :: pub_num_hub_atoms_on_node(:)  ! Num Hubbard atoms on node
  integer, public, allocatable :: pub_node_of_hub_atom(:)       ! Node of each Hubbard atom
  ! ddor: List of Hubbard atoms on all nodes
  integer, public, allocatable, dimension(:,:) :: pub_hub_atoms_on_node
  ! ddor: An array hub_hub_nat long, ordered as the input file, 
  !       which contains the original atom number of the Hubbard atom.
  integer, public, allocatable, dimension(:) :: pub_hub_atom_orig_atom

  ! Node where atom is to be found:
  ! List of atoms on all nodes
  integer, public, allocatable, dimension(:,:) :: pub_atoms_on_node
  ! cks: the new order of the atoms after distribution to nodes
  ! cks: i.e first come all atoms of node0, then all atoms of node1, etc.
  integer, public, allocatable, dimension(:) :: pub_orig_atom
  ! pdh: the opposite of pub_orig_atom
  integer, public, allocatable, dimension(:) :: pub_distr_atom
  ! List of atoms on this node
  type(element), public, allocatable, dimension(:) :: pub_elements_on_node

  ! Private variables

  logical :: atoms_distributed = .false.               ! Flag whether atoms
                                                       ! have been distributed

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_distr_atoms(elements)

    !=========================================================================!
    ! This subroutine sorts a list of elements according to a Hilbert space-  !
    ! filling curve so as to maximise the overlap between spheres on the same !
    ! node, and then distributes the atoms across the nodes.                  !
    !-------------------------------------------------------------------------!
    ! Arguments:							      !
    !   elements (input) : The list of elements read from the input file      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:				              !
    !   pub_max_atoms_on_node  : Maximum number of atoms on any node          !
    !   pub_num_atoms_on_node  : Number of atoms on each node                 !
    !   pub_atoms_on_node      : A list of atoms on each node                 !
    !   pub_elements_on_node   : A list of elements on this node              !
    !   pub_first_atom_on_node : First atom on each node                      !
    !   pub_first_ngwf_on_node : First NGWF on each node                      !
    !-------------------------------------------------------------------------!
    ! Modules used:						              !
    !   constants       : For pi                                              !
    !   geometry        : For the point type and dot product                  !
    !   simulation_cell : For the pub_cell type    			      !
    !   ion             : For the element type                                !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:					              !
    !   fcoord   : list of fractional coordinates of the atoms wrt lattice    !
    !   gcoord   : list of atomic coordinates on a coarse grid                !
    !   graycode : the Gray codes for the atoms used for sorting according to !
    !              a space-filling curve                                      !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:						      !
    !   pub_cell%nat > 0                                                      !
    !   size(elements) == pub_cell%nat                                        !
    !   pub_comms_initialised == .true.                                       !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03	                                      !
    ! Modified by Quintin Hill to use pub_cell on 17/10/2008.                 !
    !=========================================================================!

    use comms, only: comms_abort, pub_comms_initialised, pub_on_root, &
        pub_my_node_id, pub_total_num_nodes
    use constants, only : PI, LONG
    use geometry, only: operator(.DOT.)
    use hubbard_init, only: h_species
    use ion, only: element
    use rundat, only: use_space_filling_curve, pub_output_detail, pub_hubbard
    use simulation_cell, only : pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort

    implicit none

    ! Arguments
    type(element), intent(in) :: elements(:)  ! Unsorted atomic information

    ! Local variables
    integer :: iat                            ! Atom loop counter
    integer :: hub_species
    integer :: orig_atom_count
    integer :: num
    integer :: ibit                           ! Bit loop counter
    integer :: node                           ! Node counter
    integer :: node_iat
    integer :: ngrid(3),log2grid(3),abc       ! Grid sizes and their log_2s
    integer :: global_atom                    ! Counter for atom in system
    integer :: hub_atom, find_hub_atom        ! ddor: Counter used for distributing Hubbard atoms
    integer :: atom_count                     ! Counter for atoms
    integer, allocatable :: gcoord(:,:)       ! Coordinates in terms of grid
    integer, allocatable, dimension(:) :: sort_index ! Sorting index for atoms
    integer, allocatable, dimension(:) :: nfuncs     ! Number of functions/atom
    integer, allocatable, dimension(:) :: num_funcs_on_node ! Num funcs/node
    integer(kind=LONG), allocatable :: graycode(:) ! Gray code for atoms

    real(kind=DP) :: recip_twopi              ! 1 / (2 pi)
    real(kind=DP) :: av_funcs_on_node         ! Average number of NGWFs per node
    real(kind=DP) :: excess_funcs
    real(kind=DP), allocatable :: fcoord(:,:) ! Atomic coordinates in fractions
                                              ! of lattice vectors
!CW
    logical             :: check
    integer,allocatable :: dimer_table(:,:)
    integer             :: k1,k2,k3
!END CW


!CW
!put dimers on same node for DMFT
INQUIRE(FILE='mask_dimer',EXIST=check)
if(check)then
 open(unit=22,file='mask_dimer')
 k1=0
 do
  read(22,*,end=21)
  k1=k1+1
 enddo
 21 continue
 rewind(22)
 if(allocated(dimer_table)) deallocate(dimer_table)
 allocate(dimer_table(k1,2))
 do k2=1,k1
  read(22,*) dimer_table(k2,1:2)
 enddo
 close(22)
endif
!END CW


    ! Check arguments
    if (pub_cell%nat <= 0) then
       if (pub_on_root) write(stdout,*) 'Error in parallel_strategy_distr_atoms: &
            &no atoms in cell!'
       call comms_abort
    end if
    if (size(elements) /= pub_cell%nat) then
       if (pub_on_root) write(stdout,*) 'Error in parallel_strategy_distr_atoms: &
            &incompatible nat'
       call comms_abort
    end if

    ! Check comms module intialised
    if (.not. pub_comms_initialised) stop &
         'Error in parallel_strategy_distr_atoms: &
         &comms module not initialised'

    ! Reallocate module arrays if necessary
    call internal_allocate_mod_1

    ! cks: the "use_space_filling_curve" is a logical that comes from the input
    if (pub_on_root .and. pub_output_detail == VERBOSE) then
         write(stdout,'(/a)',advance='no') '... '
         if (.not. use_space_filling_curve) write(stdout,'(a)',advance='no') 'not '
         write(stdout,'(a)') 'using space-filling curve'
      end if

    ! Allocate work arrays (and initialise sort index for atoms)
    call internal_allocate_work(use_space_filling_curve)

    !---------------------------------------------!
    ! Sort atoms according to space-filling curve !
    !  -- only if required --                     !
    !---------------------------------------------!

    if (use_space_filling_curve) then

       ! Convert Cartesian coordinates from cell into fractional coordinates in
       ! the interval [0,1)
       do iat=1,pub_cell%nat
          fcoord(1,iat) = elements(iat)%centre .dot. pub_cell%b1
          fcoord(2,iat) = elements(iat)%centre .dot. pub_cell%b2
          fcoord(3,iat) = elements(iat)%centre .dot. pub_cell%b3
       end do
       recip_twopi = 0.5_DP / pi
       fcoord = fcoord * recip_twopi
       fcoord = modulo(fcoord,1.0_DP)

       ! Set up a grid which will be used to 'coarsen' the atomic coordinates
       ! which can then be used in the sorting procedure - for simplicity choose
       ! the nearest power of 2 above the simulation cell coarse grid
       ngrid(1) = pub_cell%total_pt1
       ngrid(2) = pub_cell%total_pt2
       ngrid(3) = pub_cell%total_pt3
       do abc=1,3
          log2grid(abc) = int(log(real(ngrid(abc),kind=DP))/log(2.0_DP))+1
          ngrid(abc) = 2**log2grid(abc)
       end do

       ! 'Coarsen' the fractional coordinates onto this grid
       do iat=1,pub_cell%nat
          do abc=1,3
             gcoord(abc,iat) = int(fcoord(abc,iat)*ngrid(abc))
          end do
       end do

       ! Construct Gray code for each atom (interlace binary digits from the
       ! three dimensions) - use long integers to allow grid sizes > 1024
       graycode = 0_LONG

       do iat=1,pub_cell%nat   ! Loop over all atoms in cell

          do ibit=0,log2grid(1)-1
             if (btest(gcoord(1,iat),ibit)) graycode(iat) = &
                  ior(graycode(iat),ishft(4_LONG,int(3*ibit,kind=LONG)))
          end do
          do ibit=0,log2grid(2)-1
             if (btest(gcoord(2,iat),ibit)) graycode(iat) = &
                  ior(graycode(iat),ishft(2_LONG,int(3*ibit,kind=LONG)))
          end do
          do ibit=0,log2grid(3)-1
             if (btest(gcoord(3,iat),ibit)) graycode(iat) = &
                  ior(graycode(iat),ishft(1_LONG,int(3*ibit,kind=LONG)))
          end do


       end do  ! End of loop over all atoms

       ! Sort atoms by their Gray codes
       call utils_heapsort(pub_cell%nat,graycode,sort_index)

    end if

    !-----------------------------------------------------------------------!
    ! Distribute atoms across nodes, balancing number of funcs on each node !
    !-----------------------------------------------------------------------!


    ! Fill nfuncs array with number of functions on each atom
    ! Can change this value to distribute by other types of function
    do iat=1,pub_cell%nat
       nfuncs(iat) = elements(iat)%nfunctions
    end do

    num = sum(nfuncs(:))

    ! Calculate expected number of functions per node
    av_funcs_on_node = num / real(pub_total_num_nodes,kind=DP)

    ! Work out number of atoms on each node - aim is to balance the number
    ! of NGWFs per node to balance the work
    node = 0                  ! Stores current node being filled up
    num_funcs_on_node = 0     ! Counts number of functions on each node
    pub_num_atoms_on_node = 0 ! Counts number of atoms on each node

    ! Put first atom on node 0
    pub_num_atoms_on_node(0) = 1
    num_funcs_on_node(0) = nfuncs(sort_index(1))

    ! Now loop over the remaining atoms
    do iat=2,pub_cell%nat

       ! Calculate what this node's surplus/deficit of NGWFs is
       excess_funcs = num_funcs_on_node(node)-av_funcs_on_node

       ! If it's a surplus, move on to the next node (unless we've run out)
       if (abs(excess_funcs) < abs(excess_funcs+nfuncs(sort_index(iat))) &
            .and. node < pub_total_num_nodes-1) node = node+1

       ! Put the atom on the current node
       num_funcs_on_node(node) = num_funcs_on_node(node) + &
            nfuncs(sort_index(iat))
       pub_num_atoms_on_node(node) = pub_num_atoms_on_node(node) + 1

       ! Update the expected average to take history into account:
       ! remaining number of NGWFs / remaining number of nodes
       if (node > 0) av_funcs_on_node = &
            (num-sum(num_funcs_on_node(0:node-1))) / &
            real(pub_total_num_nodes-node,kind=DP)

    end do  ! End of loop over atoms

    ! Find maximum number of atoms on any one node
    pub_max_atoms_on_node = maxval(pub_num_atoms_on_node)

    ! Initialise array indicating first atom on each node
    global_atom = 1
    do node=0,pub_total_num_nodes-1
       pub_first_atom_on_node(node) = global_atom
       global_atom = global_atom + pub_num_atoms_on_node(node)
    end do
    pub_first_atom_on_node(pub_total_num_nodes) = pub_cell%nat + 1

    ! Reallocate module arrays which depend upon pub_max_atoms_on_node
    ! if necessary
    call internal_allocate_mod_2

    ! Make up list of atoms on each node
    node = 0
    node_iat = 0
    do iat=1,pub_cell%nat
       node_iat = node_iat + 1

       pub_node_of_atom(sort_index(iat)) = node
       pub_atoms_on_node(node_iat,node) = sort_index(iat)

       ! Advance to next node if required
       if (node_iat==pub_num_atoms_on_node(node)) then
          node = node + 1
          node_iat = 0
       end if

    end do

    ! cks: Global order of atoms as distributed on nodes.
    ! cks: This is the order in which they will be arranged in matrices, etc.
    atom_count =0
    do node=0,pub_total_num_nodes-1
       do iat=1,pub_num_atoms_on_node(node)
          atom_count = atom_count+1

          orig_atom_count           = pub_atoms_on_node(iat,node)
          pub_orig_atom(atom_count) = orig_atom_count
          pub_distr_atom(orig_atom_count) = atom_count


       end do
    end do

    ! Make up list of atoms on this node
    do iat=1,pub_num_atoms_on_node(pub_my_node_id)
       pub_elements_on_node(iat) = &
           elements(pub_atoms_on_node(iat,pub_my_node_id))
    end do

    ! ddor: Distribute Hubbard DFT+U atoms and projectors over nodes if necessary,
    !       Hubbard atoms are distributed as per their host atoms,
    !       and it is necessary to have pub_orig_atom at this stage.
    if (pub_hubbard) then
       
       hub_atom = 0
       do hub_species=1,pub_cell%num_hub_species
          do iat=1,pub_cell%nat
             if (elements(iat)%species_id == &
                  h_species(hub_species)%hub_species) then
                hub_atom = hub_atom + 1
                pub_hub_atom_orig_atom(hub_atom) = iat
             end if
          end do
       end do

       pub_num_hub_atoms_on_node = 0

       do node=0,pub_total_num_nodes-1

          ! ddor: Loop over Hubbard atoms
          do hub_atom = 1,pub_cell%nat_hub
             ! Loop over atoms on this node
             do find_hub_atom=pub_first_atom_on_node(node), &
                  pub_first_atom_on_node(node+1)-1
                ! If it's the Hubbard atom we want
                if ( pub_hub_atom_orig_atom(hub_atom) == pub_orig_atom(find_hub_atom) ) then
                   pub_num_hub_atoms_on_node(node) = pub_num_hub_atoms_on_node(node) + 1

                   pub_node_of_hub_atom(hub_atom) = node
                   pub_hub_atoms_on_node(pub_num_hub_atoms_on_node(node),node) = hub_atom
                endif
             enddo
          enddo

       end do

    endif

    !CW
    if(check) deallocate(dimer_table)
    !END CW

    ! Deallocate work arrays
    call internal_deallocate_work(use_space_filling_curve)

    ! Set flag to show atoms have been distributed
    atoms_distributed = .true.

  contains

    !--------------------------------------------------------------------------

    subroutine internal_allocate_mod_1

      !=======================================================================!
      ! This subroutine reallocates module arrays as required by the parent   !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:							      !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      ! Reorganised by Nicholas Hine, 18/05/09                                !
      !=======================================================================!

      implicit none

      ! Local variables
      integer :: ierr     ! Error flag

      ! Deallocate module arrays if necessary
      if (allocated(pub_num_atoms_on_node)) then
         if (size(pub_num_atoms_on_node) /= pub_total_num_nodes) then
            deallocate(pub_num_atoms_on_node,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','pub_num_atoms_on_node',ierr)
         end if
      end if
      if (allocated(pub_first_atom_on_node)) then
         if (size(pub_first_atom_on_node) /= pub_total_num_nodes+1) then
            deallocate(pub_first_atom_on_node,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','pub_first_atom_on_node', &
                 ierr)
         end if
      end if
      if (allocated(pub_node_of_atom)) then
         if (size(pub_node_of_atom) /= pub_cell%nat) then
            deallocate(pub_node_of_atom,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','pub_node_of_atom',ierr)
         end if
      end if
      if (allocated(pub_orig_atom)) then
         if (size(pub_orig_atom) /= pub_cell%nat) then
            deallocate(pub_orig_atom,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','pub_orig_atom',ierr)
         end if
      end if
      if (allocated(pub_distr_atom)) then
         if (size(pub_distr_atom) /= pub_cell%nat) then
            deallocate(pub_distr_atom,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','pub_distr_atom',ierr)
         end if
      end if
      if (pub_hubbard) then ! ddor
         if (allocated(pub_num_hub_atoms_on_node)) then
            if (size(pub_num_hub_atoms_on_node) /= pub_total_num_nodes) then
               deallocate(pub_num_hub_atoms_on_node,stat=ierr)
               call utils_dealloc_check('internal_allocate_mod_1 &
                    &(parallel_strategy_distr_atoms)','pub_num_hub_atoms_on_node',ierr)
            end if
         end if
         if (allocated(pub_node_of_hub_atom)) then
            if (size(pub_node_of_hub_atom) /= pub_cell%nat_hub) then
               deallocate(pub_node_of_hub_atom,stat=ierr)
               call utils_dealloc_check('internal_allocate_mod_1 &
                    &(parallel_strategy_distr_atoms)','pub_node_of_hub_atom',ierr)
            end if
         end if
         if (allocated(pub_hub_atoms_on_node)) then
            if (size(pub_hub_atoms_on_node,1) /= pub_cell%nat_hub .or. &
                 size(pub_hub_atoms_on_node,2) /= pub_total_num_nodes) then
               deallocate(pub_hub_atoms_on_node,stat=ierr)
               call utils_dealloc_check('internal_allocate_mod_2 &
                    &(parallel_strategy_distr_atoms)','pub_hub_atoms_on_node',ierr)
            end if
         end if
         if (allocated(pub_hub_atom_orig_atom)) then
            if (size(pub_hub_atom_orig_atom) /= pub_cell%nat_hub) then
               deallocate(pub_hub_atom_orig_atom,stat=ierr)
               call utils_dealloc_check('internal_allocate_mod_1 &
                   &(parallel_strategy_distr_atoms)','pub_hub_atom_orig_atom',ierr)
            end if
         end if
      end if
      ! Allocate module arrays if necessary
      if (.not.allocated(pub_num_atoms_on_node)) then
         allocate(pub_num_atoms_on_node(0:pub_total_num_nodes-1),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','pub_num_atoms_on_node',ierr)
      end if
      if (.not.allocated(pub_first_atom_on_node)) then
         allocate(pub_first_atom_on_node(0:pub_total_num_nodes),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','pub_first_atom_on_node',ierr)
      end if
      if (.not.allocated(pub_node_of_atom)) then
         allocate(pub_node_of_atom(pub_cell%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','pub_node_of_atom',ierr)
      end if
      if (.not.allocated(pub_orig_atom)) then
         allocate(pub_orig_atom(pub_cell%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','pub_orig_atom',ierr)
      end if
      if (.not.allocated(pub_distr_atom)) then
         allocate(pub_distr_atom(pub_cell%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_1 &
              &(parallel_strategy_distr_atoms)','pub_distr_atom',ierr)
      end if
      if (pub_hubbard) then ! ddor
         if (.not.allocated(pub_num_hub_atoms_on_node)) then
            allocate(pub_num_hub_atoms_on_node(0:pub_total_num_nodes-1),stat=ierr)
            call utils_alloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','pub_num_hub_atoms_on_node',ierr)
         end if
         if (.not.allocated(pub_node_of_hub_atom)) then
            allocate(pub_node_of_hub_atom(pub_cell%nat_hub),stat=ierr)
            call utils_alloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','pub_node_of_hub_atom',ierr)
         end if
         if (.not. allocated(pub_hub_atoms_on_node)) then
            allocate(pub_hub_atoms_on_node(pub_cell%nat_hub, &
                 0:pub_total_num_nodes-1),stat=ierr)
            call utils_alloc_check('internal_allocate_mod_2 &
                 &(parallel_strategy_distr_atoms)','pub_hub_atoms_on_node',ierr)
            ! ddor: initialise to something obvious
            pub_hub_atoms_on_node = -1
         end if
         if (.not.allocated(pub_hub_atom_orig_atom)) then
            allocate(pub_hub_atom_orig_atom(pub_cell%nat_hub),stat=ierr)
            call utils_alloc_check('internal_allocate_mod_1 &
                 &(parallel_strategy_distr_atoms)','pub_hub_atom_orig_atom',ierr)
         end if
      end if

    end subroutine internal_allocate_mod_1

    !--------------------------------------------------------------------------

    subroutine internal_allocate_work(sfc)

      !=======================================================================!
      ! This subroutine allocates work arrays as required by the parent       !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:							      !
      ! sfc (input) : flag whether to use space-filling curve or not          !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      implicit none

      ! Arguments
      logical, intent(in) :: sfc

      ! Local variables
      integer :: ierr     ! Error flag
      integer :: iat      ! Atom loop variable

      ! Allocate the index for the sorting of the atoms
      allocate(sort_index(pub_cell%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_work &
           &(parallel_strategy_distr_atoms)','sort_index',ierr)
      allocate(nfuncs(pub_cell%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_work &
           &(parallel_strategy_distr_atoms)','nfuncs',ierr)
      allocate(num_funcs_on_node(0:pub_total_num_nodes-1),stat=ierr)
      call utils_alloc_check('internal_allocate_work &
           &(parallel_strategy_distr_atoms)','num_funcs_on_node',ierr)

      ! Allocate work arrays for space-filling curve if required
      if (sfc) then
         allocate(gcoord(3,pub_cell%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_work &
              &(parallel_strategy_distr_atoms)','gcoord',ierr)
         allocate(graycode(pub_cell%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_work &
              &(parallel_strategy_distr_atoms)','graycode',ierr)
         allocate(fcoord(3,pub_cell%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_work &
              &(parallel_strategy_distr_atoms)','fcoord',ierr)
      else
         ! cks: if no space-filling curve is used, the
         ! cks: sorting index needs to be initialised
         do iat=1, pub_cell%nat
            sort_index(iat)=iat
         end do
      endif


    end subroutine internal_allocate_work

    !--------------------------------------------------------------------------

    subroutine internal_allocate_mod_2

      !=======================================================================!
      ! This subroutine reallocates module arrays as required by the parent   !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:							      !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      implicit none

      ! Local variables
      integer :: ierr     ! Error flag

      ! Deallocate module arrays which depend upon pub_max_atoms_on_node
      ! if necessary
      if (allocated(pub_atoms_on_node)) then
         if (size(pub_atoms_on_node,1) /= pub_max_atoms_on_node .or. &
              size(pub_atoms_on_node,2) /= pub_total_num_nodes) then
            deallocate(pub_atoms_on_node,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_2 &
                 &(parallel_strategy_distr_atoms)','pub_atoms_on_node',ierr)
         end if
      end if
      if (allocated(pub_elements_on_node)) then
         if (size(pub_elements_on_node) /= pub_max_atoms_on_node) then
            deallocate(pub_elements_on_node,stat=ierr)
            call utils_dealloc_check('internal_allocate_mod_2 &
                 &(parallel_strategy_distr_atoms)','pub_elements_on_node',ierr)
         end if
      end if

      ! Allocate module arrays which depend upon pub_max_atoms_on_node
      ! if necessary
      if (.not. allocated(pub_atoms_on_node)) then
         allocate(pub_atoms_on_node(pub_max_atoms_on_node, &
              0:pub_total_num_nodes-1),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_2 &
              &(parallel_strategy_distr_atoms)','pub_atoms_on_node',ierr)
         ! cks: initialise to something obvious
         pub_atoms_on_node = -1
      end if
      if (.not. allocated(pub_elements_on_node)) then
         allocate(pub_elements_on_node(pub_max_atoms_on_node),stat=ierr)
         call utils_alloc_check('internal_allocate_mod_2 &
              &(parallel_strategy_distr_atoms)','pub_elements_on_node',ierr)
      end if

    end subroutine internal_allocate_mod_2

    !--------------------------------------------------------------------------

    subroutine internal_deallocate_work(sfc)

      !=======================================================================!
      ! This subroutine deallocates work arrays as required by the parent     !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:							      !
      ! sfc (input) : flag whether to use space-filling curve or not          !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      implicit none

      ! Arguments
      logical, intent(in) :: sfc

      ! Local variables
      integer :: ierr     ! Error flag

      ! Deallocate work arrays required by space-filling-curve
      if (sfc) then
         deallocate(fcoord,stat=ierr)
         call utils_dealloc_check('internal_deallocate_work &
              &(parallel_strategy_distr_atoms)','fcoord',ierr)
         deallocate(graycode,stat=ierr)
         call utils_dealloc_check('internal_deallocate_work &
              &(parallel_strategy_distr_atoms)','graycode',ierr)
         deallocate(gcoord,stat=ierr)
         call utils_dealloc_check('internal_deallocate_work &
              &(parallel_strategy_distr_atoms)','gcoord',ierr)
      end if

      ! Deallocate sorting index
      deallocate(num_funcs_on_node,stat=ierr)
      call utils_dealloc_check('internal_deallocate_work &
           &(parallel_strategy_distr_atoms)','num_funcs_on_node',ierr)
      deallocate(nfuncs,stat=ierr)
      call utils_dealloc_check('internal_deallocate_work &
           &(parallel_strategy_distr_atoms)','nfuncs',ierr)
      deallocate(sort_index,stat=ierr)
      call utils_dealloc_check('internal_deallocate_work &
           &(parallel_strategy_distr_atoms)','sort_index',ierr)

    end subroutine internal_deallocate_work

  end subroutine parallel_strategy_distr_atoms

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_distr_funcs(num, nfuncs_orig, sort_index, &
       first_func_on_node, num_funcs_on_node, first_func_on_atom, &
       num_funcs_on_atom, node_of_func, atom_of_func, max_funcs_on_node, &
       max_funcs_on_atom)

    !=========================================================================!
    ! This subroutine takes a list of numbers of functions on each atom and a !
    ! sorting index and fills the arrays describing the distribution of       !
    ! functions over nodes and atoms.                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:							      !
    !   num (input)         : The total number of functions being distributed.!
    !   nfuncs_orig (input) : The number of functions on each atom in the     !
    !                         order given in the input file / elements array. !
    !   sort_index (input)  : The index of atoms after sorting over over the  !
    !                         nodes of the calculation.                       !
    !   first_func_on_node,                                                   !
    !   num_funcs_on_node,                                                    !
    !   first_func_on_atom,                                                   !
    !   num_funcs_on_atom,                                                    !
    !   node_of_func,                                                         !
    !   atom_of_func (output) : Arrays describing distribution of functions   !
    !                           over nodes and atoms (meanings as above).     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:				              !
    !   pub_first_atom_on_node : First atom on each node                      !
    !-------------------------------------------------------------------------!
    ! Modules used:						              !
    !   comms           : For pub_total_num_nodes                             !
    !   simulation_cell : For the pub_cell type    			      !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:						      !
    !   pub_cell%nat > 0                                                      !
    !   pub_first_atom_on_node to be initialised and filled                   !
    !   pub_comms_initialised == .true.                                       !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine 9th July 2009.                            !
    !=========================================================================!

    use comms, only: pub_total_num_nodes
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    integer, intent(in) :: num
    integer, intent(in) :: nfuncs_orig(1:pub_cell%nat)
    integer, intent(in) :: sort_index(1:pub_cell%nat)
    integer, intent(out) :: first_func_on_node(0:pub_total_num_nodes)
    integer, intent(out) :: num_funcs_on_node(0:pub_total_num_nodes-1)
    integer, intent(out) :: first_func_on_atom(1:pub_cell%nat)
    integer, intent(out) :: num_funcs_on_atom(1:pub_cell%nat)
    integer, intent(out) :: node_of_func(1:num)
    integer, intent(out) :: atom_of_func(1:num)
    integer, intent(out) :: max_funcs_on_node
    integer, intent(out) :: max_funcs_on_atom

    ! Locals
    integer :: node
    integer :: iat, ifunc, jfunc, nfuncs

    ! Work out number of functions on each node
    node = -1                 ! Stores current node being filled up
    num_funcs_on_node(:) = 0  ! Counts number of NGWFs on each node
    ifunc = 1

    ! Loop over the atoms in their sorted order
    do iat=1,pub_cell%nat
       nfuncs = nfuncs_orig(sort_index(iat))

       ! Check if we are on a new node
       if (iat == pub_first_atom_on_node(node+1)) then
          node = node + 1
          first_func_on_node(node) = ifunc
       end if

       ! Add these functions to the current node and atom
       num_funcs_on_node(node) = num_funcs_on_node(node) + nfuncs
       num_funcs_on_atom(iat) = nfuncs

       ! Record first function on this atom (if there are any functions)
       if (nfuncs > 0) then
          first_func_on_atom(iat) = ifunc
       else
          first_func_on_atom(iat) = 0
       end if

       ! Record the node and atom of these functions
       do jfunc=ifunc,ifunc+nfuncs-1
          node_of_func(jfunc) = node
          atom_of_func(jfunc) = iat
       end do

       ! Increment current function index
       ifunc = ifunc + nfuncs

    end do  ! End of loop over atoms

    ! Record last function + 1 at end of array
    first_func_on_node(pub_total_num_nodes) = num + 1

    ! Find maximum number of functions on any node
    max_funcs_on_node = maxval(num_funcs_on_node)

    ! Find maximum number of functions on any atom
    max_funcs_on_atom = maxval(num_funcs_on_atom)

  end subroutine parallel_strategy_distr_funcs

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_check_atoms(elements)

    use comms, only: comms_abort, comms_reduce, pub_my_node_id, pub_on_root
    use constants, only: DP, PI, stdout
    use geometry, only: operator(.DOT.), operator(*), operator(+), &
        geometry_distance
    use simulation_cell, only : pub_cell
    use ion, only: element
    use rundat, only: pub_check_atoms
    use utils, only: utils_flush

    !====================================================================!
    ! Make sure the positions of all atoms are inside the simulation     !
    ! cell and there no nearly-overlapping atoms.                        !
    !--------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000 as                        !
    ! services_check_atom_positions and revised on 30/5/2001.            !
    ! Modified by Chris-Kriton Skylaris on 12/10/2004 to also check      !
    ! for nearly-overlapping atoms.                                      !
    ! Modified by Nick Hine on 06/09/2008 to only check the atoms on     !
    ! this node, and moved to parallel_strategy_mod.F90                  !
    ! Modified by Nick Hine on 30/10/2008 to report which atoms are      !
    ! outside the unit cell in a more informative way.                   !
    ! Modified to prevent hanging when only some of the nodes detect an  !
    ! atom outside of simulation cell or two atoms too close, by Nick    !
    ! Hine on 18/09/2009.                                                !
    !====================================================================!

    implicit none

    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! cks: internal variables
    integer :: atom, orig_atom, in_planes
    integer :: row, local_row
    integer :: col, per_col
    integer :: per_col_start, per_col_end
    integer :: loc1
    integer :: loc2
    integer :: loc3
    integer :: abort_count
    real(kind =DP) :: direction_a1_a2_first, direction_a1_a2_second, &
         direction_a2_a3_first, direction_a2_a3_second, &
         direction_a3_a1_first, direction_a3_a1_second
    real(kind =DP) :: row_col_dist
#ifdef ACCELRYS
    real(kind=DP), parameter :: tolerance = 1.0e-10_DP
#else
    real(kind=DP), parameter :: tolerance = 0.0_DP
#endif

    abort_count = 0

    ! cks: make sure each atom is inside the simulation cell
    ! nmdh: only check atoms on this node
    do atom=1,pub_num_atoms_on_node(pub_my_node_id)

       ! cks: find direction of vector from atom centre to first and
       !      second a1-a2 plane
       direction_a1_a2_first= &
            -(pub_cell%b3.DOT.pub_elements_on_node(atom)%centre)
       direction_a1_a2_second=direction_a1_a2_first +2.0_DP*PI


       ! cks: find direction of vector from atom centre to first and
       !      second a2-a3 plane
       direction_a2_a3_first= &
            -(pub_cell%b1.DOT.pub_elements_on_node(atom)%centre)
       direction_a2_a3_second=direction_a2_a3_first +2.0_DP*PI


       ! cks: find direction of vector from atom centre to first and
       !      second a3-a1 plane
       direction_a3_a1_first= &
            -(pub_cell%b2.DOT.pub_elements_on_node(atom)%centre)
       direction_a3_a1_second=direction_a3_a1_first +2.0_DP*PI


       ! cks: now find if atom centre is between all three pairs of parallel
       !      planes that define the simulation cell
       in_planes=0
       if (direction_a1_a2_first*direction_a1_a2_second.le.tolerance) &
            in_planes=in_planes+1
       if (direction_a2_a3_first*direction_a2_a3_second.le.tolerance) &
            in_planes=in_planes+1
       if (direction_a3_a1_first*direction_a3_a1_second.le.tolerance) &
            in_planes=in_planes+1

       if (in_planes /= 3) then
          orig_atom = &
               pub_orig_atom(atom + pub_first_atom_on_node(pub_my_node_id) - 1)
          write(stdout,'(a,i6,a)') &
               'Error in services_check_atom_positions: atom number', &
               orig_atom,' is not in simulation cell.'
          write(stdout,'(a,3f24.12)') pub_elements_on_node(atom)%symbol, &
               pub_elements_on_node(atom)%centre%x, &
               pub_elements_on_node(atom)%centre%y, &
               pub_elements_on_node(atom)%centre%z
          call utils_flush
          abort_count = abort_count + 1
       end if

    enddo

    ! Abort if any pairs of atoms are too close
    call comms_reduce('SUM',abort_count)
    if (abort_count > 0) then
       if (pub_on_root) write(stdout,*) &
            'ERROR in parallel_strategy_check_atoms: ', abort_count, &
            ' atoms are not in simulation cell'
       call comms_abort
    end if

    abort_count = 0

    if (pub_check_atoms) then

       ! cks: check that there are no atoms that are almost overlapping
       ! ndmh: share checking equally over nodes:
       ! ndmh:   only check rows on this column, and check same number
       ! ndmh:   of cols regardless of row, by modulo of periodic col
       do local_row=1,pub_num_atoms_on_node(pub_my_node_id)
          row = local_row + pub_first_atom_on_node(pub_my_node_id) - 1

          ! nmdh: find range of cols to check for this row (periodic)
          per_col_start = row + 1
          per_col_end = row + pub_cell%nat/2
          ! ndmh: avoid double counting for even-valued nat
          if (per_col_end > pub_cell%nat) then
             per_col_end = per_col_end - 1 + modulo(pub_cell%nat,2)
          end if
          do per_col=per_col_start,per_col_end

             ! ndmh: find real value of col
             col = modulo(per_col - 1,pub_cell%nat) + 1

             ! cks: initialise for current_pair
             row_col_dist =huge(1.0_DP)

             ! cks: 27 cases need to be considered for each distinct
             ! cks: pair of atoms because of periodic boundary conditions.
             do loc1 =-1, 1
                do loc2 =-1, 1
                   do loc3 =-1, 1

                      ! cks: accumulate minimum distance
                      row_col_dist=min(geometry_distance(elements(row)%centre, &
                           elements(col)%centre &
                           + real(loc1, kind =DP)*pub_cell%a1 &
                           + real(loc2, kind =DP)*pub_cell%a2 &
                           + real(loc3, kind =DP)*pub_cell%a3 ), row_col_dist)

                      if (row_col_dist <= 0.5_DP) then

                         write(stdout,*)'ERROR! The distance between atoms', &
                              row,' and ',col
                         write(stdout,*)'is ', row_col_dist, ' bohr. '
                         abort_count = abort_count + 1

                      endif

                   enddo
                enddo
             enddo

          enddo
       enddo

       ! Abort if any pairs of atoms are too close
       call comms_reduce('SUM',abort_count)
       if (abort_count > 0) then
          if (pub_on_root) write(stdout,*) &
              'ERROR in parallel_strategy_check_atoms: ',  abort_count, &
              ' pairs of atoms are too close'
          call comms_abort
       end if

    end if


  end subroutine parallel_strategy_check_atoms

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_exit

    !=========================================================================!
    ! This subroutine frees up any memory used by the parallel strategy       !
    ! module.                                                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:							      !
    !   None                						      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:				              !
    !   pub_num_overlaps      : the number of overlaps for each atom          !
    !   pub_overlap_list      : the list of atoms overlapping each atom       !
    !   pub_num_atoms_on_node : the number of atoms on each node              !
    !   pub_node_of_atom      : the node on which each atom is held           !
    !   pub_atoms_on_node     : list of atoms on each node                    !
    !   pub_elements_on_node  : the elements on this node                     !
    !-------------------------------------------------------------------------!
    ! Modules used:						              !
    !   utils                						      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:					              !
    !   None                						      !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:						      !
    !   None                						      !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/07/03                                       !
    ! Reorganised by Nicholas Hine, 18/05/09                                  !
    !=========================================================================!

    use rundat, only: pub_fine_is_dbl
    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr ! Error flag

    ! Deallocate number of overlaps
    if (allocated(pub_num_overlaps)) then
       deallocate(pub_num_overlaps,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_num_overlaps',ierr)
    end if

    ! Deallocate overlap list
    if (allocated(pub_overlap_list)) then
       deallocate(pub_overlap_list,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_overlap_list',ierr)
    end if

    ! Deallocate list of atoms on this node
    if (allocated(pub_atoms_on_node)) then
       deallocate(pub_atoms_on_node,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_atoms_on_node',ierr)
    end if
    ! Deallocate list of elements on this node
    if (allocated(pub_elements_on_node)) then
       deallocate(pub_elements_on_node,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_elements_on_node',ierr)
    end if

    !ddor: Deallocate information on DFT+U atom distribution
    if (allocated(pub_hub_atom_orig_atom)) then
       deallocate(pub_hub_atom_orig_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_hub_atom_orig_atom',ierr)
    end if
    if (allocated(pub_hub_atoms_on_node)) then
       deallocate(pub_hub_atoms_on_node,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_hub_atoms_on_node',ierr)
    end if
    if (allocated(pub_node_of_hub_atom)) then
       deallocate(pub_node_of_hub_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_node_of_hub_atom',ierr)
    end if
    if (allocated(pub_num_hub_atoms_on_node)) then
       deallocate(pub_num_hub_atoms_on_node,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_num_hub_atoms_on_node',ierr)
    end if

    ! cks: Deallocate pub_orig_atom and pub_distr_atom on this node
    if (allocated(pub_distr_atom)) then
       deallocate(pub_distr_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_distr_atom',ierr)
    end if
    if (allocated(pub_orig_atom)) then
       deallocate(pub_orig_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_orig_atom',ierr)
    end if
    ! Deallocate information on atom distribution
    if (allocated(pub_node_of_atom)) then
       deallocate(pub_node_of_atom,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_node_of_atom',ierr)
    end if
    if (allocated(pub_first_atom_on_node)) then
       deallocate(pub_first_atom_on_node,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_first_atom_on_node',ierr)
    end if
    if (allocated(pub_num_atoms_on_node)) then
       deallocate(pub_num_atoms_on_node,stat=ierr)
       call utils_dealloc_check('parallel_strategy_exit', &
            'pub_num_atoms_on_node',ierr)
    end if

  end subroutine parallel_strategy_exit

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine parallel_strategy_list_overlaps(elements,mode1,mode2,radius)

    !=========================================================================!
    ! This subroutine calculates the list of overlaps between atoms, either   !
    ! according to the region/core radii or according to an optional fixed    !
    ! radius.                                                                 !
    ! The form of overlap list required is determined by the two mode         !
    ! arguments which can be any of the following options:                    !
    !   'C' - use nonlocal pseudopotential Core radii                         !
    !   'R' - use NGWF Region radii                                           !
    !   'A' - use Conduction NGWF Region radii                                !
    !   'F' - use Fixed radius                                                !
    ! mode1 determines the radius used for the primary atom (the index in the !
    ! module arrays) whereas mode2 determines the radius for the secondary    !
    ! atom (listed as overlapping with the primary atom in the list).         !
    ! A cell list algorithm is used to obtain linear scaling for large        !
    ! systems - these cells are referred to as subcells below to avoid        !
    ! confusion with the unit cell.                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:							      !
    !   elements (input) : A list of elements with positions and radii        !
    !   mode1 (input)    : Mode for primary atom (see above)                  !
    !   mode2 (input)    : Mode for secondary atom (see above)                !
    !   radius (input)   : The optional fixed radius to be used               !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:				              !
    !   pub_max_overlaps : The maximum number of overlaps of any atom         !
    !   pub_num_overlaps : The number of overlaps of each atom                !
    !   pub_overlap_list : The list of overlaps between atoms                 !
    !-------------------------------------------------------------------------!
    ! Modules used:						              !
    !   constants       : For pi                                              !
    !   geometry        : For the point type and dot product                  !
    !   simulation_cell : For the pub_cell type    			      !
    !   ion             : For the element type                                !
    !   utils           : For checking memory allocation/deallocation         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:					              !
    !   fcoord   : list of fractional coordinates of the atoms wrt lattice    !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:						      !
    !   pub_cell%nat > 0                                                      !
    !   size(elements) == pub_cell%nat                                        !
    !   mode1 and mode2 must be one of 'C','R' or 'F'                         !
    !   If either mode1 or mode2 == 'F', radius must be present               !
    !   Sphere radii should not be larger than unit cell                      !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03	                                      !
    ! Minor modification to reduce size of arrays by Nick Hine 13/03/08       !
    ! Modified by Quintin Hill to use pub_cell on 17/10/2008.                 !
    !=========================================================================!

    use comms, only: comms_abort, comms_bcast, comms_reduce, pub_on_root,&
        pub_my_node_id, pub_total_num_nodes
    use constants, only : PI, TWO_PI
    use geometry, only: operator(.DOT.)
    use ion, only: element
    use simulation_cell, only : pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(element), intent(in) :: elements(:)      ! List of atoms
    character, intent(in) :: mode1                ! Mode for primary atom
    character, intent(in) :: mode2                ! Mode for secondary atom
    real(kind=DP), intent(in), optional :: radius ! Optional fixed radius

    ! Local variables
    logical, parameter :: periodic(3) = (/.true.,.true.,.true./) ! Apply PBC
    integer :: ierr                           ! Error flag
    integer :: iat,jat                        ! Atom loop counters
    integer :: idim                           ! Dimension loop counter
    integer :: node                           ! Node counter
    integer :: isc,iscd(3),isc23              ! Subcell loop counters/labels
    integer :: jsc,jscd(3)                    ! Subcell loop counters/labels
    integer :: ksc,kscd(3),ksc23              ! Subcell loop counters/labels
    integer :: iatsc,jatsc                    ! Atom in subcell loop counters
    integer :: num_subcells(3)                ! Number of subcells in each dim
    integer :: total_num_subcells             ! Total number of subcells
    integer :: max_atoms_subcell              ! Max num atoms in any subcell
    integer :: max_my_overlaps                ! Max num overlaps on any node
    integer :: novat                          ! Total num overlaps for atom
    integer :: nyrat                          ! Num overlaps for atom on node
    integer :: iovlap                         ! Overlap loop counter
    real(kind=DP) :: recip_twopi              ! 1 / (2 pi)
    real(kind=DP) :: subcell_size             ! Minimum subcell size
    real(kind=DP) :: spacing(3)               ! Plane spacing of unit cell faces
    real(kind=DP) :: a1(3),a2(3),a3(3)        ! Local copy of lattice vectors
    real(kind=DP) :: fextra(3)                ! Wrap-around for PBC
    real(kind=DP) :: fdiff(3)                 ! Fractional coordinate diffs
    real(kind=DP) :: adiff(3)                 ! Absolute Cartesian diffs
    real(kind=DP) :: dist_sq                  ! Squared distance between atoms
    real(kind=DP) :: cutoff                   ! Cutoff distance
    integer, allocatable :: atom_subcell(:)   ! Which subcell an atom belongs to
    integer, allocatable :: num_atoms_subcell(:) ! Num atoms in each subcell
    integer, allocatable :: subcell_list(:,:) ! Subcell list of atoms
    integer, allocatable :: num_my_overlaps(:)! Num overlaps on this node
    integer, allocatable :: num_yr_overlaps(:)! Num overlaps on other node
    integer, allocatable :: my_pub_overlap_list(:,:) ! Overlaps on this node
    integer, allocatable :: yr_pub_overlap_list(:,:) ! Overlaps on other node
    logical, allocatable :: overlapped(:)     ! Flags to avoid multiple overlaps
    real(kind=DP), allocatable :: fcoord(:,:) ! Atomic coordinates in fractions
                                              ! of lattice vectors
    real(kind=DP), allocatable :: radii(:,:)  ! Radii to be used for primary and
                                              ! secondary atoms
#ifdef ACCELRYS
    real(kind=DP), parameter :: tolerance = 1.0e-10_DP
#endif

    ! Check arguments
    call internal_check_args(present(radius))

    ! Allocate module and work arrays
    call internal_allocate_1

    ! Make up lists of radii for primary and secondary atoms according to
    ! the desired modes
    select case (mode1)
    case('C','c')
       do iat=1,pub_cell%nat
          radii(1,iat) = elements(iat)%max_core_radius
       end do
    case('R','r')
       do iat=1,pub_cell%nat
          radii(1,iat) = elements(iat)%radius
       end do
    case('A','a')
       do iat=1,pub_cell%nat
          radii(1,iat) = elements(iat)%radius_cond
       end do
    case('F','f')
       radii(1,:) = radius
    case('L','l')
       do iat=1,pub_cell%nat
          radii(1,iat) = elements(iat)%max_core_wf_radius
       end do
    end select
    select case (mode2)
    case('C','c')
       do iat=1,pub_cell%nat
          radii(2,iat) = elements(iat)%max_core_radius
       end do
    case('R','r')
       do iat=1,pub_cell%nat
          radii(2,iat) = elements(iat)%radius
       end do
    case('A','a')
       do iat=1,pub_cell%nat
          radii(2,iat) = elements(iat)%radius_cond
       end do
    case('F','f')
       radii(2,:) = radius
    case('L','l')
       do iat=1,pub_cell%nat
          radii(2,iat) = elements(iat)%max_core_wf_radius
       end do
    end select

    ! Convert Cartesian coordinates from cell into fractional coordinates in
    ! the interval [0,1)
    do iat=1,pub_cell%nat
       fcoord(1,iat) = elements(iat)%centre .dot. pub_cell%b1
       fcoord(2,iat) = elements(iat)%centre .dot. pub_cell%b2
       fcoord(3,iat) = elements(iat)%centre .dot. pub_cell%b3
    end do
#ifdef ACCELRYS
    do iat=1,pub_cell%nat
       do jat = 1,3
           if (fcoord(jat,iat) < 0.0_DP .and. fcoord(jat,iat) > (-tolerance)) then
              fcoord(jat,iat) = 0.0_DP
           elseif (fcoord(jat,iat) > 1.0_DP .and. fcoord(jat,iat) < 1.0_DP+tolerance) then
              fcoord(jat,iat) = 1.0_DP
           end if
       end do
    end do
#endif

    recip_twopi = 0.5_DP / pi
    fcoord = fcoord * recip_twopi
    fcoord = modulo(fcoord,1.0_DP)

    ! Make local copies of lattice vectors
    a1(1) = pub_cell%a1%x ; a1(2) = pub_cell%a1%y ; a1(3) = pub_cell%a1%z
    a2(1) = pub_cell%a2%x ; a2(2) = pub_cell%a2%y ; a2(3) = pub_cell%a2%z
    a3(1) = pub_cell%a3%x ; a3(2) = pub_cell%a3%y ; a3(3) = pub_cell%a3%z

    ! Figure out the maximum potential cutoff separation - this defines the
    ! minimum subcell size required
    subcell_size = maxval(radii(1,:)) + maxval(radii(2,:))

    ! Avoid subcell sizes which are unphysically small
    subcell_size = max(subcell_size,1.0_DP)

    ! Recall that the magnitude of reciprocal lattice vectors is related
    ! to plane spacings in real-space
    spacing(1) = TWO_PI / sqrt(pub_cell%b1 .dot. pub_cell%b1)
    spacing(2) = TWO_PI / sqrt(pub_cell%b2 .dot. pub_cell%b2)
    spacing(3) = TWO_PI / sqrt(pub_cell%b3 .dot. pub_cell%b3)
    if (.not. present(radius) .and. subcell_size > minval(spacing)) then
       if (pub_on_root) write(stdout,*) 'Error in parallel_strategy_list_overlaps: &
            &spheres exceed cell size.'
       call comms_abort
    end if

    ! Work out how many of our cells will fit into the unit cell
    ! Do not allow more than 1024 subcells in any direction to avoid
    ! integer overflow
    do idim=1,3
       num_subcells(idim) = min(max(int(spacing(idim) / subcell_size),1),1024)
       if (num_subcells(idim) == 2) num_subcells(idim) = 1
    end do
    total_num_subcells = num_subcells(1)*num_subcells(2)*num_subcells(3)

    ! Work out which subcell each atom belongs to - each subcell is labelled
    ! by a number from 0 to total_num_subcells-1
    do iat=1,pub_cell%nat
       do idim=1,3
          iscd(idim) = int(num_subcells(idim) * fcoord(idim,iat))
       end do
       atom_subcell(iat) = iscd(1)+num_subcells(1)*(iscd(2)+ &
            num_subcells(2)*iscd(3))
    end do

    ! Count the number of atoms in each subcell and the maximum
    allocate(num_atoms_subcell(0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell',ierr)
    num_atoms_subcell = 0
    do iat=1,pub_cell%nat
       isc = atom_subcell(iat)
       num_atoms_subcell(isc) = num_atoms_subcell(isc)+1
    end do
    max_atoms_subcell = maxval(num_atoms_subcell)

    ! Make up the cell-list i.e. the list of atoms in each subcell
    allocate(subcell_list(max_atoms_subcell,0:total_num_subcells-1),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'subcell_list',ierr)
    num_atoms_subcell = 0
    do iat=1,pub_cell%nat
       isc = atom_subcell(iat)
       num_atoms_subcell(isc) = num_atoms_subcell(isc)+1
       subcell_list(num_atoms_subcell(isc),isc) = iat
    end do

    ! Count the number of overlaps for each atom
    num_my_overlaps = 0
    ! Loop over all subcells i, distributed across nodes
    do isc=pub_my_node_id,total_num_subcells-1,pub_total_num_nodes
       iscd(1) = mod(isc,num_subcells(1))    ! Subcell i a1 location
       isc23 = isc/num_subcells(1)
       iscd(2) = mod(isc23,num_subcells(2))  ! Subcell i a2 location
       iscd(3) = isc23/num_subcells(2)       ! Subcell i a3 location
       ! Now loop over this and 26 neighbouring subcells j
       subcellj1: do ksc=0,26
          kscd(1) = mod(ksc,3)-1    ! Relative location of subcell j
          ksc23 = ksc/3             ! wrt subcell i (i.e. -1, 0, 1 in
          kscd(2) = mod(ksc23,3)-1  ! each dimension)
          kscd(3) = (ksc23/3)-1
          ! Calculate absolute location of subcell j
          jscd = iscd+kscd
          ! Apply appropriate boundary conditions
          fextra = 0.0_DP
          do idim=1,3 ! for each dimension
             if (jscd(idim) == -1) then
                if (periodic(idim)) then
                   jscd(idim) = num_subcells(idim)-1
                   fextra(idim) = -1.0_DP
                else
                   cycle subcellj1
                end if
             end if
             if (jscd(idim) == num_subcells(idim)) then
                if (periodic(idim)) then
                   jscd(idim) = 0
                   fextra(idim) = 1.0_DP
                else
                   cycle subcellj1
                end if
             end if
          end do
          ! Calculate label for subcell j
          jsc = jscd(1)+num_subcells(1)*(jscd(2)+num_subcells(2)*jscd(3))
          ! Loop over atoms in subcell i
          do iatsc=1,num_atoms_subcell(isc)
             iat = subcell_list(iatsc,isc)
             if (radii(1,iat) == 0.0_DP) cycle ! radius 0 means nothing there
             ! Loop over atoms in subcell j
             do jatsc=1,num_atoms_subcell(jsc)
                jat = subcell_list(jatsc,jsc)
                if (radii(2,jat) == 0.0_DP) cycle ! radius 0 means nothing there
                fdiff(:) = fcoord(:,jat)+fextra(:)-fcoord(:,iat)
                adiff(1) = fdiff(1)*a1(1)+fdiff(2)*a2(1)+fdiff(3)*a3(1)
                adiff(2) = fdiff(1)*a1(2)+fdiff(2)*a2(2)+fdiff(3)*a3(2)
                adiff(3) = fdiff(1)*a1(3)+fdiff(2)*a2(3)+fdiff(3)*a3(3)
                dist_sq = adiff(1)*adiff(1)+adiff(2)*adiff(2)+adiff(3)*adiff(3)
                cutoff = radii(1,iat) + radii(2,jat)
                if (dist_sq <= cutoff*cutoff) then
                   num_my_overlaps(iat) = num_my_overlaps(iat)+1
                end if
             end do ! Loop over atoms in subcell j
          end do ! Loop over atoms in subcell i
       end do subcellj1 ! Loop over subcell j
    end do ! Loop over subcell i

    ! Collect results in from all nodes
    max_my_overlaps = maxval(num_my_overlaps)
    call comms_reduce('MAX',max_my_overlaps)

    ! ndmh: Number of overlaps of one atom will never be greater than pub_cell%nat
    max_my_overlaps = min(max_my_overlaps, pub_cell%nat)

    ! Allocate workspace my_pub_overlap_list and yr_pub_overlap_list
    allocate(my_pub_overlap_list(max_my_overlaps,pub_cell%nat),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'my_pub_overlap_list',ierr)
    allocate(yr_pub_overlap_list(max_my_overlaps,pub_cell%nat),stat=ierr)
    call utils_alloc_check('parallel_strategy_list_overlaps', &
         'yr_pub_overlap_list',ierr)

    ! Make up the overlap list for this node, this time counting multiple
    ! overlaps only once
    num_my_overlaps = 0
    overlapped = .false.
    ! Loop over all subcells i, distributed across nodes
    do isc=pub_my_node_id,total_num_subcells-1,pub_total_num_nodes
       iscd(1) = mod(isc,num_subcells(1))    ! Subcell i a1 location
       isc23 = isc/num_subcells(1)
       iscd(2) = mod(isc23,num_subcells(2))  ! Subcell i a2 location
       iscd(3) = isc23/num_subcells(2)       ! Subcell i a3 location
       ! Loop over atoms in subcell i
       do iatsc=1,num_atoms_subcell(isc)
          iat = subcell_list(iatsc,isc)
          if (radii(1,iat) == 0.0_DP) cycle ! radius 0 means nothing there
          ! Now loop over this subcell and 26 neighbouring subcells j
          subcellj2: do ksc=0,26
             kscd(1) = mod(ksc,3)-1    ! Relative location of subcell j
             ksc23 = ksc/3             ! wrt subcell i (i.e. -1, 0, 1 in
             kscd(2) = mod(ksc23,3)-1  ! each dimension)
             kscd(3) = (ksc23/3)-1
             ! Calculate absolute location of subcell j
             jscd = iscd+kscd
             ! Apply appropriate boundary conditions
             fextra = 0.0_DP
             do idim=1,3 ! for each dimension
                if (jscd(idim) == -1) then
                   if (periodic(idim)) then
                      jscd(idim) = num_subcells(idim)-1
                      fextra(idim) = -1.0_DP
                   else
                      cycle subcellj2
                   end if
                end if
                if (jscd(idim) == num_subcells(idim)) then
                   if (periodic(idim)) then
                      jscd(idim) = 0
                      fextra(idim) = 1.0_DP
                   else
                      cycle subcellj2
                   end if
                end if
             end do
             ! Calculate label for subcell j
             jsc = jscd(1)+num_subcells(1)*(jscd(2)+num_subcells(2)*jscd(3))
             ! Loop over atoms in subcell j
             do jatsc=1,num_atoms_subcell(jsc)
                jat = subcell_list(jatsc,jsc)
                if (radii(2,jat) == 0.0_DP) cycle ! radius 0 means nothing there
                fdiff(:) = fcoord(:,jat)+fextra(:)-fcoord(:,iat)
                adiff(1) = fdiff(1)*a1(1)+fdiff(2)*a2(1)+fdiff(3)*a3(1)
                adiff(2) = fdiff(1)*a1(2)+fdiff(2)*a2(2)+fdiff(3)*a3(2)
                adiff(3) = fdiff(1)*a1(3)+fdiff(2)*a2(3)+fdiff(3)*a3(3)
                dist_sq = adiff(1)*adiff(1)+adiff(2)*adiff(2)+adiff(3)*adiff(3)
                cutoff = radii(1,iat) + radii(2,jat)
                if (dist_sq <= cutoff*cutoff) then
                   if (.not. overlapped(jat)) then
                      num_my_overlaps(iat) = num_my_overlaps(iat)+1
                      my_pub_overlap_list(num_my_overlaps(iat),iat) = jat
                      overlapped(jat) = .true.
                   end if
                end if
             end do ! Loop over atoms in subcell j
          end do subcellj2 ! Loop over subcell j
          ! Restore set of overlap flags
          do iovlap=1,num_my_overlaps(iat)
             overlapped(my_pub_overlap_list(iovlap,iat)) = .false.
          end do
       end do ! Loop over atoms in subcell i
    end do ! Loop over subcell i

    ! Collect results in from all nodes
    max_my_overlaps = maxval(num_my_overlaps)
    call comms_reduce('MAX',max_my_overlaps)

    pub_num_overlaps = num_my_overlaps
    call comms_reduce('SUM',pub_num_overlaps)
    pub_max_overlaps = maxval(pub_num_overlaps)

!CW
    pub_max_overlaps=max(pub_max_overlaps,pub_cell%nat_hub)
!END CW

    ! Allocate module array pub_overlap_list
    ierr = 0
    if (allocated(pub_overlap_list)) then
       if (size(pub_overlap_list,1) /= pub_max_overlaps .or. &
            size(pub_overlap_list,2) /= pub_cell%nat) then
          deallocate(pub_overlap_list,stat=ierr)
          call utils_dealloc_check('parallel_strategy_list_overlaps', &
               'pub_overlap_list',ierr)
          allocate(pub_overlap_list(pub_max_overlaps,pub_cell%nat),stat=ierr)
          call utils_alloc_check('parallel_strategy_list_overlaps', &
               'pub_overlap_list',ierr)
       end if
    else
       allocate(pub_overlap_list(pub_max_overlaps,pub_cell%nat),stat=ierr)
       call utils_alloc_check('parallel_strategy_list_overlaps', &
            'pub_overlap_list',ierr)
    end if


    ! Make up global overlap list, looping over nodes to gather up information
    pub_num_overlaps = 0
    pub_overlap_list = 0
    do node=0,pub_total_num_nodes-1
       if (pub_my_node_id == node) then
          num_yr_overlaps = num_my_overlaps
          yr_pub_overlap_list = my_pub_overlap_list
       end if
       call comms_bcast(node,num_yr_overlaps)
       call comms_bcast(node,yr_pub_overlap_list)
       do iat=1,pub_cell%nat
          nyrat = num_yr_overlaps(iat)
          novat = pub_num_overlaps(iat)
          pub_overlap_list(novat+1:novat+nyrat,iat) = &
               yr_pub_overlap_list(1:nyrat,iat)
          pub_num_overlaps(iat) = novat+nyrat
       end do
    end do

#ifdef DEBUG
    ! Verify final result
    call internal_verify
#endif

    ! Deallocate work arrays
    deallocate(yr_pub_overlap_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'yr_pub_overlap_list',ierr)
    deallocate(my_pub_overlap_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'my_pub_overlap_list',ierr)
    deallocate(subcell_list,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'subcell_list',ierr)
    deallocate(num_atoms_subcell,stat=ierr)
    call utils_dealloc_check('parallel_strategy_list_overlaps', &
         'num_atoms_subcell',ierr)

    call internal_deallocate_1

  contains

    !--------------------------------------------------------------------------
#ifdef DEBUG
    subroutine internal_verify

      !=======================================================================!
      ! This subroutine verifies the overlap list generated by the parent     !
      ! subroutine for consistency                                            !
      !-----------------------------------------------------------------------!
      ! Arguments:							      !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 16/03/04                                     !
      !=======================================================================!

      implicit none

      ! Local variables
      integer :: iat,jat   ! Atom counters
      integer :: novlaps   ! Number of overlaps
      integer :: iovlap    ! Overlap counter
      integer :: jovlap    ! Overlap counter
      logical :: found     ! Found flag

      ! Loop over all atoms in cell iat
      do iat=1,pub_cell%nat

         ! Check number of overlaps
         novlaps = pub_num_overlaps(iat)
         if (novlaps < 1 .and. radii(1,iat) > 0.0_DP .and. &
              radii(2,iat) > 0.0_DP) then
            if (pub_on_root) write(stdout,'(a/a,i6)') &
                 'Error verifying overlap list (parallel_strategy_mod.F90):', &
                 '  no overlaps detected for atom ',iat
            call comms_abort
         end if
         if (novlaps > pub_cell%nat) then
            if (pub_on_root) write(stdout,'(a/a,i6)') &
                 'Error verifying overlap list (parallel_strategy_mod.F90):', &
                 '  too many overlaps detected for atom ',iat
            call comms_abort
         end if

         ! Loop over overlapping atoms jat
         do iovlap=1,novlaps
            jat = pub_overlap_list(iovlap,iat)

            ! Check that this atom only appears once in the overlap list
            do jovlap=iovlap+1,novlaps
               if (pub_overlap_list(jovlap,iat) == jat) then
                  if (pub_on_root) write(stdout,'(a/2(a,i6))') &
                       'Error verifying overlap list &
                       &(parallel_strategy_mod.F90):', &
                       '  repeated overlap for atom ',iat,' with atom ',jat
                  call comms_abort
               end if
            end do

            ! Check that cross-reference exists (if appropriate)
            if (mode1 == mode2) then

               ! Loop over atoms which overlap atom jat
               found = .false.
               do jovlap=1,pub_num_overlaps(jat)
                  if (pub_overlap_list(jovlap,jat) == iat) then
                     found = .true.
                     exit
                  end if
               end do
               if (.not. found) then
                  if (pub_on_root) write(stdout,'(a/2(a,i6))') &
                       'Error verifying overlap list &
                       &(parallel_strategy_mod.F90):', &
                       '  no cross-reference overlap for atom ',iat, &
                       ' with atom ',jat
                  call comms_abort
               end if

            end if

         end do  ! jat

      end do  ! iat

    end subroutine internal_verify
#endif
    !--------------------------------------------------------------------------

    subroutine internal_check_args(have_radius)

      !=======================================================================!
      ! This subroutine checks the arguments to the parent subroutine         !
      !-----------------------------------------------------------------------!
      ! Arguments:							      !
      ! have_radius (input) : flag indicating presence of optional argument   !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      implicit none

      ! Arguments
      logical, intent(in) :: have_radius

      if (pub_cell%nat <= 0) then
         if (pub_on_root) write(stdout,*) &
              'Error in parallel_strategy_list_overlaps: no atoms in cell!'
         call comms_abort
      end if
      if (size(elements) /= pub_cell%nat) then
         if (pub_on_root) write(stdout,*) &
              'Error in parallel_strategy_list_overlaps: incompatible nat'
         call comms_abort
      end if
      if (mode1 /= 'C' .and. mode1 /= 'c' .and. mode1 /= 'R' .and. &
           mode1 /= 'r' .and. mode1 /= 'F' .and. mode1 /= 'f' .and. &
           mode1 /= 'A' .and. mode1 /= 'a' .and. mode1 /= 'L' .and. &
           mode1 /= 'l') then
         if (pub_on_root) write(stdout,'(4a)') &
              'Error in parallel_strategy_list_overlaps: unknown mode1 "', &
              mode1,'"'
         call comms_abort
      end if
      if (mode2 /= 'C' .and. mode2 /= 'c' .and. mode2 /= 'R' .and. &
           mode2 /= 'r' .and. mode2 /= 'F' .and. mode2 /= 'f' .and. &
           mode2 /= 'A' .and. mode2 /= 'a' .and. mode2 /= 'L' .and. &
           mode2 /= 'l') then
         if (pub_on_root) write(stdout,'(4a)') &
              'Error in parallel_strategy_list_overlaps: unknown mode2 "', &
              mode2,'"'
         call comms_abort
      end if
      if ((mode1 == 'F' .or. mode1 == 'f' .or. mode2 == 'F' .or. &
           mode2 == 'f') .and. (.not. have_radius)) then
         if (pub_on_root) write(stdout,'(a)') &
              'Error in parallel_strategy_list_overlaps: &
              &radius must be present for fixed mode'
         call comms_abort
      end if

    end subroutine internal_check_args

    !--------------------------------------------------------------------------

    subroutine internal_allocate_1

      !=======================================================================!
      ! This subroutine (de)allocates arrays as required by the parent        !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:							      !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Written by Peter Haynes, 20/11/03                                     !
      !=======================================================================!

      implicit none

      ! Allocate module array pub_num_overlaps
      ierr = 0
      if (allocated(pub_num_overlaps)) then
         if (size(pub_num_overlaps) /= pub_cell%nat) then
            deallocate(pub_num_overlaps,stat=ierr)
            call utils_dealloc_check('internal_allocate_1 &
                 &(parallel_strategy_list_overlaps)','pub_num_overlaps',ierr)
            allocate(pub_num_overlaps(pub_cell%nat),stat=ierr)
            call utils_alloc_check('internal_allocate_1 &
                 &(parallel_strategy_list_overlaps)','pub_num_overlaps',ierr)
         end if
      else
         allocate(pub_num_overlaps(pub_cell%nat),stat=ierr)
         call utils_alloc_check('internal_allocate_1 &
              &(parallel_strategy_list_overlaps)','pub_num_overlaps',ierr)
      end if


      ! Allocate work arrays
      allocate(radii(2,pub_cell%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','radii',ierr)
      allocate(atom_subcell(pub_cell%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','atom_subcell',ierr)
      allocate(fcoord(3,pub_cell%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','fcoord',ierr)
      allocate(num_my_overlaps(pub_cell%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','num_my_overlaps',ierr)
      allocate(num_yr_overlaps(pub_cell%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','num_yr_overlaps',ierr)
      allocate(overlapped(pub_cell%nat),stat=ierr)
      call utils_alloc_check('internal_allocate_1 &
           &(parallel_strategy_list_overlaps)','overlapped',ierr)

   end subroutine internal_allocate_1

    subroutine internal_deallocate_1

      !=======================================================================!
      ! This subroutine deallocates arrays as required by the parent          !
      ! subroutine.                                                           !
      !-----------------------------------------------------------------------!
      ! Arguments:							      !
      ! None                                                                  !
      !-----------------------------------------------------------------------!
      ! Added by Nicholas Hine, 18/05/09 to ensure utils_dealloc_check marks  !
      ! the right arrays as having been deallocated when doing array checking !
      !=======================================================================!

      implicit none

      ! Allocate internal arrays
      deallocate(overlapped,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'overlapped',ierr)
      deallocate(num_yr_overlaps,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'num_yr_overlaps',ierr)
      deallocate(num_my_overlaps,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'num_my_overlaps',ierr)
      deallocate(fcoord,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'fcoord',ierr)
      deallocate(atom_subcell,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'atom_subcell',ierr)
      deallocate(radii,stat=ierr)
      call utils_dealloc_check('internal_deallocate_1 &
           &(parallel_strategy_list_overlaps)', 'radii',ierr)

   end subroutine internal_deallocate_1

  end subroutine parallel_strategy_list_overlaps

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module parallel_strategy
