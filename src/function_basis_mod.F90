! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by                              !
!                                                                             !
!   Nicholas D.M. Hine                                                        !
!                                                                             !
!   in July 2009.                                                             !
!                                                                             !
!   Based on in part on previous code by                                      !
!                                                                             !
!   Chris-Kriton Skylaris, Arash A. Mostofi, Peter D. Haynes                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module function_basis

   use basis, only: FUNCTION_TIGHT_BOX, SPHERE
   use constants, only: DP, stdout


  implicit none

  private

  ! ndmh: This structure contains information describing the size, spheres,
  ! ndmh: tightboxes and parallel distribution of a given set of functions,
  ! ndmh: such as: the set of NGWFs, the nonlocal pseudopotential
  ! ndmh: projectors, the PAW partial waves, a DMFT correlated subspace, ...
  type FUNC_BASIS

     ! ndmh: the total number of functions in this set
     integer :: num

     ! ndmh: the number of functions on this node
     integer :: node_num

     ! ndmh: the number of functions on each node (0:pub_total_num_nodes-1)
     integer, allocatable :: num_on_node(:)

     ! ndmh: the first function on each node (0:pub_total_num_nodes)
     integer, allocatable :: first_on_node(:)

     ! ndmh: the number of functions on each atom (1:pub_cell%nat)
     integer, allocatable :: num_on_atom(:)

     ! ndmh: the first function on each atom (1:pub_cell%nat)
     integer, allocatable :: first_on_atom(:)

     ! ndmh: the node of each function (1:num)
     integer, allocatable :: node_of_func(:)

     ! ndmh: the atom of each function (1:num)
     integer, allocatable :: atom_of_func(:)

     ! ndmh: the maximum number of functions on any node
     integer :: max_on_node

     ! ndmh: the maximum number of functions on any atom
     integer :: max_on_atom

     ! ndmh: the number of ppds describing the functions on this node
     integer :: n_ppds

     ! ndmh: the number of ppds in all spheres over all nodes (1:num)
     integer, allocatable :: n_ppds_sphere(:)

     ! ndmh: the maximum number of ppds in any sphere of these functions
     integer :: max_n_ppds_sphere

     ! ndmh: required size of receive buffer for functions from other nodes
     integer :: func_on_grid_buffer_size

     ! ndmh: maximum sizes of tightboxes (if allocated)
     integer :: maxtight_pts1
     integer :: maxtight_pts2
     integer :: maxtight_pts3

     ! ndmh: the spheres describing the functions on this node (1:node_num)
     type(SPHERE), allocatable :: spheres(:)

     ! ndmh: the tightboxes describing the functions on this node (1:node_num)
     ! ndmh: NOT ALWAYS ALLOCATED FOR ALL FUNCTION TYPES
     type(FUNCTION_TIGHT_BOX), allocatable :: tight_boxes(:)

     ! ndmh: the tightboxes describing all the functions (1:num)
     ! ndmh: NOT ALWAYS ALLOCATED FOR ALL FUNCTION TYPES
     type(FUNCTION_TIGHT_BOX), allocatable :: all_tbs(:)

     ! ndmh: string identifying this function basis (eg 'ngwfs')
     character(len=20) :: name

  end type FUNC_BASIS

  public :: FUNC_BASIS

  ! cks
  type(SPHERE), public :: pub_buffer_sphere

   ! ndmh: special function request values and tags
   integer, parameter :: FUNCS_WAITING = -1000
   integer, parameter :: FUNCS_DONE = -2000
   integer, parameter :: req_tag = 10000000
   integer, parameter :: req_buffer_size = 40
   integer, parameter :: probe_frequency = 4
   integer, parameter :: num_buffers = 1

  ! ndmh: buffer to store requested indices
  integer :: req_buffer(req_buffer_size)
  integer :: req_buffer_index
  integer :: probe_count

  ! ndmh: public subroutines

  ! ndmh: init/exit routines
  public :: function_basis_allocate
  public :: function_basis_distribute
  public :: function_basis_deallocate
  public :: function_basis_init_spheres
  public :: function_basis_copy_spheres
  public :: function_basis_init_tight_boxes
  public :: function_basis_gath_all_tbs
  public :: function_basis_init_uni_tb
  public :: function_basis_exit_uni_tb
  public :: function_basis_est_num_psincs
  public :: function_basis_estimate_size

  ! ndmh: row sums / integrals routines
  public :: function_basis_batch_row_plan
  public :: function_basis_sum_fftbox_batch
  public :: function_basis_sum_ppd_funcs

  ! ndmh: ppd_to_tightbox (parallelised)
  public :: function_basis_ppds_to_tightbox
  public :: function_basis_tightbox_to_ppds
  public :: function_basis_ppds_to_sph_waves
  public :: function_basis_sph_waves_to_ppds

  ! ndmh: function communications routines
  public :: function_basis_init_requests
  public :: function_basis_request
  public :: function_basis_await_requests
  public :: function_basis_respond_to_reqs
  public :: function_basis_send
  public :: function_basis_recv

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_allocate(fbasis,num,name)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the function basis   !
    ! type whose size does not depend on parallel strategy determinations     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fbasis (inout)            : function basis type describing these funcs  !
    ! num (input)               : total number of functions                   !
    ! name (input)              : identifying string for this function basis  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 10/07/2009                                  !
    !=========================================================================!

    use comms, only: pub_total_num_nodes
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    integer, intent(in) :: num
    character(len=*) :: name

    ! Local Variables
    integer :: ierr

    fbasis%num = num
    fbasis%name = name

    allocate(fbasis%num_on_node(0:pub_total_num_nodes-1), stat=ierr)
    call utils_alloc_check('function_basis_allocate_'//fbasis%name, &
         'fbasis%num_on_node',ierr)

    allocate(fbasis%first_on_node(0:pub_total_num_nodes),stat=ierr)
    call utils_alloc_check('function_basis_allocate_'//fbasis%name, &
         'fbasis%first_on_node',ierr)

    allocate(fbasis%num_on_atom(1:pub_cell%nat),stat=ierr)
    call utils_alloc_check('function_basis_allocate_'//fbasis%name, &
         'fbasis%num_on_atom',ierr)

    allocate(fbasis%first_on_atom(1:pub_cell%nat),stat=ierr)
    call utils_alloc_check('function_basis_allocate_'//fbasis%name, &
         'fbasis%first_on_atom',ierr)

    allocate(fbasis%node_of_func(1:num),stat=ierr)
    call utils_alloc_check('function_basis_allocate_'//fbasis%name, &
         'fbasis%node_of_func',ierr)

    allocate(fbasis%atom_of_func(1:num),stat=ierr)
    call utils_alloc_check('function_basis_allocate_'//fbasis%name, &
         'fbasis%atom_of_func',ierr)

  end subroutine function_basis_allocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_distribute(fbasis,elements)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the function basis   !
    ! type whose size does not depend on parallel strategy determinations     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fbasis (inout)            : function basis type describing these funcs  !
    ! num (input)               : total number of functions                   !
    ! name (input)              : identifying string for this function basis  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 10/07/2009                                  !
    !=========================================================================!

    use comms, only: comms_abort, comms_reduce, pub_my_node_id
    use constants, only: stdout
    use hubbard_init, only: h_species
    use ion, only: ELEMENT
    use parallel_strategy, only: parallel_strategy_distr_funcs, pub_orig_atom, &
         pub_elements_on_node
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Local Variables
    integer :: ierr, iat, hsp
    integer, allocatable :: nfuncs(:) ! Number of funcs per atom (input order)

    allocate(nfuncs(pub_cell%nat),stat=ierr)
    call utils_alloc_check('function_basis_distribute_'//fbasis%name, &
         'nfuncs',ierr)

    ! ndmh: find number of functions on each atom
    do iat=1,pub_cell%nat
       if (fbasis%name=='ngwfs'.or.fbasis%name=='tmp_ngwfs') then
          nfuncs(iat) = elements(iat)%nfunctions
       else if (fbasis%name=='projs') then
          nfuncs(iat) = elements(iat)%nprojectors
       else if (fbasis%name=='pawpws') then
          nfuncs(iat) = elements(iat)%npawpws
       else if (fbasis%name=='hub_projs') then
          nfuncs(iat) = 0
          do hsp=1,pub_cell%num_hub_species
             if (h_species(hsp)%hub_species==elements(iat)%species_id) then
                nfuncs(iat) = 2 * h_species(hsp)%hub_ang_mom + 1
             end if
          end do
       else if (fbasis%name=='ngwfs_cond'.or.fbasis%name=='tmp_ngwfs_cond') then
          nfuncs(iat) = elements(iat)%nfunctions_cond
       else if (fbasis%name=='ngwfs_joint') then
          nfuncs(iat) = elements(iat)%nfunctions + elements(iat)%nfunctions_cond
       else if (fbasis%name=='ngwfs_aux') then
          nfuncs(iat) = elements(iat)%nfunctions_aux
       else if (fbasis%name=='corewfs') then
          nfuncs(iat) = elements(iat)%ncorewfs
       else
          write(stdout,'(2a)') 'Error in function_basis_distribute: &
               &unrecognised function basis identifier: ',fbasis%name
          call comms_abort
       end if
    end do

    call parallel_strategy_distr_funcs(fbasis%num, nfuncs, pub_orig_atom, &
         fbasis%first_on_node, fbasis%num_on_node, &
         fbasis%first_on_atom, fbasis%num_on_atom, &
         fbasis%node_of_func, fbasis%atom_of_func, fbasis%max_on_node, &
         fbasis%max_on_atom)

    fbasis%node_num = fbasis%num_on_node(pub_my_node_id)

    deallocate(nfuncs,stat=ierr)
    call utils_dealloc_check('function_basis_distribute_'//fbasis%name, &
         'nfuncs',ierr)

  end subroutine function_basis_distribute


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_deallocate(fbasis)

    !=========================================================================!
    ! This subroutine deallocates memory for the arrays in the function basis !
    ! type whose size does not depend on parallel strategy determinations     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! fbasis (inout)            : function basis type describing these funcs  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 10/07/2009                                  !
    !=========================================================================!

    use basis, only: basis_sphere_deallocate
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis

    ! Local Variables
    integer :: ierr, ifunc

    if (allocated(fbasis%all_tbs)) then
       deallocate(fbasis%all_tbs,stat=ierr)
       call utils_dealloc_check('function_basis_gath_all_tbs_'//fbasis%name, &
            'fbasis%all_tbs',ierr)
    end if

    if (allocated(fbasis%tight_boxes)) then
       deallocate(fbasis%tight_boxes,stat=ierr)
       call utils_dealloc_check('function_basis_init_tight_boxes_'//&
            fbasis%name,'fbasis%tight_boxes',ierr)
    end if

    if (fbasis%name=='ngwfs') then
       call basis_sphere_deallocate(pub_buffer_sphere)
    end if

    if (allocated(fbasis%n_ppds_sphere)) then
       deallocate(fbasis%n_ppds_sphere,stat=ierr)
       call utils_dealloc_check('function_basis_init_spheres_'//fbasis%name, &
            'fbasis%n_ppds_sphere',ierr)
    end if

    if (allocated(fbasis%spheres)) then
       do ifunc=1,fbasis%node_num
          call basis_sphere_deallocate(fbasis%spheres(ifunc))
       end do
       deallocate(fbasis%spheres,stat=ierr)
       call utils_dealloc_check('function_basis_init_spheres_'//fbasis%name, &
            'fbasis%spheres',ierr)
    end if

    deallocate(fbasis%atom_of_func,stat=ierr)
    call utils_dealloc_check('function_basis_deallocate_'//fbasis%name, &
         'fbasis%atom_of_func',ierr)

    deallocate(fbasis%node_of_func,stat=ierr)
    call utils_dealloc_check('function_basis_deallocate_'//fbasis%name, &
         'fbasis%node_of_func',ierr)

    deallocate(fbasis%first_on_atom,stat=ierr)
    call utils_dealloc_check('function_basis_deallocate_'//fbasis%name, &
         'fbasis%first_on_atom',ierr)

    deallocate(fbasis%num_on_atom,stat=ierr)
    call utils_dealloc_check('function_basis_deallocate_'//fbasis%name, &
         'fbasis%num_on_atom',ierr)

    deallocate(fbasis%first_on_node,stat=ierr)
    call utils_dealloc_check('function_basis_deallocate_'//fbasis%name, &
         'fbasis%first_on_node',ierr)

    deallocate(fbasis%num_on_node, stat=ierr)
    call utils_dealloc_check('function_basis_deallocate_'//fbasis%name, &
         'fbasis%num_on_node',ierr)

  end subroutine function_basis_deallocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_copy_spheres(fbasis,fbasis_src1,fbasis_src2)

    !========================================================================!
    ! This subroutine copies the spheres array describing the functions of a !
    ! function basis to the spheres array of another function basis.         !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! fbasis (inout)            : Function basis to which the spheres are    !
    !                             to be copied.                              !
    ! fbasis_src (inout)        : Function basis from which the spheres are  !
    !                             to be copied.                              !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/07/10.                                  !
    !========================================================================!

    use basis, only: basis_copy_sphere
    use comms, only: comms_bcast, pub_my_node_id, pub_total_num_nodes
    use ion, only: element
    use parallel_strategy, only: pub_first_atom_on_node, pub_num_atoms_on_node
    use simulation_cell, only : pub_cell
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(FUNC_BASIS), intent(in) :: fbasis_src1
    type(FUNC_BASIS), intent(in), optional :: fbasis_src2

    ! Local Variables
    integer :: ifunc_src1, ifunc_src2, ifunc
    integer :: global_ifunc, node
    integer :: count, count_src1, count_src2
    integer :: current_offset
    integer :: iat, loc_iat
    integer :: ierr
    character(len=10) :: iat_str

    allocate(fbasis%spheres(fbasis%node_num),stat=ierr)
    call utils_alloc_check('function_basis_init_spheres_'//fbasis%name, &
         'fbasis%spheres',ierr)

    allocate(fbasis%n_ppds_sphere(1:fbasis%num),stat=ierr)
    call utils_alloc_check('function_basis_init_spheres_'//fbasis%name, &
         'fbasis%n_ppds_sphere',ierr)

    ! ndmh: loop over all atoms on node and copy the spheres of each one
    count = 0
    count_src1 = 0
    count_src2 = 0
    current_offset = 1
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1

       if (.not.present(fbasis_src2)) then
          if (fbasis_src1%num_on_atom(iat)/=fbasis%num_on_atom(iat)) then
             write(iat_str,'(i10)') iat
             iat_str = adjustl(trim(iat_str))
             call utils_abort('Error in function_basis_init_spheres: &
                  &Mismatching sphere counts on atom '//iat_str)
          end if
       else
          if ((fbasis_src1%num_on_atom(iat)+fbasis_src2%num_on_atom(iat)) &
               /=fbasis%num_on_atom(iat)) then
             write(iat_str,'(i10)') iat
             iat_str = adjustl(trim(iat_str))
             call utils_abort('Error in function_basis_init_spheres: &
                  &Mismatching sphere counts on atom '//iat_str)
          end if
       end if

       do ifunc_src1=1,fbasis_src1%num_on_atom(iat)
          count_src1 = count_src1 + 1
          count = count + 1

          ! ndmh: copy the sphere from one basis to the other
          call basis_copy_sphere(fbasis%spheres(count),&
               fbasis_src1%spheres(count_src1),current_offset)
          current_offset = current_offset + &
               pub_cell%n_pts*fbasis%spheres(count)%n_ppds_sphere
       end do

       if (present(fbasis_src2)) then
          do ifunc_src2=1,fbasis_src2%num_on_atom(iat)
             count_src2 = count_src2 + 1
             count = count + 1

             ! ndmh: copy the sphere from one basis to the other
             call basis_copy_sphere(fbasis%spheres(count),&
                  fbasis_src2%spheres(count_src2),current_offset)
             current_offset = current_offset + &
                  pub_cell%n_pts*fbasis%spheres(count)%n_ppds_sphere
          end do
       end if

    end do

    ! ndmh: determine the total number of ppds belonging to the spheres on this
    ! ndmh: node
    fbasis%n_ppds = 0
    do ifunc=1,fbasis%num_on_node(pub_my_node_id)
       fbasis%n_ppds = fbasis%n_ppds + fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    ! ndmh: fill this node's entries in global n_ppds_sphere list
    do ifunc=1,fbasis%node_num
       global_ifunc = ifunc + fbasis%first_on_node(pub_my_node_id) - 1
       fbasis%n_ppds_sphere(global_ifunc) = fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    ! ndmh: collate numbers of ppds in all spheres from all nodes
    do node=0,pub_total_num_nodes-1
       if (fbasis%num_on_node(node) > 0) then
          call comms_bcast(node,fbasis%n_ppds_sphere( &
               fbasis%first_on_node(node)),fbasis%num_on_node(node))
       end if
    end do

    fbasis%max_n_ppds_sphere = maxval(fbasis%n_ppds_sphere)
    fbasis%func_on_grid_buffer_size = fbasis%max_n_ppds_sphere * &
         pub_cell%n_pts * num_buffers

  end subroutine function_basis_copy_spheres


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_init_spheres(fbasis,elements_on_node,extra_radius)

    !========================================================================!
    ! This subroutine initialises the spheres describing the functions on    !
    ! atoms on this node. It also returns the number of ppds needed for the  !
    ! storage of the functions of the current node.                          !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! fbasis (inout)            : Function basis for which the spheres are   !
    !                             to be initialised.                         !
    ! elements_on_node (input)  : Array of element structures of my_node_id. !
    !------------------------------------------------------------------------!
    ! Written as basis_init_ngwf_spheres by Chris-Kriton Skylaris in 2000.   !
    ! Modified by Chris-Kriton Skylaris on 26/8/2003 so that it works with   !
    ! the parallel version (ONETEP).                                         !
    ! Modified by Nicholas Hine on 15/09/2008 to copy initialisation of      !
    ! previous sphere when initialising many projectors per atom with same   !
    ! radius                                                                 !
    ! Generalised for function_basis module by Nicholas Hine on 10/07/2009   !
    !========================================================================!

    use basis, only: basis_copy_sphere, basis_initialise_sphere
    use comms, only: comms_abort, comms_bcast, pub_my_node_id, &
         pub_total_num_nodes
    use constants, only: stdout
    use ion, only: element
    use rundat, only: pub_cond_calculate
    use parallel_strategy, only: pub_first_atom_on_node, pub_num_atoms_on_node
    use simulation_cell, only : pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(ELEMENT), intent(in) :: elements_on_node( &
         pub_num_atoms_on_node(pub_my_node_id))
    real(kind=DP), intent(in), optional :: extra_radius

    ! Local Variables
    integer, allocatable :: ppd_list(:) ! temporary ppd index list
    integer, allocatable :: ppd_loc(:)  ! temporary ppd loc list
    real(kind=DP) :: radius
    integer :: n_ppds_buffer_sphere
    integer :: ifunc, global_ifunc, node
    integer :: count, current_offset, iat, loc_iat
    integer :: ierr

    allocate(fbasis%spheres(fbasis%node_num),stat=ierr)
    call utils_alloc_check('function_basis_init_spheres_'//fbasis%name, &
         'fbasis%spheres',ierr)

    allocate(fbasis%n_ppds_sphere(1:fbasis%num),stat=ierr)
    call utils_alloc_check('function_basis_init_spheres_'//fbasis%name, &
         'fbasis%n_ppds_sphere',ierr)

    ! ndmh: allocate temporary ppd location and number arrays
    allocate(ppd_loc(pub_cell%n_ppds),stat=ierr)
    call utils_alloc_check('function_basis_init_spheres','ppd_loc',ierr)
    allocate(ppd_list(pub_cell%n_ppds),stat=ierr)
    call utils_alloc_check('function_basis_init_spheres','ppd_list',ierr)

    ! ndmh: loop over all atoms and initialise the spheres of each one
    count = 0
    current_offset = 1
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1

       if (fbasis%name=='ngwfs'.or.fbasis%name=='tmp_ngwfs') then
          radius = elements_on_node(loc_iat)%radius
       else if (fbasis%name=='projs') then
          radius = elements_on_node(loc_iat)%max_core_radius
       else if (fbasis%name=='pawpws') then
          radius = elements_on_node(loc_iat)%max_core_radius
       else if (fbasis%name=='corewfs') then !lr408
          radius = elements_on_node(loc_iat)%max_core_wf_radius
       else if (fbasis%name=='hub_projs') then
          radius = elements_on_node(loc_iat)%radius
       else if (fbasis%name=='ngwfs_cond'.or.fbasis%name=='tmp_ngwfs_cond') then
          radius = elements_on_node(loc_iat)%radius_cond
       else if (fbasis%name=='ngwfs_aux') then
          radius = elements_on_node(loc_iat)%radius
       else
          call utils_abort('Error in function_basis_init_spheres: &
               &unrecognised function basis identifier: '//fbasis%name)
       end if

       if (present(extra_radius)) radius = radius + extra_radius

       do ifunc=1,fbasis%num_on_atom(iat)

          count = count + 1

          ! ndmh: copy the sphere if we are initialising more than 1 per atom
          if (ifunc > 1) then
             call basis_copy_sphere(fbasis%spheres(count),&
                  fbasis%spheres(count-1),current_offset)
          else
             call basis_initialise_sphere(fbasis%spheres(count), &
                  elements_on_node(loc_iat)%centre, radius, current_offset, &
                  ppd_list, ppd_loc)
          end if

          current_offset = current_offset + &
               pub_cell%n_pts*fbasis%spheres(count)%n_ppds_sphere

          if (present(extra_radius)) fbasis%spheres(count)%radius = &
               fbasis%spheres(count)%radius - extra_radius

       end do

    end do

    ! ndmh: deallocate temporary ppd location and number arrays
    deallocate(ppd_loc,stat=ierr)
    call utils_dealloc_check('function_basis_init_spheres','ppd_loc',ierr)
    deallocate(ppd_list,stat=ierr)
    call utils_dealloc_check('function_basis_init_spheres','ppd_list',ierr)

    ! ndmh: determine the total number of ppds belonging to the spheres on this
    ! ndmh: node
    fbasis%n_ppds = 0
    do ifunc=1,fbasis%num_on_node(pub_my_node_id)
       fbasis%n_ppds = fbasis%n_ppds + fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    ! ndmh: fill this node's entries in global n_ppds_sphere list
    do ifunc=1,fbasis%node_num
       global_ifunc = ifunc + fbasis%first_on_node(pub_my_node_id) - 1
       fbasis%n_ppds_sphere(global_ifunc) = fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    ! ndmh: collate numbers of ppds in all spheres from all nodes
    do node=0,pub_total_num_nodes-1
       if (fbasis%num_on_node(node) > 0) then
          call comms_bcast(node,fbasis%n_ppds_sphere( &
               fbasis%first_on_node(node)),fbasis%num_on_node(node))
       end if
    end do

    fbasis%max_n_ppds_sphere = maxval(fbasis%n_ppds_sphere)
    fbasis%func_on_grid_buffer_size = fbasis%max_n_ppds_sphere * &
         pub_cell%n_pts * num_buffers

    ! ndmh: allocate pub_buffer_sphere
    pub_buffer_sphere%centre%x = 0.0_DP
    pub_buffer_sphere%centre%y = 0.0_DP
    pub_buffer_sphere%centre%z = 0.0_DP
    pub_buffer_sphere%radius = 0.0_DP
    pub_buffer_sphere%offset = 1
    n_ppds_buffer_sphere = 0

    ! ndmh: deallocate the buffer sphere ppd_list if already allocated
    if (allocated(pub_buffer_sphere%ppd_list)) then
       n_ppds_buffer_sphere = size(pub_buffer_sphere%ppd_list,2)
       deallocate(pub_buffer_sphere%ppd_list,stat=ierr)
       call utils_dealloc_check('function_basis_init_spheres',&
            'current_sphere%ppd_list',ierr)
    end if

    ! ndmh: take whichever is bigger of old value (if there was one) and new one
    n_ppds_buffer_sphere = max(n_ppds_buffer_sphere,fbasis%max_n_ppds_sphere)

    ! ndmh: allocate the buffer sphere to this size
    pub_buffer_sphere%n_ppds_sphere = n_ppds_buffer_sphere
    allocate(pub_buffer_sphere%ppd_list(2,pub_buffer_sphere%n_ppds_sphere), &
         stat=ierr)
    call utils_alloc_check('function_basis_init_spheres',&
         'current_sphere%ppd_list',ierr)

  end subroutine function_basis_init_spheres


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_init_tight_boxes(fbasis,ngwf_basis)

    !========================================================================!
    ! This subroutine initialises the TIGHT BOXES for the functions on atoms !
    ! that belong to the current node.                                       !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! fbasis (inout) : Function basis for which the tightboxes are to be     !
    !                  initialised.                                          !
    !------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000 and revised on 30/5/2001.     !
    ! Modified by Chris-Kriton Skylaris on 21/8/2003 so that it works with   !
    ! the parallel version (ONETEP).                                         !
    ! Modified by Chris-Kriton Skylaris on 17/08/2007 to work when the       !
    ! FFT box coincides with the simulation cell.                            !
    ! Adapted from ngwfs_init_tight_boxes by Nicholas Hine on 10/07/2009.    !
    !========================================================================!

    use comms, only: comms_reduce, pub_my_node_id
    use geometry,  only: POINT, geometry_distance, local_displacement
    use parallel_strategy, only : pub_num_atoms_on_node, pub_first_atom_on_node
    use simulation_cell, only : pub_cell, pub_fftbox
    use utils, only: utils_alloc_check

    implicit none

    type(FUNC_BASIS), intent(inout) :: fbasis
    type(FUNC_BASIS), optional, intent(in) :: ngwf_basis

    ! Local Variables
    real(kind=DP) :: local_1,local_2,local_3,local_distance

    integer  :: loc_1,loc_2,loc_3,ppd_loc
    integer  :: a1_neighbour,a2_neighbour,a3_neighbour
    integer  :: start1,start2,start3
    integer  :: finish1,finish2,finish3
    integer  :: ppd_a1,ppd_a2,ppd_a3,map_a1a2
    integer  :: ppd_loc1 ! ppd a1-location is always in simcell zero
    integer  :: ppd_loc2
    integer  :: ppd_loc3
    integer  :: points1
    integer  :: points2
    integer  :: points3
    integer  :: ppd_diff,ifunc,ppd,ppd_count
    integer  :: ierr
    integer  :: hubbard_atom, hubbard_proj !ddor

    type(POINT) :: current_point

    allocate(fbasis%tight_boxes(fbasis%node_num),stat=ierr)
    call utils_alloc_check('function_basis_init_tight_boxes_'//fbasis%name, &
         'fbasis%tight_boxes',ierr)

    ! ddor: Use ngwf_basis%tight_boxes on this node for the Hubbard projectors
    if (fbasis%name=='hub_projs') then

       do hubbard_atom = pub_first_atom_on_node(pub_my_node_id), &
            pub_first_atom_on_node(pub_my_node_id+1) - 1

          do hubbard_proj = fbasis%first_on_atom(hubbard_atom) - &
               fbasis%first_on_node(pub_my_node_id) + 1, &
               fbasis%first_on_atom(hubbard_atom) + &
               fbasis%num_on_atom(hubbard_atom) - &
               fbasis%first_on_node(pub_my_node_id)
             ! ddor: Use tightbox of first NGWF on a Hubbard atom
             ifunc = ngwf_basis%first_on_atom(hubbard_atom) - &
                  ngwf_basis%first_on_node(pub_my_node_id) + 1
             fbasis%tight_boxes(hubbard_proj) = ngwf_basis%tight_boxes(ifunc)

          enddo

       enddo

       return

    end if

    ! ndmh: normal version
    do ifunc=1,fbasis%node_num

       ! cks: initialise the padding of the box
       fbasis%tight_boxes(ifunc)%pad1= 0
       fbasis%tight_boxes(ifunc)%pad2= 0
       fbasis%tight_boxes(ifunc)%pad3= 0

       ! cks: initialise border ppd coordinates
       fbasis%tight_boxes(ifunc)%start_ppds1 = 2*pub_cell%n_ppds_a1
       fbasis%tight_boxes(ifunc)%finish_ppds1 = -pub_cell%n_ppds_a1

       fbasis%tight_boxes(ifunc)%start_ppds2 = 2*pub_cell%n_ppds_a2
       fbasis%tight_boxes(ifunc)%finish_ppds2 = -pub_cell%n_ppds_a2

       fbasis%tight_boxes(ifunc)%start_ppds3 = 2*pub_cell%n_ppds_a3
       fbasis%tight_boxes(ifunc)%finish_ppds3 = -pub_cell%n_ppds_a3

       ! cks: initialise the end points of the border ppds
       start1  =pub_cell%n_pt1
       finish1 =1
       start2  =pub_cell%n_pt2
       finish2 =1
       start3  =pub_cell%n_pt3
       finish3 =1

       ! ndmh: loop over all ppds in the sphere
       do ppd_count=1,fbasis%spheres(ifunc)%n_ppds_sphere

          ppd = fbasis%spheres(ifunc)%ppd_list(1,ppd_count)
          ppd_loc = fbasis%spheres(ifunc)%ppd_list(2,ppd_count)

          ! ndmh: find ppd coords of this ppd
          ppd_a3 = (ppd-1) / (pub_cell%n_ppds_a1*pub_cell%n_ppds_a2) + 1
          map_a1a2 = ppd - (pub_cell%n_ppds_a1*pub_cell%n_ppds_a2) * (ppd_a3 - 1)
          ppd_a2 = (map_a1a2 - 1) / pub_cell%n_ppds_a1 + 1
          ppd_a1 = map_a1a2 - pub_cell%n_ppds_a1 * (ppd_a2-1)

          ! cks: find on which periodic cell the ppd is located
          a1_neighbour =nint(real(ppd_loc,DP)/9.0_DP)
          a2_neighbour =nint(real(ppd_loc-9*a1_neighbour,DP)/3.0_DP )
          a3_neighbour =ppd_loc-9*a1_neighbour-3*a2_neighbour

          ! cks: express the position of the ppd in terms of its periodic
          !      image in another cell (in integer ppd coordinates)
          !      leaving the function centre where it is.
          ppd_loc1 =ppd_a1 -a1_neighbour*pub_cell%n_ppds_a1
          ppd_loc2 =ppd_a2 -a2_neighbour*pub_cell%n_ppds_a2
          ppd_loc3 =ppd_a3 -a3_neighbour*pub_cell%n_ppds_a3

          ! cks: redefine border ppd integer coordinates as necessary,
          !      re-initialising their end points whenever a border ppd
          !      is redefined
          if (ppd_loc1 < fbasis%tight_boxes(ifunc)%start_ppds1) then
             fbasis%tight_boxes(ifunc)%start_ppds1 =ppd_loc1
          endif
          if (ppd_loc2 < fbasis%tight_boxes(ifunc)%start_ppds2) then
             fbasis%tight_boxes(ifunc)%start_ppds2 =ppd_loc2
          endif
          if (ppd_loc3 < fbasis%tight_boxes(ifunc)%start_ppds3) then
             fbasis%tight_boxes(ifunc)%start_ppds3 =ppd_loc3
          endif
          if (ppd_loc1 > fbasis%tight_boxes(ifunc)%finish_ppds1) then
             fbasis%tight_boxes(ifunc)%finish_ppds1 =ppd_loc1
          endif
          if (ppd_loc2 > fbasis%tight_boxes(ifunc)%finish_ppds2) then
             fbasis%tight_boxes(ifunc)%finish_ppds2 =ppd_loc2
          endif
          if (ppd_loc3 > fbasis%tight_boxes(ifunc)%finish_ppds3) then
             fbasis%tight_boxes(ifunc)%finish_ppds3 =ppd_loc3
          endif

       end do

       ! cks: loop over all ppds in the tightbox
       ppd=0
       do ppd_a3=fbasis%tight_boxes(ifunc)%start_ppds3, &
            fbasis%tight_boxes(ifunc)%finish_ppds3
          points3=(ppd_a3-1)*pub_cell%n_pt3

          do ppd_a2=fbasis%tight_boxes(ifunc)%start_ppds2, &
               fbasis%tight_boxes(ifunc)%finish_ppds2
             points2 =(ppd_a2-1)*pub_cell%n_pt2

             do ppd_a1=fbasis%tight_boxes(ifunc)%start_ppds1, &
                  fbasis%tight_boxes(ifunc)%finish_ppds1
                points1 =(ppd_a1-1)*pub_cell%n_pt1


                ! cks: loop over all the points of ppd
                do loc_3=0,pub_cell%n_pt3-1
                   local_3 = real(loc_3+points3, DP)*pub_cell%d3

                   do loc_2=0,pub_cell%n_pt2-1
                      local_2 = real(loc_2+points2, DP)*pub_cell%d2

                      do loc_1=0,pub_cell%n_pt1-1
                         local_1 = real(loc_1+points1, DP)*pub_cell%d1

                         current_point = local_displacement( &
                              pub_cell%a1_unit, pub_cell%a2_unit, &
                              pub_cell%a3_unit, local_1, local_2, local_3)

                         local_distance =geometry_distance( &
                              fbasis%spheres(ifunc)%centre, current_point)

                         ! cks: set the limits, in the points of the current ppd,
                         ! of the parellelepiped that exactly inscribes the
                         ! function region into it by testing every single point
                         ! of the ppd that belongs to the function region.
                         if (  local_distance.le. &
                              ( fbasis%spheres(ifunc)%radius )  ) then

                            if ( (loc_1+1.lt.start1).and. &
                                 (ppd_a1.eq.fbasis%tight_boxes(ifunc)%start_ppds1)) &
                                 start1=loc_1+1
                            if ( (loc_1+1.gt.finish1).and.&
                                 (ppd_a1.eq.fbasis%tight_boxes(ifunc)%finish_ppds1) ) &
                                 finish1=loc_1+1
                            if ( (loc_2+1.lt.start2).and.&
                                 (ppd_a2.eq.fbasis%tight_boxes(ifunc)%start_ppds2) ) &
                                 start2=loc_2+1
                            if ( (loc_2+1.gt.finish2).and.&
                                 (ppd_a2.eq.fbasis%tight_boxes(ifunc)%finish_ppds2) ) &
                                 finish2=loc_2+1
                            if ( (loc_3+1.lt.start3).and.&
                                 (ppd_a3.eq.fbasis%tight_boxes(ifunc)%start_ppds3) ) &
                                 start3=loc_3+1
                            if ( (loc_3+1.gt.finish3).and.&
                                 (ppd_a3.eq.fbasis%tight_boxes(ifunc)%finish_ppds3) ) &
                                 finish3=loc_3+1


                         endif

                      enddo  ! loc_1
                   enddo  ! loc_2
                enddo  ! loc_3

             enddo  ! ppd_a1
          enddo  ! ppd_a2
       enddo  ! ppd_a3

       fbasis%tight_boxes(ifunc)%START_PTS1  = start1
       fbasis%tight_boxes(ifunc)%FINISH_PTS1 = finish1
       fbasis%tight_boxes(ifunc)%START_PTS2  = start2
       fbasis%tight_boxes(ifunc)%FINISH_PTS2 = finish2
       fbasis%tight_boxes(ifunc)%START_PTS3  = start3
       fbasis%tight_boxes(ifunc)%FINISH_PTS3 = finish3

       ! cks: now finish box specification by calculating the number of (grid)
       !      points in every lattice vector direction.
       if (pub_fftbox%coin1) then
          fbasis%tight_boxes(ifunc)%start_ppds1  = 1
          fbasis%tight_boxes(ifunc)%finish_ppds1 = pub_cell%n_ppds_a1
          fbasis%tight_boxes(ifunc)%start_pts1   = 1
          fbasis%tight_boxes(ifunc)%finish_pts1  = pub_cell%n_pt1
          fbasis%tight_boxes(ifunc)%tight_pts1   = pub_cell%total_pt1
       else
          ppd_diff = fbasis%tight_boxes(ifunc)%finish_ppds1 &
               - fbasis%tight_boxes(ifunc)%start_ppds1
          fbasis%tight_boxes(ifunc)%tight_pts1 = &
               fbasis%tight_boxes(ifunc)%finish_pts1 &
               - fbasis%tight_boxes(ifunc)%start_pts1 &
               + 1 + ppd_diff*pub_cell%n_pt1 + 2*fbasis%tight_boxes(ifunc)%pad1
       endif


       if (pub_fftbox%coin2) then
          fbasis%tight_boxes(ifunc)%start_ppds2  = 1
          fbasis%tight_boxes(ifunc)%finish_ppds2 = pub_cell%n_ppds_a2
          fbasis%tight_boxes(ifunc)%start_pts2   = 1
          fbasis%tight_boxes(ifunc)%finish_pts2  = pub_cell%n_pt2
          fbasis%tight_boxes(ifunc)%tight_pts2   = pub_cell%total_pt2
       else
          ppd_diff = fbasis%tight_boxes(ifunc)%finish_ppds2 &
               - fbasis%tight_boxes(ifunc)%start_ppds2
          fbasis%tight_boxes(ifunc)%tight_pts2 = &
               fbasis%tight_boxes(ifunc)%finish_pts2 &
               -fbasis%tight_boxes(ifunc)%start_pts2 &
               + 1 + ppd_diff*pub_cell%n_pt2 + 2*fbasis%tight_boxes(ifunc)%pad2
       endif


       if (pub_fftbox%coin3) then
          fbasis%tight_boxes(ifunc)%start_ppds3  = 1
          fbasis%tight_boxes(ifunc)%finish_ppds3 = pub_cell%n_ppds_a3
          fbasis%tight_boxes(ifunc)%start_pts3   = 1
          fbasis%tight_boxes(ifunc)%finish_pts3  = pub_cell%n_pt3
          fbasis%tight_boxes(ifunc)%tight_pts3   = pub_cell%total_pt3
       else
          ppd_diff = fbasis%tight_boxes(ifunc)%finish_ppds3 &
               - fbasis%tight_boxes(ifunc)%start_ppds3
          fbasis%tight_boxes(ifunc)%tight_pts3 = &
               fbasis%tight_boxes(ifunc)%finish_pts3 &
               - fbasis%tight_boxes(ifunc)%start_pts3 &
               + 1 + ppd_diff*pub_cell%n_pt3 + 2*fbasis%tight_boxes(ifunc)%pad3
       endif

    enddo

    ! ndmh: find maximum size of tightboxes across all functions in set
    fbasis%maxtight_pts1 = maxval(fbasis%tight_boxes(:)%tight_pts1)
    fbasis%maxtight_pts2 = maxval(fbasis%tight_boxes(:)%tight_pts2)
    fbasis%maxtight_pts3 = maxval(fbasis%tight_boxes(:)%tight_pts3)
    call comms_reduce('MAX',fbasis%maxtight_pts1)
    call comms_reduce('MAX',fbasis%maxtight_pts2)
    call comms_reduce('MAX',fbasis%maxtight_pts3)

  end subroutine function_basis_init_tight_boxes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_gath_all_tbs(fbasis, ngwf_basis)   ! input

    !=================================================================!
    ! This subroutine initialises the fbasis%all_tbs array that holds !
    ! together all tightboxes of all nodes for these functions.       !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 24/03/2004.                 !
    ! Modified by Peter Haynes on 26/5/2005.                          !
    ! Moved to function_basis_mod by Nicholas Hine, and adapted       !
    ! to take FUNC_BASIS argument 15/07/2009.                         !
    !=================================================================!

    use comms, only : comms_bcast, pub_my_node_id, pub_total_num_nodes
    use parallel_strategy, only : pub_first_atom_on_node, pub_num_atoms_on_node
    use utils, only: utils_dealloc_check, utils_alloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: fbasis
    type(FUNC_BASIS), optional, intent(in) :: ngwf_basis

    ! Local variables
    integer :: local_func             ! Index of function on local node
    integer :: global_func            ! Global index of function
    integer :: node                   ! Node counter
    integer :: ierr                   ! Error flag
    integer, allocatable :: buf(:,:)  ! Buffer for comms
    integer :: atom, proj             ! ddor: Used for Hubbard projector tbs

    ! pdh: allocate memory for all tight boxes
    if (allocated(fbasis%all_tbs)) then
       if (size(fbasis%all_tbs) /= fbasis%num) then
          deallocate(fbasis%all_tbs, stat=ierr)
          call utils_dealloc_check('function_basis_gath_all_tbs_'//fbasis%name,&
               'fbasis%all_tbs', ierr)
       end if
    end if

    if (.not. allocated(fbasis%all_tbs)) then
       allocate(fbasis%all_tbs(fbasis%num),stat=ierr)
       call utils_alloc_check('function_basis_gath_all_tbs_'//fbasis%name, &
            'fbasis%all_tbs',ierr)
    end if

    !ddor: Use ngwf_basis%all_tbs for the Hubbard projectors
    if (fbasis%name=='hub_projs') then

       ! ddor: loop over nodes
       do node=0,pub_total_num_nodes-1

          do atom = pub_first_atom_on_node(node), &
               pub_first_atom_on_node(node) + pub_num_atoms_on_node(node) - 1

             do proj = fbasis%first_on_atom(atom), &
                  fbasis%first_on_atom(atom) + fbasis%num_on_atom(atom) - 1

                ! ddor: Use tightbox of first NGWF on a Hubbard atom
                fbasis%all_tbs(proj) = &
                     &ngwf_basis%all_tbs( ngwf_basis%first_on_atom(atom) )

             enddo

          enddo

       enddo

    else

       ! pdh: allocate comms buffer
       allocate(buf(18,maxval(fbasis%num_on_node)),stat=ierr)
       call utils_alloc_check('function_basis_gath_all_tbs','buf',ierr)

       ! pdh: loop over nodes
       do node=0,pub_total_num_nodes-1

          ! pdh: pack up tight boxes on local_node into buf
          if (node == pub_my_node_id) then
             do local_func=1,fbasis%num_on_node(node)
                buf( 1,local_func) = fbasis%tight_boxes(local_func)%pad1
                buf( 2,local_func) = fbasis%tight_boxes(local_func)%pad2
                buf( 3,local_func) = fbasis%tight_boxes(local_func)%pad3
                buf( 4,local_func) = fbasis%tight_boxes(local_func)%start_ppds1
                buf( 5,local_func) = fbasis%tight_boxes(local_func)%start_ppds2
                buf( 6,local_func) = fbasis%tight_boxes(local_func)%start_ppds3
                buf( 7,local_func) = fbasis%tight_boxes(local_func)%finish_ppds1
                buf( 8,local_func) = fbasis%tight_boxes(local_func)%finish_ppds2
                buf( 9,local_func) = fbasis%tight_boxes(local_func)%finish_ppds3
                buf(10,local_func) = fbasis%tight_boxes(local_func)%start_pts1
                buf(11,local_func) = fbasis%tight_boxes(local_func)%start_pts2
                buf(12,local_func) = fbasis%tight_boxes(local_func)%start_pts3
                buf(13,local_func) = fbasis%tight_boxes(local_func)%finish_pts1
                buf(14,local_func) = fbasis%tight_boxes(local_func)%finish_pts2
                buf(15,local_func) = fbasis%tight_boxes(local_func)%finish_pts3
                buf(16,local_func) = fbasis%tight_boxes(local_func)%tight_pts1
                buf(17,local_func) = fbasis%tight_boxes(local_func)%tight_pts2
                buf(18,local_func) = fbasis%tight_boxes(local_func)%tight_pts3
             end do
          end if

          ! pdh: broadcast tight-boxes on local_node
          call comms_bcast(node,buf,fbasis%num_on_node(node)*18)

          ! pdh: loop over funcs
          local_func = 0
          do global_func=1,fbasis%num

             ! pdh: copy into global position
             if (node == fbasis%node_of_func(global_func)) then
                local_func = local_func + 1

                fbasis%all_tbs(global_func)%pad1         = buf( 1,local_func)
                fbasis%all_tbs(global_func)%pad2         = buf( 2,local_func)
                fbasis%all_tbs(global_func)%pad3         = buf( 3,local_func)
                fbasis%all_tbs(global_func)%start_ppds1  = buf( 4,local_func)
                fbasis%all_tbs(global_func)%start_ppds2  = buf( 5,local_func)
                fbasis%all_tbs(global_func)%start_ppds3  = buf( 6,local_func)
                fbasis%all_tbs(global_func)%finish_ppds1 = buf( 7,local_func)
                fbasis%all_tbs(global_func)%finish_ppds2 = buf( 8,local_func)
                fbasis%all_tbs(global_func)%finish_ppds3 = buf( 9,local_func)
                fbasis%all_tbs(global_func)%start_pts1   = buf(10,local_func)
                fbasis%all_tbs(global_func)%start_pts2   = buf(11,local_func)
                fbasis%all_tbs(global_func)%start_pts3   = buf(12,local_func)
                fbasis%all_tbs(global_func)%finish_pts1  = buf(13,local_func)
                fbasis%all_tbs(global_func)%finish_pts2  = buf(14,local_func)
                fbasis%all_tbs(global_func)%finish_pts3  = buf(15,local_func)
                fbasis%all_tbs(global_func)%tight_pts1   = buf(16,local_func)
                fbasis%all_tbs(global_func)%tight_pts2   = buf(17,local_func)
                fbasis%all_tbs(global_func)%tight_pts3   = buf(18,local_func)
             end if
          end do

       end do   ! loop over nodes

       ! Deallocate comms buffer
       deallocate(buf,stat=ierr)
       call utils_dealloc_check('function_basis_gath_all_tbs','buf',ierr)

    endif

  end subroutine function_basis_gath_all_tbs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_init_uni_tb(ngwf_basis,ngwf_basis2) ! input

    !=================================================================!
    ! This subroutine finds the size of the universal NGWF tightbox   !
    ! and initialises its reciprocal space grid if required.          !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/07/2009 from bits of the         !
    ! routine basis_gath_all_tbs written by Chris-Kriton Skylaris,    !
    ! and modified by Alvaro Ruiz-Serrano.                            !
    !=================================================================!

    use geometry, only: point, operator(*)
    use rundat, only: pub_tightbox_fft_coarse, pub_tightbox_fft_fine, &
         pub_usehfx
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3, pub_tb_recip_grid, pub_tb_recip_grid_fine
    use utils, only: utils_alloc_check

    implicit none

    ! Argument
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    ! lr408: Optional argument for the case of conduction states to ensure universal
    ! lr408: NGWF tightbox is big enough for both valence and conduction NGWFs
    type(FUNC_BASIS), optional, intent(in) :: ngwf_basis2

    ! ars : the universal tightbox will have an odd number of points
    ! ars : so it can be used to perform FFT without problems
    pub_maxtight_pts1 = ngwf_basis%maxtight_pts1
    pub_maxtight_pts2 = ngwf_basis%maxtight_pts2
    pub_maxtight_pts3 = ngwf_basis%maxtight_pts3


    ! lr408: Check maximum size needed to fit both conduction and valence NGWFs
    if (present(ngwf_basis2)) then

       pub_maxtight_pts1 = max(ngwf_basis2%maxtight_pts1,pub_maxtight_pts1)
       pub_maxtight_pts2 = max(ngwf_basis2%maxtight_pts2,pub_maxtight_pts2)
       pub_maxtight_pts3 = max(ngwf_basis2%maxtight_pts3,pub_maxtight_pts3)

    end if

    if (pub_tightbox_fft_coarse.or.pub_tightbox_fft_fine) then

       if (mod(pub_maxtight_pts1,2).eq.0) &
            pub_maxtight_pts1 = pub_maxtight_pts1 + 1

       if (mod(pub_maxtight_pts2,2).eq.0) &
            pub_maxtight_pts2 = pub_maxtight_pts2 + 1

       if (mod(pub_maxtight_pts3,2).eq.0) &
            pub_maxtight_pts3 = pub_maxtight_pts3 + 1

    endif

    call internal_recip_grid_init

  contains

    subroutine internal_recip_grid_init

      !=====================================================!
      ! This subroutine initialises pub_tb_recip_grid and   !
      ! pub_tb_recip_grid_fine if necessary. These are used !
      ! to create spherical waves in tightboxes.            !
      !-----------------------------------------------------!
      ! Based on simulation_cell_fftbox_init.               !
      ! Written by Quintin Hill on 20/03/2009.              !
      !=====================================================!

      implicit none

      type(point)   :: tb_b1, tb_b2, tb_b3 ! Reciprocal lattice vectors
      real(kind=DP) :: b1(3), b2(3), b3(3) ! Reciprocal lattice vectors
      real(kind=DP) :: g3(3),g23(3),g(3)
      real(kind=DP) :: gsq
      integer :: i1,i2,i3
      integer :: k1,k2,k3
      integer :: n1half,n2half,n3half
      integer :: ierr

      tb_b1 =( real(pub_cell%total_pt1, kind=DP)/&
           real(pub_maxtight_pts1, kind=DP) ) * pub_cell%b1
      tb_b2 =( real(pub_cell%total_pt2, kind=DP)/&
           real(pub_maxtight_pts2, kind=DP) ) * pub_cell%b2
      tb_b3 =( real(pub_cell%total_pt3, kind=DP)/&
           real(pub_maxtight_pts3, kind=DP) )  * pub_cell%b3

      ! local copies of FFT box reciprocal lattice vectors
      b1(1) = tb_b1%x ; b1(2) = tb_b1%y ; b1(3) = tb_b1%z
      b2(1) = tb_b2%x ; b2(2) = tb_b2%y ; b2(3) = tb_b2%z
      b3(1) = tb_b3%x ; b3(2) = tb_b3%y ; b3(3) = tb_b3%z

      if (pub_tightbox_fft_coarse) then

         allocate(pub_tb_recip_grid(5,pub_maxtight_pts1, &
              pub_maxtight_pts2,pub_maxtight_pts3),stat=ierr)
         call utils_alloc_check('function_basis_init_uni_tb', &
              'pub_tb_recip_grid', ierr)

         !qoh:  loop over tightbox reciprocal grid
         pub_tb_recip_grid = 0.0_DP
         n1half = pub_maxtight_pts1/2+1
         n2half = pub_maxtight_pts2/2+1
         n3half = pub_maxtight_pts3/2+1
         do i3=1,pub_maxtight_pts3
            if (i3 > n3half) then
               k3 = i3 - pub_maxtight_pts3 - 1
            else
               k3 = i3 - 1
            end if
            g3 = k3 * b3
            do i2=1,pub_maxtight_pts2
               if (i2 > n2half) then
                  k2 = i2 - pub_maxtight_pts2 - 1
               else
                  k2 = i2 - 1
               end if
               g23 = g3 + k2 * b2
               do i1=1,pub_maxtight_pts1
                  if (i1 > n1half) then
                     k1 = i1 - pub_maxtight_pts1 - 1
                  else
                     k1 = i1 - 1
                  end if
                  g = g23 + k1 * b1
                  pub_tb_recip_grid(1:3,i1,i2,i3) = g
                  gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
                  pub_tb_recip_grid(4,i1,i2,i3) = sqrt(gsq)
                  pub_tb_recip_grid(5,i1,i2,i3) = 0.5_DP * gsq
               end do
            end do
         end do

      end if

      if (pub_tightbox_fft_fine .and. pub_usehfx) then

         allocate(pub_tb_recip_grid_fine(5,pub_maxtight_pts1*2,&
              pub_maxtight_pts2*2,pub_maxtight_pts3*2),stat=ierr)
         call utils_alloc_check('function_basis_init_uni_tb', &
              'pub_tb_recip_grid_fine',ierr)
         pub_tb_recip_grid_fine = 0.0_DP

         n1half = pub_maxtight_pts1+1
         n2half = pub_maxtight_pts2+1
         n3half = pub_maxtight_pts3+1

         !qoh:  loop over tightbox reciprocal grid
         do i3=1,pub_maxtight_pts3*2
            if (i3 > n3half) then
               k3 = i3 - pub_maxtight_pts3*2 - 1
            else
               k3 = i3 - 1
            end if
            g3 = k3 * b3
            do i2=1,pub_maxtight_pts2*2
               if (i2 > n2half) then
                  k2 = i2 - pub_maxtight_pts2*2 - 1
               else
                  k2 = i2 - 1
               end if
               g23 = g3 + k2 * b2
               do i1=1,pub_maxtight_pts1*2
                  if (i1 > n1half) then
                     k1 = i1 - pub_maxtight_pts1*2 - 1
                  else
                     k1 = i1 - 1
                  end if
                  g = g23 + k1 * b1
                  pub_tb_recip_grid_fine(1:3,i1,i2,i3) = g
                  gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
                  pub_tb_recip_grid_fine(4,i1,i2,i3) = sqrt(gsq)
                  pub_tb_recip_grid_fine(5,i1,i2,i3) = 0.5_DP * gsq
               end do
            end do
         end do
      end if

    end subroutine internal_recip_grid_init

  end subroutine function_basis_init_uni_tb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_exit_uni_tb

    !=================================================================!
    ! This subroutine deallocates arrays relating to the universal    !
    ! NGWF tightbox.                                                  !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/07/2009 from bits of the         !
    ! routine basis_exit by Chris-Kriton Skylaris and Quintin Hill.   !
    !=================================================================!

    use geometry, only: point, operator(*)
    use simulation_cell, only: pub_tb_recip_grid, pub_tb_recip_grid_fine
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    integer :: ierr

    ! qoh: Deallocate tightbox reciprocal grids
    if (allocated(pub_tb_recip_grid)) then
       deallocate(pub_tb_recip_grid,stat=ierr)
       call utils_dealloc_check('function_basis_exit_uni_tb', &
            'pub_tb_recip_grid', ierr)
    end if

    if (allocated(pub_tb_recip_grid_fine)) then
       deallocate(pub_tb_recip_grid_fine,stat=ierr)
       call utils_dealloc_check('function_basis_exit_uni_tb', &
            'pub_tb_recip_grid_fine', ierr)
    end if

  end subroutine function_basis_exit_uni_tb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function function_basis_est_num_psincs(fbasis)

    !========================================================================!
    ! This function returns an estimate of the total number of psincs within !
    ! all the spheres of a given function basis, based on the sum of their   !
    ! volumes.                                                               !
    !------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 16/12/2004.                        !
    ! Moved to function_basis by Nicholas Hine on 28/04/2010.                !
    !========================================================================!

    use comms, only: comms_reduce
    use constants, only: DP, PI
    use simulation_cell, only: pub_cell

    implicit none

    ! ndmh: Arguments
    type(FUNC_BASIS) :: fbasis

    ! cks: Local Variables
    integer :: row ! loop counter
    real(kind =DP) :: vol  ! total volume of all NGWF spheres

    ! cks: calculate the sum of all spheres' volumes
    vol = 0.0_DP
    do row=1,fbasis%node_num
       vol = vol +(4.0_DP/3.0_DP)*PI*(fbasis%spheres(row)%radius**3)
    enddo
    call comms_reduce('SUM', vol)

    function_basis_est_num_psincs = nint(vol/pub_cell%weight)

  end function function_basis_est_num_psincs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_estimate_size(fbasis,local_size,global_size)

    !=========================================================================!
    ! This subroutine estimates the total size of an allocated function_basis.!
    ! Test/Prototype for eventual widespread memory-usage-estimation code.    !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !    fbasis       (in): Function basis type whose size is to be estimated.!
    !    local_size  (out): Size on this node.                                !
    !    global_size (out): Size on all nodes.                                !
    !-------------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine in November 2009.          !
    !=========================================================================!

    use comms, only: comms_reduce
    use constants, only: int_size, real_size
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    integer, intent(out) :: local_size, global_size

    ! Local Variables
    integer :: distr_array_size
    integer :: spheres_array_size
    integer :: tightbox_array_size
    integer :: ifunc

    ! ndmh: size of basic information (ints and characters)
    distr_array_size = 7*int_size + 20

    ! ndmh: size of allocated arrays
    distr_array_size = distr_array_size + int_size * &
         ( size(fbasis%num_on_node) &
         + size(fbasis%first_on_node) &
         + size(fbasis%num_on_atom) &
         + size(fbasis%first_on_atom) &
         + size(fbasis%node_of_func) &
         + size(fbasis%atom_of_func) &
         + size(fbasis%n_ppds_sphere))

    ! ndmh: size of spheres array
    spheres_array_size = 0
    do ifunc=1,fbasis%node_num
       spheres_array_size = spheres_array_size + real_size + 2*int_size + &
            int_size*2*fbasis%spheres(ifunc)%n_ppds_sphere
    end do

    ! ndmh: size of tightboxes arrays
    tightbox_array_size = 0
    if (allocated(fbasis%tight_boxes)) then
       tightbox_array_size = tightbox_array_size + 18*int_size*fbasis%node_num
    end if
    if (allocated(fbasis%all_tbs)) then
       tightbox_array_size = tightbox_array_size + 18*int_size*fbasis%num
    end if

    local_size = distr_array_size + spheres_array_size + tightbox_array_size

    ! ndmh: find global size summed over all nodes
    global_size = local_size
    call comms_reduce('SUM',global_size)

  end subroutine function_basis_estimate_size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_tightbox_to_ppds(tb_in, &
       read_tb_orig1,read_tb_orig2,read_tb_orig3, &
       ifunc,tb_node_id,tb_n1,tb_n2,tb_n3,funcs_on_grid,fbasis, &
       tb_start1,tb_start2,tb_start3,fftbox_complex,fftbox_complex_shifted, &
       fftbox_real,sendbuf,recvbuf)

    !=========================================================================!
    ! This subroutine fills the ppds of a given function (which can be on     !
    ! any node) with a tightbox from (potentially) another node. The          !
    ! arguments specify the original location of the function within the      !
    ! input tightbox, which may differ from the new location of the function  !
    ! with respect to the grid.                                               !
    ! In such cases, the FFTbox is shifted by the Fourier Transform/Phase     !
    ! Shift/Fourier Transform. This routine is intended as part of the NGWF   !
    ! IO system, and as such it effectively serialises the comms. It should   !
    ! NOT be used as part of the main code, as it would be very slow and      !
    ! scale poorly with processor count as currently implemented.             !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !  tb_n1,tb_n2,tb_n3 (input)         : Tightbox sizes in 1,2,3-directions !
    !  ifunc (input)                     : Index of the function required     !
    !  tb_node_id (input)                : Node with the tightbox to extract  !
    !  fbasis (input)                    : Function basis describing funcs    !
    !  funcs_on_grid (input)             : PPD data of functions              !
    !  sendbuf, recvbuf (input/output)   : Buffers for send/recv of tb origin !
    !  tb_in (input/output)              : Origin tightbox array              !
    !  tb_start1,2,3 (in)                : Origin of a tightbox within FFTbox !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2010, reusing bits of restart_mod  !
    ! routines (eg restart_ngwfs_tightbox_input) originally written by        !
    ! Chris-Kriton Skylaris in March 2004.                                    !
    !=========================================================================!

    use basis, only: basis_copy_tightbox_to_fftbox, basis_clean_function, &
         basis_extract_function_from_box, basis_function_origin_wrt_tb, &
         basis_phase_on_fftbox_recip
    use comms, only: comms_recv, comms_send, pub_my_node_id
    use constants, only: stdout
    use fourier, only: fourier_apply_box
    use parallel_strategy, only: pub_node_of_atom
    use simulation_cell, only: pub_cell, pub_fftbox

    implicit none

    ! Arguments
    integer, intent(in)           :: tb_n1,tb_n2,tb_n3
    integer, intent(in)           :: tb_start1,tb_start2,tb_start3
    integer, intent(in)           :: ifunc
    integer, intent(in)           :: tb_node_id
    type(FUNC_BASIS), intent(in)  :: fbasis
    real(kind=DP), intent(inout)  :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(inout)  :: sendbuf(3), recvbuf(3)
    real(kind=DP), intent(inout)  :: tb_in(tb_n1,tb_n2,tb_n3)
    real(kind=DP), intent(inout)  :: read_tb_orig1, read_tb_orig2, read_tb_orig3
    real(kind=DP), intent(out)    :: fftbox_real(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    complex(kind=DP), intent(out) :: fftbox_complex(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    complex(kind=DP), intent(out) :: fftbox_complex_shifted( &
         pub_fftbox%total_ld1,pub_fftbox%total_ld2,pub_fftbox%total_pt3)

    ! Locals
    integer :: loc_ifunc
    integer :: fft_n1, fft_n2, fft_ld1, fft_ld2, fft_n3
    real(kind=DP) :: tb_orig1, tb_orig2, tb_orig3

    fft_ld1 = pub_fftbox%total_ld1
    fft_ld2 = pub_fftbox%total_ld2
    fft_n1 = pub_fftbox%total_pt1
    fft_n2 = pub_fftbox%total_pt2
    fft_n3 = pub_fftbox%total_pt3

    ! ndmh: local node needs the ppds of this function, so must receive them
    ! ndmh: into its tightbox buffer (and receive the read offset from the
    ! ndmh: origin) and extract the ppds.
    if (pub_my_node_id==fbasis%node_of_func(ifunc)) then

       loc_ifunc = ifunc - fbasis%first_on_node(pub_my_node_id) + 1

       ! ndmh: recv tightbox data, and origin of atom in read tightbox,
       ! ndmh: if node holding tb is not local node
       if (.not. (pub_my_node_id==tb_node_id)) then
          call comms_recv(tb_node_id,recvbuf,tag=ifunc)
          read_tb_orig1 = recvbuf(1)
          read_tb_orig2 = recvbuf(2)
          read_tb_orig3 = recvbuf(3)
          call comms_recv(tb_node_id,tb_in,tag=ifunc+fbasis%num)
       end if

       ! cks: find origin of NGWF wrt to tightbox in terms of number of
       ! cks: grid points in each lattice vector direction
       call basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3, &
            fbasis%spheres(loc_ifunc)%centre,fbasis%tight_boxes(loc_ifunc))

       ! ndmh: check if the origin has changed since last time
       if ((read_tb_orig1/=tb_orig1) .or. &
            (read_tb_orig2/=tb_orig2) .or. &
            (read_tb_orig3/=tb_orig3)) then

          ! cks: apply phase factor to translate NGWF from beginning of
          ! cks: coordinates in universal tightbox to the atomic centre
          ! cks: of NGWF
          !============== PHASE SHIFT ======================================

          ! cks: initialise complex fftbox
          fftbox_complex = (0.0_DP,0.0_DP)

          ! ndmh: deposit tightbox at (tb_start1,tb_start2,tb_start3)
          ! ndmh: within complex FFTbox
          fftbox_complex(tb_start1:tb_start1+tb_n1-1, &
               tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
               cmplx(tb_in(1:tb_n1,1:tb_n2,1:tb_n3),0.0_DP,kind=DP)

          ! ndmh: in-place FFT to reciprocal space
          call fourier_apply_box('Coarse','Forward',fftbox_complex)

          ! ndmh: apply phase factor to move function by required shift
          call basis_phase_on_fftbox_recip(fftbox_complex_shifted, &
               fftbox_complex, fft_n1, fft_n2, fft_n3, fft_ld1, fft_ld2, &
               read_tb_orig1-tb_orig1, read_tb_orig2-tb_orig2, &
               read_tb_orig3-tb_orig3)

          ! cks: g=0 element must be real
          fftbox_complex_shifted(1,1,1) = &
               cmplx(real(fftbox_complex_shifted(1,1,1),kind=DP),0.0_DP,kind=DP)

          ! ndmh: in-place FFT back to real space
          call fourier_apply_box('Coarse','Backward',fftbox_complex_shifted)

          ! ndmh: transfer to real-valued FFTbox for extraction
          fftbox_real = real(fftbox_complex_shifted,kind=DP)

          !=========== END PHASE SHIFT ======================================

       else ! ndmh: no need to shift function

          fftbox_real = 0.0_DP

          ! ndmh: deposit tightbox at (tb_start1,tb_start2,tb_start3)
          ! ndmh: within real FFTbox
          fftbox_real(tb_start1:tb_start1+tb_n1-1, &
               tb_start2:tb_start2+tb_n2-1,tb_start3:tb_start3+tb_n3-1) = &
               tb_in(1:tb_n1,1:tb_n2,1:tb_n3)

       end if

       ! ndmh: extract ppds of function from FFTbox
       call basis_extract_function_from_box(funcs_on_grid,fft_ld1,fft_ld2, &
            fft_n3, fftbox_real, fbasis%spheres(loc_ifunc), &
            fbasis%tight_boxes(loc_ifunc), tb_start1, tb_start2, tb_start3, &
            fbasis%spheres(loc_ifunc)%offset)

       ! ndmh: zero points outside localisation sphere
       call basis_clean_function(funcs_on_grid, fbasis%spheres(loc_ifunc), &
            fbasis%n_ppds)

       ! ndmh: local node does not have the ppds of this function. Do nothing
       ! ndmh: unless this is the sending node, in which case, send the
       ! ndmh: tightbox and its position within the tightbox as written.
    else

       if (pub_my_node_id==tb_node_id) then
          ! ndmh: send the read function to its node
          sendbuf(1) = read_tb_orig1
          sendbuf(2) = read_tb_orig2
          sendbuf(3) = read_tb_orig3
          call comms_send(fbasis%node_of_func(ifunc),sendbuf,tag=ifunc)
          call comms_send(fbasis%node_of_func(ifunc),tb_in,tag=ifunc+fbasis%num)
       endif

    endif

  end subroutine function_basis_tightbox_to_ppds


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_ppds_to_tightbox(tb_out,tb_orig1,tb_orig2,tb_orig3, &
       ifunc,tb_node_id,tb_n1,tb_n2,tb_n3,funcs_on_grid,fbasis,sendbuf,recvbuf)

    !=========================================================================!
    ! This subroutine places the ppds of a given function (which can be on    !
    ! any node) into a tightbox on another node. This routine is intended as  !
    ! part of the NGWF IO system, and as such it effectively serialises the   !
    ! comms. It should NOT be used as part of the main code, as that would    !
    ! be very slow in such a context.                                         !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !  tb_n1,tb_n2,tb_n3 (input)         : Tightbox sizes in 1,2,3-directions !
    !  ifunc (input)                     : Index of the function required     !
    !  tb_node_id (input)                : Node of the destination tightbox   !
    !  fbasis (input)                    : Function basis describing funcs    !
    !  funcs_on_grid (input)             : PPD data of functions              !
    !  sendbuf, recvbuf (input/output)   : Buffers for send/recv of tb origin !
    !  tb_out (output)                   : Destination tightbox array         !
    !  tb_orig1, tb_orig2, tb_orig3 (out): Origin of atom within tightbox     !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2010, reusing bits of restart_mod  !
    ! routines (eg restart_ngwfs_tightbox_output) by Chris-Kriton Skylaris.   !
    !=========================================================================!

    use basis, only: basis_copy_function_to_box, &
         basis_function_origin_wrt_tb
    use comms, only: comms_recv, comms_send, pub_my_node_id
    use parallel_strategy, only: pub_node_of_atom
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    integer, intent(in)          :: tb_n1,tb_n2,tb_n3
    integer, intent(in)          :: ifunc
    integer, intent(in)          :: tb_node_id
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(in)    :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(inout) :: sendbuf(3), recvbuf(3)
    real(kind=DP), intent(out)   :: tb_out(tb_n1,tb_n2,tb_n3)
    real(kind=DP), intent(out)   :: tb_orig1, tb_orig2, tb_orig3

    ! Locals
    integer :: loc_ifunc

    ! ndmh: local node has the ppds of this function, so must copy them to its
    ! ndmh: tightbox buffer (and calculate the offset from the origin) and
    ! ndmh: send it (unless it is also the node that is receiving it).
    if (pub_my_node_id==fbasis%node_of_func(ifunc)) then

       loc_ifunc = ifunc - fbasis%first_on_node(pub_my_node_id) + 1

       ! Initialise to zero
       tb_out = 0.0_DP

       ! Put function in tightbox
       call basis_copy_function_to_box(tb_out,tb_n1,tb_n2,tb_n3,1,1,1, &
            fbasis%tight_boxes(loc_ifunc), funcs_on_grid, &
            fbasis%spheres(loc_ifunc))

       ! cks: find origin of NGWF wrt to tightbox in terms of number of
       ! cks: grid points in each lattice vector direction
       call basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3,  &
            fbasis%spheres(loc_ifunc)%centre, &
            fbasis%tight_boxes(loc_ifunc))

       ! ndmh: send, if destination node is not local
       if (.not. (pub_my_node_id==tb_node_id)) then
          sendbuf(1) = tb_orig1
          sendbuf(2) = tb_orig2
          sendbuf(3) = tb_orig3
          call comms_send(tb_node_id,sendbuf,tag=ifunc)
          call comms_send(tb_node_id,tb_out,tag=ifunc+fbasis%num)
       end if

       ! ndmh: local node does not have the ppds of this function. Do nothing
       ! ndmh: unless this is the receiving node (in which case, receive the
       ! ndmh: tightbox.
    else

       if (pub_my_node_id==tb_node_id) then
          ! ndmh: receive the function from its node
          call comms_recv(fbasis%node_of_func(ifunc),recvbuf,tag=ifunc)
          tb_orig1 = recvbuf(1)
          tb_orig2 = recvbuf(2)
          tb_orig3 = recvbuf(3)
          call comms_recv(fbasis%node_of_func(ifunc),tb_out, &
               tag=ifunc+fbasis%num)
       endif

    endif

  end subroutine function_basis_ppds_to_tightbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_batch_row_plan(num_plan_steps,my_plan,idx_len, &
       sparse_idx,mat,local_start,local_end,scheme,reqs)

    !=========================================================================!
    ! This subroutine creates a list of row, column pairs of funcs that need  !
    ! to be calculated on the local node. This can then be used to execute    !
    ! row sums such as density_batch_row_sums and ngwf_gradient_fb_sums in a  !
    ! maximally efficient manner                                              !
    !-------------------------------------------------------------------------!
    !  Arguments:                                                             !
    !    num_plan_steps (inout) : on input: number of steps in my_plan array  !
    !                             on output: number of steps actually needed  !
    !    plan(1,:,:)    (out) : rows of row/col pairs this node will          !
    !                           need to calculate                             !
    !    plan(2,:,:)    (out) : cols of row/col pairs this node will          !
    !                           need to calculate                             !
    !    idx_len        (in)  : length of sparse index sparse_idx             !
    !    sparse_idx     (in)  : index of the matrix being calculated          !
    !    local_start    (in)  : start of col functions in this batch          !
    !    local_start    (in)  : end of col functions in this batch            !
    !    scheme         (in)  : string identifying which scheme to use        !
    !                           'FULL' = calculate all elements               !
    !                           'LOWER' = calculate lower triangle only       !
    !                           'ALTERNATE' = calculate alternating elements  !
    !                           from upper and lower triangles                !
    !    reqs           (out) : array of logical flags to store whether this  !
    !                           node will need to request anything from each  !
    !                           of the other nodes                            !
    !-------------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine in May 2008 based on code  !
    ! from the old integrals_brappd_ketfftbox by Chris-Kriton Skylaris        !
    !=========================================================================!

    use comms, only: comms_abort, pub_my_node_id, pub_on_root, &
         pub_total_num_nodes
    use constants, only: stdout
    use parallel_strategy, only: pub_first_atom_on_node
    use sparse, only: SPAM3, sparse_first_elem_on_node, sparse_atom_of_elem, &
         pattern_full, pattern_alternate, pattern_lower
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! Arguments
    integer, intent(inout) :: num_plan_steps
    integer, intent(out) :: my_plan(2,num_plan_steps)
    integer, intent(in) :: idx_len
    integer, intent(in) :: sparse_idx(idx_len)
    type(SPAM3),intent(in) :: mat
    integer, intent(in) :: local_start, local_end
    character(*), intent(in) :: scheme
    logical, intent(out) :: reqs(0:pub_total_num_nodes-1)

    ! Locals
    integer :: idx                   ! sparse matrix index counter
    integer :: local_col_atom        ! atom of local_col
    integer :: loc_local_col_atom    ! local index of above atom
    integer :: pattern               ! pattern determining elements calculated
    integer :: step                  ! node-block loop counter
    integer :: plan_step             ! position counter for plan array
    integer :: recv_node             ! node to receive data from
    integer :: recv_row              ! the function being received
    integer :: recv_row_atom         ! the atom of the function being received
    integer :: col                   ! global index of column function
    integer :: local_col             ! local index of column function
    integer :: first_col             ! first column on this node
    integer :: first_row             ! first row on recvnode
    integer :: last_row              ! last row on recvnode

    ! ndmh: initialise plan variables
    my_plan = 0
    plan_step = 1
    reqs = .false.

    ! ndmh: Set calculation pattern
    select case (scheme)
    case ('FULL','full','ASYM','asym')
       pattern = pattern_full
    case ('LOWER','lower')
       pattern = pattern_lower
    case ('ALT','alt','ALTERNATE','alternate')
       pattern = pattern_alternate
    case default
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in function_basis_batch_row_plan: calculation pattern "',&
            trim(scheme), '" not recognised'
       call comms_abort
    end select

    ! ndmh: Loop over segments of the matrix, starting from
    ! ndmh: the diagonal segments
    first_col = sparse_first_elem_on_node(pub_my_node_id,mat,'C')
    do step=0,pub_total_num_nodes-1

       ! ndmh: find node to check overlaps with on this step
       recv_node = modulo(pub_my_node_id + step, pub_total_num_nodes)

       ! Loop over the rows on recv_node
       first_row = sparse_first_elem_on_node(recv_node,mat,'R')
       last_row = sparse_first_elem_on_node(recv_node+1,mat,'R') - 1
       do recv_row=first_row,last_row
          recv_row_atom = sparse_atom_of_elem(recv_row,mat,'R')

          ! ndmh: loop over cols on this node to see whether they
          ! ndmh: overlap recv_row
          do local_col=local_start,local_end
             col = local_col + first_col - 1
             local_col_atom = sparse_atom_of_elem(col,mat,'C')
             loc_local_col_atom = local_col_atom - &
                  pub_first_atom_on_node(pub_my_node_id) + 1
             do idx=sparse_idx(loc_local_col_atom), &
                  sparse_idx(loc_local_col_atom+1)-1

                ! ndmh: add to plan only if sparse element exists...
                if (sparse_idx(idx) == recv_row_atom) then

                   ! ndmh: ... and is in the pattern being used
                   if (internal_in_pattern(recv_row,col)) then
                      my_plan(1,plan_step) = recv_row
                      my_plan(2,plan_step) = col
                      reqs(recv_node) = .true.
                      plan_step = plan_step + 1
                   end if
                end if

             end do
          end do
       end do

    end do

    ! ndmh: return index of last entry in plan on this node
    num_plan_steps = plan_step - 1

  contains

    !=========================================================================!
    ! Checks if this node needs to receive and use the data for this col, row !
    ! combination from another node, based on the variable pattern.           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! row (input) : row function to consider receiving in full matrix         !
    ! col (input) : column function to consider receiving in full matrix      !
    !-------------------------------------------------------------------------!
    ! Written by Nick Hine, April 2008                                        !
    !=========================================================================!

    logical function internal_in_pattern(row,col)

      integer, intent(in) :: col
      integer, intent(in) :: row

      select case (pattern)
      case (pattern_full)
         internal_in_pattern = .true.
      case (pattern_lower)
         if (col<=row) then
            internal_in_pattern = .true.
         else
            internal_in_pattern = .false.
         end if
      case (pattern_alternate)
         if (((col<row).and.(mod(col+row,2)==1)).or. &
              ((col>row).and.(mod(col+row,2)==0)).or. &
              (col==row)) then
            internal_in_pattern = .true.
         else
            internal_in_pattern = .false.
         end if
      case default
         internal_in_pattern = .false.
         if (pub_on_root) then
            write(stdout,'(a)')'Invalid pattern supplied to &
                 &internal_in_pattern (function_basis_batch_row_plan)'
            call comms_abort
         end if
      end select

    end function internal_in_pattern

  end subroutine function_basis_batch_row_plan


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_sum_fftbox_batch(fftbox_batch, &
       funcs_on_grid, row_basis, col_basis, batch_size, local_start, &
       local_end, mat_idx, idx_len, coeff_mat, mmat, nmat, &
       row_on_grid_buffer, common_fac)

    !=========================================================================!
    ! This subroutine calculates sums in fftboxes of functions multiplied by  !
    ! coefficients in the form of matrices. The fftbox data is A_\alpha(r)    !
    ! corresponding to a given column function g_\alpha (where the set of     !
    ! column functions is not necessarily the same as the row functions).     !
    ! The fftbox data A_\alpha(r) is constructed from the sum of functions    !
    ! f_\beta multiplied by coefficients M_\beta\alpha, for functions f_\beta !
    ! which overlap \alpha (ie have a nonzero matrix element S_{\beta\alpha}  !
    ! The sum is given by                                                     !
    ! A_\alpha(r) = \sum_beta M_{\beta\alpha} f_\beta(r)                      !
    ! The batch contains a subset of the functions f_\alpha that belong to    !
    ! the node pub_my_node_id. Increasing the size of the batch               !
    ! increases parallel efficiency by reducing communication but also        !
    ! increases the memory use per processor.                                 !
    ! Both A_\alpha and M may be arrays, running from mmat to nmat - the sum  !
    ! will be taken over this range of arrays.                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! row_fftbox_batch (input/output) : FFTboxes containing results A(r)      !
    ! funcs_on_grid (input) : funcs of pub_my_node_id in ppd representation   !
    ! row_basis (input) : Function basis type for the functions being summed  !
    ! col_basis (input) : Function basis type for the matrix columns          !
    ! batch_size (input) : Number of fftboxes (phi_a functions) in each batch !
    ! local_start (input): Number of the first function of the current batch  !
    !    in the counting scheme of all functions of pub_my_node_id.           !
    ! local_end (input)  : Number of the last function of the current batch   !
    !    in the counting scheme of all functions of pub_my_node_id.           !
    ! mat_idx (input) : Overlap matrix index between row funcs and col funcs  !
    ! idx_len (input) : Length of overlap matrix index                        !
    ! coeff_mat (input) : Coefficients matrix M_{\alpha\beta} in SPAM3 format.!
    ! mmat (input) : First entry in the array of A's and M's to use.          !
    ! nmat (input) : Last entry in the array of A's and M's to use.           !
    !-------------------------------------------------------------------------!
    ! This subroutine was written by Chris-Kriton Skylaris on 18/9/2003       !
    ! and is capable of running on parallel computers with an                 !
    ! arbitary number of processors.                                          !
    ! Rewritten for SPAM3, moved to function_basis and made able to support   !
    ! different function sets by Nicholas Hine, April-May 2009.               !
    ! Modified to remove integer buffers and consolidate MPI usage            !
    ! by Nicholas Hine, June 2009.                                            !
    ! Modified to allow different function basis for rows and colums by       !
    ! David O'Regan, September 2009.                                          !
    !=========================================================================!

    use basis, only: basis_location_func_wrt_cell, &
         basis_ket_start_wrt_fftbox, basis_location_fb_wrt_box, &
         basis_add_function_to_box
    use comms, only: comms_barrier, comms_free, pub_my_node_id, &
         pub_total_num_nodes
    use constants, only: DP, stdout
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_get_element, sparse_node_of_elem, &
         sparse_first_elem_on_node, sparse_num_elems_on_node, &
         sparse_node_num_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size
    integer, intent(in) :: nmat
    real(kind=DP), intent(inout) :: fftbox_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3, pub_cell%num_spins, &
         nmat, batch_size)
    integer, intent(in) :: local_start, local_end
    integer, intent(in) :: mmat ! ddor: Starting index
    integer, intent(in) :: idx_len
    integer, intent(in) :: mat_idx(idx_len)
    type(FUNC_BASIS), intent(in) :: row_basis
    type(FUNC_BASIS), intent(in) :: col_basis
    type(SPAM3), intent(in) :: coeff_mat(pub_cell%num_spins,nmat)
    real(kind=DP), intent(in)  :: funcs_on_grid(row_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP) :: row_on_grid_buffer(row_basis%func_on_grid_buffer_size)
    real(kind=DP), intent(in) :: common_fac

    ! Local Variables
    integer :: local_row, local_col
    integer :: col
    integer :: recv_row
    integer :: recv_node
    integer :: batch_count
    integer :: col_start1, col_start2, col_start3
    integer :: row_start1, row_start2, row_start3
    integer :: col_cell_start1, col_cell_start2, col_cell_start3
    integer :: row_cell_start1, row_cell_start2, row_cell_start3
    integer :: is
    integer :: imat
    integer :: ierr             ! Error flag
    logical :: remote_row
    integer :: n_row_ppds
    real(kind=DP) :: coeff(pub_cell%num_spins,nmat)
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
    integer, allocatable :: plan(:,:)         ! Plan for this node
    integer, allocatable :: func_requests(:)
    integer, allocatable :: request_handles(:)

    ! ndmh: allocate workspace
    plan_steps = sparse_node_num_element(coeff_mat(1,mmat))
    allocate(plan(2,plan_steps),stat=ierr)
    call utils_alloc_check('function_basis_sum_fftbox_batch','plan',ierr)
    allocate(func_requests(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('function_basis_sum_fftbox_batch','func_requests',ierr)
    allocate(request_handles(0:pub_total_num_nodes+2),stat=ierr)
    call utils_alloc_check('function_basis_sum_fftbox_batch','request_handles',ierr)
    allocate(reqs_in(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('function_basis_sum_fftbox_batch','reqs_in',ierr)
    allocate(reqs_out(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('function_basis_sum_fftbox_batch','reqs_out',ierr)

    ! ddor: changed from coeff_mat(1,1) to coeff_mat(1,mmat) throughout

    ! ndmh: create a plan from the matrix index
    call function_basis_batch_row_plan(plan_steps,plan,idx_len,mat_idx, &
         coeff_mat(1,mmat),local_start,local_end,'FULL',reqs_in)

    ! ndmh: initializations
    prev_recv_row = -1
    prev_req_row = -1
    local_row = -1
    remote_row = .false.

    ! ndmh: initialisation of send request receive operations
    call function_basis_init_requests(func_requests,request_handles,reqs_in,reqs_out)

    ! cks: col funtions are by definition placed in the centre of the
    ! cks: fftbox and stay there, or they are left where they are in the
    ! cks: simulation cell (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(col_start1, col_start2, col_start3, &
         pub_fftbox%total_pt1, pub_fftbox%total_pt2, pub_fftbox%total_pt3)

    ! ndmh: loop over the steps of the plan
    do plan_step=1-lookahead,plan_steps

       ! ndmh: check if we need to request a function to be sent to this node
       if (plan_step+lookahead<=plan_steps) then
          req_row = plan(1,plan_step+lookahead)
          req_node = sparse_node_of_elem(req_row,coeff_mat(1,mmat),'R')
          ! ndmh: if this function is not local to this node and we do not already
          ! ndmh: have it, then send a request for it to req_node
          if ((req_node /= pub_my_node_id) .and. (req_row /= prev_req_row)) then
             call function_basis_request(req_node,req_row)
             prev_req_row = req_row
          end if
       end if

       ! ndmh: respond to any send requests made of this node
       call function_basis_respond_to_reqs(func_requests,request_handles, &
            row_basis,funcs_on_grid)

       ! ndmh: cycle if still pre-start
       if (plan_step < 1) cycle

       ! ndmh: find the row and col of the matrix element to calculate
       recv_row = plan(1,plan_step)
       col = plan(2,plan_step)

       ! ndmh: receive or copy the ppd list and offset into the buffer sphere
       if (prev_recv_row /= recv_row) then

          ! ndmh: find the node the row function is stored on
          recv_node = sparse_node_of_elem(recv_row,coeff_mat(1,mmat),'R')

          if (recv_node==pub_my_node_id) then
             local_row = recv_row - sparse_first_elem_on_node( &
                  pub_my_node_id,coeff_mat(1,mmat),'R') + 1
             n_row_ppds = row_basis%spheres(local_row)%n_ppds_sphere
             pub_buffer_sphere%n_ppds_sphere = n_row_ppds
             pub_buffer_sphere%ppd_list(:,1:n_row_ppds) = &
                  row_basis%spheres(local_row)%ppd_list(:,1:n_row_ppds)
             pub_buffer_sphere%offset = row_basis%spheres(local_row)%offset
             remote_row = .false.
          else
             call function_basis_recv(recv_node,recv_row,func_requests, &
                  request_handles,pub_buffer_sphere,row_on_grid_buffer, &
                  row_basis,funcs_on_grid)
             remote_row = .true.
          end if
          prev_recv_row = recv_row

          call basis_location_func_wrt_cell(row_cell_start1, &
               row_cell_start2, row_cell_start3, row_basis%all_tbs(recv_row))
       end if

       ! ndmh: find coefficients multiplying this function for each fftbox
       do imat=mmat,nmat
          do is=1,pub_cell%num_spins
             call sparse_get_element(coeff(is,imat),coeff_mat(is,imat), &
                  recv_row,col)
             coeff(is,imat) = coeff(is,imat) * common_fac
          end do
       end do

       ! ndmh: find local column and position of column within this batch
       local_col = col - sparse_first_elem_on_node(pub_my_node_id, &
            coeff_mat(1,mmat),'C') + 1
       batch_count = local_col - local_start + 1

       ! ndmh: Find position of column function wrt simulation cell
       call basis_location_func_wrt_cell(col_cell_start1, &
            col_cell_start2, col_cell_start3, col_basis%all_tbs(col))

       ! ndmh: find where to deposit row function in box
       call basis_location_fb_wrt_box(row_start1,row_start2,row_start3, & ! out
            col_start1,col_start2,col_start3, &                           ! in
            col_cell_start1,col_cell_start2,col_cell_start3, &            ! in
            row_cell_start1,row_cell_start2,row_cell_start3, &            ! in
            pub_cell%total_pt1,pub_cell%total_pt2,pub_cell%total_pt3)     ! in

       ! ndmh: add function times coeff into the fftboxes of the batch
       do imat=mmat,nmat
          do is=1,pub_cell%num_spins
             if (remote_row) then
                call basis_add_function_to_box( &
                     fftbox_batch(:,:,:,is,imat,batch_count), &
                     pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                     pub_fftbox%total_pt3, &
                     row_start1,row_start2,row_start3, &
                     row_basis%all_tbs(recv_row),row_on_grid_buffer, &
                     pub_buffer_sphere,coeff(is,imat))
             else
                call basis_add_function_to_box( &
                     fftbox_batch(:,:,:,is,imat,batch_count), &
                     pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                     pub_fftbox%total_pt3, &
                     row_start1, row_start2,row_start3, &
                     row_basis%all_tbs(recv_row),funcs_on_grid, &
                     row_basis%spheres(local_row),coeff(is,imat))
             end if
          end do  ! is
       end do  ! imat

    end do  ! plan_step

    ! ndmh: send signal indicating completion to all nodes
    do req_node=0,pub_total_num_nodes-1
       if (reqs_in(req_node)) then
          call function_basis_request(req_node,-2000)
       end if
    end do

    call function_basis_await_requests(func_requests,request_handles, &
         row_basis,funcs_on_grid)

    ! ndmh: synchronise again so that non-blocking sends are completed and all
    ! ndmh: handles restored before deallocation
    call comms_free
    call comms_barrier

    ! ndmh: deallocate local workspace
    deallocate(reqs_out,stat=ierr)
    call utils_dealloc_check('function_basis_sum_fftbox_batch','reqs_out',ierr)
    deallocate(reqs_in,stat=ierr)
    call utils_dealloc_check('function_basis_sum_fftbox_batch','reqs_in',ierr)
    deallocate(request_handles,stat=ierr)
    call utils_dealloc_check('function_basis_sum_fftbox_batch','request_handles',ierr)
    deallocate(func_requests,stat=ierr)
    call utils_dealloc_check('function_basis_sum_fftbox_batch','func_requests',ierr)
    deallocate(plan,stat=ierr)
    call utils_dealloc_check('function_basis_sum_fftbox_batch','plan',ierr)

  end subroutine function_basis_sum_fftbox_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_sum_ppd_funcs(sum_on_grid, &        ! inout
       sum_basis, coeff_mat, mmat, nmat, index_mat, &           ! input
       funcs_on_grid, fun_basis)                                ! input

    !=========================================================================!
    ! This subroutine calculates sums in ppds of functions in ppds multiplied !
    ! by coefficients in the form of SPAM3 matrices, where only the functions !
    ! for which there is a nonzero element of an overlap matrix are included  !
    ! in the summation.                                                       !
    ! The sum is given by                                                     !
    ! |sum_a> = \sum_b M_ab |f_b>     for all b where S_ab /= 0               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! sum_on_grid (inout) : sum funcs of pub_my_node_id in ppd representation !
    ! funcs_on_grid (input) : funcs of pub_my_node_id in ppd representation   !
    ! fbasis (input) : Function basis type for the functions being summed     !
    ! coeff_mat (input) : coefficients matrix M_{\alpha\beta} in SPAM3 format.!
    ! index_mat (input) : index matrix S_{\alpha\beta} in SPAM3 format.       !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in May 2010, reusing bits of the routine       !
    ! integrals_brappd_ketppd (originally integrals_brackets).                !
    !=========================================================================!

    use simulation_cell, only: pub_cell
    use comms, only: comms_barrier, comms_free, comms_reduce, pub_my_node_id, &
         pub_total_num_nodes
    use sparse, only: SPAM3, sparse_get_element, sparse_node_of_elem, &
         sparse_first_elem_on_node, sparse_num_elems_on_node, &
         sparse_node_num_element, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: mmat, nmat
    type(FUNC_BASIS), intent(in) :: fun_basis
    type(FUNC_BASIS), intent(in) :: sum_basis
    real(kind=DP), intent(inout) :: sum_on_grid(sum_basis%n_ppds * &
         pub_cell%n_pts,nmat)
    type(SPAM3), intent(in) :: coeff_mat(nmat)
    type(SPAM3), intent(in) :: index_mat
    real(kind=DP), intent(in) :: funcs_on_grid(fun_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    integer :: col, recv_row
    integer :: local_col, local_row
    integer :: recv_node
    integer :: n_ppds
    logical :: remote_row                ! pdh: flag to avoid pointers
    real(kind=DP) :: coeff_mat_el        ! matrix element
    integer :: ierr                      ! ndmh: error flag
    integer :: idx_len
    integer :: imat
    integer,allocatable :: index_mat_idx(:)
    real(kind=DP),allocatable :: func_on_grid_buffer(:)
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

    call timer_clock('function_basis_sum_ppd_funcs',1)

    ! cks: synchronize PEs
    call comms_barrier

    ! ndmh: allocate workspace
    plan_steps = sparse_node_num_element(index_mat)
    allocate(plan(2,plan_steps),stat=ierr)
    call utils_alloc_check('function_basis_sum_ppd_funcs','plan',ierr)
    idx_len = sparse_index_length(index_mat)
    allocate(index_mat_idx(idx_len),stat=ierr)
    call utils_alloc_check('function_basis_sum_ppd_funcs','index_mat_idx',ierr)
    allocate(func_requests(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('function_basis_sum_ppd_funcs','func_requests',ierr)
    allocate(request_handles(0:pub_total_num_nodes+2),stat=ierr)
    call utils_alloc_check('function_basis_sum_ppd_funcs','request_handles',ierr)
    allocate(reqs_in(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('function_basis_sum_ppd_funcs','reqs_in',ierr)
    allocate(reqs_out(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('function_basis_sum_ppd_funcs','reqs_out',ierr)
    allocate(func_on_grid_buffer(fun_basis%max_n_ppds_sphere*pub_cell%n_pts), &
         stat=ierr)
    call utils_alloc_check('function_basis_sum_ppd_funcs', &
         'func_on_grid_buffer',ierr)

    ! ndmh: get the matrix index
    call sparse_generate_index(index_mat_idx,index_mat)

    ! ndmh: create a plan from the matrix index
    local_col = sparse_num_elems_on_node(pub_my_node_id,coeff_mat(mmat),'C')
    call function_basis_batch_row_plan(plan_steps,plan,idx_len,index_mat_idx, &
         coeff_mat(mmat),1,local_col,'FULL',reqs_in)

    ! ndmh: initializations
    prev_recv_row = -1
    prev_req_row = -1
    local_row = -1
    remote_row = .false.

    ! ndmh: initialisation of send request receive operations
    call function_basis_init_requests(func_requests,request_handles,reqs_in,reqs_out)

    do plan_step=1,plan_steps

       ! ndmh: check if we need to request a function to be sent to this node
       if (plan_step+lookahead<=plan_steps) then
          req_row = plan(1,plan_step+lookahead)
          req_node = sparse_node_of_elem(req_row,coeff_mat(mmat),'R')
          ! ndmh: if this bra is not local to this node and we do not already
          ! ndmh: have it, then send a request for it to req_node
          if ((req_node /= pub_my_node_id) .and. (req_row /= prev_req_row)) then
             call function_basis_request(req_node,req_row)
             prev_req_row = req_row
          end if
       end if

       ! ndmh: respond to any send requests made of this node
       call function_basis_respond_to_reqs(func_requests,request_handles, &
            fun_basis,funcs_on_grid)

       ! ndmh: find the row and col of the matrix element to calculate
       recv_row = plan(1,plan_step)
       col = plan(2,plan_step)

       ! ndmh: find the node the row function is stored on
       recv_node = sparse_node_of_elem(recv_row,coeff_mat(mmat),'R')

       ! ndmh: receive or copy the ppd list and offset into the buffer sphere
       if (prev_recv_row /= recv_row) then
          if (recv_node==pub_my_node_id) then
             local_row = recv_row - sparse_first_elem_on_node(pub_my_node_id, &
                  coeff_mat(mmat),'R') + 1
             n_ppds = fun_basis%spheres(local_row)%n_ppds_sphere
             pub_buffer_sphere%n_ppds_sphere = &
                  fun_basis%spheres(local_row)%n_ppds_sphere
             pub_buffer_sphere%ppd_list(:,1:n_ppds) = &
                  fun_basis%spheres(local_row)%ppd_list(:,1:n_ppds)
             pub_buffer_sphere%offset = fun_basis%spheres(local_row)%offset
             remote_row = .false.
          else
             call function_basis_recv(recv_node,recv_row,func_requests, &
                  request_handles,pub_buffer_sphere,func_on_grid_buffer, &
                  fun_basis,funcs_on_grid)
             remote_row = .true.
          end if
          prev_recv_row = recv_row
       end if

       ! Find local index of col on this node
       local_col = col - &
            sparse_first_elem_on_node(pub_my_node_id,coeff_mat(mmat),'C') + 1


       do imat=mmat,nmat

          ! ndmh: get matrix element from SPAM3
          call sparse_get_element(coeff_mat_el,coeff_mat(imat), &
               recv_row,col)
          ! ndmh: add coeff_mat matrix element times function to result
          if (remote_row) then
             call internal_add_ppds(sum_on_grid(:,imat), &
                  func_on_grid_buffer, sum_basis%n_ppds, &
                  pub_buffer_sphere%n_ppds_sphere,sum_basis%spheres(local_col),&
                  pub_buffer_sphere, coeff_mat_el)
          else
             call internal_add_ppds(sum_on_grid(:,imat), &
                  funcs_on_grid, sum_basis%n_ppds, fun_basis%n_ppds, &
                  sum_basis%spheres(local_col),fun_basis%spheres(local_row), &
                  coeff_mat_el)
          end if

       end do

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
         fun_basis,funcs_on_grid)

    ! ndmh: synchronise again so that non-blocking sends are completed and all
    ! ndmh: handles restored before deallocation
    call comms_free
    call comms_barrier

    ! ndmh: deallocate workspace
    deallocate(func_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('function_basis_sum_ppd_funcs','func_on_grid_buffer', &
         ierr)
    deallocate(reqs_out,stat=ierr)
    call utils_dealloc_check('function_basis_sum_ppd_funcs','reqs_out',ierr)
    deallocate(reqs_in,stat=ierr)
    call utils_dealloc_check('function_basis_sum_ppd_funcs','reqs_in',ierr)
    deallocate(request_handles,stat=ierr)
    call utils_dealloc_check('function_basis_sum_ppd_funcs','request_handles',ierr)
    deallocate(func_requests,stat=ierr)
    call utils_dealloc_check('function_basis_sum_ppd_funcs','func_requests',ierr)
    deallocate(index_mat_idx,stat=ierr)
    call utils_dealloc_check('function_basis_sum_ppd_funcs','index_mat_idx',ierr)
    deallocate(plan,stat=ierr)
    call utils_dealloc_check('function_basis_sum_ppd_funcs','plan',ierr)

    call timer_clock('function_basis_sum_ppd_funcs',2)

  contains

    ! ndmh: copy of integrals_bra_dot_ket_ppds, turned around logically
    subroutine internal_add_ppds(sumfunc, func, n_sum_ppds, n_func_ppds, &
         sum_sphere, func_sphere, coeff)

      implicit none

      ! Arguments
      integer, intent(in) :: n_sum_ppds, n_func_ppds
      real(kind=DP), intent(inout) :: sumfunc(n_sum_ppds*pub_cell%n_pts)
      real(kind=DP), intent(in) :: func(n_func_ppds*pub_cell%n_pts)
      type(SPHERE), intent(in)  :: sum_sphere
      type(SPHERE), intent(in)  :: func_sphere
      real(kind=DP), intent(in) :: coeff

      ! Local Variables
      integer :: isum, ifunc
      integer :: sum_ppd, func_ppd
      integer :: sum_start, func_start
      integer :: i                       ! loop counter

      ! Start at the beginning of the sum sphere ppd list
      ifunc = 1
      func_ppd = func_sphere%ppd_list(1,ifunc)

      ! ndmh: calculate sum = sum + alpha*func
      do isum=1,sum_sphere%n_ppds_sphere
         ! ndmh: find ppd number and start position for isum
         sum_ppd = sum_sphere%ppd_list(1,isum)
         sum_start = sum_sphere%offset + (isum-1)*pub_cell%n_pts
         do
            ! ndmh: keep moving on while func_ppd is less than sum_ppd
            if (func_ppd < sum_ppd) then
               ifunc = ifunc + 1
               if (ifunc > func_sphere%n_ppds_sphere) exit
               func_ppd = func_sphere%ppd_list(1,ifunc)
               cycle
            end if

            ! ndmh: ppd numbers match, so add this ppd of func to sum
            if (func_ppd == sum_ppd) then
               func_start = func_sphere%offset + (ifunc-1)*pub_cell%n_pts
               do i=0,pub_cell%n_pts-1
                  sumfunc(sum_start+i) = sumfunc(sum_start+i) + &
                       coeff*func(func_start+i)
               end do
            end if

            ! ndmh: move on to next sum ppd as soon as we have found a match
            ! ndmh: or moved beyond sum_ppd
            if (func_ppd >= sum_ppd) exit

         end do
      end do

    end subroutine internal_add_ppds

  end subroutine function_basis_sum_ppd_funcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_init_requests(requests,handles,reqs_in,reqs_out)

    use comms, only: comms_alltoall, comms_irecv, pub_my_node_id, &
         pub_null_handle, pub_total_num_nodes

    implicit none

    ! Arguments
    integer,intent(out) :: requests(0:pub_total_num_nodes-1)
    integer,intent(out) :: handles(0:pub_total_num_nodes+2)
    logical,intent(inout) :: reqs_in(0:pub_total_num_nodes-1)
    logical,intent(out) :: reqs_out(0:pub_total_num_nodes-1)

    ! Locals
    integer :: node

    ! ndmh: initialisation
    probe_count = 0
    reqs_in(pub_my_node_id) = .false.
    reqs_out = reqs_in
    call comms_alltoall(reqs_out,1)

    ! ndmh: loop over all nodes
    do node=0,pub_total_num_nodes-1
       if (reqs_out(node)) then
          ! ndmh: start asynchronous receive for request from this node
          requests(node) = FUNCS_WAITING
          call comms_irecv(node,requests(node),1,tag=req_tag+node, &
               handle=handles(node))
       else
          requests(node) = FUNCS_DONE
          handles(node) = pub_null_handle
       end if
    end do

    req_buffer_index = 1

  end subroutine function_basis_init_requests


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_request(req_node,req_row)

    use comms, only: comms_send, pub_my_node_id

    implicit none

    ! Arguments
    integer,intent(in) :: req_node
    integer,intent(in) :: req_row

    ! Store requested row index in buffer to asynchronously send
    req_buffer(req_buffer_index) = req_row

    ! ndmh: send request for function required
    call comms_send(req_node,req_buffer(req_buffer_index),1, &
         tag=req_tag+pub_my_node_id)

    ! ndmh: start asynchronous receive operations for the incoming ngwf
    ! ndmh: kept in case multiple buffers get implemented
    !call comms_irecv(recv_node,buffer_sphere%ppd_list(:,:),2*n_ppds, &
    !     tag=recv_row*2+0,handle=handles(pub_total_num_nodes+0))
    !call comms_irecv(recv_node,func_on_grid_buffer,n_pts, &
    !     tag=recv_row*2+1,handle=handles(pub_total_num_nodes+1))

    ! ndmh: advance buffer index and cycle if filled
    req_buffer_index = req_buffer_index + 1
    if (req_buffer_index > req_buffer_size) req_buffer_index = 1

  end subroutine function_basis_request


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_await_requests(requests,handles,fbasis, &
       funcs_on_grid)

    use comms, only: comms_waitany, comms_wait, pub_total_num_nodes
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    integer,intent(inout)      :: requests(0:pub_total_num_nodes-1)
    integer,intent(inout)      :: handles(0:pub_total_num_nodes+2)
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(in)  :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)

    ! Locals
    integer :: node

    do
       if (all(requests(:)==FUNCS_DONE)) exit

       ! ndmh: respond to any outstanding requests
       call function_basis_respond_to_reqs(requests,handles,fbasis, &
            funcs_on_grid)

       ! ndmh: if all other nodes have not yet finished, wait for a new request
       if (any(requests(:)==FUNCS_WAITING)) then
          call comms_waitany(pub_total_num_nodes,handles)
       end if
    end do

    ! ndmh: ensure all nodes have completed their final receive operation
    do node=0,pub_total_num_nodes-1
       call comms_wait(handles(node))
    end do

  end subroutine function_basis_await_requests


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_respond_to_reqs(requests,handles,fbasis, &
       funcs_on_grid)

    use basis, only: SPHERE
    use comms, only: comms_irecv, comms_wait, comms_probe, pub_my_node_id, &
         pub_total_num_nodes
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    integer,intent(inout) :: requests(0:pub_total_num_nodes-1)
    integer,intent(inout) :: handles(0:pub_total_num_nodes+3)
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(in)  :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)

    ! Locals
    integer :: node
    logical :: result
    logical :: can_exit

    ! Every probe_frequency times this routine is called, we must enter
    ! the MPI library to ensure any ongoing sends complete even if no
    ! requests are received.
    probe_count = probe_count + 1
    if (probe_count >= probe_frequency) then
       call comms_probe(result,pub_my_node_id)
       probe_count = 0
    end if

    ! ndmh: do not continue until there are no outstanding requests
    can_exit = .true.
    do
       ! ndmh: loop over nodes checking for requests
       do node=0,pub_total_num_nodes-1

          if (requests(node)<0) cycle
          can_exit = .false.

          ! ndmh: ensure receive has completed
          call comms_wait(handles(node))

          ! ndmh: send requested ngwf to node
          call function_basis_send(node,requests(node),fbasis,funcs_on_grid)

          ! ndmh: re-initialise the request receive operation
          requests(node) = FUNCS_WAITING
          call comms_irecv(node,requests(node),1,tag=req_tag+node, &
               handle=handles(node))
          call comms_probe(result,node)

       end do
       if (can_exit) then
          exit
       else
          can_exit = .true.
       end if
    end do

  end subroutine function_basis_respond_to_reqs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_send(send_node,send_row,fbasis,funcs_on_grid)

    use basis, only: SPHERE
    use comms, only: comms_send, pub_my_node_id
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    integer,intent(in) :: send_node
    integer,intent(in) :: send_row
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(in)  :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)

    ! Locals
    integer :: local_row
    integer :: n_ppds,n_pts,start_pt

    local_row = send_row - fbasis%first_on_node(pub_my_node_id) + 1
    n_ppds = fbasis%n_ppds_sphere(send_row)
    n_pts = n_ppds*pub_cell%n_pts
    start_pt = fbasis%spheres(local_row)%offset

    ! ndmh: send sphere information and ppd data for this function
    call comms_send(send_node,fbasis%spheres(local_row)%ppd_list(:,:), &
         2*n_ppds,tag=send_row*2+0)
    call comms_send(send_node,funcs_on_grid(start_pt:start_pt+n_pts-1), &
         n_pts,tag=send_row*2+1)

  end subroutine function_basis_send


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine function_basis_recv(recv_node,recv_row,requests,handles, &
       buffer_sphere,func_on_grid_buffer,fbasis,funcs_on_grid)

    use basis, only: SPHERE
    use comms, only: comms_irecv, comms_waitany, comms_wait, &
         pub_null_handle, pub_total_num_nodes
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    integer,intent(in) :: recv_node
    integer,intent(in) :: recv_row
    integer,intent(inout) :: requests(0:pub_total_num_nodes-1)
    integer,intent(inout) :: handles(0:pub_total_num_nodes+2)
    type(SPHERE),intent(inout) :: buffer_sphere
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP),intent(out) :: func_on_grid_buffer( &
         fbasis%max_n_ppds_sphere*pub_cell%n_pts)
    real(kind=DP), intent(in)  :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)

    ! Locals
    integer :: n_ppds,n_pts
    integer :: count

    ! ndmh: find sizes of buffers
    n_ppds = fbasis%n_ppds_sphere(recv_row)
    n_pts = n_ppds*pub_cell%n_pts
    buffer_sphere%n_ppds_sphere = n_ppds
    buffer_sphere%offset = 1

    ! ndmh: start asynchronous receive operations for the incoming ngwf
    call comms_irecv(recv_node,buffer_sphere%ppd_list(:,:),2*n_ppds, &
         tag=recv_row*2+0,handle=handles(pub_total_num_nodes+0))
    call comms_irecv(recv_node,func_on_grid_buffer,n_pts, &
         tag=recv_row*2+1,handle=handles(pub_total_num_nodes+1))
    count = 0
    do
       count = count + 1

       if (handles(pub_total_num_nodes+1)==pub_null_handle) then
          ! call comms_wait to ensure receive has completed before exiting
          call comms_wait(handles(pub_total_num_nodes+0))
          call comms_wait(handles(pub_total_num_nodes+1))
          exit
       end if

       if (any(requests(:)>0)) then
          ! ndmh: send any functions required by other nodes, if there are
          ! ndmh: outstanding requests
          call function_basis_respond_to_reqs(requests,handles,fbasis, &
               funcs_on_grid)
       else
          ! ndmh: wait for either an incoming request, or one of the two recv
          ! ndmh: operations to complete
          call comms_waitany(pub_total_num_nodes+2,handles)
       end if

    end do

  end subroutine function_basis_recv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine function_basis_ppds_to_sph_waves(sw_coeffs, maxln, maxn, &
       ifunc,sw_node_id,funcs_on_grid,fbasis)

    !=======================================================================!
    ! Converts a function in PPD representation to a linear combination of  !
    ! spherical waves. The result is sent to the root node for writing on a !
    ! file                                                                  !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Originally written by Alvaro Ruiz Serrano in January 2009 based on    !
    ! previous code by Mark Robinson.                                       !
    ! Revisited by Alvaro Ruiz Serrano in November 2011 for re-structuring  !
    ! and bug-fixing inspired on function_basis_ppds_to_tightbox.           !
    !=======================================================================!

    use basis, only: basis_location_func_wrt_cell, basis_copy_function_to_box,&
         basis_func_centre_wrt_fftbox, SPHERE, FUNCTION_TIGHT_BOX
    use comms, only: comms_recv, comms_send, pub_my_node_id
    use constants, only: DP, stdout
    use geometry, only: POINT
    use ion, only: element
    use rundat, only: pub_write_max_l
    use simulation_cell, only: pub_cell, pub_tb_recip_grid, &
         pub_maxtight_pts1, pub_maxtight_pts2, pub_maxtight_pts3
    use spherical_wave, only: sw_recp_generate_in_tb, sw_bessel_zeros
    use utils, only: utils_alloc_check, utils_dealloc_check
    use timer, only: timer_clock


    implicit none

    ! ars: arguments
    real(kind=DP), intent(inout) :: sw_coeffs(0:pub_write_max_l,-pub_write_max_l:pub_write_max_l,1:maxn)
    integer,       intent(in   ) :: maxln(0:pub_write_max_l)
    integer,       intent(in   ) :: maxn
    integer, intent(in)          :: ifunc
    integer, intent(in)          :: sw_node_id
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(in)    :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)


    ! ars: local variables
    integer :: loc_ifunc, tb_n1, tb_n2, tb_n3, start1, start2, start3
    integer :: ll, mm, nn, xx, yy, zz, ierr
    real(kind=DP) :: radius, qnl, sw_dot_func, sw_dot_sw
    real(kind=DP), allocatable :: sw_tb(:,:,:), func_tb(:,:,:)
    complex(kind=DP), allocatable :: sw_work(:,:,:)
    type(POINT) :: centre, func_centre

    ! ars: initialize swcoeff
    sw_coeffs(:,:,:)=0.0_DP


    ! ndmh: local node has the ppds of this function, so must copy them to its
    ! ndmh: tightbox buffer (and calculate the offset from the origin) and
    ! ndmh: send it (unless it is also the node that is receiving it).
    if (pub_my_node_id==fbasis%node_of_func(ifunc)) then

       ! ars: find local function, its radius, centre and tighbox
       loc_ifunc = ifunc - fbasis%first_on_node(pub_my_node_id) + 1
       centre = fbasis%spheres(loc_ifunc)%centre
       radius = fbasis%spheres(loc_ifunc)%radius
       tb_n1 = pub_maxtight_pts1
       tb_n2 = pub_maxtight_pts2
       tb_n3 = pub_maxtight_pts3

       ! ars: allocate workspace
       allocate(sw_tb(tb_n1,tb_n2,tb_n3), stat=ierr)
       call utils_alloc_check('function_basis_ppds_to_sph_waves','sw_tb',ierr)
       allocate(sw_work(tb_n1,tb_n2,tb_n3), stat=ierr)
       call utils_alloc_check('function_basis_ppds_to_sph_waves','sw_work',ierr)
       allocate(func_tb(tb_n1,tb_n2,tb_n3), stat=ierr)
       call utils_alloc_check('function_basis_ppds_to_sph_waves','func_tb',ierr)


       ! ars : put functions in tightbox
       call basis_copy_function_to_box(func_tb, tb_n1,tb_n2,tb_n3,&
            1, 1, 1, fbasis%tight_boxes(loc_ifunc), funcs_on_grid, fbasis%spheres(loc_ifunc))

       ! ars : find the centre of the function wrt the tightbox
       call basis_location_func_wrt_cell(start1,start2,start3, fbasis%tight_boxes(loc_ifunc))

       func_centre=basis_func_centre_wrt_fftbox(centre, &
            1,1,1,& ! ars : position of the tightbox wrt FFTbox
            start1, start2, start3) ! ars : position of the FFTbox wrt cell


       ! loop over angular momentum l
       do ll=0,pub_write_max_l
          ! loop over azimuthal angular momentum m
          do mm=-ll,+ll
             ! loop over n
             do nn=1,maxln(ll)

                ! ars : calculate q_nl
                qnl = sw_bessel_zeros(nn,ll)/radius

                ! ars : generate the SW in the tightbox
                ! ars : the centre of the function and the spherical wave
                ! ars : must coincide in order to calculate the overlap
                call sw_recp_generate_in_tb(sw_tb,sw_work,ll,mm,&
                     qnl,radius, func_centre, tb_n1,tb_n2,tb_n3)


                ! ars: calculate coefficients
                sw_dot_func=0.0_DP
                sw_dot_sw=0.0_DP
                do zz=1, tb_n3
                   do yy=1, tb_n2
                      do xx=1, tb_n1
                         sw_dot_func = sw_dot_func + sw_tb(xx,yy,zz)*func_tb(xx,yy,zz)
                         sw_dot_sw = sw_dot_sw + sw_tb(xx,yy,zz)*sw_tb(xx,yy,zz)
                      enddo
                   enddo
                enddo
                sw_coeffs(ll,mm,nn)=sw_dot_func/sw_dot_sw

             enddo  ! end n loop
          enddo     ! end m loop
       enddo        ! end l loop

       ! ndmh: send, if destination node is not local
       if (.not. (pub_my_node_id==sw_node_id)) then
          call comms_send(sw_node_id,sw_coeffs,tag=ifunc+fbasis%num)
       end if


       ! ars: deallocate workspace
       deallocate(sw_tb,stat=ierr)
       call utils_dealloc_check('function_basis_ppds_to_sph_waves','sw_tb',ierr)
       deallocate(sw_work,stat=ierr)
       call utils_dealloc_check('function_basis_ppds_to_sph_waves','sw_work',ierr)
       deallocate(func_tb,stat=ierr)
       call utils_dealloc_check('function_basis_ppds_to_sph_waves','func_tb',ierr)


       ! ndmh: local node does not have the ppds of this function. Do nothing
       ! ndmh: unless this is the receiving node (in which case, receive the
       ! ndmh: tightbox.
    else

       if (pub_my_node_id==sw_node_id) then
          ! ndmh: receive the function from its node
          call comms_recv(fbasis%node_of_func(ifunc),sw_coeffs,tag=ifunc+fbasis%num)
       endif

    end if

  end subroutine function_basis_ppds_to_sph_waves


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine function_basis_sph_waves_to_ppds(fbasis,funcs_on_grid,sw_coeffs,&
       maxln,maxl,maxn,ifunc,sw_node_id)

    !=======================================================================!
    ! Converts a function in spherical wave representation to PPDs.         !
    ! The result is sent to the node that owns the function for continuing  !
    ! parallel calculation.                                                 !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Originally written by Alvaro Ruiz Serrano in January 2009.            !
    ! Revisited by Alvaro Ruiz Serrano in November 2011 for re-structuring  !
    ! and bug-fixing inspired on function_basis_tightbox_to_ppds.           !
    !=======================================================================!

    use basis, only: basis_location_func_wrt_cell, basis_copy_function_to_box,&
         basis_func_centre_wrt_fftbox, SPHERE, FUNCTION_TIGHT_BOX,&
         basis_put_tightbox_in_fftbox,basis_extract_function_from_box,&
         basis_clean_function
    use comms, only: comms_recv, comms_send, pub_my_node_id
    use constants, only: DP, stdout
    use geometry, only: POINT
    use ion, only: element
    use simulation_cell, only: pub_cell, pub_tb_recip_grid, pub_fftbox,&
         pub_maxtight_pts1, pub_maxtight_pts2, pub_maxtight_pts3
    use spherical_wave, only: sw_recp_generate_in_tb, sw_bessel_zeros
    use utils, only: utils_alloc_check, utils_dealloc_check
    use timer, only: timer_clock

    implicit none

    ! ars: arguments
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind=DP), intent(inout) :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(inout) :: sw_coeffs(0:maxl,-maxl:maxl,1:maxn)
    integer,       intent(in   ) :: maxln(0:maxl)
    integer,       intent(in   ) :: maxl
    integer,       intent(in   ) :: maxn
    integer, intent(in)          :: ifunc
    integer, intent(in)          :: sw_node_id


    ! ars: local variables
    integer :: loc_ifunc, tb_n1, tb_n2, tb_n3, start1, start2, start3
    integer :: ll, mm, nn, ierr
    integer :: offset, lastp, npoints
    real(kind=DP) :: radius, qnl
    ! ars: << SW FFTbox and tightbox >>
    real(kind=DP), allocatable    :: func_fftbox(:,:,:)
    real(kind=DP), allocatable    :: func_tb(:,:,:)
    complex(kind=DP), allocatable :: sw_work(:,:,:)
    real(kind=DP), allocatable    :: sw_tb(:,:,:)
    type(POINT) :: centre, func_centre

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~ Generate NGWFs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !################## Generate NGWFs in reciprocal space #####################!


    ! ars : 1st : The SW are generated in the universal tightbox (sw_tb)
    ! ars : 2nd : Generate the NGWFs also in the universal tightbox (ngwfs_tb)
    ! ars : 3rd : Put the NGWFs in a FFTbox in an arbitrary position. In this
    ! ars :       case, the universal tightbox is placed in the bottom-left
    ! ars :       corner of the FFTbox ([1,1,1]). The FFTbox is only used to
    ! ars :       extract ppds, so the position of the universal tightbox
    ! ars :       wrt to the FFTbox is irrelevant.
    ! ars : 4th : Extract ppds from the FFTbox to generate ngwfs_on_grid


    ! ndmh: local node needs the ppds of this function, so must receive them
    ! ndmh: into its tightbox buffer (and receive the read offset from the
    ! ndmh: origin) and extract the ppds.
    if (pub_my_node_id==fbasis%node_of_func(ifunc)) then

       ! ars: find local function, its radius, centre and tighbox
       loc_ifunc = ifunc - fbasis%first_on_node(pub_my_node_id) + 1
       centre = fbasis%spheres(loc_ifunc)%centre
       radius = fbasis%spheres(loc_ifunc)%radius
       offset = fbasis%spheres(loc_ifunc)%offset
       npoints = fbasis%spheres(loc_ifunc)%n_ppds_sphere * pub_cell%n_pts
       lastp = offset + npoints -1
       tb_n1 = pub_maxtight_pts1
       tb_n2 = pub_maxtight_pts2
       tb_n3 = pub_maxtight_pts3

       ! ndmh: recv tightbox data, and origin of atom in read tightbox,
       ! ndmh: if node holding tb is not local node
       if (.not.(pub_my_node_id==sw_node_id)) then
          call comms_recv(sw_node_id,sw_coeffs,tag=ifunc+fbasis%num)
       end if

       ! ars : allocate workspace
       allocate(func_fftbox(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
            pub_fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('restart_ngwfs_swcoeff_input','ngwfs_fftbox',ierr)
       allocate(func_tb(tb_n1, tb_n2, tb_n3),stat=ierr)
       call utils_alloc_check('restart_ngwfs_swcoeff_input','ngwfs_tb',ierr)
       allocate(sw_work(tb_n1, tb_n2, tb_n3),stat=ierr)
       call utils_alloc_check('restart_ngwfs_swcoeff_input','sw_work',ierr)
       allocate(sw_tb(tb_n1, tb_n2, tb_n3),stat=ierr)
       call utils_alloc_check('restart_ngwfs_swcoeff_input','sw_tb',ierr)

       ! ars : initialisation of boxes
       func_fftbox=0.0_DP
       func_tb=0.0_DP

       ! ars : the SW will be created in the NGWFs position wrt to the
       ! ars : universal tightbox
       call basis_location_func_wrt_cell( &
            start1,start2,start3, fbasis%tight_boxes(loc_ifunc))

       ! ars : vector that points from the origin of the universal tightbox
       ! ars : to the centre of the localization sphere
       func_centre=basis_func_centre_wrt_fftbox( &
            fbasis%spheres(loc_ifunc)%centre,&
            1,1,1,& ! ars : position of the tightbox wrt the fftbox
                    ! ars : (arbitrarily placed in [1,1,1])
            start1, start2, start3) ! ars : position of the FFTbox wrt cell

       ! ars : loop over angular momentum l
       do ll=0,maxl

          ! ars : loop over azimuthal angular momentum m
          do mm=-ll,+ll

             ! ars : loop over n
             do nn=1,maxln(ll)

                qnl=sw_bessel_zeros(nn,ll)/radius

                ! ars : generate the SW in the universal tightbox
                call sw_recp_generate_in_tb(sw_tb,sw_work,ll,mm,qnl,radius,&
                     func_centre,tb_n1, tb_n2, tb_n3)

                ! ars : generate the ngwfs in the universal tightbox
                ! ars : by a linear combination of spherical waves
                func_tb = func_tb + sw_tb*sw_coeffs(ll,mm,nn)

             enddo  ! end n loop
          enddo     ! end m loop
       enddo        ! end l loop


       ! ars : put the NGWFs in the FFTbox
       call basis_put_tightbox_in_fftbox(func_fftbox, &   ! input/output
            1, 1, 1, & ! ars : position of the tightbox wrt the fftbox
            func_tb, tb_n1, tb_n2, tb_n3, 1.0_DP)  ! ars : factor = 1. No scaling.

       ! ars : extract ppds from FFTbox and fill ngwfs_on_grid
       call basis_extract_function_from_box(&
            funcs_on_grid(offset:lastp),pub_fftbox%total_ld1, &
            pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            func_fftbox,fbasis%spheres(loc_ifunc),&
            fbasis%tight_boxes(loc_ifunc),&
            1,1,1, 1)



       ! ars : deallocate workspace
       deallocate(sw_tb,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_swcoeff_input','sw_tb',ierr)
       deallocate(sw_work,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_swcoeff_input','sw_work',ierr)
       deallocate(func_tb,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_swcoeff_input','ngwfs_tb',ierr)
       deallocate(func_fftbox,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_swcoeff_input','ngwfs_fftbox',ierr)

       ! ars : shave functions - localise within sphere
       call basis_clean_function(funcs_on_grid,&
            fbasis%spheres(loc_ifunc), fbasis%n_ppds)


       ! ndmh: local node does not have the ppds of this function. Do nothing
       ! ndmh: unless this is the sending node, in which case, send the
       ! ndmh: tightbox and its position within the tightbox as written.
    else

       if (pub_my_node_id==sw_node_id) then
          ! ndmh: send the read function to its node
          call comms_send(fbasis%node_of_func(ifunc),sw_coeffs,tag=ifunc+fbasis%num)
       endif

    endif

    !#################### End generate NGWFs in reciprocal space ###############!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End generate NGWFs ~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  end subroutine function_basis_sph_waves_to_ppds

end module function_basis
