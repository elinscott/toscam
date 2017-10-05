! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                NGWF DIIS optimisation module                   !
!                                                                !
! This module performs the optimisation of the electronic energy !
! with respect to the NGWFs by Direct Inversion in an Iterative  !
! Subspace (DIIS). At each DIIS step, the NGWFs residual vector  !
! (here taken as the NGWFS gradients) are minimised within the   !
! current subspace.                                              !
!----------------------------------------------------------------!
! This module was created by Simon M.-M. Dubois on 28/10/2010    !
! out of parts of the ngwf_cg_mod.                               !
!----------------------------------------------------------------!
! WARNING : this module is still under development/testing.      !
!================================================================!

module ngwf_diis

  implicit none

  private

  public :: ngwf_diis_optimise

contains

  subroutine ngwf_diis_optimise(total_energy, converged, &
       ham, denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &
       hub, elements, ewald_energy, localpseudo_fine, &
       core_density_fine, lhxc_fine)

    !==========================================================================!
    ! This subroutine minimises the total energy with respect to the NGWF      !
    ! expansion coefficients and the density kernel elements subject to the    !
    ! constraint that the density matrix is idempotent and integrates to       !
    ! the correct number of electrons.                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! total_energy  (output)  : the total energy                               !
    ! converged     (output)  : whether the total energy was converged         !
    ! rep           (inout)   : NGWF Representation (functions and matrices)   !
    ! denskern      (inout)   : density kernel for only alpha (or beta)        !
    !                           electrons                                      !
    ! elements         (input): all elements in the order they appear in       !
    !                           the input file                                 !
    ! ngwf_basis       (input): Function basis type describing the NGWFs       !
    ! proj_basis       (input): Function basis type describing the nonlocal    !
    !                           pseudopotential projectors                     !
    ! hub_proj_basis   (input): Function basis type describing the Hubbard     !
    !                           projectors                                     !
    ! ewald_energy     (input): Ewald energy                                   !
    ! localpseudo_fine (input): local pseudopotential on fine grid             !
    ! core_density_fine(input): core density on fine grid for NLCC             !
    ! lhxc_fine        (input): Local-Hartree-Exhange-Correlation potential    !
    !--------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois in July 2010.                              !
    !==========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: comms_abort, pub_on_root, pub_root_node_id, comms_bcast, &
         comms_reduce, pub_my_node_id, comms_barrier
    use constants, only: DP, verbose, max_spins, normal, stdout, brief, stderr
    use electronic, only: electronic_energy, electronic_lagrangian
    use function_basis, only: FUNC_BASIS, function_basis_est_num_psincs
    use geometry, only: point
    use hamiltonian, only: hamiltonian_build_matrix, &
         hamiltonian_dens_indep_matrices, hamiltonian_energy_components
    use hubbard_build, only: HUBBARD_MODEL, hubbard_energy_info
    use ion, only: element
    use ngwf_gradient, only: ngwf_gradient_lnv, ngwf_gradient_exit
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use restart, only: restart_kernel_write, restart_ngwfs_tightbox_output, &
         restart_sph_waves_output, restart_ngwfs_tightbox_store, &
         restart_kernel_store, store_tightbox_ngwfs, store_denskern, &
         retrieve_tightbox_ngwfs, retrieve_denskern
    use rundat, only: pub_kernel_update, write_tightbox_ngwfs, ngwf_cg_type, &
         pub_output_detail, pub_nnho, precond_real, precond_recip, &
         minit_lnv, ngwf_threshold_orig, maxit_ngwf_cg, maxit_ngwf_diis, &
         pub_usehfx, pub_write_ngwf_plot, pub_cube_format, pub_dx_format, &
         pub_grd_format, pub_hubbard, print_qc, k_zero, lnv_threshold_orig, &
         maxit_lnv, pub_write_sw_ngwfs, task, write_denskern, pub_elec_cg_max, &
         pub_write_converged_dk_ngwfs, pub_paw, ngwf_diis_size, &
         ngwf_diis_max_step, ngwf_diis_reset, pub_devel_code
    use integrals, only : integrals_brappd_ketppd
    use services, only: services_flush, services_polak_cg_coeff, &
         services_cubic_fit_minimum, services_line_search_parabola
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_trace, &
         sparse_copy, sparse_product
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_ngwfs
    use wrappers, only : wrappers_ddot

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(ELEMENT), intent(in)    :: elements(pub_cell%nat)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(SPAM3), intent(inout)   :: denskern(pub_cell%num_spins)
    type(NGWF_REP), intent(inout) :: rep
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(in) :: ewald_energy
    real(kind=DP), intent(in) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(out)   :: total_energy
    logical, intent(out) :: converged  ! pdh: convergence flag
    type(NGWF_HAM), intent(inout) :: ham

    ! Local Variables
    type(SPAM3), allocatable :: pur_denskern(:)

    real(kind=DP), allocatable, dimension(:) :: cov_grad_on_grid
    real(kind=DP), allocatable, dimension(:) :: contra_grad_on_grid
    real(kind=DP), allocatable, dimension(:) :: direction_on_grid
    real(kind=DP), allocatable, dimension(:) :: contra_direction_on_grid
    real(kind=DP), allocatable, dimension(:) :: prev_direction_on_grid
    real(kind=DP), allocatable, dimension(:) :: prev_contra_grad_on_grid
    real(kind=DP), allocatable, dimension(:) :: prev_contra_direction_on_grid
    real(kind=DP), allocatable, dimension(:) :: trial_ngwfs_on_grid

    type(SPAM3), allocatable :: denskern_at_point(:) ! ddor
    type(SPAM3), allocatable :: pur_denskern_at_point(:) ! ddor

    real(kind=DP), allocatable, dimension(:,:) :: ngwfs_history
    real(kind=DP), allocatable, dimension(:,:) :: cov_grad_history
    real(kind=DP), allocatable, dimension(:,:) :: contra_grad_history

    type(SPAM3), allocatable, dimension(:,:)     :: diis_overlap
    type(SPAM3), allocatable, dimension(:,:)     :: diis_grad_overlap
    real(kind=DP), allocatable, dimension(:,:)   :: diis_grad
    real(kind=DP), allocatable, dimension(:,:)   :: diis_ngwfs

    real(kind=DP), allocatable, dimension(:)   :: diis_cov_grad_on_grid
    real(kind=DP), allocatable, dimension(:)   :: diis_contra_grad_on_grid
    real(kind=DP), allocatable, dimension(:)   :: diis_ngwfs_on_grid
    real(kind=DP), allocatable, dimension(:)   :: start_ngwfs_on_grid

    real(kind=DP), allocatable, dimension(:,:) :: diis_matrix_a
    real(kind=DP), allocatable, dimension(:,:) :: diis_matrix_b
    real(kind=DP), allocatable, dimension(:)   :: diis_coef
    real(kind=DP), allocatable, dimension(:)   :: diis_eig
    real(kind=DP), allocatable, dimension(:)   :: work

    real(kind=DP) :: last_n_energies(3) ! space to store the 3 most recent energies
    real(kind=DP) :: mu(max_spins)
    real(kind=DP) :: rms_gradient
    real(kind=DP) :: previous_rms_gradient ! RMS NGWF grad of previous iteration
    real(kind=DP) :: line_search_coeff
    real(kind=DP) :: lnv_threshold
    real(kind=DP) :: ngwf_threshold
    real(kind=DP) :: F0,F1,F2,Fdiis,G_init,trial_length, rms0, rmsdiis
    real(kind=DP) :: quadratic_coeff,cubic_coeff,rejected_quadratic_coeff
    real(kind=DP) :: predicted_functional
    real(kind=DP) :: previous_g_dot_g
    real(kind=DP) :: current_g_dot_g
    real(kind=DP) :: cg_coeff

    character(len=80), allocatable, dimension(:) :: summary_lines
    integer, allocatable, dimension(:)           :: ipiv
    integer, allocatable, dimension(:)           :: history_idx
    integer :: min_eigenvalue
    integer :: current_idx
    integer :: subspace_size
    integer :: isub, jsub        ! index for the diis-subspace dimensions
    integer :: iteration         ! current iteration
    integer :: est_num_psincs    ! estimate number of psincs in all NGWF spheres
    integer :: current_maxit_lnv ! number of lnv iterations to do for current NGWF step
    integer :: cg_count          ! current number of steps since CG reset
    integer :: is         ! pdh: spin loop counter
    integer :: ierr       ! error flag
    integer :: minit      ! qoh: Minimum number of DIIS iterations
    logical :: trial2     ! ndmh: flag to perform second trial step
    logical :: retrial1   ! ndmh: flag to perform repeat first trial step
    logical :: reversing  ! ndmh: line search is going uphill
    logical :: line_search_success ! ndmh: line search fit success flag
    logical :: use_diis
    logical :: diis_coef_success
    logical :: diis_test_ngwfs
    integer :: ngwf_diis_type
    real(kind=DP) :: ngwf_diis_damping
    logical :: ngwf_diis_upd_sub
    logical :: ngwf_diis_xtpol_gradient
    logical :: ngwf_diis_check_grad
    logical :: ngwf_diis_check_ener
    logical :: ngwf_diis_kernel_upd
    logical :: ngwf_diis_use_cg_dir
    logical :: ngwf_diis_store_cg_dir


    ! cks: FOR PLOTTING FUNCTIONS
    character(len=64) :: fun_name, iter_number
    character(len=200) :: diis_devel_code
    integer            :: start_pos, stop_pos, test_pos

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering ngwf_diis_optimise'
#endif

    ! Flush output
    call services_flush

    ! Start timer
    call timer_clock('ngwf_diis_optimise',1)

    ngwf_diis_type = 3
    ngwf_diis_damping = 0.0_DP
    ngwf_diis_upd_sub = .false.
    ngwf_diis_xtpol_gradient = .false.
    ngwf_diis_check_grad = .true.
    ngwf_diis_check_ener = .false.
    ngwf_diis_kernel_upd = .false.
    ngwf_diis_use_cg_dir = .false.
    ngwf_diis_store_cg_dir = .false.

    ! Check local copy of devel_code string and act accordingly
    if (pub_on_root) then
       diis_devel_code=pub_devel_code
       if (len_trim(diis_devel_code)>0) then
          start_pos=index(diis_devel_code,'DIIS:')
          stop_pos=index(diis_devel_code,':DIIS')
          if (stop_pos<=0) stop_pos=len_trim(diis_devel_code) !missing end so go to end of string
          if (start_pos>0) then

             !look for a setting of diis type
             test_pos=index(diis_devel_code,'TYPE=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len('TYPE=')
                read(diis_devel_code(test_pos:test_pos+ &
                     & index(diis_devel_code(test_pos:stop_pos),':')-2),*) ngwf_diis_type
             end if
             !look for a setting of damping parameter
             test_pos=index(diis_devel_code,'DAMP=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len('DAMP=')
                read(diis_devel_code(test_pos:test_pos+ &
                     & index(diis_devel_code(test_pos:stop_pos),':')-2),*) ngwf_diis_damping
             end if
             !look for diis_update_subspace flag
             if (index(diis_devel_code(start_pos:stop_pos),'UPD_SUB=T')>0) then
                ngwf_diis_upd_sub=.true.
             else if (index(diis_devel_code(start_pos:stop_pos),'UPD_SUB=F')>0) then
                ngwf_diis_upd_sub=.false.
             end if
             !look for diis_extrapolate_gradient flag
             if (index(diis_devel_code(start_pos:stop_pos),'XTPOL_GRAD=T')>0) then
                ngwf_diis_xtpol_gradient=.true.
             else if (index(diis_devel_code(start_pos:stop_pos),'XTPOL_GRAD=F')>0) then
                ngwf_diis_xtpol_gradient=.false.
             end if
             !look for diis_gradient_check flag
             if (index(diis_devel_code(start_pos:stop_pos),'CHECK_GRAD=T')>0) then
                ngwf_diis_check_grad=.true.
             else if (index(diis_devel_code(start_pos:stop_pos),'CHECK_GRAD=F')>0) then
                ngwf_diis_check_grad=.false.
             end if
             !look for diis_energy_check flag
             if (index(diis_devel_code(start_pos:stop_pos),'CHECK_ENER=T')>0) then
                ngwf_diis_check_ener=.true.
             else if (index(diis_devel_code(start_pos:stop_pos),'CHECK_ENER=F')>0) then
                ngwf_diis_check_ener=.false.
             end if
             !look for diis_kernel_update flag
             if (index(diis_devel_code(start_pos:stop_pos),'UPD_KERN=T')>0) then
                ngwf_diis_kernel_upd=.true.
             else if (index(diis_devel_code(start_pos:stop_pos),'UPD_KERN=F')>0) then
                ngwf_diis_kernel_upd=.false.
             end if
             !look for diis_use_cg_directions flag
             if (index(diis_devel_code(start_pos:stop_pos),'USE_CG_DIR=T')>0) then
                ngwf_diis_use_cg_dir=.true.
             else if (index(diis_devel_code(start_pos:stop_pos),'USE_CG_DIR=F')>0) then
                ngwf_diis_use_cg_dir=.false.
             end if
             !look for diis_store_cg_directions flag
             if (index(diis_devel_code(start_pos:stop_pos),'STORE_CG_DIR=T')>0) then
                ngwf_diis_store_cg_dir=.true.
             else if (index(diis_devel_code(start_pos:stop_pos),'STORE_CG_DIR=F')>0) then
                ngwf_diis_store_cg_dir=.false.
             end if

             if (pub_on_root.and.(pub_output_detail>=NORMAL)) then
                write(stdout,'(/a)') 'Diis: Processing of the devel_code string...'
                write(stdout,'(a,i2)') 'Diis: type = ', ngwf_diis_type
                write(stdout,'(a,f6.4)') 'Diis: damping parameter = ', ngwf_diis_damping
                write(stdout,'(a,l1)') 'Diis: xtpol_gradient = ', ngwf_diis_xtpol_gradient
                write(stdout,'(a,l1)') 'Diis: upd_sub = ', ngwf_diis_upd_sub
                write(stdout,'(a,l1)') 'Diis: check_grad = ', ngwf_diis_check_grad
                write(stdout,'(a,l1)') 'Diis: check_ener = ', ngwf_diis_check_ener
                write(stdout,'(a,l1)') 'Diis: kernel_upd = ', ngwf_diis_kernel_upd
                write(stdout,'(a,l1)') 'Diis: use_cg_dir = ', ngwf_diis_use_cg_dir
                write(stdout,'(a,l1)') 'Diis: store_cg_dir = ', ngwf_diis_store_cg_dir
             end if
         end if
       end if
    end if
    call comms_bcast(pub_root_node_id,ngwf_diis_type)
    call comms_bcast(pub_root_node_id,ngwf_diis_damping)
    call comms_bcast(pub_root_node_id,ngwf_diis_xtpol_gradient)
    call comms_bcast(pub_root_node_id,ngwf_diis_upd_sub)
    call comms_bcast(pub_root_node_id,ngwf_diis_check_grad)
    call comms_bcast(pub_root_node_id,ngwf_diis_check_ener)
    call comms_bcast(pub_root_node_id,ngwf_diis_kernel_upd)
    call comms_bcast(pub_root_node_id,ngwf_diis_use_cg_dir)
    call comms_bcast(pub_root_node_id,ngwf_diis_store_cg_dir)

    if (ngwf_diis_reset.lt.ngwf_diis_size) ngwf_diis_reset = ngwf_diis_size



    ! cks: write initial NGWFs in plotting formats
    if (pub_write_ngwf_plot .and. &
         (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
       call visual_ngwfs(rep%ngwfs_on_grid, ngwf_basis, 'initial', elements)
    endif

    ! pdh: allocate local sparse matrices
    allocate(pur_denskern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','pur_denskern',ierr)
    do is=1,pub_cell%num_spins
       call sparse_create(pur_denskern(is),denskern(is))
    end do
    allocate(denskern_at_point(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','denskern_at_point',ierr)
    do is=1, pub_cell%num_spins
       call sparse_create(denskern_at_point(is),denskern(is))
    end do
    allocate(pur_denskern_at_point(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','pur_denskern_at_point',ierr)
    do is=1, pub_cell%num_spins
       call sparse_create(pur_denskern_at_point(is),pur_denskern(is))
    end do

    ! Allocate workspace
    allocate(diis_ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise', &
         'diis_ngwfs_on_grid',ierr)
    allocate(start_ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise', &
         'start_ngwfs_on_grid',ierr)
    allocate(contra_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise', &
         'contra_grad_on_grid',ierr)
    allocate(cov_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise', &
         'cov_grad_on_grid',ierr)
    allocate(direction_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise', &
         'direction_on_grid',ierr)
    allocate(contra_direction_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise', &
         'contra_direction_on_grid',ierr)
    allocate(trial_ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise', &
         'trial_ngwfs_on_grid',ierr)
    allocate(prev_direction_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), &
            stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','prev_direction_on_grid',ierr)
    allocate(prev_contra_direction_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), &
            stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','prev_contra_direction_on_grid',ierr)

    if (ngwf_cg_type == 'NGWF_POLAK') then
       allocate(prev_contra_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), &
            stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','prev_contra_grad_on_grid', &
            ierr)
    end if

    allocate(summary_lines(max(maxit_ngwf_diis+1,2)),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','summary_lines',ierr)

    allocate(cov_grad_history(ngwf_basis%n_ppds*pub_cell%n_pts,ngwf_diis_size), &
         stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','cov_grad_history',ierr)
    allocate(contra_grad_history(ngwf_basis%n_ppds*pub_cell%n_pts,ngwf_diis_size), &
         stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','contra_grad_history',ierr)
    allocate(ngwfs_history(ngwf_basis%n_ppds*pub_cell%n_pts,ngwf_diis_size), &
         stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','ngwfs_history',ierr)
    allocate(diis_cov_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','diis_cov_grad_on_grid',ierr)
    allocate(diis_contra_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','diis_contra_grad_on_grid',ierr)
    allocate(history_idx(ngwf_diis_size),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','history_idx',ierr)
    allocate(diis_coef(ngwf_diis_size),stat=ierr)
    call utils_alloc_check('ngwf_diis_optimise','diis_coef',ierr)

    if (ngwf_diis_type == 1) then
       allocate(diis_matrix_a(ngwf_diis_size+1,ngwf_diis_size+1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_matrix_a',ierr)
       allocate(diis_matrix_b(ngwf_diis_size+1,1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_matrix_b',ierr)
       allocate(ipiv(ngwf_diis_size+1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','ipiv',ierr)
       allocate(diis_grad(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_grad',ierr)
       allocate(diis_eig(1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_eig',ierr)
       allocate(work(1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','work',ierr)
       allocate(diis_ngwfs(1,1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_ngwfs',ierr)
       allocate(diis_overlap(1,1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_overlap',ierr)
       allocate(diis_grad_overlap(1,1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_grad_overlap',ierr)

    elseif (ngwf_diis_type == 2) then
       allocate(diis_matrix_a(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_matrix_a',ierr)
       allocate(diis_matrix_b(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_matrix_b',ierr)
       allocate(diis_eig(ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_eig',ierr)
       allocate(work(8*ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','work',ierr)
       allocate(ipiv(1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','ipiv',ierr)
       allocate(diis_grad(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_grad',ierr)
       allocate(diis_ngwfs(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_ngwfs',ierr)
       allocate(diis_overlap(1,1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_overlap',ierr)
       allocate(diis_grad_overlap(1,1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_grad_overlap',ierr)

    elseif (ngwf_diis_type == 3) then
       allocate(diis_matrix_a(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_matrix_a',ierr)
       allocate(diis_matrix_b(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_matrix_b',ierr)
       allocate(diis_eig(ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_eig',ierr)
       allocate(work(8*ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','work',ierr)
       allocate(ipiv(1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','ipiv',ierr)
       allocate(diis_overlap(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_overlap',ierr)
       allocate(diis_grad_overlap(ngwf_diis_size,ngwf_diis_size),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_grad_overlap',ierr)
       allocate(diis_grad(1,1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_grad',ierr)
       allocate(diis_ngwfs(1,1),stat=ierr)
       call utils_alloc_check('ngwf_diis_optimise','diis_ngwfs',ierr)

    endif

    ! Parameter initialisations
    ngwf_threshold     = ngwf_threshold_orig
    lnv_threshold      = lnv_threshold_orig
    est_num_psincs     = function_basis_est_num_psincs(ngwf_basis)

    ! Variable intialisations
    current_maxit_lnv     = minit_lnv
    previous_rms_gradient = huge(1.0_DP)
    line_search_coeff  = 0.15_DP
    trial_length       = 0.1_DP
    F0=0.0_DP ; F1=0.0_DP ; F2=0.0_DP ; Fdiis=0.0_DP
    G_init=0.0_DP ; diis_coef(:)=0.0_DP ; cg_coeff=0.0_DP
    rms_gradient=1.0_DP ; mu=0.0_DP ; total_energy=0.0_DP ; cg_count=0
    predicted_functional =0.0_DP
    line_search_success = .true.
    trial2 = .false.
    retrial1 = .false.
    reversing = .false.
    last_n_energies(:) = huge(1.0_DP)
    quadratic_coeff = 0.0_DP
    rejected_quadratic_coeff = 0.0_DP
    cubic_coeff = 0.0_DP
    converged = .false.
    subspace_size = 0
    history_idx(:) = 0
    current_idx = 0
    previous_rms_gradient = huge(1.0_DP)
    current_g_dot_g  = 0.0_DP
    previous_g_dot_g = 0.0_DP
    prev_direction_on_grid = 0.0_DP
    prev_contra_direction_on_grid = 0.0_DP


    ! cks: prepare to output summary of NGWF optimisation
    if (pub_on_root) write(summary_lines(1),'(a80)') '|ITER|    RMS GRADIENT   &
         &|     TOTAL ENERGY    |   step   |     Epredicted  '
    summary_lines(1) = adjustl(summary_lines(1))

    ! cks: First message for brief output level
    if (pub_on_root .and. pub_output_detail == BRIEF) then
       write(stdout,'(/a)')'########################################&
            &########################################'
       write(stdout,'(a)')'##################### &
            & NGWF self-consistent optimisation  ######################'
       write(stdout,'(a)')'########################################&
            &########################################'
       write(stdout,'(a80)') summary_lines(1)
    end if

    ! qoh: Allow blank calculations for timings purposes if maxit_ngwf_diis < 0
    minit = 1
    if (maxit_ngwf_diis < 0) minit = 0


    !#################################################################
    !####################  NGWF iteration loop  ######################
    !#################################################################

    iteration = 0 ! ndmh: prevent unitialised variable when maxit_ngwf_diis = 0
    ITERATION_LOOP : do iteration=1,max(maxit_ngwf_diis,minit)

       !==============================================================
       !========= Iteration header ===================================

       if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
            maxit_ngwf_diis > 0) then
          write(stdout,'(/a)')'########################################&
               &########################################'
          write(stdout,'(a,i4.3,a)')'########################## &
               &NGWF DIIS iteration ',iteration,' ############################'
          write(stdout,'(a)')'########################################&
               &########################################'
          if (precond_recip) then
             write(stdout,'(a,f7.4,a)') &
                  '****** Reciprocal space K.E. preconditioning with &
                  &k_zero = ', k_zero, ' a0^-1 *******'
          end if
          if (precond_real) then
             write(stdout,'(a,f7.4,a)') &
                  '************ Real space K.E. preconditioning with &
                  &k_zero = ', k_zero, ' a0^-1 *******'
          end if
       end if
       !call services_flush

       !========= End iteration header ===============================
       !==============================================================

       !==============================================================
       !========= Optimise density kernel ============================


       ! cks: When iteration=1 these matrices are already initialised
       if (iteration > 1) call hamiltonian_dens_indep_matrices(rep, &
            ngwf_basis, proj_basis, hub_proj_basis, hub)

       total_energy = electronic_energy(denskern, pur_denskern, ham, &
            lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
            localpseudo_fine, core_density_fine, ewald_energy, elements, &
            lnv_threshold, current_maxit_lnv, .true.)

       last_n_energies(1) = last_n_energies(2)
       last_n_energies(2) = last_n_energies(3)
       last_n_energies(3) = total_energy

       F0 = electronic_lagrangian(total_energy, rep%overlap, denskern, &
            ham%ham, mu, rep%n_occ)

       !========= End optimise density kernel ========================
       !==============================================================

       if (pub_nnho) then
          call internal_hybridize_ngwfs
       endif

       !==============================================================
       !========= NGWF Gradient ======================================

       call hamiltonian_build_matrix(ham, rep)

       ! pdh: exit if no NGWF optimisation required
       if (maxit_ngwf_diis == 0) then
          converged = .false.
          exit
       end if

       call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, &   ! out
            denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &     ! in
            lhxc_fine, ham, hub, mu, elements)     ! in

       rms0 = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
                                    contra_grad_on_grid,1,cov_grad_on_grid,1)
       call comms_reduce('SUM', rms0)
       rms0 = sqrt(abs(rms0) / real(est_num_psincs, kind=DP))

       !========= End NGWF Gradient ==================================
       !==============================================================

       !==============================================================
       !========= Test convergence ===================================

       converged = internal_test_convergence()
       call comms_bcast(pub_root_node_id,converged)
       if (converged) exit

       !========= End test convergence ===============================
       !==============================================================


       ! Write out energy components
       if ((mod(iteration, 5) == 0 .and. pub_output_detail == VERBOSE) .or. &
            iteration == maxit_ngwf_cg) &
            call hamiltonian_energy_components( &
                 pur_denskern, rep, localpseudo_fine, core_density_fine, &
                 ngwf_basis, hub_proj_basis, hub, ewald_energy, ham%hfexchange)

       ! Write NGWFs of current iteration in plotting formats
       if (pub_write_ngwf_plot .and. &
            (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
          write(iter_number,*) iteration
          write(fun_name,'(a64)') 'iteration'// &
               trim(adjustl(iter_number))
          call visual_ngwfs(rep%ngwfs_on_grid, ngwf_basis, fun_name, elements)
       end if


       !==============================================================
       !========= Functional at initial point ========================

       !rms0 = rms_gradient
       !F0 = electronic_lagrangian(total_energy, rep%overlap, denskern, &
       !     mu, rep%n_occ)

       call internal_find_direction

       !========= End functional at initial point ====================
       !==============================================================

       !##############################################################
       !######### Update DIIS history ################################

       if (pub_on_root .and. (pub_output_detail >= NORMAL) ) write(stdout,'(/a)') &
            '--------------------------------- RM-DIIS step --&
            &-------------------------------'

       if (ngwf_diis_store_cg_dir) then
          call internal_update_diis_subspace(rep%ngwfs_on_grid,&
                                direction_on_grid,contra_direction_on_grid)
       else
          call internal_update_diis_subspace(rep%ngwfs_on_grid,&
                                cov_grad_on_grid,contra_grad_on_grid)
       endif

       !######### End update DIIS history ############################
       !##############################################################


       !##############################################################
       !######### Compute DIIS ngwfs #################################

       if (subspace_size .gt.1) then
          start_ngwfs_on_grid = rep%ngwfs_on_grid
          call internal_diis_ngwfs
       endif

       !######### End compute DIIS ngwfs #############################
       !##############################################################


       !##############################################################
       !######### Functional at DIIS #################################

       if (subspace_size .gt.1 .and. diis_test_ngwfs) then

          ! Functional at DIIS NGWFs
          rep%ngwfs_on_grid = diis_ngwfs_on_grid
          call hamiltonian_dens_indep_matrices(rep, ngwf_basis, &
               proj_basis, hub_proj_basis, hub)

          if (ngwf_diis_kernel_upd .or. pub_kernel_update) then
             do is=1,pub_cell%num_spins
                call sparse_copy(denskern_at_point(is),denskern(is))
                call sparse_copy(pur_denskern_at_point(is),pur_denskern(is))
             enddo

             Fdiis = electronic_energy(denskern, pur_denskern, ham, &
                  lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
                  localpseudo_fine, core_density_fine, ewald_energy, elements, &
                  lnv_threshold, current_maxit_lnv, .true.)
          else
             Fdiis = electronic_energy(denskern, pur_denskern, ham, &
                  lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
                  localpseudo_fine, core_density_fine, ewald_energy, elements, &
                  lnv_threshold, current_maxit_lnv, .false.)
          endif

          Fdiis = electronic_lagrangian(Fdiis,rep%overlap,denskern,ham%ham,mu,rep%n_occ)

          if (pub_on_root .and. (pub_output_detail >= NORMAL) ) &
             write(stdout,'(a,f12.8,a,f12.8,a)') 'DIIS predicted functional :  Fdiis   = ', Fdiis,   ', (F0   = ', F0, ')'


       endif

       !######### End functional at DIIS #############################
       !##############################################################



       !##############################################################
       !######### Gradients at DIIS ##################################

       if (subspace_size .gt.1 .and. diis_test_ngwfs) then

          ! Gradients at DIIS NGWFS
          if (ngwf_diis_xtpol_gradient) then
             diis_cov_grad_on_grid(:) = 0.0_dp
             diis_contra_grad_on_grid(:) = 0.0_dp
             do isub=1,subspace_size
                diis_cov_grad_on_grid = diis_cov_grad_on_grid + &
                            diis_coef(isub)*cov_grad_history(:,isub)
                diis_contra_grad_on_grid = diis_contra_grad_on_grid + &
                            diis_coef(isub)*contra_grad_history(:,isub)
             enddo
          else
             call hamiltonian_build_matrix(ham, rep)
             call ngwf_gradient_lnv(diis_contra_grad_on_grid, &
                  diis_cov_grad_on_grid, & ! out
                  denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &  ! in
                  lhxc_fine, ham, hub, mu, elements)        ! in
          endif

          rmsdiis = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
                                       diis_contra_grad_on_grid,1,diis_cov_grad_on_grid,1)
          call comms_reduce('SUM', rmsdiis)
          rmsdiis = sqrt(abs(rmsdiis) / real(est_num_psincs, kind=DP))

          if (pub_on_root .and. (pub_output_detail >= NORMAL) ) &
             write(stdout,'(a,f12.8,a,f12.8,a)') 'DIIS predicted RMS :         RMSdiis = ', rmsdiis, ', (RMS0 = ', rms0, ')'

       endif


       !######### End gradients at DIIS ##############################
       !##############################################################


       !##############################################################
       !######### Test DIIS ##########################################

       use_diis = .false.
       if (subspace_size .gt.1 .and. diis_test_ngwfs) then
          use_diis = .true.
          if (ngwf_diis_check_ener .and. Fdiis .gt. F0) use_diis = .false.
          if (ngwf_diis_check_grad .and. rmsdiis .gt. rms0) use_diis = .false.
       endif

       if (subspace_size .gt.1 .and. use_diis) then
          start_ngwfs_on_grid = rep%ngwfs_on_grid
          cov_grad_on_grid(:) = diis_cov_grad_on_grid(:)
          contra_grad_on_grid(:) = diis_contra_grad_on_grid(:)
          F0 = Fdiis
          rms0 = rmsdiis
          call internal_find_direction
          if (pub_on_root) then
             write(stdout,'(a)') 'DIIS ngwfs accepted ! '
          endif

       elseif (subspace_size .gt.1) then
          rep%ngwfs_on_grid = start_ngwfs_on_grid
          if (ngwf_diis_kernel_upd .or. pub_kernel_update) then
             do is=1,pub_cell%num_spins
                call sparse_copy(denskern(is),denskern_at_point(is))
                call sparse_copy(pur_denskern(is),pur_denskern_at_point(is))
             enddo
          endif
          if (pub_on_root) then
             write(stdout,'(a)') 'DIIS ngwfs rejected ! '
          endif

       endif
       if (pub_on_root .and. (pub_output_detail >= NORMAL) ) write(stdout,'(/a)') &
            '--------------------------------- RM-DIIS step --&
            &-------------------------------'

       if (pub_on_root .and. (pub_output_detail >= NORMAL) ) write(stdout,'(/a)') &
            '------------------------------- End RM-DIIS step --&
            &-----------------------------'



       !######### End test DIIS ######################################
       !##############################################################


       !==============================================================
       !========= Line search by fitting parabola ====================

       if (pub_on_root .and. (pub_output_detail >= NORMAL) ) write(stdout,'(/a)') &
            '------------------------------- NGWF line search &
            &-------------------------------'

       prev_direction_on_grid = direction_on_grid
       prev_contra_direction_on_grid = contra_direction_on_grid
       previous_g_dot_g = current_g_dot_g

       if (ngwf_cg_type == 'NGWF_POLAK') &
            prev_contra_grad_on_grid = contra_grad_on_grid

       ! Functional at trial NGWFs
       start_ngwfs_on_grid = rep%ngwfs_on_grid
       rep%ngwfs_on_grid = start_ngwfs_on_grid + trial_length*direction_on_grid

       call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
            hub_proj_basis, hub)

       F1 = electronic_energy(denskern, pur_denskern, ham, &
            lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
            localpseudo_fine, core_density_fine, ewald_energy, elements, &
            lnv_threshold, current_maxit_lnv, pub_kernel_update)

       F1 = electronic_lagrangian(F1,rep%overlap,denskern,ham%ham,mu,rep%n_occ)

       ! Search by fitting parabola
       call services_line_search_parabola(&
            quadratic_coeff, predicted_functional,line_search_success, & !output
            G_init, F0, F1, trial_length, ngwf_diis_max_step)              !input

       line_search_coeff = quadratic_coeff

       !========= End line search by fitting parabola ================
       !==============================================================


       !==============================================================
       !========= Advanced line search initialization ================

       if (G_init * quadratic_coeff > 0.0_DP) then
          trial2 = .true.
          trial_ngwfs_on_grid = start_ngwfs_on_grid + &
               2.0_DP * trial_length * direction_on_grid
       end if

       ! ndmh: protection against bad line search results: redo trial step if
       ! ndmh: line search result is too much bigger or smaller than trial step
       ! ndmh: or if we are searching 'uphill' and went past trial step position
       if ( ((line_search_coeff < 0.05_DP*trial_length) .or. &
            (line_search_coeff > 20.0_DP*trial_length) .or. &
            (.not. line_search_success) .or. &
            (reversing .and. (line_search_coeff > trial_length))) &
            .and. .not. trial2 ) then
          retrial1 = .true.
       end if

       if (retrial1) then
          rejected_quadratic_coeff = quadratic_coeff
          trial_ngwfs_on_grid = start_ngwfs_on_grid + &
               quadratic_coeff * direction_on_grid
       end if

       !========= End advanced line search initialization ============
       !==============================================================


       !==============================================================
       !========= Advanced line search ===============================


       if (trial2 .or. retrial1) then

          ! Functional at trial length 2
          call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
               hub_proj_basis, hub)

          F2 = electronic_energy(denskern, pur_denskern, ham, &
               lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
               localpseudo_fine, core_density_fine, ewald_energy, elements, &
               lnv_threshold, current_maxit_lnv, pub_kernel_update)

          F2 = electronic_lagrangian(F2,rep%overlap,denskern,ham%ham,mu,rep%n_occ)

          ! Cubic fit
          if (trial2) then
             call services_cubic_fit_minimum( &
                  cubic_coeff, predicted_functional, line_search_success, & !output
                  F0, F1, F2, G_init, trial_length, 2.0_DP*trial_length, &  !input
                  ngwf_diis_max_step)                                         !input
             line_search_coeff = cubic_coeff
          end if

          ! Quadratic fit at new trial length
          if (retrial1) then
             call services_line_search_parabola(&
                  quadratic_coeff, predicted_functional, &               !output
                  line_search_success, G_init, F0, F2, &                 !input
                  line_search_coeff, ngwf_diis_max_step)                   !input
             line_search_coeff = quadratic_coeff
          end if

       else
          ! ndmh: no second trial step
          F2 = 0.0_DP
       end if

       !========= End advanced line search ===========================
       !==============================================================

       ! ndmh: reset flags
       trial2 = .false.
       retrial1 = .false.
       reversing = .false.
       F2 = 0.0_DP


       !==============================================================
       !========= New values of NGWFs ================================

       ! cks: set the new values of the NGWFs
       rep%ngwfs_on_grid = start_ngwfs_on_grid + &
            line_search_coeff * direction_on_grid

       !========= End new values NGWFs ===============================
       !==============================================================

       if (pub_on_root .and. pub_output_detail >= NORMAL) &
            call internal_print_search_info


       ! cks: write new NGWFs to file if required
       if (write_tightbox_ngwfs.and.(.not.pub_write_converged_dk_ngwfs)) then
          ! cks: in universal tightbox representation
          call restart_ngwfs_tightbox_output(rep%ngwfs_on_grid,ngwf_basis,&
               elements, 'tightbox_ngwfs')
       endif

       ! cks: write line with summary info for curent iteration
       if (pub_on_root) then
          if(abs(total_energy)<100000_DP) then
             write(summary_lines(iteration+1),'(i4,f21.14,f22.14,f11.6,f22.14)')&
                  iteration,rms_gradient,total_energy,line_search_coeff, &
                  predicted_functional
          else
             write(summary_lines(iteration+1),'(i4,f21.14,f22.12,f11.6,f22.12)')&
                  iteration,rms_gradient,total_energy,line_search_coeff, &
                  predicted_functional
          end if
          if (pub_output_detail == BRIEF) then
             write(stdout,'(a80)') summary_lines(iteration+1)
          end if
       end if

    enddo ITERATION_LOOP

    !#################################################################
    !#################  End of NGWF iteration loop  ##################
    !#################################################################


    ! write NGWFs if required and minimisation ends in first iteration
    ! ndmh: or if writing final converged results only
    if (write_tightbox_ngwfs .and. ((iteration.eq.1).or. &
        pub_write_converged_dk_ngwfs)) then
       call restart_ngwfs_tightbox_output(rep%ngwfs_on_grid,ngwf_basis,&
            elements, 'tightbox_ngwfs')
    endif

    ! smmd: store NGWFs for extrapolation at subsequent MD/geom steps
    if (store_tightbox_ngwfs) &
       call restart_ngwfs_tightbox_store(rep%ngwfs_on_grid,ngwf_basis,&
            elements)

    ! ndmh: write denskern here if we have not been writing it as we go along
    if (write_denskern .and. pub_write_converged_dk_ngwfs) &
         call restart_kernel_write(denskern)

    ! smmd: store denskern for extrapolation at subsequent MD/geom steps
    if (store_denskern) &
       call restart_kernel_store(denskern)

    ! ars: write in spherical waves representation
    if (pub_write_sw_ngwfs) then! .and. iteration.eq.1) then
       call restart_sph_waves_output(rep%ngwfs_on_grid,ngwf_basis,elements,"sw_ngwfs")
    endif

    if (.not. converged) then
       ! cks: reset iteration number just for storing final line of calculation
       !      summary
       iteration = iteration - 1
       ! cks: print warning that calculation failed to converge
       if (pub_on_root) write(stdout,'(a,i4,a)') &
            'WARNING: maximum number of NGWF DIIS iterations (',maxit_ngwf_diis, &
            ') exceeded!'
    end if

    ! cks: print calculation summary
    call internal_calculation_summary

    ! cks: print quality control information
    if (print_qc) call internal_qc_output

    ! ddor: Write out DFT+U occupancies if it hasn't already been done
    if (pub_hubbard.and.pub_on_root.and.(pub_output_detail .ne. VERBOSE)) then
       call hubbard_energy_info(hub,hub_proj_basis)
    endif

    ! clean up various routines
    call ngwf_gradient_exit


    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(pur_denskern(is))
    end do
    deallocate(pur_denskern,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'pur_denskern',ierr)

    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(denskern_at_point(is))
    end do
    deallocate(denskern_at_point,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'denskern_at_point',ierr)

    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(pur_denskern_at_point(is))
    end do
    deallocate(pur_denskern_at_point,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'pur_denskern_at_point',ierr)

    deallocate(diis_ngwfs_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_ngwfs_on_grid', ierr)
    deallocate(start_ngwfs_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'start_ngwfs_on_grid',ierr)
    deallocate(contra_grad_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'contra_grad_on_grid',ierr)
    deallocate(cov_grad_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'cov_grad_on_grid',ierr)
    deallocate(direction_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'direction_on_grid',ierr)
    deallocate(contra_direction_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'contra_direction_on_grid',ierr)
    deallocate(trial_ngwfs_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'trial_ngwfs_on_grid',ierr)
    deallocate(prev_direction_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'prev_direction_on_grid',ierr)
    deallocate(prev_contra_direction_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'prev_contra_direction_on_grid',ierr)

    if (ngwf_cg_type == 'NGWF_POLAK') then
       deallocate(prev_contra_grad_on_grid,stat=ierr)
       call utils_dealloc_check('ngwf_diis_optimise', &
            'prev_contra_grad_on_grid', ierr)
    end if

    deallocate(summary_lines,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'summary_lines',ierr)
    deallocate(cov_grad_history,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'cov_grad_history',ierr)
    deallocate(contra_grad_history,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'contra_grad_history',ierr)
    deallocate(ngwfs_history,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'ngwfs_history',ierr)
    deallocate(diis_cov_grad_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_cov_grad_on_grid',ierr)
    deallocate(diis_contra_grad_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_contra_grad_on_grid',ierr)
    deallocate(history_idx,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'history_idx',ierr)
    deallocate(diis_coef,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_coef',ierr)
    deallocate(diis_matrix_a,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_matrix_a',ierr)
    deallocate(diis_matrix_b,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_matrix_b',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'ipiv',ierr)
    deallocate(diis_grad,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_grad',ierr)
    deallocate(diis_eig,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_eig',ierr)
    deallocate(work,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'work',ierr)
    deallocate(diis_ngwfs,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_ngwfs',ierr)
    deallocate(diis_overlap,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_overlap',ierr)
    deallocate(diis_grad_overlap,stat=ierr)
    call utils_dealloc_check('ngwf_diis_optimise', &
         'diis_grad_overlap',ierr)

    call timer_clock('ngwf_diis_optimise',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving ngwf_diis_optimise'
#endif

    ! Flush output
    call services_flush

    return

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_update_diis_subspace(new_ngwfs,new_cov_res,new_contra_res)

      !==============================================================!
      ! Update subspace size, history_index as well as history datas !
      !--------------------------------------------------------------!
      ! Written by Simon M.-M. Dubois on September 2010              !
      !==============================================================!

      implicit none

      ! Arguments
      real(kind=DP), intent(in) :: new_ngwfs(ngwf_basis%n_ppds * pub_cell%n_pts)
      real(kind=DP), intent(in) :: new_cov_res(ngwf_basis%n_ppds * pub_cell%n_pts)
      real(kind=DP), intent(in) :: new_contra_res(ngwf_basis%n_ppds * pub_cell%n_pts)

      !==============================================================
      !========= Update subspace size and history index =============

      ! If required reset the diis subspace
      if (mod(iteration,ngwf_diis_reset+1)==0) then
         subspace_size = 0
         history_idx(:) = 0
      endif

      ! Update history index
      do isub = 1,min(subspace_size,ngwf_diis_size)
         history_idx(isub) = history_idx(isub) + 1
      enddo

      ! Increase the size of the diis subspace if required
      if (subspace_size.lt.ngwf_diis_size) then
         subspace_size = subspace_size + 1
         current_idx  = subspace_size
      else
         current_idx = maxloc(history_idx(1:ngwf_diis_size),1)
      endif
      history_idx(current_idx) = 1

      if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
         write(stdout,'(a,i3.3)') 'Diis subspace size ', subspace_size
         write(stdout,'(a,i3.3)') 'Diis current index ', current_idx
      endif

      !==============================================================
      !========= Update DIIS history ================================

      cov_grad_history(:,current_idx) = new_cov_res(:)
      contra_grad_history(:,current_idx) = new_contra_res(:)
      ngwfs_history(:,current_idx) = new_ngwfs(:)

      if (ngwf_diis_type == 1) then

         !===== Update diis_grad ====================================
         diis_grad(current_idx,current_idx) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
             new_contra_res(:),1,new_cov_res(:),1)
         call comms_reduce('SUM',diis_grad(current_idx,current_idx))
         do isub = 1, min(subspace_size,ngwf_diis_size)
            if (isub == current_idx) cycle
            diis_grad(current_idx,isub) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
                new_contra_res,1,cov_grad_history(:,isub),1)
            call comms_reduce('SUM', diis_grad(current_idx,isub))
            diis_grad(isub,current_idx) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
                contra_grad_history(:,isub),1,new_cov_res,1)
            call comms_reduce('SUM', diis_grad(isub,current_idx))
         enddo


      elseif (ngwf_diis_type == 2) then

         !===== Update diis_grad ====================================
         diis_grad(current_idx,current_idx) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
             new_cov_res(:),1,new_cov_res(:),1)
         call comms_reduce('SUM', diis_grad(current_idx,current_idx))
         do isub = 1, min(subspace_size,ngwf_diis_size)
            if (isub .gt. current_idx) then
               diis_grad(current_idx,isub) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
                     new_cov_res,1,cov_grad_history(:,isub),1)
               call comms_reduce('SUM', diis_grad(current_idx,isub))
            elseif (isub .lt. current_idx) then
               diis_grad(isub,current_idx) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
                   cov_grad_history(:,isub),1,new_cov_res,1)
               call comms_reduce('SUM', diis_grad(isub,current_idx))
            endif
         enddo

         !===== Update diis_ngwfs ====================================
         diis_ngwfs(current_idx,current_idx) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
             new_ngwfs,1,new_ngwfs,1)
         call comms_reduce('SUM', diis_ngwfs(current_idx,current_idx))
         do isub = 1, min(subspace_size,ngwf_diis_size)
            if (isub .gt. current_idx) then
               diis_ngwfs(current_idx,isub) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
                   new_ngwfs,1,ngwfs_history(:,isub),1)
               call comms_reduce('SUM', diis_ngwfs(current_idx,isub))
            elseif (isub .lt. current_idx) then
               diis_ngwfs(isub,current_idx) = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
                   ngwfs_history(:,isub),1,new_ngwfs,1)
               call comms_reduce('SUM', diis_ngwfs(isub,current_idx))
            endif
         enddo

      elseif (ngwf_diis_type == 3) then

         !===== Update diis_overlap =================================
         if (iteration .le. subspace_size) then
            call sparse_create(diis_overlap(current_idx,current_idx),rep%overlap)
         endif
         call sparse_copy(diis_overlap(current_idx,current_idx),rep%overlap)
         do isub = 1, min(subspace_size,ngwf_diis_size)
            if (isub == current_idx) cycle
            if (iteration .le. subspace_size) then
               call sparse_create(diis_overlap(current_idx,isub),rep%overlap)
               call sparse_create(diis_overlap(isub,current_idx),rep%overlap)
            endif
            if (isub .gt. current_idx) then
               call integrals_brappd_ketppd(diis_overlap(current_idx,isub),new_ngwfs,ngwf_basis,ngwfs_history(:,isub),ngwf_basis)
            elseif (isub .lt. current_idx) then
               call integrals_brappd_ketppd(diis_overlap(isub,current_idx),ngwfs_history(:,isub),ngwf_basis,new_ngwfs,ngwf_basis)
            endif
         enddo

         !===== Update diis_grad_overlap ============================
         if (iteration .le. subspace_size) then
            call sparse_create(diis_grad_overlap(current_idx,current_idx),rep%overlap)
         endif
         call integrals_brappd_ketppd(diis_grad_overlap(current_idx,current_idx),new_cov_res,ngwf_basis,new_cov_res,ngwf_basis)
         do isub = 1, min(subspace_size,ngwf_diis_size)
            if (isub == current_idx) cycle
            if (iteration .le. subspace_size) then
               call sparse_create(diis_grad_overlap(current_idx,isub),rep%overlap)
               call sparse_create(diis_grad_overlap(isub,current_idx),rep%overlap)
            endif
            if (isub .gt. current_idx) then
               call integrals_brappd_ketppd(diis_grad_overlap(current_idx,isub), &
                    new_cov_res,ngwf_basis,cov_grad_history(:,isub),ngwf_basis)
            elseif (isub .lt. current_idx) then
               call integrals_brappd_ketppd(diis_grad_overlap(isub,current_idx), &
                    cov_grad_history(:,isub),ngwf_basis,new_cov_res,ngwf_basis)
            endif
         enddo

      endif

    end subroutine internal_update_diis_subspace

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_diis_ngwfs

      !==============================================================!
      ! Update subspace size, history_index as well as history datas !
      !--------------------------------------------------------------!
      ! Written by Simon M.-M. Dubois on September 2010              !
      !==============================================================!

      implicit none

      ! Local variables
      real(kind=DP) :: norm_ngwfs
      real(kind=DP) :: diis_norm_ngwfs
      character(len=20) :: fmt,tmp

      norm_ngwfs = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts,rep%ngwfs_on_grid,1,rep%ngwfs_on_grid,1)
      call comms_reduce('SUM', norm_ngwfs)

      !==============================================================
      !===== Linear DIIS ============================================
      !==============================================================
      if (ngwf_diis_type == 1) then

         ! Initialize the DIIS matrix A
         do isub=1,subspace_size
            diis_matrix_a(isub,1:subspace_size) = diis_grad(isub,1:subspace_size)/diis_grad(isub,isub)
            diis_matrix_a(isub,isub) = diis_matrix_a(isub,isub)*(1+ngwf_diis_damping)
            diis_matrix_a(isub,subspace_size+1) = -1.0_dp / diis_grad(isub,isub)
         enddo
         diis_matrix_a(subspace_size+1,:) = 1.0_dp
         diis_matrix_a(subspace_size+1,subspace_size+1) = 0.0_dp

         ! Initialize the DIIS matrix B
         diis_matrix_b(1:subspace_size,1) = 0.0_dp
         diis_matrix_b(subspace_size+1,1) = 1.0_dp

         ! Compute the DIIS coefficients
         call comms_barrier
         call dgesv(subspace_size+1,1,diis_matrix_a(1:subspace_size+1,1:subspace_size+1), &
                    subspace_size+1,ipiv,diis_matrix_b(1:subspace_size+1,1), subspace_size+1, &
                    ierr)
         if (ierr.ne.0) then
            write(stderr,'(a,i6)') 'Error in ngwf_diis_optimise: &
                 &dgesv failed with code ',ierr
            diis_coef_success = .false.
            diis_test_ngwfs = .false.
         else
            diis_coef_success = .true.
            diis_test_ngwfs = .true.
         endif

         if (diis_coef_success) then
            diis_coef(1:subspace_size) = diis_matrix_b(1:subspace_size,1)

            if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
               write(stdout,'(a)') 'DIIS coefficents : '
               write(tmp,'(i5)') subspace_size
               write(fmt,'(a,a,a)') '(',trim(adjustl(tmp)),'(f8.4,x),f8.4)'
               write(stdout,fmt) diis_coef(1:subspace_size)
            endif

            ! Perform the DIIS step
            diis_ngwfs_on_grid(:) = 0.0_dp
            do isub=1,subspace_size
               diis_ngwfs_on_grid = diis_ngwfs_on_grid + diis_coef(isub)*ngwfs_history(:,isub)
            enddo
         endif

      !==============================================================
      !===== Eigenvalue DIIS ========================================
      !==============================================================
      elseif (ngwf_diis_type == 2) then

         ! Initialize the DIIS matrices
         diis_matrix_b(1:subspace_size,1:subspace_size) = diis_grad(1:subspace_size,1:subspace_size)
         diis_matrix_a(1:subspace_size,1:subspace_size) = diis_ngwfs(1:subspace_size,1:subspace_size)

         call comms_barrier
         call dsygv(1,'V','U',subspace_size,diis_matrix_b(1:subspace_size,1:subspace_size), &
                    subspace_size,diis_matrix_a(1:subspace_size,1:subspace_size), subspace_size, &
                    diis_eig(1:subspace_size),work,8*subspace_size,ierr)

         if (ierr.ne.0) then
            write(stderr,'(a,i6)') 'Error in ngwf_diis_optimise: &
                 &dsygv failed with code ',ierr
            diis_coef_success = .false.
            diis_test_ngwfs = .false.
         else
            diis_coef_success = .true.
            diis_test_ngwfs = .true.
         endif

         if (diis_coef_success) then
            min_eigenvalue = minloc(abs(diis_eig(1:subspace_size)),1)
            diis_coef(1:subspace_size) = &
                 sign(1.0_dp,diis_matrix_b(current_idx,min_eigenvalue)) * &
                 diis_matrix_b(1:subspace_size,min_eigenvalue)

            ! Perform the DIIS step
            diis_ngwfs_on_grid(:) = 0.0_dp
            do isub=1,subspace_size
               diis_ngwfs_on_grid = diis_ngwfs_on_grid + diis_coef(isub)*ngwfs_history(:,isub)
            enddo

            ! Compute the norm
            diis_norm_ngwfs = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts,diis_ngwfs_on_grid,1,diis_ngwfs_on_grid,1)
            call comms_reduce('SUM', diis_norm_ngwfs)

            ! Normalise the ngwfs
            diis_ngwfs_on_grid = diis_ngwfs_on_grid * sqrt(norm_ngwfs/diis_norm_ngwfs)
            diis_coef(1:subspace_size) = diis_coef(1:subspace_size) * sqrt(norm_ngwfs/diis_norm_ngwfs)

            if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
               write(stdout,'(a)') 'DIIS coefficents : '
               write(tmp,'(i5)') subspace_size
               write(fmt,'(a,a,a)') '(',trim(adjustl(tmp)),'(f8.4,x),f8.4)'
               write(stdout,fmt)  diis_coef(1:subspace_size)
            endif

         endif

      !==============================================================
      !===== Eigenvalue DIIS ========================================
      !==============================================================
      elseif (ngwf_diis_type == 3) then

         ! Initialize the DIIS matrices
         do isub = 1, min(subspace_size,ngwf_diis_size)
            do jsub = 1, min(subspace_size,ngwf_diis_size)
               diis_matrix_b(isub,jsub) = sparse_trace(rep%inv_overlap,diis_grad_overlap(isub,jsub))
               diis_matrix_a(isub,jsub) = sparse_trace(rep%inv_overlap,diis_overlap(isub,jsub))
            enddo
         enddo

         call dsygv(1,'V','U',subspace_size,diis_matrix_b(1:subspace_size,1:subspace_size), &
                    subspace_size,diis_matrix_a(1:subspace_size,1:subspace_size), subspace_size, &
                    diis_eig(1:subspace_size),work,8*subspace_size,ierr)

         if (ierr.ne.0) then
            write(stderr,'(a,i6)') 'Error in ngwf_diis_optimise: &
                 &dsygv failed with code ',ierr
            diis_coef_success = .false.
            diis_test_ngwfs = .false.
         else
            diis_coef_success = .true.
            diis_test_ngwfs = .true.
         endif

         if (diis_coef_success) then
            min_eigenvalue = minloc(abs(diis_eig(1:subspace_size)),1)
            diis_coef(1:subspace_size) = &
                 sign(1.0_dp,diis_matrix_b(current_idx,min_eigenvalue)) * &
                 diis_matrix_b(1:subspace_size,min_eigenvalue)

            ! Perform the DIIS step
            diis_ngwfs_on_grid(:) = 0.0_dp
            do isub=1,subspace_size
               diis_ngwfs_on_grid = diis_ngwfs_on_grid + diis_coef(isub)*ngwfs_history(:,isub)
            enddo

            ! Compute the norm
            diis_norm_ngwfs = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts,diis_ngwfs_on_grid,1,diis_ngwfs_on_grid,1)
            call comms_reduce('SUM', diis_norm_ngwfs)

            ! Normalise the ngwfs
            diis_ngwfs_on_grid = diis_ngwfs_on_grid * sqrt(norm_ngwfs/diis_norm_ngwfs)
            diis_coef(1:subspace_size) = diis_coef(1:subspace_size) * sqrt(norm_ngwfs/diis_norm_ngwfs)

            if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
               write(stdout,'(a)') 'DIIS coefficents : '
               write(tmp,'(i5)') subspace_size
               write(fmt,'(a,a,a)') '(',trim(adjustl(tmp)),'(f8.4,x),f8.4)'
               write(stdout,fmt)  diis_coef(1:subspace_size)
            endif


         endif

      endif

    end subroutine internal_diis_ngwfs

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_transform_kernel(new_ngwfs_on_grid,old_ngwfs_on_grid)

      !==============================================================!
      ! Transforms the density kernel to its representation in terms !
      ! of new NGWFs                                                 !
      !--------------------------------------------------------------!
      ! Written by Nicholas Hine on 22/10/09                         !
      !==============================================================!

      use integrals, only: integrals_brappd_ketppd
      use sparse, only: SPAM3, sparse_create,  sparse_transpose, &
           sparse_product, sparse_destroy

      implicit none

      ! Arguments
      real(kind=DP), intent(in) :: new_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      real(kind=DP), intent(in) :: old_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)

      ! Local Variables
      type(SPAM3) :: overlap_old_new
      type(SPAM3) :: sk,ks,ksk

      ! Create temporary matrix structures
      call sparse_create(overlap_old_new,rep%overlap)
      call sparse_create(sk,rep%overlap,denskern(1))
      call sparse_create(ks,denskern(1),rep%overlap)
      call sparse_create(ksk,ks,denskern(1))

      ! Calculate overlap of old NGWFs with new NGWFs
      call integrals_brappd_ketppd(overlap_old_new, &
           old_ngwfs_on_grid,ngwf_basis,new_ngwfs_on_grid,ngwf_basis)

      ! Calculate sk_a^b = <f_a|f'_c>.(S^-1)^cb
      call sparse_product(sk,overlap_old_new,rep%inv_overlap)
      ! Transpose to find ks^a_b = (S^-1)^ac.<f'_c|f_b>
      call sparse_transpose(ks,sk)

      ! Calculate density kernel transformed to new NGWF basis
      ! K'^ah = (S^-1)^ab.<f'_b|f_g> K^ge <f_e|f'_z>.(S^-1)^zh
      do is=1,pub_cell%num_spins
         call sparse_product(ksk,ks,denskern(is))
         call sparse_product(denskern(is),ksk,sk)
      end do

      ! Destroy temporary matrix structures
      call sparse_destroy(ksk)
      call sparse_destroy(ks)
      call sparse_destroy(sk)
      call sparse_destroy(overlap_old_new)

    end subroutine internal_transform_kernel

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_hybridize_ngwfs

      !================================================!
      ! Construct atomic hybrids from current NGWFs.   !
      !------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 17/02/2006 !
      !================================================!

      use nho, only: nho_generate
      use projectors, only: projectors_func_ovlp_box
      use pseudopotentials, only: nlps_projectors
      use rundat, only: pub_any_nl_proj
      use sparse, only: SPAM3, sparse_create, sparse_scale, sparse_axpy,&
           sparse_transpose, sparse_product, sparse_destroy
      implicit none

      type(SPAM3) :: ao_to_hy
      type(SPAM3) :: trans_ao_to_hy
      type(SPAM3) :: hy_to_ao
      type(SPAM3) :: trans_hy_to_ao
      type(SPAM3) :: k_hao
      type(SPAM3) :: h_aoh
      type(SPAM3) :: kernel

      ao_to_hy%structure ='D'
      call sparse_create(ao_to_hy)
      trans_ao_to_hy%structure ='D'
      call sparse_create(trans_ao_to_hy)

      hy_to_ao%structure ='D'
      call sparse_create(hy_to_ao)
      trans_hy_to_ao%structure ='D'
      call sparse_create(trans_hy_to_ao)

      call sparse_create(k_hao, denskern(1), hy_to_ao)
      call sparse_create(h_aoh, ham%ham(1), ao_to_hy)

      ! pdh: take average of up and down kernels for spin polarisation
      call sparse_create(kernel, denskern(1))
      do is=1,pub_cell%num_spins
         call sparse_axpy(kernel,denskern(is),1.0_DP)
      end do
      call sparse_scale(kernel,1.0_DP/pub_cell%num_spins)

      ! cks: convert the NGWFs to NHOs
      call nho_generate(ao_to_hy, hy_to_ao, &  ! output
           rep%ngwfs_on_grid, prev_direction_on_grid, &   ! in-out
           kernel, rep%overlap, elements, ngwf_basis)     ! input

      call sparse_transpose(trans_ao_to_hy, ao_to_hy)
      call sparse_transpose(trans_hy_to_ao, hy_to_ao)

      ! cks: transform denskern
      do is=1,pub_cell%num_spins
         call sparse_product(k_hao, denskern(is), trans_hy_to_ao)
         call sparse_product(denskern(is), hy_to_ao, k_hao)
      end do
      ! ndmh: mark sk and ks workspaces in kernel mod as invalid
      !call kernel_workspace_invalidate()

      ! cks: transform pur_denskern
      do is=1,pub_cell%num_spins
         call sparse_product(k_hao, pur_denskern(is), trans_hy_to_ao)
         call sparse_product(pur_denskern(is), hy_to_ao, k_hao)
      end do

      ! cks: transform inverse overlap
      call sparse_product(k_hao, rep%inv_overlap, trans_hy_to_ao)
      call sparse_product(rep%inv_overlap, hy_to_ao, k_hao)

      ! cks: transform overlap
      call sparse_product(h_aoh, rep%overlap, ao_to_hy)
      call sparse_product(rep%overlap, trans_ao_to_hy, h_aoh)

      ! cks: transform kinetic
      call sparse_product(h_aoh, rep%kinet, ao_to_hy)
      call sparse_product(rep%kinet, trans_ao_to_hy, h_aoh)

      ! cks: transform nonlocpot
      if (pub_any_nl_proj) then
         call sparse_product(h_aoh, rep%nonlocpot(1), ao_to_hy)
         call sparse_product(rep%nonlocpot(1), trans_ao_to_hy, h_aoh)
      end if
      if (pub_paw) then
         do is=1,pub_cell%num_spins
            call sparse_product(h_aoh, ham%nonlocpot(is), ao_to_hy)
            call sparse_product(ham%nonlocpot(is), trans_ao_to_hy, h_aoh)
         end do
      end if

      ! cks: transform lhxc
      do is=1,pub_cell%num_spins
         call sparse_product(h_aoh, ham%lhxc(is), ao_to_hy)
         call sparse_product(ham%lhxc(is), trans_ao_to_hy, h_aoh)
      end do

      ! cks: temporary measure - recompute sp_overlap matrix
      if (pub_any_nl_proj) then
         call projectors_func_ovlp_box(rep%sp_overlap, &
              rep%ngwfs_on_grid,ngwf_basis,proj_basis,nlps_projectors)
      end if

      if (pub_hubbard) then
         ! ddor: recompute the DFT+U on-site ngwf-projector overlap matrix
         !       since the NGWFs have been modified.
         call projectors_func_ovlp_box(rep%hub_overlap, &                 ! out
              rep%ngwfs_on_grid,ngwf_basis,hub_proj_basis,hub%projectors) ! in
      endif

      call sparse_destroy(kernel)
      call sparse_destroy(h_aoh)
      call sparse_destroy(ao_to_hy)
      call sparse_destroy(trans_ao_to_hy)
      call sparse_destroy(hy_to_ao)
      call sparse_destroy(trans_hy_to_ao)
      call sparse_destroy(k_hao)

    end subroutine internal_hybridize_ngwfs

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_qc_output

      !====================================================!
      ! This subroutine prints out quality control info in !
      ! a form that can be easily accessed and compared.   !
      !----------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 12/03/2005.    !
      !====================================================!

      use comms, only: comms_barrier
      use constants, only: stdout

      if (pub_on_root) then
         write(stdout,'(a1)')' '
         write(stdout, '(a30,   i9)')'<QC>        [NGWF iterations]:', &
              iteration
         write(stdout, '(a30,f18.8)')'<QC>           [total_energy]:', &
              total_energy
         write(stdout, '(a30,f18.8)')'<QC>                     [F0]:', &
              F0
         write(stdout, '(a30,f18.8)')'<QC>                  [Fdiis]:', &
              Fdiis
         write(stdout, '(a30,f18.8)')'<QC>                     [F1]:', &
              F1
         write(stdout, '(a30,f18.8)')'<QC>                     [F2]:', &
              F2
         write(stdout, '(a30,f18.8)')'<QC>   [predicted_functional]:', &
              predicted_functional
         write(stdout,'(a30,f22.12)')'<QC>           [rms_gradient]:', &
              rms_gradient
         write(stdout,'(a30,f22.12)')'<QC>                 [G_init]:', &
              G_init
         write(stdout, '(a30,f16.6)')'<QC>        [quadratic_coeff]:', &
              quadratic_coeff
         write(stdout, '(a30,f16.6)')'<QC>            [cubic_coeff]:', &
              cubic_coeff
         write(stdout, '(a30,f16.6)')'<QC>           [trial_length]:', &
              trial_length
      endif

      call comms_barrier


    end subroutine internal_qc_output


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    logical function internal_test_convergence()

      use hubbard_build, only: hubbard_spin_splitting_zero
      use rundat, only: delta_e_conv, maxit_pen, pub_hub_ngwf_spin_thr, &
           pen_param, pub_write_ngwf_grad_plot
      use sparse, only: sparse_is_dense
      use wrappers, only : wrappers_ddot

      implicit none

      logical :: energy_gain_converged
      logical :: ngwf_grad_converged

      ! ddor: Used in the case of DFT+U with self-consistent projectors
      character(len=64) :: our_name, hub_proj_iter_number


      ! calculate rms gradient
      ! cks: parallel calculation of NGWF rms_gradient
      rms_gradient = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
           contra_grad_on_grid,1,cov_grad_on_grid,1)
      call comms_reduce('SUM', rms_gradient)

      ! cks: store for Fletcher-Reeves calculation of CG coeff
      current_g_dot_g = rms_gradient
      rms_gradient = sqrt(abs(rms_gradient) / real(est_num_psincs, kind=DP))

      ! cks: set number of lnv iterations for next NGWF step
      if (rms_gradient > 0.9_DP*previous_rms_gradient) then
         current_maxit_lnv = current_maxit_lnv +1
      end if

      if (current_maxit_lnv > maxit_lnv) current_maxit_lnv =maxit_lnv
      ! cks: store current rms grad for next NGWF step
      previous_rms_gradient = rms_gradient

      ! cks: Test for if allowed to use energy gain as convergence criterion
      energy_gain_converged = .false.
      energy_gain_converged = delta_e_conv .and. &
           ( last_n_energies(3) > last_n_energies(2)) .and. &
           ( last_n_energies(2) > last_n_energies(1))

      ! ndmh: write gradient on grid to file for visualisation
      if (pub_write_ngwf_grad_plot .and. &
           (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
         call visual_ngwfs(cov_grad_on_grid, ngwf_basis, &
              'ngwf_cov_grad', elements)
         call visual_ngwfs(contra_grad_on_grid, ngwf_basis, &
              'ngwf_contra_grad', elements)
      end if

      ! cks: NGWF RMS gradient convergence criterion
      ngwf_grad_converged = (rms_gradient < ngwf_threshold)

      if ( ngwf_grad_converged .or. energy_gain_converged ) then

         ! cks: print details only when output is not set to brief
         if (pub_output_detail >= NORMAL) then

            if (pub_on_root) then
               write(stdout,'(/a)')'           ............................&
                    &............................'
               write(stdout,'(a)')'           | *** Wannier-like function &
                    &optimisation converged *** |'
               write(stdout,'(a,f19.14,a)')'           | RMS NGWF gradient = ',&
                    rms_gradient,'              |'
               write(stdout,'(a)')'           | Criteria satisfied: &
                    &                                 |'
               if (ngwf_grad_converged) write(stdout,'(a)') &
                    '           | -> RMS NGWF gradient lower than set threshold.       |'
               if (energy_gain_converged) then
                  write(stdout,'(a)') &
                       '           | -> Maximum degree of convergence for applied level   |'
                  write(stdout,'(a)') &
                       '           |    of density kernel truncation has been reached.    | '
               endif
               write(stdout,'(a)') '           ===========================&
                    &============================='
            end if
         else
            if (pub_on_root.and.energy_gain_converged) then
               write(stdout,'(a)') &
                    'Maximum degree of convergence for applied level of &
                    &kernel truncation reached.'
            end if
         end if

         ! ndmh: if there is no kernel truncation but energy gain convergence
         ! ndmh: has been set due to rising energy, something must be wrong
         ! ndmh: with the kernel (unless unreasonable demands on NGWF
         ! ndmh: convergence have been requested), so warn user
         if (pub_on_root.and.energy_gain_converged.and. &
              (sparse_is_dense(denskern(1)))) then
            write(stdout,'(a,f8.4,a)')  'WARNING: No kernel truncation, &
                 &yet energy has risen on last 3 iterations.'
            write(stdout,'(a)') 'WARNING: Either kernel occupation numbers &
                 &may be unreliable, or'
            write(stdout,'(a)') 'WARNING: requested NGWF gradient tolerance &
                 &may be unachievable.'
            write(stdout,'(a)') 'WARNING: Check eigenvalues with properties &
                 &calculation if feasible,'
            write(stdout,'(a,f8.4,a,i3)') 'WARNING: else restart calculation &
                 &with pen_param >=', 2.0_DP*pen_param,' and maxit_pen >=', &
                 maxit_pen*2
         end if

         call hamiltonian_energy_components(  &
              pur_denskern, rep, localpseudo_fine, core_density_fine, &
              ngwf_basis, hub_proj_basis, hub, ewald_energy, ham%hfexchange)

         ! cks: write final NGWFs in plotting formats
         if (pub_write_ngwf_plot .and. &
              ( pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then

          ! ddor: If carrying out a DFT+U calculation with self-consistent
          ! ddor: projectors, then write out the final NGWFs on each projector
          ! ddor: optimisation iteration.
          if ( task == 'HUBBARDSCF' ) then
             write(hub_proj_iter_number,*) hub%consistency_iteration
             write(our_name,'(a64)') 'hub_consistency_iteration'// &
             trim(adjustl(hub_proj_iter_number))
             call visual_ngwfs(rep%ngwfs_on_grid,ngwf_basis,our_name, elements)
          else
             call visual_ngwfs(rep%ngwfs_on_grid,ngwf_basis,'final', elements)
          endif

         end if

         internal_test_convergence = .true.
      else
         internal_test_convergence = .false.

         !ddor: Return DFT+U spin-splitting to zero if we have
         !      proceeded far enough in the calculation.
         if ((pub_hubbard) .and. (rms_gradient .lt. pub_hub_ngwf_spin_thr)) then
            call hubbard_spin_splitting_zero(hub,hub_proj_basis)
         endif

      end if

    end function internal_test_convergence


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    ! cks: print a concise summary of the calculation
    subroutine internal_calculation_summary

      implicit none

      integer :: sumrow
      integer :: lastrow

      lastrow = max(iteration+1,2)

      if (pub_on_root) then
         ! ndmh: adapt number of digits in total energy to suit its magnitude
         if(abs(total_energy)<100000_DP) then
            write(summary_lines(lastrow),'(i4, f21.14, f22.14, a)') &
                 iteration, rms_gradient, total_energy, '  <-- DIIS'
         else
            write(summary_lines(lastrow),'(i4, f21.14, f22.12, a)') &
                 iteration, rms_gradient, total_energy, '  <-- DIIS'
         end if

         if (pub_output_detail >= NORMAL) then
            write(stdout,'(/20x,a)') '<<<<< CALCULATION SUMMARY >>>>>'
            do sumrow=1,lastrow
               write(stdout,'(a80)') summary_lines(sumrow)
            end do
         else
            write(stdout,'(a80)') summary_lines(iteration+1)
         endif

      end if

    end subroutine internal_calculation_summary

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_find_direction

      use wrappers, only : wrappers_ddot

      implicit none

      if (.not.line_search_success) then
         trial_length = trial_length * 0.1_DP
      else if (line_search_coeff > 0.0_DP) then
         trial_length = &
              max(sqrt(trial_length * line_search_coeff),epsilon(1.0_DP))
      end if

#ifdef FD
      trial_length = 0.0001_DP
#endif


      ! calculate CG coefficient
      if ( iteration > 1 .and. ngwf_diis_use_cg_dir) then
         if ((cg_count >= pub_elec_cg_max) .or. (.not.line_search_success)) then
            ! cks: reset CG after "cg_max" steps
            ! ndmh: or after a fitting failure
            cg_coeff = 0.0_DP
            cg_count = 0
            if (pub_on_root  .and. (.not.line_search_success) .and. &
                 (pub_output_detail >= NORMAL) ) write(stdout,'(a)') &
                 'NOTE: Resetting NGWF CG directions'
         else
            ! cks <<< FLETCHER >>>
            if (ngwf_cg_type == 'NGWF_FLETCHER') then

               ! cks: original Fletcher-Reeves formula (cheaper in memory)
               if (abs(previous_g_dot_g) > epsilon(1.0_DP)) then
                  cg_coeff = current_g_dot_g / previous_g_dot_g
                  ! cks: protection from crazy coefficients
                  if (abs(cg_coeff) > 2.0_DP ) then
                     if (pub_on_root  .and.(pub_output_detail >= NORMAL) ) &
                          write(stdout,'(a,f8.4,a)') &
                          'WARNING: NGWF Fletcher-Reeves CG coeff too large (',&
                          cg_coeff, ') - setting to zero'
                     cg_coeff = 0.0_DP
                     cg_count = 0
                  end if
               else
                  if (pub_on_root .and.(pub_output_detail >= NORMAL) ) &
                       write(stdout,*) ' NOTE: previous_g_dot_g=', &
                       previous_g_dot_g, &
                       'CG coeffient set to zero in ngwf_diis_optimise'
                  cg_coeff = 0.0_DP
                  cg_count = 0
               end if

               ! cks: <<< POLAK >>>
            else if (ngwf_cg_type == 'NGWF_POLAK') then

               ! POLAK FORMULA
               cg_coeff = services_polak_cg_coeff(prev_direction_on_grid, &
                    cov_grad_on_grid,contra_grad_on_grid, &
                    prev_contra_grad_on_grid, ngwf_basis%n_ppds*pub_cell%n_pts)

               ! FLETCHER-REEVES FORMULA
               !        cg_coeff=services_fr2_cg_coeff(prev_direction_on_grid, &
               !             cov_grad_on_grid,contra_grad_on_grid, &
               !             prev_contra_grad_on_grid,n_ngwf_ppds*cell%n_pts)

            end if


            ! cks: re-initialise the periodic reset process if cg_coeff was zero
            ! cks: otherwise increase cg_count
            if (cg_coeff == 0.0_DP) then
               cg_count = 0
            else
               cg_count = cg_count + 1
            end if

         end if
      else
         cg_coeff = 0.0_DP
      endif

      ! Find search direction
      if (pub_elec_cg_max > 0) then
         direction_on_grid = -cov_grad_on_grid + cg_coeff * &
              prev_direction_on_grid
         contra_direction_on_grid = -cov_grad_on_grid + cg_coeff * &
              prev_contra_direction_on_grid
      else
         ! cks:steepest descents
         direction_on_grid = -cov_grad_on_grid
         contra_direction_on_grid = -contra_grad_on_grid
      endif


      ! Slope of energy in search direction
      G_init = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
           contra_grad_on_grid,1,direction_on_grid,1)

      ! cks: collect the work of each node
      call comms_reduce('SUM', G_init)

      ! take action in case of positive slope along search direction
      if (G_init > 0.0_DP) then

         if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
            write(stdout,'(a,e16.6)') &
                 'WARNING: slope along search direction is positive:', G_init
            write(stdout,'(a)') '         Resetting conjugate gradients!'
         end if
         direction_on_grid = -cov_grad_on_grid
         contra_direction_on_grid = -contra_grad_on_grid

         cg_count = 0
         G_init = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
              contra_grad_on_grid,1,direction_on_grid,1)

         ! cks: collect the work of each node
         call comms_reduce('SUM', G_init)


         if (G_init > 0.0_DP) then
            if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
               write(stdout,'(a)') 'WARNING: slope along search direction is still &
                    &positive.'
               write(stdout,'(a)') '         Reversing search direction!!'
            end if
            direction_on_grid = -direction_on_grid
            contra_direction_on_grid = -contra_direction_on_grid
            G_init = -G_init
            ! ndmh: if searching 'uphill', always re-check final step
            ! ndmh: before accepting, to avoid accepting very bad steps.
            reversing = .true.
         end if

      end if

      call services_flush

    end subroutine internal_find_direction


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_print_search_info

      implicit none

      write(stdout,'(a,f22.14)')   'RMS gradient                = ',rms_gradient
      write(stdout,'(a,f22.6)')    'Trial step length           = ',trial_length
      write(stdout,'(a,f22.11)')   'Gradient along search dir.  = ',G_init
#ifdef FD
      write(stdout,'(a,f22.11)')   'Gradient by FD              = ',&
           (F1-F0)/trial_length
#endif
      write(stdout,'(a,f22.14)')   'Functional at step 0        = ',F0
      write(stdout,'(a,f22.14)')   'Functional at rmm DIIS      = ',Fdiis
      write(stdout,'(a,f22.14)')   'Functional at step 1        = ',F1
      if (trial2) then
         write(stdout,'(a,f22.14)')'Functional at step 2        = ',F2
         write(stdout,'(a,f22.14)')'Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'Rejected quadratic step     = ', &
              quadratic_coeff
         write(stdout,'(a,f22.6)') 'Selected cubic step         = ',cubic_coeff
      else if (retrial1) then
         write(stdout,'(a,f22.6)') 'Rejected quadratic step     = ', &
              rejected_quadratic_coeff
         write(stdout,'(a,f22.14)')'Functional at new step 1    = ',F2
         write(stdout,'(a,f22.14)')'Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'Selected quadratic step     = ', &
              quadratic_coeff
      else
         write(stdout,'(a,f22.14)')'Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'Selected quadratic step     = ', &
              quadratic_coeff
      endif
      write(stdout,'(a)')'--------------------------- NGWF line search &
           &finished --------------------------'
      write(stdout,'(a)')'                                               '

    end subroutine internal_print_search_info


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine ngwf_diis_optimise

end module ngwf_diis
