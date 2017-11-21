! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!            NGWF Conjugate Gradients optimisation module        !
!                                                                !
! This module performs the optimisation of the electronic energy !
! with respect to the NGWFs by Conjugate Gradients Minimisation. !
!----------------------------------------------------------------!
! This module was created by Nicholas Hine on 28/04/2010 out of  !
! part of the previous version of electronic_mod.                !
! Originally written by Chris-Kriton Skylaris in 2000.           !
! Subsequent modifications by Chris-Kriton Skylaris,             !
! Arash A. Mostofi, Peter D. Haynes and Nicholas D.M. Hine.      !
! Modified by Nicholas Hine to use the NGWF_REP container for    !
! NGWF-related data and matrices, October 2010.                  !
!================================================================!


module ngwf_cg

  implicit none

  private

  public :: ngwf_cg_optimise

contains

  subroutine ngwf_cg_optimise(total_energy, converged, &              !   out
       ham, denskern, rep, ngwf_nonsc_forces, lhxc_fine, &            ! inout
       ngwf_basis, proj_basis, hub_proj_basis, hub, &                 ! in
       elements, ewald_energy, localpseudo_fine, core_density_fine, & ! in
       val_rep, val_ngwf_basis, val_dkn, val_ham) ! lr408: optional cond args


    !==========================================================================!
    ! This subroutine minimises the total energy with respect to the NGWF      !
    ! expansion coefficients and the density kernel elements subject to the    !
    ! constraint that the density matrix is idempotent and integrates to       !
    ! the correct number of electrons.                                         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! total_energy    (output): the total energy                               !
    ! converged       (output): whether the total energy was converged         !
    ! rep              (inout): NGWF Representation (functions and matrices)   !
    ! denskern         (inout): density kernel for only alpha (or beta)        !
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
    ! lhxc_fine        (inout): Local-Hartree-Exhange-Correlation potential    !
    ! ngwf_nonsc_forces(inout): Outer loop non self-consistent force correction!
    ! For Conduction NGWF optimisation only:                                  !
    ! val_rep          (input): Valence NGWF representation (optional)         !
    ! val_ngwf_basis   (input): Valence NGWF basis (optional)                  !
    ! val_dkn          (input): Valence density kernel (optional)              !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2001.                           !
    ! Modified by Arash Mostofi.                                               !
    ! Modified by Chris-Kriton Skylaris in June 2003.                          !
    ! Modified by Chris-Kriton Skylaris October-December 2003 so that          !
    ! it runs in parallel.                                                     !
    ! Modified by Chris-Kriton Skylaris on 22/02/2004 to reduce                !
    ! memory requirements.                                                     !
    ! Modified by Peter Haynes on 1/07/2004 to use fourier parallelisation     !
    ! Fix of bug in calculation summary by Chris-Kriton Skylaris on 16/07/2004.!
    ! Addition of qc printout and improvements in tests for energy convergence !
    ! by Chris-Kriton Skylaris on 12/03/2005.                                  !
    ! Modification to increase number of lnv iterations when NGWF convergence  !
    ! stagnates by Chris-Kriton Skylaris on 22/03/2005.                        !
    ! Modified for spin polarisation by Peter Haynes, July 2006                !
    ! Modified for Nonlinear Core Corrections by Nicholas Hine, January 2009   !
    ! Modified by David O'Regan for DFT+U, April 2009                          !
    ! Adapted for SPAM3 and function basis by Nicholas Hine, May-July 2009.    !
    ! Renamed from electronic_energy_minimise_tc to ngwf_cg_optimise and moved !
    ! to new module by Nicholas Hine on 28/04/2010.                            !
    ! Modifications for NGWF_REP and NGWF_HAM types by Nicholas Hine in        !
    ! October 2010.                                                            !
    ! Modified by Laura Ratcliff for conduction calculations, Oct 2010.        !
    ! Modified for Christoffel Kernel updates by David O'Regan and Nicholas    !
    ! Hine in November 2010.                                                   !
    ! New convergence criteria by Nicholas Hine, 11/11/2010.                   !
    ! Non self-consistent forces correction by Alvaro Ruiz Serrano, 19/11/2010.!
    ! Return flag indicating convergence added by Nicholas Hine 20/04/2011.    !
    !==========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: comms_abort, pub_on_root, pub_root_node_id, comms_bcast, &
         comms_reduce
    use constants, only: DP, verbose, max_spins, normal, stdout, brief
    use electronic, only: electronic_energy, electronic_lagrangian
    use function_basis, only: FUNC_BASIS, function_basis_est_num_psincs
    use geometry, only: point
    use hamiltonian, only: hamiltonian_build_matrix, &
         hamiltonian_dens_indep_matrices, hamiltonian_energy_components, &
         hamiltonian_dens_dep_nonsc
    use hubbard_build, only: HUBBARD_MODEL, hubbard_energy_info
    use ion, only: element
    use kernel, only: kernel_normalise
    use ngwf_gradient, only: ngwf_gradient_lnv
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_ham_create, &
         ngwf_ham_destroy
    use nonsc_forces, only: nonsc_forces_ngwfs_calc
    use restart, only: restart_kernel_write, restart_ngwfs_tightbox_output, &
         restart_sph_waves_output, restart_kernel_store, &
         restart_ngwfs_tightbox_store, store_tightbox_ngwfs, store_denskern, &
         retrieve_tightbox_ngwfs, retrieve_denskern
    use rundat, only: pub_kernel_update, write_tightbox_ngwfs, ngwf_cg_type, &
         pub_output_detail, pub_nnho, precond_real, precond_recip, &
         maxit_ngwf_cg, minit_lnv, ngwf_threshold_orig, pub_elec_cg_max, &
         pub_usehfx, pub_write_ngwf_plot, pub_cube_format, pub_dx_format, &
         pub_grd_format, pub_hubbard, print_qc, k_zero, lnv_threshold_orig, &
         maxit_lnv, pub_write_sw_ngwfs, task, ngwf_cg_max_step, &
         write_denskern, pub_write_converged_dk_ngwfs, pub_aug, &
         pub_cond_calculate, pub_kernel_christoffel_update, pub_elec_force_tol,&
         pub_write_forces, pub_nonsc_forces, pub_rootname, pub_devel_code, &
         pub_write_ngwf_radial, pub_write_ngwf_grad_radial
!CW
    use rundat,only : pub_dmft_spoil_kernel,pub_dmft_fully_sc,pub_dmft_fully_sc_h
!END CW
    use services, only: services_flush, services_polak_cg_coeff, &
         services_cubic_fit_minimum, services_line_search_parabola
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale, &
         sparse_copy, sparse_max_abs_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_ngwfs, visual_ngwfs_radial
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
    type(NGWF_HAM), intent(inout) :: ham
    ! ars: non self-consistent forces correction due to NGWF loop
    real(kind=DP), intent(inout) :: ngwf_nonsc_forces(:,:)
    logical, intent(out) :: converged  ! pdh: convergence flag

    ! lr408: Optional arguments needed for conduction calculation
    type(NGWF_REP), optional, intent(in) :: val_rep
    type(NGWF_HAM), optional, intent(in) :: val_ham
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis
    type(SPAM3), optional, intent(in)      :: val_dkn(pub_cell%num_spins)

    ! Local Variables
    type(SPAM3), allocatable :: pur_denskern(:)
    type(SPAM3), allocatable :: start_denskern(:)
    type(SPAM3) :: start_sp_overlap
    real(kind=DP), allocatable, dimension(:) :: cov_grad_on_grid
    real(kind=DP), allocatable, dimension(:) :: prev_contra_grad_on_grid
    real(kind=DP), allocatable, dimension(:) :: contra_grad_on_grid
    real(kind=DP), allocatable, dimension(:) :: direction_on_grid
    real(kind=DP), allocatable, dimension(:) :: start_ngwfs_on_grid
    real(kind=DP), allocatable, dimension(:) :: prev_direction_on_grid
    real(kind=DP), allocatable :: total_forces(:,:,:)
    real(kind=DP) :: last_n_energies(3) ! space to store the 3 most recent energies
    real(kind=DP) :: mu(max_spins)
    real(kind=DP) :: rms_gradient
    real(kind=DP) :: max_gradient
    real(kind=DP) :: previous_rms_gradient ! RMS NGWF grad of previous iteration
    real(kind=DP) :: line_search_coeff
    real(kind=DP) :: lnv_threshold
    real(kind=DP) :: ngwf_threshold
    real(kind=DP) :: F0,F1,F2,G_init,trial_length
    real(kind=DP) :: cg_coeff
    real(kind=DP) :: quadratic_coeff,cubic_coeff,rejected_quadratic_coeff
    real(kind=DP) :: predicted_functional
    real(kind=DP) :: previous_g_dot_g
    real(kind=DP) :: current_g_dot_g
    real(kind=DP) :: previous_dir_dot_g
    character(len=80), allocatable, dimension(:) :: summary_lines
    integer :: iteration         ! current iteration
    integer :: cg_count          ! current number of steps since CG reset
    integer :: est_num_psincs    ! estimate number of psincs in all NGWF spheres
    integer :: current_maxit_lnv ! number of lnv iterations to do for current NGWF step
    integer :: is         ! pdh: spin loop counter
    integer :: ierr       ! error flag
    integer :: minit      ! qoh: Minimum number of CG iterations
    logical :: trial2     ! ndmh: flag to perform second trial step
    logical :: retrial1   ! ndmh: flag to perform repeat first trial step
    logical :: reversing  ! ndmh: line search is going uphill
    logical :: line_search_success ! ndmh: line search fit success flag
    logical :: check_conjugacy     ! ndmh: flag for doing cg conjugacy check
    logical :: updated_shift ! lr408: Flag needed for conduction calculations

    ! ndmh: FD variables
    integer, parameter :: nfd=25
    integer :: ifd
    real(kind=DP) :: fd_trial_step(nfd)=(/0.000001_DP,0.00001_DP,0.0001_DP, &
         0.001_DP,0.01_DP,0.02_DP,0.04_DP,0.06_DP,0.08_DP,0.1_DP, &
         0.15_DP,0.2_DP,0.3_DP,0.4_DP,0.5_DP,0.6_DP,0.7_DP,0.8_DP,0.9_DP, &
         1.0_DP,1.2_DP,1.4_DP,1.6_DP,1.8_DP,2.0_DP/)
    real(kind=DP) :: FFD(nfd)

    ! cks: FOR PLOTTING FUNCTIONS
    character(len=64) :: fun_name, iter_number

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering ngwf_cg_optimise'
#endif

    ! Flush output
    call services_flush

    ! Start timer
    call timer_clock('ngwf_cg_optimise',1)

    ! cks: write initial NGWFs in plotting formats
    if (pub_write_ngwf_plot .and. &
         (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
       call visual_ngwfs(rep%ngwfs_on_grid, ngwf_basis, 'initial', elements)
    endif

    ! smmd: write initial NGWFs radial distribution
    if (pub_write_ngwf_radial.ge.1) then
       call visual_ngwfs_radial(rep%ngwfs_on_grid, ngwf_basis, 'initial','ngwf')
    endif


    ! pdh: allocate local sparse matrices
    allocate(pur_denskern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('ngwf_cg_optimise','pur_denskern',ierr)
    do is=1,pub_cell%num_spins
       call sparse_create(pur_denskern(is),denskern(is))
    end do

    ! Allocate workspace
    allocate(contra_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_cg_optimise', &
         'contra_grad_on_grid',ierr)
    allocate(direction_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_cg_optimise', &
         'direction_on_grid',ierr)
    allocate(start_ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_cg_optimise', &
         'start_ngwfs_on_grid',ierr)

    if (pub_elec_cg_max > 0) then
       ! cks: allocate properly only when not doing steepest descents
       allocate(prev_direction_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), &
            stat=ierr)
    else
       ! cks: otherwise allocate only token memory to keep compiler happy
       allocate(prev_direction_on_grid(1),stat=ierr)
    end if
    call utils_alloc_check('ngwf_cg_optimise','prev_direction_on_grid',ierr)

    allocate(cov_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('ngwf_cg_optimise','cov_grad_on_grid',ierr)

    if (ngwf_cg_type == 'NGWF_POLAK') then
       allocate(prev_contra_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), &
            stat=ierr)
       call utils_alloc_check('ngwf_cg_optimise','prev_contra_grad_on_grid', &
            ierr)
    end if
    if (pub_kernel_christoffel_update) then
       allocate(start_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('ngwf_cg_optimise','start_denskern',ierr)
       do is=1,pub_cell%num_spins
          call sparse_create(start_denskern(is),denskern(is))
       end do
       if (pub_aug) call sparse_create(start_sp_overlap,rep%sp_overlap)
    end if
    allocate(summary_lines(max(maxit_ngwf_cg+1,2)),stat=ierr)
    call utils_alloc_check('ngwf_cg_optimise','summary_lines',ierr)

    ! ndmh: for force convergence testing
    if (pub_elec_force_tol > 0.0_DP) then
       allocate(total_forces(3,pub_cell%nat,3),stat=ierr)
       call utils_alloc_check('ngwf_cg_optimise','total_forces',ierr)
    end if

    ! cks: <<< parameter initialisations >>>
    ngwf_threshold     = ngwf_threshold_orig
    lnv_threshold      = lnv_threshold_orig
    est_num_psincs     = function_basis_est_num_psincs(ngwf_basis)

    ! cks: <<< variable intialisations >>>
    current_maxit_lnv     = minit_lnv
    previous_rms_gradient = huge(1.0_DP)
    max_gradient          = huge(1.0_DP)
    line_search_coeff  = 0.15_DP
    trial_length       = 0.1_DP
    F0=0.0_DP ; F1=0.0_DP ; F2=0.0_DP ; G_init=0.0_DP ; cg_coeff=0.0_DP
    rms_gradient=1.0_DP ; mu=0.0_DP ; total_energy=0.0_DP ; cg_count=0
    predicted_functional =0.0_DP
    current_g_dot_g  =0.0_DP
    previous_g_dot_g =0.0_DP
    line_search_success = .true.
    trial2 = .false.
    retrial1 = .false.
    reversing = .false.
    prev_direction_on_grid = 0.0_DP
    last_n_energies(:) = huge(1.0_DP)
    quadratic_coeff = 0.0_DP
    rejected_quadratic_coeff = 0.0_DP
    cubic_coeff = 0.0_DP
    converged = .false.

    ! lr408: Default is false
    updated_shift = .false.

    ! ars: nonsc_forces zero by default
    ngwf_nonsc_forces(:,:) = 0.0_DP

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

    ! qoh: Allow blank calculations for timings purposes if maxit_ngwf_cg < 0
    minit = 1
    if (maxit_ngwf_cg < 0) minit = 0

    ! cks: NGWF iterations loop
    iteration = 0 ! ndmh: prevent unitialised variable when maxit_ngwf_cg = 0

    do iteration=1,max(maxit_ngwf_cg,minit)

       ! cks: ++++++++++ ITERATION HEADER +++++++++++++++++++++++++++++
       if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
            maxit_ngwf_cg > 0) then
          write(stdout,'(/a)')'########################################&
               &########################################'
          if (pub_cond_calculate) then
             write(stdout,'(a,i4.3,a)')'######################### COND NGWF &
                  &CG iteration ',iteration,' ##########################'
          else
             write(stdout,'(a,i4.3,a)')'########################### NGWF &
                  &CG iteration ',iteration,' #############################'
          end if
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
       call services_flush
       ! cks: ++++++ END ITERATION HEADER +++++++++++++++++++++++++++++



       ! lr408: Update conduction Hamiltonian
       if (pub_cond_calculate) call hamiltonian_dens_dep_nonsc(ham,rep, &
            ngwf_basis,lhxc_fine,hub,val_rep,val_ham,val_dkn,updated_shift)

!CW
!CW COMMENT : here update density kernel
!END CW
       ! cks: +++++++++++++++++++++++ OPTIMISE DENSITY KERNEL ++++++++++++++++++
           total_energy = electronic_energy(denskern, pur_denskern, ham, &
           lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
           localpseudo_fine, core_density_fine, ewald_energy, elements, &
           lnv_threshold, current_maxit_lnv, kernel_update=.true.)
       ! cks: ++++++++++++++++++++ END OPTIMISE DENSITY KERNEL +++++++++++++++++


       ! cks: =========== HYBRID NGWFs==========================================
       if (pub_nnho) then
          call internal_hybridize_ngwfs
       endif
       ! cks: ======= END HYBRID NGWFs==========================================


       ! lr408: Update Hamiltonian (in conduction case this is unnecessary)
!CW
       if (.not. pub_cond_calculate.and..not.pub_dmft_spoil_kernel) then
          if(.not.pub_dmft_fully_sc.or.(pub_dmft_fully_sc.and.pub_dmft_fully_sc_h)) call hamiltonian_build_matrix(ham, rep)
       end if
!END CW

       ! ars: calculate gradient if nonsc forces required
       if ((maxit_ngwf_cg == 0).and..not.pub_nonsc_forces) then
          converged = .false.
          exit
       end if

       ! cks: ~~~~~~~~~~~~~~~~~~~~~~~ NGWF GRADIENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if (.not.pub_cond_calculate) then
          call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, & ! out
               denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &   ! in
               lhxc_fine, ham, hub, mu, elements)                         ! in
       else
          ! lr408: Need extra arguments for conduction calculation
          call ngwf_gradient_lnv(contra_grad_on_grid, cov_grad_on_grid, & ! out
               denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &   ! in
               lhxc_fine, ham, hub, mu, elements, & ! in
               val_dkn, val_rep, val_ham%ham, val_ngwf_basis,ham%cond_shift)! in
       end if
       ! cks: ~~~~~~~~~~~~~~~~~~~~ END NGWF GRADIENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~


       ! ndmh: check conjugacy condition on search direction
       check_conjugacy = .false.
       if (check_conjugacy) then
          previous_dir_dot_g = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
               contra_grad_on_grid,1,prev_direction_on_grid,1)
          call comms_reduce('SUM',previous_dir_dot_g)
       end if

       ! cks: ###################### TEST CONVERGENCE ##########################
       converged = internal_test_convergence()
       call comms_bcast(pub_root_node_id,converged)
       if (converged) exit
       ! cks: #################### END TEST CONVERGENCE ########################

       ! pdh: exit if no NGWF optimisation required
       if (maxit_ngwf_cg == 0) then
          converged = .false.
          exit
       end if


       ! cks: ===== write out different energy components every so often =======
       ! lr408: Don't print out components if this is a conduction calculation
!CW
       if (((mod(iteration, 5) == 0 .and. pub_output_detail == VERBOSE) .or. &
            iteration == maxit_ngwf_cg).and. (.not.pub_cond_calculate) ) &
            call hamiltonian_energy_components( &
            pur_denskern, rep, localpseudo_fine, core_density_fine, &
            ngwf_basis, hub_proj_basis, hub, ewald_energy, ham%hfexchange)
!END CW
       ! cks: ==================================================================


       ! cks: **************** FUNCTIONAL AT INITIAL POINT *********************
       F0 = electronic_lagrangian(total_energy, rep%overlap, denskern, &
            ham%ham, mu, rep%n_occ)
       ! cks: ************** END FUNCTIONAL AT INITIAL POINT *******************


       ! cks: ************************ LINE SEARCH *****************************
       ! the density kernel is not updated during the line search

       call internal_find_direction

       ! cks: store direction if doing conjugate gradients
       if (pub_elec_cg_max > 0) then
          prev_direction_on_grid = direction_on_grid
       endif

       if (ngwf_cg_type == 'NGWF_POLAK') &
            prev_contra_grad_on_grid = contra_grad_on_grid

       previous_g_dot_g = current_g_dot_g

       start_ngwfs_on_grid = rep%ngwfs_on_grid

       ! ndmh: if doing Christoffel updates to the kernel, save the denskern at
       ! the initial point now
       if (pub_kernel_christoffel_update) then
          do is=1,pub_cell%num_spins
             call sparse_copy(start_denskern(is),denskern(is))
             if (pub_aug) call sparse_copy(start_sp_overlap,rep%sp_overlap)
          end do
       end if

       if (index(pub_devel_code,'NGWF_FD')>0) then
          do ifd=1,nfd
             rep%ngwfs_on_grid = start_ngwfs_on_grid + &
                  fd_trial_step(ifd)*direction_on_grid

             ! ndmh: %%%%%%%%%%%%% FUNCTIONAL AT FD STEP (DEBUG) %%%%%%%%%%%%%%%
             ! ndmh: calculate matrix elements with NGWFs at FD step
             call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
                  hub_proj_basis, hub, val_rep, val_ngwf_basis) ! 2x optional args

             if (pub_kernel_christoffel_update) &
                  call internal_kernel_christoffel( &
                  rep%ngwfs_on_grid,start_ngwfs_on_grid,denskern, &
                  start_denskern,start_sp_overlap)

             ! lr408: Update conduction Hamiltonian
             if (pub_cond_calculate) call hamiltonian_dens_dep_nonsc(ham,rep, &
                  ngwf_basis,lhxc_fine,hub,val_rep,val_ham,val_dkn, &
                  updated_shift)
             FFD(ifd) = electronic_energy(denskern, pur_denskern, ham, &
                  lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
                  localpseudo_fine, core_density_fine, ewald_energy, elements, &
                  lnv_threshold, current_maxit_lnv, pub_kernel_update)

             FFD(ifd) = electronic_lagrangian(FFD(ifd),rep%overlap,denskern, &
                  ham%ham,mu,rep%n_occ)

             if ((pub_output_detail>=VERBOSE).and.(.not.pub_cond_calculate)) then
                call hamiltonian_energy_components(  &
                     pur_denskern, rep, localpseudo_fine, core_density_fine, &
                     ngwf_basis, hub_proj_basis, hub, ewald_energy, &
                     ham%hfexchange)
             end if
             ! ndmh: %%%%%%%%%%%%% FUNCTIONAL AT FD STEP (DEBUG) %%%%%%%%%%%%%%%
          end do
       end if

       rep%ngwfs_on_grid = start_ngwfs_on_grid + trial_length*direction_on_grid

       ! cks: %%%%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 1 %%%%%%%%%%%%%%%%%%
       ! cks: calculate matrix elements with trial NGWFs
       call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
            hub_proj_basis, hub, val_rep, val_ngwf_basis) ! 2x optional args

       if (pub_kernel_christoffel_update) call internal_kernel_christoffel( &
            rep%ngwfs_on_grid,start_ngwfs_on_grid,denskern,start_denskern, &
            start_sp_overlap)

       ! lr408: Update conduction Hamiltonian
       if (pub_cond_calculate) call hamiltonian_dens_dep_nonsc(ham,rep, &
            ngwf_basis,lhxc_fine,hub,val_rep,val_ham,val_dkn,updated_shift)
       F1 = electronic_energy(denskern, pur_denskern, ham, &
            lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
            localpseudo_fine, core_density_fine, ewald_energy, elements, &
            lnv_threshold, current_maxit_lnv, pub_kernel_update)

       F1 = electronic_lagrangian(F1,rep%overlap,denskern,ham%ham,mu,rep%n_occ)

       ! cks: %%%%%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 1 %%%%%%%%%%%%%%%%


       ! cks: ===================== SEARCH BY FITTING PARABOLA =================
       call services_line_search_parabola(&
            quadratic_coeff, predicted_functional,line_search_success, & !output
            G_init, F0, F1, trial_length, ngwf_cg_max_step)              !input

       line_search_coeff = quadratic_coeff

       ! cks: ================= END SEARCH BY FITTING PARABOLA =================


       ! cks: CUBIC CUBIC CUBIC ----- SEARCH BY FITTING CUBIC ---- CUBIC CUBIC
       if (G_init * quadratic_coeff > 0.0_DP) then

          trial2 = .true.

          rep%ngwfs_on_grid = start_ngwfs_on_grid + &
               2.0_DP * trial_length * direction_on_grid
       end if

       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP
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
          rep%ngwfs_on_grid = start_ngwfs_on_grid + &
               quadratic_coeff * direction_on_grid
       end if
       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP

       if (trial2 .or. retrial1) then

          ! cks: %%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%

          call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
               hub_proj_basis, hub, val_rep, val_ngwf_basis) ! 2x optional args

          if (pub_kernel_christoffel_update) call internal_kernel_christoffel( &
               rep%ngwfs_on_grid,start_ngwfs_on_grid,denskern,start_denskern, &
               start_sp_overlap)

          ! lr408: Update conduction Hamiltonian

          if (pub_cond_calculate) call hamiltonian_dens_dep_nonsc(ham,rep, &
               ngwf_basis,lhxc_fine,hub,val_rep,val_ham,val_dkn,updated_shift)

          F2 = electronic_energy(denskern, pur_denskern, ham, &
               lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
               localpseudo_fine, core_density_fine, ewald_energy, elements, &
               lnv_threshold, current_maxit_lnv, pub_kernel_update)


          F2 = electronic_lagrangian(F2,rep%overlap,denskern,ham%ham,mu,rep%n_occ)

          ! cks: %%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%

          if (trial2) then
             ! Cubic Fit
             call services_cubic_fit_minimum( &
                  cubic_coeff, predicted_functional, line_search_success, & !output
                  F0, F1, F2, G_init, trial_length, 2.0_DP*trial_length, &  !input
                  ngwf_cg_max_step)                                         !input

             line_search_coeff = cubic_coeff
          end if

          if (retrial1) then
             ! ndmh: quadratic fit at new trial length
             call services_line_search_parabola(&
                  quadratic_coeff, predicted_functional, &               !output
                  line_search_success, G_init, F0, F2, &                 !input
                  line_search_coeff, ngwf_cg_max_step)                   !input

             line_search_coeff = quadratic_coeff
          end if

          ! AAM: quadratic fit using two trial steps
          !          quadratic_coeff=services_parabolic_step( &
          !               F0,F1,F2,trial_length,2.0_DP*trial_length)
          !          line_search_coeff=quadratic_coeff
       else
          ! ndmh: no second trial step
          F2 = 0.0_DP
       end if
       ! cks: CUBIC CUBIC --------- END SEARCH BY FITTING CUBIC -------- CUBIC


       ! cks: set the new values of the NGWFs
       rep%ngwfs_on_grid = start_ngwfs_on_grid + &
            line_search_coeff * direction_on_grid

       ! ndmh: re-calculate density-independent terms in hamiltonian
       call hamiltonian_dens_indep_matrices(rep, &
            ngwf_basis, proj_basis, hub_proj_basis, hub, &
            val_rep, val_ngwf_basis)

       if (pub_kernel_christoffel_update) call internal_kernel_christoffel( &
            rep%ngwfs_on_grid,start_ngwfs_on_grid,denskern,start_denskern, &
            start_sp_overlap)

       if (pub_on_root .and. pub_output_detail >= NORMAL) &
            call internal_print_search_info

       ! cks: write new NGWFs to file if required
       if (write_tightbox_ngwfs.and.(.not.pub_write_converged_dk_ngwfs)) then
          ! cks: in universal tightbox representation
          call restart_ngwfs_tightbox_output(rep%ngwfs_on_grid,ngwf_basis,&
               elements, 'tightbox_'//ngwf_basis%name)
       endif
       ! cks: ************************* END LINE SEARCH ***********************

       ! ndmh: reset flags
       trial2 = .false.
       retrial1 = .false.
       reversing = .false.
       F2 = 0.0_DP

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


    end do


    ! ars: calculate ngwf_nonsc_forces if required
    if (pub_nonsc_forces .and. (pub_write_forces .or. &
         task .eq. 'GEOMETRYOPTIMIZATION' .or. &
         task .eq. 'TRANSITIONSTATESEARCH'.or. &
         task .eq. 'MOLECULARDYNAMICS' .or. &
         task .eq. 'PHONON')) then
       call nonsc_forces_ngwfs_calc(ngwf_nonsc_forces,&
            rep%ngwfs_on_grid, contra_grad_on_grid, ngwf_basis)
    end if

    ! write NGWFs if required and minimisation ends in first iteration
    ! ndmh: or if writing final converged results only
    if (write_tightbox_ngwfs .and. ((iteration.eq.1).or. &
         pub_write_converged_dk_ngwfs)) then
       call restart_ngwfs_tightbox_output(rep%ngwfs_on_grid,ngwf_basis,&
            elements, 'tightbox_'//ngwf_basis%name)
    endif

    ! smmd: store NGWFs for extrapolation at subsequent MD/geom steps
    if (store_tightbox_ngwfs) &
       call restart_ngwfs_tightbox_store(rep%ngwfs_on_grid,ngwf_basis,&
            elements)

    ! ndmh: write denskern here if we have not been writing it as we go along
    if (write_denskern .and. pub_write_converged_dk_ngwfs) &
         call restart_kernel_write(denskern,write_cond=pub_cond_calculate)

    ! smmd: store denskern for extrapolation at subsequent MD/geom steps
    if (store_denskern) &
       call restart_kernel_store(denskern)

    ! ars: write in spherical waves representation
    if (pub_write_sw_ngwfs) then! .and. iteration.eq.1) then
       call restart_sph_waves_output(rep%ngwfs_on_grid,ngwf_basis,elements, &
            'sw_ngwfs')
    endif

    if (.not. converged) then
       ! cks: reset iteration number just for storing final line of calculation
       !      summary
       iteration = iteration - 1
       ! cks: print warning that calculation failed to converge
       if (pub_on_root) write(stdout,'(a,i4,a)') &
            'WARNING: maximum number of NGWF CG iterations (',maxit_ngwf_cg, &
            ') exceeded!'
    end if

    ! cks: print calculation summary
    call internal_calculation_summary

    ! cks: print quality control information
    if (print_qc) call internal_qc_output

    ! ddor: Write out DFT+U occupancies if it hasn't already been done
    if (pub_hubbard.and.(pub_output_detail .ne. VERBOSE)) then
       call hubbard_energy_info(hub,hub_proj_basis)
    endif

    ! ndmh: for force convergence testing
    if (pub_elec_force_tol > 0.0_DP) then
       deallocate(total_forces,stat=ierr)
       call utils_dealloc_check('ngwf_cg_optimise','total_forces',ierr)
    end if
    deallocate(summary_lines,stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','summary_lines',ierr)
    if (pub_kernel_christoffel_update) then
       if (pub_aug) call sparse_destroy(start_sp_overlap)
       do is=pub_cell%num_spins,1,-1
          call sparse_destroy(start_denskern(is))
       end do
       deallocate(start_denskern,stat=ierr)
       call utils_dealloc_check('ngwf_cg_optimise','start_denskern',ierr)
    end if
    if (ngwf_cg_type == 'NGWF_POLAK') then
       deallocate(prev_contra_grad_on_grid,stat=ierr)
       call utils_dealloc_check('ngwf_cg_optimise','prev_contra_grad_on_grid', &
            ierr)
    end if
    deallocate(cov_grad_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','cov_grad_on_grid',ierr)
    deallocate(prev_direction_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','prev_direction_on_grid',ierr)
    deallocate(start_ngwfs_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','start_ngwfs_on_grid',ierr)
    deallocate(direction_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','direction_on_grid',ierr)
    deallocate(contra_grad_on_grid,stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','contra_grad_on_grid',ierr)

    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(pur_denskern(is))
    end do
    deallocate(pur_denskern,stat=ierr)
    call utils_dealloc_check('ngwf_cg_optimise','pur_denskern',ierr)

    call timer_clock('ngwf_cg_optimise',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving ngwf_cg_optimise'
#endif

    ! Flush output
    call services_flush

    return

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_transform_kernel(new_ngwfs_on_grid,&
         old_ngwfs_on_grid, old_sp_overlap)

      !==============================================================!
      ! Transforms the density kernel to its representation in terms !
      ! of new NGWFs                                                 !
      !--------------------------------------------------------------!
      ! Written by Nicholas Hine on 22/10/09                         !
      ! PAW added by David O'Regan 29/3/11                           !
      !==============================================================!

      use augmentation, only: augmentation_overlap
      use integrals, only: integrals_brappd_ketppd
      use rundat, only: pub_aug
      use sparse, only: SPAM3, sparse_create,  sparse_transpose, &
           sparse_product, sparse_destroy, sparse_copy, sparse_axpy

      implicit none

      ! Arguments
      real(kind=DP), intent(in) :: new_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      real(kind=DP), intent(in) :: old_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      type(SPAM3), intent(in) :: old_sp_overlap

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

      if (pub_aug) then
         ! ddor: Calculate the augmentation of the overlap matrix
         call augmentation_overlap(overlap_old_new,old_sp_overlap, &
              rep%sp_overlap)
      end if

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

    subroutine internal_kernel_christoffel(new_ngwfs_on_grid, &
         old_ngwfs_on_grid,new_denskern,old_denskern,old_sp_overlap)

      !==============================================================!
      ! Updates the density kernel with terms from the Christoffel   !
      ! symbols which preserve the completeness of the basis.        !
      !--------------------------------------------------------------!
      ! Written by David O'Regan in June 2010                        !
      ! Modifications by Nicholas Hine, November 2010.               !
      ! PAW added by David O'Regan 29/3/11                           !
      !==============================================================!

      use augmentation, only: augmentation_overlap
      use integrals, only: integrals_brappd_ketppd
      use rundat, only: pub_aug
      use sparse, only: SPAM3, sparse_create,  sparse_transpose, &
           sparse_product, sparse_destroy, sparse_axpy

      implicit none

      ! Arguments
      real(kind=DP), intent(in) :: new_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      real(kind=DP), intent(in) :: old_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      type(SPAM3), intent(inout) :: new_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: old_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: old_sp_overlap

      ! Local Variables
      type(SPAM3) :: step, tc_step, step_sp_overlap
      type(SPAM3) :: ktmp
      real(kind=DP), allocatable :: step_on_grid(:)

      allocate(step_on_grid(ngwf_basis%n_ppds * pub_cell%n_pts),stat=ierr)
      call utils_alloc_check('internal_kernel_christoffel',&
           &'step_on_grid',ierr)

      ! Create temporary matrix structures
      call sparse_create(step,rep%overlap)
      call sparse_create(tc_step,rep%overlap,rep%inv_overlap)
      call sparse_create(ktmp,old_denskern(1))

      ! Calculate change |dphi> in NGWFs
      step_on_grid = new_ngwfs_on_grid - old_ngwfs_on_grid

      ! Calculate overlap of the change with the old NGWFs
      call integrals_brappd_ketppd(step, &
           step_on_grid,ngwf_basis,old_ngwfs_on_grid,ngwf_basis)

      ! ndmh: Apply the augmentated overlap operator to this matrix
      if (pub_aug) then

         ! ndmh: <proj_i|dphi_a> = <proj_i|new_phi_a> - <proj_i|old_phi_a>
         call sparse_create(step_sp_overlap,old_sp_overlap)
         call sparse_copy(step_sp_overlap,rep%sp_overlap)
         call sparse_axpy(step_sp_overlap,old_sp_overlap,-1.0_DP)

        ! ndmh: Calculate the augmentation of the "step" matrix
         call augmentation_overlap(step,old_sp_overlap,step_sp_overlap)

         ! ndmh: Clean up <proj_i|dphi_a> (no longer needed)
         call sparse_destroy(step_sp_overlap)

      end if

      ! Raise an index in the step-NGWF overlap
      call sparse_product(tc_step,step,rep%inv_overlap)

      do is=1,pub_cell%num_spins
         ! Compute -K <g|phi> S^^
         call sparse_product(new_denskern(is),old_denskern(is),tc_step)
         call sparse_scale(new_denskern(is),-1.0_DP)
         ! Transpose to get -S^^ <phi|g> K
         call sparse_transpose(ktmp,new_denskern(is))
         ! Add this to what we got before
         call sparse_axpy(new_denskern(is),ktmp,1.0_DP)
         ! Add the previous kernel to the calculated change to get new kernel
         call sparse_axpy(new_denskern(is),old_denskern(is),1.0_DP)
      enddo

      ! Destroy temporary matrix structures
      call sparse_destroy(ktmp)
      call sparse_destroy(tc_step)
      call sparse_destroy(step)

      deallocate(step_on_grid,stat=ierr)
      call utils_dealloc_check('internal_kernel_christoffel',&
           &'step_on_grid',ierr)

    end subroutine internal_kernel_christoffel

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

      call sparse_create(k_hao, denskern(1))
      call sparse_create(h_aoh, ham%ham(1))

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
      if (pub_aug) then
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
         write(stdout, '(a30,f16.6)')'<QC>               [cg_coeff]:', &
              cg_coeff
         write(stdout, '(a30,f16.6)')'<QC>           [trial_length]:', &
              trial_length
      endif

      call comms_barrier


    end subroutine internal_qc_output


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    logical function internal_test_convergence()

      use forces, only: forces_calculate
      use hubbard_build, only: hubbard_spin_splitting_zero
      use rundat, only: delta_e_conv, maxit_pen, pub_hub_ngwf_spin_thr, &
           pen_param, pub_write_ngwf_grad_plot, pub_elec_energy_tol, &
           pub_ngwf_max_grad, pub_write_forces
      use sparse, only: sparse_is_dense
      use wrappers, only : wrappers_ddot

      implicit none

      logical :: energy_increasing
      logical :: energy_change_converged
      logical :: ngwf_grad_converged
      logical :: ngwf_max_grad_converged
      logical :: forces_converged
      logical :: all_converged
      logical :: orig_write_forces
      logical :: test_e_conv, test_f_conv
      logical :: test_rms_grad, test_max_grad
      real(kind=DP) :: max_force_diff
      real(kind=DP) :: max_energy_diff
      real(kind=DP) :: diff12, diff23
      integer :: iat
      integer :: ipt

      ! ddor: Used in the case of DFT+U with self-consistent projectors
      character(len=64) :: our_name, hub_proj_iter_number

      ! ndmh: If required, recalculate forces to see if they are converged
      forces_converged = .false.
!CW
      test_f_conv = (pub_elec_force_tol > 0.0_DP).and.(iteration>1).and. &
           (.not.pub_cond_calculate).and..not.pub_dmft_spoil_kernel
      if ((pub_elec_force_tol > 0.0_DP).and.(.not.pub_cond_calculate).and..not.pub_dmft_spoil_kernel) then
!END CW
         ! ndmh: Update storage of most recent three forces
         total_forces(:,:,1) = total_forces(:,:,2)
         total_forces(:,:,2) = total_forces(:,:,3)

         ! ndmh: temporarily override write_forces (we do not want to display
         ! ndmh: the forces here every iteration)
         orig_write_forces = pub_write_forces
         pub_write_forces = .false.

         ! ars: calculate nonsc_forces if required
         if (pub_nonsc_forces) then
            call nonsc_forces_ngwfs_calc(ngwf_nonsc_forces,&
                 rep%ngwfs_on_grid, contra_grad_on_grid, ngwf_basis)
         end if

         call forces_calculate(total_forces(:,:,3),denskern,ham,lhxc_fine, &
              rep,ngwf_basis,proj_basis,hub_proj_basis,hub,localpseudo_fine, &
              core_density_fine,elements, ngwf_nonsc_forces)
         pub_write_forces = orig_write_forces

         diff12 = 0.0_DP
         diff23 = 0.0_DP
         ! Loop over atoms to find greatest change in total_forces
         do iat=1,pub_cell%nat
            diff12 = max(diff12, &
                 sqrt(sum((total_forces(:,iat,1)-total_forces(:,iat,2))**2)))
            diff23 = max(diff23, &
                 sqrt(sum((total_forces(:,iat,2)-total_forces(:,iat,3))**2)))
         end do
         if (iteration<2) diff12 = 0.0_DP
         max_force_diff = max(diff12,diff23)

         if ((diff23<pub_elec_force_tol).and.(diff12<pub_elec_force_tol)) then
            forces_converged = .true.
         end if

      end if

      ! cks: Update storage of most recent three energies
      last_n_energies(1) = last_n_energies(2)
      last_n_energies(2) = last_n_energies(3)
      last_n_energies(3) = total_energy

      ! ndmh: Check energy change per atom vs tolerance
      energy_change_converged = .false.
      test_e_conv = (pub_elec_energy_tol > 0.0_DP).and.(iteration>1)
      if (test_e_conv) then
         diff12 = abs(last_n_energies(1)-last_n_energies(2)) &
              / real(pub_cell%nat,kind=DP)
         diff23 = abs(last_n_energies(2)-last_n_energies(3)) &
              / real(pub_cell%nat,kind=DP)
         max_energy_diff = max(diff12,diff23)
         if (max_energy_diff<pub_elec_energy_tol) then
            energy_change_converged = .true.
         end if
      end if

      ! calculate rms gradient
      ! cks: parallel calculation of NGWF rms_gradient
      rms_gradient = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
           contra_grad_on_grid,1,cov_grad_on_grid,1)
      call comms_reduce('SUM', rms_gradient)

      ! cks: store for Fletcher-Reeves calculation of CG coeff
      current_g_dot_g = rms_gradient
      rms_gradient = sqrt(abs(rms_gradient) / real(est_num_psincs, kind=DP))

      ! ndmh: Calculate max value of <g^a|g_a> over all NGWFs on all nodes
      max_gradient = 0.0_DP
      test_max_grad = (pub_ngwf_max_grad > 0.0_DP)
      if (test_max_grad) then
         do ipt=1,ngwf_basis%n_ppds*pub_cell%n_pts
            max_gradient = max(max_gradient, &
                 sqrt(abs(contra_grad_on_grid(ipt)*cov_grad_on_grid(ipt))))
         end do
         call comms_reduce('MAX', max_gradient)
      end if

      ! cks: Test for if allowed to use energy gain as convergence criterion
      energy_increasing = .false.
      energy_increasing = delta_e_conv .and. &
           ( last_n_energies(3) > last_n_energies(2)) .and. &
           ( last_n_energies(2) > last_n_energies(1))

      ! cks: NGWF RMS gradient convergence criterion
      test_rms_grad = (ngwf_threshold > 0.0_DP)
      if (test_rms_grad) then
         ngwf_grad_converged = (rms_gradient < ngwf_threshold)
      else
         ngwf_grad_converged = .false.
      end if

      ! ndmh: NGWF MAX gradient convergence criterion
      ngwf_max_grad_converged = .false.
      if (test_max_grad) then
         ngwf_max_grad_converged = (max_gradient < pub_ngwf_max_grad)
      end if

      ! ndmh: Check all relevant criteria for convergence
      all_converged = .true.
      if ((pub_elec_force_tol > 0.0_DP).and.(iteration<2)) &
           all_converged = .false.
      if ((pub_elec_energy_tol > 0.0_DP).and.(iteration<2)) &
           all_converged = .false.
      if (test_f_conv) &
           all_converged = all_converged.and.forces_converged
      if (test_e_conv) &
           all_converged = all_converged.and.energy_change_converged
      if (test_rms_grad) &
           all_converged = all_converged.and.ngwf_grad_converged
      if (test_max_grad) &
           all_converged = all_converged.and.ngwf_max_grad_converged
      if (delta_e_conv) &
           all_converged = all_converged.or.energy_increasing

      if (all_converged) then

         ! cks: print details only when output is not set to brief
         if (pub_output_detail >= NORMAL) then

            if (pub_on_root) then
               write(stdout,'(/a)')'           ............................&
                    &............................'
               write(stdout,'(a)')'           | *** NGWF optimisation &
                    &converged ***                  |'
               if (test_rms_grad) write(stdout,'(a,f19.14,a)') &
                    '           | RMS NGWF gradient = ',rms_gradient,'              |'
               if (test_max_grad) write(stdout,'(a,f19.14,a)') &
                    '           | MAX NGWF gradient = ',max_gradient,'              |'
               if (test_f_conv) write(stdout,'(a,f17.14,a)') &
                    '           | Maximum change in forces = ',max_force_diff,' au      |'
               if (test_e_conv) write(stdout,'(a,f17.14,a)') &
                    '           | Maximum change in energy = ',max_energy_diff,' / atom  |'
               write(stdout,'(a)')'           | Criteria satisfied: &
                    &                                 |'
               if (test_rms_grad.and.ngwf_grad_converged) write(stdout,'(a)') &
                    '           | -> RMS NGWF gradient lower than set threshold.       |'
               if (test_max_grad.and.ngwf_max_grad_converged) write(stdout,'(a)') &
                    '           | -> MAX NGWF gradient lower than set threshold.       |'
               if (test_e_conv.and.energy_change_converged) write(stdout,'(a)') &
                    '           | -> Energy change per atom lower than set threshold.  |'
               if (test_f_conv.and.forces_converged) write(stdout,'(a)') &
                    '           | -> Maximum force change lower than set threshold.    |'
               if (energy_increasing) then
                  write(stdout,'(a)') &
                       '           | -> Maximum degree of convergence for applied level   |'
                  write(stdout,'(a)') &
                       '           |    of density kernel truncation has been reached.    | '
               endif
               write(stdout,'(a)') '           ===========================&
                    &============================='
            end if
         end if

         ! ndmh: if there is no kernel truncation but energy gain convergence
         ! ndmh: has been set due to rising energy, something must be wrong
         ! ndmh: with the kernel (unless unreasonable demands on NGWF
         ! ndmh: convergence have been requested), so warn user
         if (pub_on_root.and.energy_increasing.and. &
              (sparse_is_dense(denskern(1)))) then
            write(stdout,'(a,f8.4,a)')  'WARNING: No kernel truncation, &
                 &yet energy has risen on last 2 iterations.'
            write(stdout,'(a,f8.4,a)')  'WARNING: This most likely indicates, &
                 &a problem with the density kernel.'
            write(stdout,'(a)') 'WARNING: Either kernel occupation numbers &
                 &may be unreliable, or'
            write(stdout,'(a)') 'WARNING: requested NGWF gradient tolerance &
                 &may be unachievable.'
         end if

         if (.not. pub_cond_calculate) then
            call hamiltonian_energy_components(  &
                 pur_denskern, rep, localpseudo_fine, core_density_fine, &
                 ngwf_basis, hub_proj_basis, hub, ewald_energy, ham%hfexchange)
         end if
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

         ! smmd: Write NGWFs radial distribution for final iteration
         if (pub_write_ngwf_radial.ge.1) then
            call visual_ngwfs_radial(rep%ngwfs_on_grid, ngwf_basis, 'final', 'ngwf')
         end if

         ! smmd: Write covariant gradient radial distribution for final iteration
         if (pub_write_ngwf_grad_radial.ge.1) then
            call visual_ngwfs_radial(cov_grad_on_grid, ngwf_basis, 'final', 'cov_grad')
            call visual_ngwfs_radial(contra_grad_on_grid, ngwf_basis, 'final', 'contra_grad') ! ars
         end if

         internal_test_convergence = .true.
      else

         ! ndmh: print details only when output is set to verbose
         if (pub_output_detail >= VERBOSE) then

            if (pub_on_root) then
               write(stdout,'(/a)')'           ............................&
                    &............................'
               write(stdout,'(a)')'           | *** NGWF optimisation &
                    &not converged ***              |'
               if (test_rms_grad) write(stdout,'(a,f19.14,a)') &
                    '           | RMS NGWF gradient = ',rms_gradient,'              |'
               if (test_max_grad) write(stdout,'(a,f19.14,a)') &
                    '           | MAX NGWF gradient = ',max_gradient,'              |'
               if (test_f_conv) write(stdout,'(a,f17.14,a)') &
                    '           | Maximum change in forces = ',max_force_diff,' au      |'
               if (test_e_conv) write(stdout,'(a,f17.14,a)') &
                    '           | Maximum change in energy = ',max_energy_diff,' / atom  |'
               write(stdout,'(a)')'           | Criteria not satisfied: &
                    &                             |'
               if (test_rms_grad.and.(.not.ngwf_grad_converged)) &
                    write(stdout,'(a)') &
                    '           | -> RMS NGWF gradient higher than set threshold.      |'
               if (test_max_grad.and.(.not.ngwf_max_grad_converged)) &
                    write(stdout,'(a)') &
                    '           | -> MAX NGWF gradient higher than set threshold.      |'
               if (test_e_conv.and.(.not.energy_change_converged)) &
                    write(stdout,'(a)') &
                    '           | -> Energy change per atom higher than set threshold. |'
               if (test_f_conv.and.(.not.forces_converged)) &
                    write(stdout,'(a)') &
                    '           | -> Maximum force change higher than set threshold.   |'
               write(stdout,'(a)') '           ===========================&
                    &============================='
            end if

            if (pub_on_root.and.energy_increasing) then
               write(stdout,'(a)') 'Energy is increasing: Maximum degree of &
                    &convergence may have been reached'
               write(stdout,'(a)') 'for applied level of kernel truncation.'
            end if
         end if

         internal_test_convergence = .false.
      end if

      ! cks: set number of lnv iterations for next NGWF step
      if (rms_gradient > 0.9_DP*previous_rms_gradient) then
         current_maxit_lnv = current_maxit_lnv +1
      end if

      if (current_maxit_lnv > maxit_lnv) current_maxit_lnv = maxit_lnv
      ! cks: store current rms grad for next NGWF step
      previous_rms_gradient = rms_gradient

      if (.not.all_converged) then
         !ddor: Return DFT+U spin-splitting to zero if we have
         !      proceeded far enough in the calculation.
         if ((pub_hubbard) .and. (rms_gradient .lt. pub_hub_ngwf_spin_thr)) then
            call hubbard_spin_splitting_zero(hub,hub_proj_basis)
         endif
      end if

      ! cks: Write NGWFs of current iteration in plotting formats
      if (pub_write_ngwf_plot .and. &
           (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
         write(iter_number,*) iteration
         write(fun_name,'(a64)') 'iteration'// &
              trim(adjustl(iter_number))
         call visual_ngwfs(rep%ngwfs_on_grid, ngwf_basis, fun_name, elements)
      end if

      ! ndmh: write gradient on grid to file for visualisation
      if (pub_write_ngwf_grad_plot .and. &
           (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
         call visual_ngwfs(cov_grad_on_grid, ngwf_basis, &
              'ngwf_cov_grad', elements)
         call visual_ngwfs(contra_grad_on_grid, ngwf_basis, &
              'ngwf_contra_grad', elements)
      end if

      ! smmd: Write NGWFs radial distribution for current iteration
      if (pub_write_ngwf_radial.ge.2) then
         write(iter_number,*) iteration
         write(fun_name,'(a64)') 'iter'//trim(adjustl(iter_number))
         call visual_ngwfs_radial(rep%ngwfs_on_grid, ngwf_basis, fun_name, 'ngwf')
      end if

      ! smmd: Write covariant gradient radial distribution for current iteration
      if (pub_write_ngwf_grad_radial.ge.2) then
         write(iter_number,*) iteration
         write(fun_name,'(a64)') 'iter'//trim(adjustl(iter_number))
         call visual_ngwfs_radial(cov_grad_on_grid, ngwf_basis, fun_name, 'cov_grad')
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
                 iteration, rms_gradient, total_energy, '  <-- CG'
         else
            write(summary_lines(lastrow),'(i4, f21.14, f22.12, a)') &
                 iteration, rms_gradient, total_energy, '  <-- CG'
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
         trial_length = trial_length * 0.5_DP
      else if (line_search_coeff > 0.0_DP) then
         trial_length = &
              max(sqrt(trial_length * line_search_coeff),epsilon(1.0_DP))
      end if
      trial_length = max(trial_length,0.0001_DP)

      if (pub_on_root .and. (pub_output_detail >= NORMAL) ) write(stdout,'(/a)') &
           '------------------------------- NGWF line search &
           &-------------------------------'


      ! calculate CG coefficient
      if ( iteration > 1 ) then
         if ((cg_count >= pub_elec_cg_max) .or. (.not.line_search_success) &
              .or. (updated_shift)) then
            ! cks: reset CG after "cg_max" steps
            ! ndmh: or after a fitting failure
            ! lr408: or if conduction shift has just been updated
            cg_coeff = 0.0_DP
            cg_count = 0
            if (pub_on_root  .and. (.not.line_search_success) .and. &
                 (pub_output_detail >= NORMAL) ) write(stdout,'(a)') &
                 'Resetting NGWF CG'
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
                       write(stdout,*) ' previous_g_dot_g=', &
                       previous_g_dot_g, &
                       'CG coeffient set to zero in ngwf_cg_optimise'
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

      ! cks: HACK
      !      cg_coeff = 0.0_DP      ! STEEPEST DESCENTS ONLY
      !      write(stdout,*)'Steepest descents enforced!'

      ! Find search direction
      if (pub_elec_cg_max > 0) then
         direction_on_grid = -cov_grad_on_grid + cg_coeff * &
              prev_direction_on_grid
      else
         ! cks:steepest descents
         direction_on_grid = -cov_grad_on_grid
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

      use rundat, only: pub_ngwf_max_grad

      implicit none

      write(stdout,'(a,f22.14)')   'RMS gradient                = ',rms_gradient
      if (pub_ngwf_max_grad>0.0_DP) write(stdout,'(a,f22.14)') &
           'MAX gradient                = ',max_gradient
      if (index(pub_devel_code,'NGWF_FD')>0) then
         do ifd=1,nfd
            write(stdout,'(a,f22.6)')    'FD trial step length        = ', &
                 fd_trial_step(ifd)
         end do
      end if
      write(stdout,'(a,f22.6)')    'Trial step length           = ',trial_length
      write(stdout,'(a,f22.11)')   'Gradient along search dir.  = ',G_init
      if (index(pub_devel_code,'NGWF_FD')>0) then
         do ifd=1,nfd
            write(stdout,'(a,f22.11)')   'Gradient by FD              = ',&
                 (FFD(ifd)-F0)/fd_trial_step(ifd)
         end do
      end if
      if (check_conjugacy) then
         write(stdout,'(a,f22.11)')   'Conjugacy Test              = ', &
              previous_dir_dot_g
      end if
      write(stdout,'(a,f22.14)')   'Functional at step 0        = ',F0
      if (index(pub_devel_code,'NGWF_FD')>0) then
         do ifd=1,nfd
            write(stdout,'(a,f22.14)')   'Functional at FD step       = ',FFD(ifd)
         end do
      end if
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
      if ( pub_elec_cg_max > 0 ) then
         write(stdout,'(a,f22.6)')    'Conjugate gradients coeff.  = ',cg_coeff
      endif
      write(stdout,'(a)')'--------------------------- NGWF line search &
           &finished --------------------------'
      write(stdout,'(a)')'                                               '

    end subroutine internal_print_search_info


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine ngwf_cg_optimise

end module ngwf_cg
