! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================== !
!                                                                   !
! Constrain-potential Conjugate Gradients optimisation module       !
!                                                                   !
! This module performs the optimisation (i.e.**maximisation**) of   !
! the electronic energy with respect to the constraining-potentials !
! for the constrained-DFT functionality of ONETEP.                  ! 
!-------------------------------------------------------------------!
! Written by Gilberto Teobaldi in December 2011 using the           !
! ngwf_cg_mod.F90 module as inspiration...                          !
!===================================================================!

module cdft_cg

  implicit none

  private

  public :: cdft_u_cg_optimise
  public :: cdft_restart_read
  public :: cdft_update_matrices

!gibo (4 humans)
! Be aware that this module contains **3** public subroutines
! They are listed above.
! In turn, the cdft_u_cg_optimise subroutine (below) contains
! several internal (private) subroutines/functions, 
! whose name start with 'initial_'
CONTAINS

  subroutine cdft_u_cg_optimise(total_energy, ngwfs_converged, cdft_converged, & !   out
       ham, denskern, rep, ngwf_nonsc_forces, lhxc_fine, &              ! inout
       ngwf_basis, proj_basis, hub_proj_basis, hub, &                   ! in
       elements, ewald_energy, localpseudo_fine, core_density_fine , &  ! in
       val_rep, val_ngwf_basis, val_dkn, val_ham) ! lr408: optional cond args

    !==========================================================================!
    ! This subroutine maximise the total energy with respect to the            !
    ! constraining-potential of the constrained-DFT mode                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! total_energy    (output): the total energy                               !
    ! ngwf_converged       (output): whether the NGWFs were converged          !
    ! cdft_converged       (output): whether the cDFT U-pot were converged     !
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
    ! For Conduction NGWF optimisation only:                                   !
    ! val_rep          (input): Valence NGWF representation (optional)         !
    ! val_ngwf_basis   (input): Valence NGWF basis (optional)                  !
    ! val_dkn          (input): Valence density kernel (optional)              !
    !--------------------------------------------------------------------------!
    ! Adapted from ngwf_cg_optimise by Gilberto Teobaldi in December 2011      !
    !==========================================================================!


    use cell_grid, only: pub_fine_grid
    use comms, only: comms_abort, pub_on_root, pub_root_node_id, pub_my_node_id, &
                     comms_bcast, comms_reduce, &
                     comms_barrier
    use constants, only: DP, verbose, normal, stdout, brief, &
        HARTREE_IN_EVS
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_indep_matrices
    use hubbard_build, only: HUBBARD_MODEL, cdft_energy_info
    use hubbard_init, only: h_species
    use ion, only: element
    use kernel, only: kernel_init_core_ham
    use ngwf_cg, only: ngwf_cg_optimise
    use ngwf_diis, only: ngwf_diis_optimise
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use palser_mano, only: palser_mano_kernel_optimise
    use parallel_strategy, only: pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node
    use restart, only: restart_kernel_write
    use rundat, only: pub_output_detail, pub_maxit_palser_mano, pub_aug, &
        pub_any_nl_proj, coreham_denskern_guess, pub_rootname, pub_devel_code,&
        print_qc, maxit_ngwf_cg, maxit_ngwf_diis, pub_nonsc_forces, &
        pub_write_converged_dk_ngwfs, write_denskern, &
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge, pub_cdft_group_spin, pub_cdft_group_charge_diff,&
        pub_cdft_group_spin_diff, pub_cdft_group_u, pub_cdft_group_diff_u, &
        maxit_cdft_u_cg, cdft_cg_type, &
        cdft_cg_threshold, pub_cdft_cg_max, cdft_cg_max_step, pub_cdft_continuation
    use services, only: services_flush, services_line_search_parabola, &
         services_cubic_fit_maximum
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_copy, sparse_axpy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
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
    real(kind=DP), intent(inout) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(inout) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(out)   :: total_energy
    type(NGWF_HAM), intent(inout) :: ham
    ! ars: non self-consistent forces correction due to NGWF loop
    real(kind=DP), intent(inout) :: ngwf_nonsc_forces(:,:)
    !logical, intent(out) :: converged  ! pdh: convergence flag
    logical, intent(out) :: ngwfs_converged  ! pdh: convergence flag
    logical, intent(out) :: cdft_converged   !gibo: cDFT convergence flag

    ! lr408: Optional arguments needed for conduction calculation
    type(NGWF_REP), optional, intent(in) :: val_rep
    type(NGWF_HAM), optional, intent(in) :: val_ham
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis
    type(SPAM3), optional, intent(in)      :: val_dkn(pub_cell%num_spins)


    ! Local Variables
    !local-array for cdft_gradient
    real(kind=DP), allocatable, dimension(:) :: local_cdft_gradient
    !local-array for previous step cdft_gradient
    real(kind=DP), allocatable, dimension(:) :: local_previous_cdft_gradient
    !local-array for U-opt direction
    real(kind=DP), allocatable, dimension(:) :: opt_direction
    !local-array for PREVIOUS U-opt direction
    real(kind=DP), allocatable, dimension(:) :: prev_opt_direction
    !local-array for (being-)optimised U-potentials
    real(kind=DP), allocatable, dimension(:) :: local_start_cdft_u
    real(kind=DP), allocatable, dimension(:) :: delta_u
    real(kind=DP) :: total_energy_current
    !real(kind=DP) :: hubbard_energy
    real(kind=DP) :: last_n_energies(3) ! space to store the 3 most recent energies
    real(kind=DP) :: rms_gradient
    real(kind=DP) :: max_gradient
    real(kind=DP) :: previous_rms_gradient ! RMS NGWF grad of previous iteration
    real(kind=DP) :: line_search_coeff
    !real(kind=DP) :: lnv_threshold
    !real(kind=DP) :: ngwf_threshold
    real(kind=DP) :: F0,F1,F2,G_init,trial_length
    real(kind=DP) :: trial_length_x2
    real(kind=DP) :: cg_coeff
    real(kind=DP) :: quadratic_coeff,cubic_coeff,rejected_quadratic_coeff
    real(kind=DP) :: predicted_functional
    real(kind=DP) :: previous_g_dot_g
    real(kind=DP) :: current_g_dot_g
    real(kind=DP) :: previous_dir_dot_g
    character(len=80), allocatable, dimension(:) :: summary_lines
    integer :: iteration         ! current iteration
    integer :: cg_count          ! current number of steps since CG reset
    integer :: is         ! pdh: spin loop counter
    integer :: ierr       ! error flag
    integer :: minit      ! qoh: Minimum number of CG iterations
    logical :: trial2     ! ndmh: flag to perform second trial step
    logical :: retrial1   ! ndmh: flag to perform repeat first trial step
    logical :: reversing  ! ndmh: line search is going uphill
    logical :: line_search_success ! ndmh: line search fit success flag
    logical :: check_conjugacy     ! ndmh: flag for doing cg conjugacy check
    !logical :: updated_shift ! lr408: Flag needed for conduction calculations

!4 GIBO this local variables are needed
    integer :: hat, species
    !integer :: node
    integer :: g_counter ! 1D U-gradient array counter
    integer :: cdft_gradient_size
    !real(kind=DP) :: rms_cdft_grad_up, rms_cdft_grad_down, rms_cdft_grad_tot
    real(kind=DP) :: cdft_threshold

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering cdft_u_cg_optimise'
#endif

    ! Flush output
    call services_flush

    ! Start timer
    call timer_clock('cdft_u_cg_optimise',1)


    ! allocate workspace for U-optimisation
    cdft_gradient_size = pub_cell%num_spins*pub_cell%nat_hub
    allocate(local_cdft_gradient(cdft_gradient_size),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise', &
         'local_cdft_gradient',ierr)
    allocate(opt_direction(cdft_gradient_size),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise', &
         'opt_direction',ierr)
    allocate(local_start_cdft_u(cdft_gradient_size),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise', &
         'local_start_cdft_u',ierr)
    allocate(delta_u(cdft_gradient_size),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise', &
         'delta_u',ierr)

    if (pub_cdft_cg_max > 0) then
     allocate(prev_opt_direction(cdft_gradient_size),stat=ierr) 
     call utils_alloc_check('cdft_u_cg_optimise', &
          'prev_opt_direction',ierr)
    endif

     if (cdft_cg_type == 'NGWF_POLAK') then
      allocate(local_previous_cdft_gradient(cdft_gradient_size), stat=ierr)
      call utils_alloc_check('cdft_u_cg_optimise', &
           'local_previous_cdft_gradient',ierr)
     endif 

    ! cks: <<< parameter initialisations >>>
    cdft_threshold     = cdft_cg_threshold

    ! cks: <<< variable intialisations >>>
    !current_maxit_lnv     = minit_lnv
    previous_rms_gradient = huge(1.0_DP)
    max_gradient          = huge(1.0_DP)
    line_search_coeff  = 0.15_DP
    trial_length       = 0.1_DP
    F0=0.0_DP ; F1=0.0_DP ; F2=0.0_DP ; G_init=0.0_DP ; cg_coeff=0.0_DP
    !rms_gradient=1.0_DP ; mu=0.0_DP ; total_energy=0.0_DP ; cg_count=0
    rms_gradient=1.0_DP ;  cg_count=0
    total_energy=0.0_DP ; total_energy_current=0.0_DP 
    predicted_functional =0.0_DP
    current_g_dot_g  =0.0_DP
    previous_g_dot_g =0.0_DP
    line_search_success = .true.
    trial2 = .false.
    retrial1 = .false.
    reversing = .false.
    prev_opt_direction(:) = 0.0_DP
    last_n_energies(:) = -huge(1.0_DP)
    quadratic_coeff = 0.0_DP
    rejected_quadratic_coeff = 0.0_DP
    cubic_coeff = 0.0_DP
    cdft_converged = .false.


    !allocate(summary_lines(max(maxit_ngwf_cg+1,2)),stat=ierr)
    allocate(summary_lines(max(maxit_cdft_u_cg+1,2)),stat=ierr)
    call utils_alloc_check('cdft_u_cg_optimise','summary_lines',ierr)


    ! if this a cDFT restart, get the latest U-potentials from (.cdft) file
    !if (pub_cdft_continuation) call cdft_restart_read()
    if (pub_cdft_continuation) call cdft_restart_read(hub)


    !gibo: prepare to output summary of cDFT-U optimisation
    if (pub_on_root) write(summary_lines(1),'(a80)') '|ITER|    RMS GRADIENT   &
         &|     TOTAL ENERGY    |   step   |     Epredicted  '
    summary_lines(1) = adjustl(summary_lines(1))

    !gibo: First message for brief output level
    if (pub_on_root .and. pub_output_detail == BRIEF) then
       write(stdout,'(/a)')'oooooooooooooooooooooooooooooooooooooooo&
            &oooooooooooooooooooooooooooooooooooooooo'
       write(stdout,'(a)')'ooooooooooooooooooooo &
            & cDFT U-optimisation  oooooooooooooooooooooooooooooooooooo'
       write(stdout,'(a)')'oooooooooooooooooooooooooooooooooooooooo&
            &oooooooooooooooooooooooooooooooooooooooo'
       write(stdout,'(a80)') summary_lines(1)
    end if

    !Allow blank calculations for timings purposes if maxit_cdft_u_cg < 0
    minit = 1
    if (maxit_cdft_u_cg < 0) minit = 0

    iteration = 0 ! ndmh: prevent unitialised variable when maxit_ngwf_cg = 0

    !gibo: U-optimisation loop: START
    OPT_LOOP: do iteration=1,max(maxit_cdft_u_cg,minit)


       ! cks: ++++++++++ START ITERATION HEADER +++++++++++++++++++++++++++++
       if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
            maxit_cdft_u_cg > 0) then
          write(stdout,'(/a)')'oooooooooooooooooooooooooooooooooooooooo&
               &oooooooooooooooooooooooooooooooooooooooo'
          write(stdout,'(a,i4.3,a)')'ooooooooooooooooooooooooooo cDFT &
                  &CG iteration ',iteration,' ooooooooooooooooooooooooooooo'
          write(stdout,'(a/)')'oooooooooooooooooooooooooooooooooooooooo&
               &oooooooooooooooooooooooooooooooooooooooo'
       end if
       call services_flush
       ! cks: ++++++ END ITERATION HEADER +++++++++++++++++++++++++++++



       !gibo: write list of the latest cDFT U-potentials into file
       call internal_cdft_restart_write()

       !gibo: ============ OPTIMISE NGWF (and DENSITY KERNEL) ===========
       if ((maxit_ngwf_cg .gt. 0).or.(maxit_ngwf_diis==0).or.pub_nonsc_forces) then
          call ngwf_cg_optimise(total_energy, ngwfs_converged, &          ! out
               ham, denskern, rep, ngwf_nonsc_forces, lhxc_fine, &        ! inout
               ngwf_basis, proj_basis, hub_proj_basis, hub, &             ! in
               elements, ewald_energy, localpseudo_fine, core_density_fine )
       elseif (maxit_ngwf_diis .gt. 0) then
          call ngwf_diis_optimise(total_energy, ngwfs_converged, ham, &
               denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, hub, &
               elements, ewald_energy, localpseudo_fine, core_density_fine, &
               lhxc_fine)
       endif

       total_energy_current = total_energy
       !gibo: ============ OPTIMISE NGWF (and DENSITY KERNEL) ===========


        !gibo: reduce hub%cdft_gradient(:,:) across the nodes
        call comms_reduce('SUM',hub%cdft_gradient)

      !======= CDFT (U) GRADIENT=============== START
      !gibo: store locally the hub%cdft_gradient
      ![previously calculated in sbrtne cdft_energy_total,
      ! see hubbard_build_mod.F90]
      g_counter = 0
      do hat = 1, pub_cell%nat_hub
       do is = 1, pub_cell%num_spins
        g_counter = g_counter + 1
        local_cdft_gradient(g_counter) = hub%cdft_gradient(is, hat)
       enddo
      enddo 
      !======= CDFT (U) GRADIENT=============== END

       ! ndmh: check conjugacy condition on search direction
       check_conjugacy = .false.
       if (check_conjugacy) then
        !gibo(4 humans): ddot is the double precision BLAS
        !                subroutine dedicated to the dot_product between two vectors...
        previous_dir_dot_g = wrappers_ddot(cdft_gradient_size,&
             local_cdft_gradient,1,prev_opt_direction,1)
!          previous_dir_dot_g = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
!               contra_grad_on_grid,1,prev_direction_on_grid,1)
!          call comms_reduce('SUM',previous_dir_dot_g)
       end if

      !======= TEST CDFT-OPT CONVERGENCE=============== START
      cdft_converged = internal_test_convergence()
      call comms_bcast(pub_root_node_id,cdft_converged)
      if (cdft_converged) then
        total_energy = total_energy_current
        EXIT OPT_LOOP ! exit U-optimisation loop
      endif
      !======= TEST CDFT-OPT CONVERGENCE=============== END


       ! pdh: exit if no cDFT U-optimisation is required
       if (maxit_cdft_u_cg == 0) then
          cdft_converged = .false.
          total_energy = total_energy_current
          EXIT OPT_LOOP ! exit U-optimisation loop
       end if

       ! cks: **************** FUNCTIONAL AT INITIAL POINT *********************
       !gibo: use total_energy [from last NGWF optimisation] for F0
       F0 = total_energy_current 
       !F0 = electronic_lagrangian(total_energy, rep%overlap, denskern, &
       !     ham%ham, mu, rep%n_occ)
       ! cks: ************** END FUNCTIONAL AT INITIAL POINT *******************

       ! cks: ************************ START LINE SEARCH ********************

       ! find new CG-(maximisation) direction
       call internal_find_direction

       ! cks: store direction if doing conjugate gradients
       !if (pub_elec_cg_max > 0) then
       if (pub_cdft_cg_max > 0) then
        prev_opt_direction(:) = opt_direction(:) 
          !prev_direction_on_grid = direction_on_grid
       endif

       !if (ngwf_cg_type == 'NGWF_POLAK') &
       if (cdft_cg_type == 'NGWF_POLAK') &
            local_previous_cdft_gradient(:) = local_cdft_gradient(:)
            !prev_contra_grad_on_grid = contra_grad_on_grid

       previous_g_dot_g = current_g_dot_g

      !gibo: make a local copy of the hub%cdft_u... potentials
      !      and change hub%cdft_u along the search direction
      if (pub_cdft_atom_charge .OR. pub_cdft_group_charge .OR. &
          pub_cdft_group_charge_diff) then

        g_counter = 0
        do hat = 1, pub_cell%nat_hub
          g_counter = g_counter + 2
          local_start_cdft_u(g_counter-1) = hub%cdft_u_charge_up(hat)
          local_start_cdft_u(g_counter)   = hub%cdft_u_charge_down(hat)
        enddo

       elseif (pub_cdft_atom_spin .OR. pub_cdft_group_spin .OR. &
               pub_cdft_group_spin_diff) then

        g_counter = 0
        do hat = 1, pub_cell%nat_hub
          g_counter = g_counter + 2
          local_start_cdft_u(g_counter-1) = hub%cdft_u_spin(hat)
          local_start_cdft_u(g_counter)   = hub%cdft_u_spin(hat)
        enddo

      endif

!!      !update cDFT U-potentials on the basis of local_start_cdft_u and delta_u
      call internal_cdft_update_u(trial_length)

      !gibo: write a list of the latest cDFT U-potentials into file
      call internal_cdft_restart_write()

!      ! update hub%up/down_matrices, which contain the
!      ! (diagonal) cDFT Hamiltonian
      !call internal_cdft_update_matrices(hub, hub_proj_basis)
      call cdft_update_matrices(hub, hub_proj_basis)

      !gibo: re-initialise the hamiltonian and density kernel,
      !      but keep the latest localpseudo_fine, core_density_fine,
      !      lhxc_fine and rep
      call internal_cdft_kernel_reinit()

       
       ! cks: %%%%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 1 %%%%%%%%%%%%%%%%%%

       !gibo: ============ OPTIMISE NGWF (and DENSITY) KERNEL ===========
       if ((maxit_ngwf_cg .gt. 0).or.(maxit_ngwf_diis==0).or.pub_nonsc_forces) then
          call ngwf_cg_optimise(total_energy, ngwfs_converged, &          ! out
               ham, denskern, rep, ngwf_nonsc_forces, lhxc_fine, &        ! inout
               ngwf_basis, proj_basis, hub_proj_basis, hub, &             ! in
               elements, ewald_energy, localpseudo_fine, core_density_fine )
       elseif (maxit_ngwf_diis .gt. 0) then
          call ngwf_diis_optimise(total_energy, ngwfs_converged, ham, &
               denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, hub, &
               elements, ewald_energy, localpseudo_fine, core_density_fine, &
               lhxc_fine)
       endif

       F1 = total_energy
       !gibo: ============ OPTIMISE NGWF (and DENSITY) KERNEL ===========

       ! cks: %%%%%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 1 %%%%%%%%%%%%%%%%


       ! cks: ===================== SEARCH BY FITTING PARABOLA =================
       call services_line_search_parabola(&
            quadratic_coeff, predicted_functional,line_search_success, & !output
            G_init, F0, F1, trial_length, cdft_cg_max_step)              !input
       !call services_line_search_parabola(&
       !     quadratic_coeff, predicted_functional,line_search_success, & !output
       !     G_init, F0, F1, trial_length, ngwf_cg_max_step)              !input

       line_search_coeff = quadratic_coeff
       ! cks: ================= END SEARCH BY FITTING PARABOLA =================

       ! cks: CUBIC CUBIC CUBIC ----- SEARCH BY FITTING CUBIC ---- CUBIC CUBIC
       !gibo: inverted because we are (hopefully) dealing with a convex parabola
       !      and the gradient is positive (G_init >0 )along a good search direction
       !      (quadratic_coeff >0). Therefore we need trial2=.true. *ONLY*
       !      if G_init * quadratic_coeff < 0. 
       !CUBIC: if (G_init * quadratic_coeff > 0.0_DP) then
       CUBIC: if (G_init * quadratic_coeff < 0.0_DP) then

          trial2 = .true.

!!      !update cDFT U-potentials on the basis of local_start_cdft_u and delta_u
        trial_length_x2 = 2.0_DP*trial_length
        call internal_cdft_update_u(trial_length_x2)
        trial_length = 0.5_DP*trial_length_x2

        ! update hub%up/down_matrices, which contain the
        ! (diagonal) cDFT Hamiltonian
        !call internal_cdft_update_matrices(hub, hub_proj_basis)
        call cdft_update_matrices(hub, hub_proj_basis)

!          rep%ngwfs_on_grid = start_ngwfs_on_grid + &
!               2.0_DP * trial_length * direction_on_grid
       end if CUBIC


       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP
       ! ndmh: protection against bad line search results: redo trial step if
       ! ndmh: line search result is too much bigger or smaller than trial step
       ! gibo: or if we are searching 'downhill' and went past trial step position
       if ( ((line_search_coeff < 0.05_DP*trial_length) .or. &
            (line_search_coeff > 20.0_DP*trial_length) .or. &
            (.not. line_search_success) .or. &
            (reversing .and. (line_search_coeff > trial_length))) &
            .and. .not. trial2 ) then

          retrial1 = .true.
       end if

       RETRIAL_1: if (retrial1) then

        rejected_quadratic_coeff = quadratic_coeff

!!!      !update cDFT U-potentials on the basis of local_start_cdft_u and delta_u
        call internal_cdft_update_u(quadratic_coeff)

        ! update hub%up/down_matrices, which contain the
        ! (diagonal) cDFT Hamiltonian
        !call internal_cdft_update_matrices(hub, hub_proj_basis)
        call cdft_update_matrices(hub, hub_proj_basis)

!          rep%ngwfs_on_grid = start_ngwfs_on_grid + &
!               quadratic_coeff * direction_on_grid
       end if RETRIAL_1
       ! ndmh: SECOND TRIAL STEP ----- SECOND TRIAL STEP ----- SECOND TRIAL STEP

       if (trial2 .or. retrial1) then

          ! cks: %%%%%%%%%%%%%%%% FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%
       !gibo: write list of the latest cDFT U-potentials into file
       call internal_cdft_restart_write()

      !gibo: re-initialise the hamiltonian and density kernel,
      !      but keep the latest localpseudo_fine, core_density_fine,
      !      lhxc_fine and rep
      call internal_cdft_kernel_reinit()

          !gibo: ============ OPTIMISE NGWF (and DENSITY) KERNEL ===========
          if ((maxit_ngwf_cg .gt. 0).or.(maxit_ngwf_diis==0).or.pub_nonsc_forces) then
             call ngwf_cg_optimise(total_energy, ngwfs_converged, &          ! out
                  ham, denskern, rep, ngwf_nonsc_forces, lhxc_fine, &        ! inout
                  ngwf_basis, proj_basis, hub_proj_basis, hub, &             ! in
                  elements, ewald_energy, localpseudo_fine, core_density_fine )
          elseif (maxit_ngwf_diis .gt. 0) then
             call ngwf_diis_optimise(total_energy, ngwfs_converged, ham, &
                  denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, hub, &
                  elements, ewald_energy, localpseudo_fine, core_density_fine, &
                  lhxc_fine)
          endif
   
          F2 = total_energy
          !gibo: ============ OPTIMISE NGWF (and DENSITY) KERNEL ===========

          ! cks: %%%%%%%%%%%%% END FUNCTIONAL AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%


          if (trial2) then
             ! Cubic Fit
             !call services_cubic_fit_minimum( &
             !     cubic_coeff, predicted_functional, line_search_success, & !output
             !     F0, F1, F2, G_init, trial_length, 2.0_DP*trial_length, &  !input
             !     cdft_cg_max_step)                                         !input

             call services_cubic_fit_MAXIMUM( &
                  cubic_coeff, predicted_functional, line_search_success, & !output
                  F0, F1, F2, G_init, trial_length, 2.0_DP*trial_length, &  !input
                  cdft_cg_max_step)                                         !input

             line_search_coeff = cubic_coeff
          end if

          if (retrial1) then
             ! ndmh: quadratic fit at new trial length
             call services_line_search_parabola(&
                  quadratic_coeff, predicted_functional, &               !output
                  line_search_success, G_init, F0, F2, &                 !input
                  line_search_coeff, cdft_cg_max_step)                   !input

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

!       ! cks: set the new values of the NGWFs
!       rep%ngwfs_on_grid = start_ngwfs_on_grid + &
!            line_search_coeff * direction_on_grid

!!      !update cDFT U-potentials on the basis of local_start_cdft_u and delta_u
      call internal_cdft_update_u(line_search_coeff)

!      ! update hub%up/down_matrices, which contain the
!      ! (diagonal) cDFT Hamiltonian
      !call internal_cdft_update_matrices(hub, hub_proj_basis)
      call cdft_update_matrices(hub, hub_proj_basis)

      !gibo: re-initialise the hamiltonian and density kernel,
      !      but keep the latest localpseudo_fine, core_density_fine,
      !      lhxc_fine and rep
      call internal_cdft_kernel_reinit()

        !if required, print out information on current line search
        if (pub_on_root .and. pub_output_detail >= NORMAL) &
            call internal_print_search_info

       ! cks: ************************* END LINE SEARCH ***********************

       ! ndmh: reset flags
       trial2 = .false.
       retrial1 = .false.
       reversing = .false.
       F2 = 0.0_DP

       ! cks: write line with summary info for curent iteration
       if (pub_on_root) then
          !if(abs(total_energy)<100000_DP) then
          if(abs(total_energy_current)<100000_DP) then
             write(summary_lines(iteration+1),'(i4,f21.14,f22.14,f11.6,f22.14)')&
                iteration,rms_gradient,total_energy_current,line_search_coeff, &
                predicted_functional
          else
             write(summary_lines(iteration+1),'(i4,f21.14,f22.12,f11.6,f22.12)')&
                iteration,rms_gradient,total_energy_current,line_search_coeff, &
                predicted_functional
          end if
          if (pub_output_detail == BRIEF) then
             write(stdout,'(a80)') summary_lines(iteration+1)
          end if
       end if

       !gibo: write intermediate summary of U-optimisation
       call internal_cdft_opt_summary


    !gibo: U-optimisation loop: END
    enddo OPT_LOOP


    if (.not. cdft_converged) then
       ! cks: reset iteration number just for storing final line of calculation
       !      summary
       iteration = iteration - 1
       ! cks: print warning that calculation failed to converge
       if (pub_on_root) write(stdout,'(a,i4,a)') &
            'WARNING: maximum number of cDFT CG iterations (',maxit_cdft_u_cg, &
            ') exceeded!'
    end if


    ! cks: print calculation summary
    call internal_calculation_summary

    ! cks: print quality control information
    if (print_qc) call internal_qc_output

    ! gibo: write out cDFT(+U) occupancies if it hasn't already been done
    if (pub_output_detail .ne. VERBOSE) then
       !call hubbard_energy_info(hub,hub_proj_basis)
       call cdft_energy_info(hub,hub_proj_basis)
    endif

    !gibo: write out atom-resolved  summary of U-optimisation
    call internal_cdft_opt_summary

    ! deallocate workspace for U-optimisation
    deallocate(local_cdft_gradient,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','local_cdft_gradient',ierr)
    deallocate(opt_direction,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','opt_direction',ierr)
    deallocate(local_start_cdft_u,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','local_start_cdft_u',ierr)
    deallocate(delta_u,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','delta_u',ierr)

    if (pub_cdft_cg_max > 0) then
     deallocate(prev_opt_direction,stat=ierr)
     call utils_dealloc_check('cdft_u_cg_optimise', &
          'previous_opt_direction',ierr)
    endif

    if (cdft_cg_type == 'NGWF_POLAK') then
     deallocate(local_previous_cdft_gradient,stat=ierr)
     call utils_dealloc_check('cdft_u_cg_optimise', &
          'local_previous_cdft_gradient',ierr)
    endif

    deallocate(summary_lines,stat=ierr)
    call utils_dealloc_check('cdft_u_cg_optimise','summary_lines',ierr)     

    call timer_clock('cdft_u_cg_optimise',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving cdft_u_cg_optimise'
#endif

    ! Flush output
    call services_flush

    return

! gibo (4 humans)
! Below, the (internal) subrtuines/functions contained by the 
! 'cdft_u_cg_optimise' subroutine
  contains

!##############################################################################
!##############################################################################

    subroutine internal_cdft_update_u(step)

    !==========================================================!
    ! This subroutine updates the cDFT U potentials.           !
    ! It also protects agains crazy steps which would turn     !
    ! acceptor-atoms into donor-atoms [U_old*U_new <0]         !
    ! and viceversa...                                         !
    !----------------------------------------------------------!
    ! Written for cDFT-OPT module by Gilberto Teobaldi         !
    ! in December 2011                                         !
    !==========================================================!

    implicit none

    real(kind=DP), intent(inout) :: step

    integer :: hat, g_counter
    real(kind=DP) :: delta_scale
    real(kind=DP), parameter :: tol_u = 1.E-9_DP
    logical :: good_step

    delta_scale = 1.0_DP

    !update hub%cdft_u_charge_up/down .AND. the (energy-)suggested cDFT step
    CDFT_TYPE: if (pub_cdft_atom_charge .OR. pub_cdft_group_charge .OR. &
        pub_cdft_group_charge_diff) then


      !protect against sign-inversion of U-potentials
      CHECK_DONOR_ACCEPTOR_charge: do

         good_step =  .TRUE.

         !scale line_search step according to delta_scale
         step = step * delta_scale

         g_counter = 0
         do hat = 1, pub_cell%nat_hub
           g_counter = g_counter + 2

           ! work out tentative change for cDFT potentials
           delta_u(g_counter-1) = step*opt_direction(g_counter-1)
           delta_u(g_counter)   = step*opt_direction(g_counter)

           ! update cDFT potentials
           hub%cdft_u_charge_up(hat) = local_start_cdft_u(g_counter-1) &
                    + delta_scale*delta_u(g_counter-1)
           hub%cdft_u_charge_down(hat) = local_start_cdft_u(g_counter) &
                    + delta_scale*delta_u(g_counter)
 
           !change of sign for potential_UP of atom hat! Reduce delta_scale
           if (hub%cdft_u_charge_up(hat)*local_start_cdft_u(g_counter-1) < &
                epsilon(1.0_DP) ) then
                delta_scale = 0.1_DP*delta_scale
                good_step = .FALSE.
                CYCLE CHECK_DONOR_ACCEPTOR_charge
           endif
 
           !change of sign for potential_DOWN of atom hat! Reduce delta_scale
           if (hub%cdft_u_charge_down(hat)*local_start_cdft_u(g_counter) < &
                epsilon(1.0_DP) ) then
                delta_scale = 0.1_DP*delta_scale
                good_step = .FALSE.
                CYCLE CHECK_DONOR_ACCEPTOR_charge
           endif
 
         enddo
 
         if (good_step) EXIT CHECK_DONOR_ACCEPTOR_charge
 
         if (delta_scale < tol_u) then
           if (pub_on_root) &
                write(stdout,'(a/a,a)') &
                'WARNING in internal_cdft_update_u:',    &
                ' prevented change of sign for cDFT potentials!!!', &
                ' Disaster may loom large...'
 
           EXIT CHECK_DONOR_ACCEPTOR_charge
         endif

      enddo CHECK_DONOR_ACCEPTOR_charge


    !update hub%cdft_u_charge_up/down .AND. the (energy-)suggested cDFT step
    elseif (pub_cdft_atom_spin .OR. pub_cdft_group_spin .OR. &
             pub_cdft_group_spin_diff) then

      !protect against sign-inversion of U-potentials
      CHECK_DONOR_ACCEPTOR_spin: do

         good_step =  .TRUE.

         !scale line_search step according to delta_scale
         step = step * delta_scale

         g_counter = 0
         do hat = 1, pub_cell%nat_hub
           g_counter = g_counter + 2

           ! work-out tentative change for cDFT potentials
           delta_u(g_counter-1) = step*opt_direction(g_counter-1)
           !delta_u(g_counter)   = step*opt_direction(g_counter)

           ! update cDFT potentials
           hub%cdft_u_spin(hat) = local_start_cdft_u(g_counter-1) &
                    + delta_scale*delta_u(g_counter-1)
           ! gibo: no need to update for spin-beta 
           !       [implictly accounted for in cdft_energy_info]
      
           !change of sign for spin-potential of atom hat! Reduce delta_scale
           if (hub%cdft_u_spin(hat)*local_start_cdft_u(g_counter-1) < &
                epsilon(1.0_DP) ) then
                delta_scale = 0.1_DP*delta_scale
                good_step = .FALSE.
                CYCLE CHECK_DONOR_ACCEPTOR_spin
           endif
      
         enddo
      
         if (good_step) EXIT CHECK_DONOR_ACCEPTOR_spin
      
         if (delta_scale < tol_u) then
           if (pub_on_root) &
                write(stdout,'(a/a,a)') &
                'WARNING in internal_cdft_update_u:',    &
                ' prevented change of sign for cDFT potentials!!!', &
                ' Disaster may loom large...'
      
           EXIT CHECK_DONOR_ACCEPTOR_spin
         endif

      enddo CHECK_DONOR_ACCEPTOR_spin

    endif CDFT_TYPE


!NEW-12.12.11
    !for group_charge/spin_[diff] runs, 
    !update global pub_cdft_group_diff_u potential
    if (pub_cdft_group_charge_diff) &
           pub_cdft_group_diff_u = ABS(hub%cdft_u_charge_up(1))

    if (pub_cdft_group_charge) &
           pub_cdft_group_u = ABS(hub%cdft_u_charge_up(1))

    if (pub_cdft_group_spin_diff) &
           pub_cdft_group_diff_u = ABS(hub%cdft_u_spin(1))

    if (pub_cdft_group_spin) &
           pub_cdft_group_u = ABS(hub%cdft_u_spin(1))

     !we do not really need this since we sum_reduce the cdft_gradient 
     !at each cDFT-step. However, I leave it for extra safety...
    if (pub_cdft_group_charge_diff .OR. pub_cdft_group_spin_diff) &
            call comms_bcast(pub_root_node_id,pub_cdft_group_diff_u)

    if (pub_cdft_group_charge .OR. pub_cdft_group_spin) &
            call comms_bcast(pub_root_node_id,pub_cdft_group_u)
!NEW-12.12.11

    end subroutine internal_cdft_update_u

!!##############################################################################
!!##############################################################################
!
!    subroutine internal_cdft_update_matrices(hub, hub_proj_basis)
!
!    !==========================================================!
!    ! This subroutine updates the cDFT[+U] SPAM3 matrices      !
!    !----------------------------------------------------------!
!    ! Written for cDFT-OPT module,and based on                 !
!    ! hubbard_build_matrices by Gilberto Teobaldi in Dec. 2011 !
!    !==========================================================!
!
!    use comms, only: pub_my_node_id, pub_on_root
!    use function_basis, only: FUNC_BASIS
!    use hubbard_init, only: h_species
!!gibo-start
!    use hubbard_build, only: HUBBARD_MODEL
!!gibo-end
!    use parallel_strategy, only:  pub_num_hub_atoms_on_node, &
!         pub_hub_atoms_on_node, pub_distr_atom, pub_first_atom_on_node
!!gibo-start-10.12.11
!!    use rundat, only: pub_aug, pub_cdft, pub_cdft_hubbard, &
!!         task, pub_hubbard_restart, pub_hubbard_atomsolve
!!gibo-end-10.12.11
!    use simulation_cell, only: pub_cell
!!gibo-start
!    !use sparse, only: SPAM3, sparse_create, &
!    !     sparse_scale, sparse_put_element
!    use sparse, only: SPAM3, sparse_put_element
!!gibo-end
!
!    implicit none
!
!    ! Arguments
!    type(HUBBARD_MODEL), intent(inout) :: hub
!    type(FUNC_BASIS), intent(in) :: hub_proj_basis
!
!    ! Local Variables
!    integer :: sp, hat_on_node, hat, theatom
!    integer :: hub_proj
!    integer :: is, ierr
!    real(kind=DP) :: u_charge_up, u_charge_down, u_spin_up, u_spin_down
!
!#ifdef DEBUG
!    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
!         &hubbard_build_matrices'
!#endif
!
!!====select CDFT-mode and iterate...
!    ! gibo: ATOM/GROUP-CHARGE-[DIFF] constrained run
!    CDFT_MODE: if ((pub_cdft_atom_charge .OR. pub_cdft_group_charge).OR.&
!                   pub_cdft_group_charge_diff) then
!!====select CDFT-mode and iterate...
!
!    ! ddor: Loop over Hubbard atoms on my node
!    cdft_atoms_1 : do hat_on_node = 1, &
!         pub_num_hub_atoms_on_node(pub_my_node_id)
!
!       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
!       sp = hub%species_number(hat)
!       theatom = pub_distr_atom(hub%orig(hat))
!
!       u_charge_up   = hub%cdft_u_charge_up(hat)
!       u_charge_down = hub%cdft_u_charge_down(hat)
!
!       ! ddor: Loop over Hubbard projectors on my node
!       hubprojs_cdft_1: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
!            hub_proj_basis%first_on_atom(theatom) + &
!            hub_proj_basis%num_on_atom(theatom) - 1
!
!          call sparse_put_element(u_charge_up,&
!               hub%up_matrix,hub_proj,hub_proj)
!          call sparse_put_element(u_charge_down,&
!               hub%down_matrix,hub_proj,hub_proj)
!
!       end do hubprojs_cdft_1
!
!    end do cdft_atoms_1
!
!
!!====select CDFT-mode and iterate...
!    ! gibo: ATOM/GROUP-SPIN-[DIFF] constrained run
!    else if ((pub_cdft_atom_spin .OR. pub_cdft_group_spin).OR.&
!                   pub_cdft_group_spin_diff) then
!!====select CDFT-mode and iterate...
!
!    ! ddor: Loop over Hubbard atoms on my node
!    cdft_atoms_2 : do hat_on_node = 1, &
!         pub_num_hub_atoms_on_node(pub_my_node_id)
!
!       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
!       sp = hub%species_number(hat)
!       theatom = pub_distr_atom(hub%orig(hat))
!
!       ! MIND THE SIGNS: Ecdft_spin = Us*[Nup -Ndown - DN_target]
!       ! thus V_up = Us, and V_down = -Us
!       u_spin_up   =  hub%cdft_u_spin(hat)
!       u_spin_down = -hub%cdft_u_spin(hat)
!
!
!       ! ddor: Loop over Hubbard projectors on my node
!       hubprojs_cdft_2: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
!            hub_proj_basis%first_on_atom(theatom) + &
!            hub_proj_basis%num_on_atom(theatom) - 1
!
!          call sparse_put_element(u_spin_up,&
!               hub%up_matrix,hub_proj,hub_proj)
!          call sparse_put_element(u_spin_down,&
!               hub%down_matrix,hub_proj,hub_proj)
!
!       end do hubprojs_cdft_2
!
!    end do cdft_atoms_2
!
!!====select CDFT-mode and iterate...
!    end if CDFT_MODE
!!====select CDFT-mode and iterate...
!
!
!#ifdef DEBUG
!    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
!         &hubbard_build_matrices'
!#endif
!    end subroutine internal_cdft_update_matrices

!##############################################################################
!##############################################################################

    logical function internal_test_convergence()


!gibo-start
!      use rundat, only: delta_e_conv, maxit_pen, pub_hub_ngwf_spin_thr, &
!           pen_param, pub_write_ngwf_grad_plot, pub_elec_energy_tol, &
!           pub_ngwf_max_grad, pub_write_forces, &
!           pub_cdft_max_grad, pub_cdft_elec_energy_tol
      use rundat, only: delta_e_conv, pub_cdft_max_grad, pub_cdft_elec_energy_tol
!gibo-end
      use sparse, only: sparse_is_dense
      use wrappers, only : wrappers_ddot

      implicit none

!4 GIBO: TO BE CHECKED
      !logical :: energy_increasing
      logical :: energy_change_converged
      !logical :: ngwf_grad_converged
      !logical :: ngwf_max_grad_converged
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
!4 GIBO: TO BE CHECKED

!gibo-5.12.11
      logical :: energy_decreasing
      logical :: cdft_grad_converged
      logical :: cdft_max_grad_converged
!gibo-5.12.11

      ! cks: Update storage of most recent three energies
      last_n_energies(1) = last_n_energies(2)
      last_n_energies(2) = last_n_energies(3)
      !last_n_energies(3) = total_energy
      last_n_energies(3) = total_energy_current

      ! ndmh: Check energy change per atom vs tolerance
      !gibo:adapted to check energy change per CDFT-atom vs tolerance
      energy_change_converged = .false.
      !test_e_conv = (pub_elec_energy_tol > 0.0_DP).and.(iteration>1)
      test_e_conv = (pub_cdft_elec_energy_tol > 0.0_DP).and.(iteration>1)
      if (test_e_conv) then
         diff12 = abs(last_n_energies(1)-last_n_energies(2)) &
              !/ real(pub_cell%nat,kind=DP)
              / real(pub_cell%nat_hub,kind=DP)
         diff23 = abs(last_n_energies(2)-last_n_energies(3)) &
              !/ real(pub_cell%nat,kind=DP)
              / real(pub_cell%nat_hub,kind=DP)
         max_energy_diff = max(diff12,diff23)
         !if (max_energy_diff<pub_elec_energy_tol) then
         if (max_energy_diff<pub_cdft_elec_energy_tol) then
            energy_change_converged = .true.
         end if
      end if

      ! calculate rms gradient
      ! cks: parallel calculation of NGWF rms_gradient
      !gibo(4 humans): ddot is the double precision BLAS
      !                subroutine dedicated to the dot_product between two vectors...
!      rms_gradient = wrappers_ddot(pub_cell%num_spins*pub_cell%nat_hub,&
      rms_gradient = wrappers_ddot(cdft_gradient_size,&
           local_cdft_gradient,1,local_cdft_gradient,1)
!gibo:hub%cdft_gradient has already sum-reduced (cdft_u_cg_optimise) so
!     no need to repeat the operation here
! 4 GIBO, we have already reduced the gradient = WE MOST LIKELY DO NOT NEED THIS
!      call comms_reduce('SUM', rms_gradient)
! 4 GIBO, we have already reduced the gradient = WE MOST LIKELY DO NOT NEED THIS
      !rms_gradient = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
      !     contra_grad_on_grid,1,cov_grad_on_grid,1)
      !call comms_reduce('SUM', rms_gradient)

      ! cks: store for Fletcher-Reeves calculation of CG coeff
      current_g_dot_g = rms_gradient
      rms_gradient = sqrt(abs(rms_gradient) / &
                     !real(pub_cell%num_spins*pub_cell%nat_hub, kind=DP))
                     real(cdft_gradient_size, kind=DP))
      !current_g_dot_g = rms_gradient
      !rms_gradient = sqrt(abs(rms_gradient) / real(est_num_psincs, kind=DP))

      ! ndmh: Calculate max value of <g^a|g_a> over all NGWFs on all nodes
      !gibo: Calculate max value of cdft_gradient
      max_gradient = 0.0_DP
      !test_max_grad = (pub_ngwf_max_grad > 0.0_DP)
      test_max_grad = (pub_cdft_max_grad > 0.0_DP)
      if (test_max_grad) then
        do g_counter = 1, cdft_gradient_size
          max_gradient = max(max_gradient, &
          sqrt(local_cdft_gradient(g_counter)*local_cdft_gradient(g_counter) ) )
        enddo
!gibo:hub%cdft_gradient has already sum-reduced (cdft_u_cg_optimise) so
!     no need to repeat the operation here
          
!         do ipt=1,ngwf_basis%n_ppds*pub_cell%n_pts
!            max_gradient = max(max_gradient, &
!                 sqrt(abs(contra_grad_on_grid(ipt)*cov_grad_on_grid(ipt))))
!         end do
!         call comms_reduce('MAX', max_gradient)
      end if

      ! cks: Test for if allowed to use energy gain as convergence criterion
      !gibo: adapted to CDFT Ecdft **OPTIMISATION**
      !energy_increasing = .false.
      energy_decreasing = .false.
      !energy_increasing = delta_e_conv .and. &
      energy_decreasing = delta_e_conv .and. &
           ( last_n_energies(3) < last_n_energies(2)) .and. &
           ( last_n_energies(2) < last_n_energies(1))

      !cDFT RMS gradient convergence criterion
      test_rms_grad = (cdft_threshold > 0.0_DP)
      if (test_rms_grad) then
         cdft_grad_converged = (rms_gradient < cdft_threshold)
      else
         cdft_grad_converged = .false.
      end if

      !cDFT MAX gradient convergence criterion
      !ngwf_max_grad_converged = .false.
      cdft_max_grad_converged = .false.
      if (test_max_grad) then
         cdft_max_grad_converged = (max_gradient < pub_cdft_max_grad)
      end if

      ! ndmh: Check all relevant criteria for convergence
      all_converged = .true.
!      if ((pub_elec_force_tol > 0.0_DP).and.(iteration<2)) &
!           all_converged = .false.
      if ((pub_cdft_elec_energy_tol > 0.0_DP).and.(iteration<2)) &
           all_converged = .false.
!      if (test_f_conv) &
!           all_converged = all_converged.and.forces_converged
      if (test_e_conv) &
           all_converged = all_converged.and.energy_change_converged
      if (test_rms_grad) &
           !all_converged = all_converged.and.ngwf_grad_converged
           all_converged = all_converged.and.cdft_grad_converged
      if (test_max_grad) &
           !all_converged = all_converged.and.ngwf_max_grad_converged
           all_converged = all_converged.and.cdft_max_grad_converged
      if (delta_e_conv) &
           !all_converged = all_converged.or.energy_increasing
           all_converged = all_converged.or.energy_decreasing

!--------------------------------
      ALLCONVERGED: if (all_converged) then

         ! cks: print details only when output is not set to brief
         if (pub_output_detail >= NORMAL) then

            if (pub_on_root) then
               write(stdout,'(/a)')'           oooooooooooooooooooooooooooo&
                    &oooooooooooooooooooooooooooo'
!               write(stdout,'(/a)')'           ............................&
!                    &............................'
               write(stdout,'(a)')'           | *** cDFT optimisation &
                    &converged ***                  |'
               if (test_rms_grad) write(stdout,'(a,f19.14,a)') &
                    !'           | RMS NGWF gradient = ',rms_gradient,'              |'
                    '           o RMS cDFT gradient = ',rms_gradient,'              o'
               if (test_max_grad) write(stdout,'(a,f19.14,a)') &
                    !'           | MAX NGWF gradient = ',max_gradient,'              |'
                    '           o MAX cDFT gradient = ',max_gradient,'              o'
!               if (test_f_conv) write(stdout,'(a,f17.14,a)') &
!                    '           | Maximum change in forces = ',max_force_diff,' au      |'
               if (test_e_conv) write(stdout,'(a,f17.14,a)') &
                    '           o Maximum change in energy = ',max_energy_diff,' / atom  o'
               write(stdout,'(a)')'           o Criteria satisfied: &
                    &                                 o'
               !if (test_rms_grad.and.ngwf_grad_converged) write(stdout,'(a)') &
               !     '           | -> RMS NGWF gradient lower than set threshold.       |'
               if (test_rms_grad.and.cdft_grad_converged) write(stdout,'(a)') &
                    '           o -> RMS cDFT gradient lower than set threshold.       o'
               !if (test_max_grad.and.ngwf_max_grad_converged) write(stdout,'(a)') &
               !     '           | -> MAX NGWF gradient lower than set threshold.       |'
               if (test_max_grad.and.cdft_max_grad_converged) write(stdout,'(a)') &
                    '           o -> MAX cDFT gradient lower than set threshold.       o'
               if (test_e_conv.and.energy_change_converged) write(stdout,'(a)') &
                    '           o -> Energy change per atom lower than set threshold.  o'
!               if (test_f_conv.and.forces_converged) write(stdout,'(a)') &
!                    '           | -> Maximum force change lower than set threshold.    |'
               !if (energy_increasing) then
               if (energy_decreasing) then
                  write(stdout,'(a)') &
                       '           o -> WARNING!!!!                                       o'
                  write(stdout,'(a)') &
                       '           o -> Energy DECREASE over the last 3 cDFT CG-steps     o'
                  write(stdout,'(a)') &
                       '           o -> Disaster may be looming large....                 o'
!                  write(stdout,'(a)') &
!                       '           | -> Much better to stop here, aborting...             |'
!                  call comms_abort
!                  write(stdout,'(a)') &
!                       '           | -> Maximum degree of convergence for applied level   |'
!                  write(stdout,'(a)') &
!                       '           |    of density kernel truncation has been reached.    | '
               endif
               write(stdout,'(a)') '           ooooooooooooooooooooooooooo&
                    &ooooooooooooooooooooooooooooo'
            end if
         end if

!gibo. NGWF (gradient) plotting calls eliminated
!      [dealt with in ngwf_cg_mod.F90]

       internal_test_convergence = .true.

      else

         ! ndmh: print details only when output is set to verbose
         if (pub_output_detail >= VERBOSE) then

            if (pub_on_root) then
               write(stdout,'(/a)')'           oooooooooooooooooooooooooooo&
                    &oooooooooooooooooooooooooooo'
!               write(stdout,'(/a)')'           ............................&
!                    &............................'
               write(stdout,'(a)')'           o *** cDFT optimisation &
                    &not converged ***              o'
               if (test_rms_grad) write(stdout,'(a,f19.14,a)') &
                    !'           | RMS NGWF gradient = ',rms_gradient,'              |'
                    '           o RMS cDFT gradient = ',rms_gradient,'              o'
               if (test_max_grad) write(stdout,'(a,f19.14,a)') &
                    !'           | MAX NGWF gradient = ',max_gradient,'              |'
                    '           o MAX cDFT gradient = ',max_gradient,'              o'
!               if (test_f_conv) write(stdout,'(a,f17.14,a)') &
!                    '           | Maximum change in forces = ',max_force_diff,' au      |'
               if (test_e_conv) write(stdout,'(a,f17.14,a)') &
                    '           o Maximum change in energy = ',max_energy_diff,' / atom  o'
               write(stdout,'(a)')'           o Criteria not satisfied: &
                    &                             o'
               !if (test_rms_grad.and.ngwf_grad_converged) write(stdout,'(a)') &
               !     '           | -> RMS NGWF gradient lower than set threshold.       |'
               if (test_rms_grad.and.(.not.cdft_grad_converged)) &
                    write(stdout,'(a)') &
                    '           o -> RMS cDFT gradient higher than set threshold.      o'
               !if (test_max_grad.and.ngwf_max_grad_converged) write(stdout,'(a)') &
               !     '           | -> MAX NGWF gradient lower than set threshold.       |'
               if (test_max_grad.and.(.not.cdft_max_grad_converged)) &
                    write(stdout,'(a)') &
                    '           o -> MAX cDFT gradient higher than set threshold.      o'
               if (test_e_conv.and.(.not.energy_change_converged)) &
                    write(stdout,'(a)') &
                    '           o -> Energy change per atom higher than set threshold. o'
!               if (test_f_conv.and.forces_converged) write(stdout,'(a)') &
!                    '           | -> Maximum force change lower than set threshold.    |'
               !if (energy_increasing) then
               if (energy_decreasing) then
                  write(stdout,'(a)') &
                       '           o -> WARNING!!!!                                       o'
                  write(stdout,'(a)') &
                       '           o -> Energy DECREASE over the last 3 cDFT CG-steps     o'
                  write(stdout,'(a)') &
                       '           o -> Disaster may be looming large....                 o'
!                  write(stdout,'(a)') &
!                       '           | -> Much better to stop here, aborting...             |'
!                  call comms_abort
!                  write(stdout,'(a)') &
!                       '           | -> Maximum degree of convergence for applied level   |'
!                  write(stdout,'(a)') &
!                       '           |    of density kernel truncation has been reached.    | '
               endif
               write(stdout,'(a)') '           ooooooooooooooooooooooooooo&
                    &ooooooooooooooooooooooooooooo'
            end if
         end if

         internal_test_convergence = .false.

      endif ALLCONVERGED
!--------------------------------

!gibo. NGWF (gradient) plotting calls eliminated
!      [dealt with in ngwf_cg_mod.F90]


    end function internal_test_convergence

!##############################################################################
!##############################################################################

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
           'ooooooooooooooooooooooooooooooo cDFT line search &
           &ooooooooooooooooooooooooooooooo'


      ! calculate CG coefficient
      if ( iteration > 1 ) then
         !if ((cg_count >= pub_elec_cg_max) .or. (.not.line_search_success) &
         !     .or. (updated_shift)) then
         if ((cg_count >= pub_cdft_cg_max) .or. (.not.line_search_success)) then
            ! cks: reset CG after "cg_max" steps
            ! ndmh: or after a fitting failure
            ! lr408: or if conduction shift has just been updated
            cg_coeff = 0.0_DP
            cg_count = 0
            if (pub_on_root  .and. (.not.line_search_success) .and. &
                 (pub_output_detail >= NORMAL) ) write(stdout,'(a)') &
                 !'Resetting NGWF CG'
                 'Resetting cDFT CG'
         else
            ! cks <<< FLETCHER >>>
            !if (ngwf_cg_type == 'NGWF_FLETCHER') then
            if (cdft_cg_type == 'NGWF_FLETCHER') then

               ! cks: original Fletcher-Reeves formula (cheaper in memory)
               if (abs(previous_g_dot_g) > epsilon(1.0_DP)) then
                  cg_coeff = current_g_dot_g / previous_g_dot_g
                  ! cks: protection from crazy coefficients
                  if (abs(cg_coeff) > 2.0_DP ) then
                     if (pub_on_root  .and.(pub_output_detail >= NORMAL) ) &
                          write(stdout,'(a,f8.4,a)') &
                          !'WARNING: NGWF Fletcher-Reeves CG coeff too large (',&
                          'WARNING: cDFT Fletcher-Reeves CG coeff too large (',&
                          cg_coeff, ') - setting to zero'
                     cg_coeff = 0.0_DP
                     cg_count = 0
                  end if
               else
                  if (pub_on_root .and.(pub_output_detail >= NORMAL) ) &
                       write(stdout,*) ' previous_g_dot_g=', &
                       previous_g_dot_g, &
                       !'CG coeffient set to zero in ngwf_cg_optimise'
                       'CG coeffient set to zero in cdft_cg_u_optimise'
                  cg_coeff = 0.0_DP
                  cg_count = 0
               end if

               ! cks: <<< POLAK >>>
            !else if (ngwf_cg_type == 'NGWF_POLAK') then
            else if (cdft_cg_type == 'NGWF_POLAK') then

               ! POLAK FORMULA
               !cg_coeff = services_polak_cg_coeff(prev_direction_on_grid, &
               !     cov_grad_on_grid,contra_grad_on_grid, &
               !     prev_contra_grad_on_grid, ngwf_basis%n_ppds*pub_cell%n_pts)

               !gibo: calculate the POLAK-RIBIERE CG coefficient
               cg_coeff = internal_cdft_polak_cg_coeff()

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
      !if (pub_elec_cg_max > 0) then
      if (pub_cdft_cg_max > 0) then
         !****WARNING****
         ! Since we want to maximise E_cDFT we need to go uphill
         ! [i.e. follow to gradient i.e. no minus the gradient
         !opt_direction(:) = -local_cdft_gradient(:) + &  ! OK for minimisation
         opt_direction(:) =  local_cdft_gradient(:) + &  ! OK for MAXIMISATION
                                 cg_coeff * prev_opt_direction(:)

         !direction_on_grid = -cov_grad_on_grid + cg_coeff * &
         !     prev_direction_on_grid
      else
         ! cks:steepest descents
         !direction_on_grid = -cov_grad_on_grid
         !opt_direction(:) = -local_cdft_gradient(:) ! OK for minimisation
         opt_direction(:) =  local_cdft_gradient(:) ! OK for MAXIMISATION
      endif


      ! Slope of energy in search direction
      G_init = wrappers_ddot(cdft_gradient_size, &
              local_cdft_gradient,1,opt_direction,1)
      !G_init = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
      !     contra_grad_on_grid,1,direction_on_grid,1)

!gibo: no need to collect anything since we are working on 
!      the sum-reduced cdft_gradient
      ! cks: collect the work of each node
      !call comms_reduce('SUM', G_init)

      ! take action in case of positive slope along search direction
      !gibo: for cDFT, take action against NEGATIVE slope along search direction
      !if (G_init > 0.0_DP) then
      if (G_init < 0.0_DP) then

         if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
            write(stdout,'(a,e16.6)') &
               'WARNING: slope along cDFT search direction is negative:', G_init
            write(stdout,'(a)') '         Resetting conjugate gradients!'
         end if
         !direction_on_grid = -cov_grad_on_grid
         !opt_direction(:) = -local_cdft_gradient(:) ! OK for minimisation
         opt_direction(:) =  local_cdft_gradient(:) ! OK for MAXIMISATION

         cg_count = 0
         G_init = wrappers_ddot(cdft_gradient_size, &
                 local_cdft_gradient,1,opt_direction,1)
         !G_init = wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
         !     contra_grad_on_grid,1,direction_on_grid,1)

!gibo: no need to collect anything since we are working on 
!      the sum-reduced cdft_gradient
         ! cks: collect the work of each node
         !call comms_reduce('SUM', G_init)

         !gibo: for cDFT, act against NEGATIVE slope along search direction
         !if (G_init > 0.0_DP) then
         if (G_init < 0.0_DP) then
            if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then
               write(stdout,'(a)') 'WARNING: slope along cDFT search direction &
                    &is still negative.'
               write(stdout,'(a)') '         Reversing search direction!!'
            end if
            !direction_on_grid = -direction_on_grid
            opt_direction(:) = -opt_direction(:)
            G_init = -G_init
            ! ndmh: if searching 'uphill', (or downhill for cDFT)
            ! ndmh: always re-check final step
            ! ndmh: before accepting, to avoid accepting very bad steps.
            reversing = .true.
         end if

      end if

      call services_flush


    end subroutine internal_find_direction

!##############################################################################
!##############################################################################

   function internal_cdft_polak_cg_coeff()
 
     !===============================================================!
     ! This subroutine calculates the conjugate gradient coefficient !
     ! according to the Polak formula:                               !
     !  b_{i+1}=g_{i}^T * [g_{i}-g{i-1}]/[G_{i-1}^T * G_{i-1}]       !
     !                                                               !
     ! for the cDFT optimisation of the U-opt (1D) array             !
     !---------------------------------------------------------------!
     ! Adapted from services_polak_cg_coeff by Gilberto Teobaldi     !
     ! in December 2011                                              !
     !===============================================================!
 
     use comms, only: pub_on_root, comms_reduce
     use constants, only: DP, stdout, NORMAL
     use rundat, only : pub_output_detail
!gibo-start
     use wrappers, only : wrappers_ddot
!gibo-end
 
     real(kind=DP)  :: internal_cdft_polak_cg_coeff
     real(kind=DP)  :: eps
     real(kind=DP)  :: denominator
 
     eps =epsilon(1.0_DP)
 
     denominator =  wrappers_ddot(cdft_gradient_size,&
              local_cdft_gradient,1,local_cdft_gradient,1)
 
     if ( abs(denominator) .gt. eps ) then
        internal_cdft_polak_cg_coeff = &
           wrappers_ddot(cdft_gradient_size,local_cdft_gradient,1, &
           (local_cdft_gradient - local_previous_cdft_gradient),1) / denominator

!       ! cks: contribution from my node
!       services_polak_cg_coeff = &
!            sum(cov_grad(1: vec_size) &
!            *(contra_grad(1: vec_size) -prev_contra_grad(1: vec_size) ) ) / denominator
!       ! cks: sum of contributions from all nodes
!       call comms_reduce('SUM', services_polak_cg_coeff)
     else
        if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
             write(stdout,'(a/a)') &
             'WARNING: zero denominator in internal_cdft_polak_cg_coeff', &
             '         Setting to zero'
        internal_cdft_polak_cg_coeff = 0.0_DP
     endif
 
     !if ( abs(services_polak_cg_coeff).gt.(2.0_DP) ) then
     if ( abs(internal_cdft_polak_cg_coeff).gt.(2.0_DP) ) then
 
        if (pub_on_root .and. pub_output_detail >= NORMAL) &
             write(stdout,'(a,f11.5,a)') &
             'WARNING: internal_cdft_polak_cg_coeff too large &
             &(',internal_cdft_polak_cg_coeff,'). Setting to zero'
        internal_cdft_polak_cg_coeff = 0.0_DP
     endif
 
   end function internal_cdft_polak_cg_coeff
 
!##############################################################################
!##############################################################################

    subroutine internal_print_search_info

      use rundat, only: pub_cdft_max_grad

      implicit none

      write(stdout,'(/a,f22.14)')   'cDFT RMS gradient                = ',rms_gradient
      if (pub_cdft_max_grad>0.0_DP) write(stdout,'(a,f22.14)') &
           'cDFT MAX gradient                = ',max_gradient
      !if (index(pub_devel_code,'NGWF_FD')>0) then
      !   do ifd=1,nfd
      !      write(stdout,'(a,f22.6)')    'FD trial step length        = ', &
      !           fd_trial_step(ifd)
      !   end do
      !end if
      write(stdout,'(a,f22.6)')    'cDFT Trial step length           = ',trial_length
      write(stdout,'(a,f22.11)')   'cDFT Gradient along search dir.  = ',G_init
      !if (index(pub_devel_code,'NGWF_FD')>0) then
      !   do ifd=1,nfd
      !      write(stdout,'(a,f22.11)')   'Gradient by FD              = ',&
      !           (FFD(ifd)-F0)/fd_trial_step(ifd)
      !   end do
      !end if
      if (check_conjugacy) then
         write(stdout,'(a,f22.11)')   'cDFT Conjugacy Test              = ', &
              previous_dir_dot_g
      end if
      write(stdout,'(a,f22.14)')   'cDFT Functional at step 0        = ',F0
      !if (index(pub_devel_code,'NGWF_FD')>0) then
      !   do ifd=1,nfd
      !      write(stdout,'(a,f22.14)')   'Functional at FD step       = ',FFD(ifd)
      !   end do
      !end if
      write(stdout,'(a,f22.14)')   'cDFT Functional at step 1        = ',F1
      if (trial2) then
         write(stdout,'(a,f22.14)')'cDFT Functional at step 2        = ',F2
         write(stdout,'(a,f22.14)')'cDFT Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'cDFT Rejected quadratic step     = ', &
              quadratic_coeff
         write(stdout,'(a,f22.6)') 'cDFT Selected cubic step         = ',cubic_coeff
      else if (retrial1) then
         write(stdout,'(a,f22.6)') 'cDFT Rejected quadratic step     = ', &
              rejected_quadratic_coeff
         write(stdout,'(a,f22.14)')'cDFT Functional at new step 1    = ',F2
         write(stdout,'(a,f22.14)')'cDFT Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'cDFT Selected quadratic step     = ', &
              quadratic_coeff
      else
         write(stdout,'(a,f22.14)')'cDFT Functional predicted        = ', &
              predicted_functional
         write(stdout,'(a,f22.6)') 'cDFT Selected quadratic step     = ', &
              quadratic_coeff
      endif
      !if ( pub_elec_cg_max > 0 ) then
      if ( pub_cdft_cg_max > 0 ) then
         write(stdout,'(a,f22.6)')    'cDFT Conjugate gradients coeff.  = ',cg_coeff
      endif
      write(stdout,'(a)')'oLooooooooooooooooooooooooo cDFT line search &
           &finished oooooooooooooooooooooooooo'
      write(stdout,'(a)')'                                               '

    end subroutine internal_print_search_info
 
!##############################################################################
!##############################################################################
 
   !gibo: print out intermediate summary of cDFT-optimisation
   !      with atom resolved break-down of the gradient
   subroutine internal_cdft_opt_summary
 
    implicit none

    real(kind=DP), parameter :: zero=0.0_DP
 
    if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
         maxit_cdft_u_cg > 0) then
 
        write(stdout,'(/a)')'oSoooooooooooooooooooooooooooooooooooooo&
             &oooooooooooooooooooooooooooooooooooooooo'
        write(stdout,'(a,i4.3,a,f24.14/)') &
                 & '    cDFT summary, CG iteration: ', iteration, &
                 !& '       Energy (Ha): ', total_energy
                 & '       Energy (Ha): ', total_energy_current
        write(stdout,'(a,f15.9)')  'TOTAL RMS cDFT-GRADIENT(UP+DOWN):', &
               rms_gradient
        write(stdout,'(a,f15.9)')  'MAX cDFT-GRADIENT(UP+DOWN):      ', &
               max_gradient
    endif
    
    if (pub_on_root .and. pub_output_detail >= VERBOSE .and. &
         maxit_cdft_u_cg > 0) then
 
        CHARGE_SPIN_CDFT: if (pub_cdft_atom_charge .OR. pub_cdft_group_charge &
            .OR. pub_cdft_group_charge_diff) then
    
           write(stdout,'(/a)') '  site    label       Uq(UP)    Uq(DOWN) &
           &      Us      dE/dUq(UP)  dE/dUq(DOWN) '
      
           do species = 1, pub_cell%num_hub_species
              do hat = 1 , pub_cell%nat_hub
      
               if (hub%species_number(hat) == species) then
                   g_counter = 2*hat
                   write(stdout,'(i5,7x,a,5(f12.6))') hat, &
                     & h_species(species)%hub_species, &
                     & hub%cdft_u_charge_up(hat)   * HARTREE_IN_EVS, &
                     & hub%cdft_u_charge_down(hat) * HARTREE_IN_EVS, &
                     & zero                                        , &
                     !& hub%cdft_gradient(1,hat), hub%cdft_gradient(2,hat)
                     & local_cdft_gradient(g_counter-1), local_cdft_gradient(g_counter)
               endif
      
              enddo
           enddo


        elseif (pub_cdft_atom_spin .OR. pub_cdft_group_spin &
               .OR. pub_cdft_group_spin_diff) then

         write(stdout,'(/a)') '  site    label       Uq(UP)    Uq(DOWN) &
         &      Us      dE/dUs(UP)  dE/dUs(DOWN) '
     
           do species = 1, pub_cell%num_hub_species
              do hat = 1 , pub_cell%nat_hub
     
               if (hub%species_number(hat) == species) then
                   g_counter = 2*hat
                   write(stdout,'(i5,7x,a,5(f12.6))') hat, &
                     & h_species(species)%hub_species, &
                     & zero                                        , &
                     & zero                                        , &
                     & hub%cdft_u_spin(hat)        * HARTREE_IN_EVS, &
                     !& hub%cdft_gradient(1,hat), hub%cdft_gradient(2,hat)
                     & local_cdft_gradient(g_counter-1), local_cdft_gradient(g_counter)
               endif
     
              enddo
           enddo

        endif CHARGE_SPIN_CDFT
 
    endif

    if (pub_on_root .and. pub_output_detail >= NORMAL .and. &
         maxit_cdft_u_cg > 0) then

        write(stdout,'(a/)')'oooooooooooooooooooooooooooooooooooooooo&
             &oooooooooooooooooooooooooooooooooooooooo'

    endif

    call services_flush
 
   end subroutine internal_cdft_opt_summary

!##############################################################################
!##############################################################################

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
           !write(stdout,'(/20x,a)') '<<<<< CALCULATION SUMMARY >>>>>'
           write(stdout,'(/15x,a)') 'ooooo CDFT OPTIMISATION SUMMARY ooooo'
           do sumrow=1,lastrow
              write(stdout,'(a80)') summary_lines(sumrow)
           end do
        else
           write(stdout,'(a80)') summary_lines(iteration+1)
        endif

     end if

   end subroutine internal_calculation_summary

!##############################################################################
!##############################################################################

    subroutine internal_qc_output

      !====================================================!
      ! This subroutine prints out quality control info in !
      ! a form that can be easily accessed and compared.   !
      !----------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 12/03/2005.    !
      ! Adapted to cDFT by Gilberto Teobaldi on 7/12/2011  !
      !====================================================!

      use comms, only: comms_barrier
      use constants, only: stdout

      if (pub_on_root) then
         write(stdout,'(a1)')' '
         !write(stdout, '(a30,   i9)')'<QC>        [NGWF iterations]:', &
         write(stdout, '(a30,   i9)')'<QC>        [cDFT iterations]:', &
              iteration
         !write(stdout, '(a30,f18.8)')'<QC>           [total_energy]:', &
              !total_energy
         write(stdout, '(a30,f18.8)')'<QC>      [cDFT total_energy]:', &
              total_energy_current
         !write(stdout, '(a30,f18.8)')'<QC>                     [F0]:', &
         write(stdout, '(a30,f18.8)')'<QC>                [cDFT F0]:', &
              F0
         !write(stdout, '(a30,f18.8)')'<QC>                     [F1]:', &
         write(stdout, '(a30,f18.8)')'<QC>                [cDFT F1]:', &
              F1
         !write(stdout, '(a30,f18.8)')'<QC>                     [F2]:', &
         write(stdout, '(a30,f18.8)')'<QC>                [cDFT F2]:', &
              F2
         !write(stdout, '(a30,f18.8)')'<QC>   [predicted_functional]:', &
         write(stdout, '(a30,f18.8)')'<QC>[cDFT predicted_functinl]:', &
              predicted_functional
         !write(stdout,'(a30,f22.12)')'<QC>           [rms_gradient]:', &
         write(stdout,'(a30,f22.12)')'<QC>      [cDFT rms_gradient]:', &
              rms_gradient
         !write(stdout,'(a30,f22.12)')'<QC>                 [G_init]:', &
         write(stdout,'(a30,f22.12)')'<QC>            [cDFT G_init]:', &
              G_init
         !write(stdout, '(a30,f16.6)')'<QC>        [quadratic_coeff]:', &
         write(stdout, '(a30,f16.6)')'<QC>   [cDFT quadratic_coeff]:', &
              quadratic_coeff
         !write(stdout, '(a30,f16.6)')'<QC>            [cubic_coeff]:', &
         write(stdout, '(a30,f16.6)')'<QC>       [cDFT cubic_coeff]:', &
              cubic_coeff
         !write(stdout, '(a30,f16.6)')'<QC>               [cg_coeff]:', &
         write(stdout, '(a30,f16.6)')'<QC>          [cDFT cg_coeff]:', &
              cg_coeff
         !write(stdout, '(a30,f16.6)')'<QC>           [trial_length]:', &
         write(stdout, '(a30,f16.6)')'<QC>      [cDFT trial_length]:', &
              trial_length
      endif

      call comms_barrier


    end subroutine internal_qc_output

!##############################################################################
!##############################################################################

    subroutine internal_cdft_restart_write()

    !==========================================================================!
    ! This subroutine write the cDFT U-potentials                              !
    ! [i.e. the hub%cdft_u_charge_up/down or hub%cdft_u_spin arrays]           !
    ! to (.cdft) file                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Gilberto Teobaldi in December 2011                            !
    !==========================================================================!

      use constants,       only: stdout
      use comms,           only: pub_on_root, comms_abort
      use rundat,          only: pub_rootname
      use simulation_cell, only: pub_cell
      use utils,           only: utils_unit, utils_open_unit_check, &
             utils_close_unit_check

      implicit none

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: species, hat, g_counter
      character(len=80) :: filename
      real(KIND=DP), parameter :: zero= 0.0_DP

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering internal_cdft_restart_write'
#endif

      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.cdft'

         open(unit=io_unit,iostat=io_stat,file=filename,status='REPLACE',&
              form='FORMATTED',position='REWIND',action='WRITE')
         call utils_open_unit_check('internal_cdft_restart_write',&
                                     filename,io_stat)

         write(io_unit,'(a/)') '  site    label       l        Z        Uq(UP) &
                &    Uq(DOWN)       Us         Nup        Ndown      Nup-Ndown'

        CHARGE_SPIN_CDFT: if (pub_cdft_atom_charge .OR. pub_cdft_group_charge &
            .OR. pub_cdft_group_charge_diff) then

          do species = 1, pub_cell%num_hub_species
            do hat = 1 , pub_cell%nat_hub

             if (hub%species_number(hat) == species) then
                 g_counter = 2*hat
                 write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,     &
                   & h_species(species)%hub_species,               &
                   & h_species(species)%hub_ang_mom,               &
                   & h_species(species)%hub_charge,                &
                   & hub%cdft_u_charge_up(hat)   * HARTREE_IN_EVS, &
                   & hub%cdft_u_charge_down(hat) * HARTREE_IN_EVS, &
                   !& hub%cdft_u_spin(hat)        * HARTREE_IN_EVS, &
                   & zero,                                         &
                   & h_species(species)%cdft_target_up,            &
                   & h_species(species)%cdft_target_down,          &
                   & h_species(species)%cdft_target_spin
             endif

            enddo
          enddo

        elseif (pub_cdft_atom_spin .OR. pub_cdft_group_spin &
               .OR. pub_cdft_group_spin_diff) then      

          do species = 1, pub_cell%num_hub_species
            do hat = 1 , pub_cell%nat_hub

             if (hub%species_number(hat) == species) then
                 g_counter = 2*hat
                 write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,     &
                   & h_species(species)%hub_species,               &
                   & h_species(species)%hub_ang_mom,               &
                   & h_species(species)%hub_charge,                &
                   & zero,                                         &
                   & zero,                                         &
                   !& hub%cdft_u_charge_up(hat)   * HARTREE_IN_EVS, &
                   !& hub%cdft_u_charge_down(hat) * HARTREE_IN_EVS, &
                   & hub%cdft_u_spin(hat)        * HARTREE_IN_EVS, &
                   & h_species(species)%cdft_target_up,            &
                   & h_species(species)%cdft_target_down,          &
                   & h_species(species)%cdft_target_spin
             endif

            enddo
          enddo

        endif CHARGE_SPIN_CDFT

         close(io_unit,iostat=io_stat)
         call utils_open_unit_check('internal_cdft_restart_write',&
                                     filename,io_stat)

      endif

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving internal_cdft_restart_write'
#endif

    end subroutine internal_cdft_restart_write

!!##############################################################################
!!##############################################################################
!
!    subroutine cdft_restart_read()
!
!    !==========================================================================!
!    ! This subroutine reads the cDFT U-potential from (.cdft) file and update  !
!    ! the hub%cdft_u_charge_up/down or hub%cdft_u_spin arrays                  !
!    !--------------------------------------------------------------------------!
!    ! Written by Gilberto Teobaldi in December 2011                            !
!    !==========================================================================!
!
!      use constants,       only: stdout
!      use comms,           only: pub_on_root, comms_abort
!!gibo-start
!      use hubbard_build,   only: HUBBARD_MODEL
!      use hubbard_init,    only: h_species
!!gibo-end
!      use rundat,          only: pub_rootname
!      use simulation_cell, only: pub_cell
!      use utils,           only: utils_unit, utils_open_unit_check, &
!             utils_close_unit_check
!
!      implicit none
!
!      ! Local variables
!      integer           :: io_unit,io_stat
!      integer           :: species, hat, g_counter, dummy_int
!      character(len=80) :: filename
!      !real(KIND=DP), parameter :: zero= 0.0_DP
!      real(KIND=DP) :: zero
!
!#ifdef DEBUG
!      if (pub_on_root) write(stdout,'(/a)') &
!         'DEBUG: Entering cdft_restart_read'
!#endif
!
!      if (pub_on_root) then
!
!         ! Find available unit specifier
!         io_unit = utils_unit()
!         write(filename,'(a,a)') trim(pub_rootname),'.cdft'
!
!         open(unit=io_unit,iostat=io_stat,file=filename,&
!              access='SEQUENTIAL',form='FORMATTED',position='REWIND', &
!              action='READ')
!
!
!         RESTART_OK: if (io_stat == 0) then
!
!!         write(io_unit,'(a/)') '  site    label       l        Z        Uq(UP) &
!!                &    Uq(DOWN)       Us         Nup        Ndown      Nup-Ndown'
!
!           read(io_unit,'(/,/)')  ! skip the 2-line header of .cdft file
!
!           CHARGE_SPIN_CDFT: if (pub_cdft_atom_charge .OR. &
!                   pub_cdft_group_charge .OR. pub_cdft_group_charge_diff) then
!   
!             do species = 1, pub_cell%num_hub_species
!               do hat = 1 , pub_cell%nat_hub
!   
!                if (hub%species_number(hat) == species) then
!                    g_counter = 2*hat
!                    !write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,     &
!                    read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
!                      & h_species(species)%hub_species,                 &
!                      & h_species(species)%hub_ang_mom,                 &
!                      & h_species(species)%hub_charge,                  &
!                      & hub%cdft_u_charge_up(hat),                      &
!                      & hub%cdft_u_charge_down(hat),                    &
!                      !& hub%cdft_u_charge_up(hat)   * HARTREE_IN_EVS,   &
!                      !& hub%cdft_u_charge_down(hat) * HARTREE_IN_EVS,   &
!                      & zero,                                           &
!                      & h_species(species)%cdft_target_up,              &
!                      & h_species(species)%cdft_target_down,            &
!                      & h_species(species)%cdft_target_spin
!
!                      ! internally convert U-potentials in Hartree
!                      hub%cdft_u_charge_up(hat)   = &
!                          hub%cdft_u_charge_up(hat)/HARTREE_IN_EVS
!                      hub%cdft_u_charge_down(hat) = &
!                          hub%cdft_u_charge_down(hat)/HARTREE_IN_EVS
!                endif
!   
!               enddo
!             enddo
!
!             !broadcast updated U-potentials
!             call comms_bcast(pub_root_node_id, hub%cdft_u_charge_up)
!             call comms_bcast(pub_root_node_id, hub%cdft_u_charge_down)
!   
!           elseif (pub_cdft_atom_spin .OR. pub_cdft_group_spin &
!                  .OR. pub_cdft_group_spin_diff) then
!   
!             do species = 1, pub_cell%num_hub_species
!               do hat = 1 , pub_cell%nat_hub
!   
!                if (hub%species_number(hat) == species) then
!                    g_counter = 2*hat
!                    !write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,     &
!                    read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
!                      & h_species(species)%hub_species,                 &
!                      & h_species(species)%hub_ang_mom,                 &
!                      & h_species(species)%hub_charge,                  &
!                      & zero,                                           &
!                      & zero,                                           &
!                      & hub%cdft_u_spin(hat),                           &
!                      !& hub%cdft_u_spin(hat)        * HARTREE_IN_EVS,   &
!                      & h_species(species)%cdft_target_up,              &
!                      & h_species(species)%cdft_target_down,            &
!                      & h_species(species)%cdft_target_spin
!
!                      ! internally convert U-potentials in Hartree
!                      hub%cdft_u_spin(hat)   = &
!                          hub%cdft_u_spin(hat)/HARTREE_IN_EVS
!                endif
!   
!               enddo
!             enddo
!
!             !broadcast updated U-potentials
!             call comms_bcast(pub_root_node_id, hub%cdft_u_spin)
!   
!           endif CHARGE_SPIN_CDFT
!
!           close(unit=io_unit,iostat=io_stat)
!           call utils_close_unit_check('cdft_restart_read',filename,io_stat)
!
!         else
!           write(stdout,'(/,3a)') 'WARNING : cdft_restart_read could &
!             &not open "', trim(filename),'"'
!           write(stdout,'(/a)') '...aborting...'
!           call comms_abort
!         endif RESTART_OK
!
!      endif
!
!#ifdef DEBUG
!      if (pub_on_root) write(stdout,'(/a)') &
!         'DEBUG: Leaving cdft_restart_read'
!#endif
!
!    end subroutine cdft_restart_read

!##############################################################################
!##############################################################################

   subroutine internal_cdft_kernel_reinit()

   implicit none

   ! STEP_1 (probably redundant: it is done at the end of ngwf_cg_optimise)
    ! cks: calculate density INDEPENDENT matrix elements with current NGWFs
    call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
         hub_proj_basis, hub, val_rep, val_ngwf_basis) ! 2x optional args

   ! STEP_2: re-initialise ham and denskern
   !         [by Palser Manolopulos Canonical Purification
    do is=1,pub_cell%num_spins
        if (pub_on_root) then
           if (is == 1) then
              write(stdout,'(a)',advance='no') &
                   'cDFT: Up spin density kernel initialisation ...'
           elseif (is == 2) then
              write(stdout,'(a)',advance='no') &
                   'cDFT: Down spin density kernel initialisation ...'
           endif
        endif

        ! ham = Hartree + locpot
        call sparse_copy(ham%ham(is),ham%lhxc(is))

        ! ham = (Hartree + locpot) + kinet
        call sparse_axpy(ham%ham(is),rep%kinet,1.0_DP)

        ! ham = (Hartree + locpot + kinet) + nonlocpot
        if (pub_any_nl_proj.and.(.not.pub_aug)) then
           call sparse_axpy(ham%ham(is),rep%nonlocpot(1),1.0_DP)
        end if
        if (pub_aug) then
           call sparse_axpy(ham%ham(is),ham%nonlocpot(is),1.0_DP)
        end if

        if (pub_maxit_palser_mano < 0) then
           ! Diagonalise the initial Hamiltonian
           call kernel_init_core_ham(denskern(is), &
                ham%ham(is), rep%overlap, rep%n_occ(is))
        else
           ! Palser-Manolopoulos Canonical Purification
           call palser_mano_kernel_optimise(denskern(is), &
                ham%ham(is), rep%overlap, rep%inv_overlap, rep%n_occ(is), &
                num_iter=pub_maxit_palser_mano)
        endif

        if (pub_on_root) write(stdout,'(a)') '... done'
    enddo

     ! cks: output density kernel to file if this is requested
     if (write_denskern.and.(.not.pub_write_converged_dk_ngwfs)) &
         ! cDFT is *not* a conduction_NGWF call therefore print valence_denskern
         call restart_kernel_write(denskern,.FALSE.) 

   end subroutine internal_cdft_kernel_reinit

!##############################################################################
!##############################################################################


  end subroutine cdft_u_cg_optimise

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    subroutine cdft_update_matrices(hub, hub_proj_basis)

    !==========================================================!
    ! This subroutine updates the cDFT[+U] SPAM3 matrices      !
    !----------------------------------------------------------!
    ! Written for cDFT-OPT module,and based on                 !
    ! hubbard_build_matrices by Gilberto Teobaldi in Dec. 2011 !
    !==========================================================!

    use comms, only: pub_my_node_id, pub_on_root
!gibo-start
   use constants, only: DP, stdout, VERBOSE
!gibo-end
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
!gibo-start
    use hubbard_build, only: HUBBARD_MODEL
!gibo-end
    use parallel_strategy, only:  pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom, pub_first_atom_on_node
!gibo-start-10.12.11
!    use rundat, only: pub_aug, pub_cdft, pub_cdft_hubbard, &
!         task, pub_hubbard_restart, pub_hubbard_atomsolve
    use rundat, only: pub_cdft, &
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge, pub_cdft_group_spin, &
        pub_cdft_group_charge_diff, pub_cdft_group_spin_diff
!gibo-end-10.12.11
    use simulation_cell, only: pub_cell
!gibo-start
    !use sparse, only: SPAM3, sparse_create, &
    !     sparse_scale, sparse_put_element
    use sparse, only: SPAM3, sparse_put_element
!gibo-end

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: sp, hat_on_node, hat, theatom
    integer :: hub_proj
    real(kind=DP) :: u_charge_up, u_charge_down, u_spin_up, u_spin_down

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_build_matrices'
#endif

!====select CDFT-mode and iterate...
    ! gibo: ATOM/GROUP-CHARGE-[DIFF] constrained run
    CDFT_MODE: if ((pub_cdft_atom_charge .OR. pub_cdft_group_charge).OR.&
                   pub_cdft_group_charge_diff) then
!====select CDFT-mode and iterate...

    ! ddor: Loop over Hubbard atoms on my node
    cdft_atoms_1 : do hat_on_node = 1, &
         pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
       sp = hub%species_number(hat)
       theatom = pub_distr_atom(hub%orig(hat))

       u_charge_up   = hub%cdft_u_charge_up(hat)
       u_charge_down = hub%cdft_u_charge_down(hat)

       ! ddor: Loop over Hubbard projectors on my node
       hubprojs_cdft_1: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
            hub_proj_basis%first_on_atom(theatom) + &
            hub_proj_basis%num_on_atom(theatom) - 1

          call sparse_put_element(u_charge_up,&
               hub%up_matrix,hub_proj,hub_proj)
          call sparse_put_element(u_charge_down,&
               hub%down_matrix,hub_proj,hub_proj)

       end do hubprojs_cdft_1

    end do cdft_atoms_1


!====select CDFT-mode and iterate...
    ! gibo: ATOM/GROUP-SPIN-[DIFF] constrained run
    else if ((pub_cdft_atom_spin .OR. pub_cdft_group_spin).OR.&
                   pub_cdft_group_spin_diff) then
!====select CDFT-mode and iterate...

    ! ddor: Loop over Hubbard atoms on my node
    cdft_atoms_2 : do hat_on_node = 1, &
         pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
       sp = hub%species_number(hat)
       theatom = pub_distr_atom(hub%orig(hat))

       ! MIND THE SIGNS: Ecdft_spin = Us*[Nup -Ndown - DN_target]
       ! thus V_up = Us, and V_down = -Us
       u_spin_up   =  hub%cdft_u_spin(hat)
       u_spin_down = -hub%cdft_u_spin(hat)


       ! ddor: Loop over Hubbard projectors on my node
       hubprojs_cdft_2: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
            hub_proj_basis%first_on_atom(theatom) + &
            hub_proj_basis%num_on_atom(theatom) - 1

          call sparse_put_element(u_spin_up,&
               hub%up_matrix,hub_proj,hub_proj)
          call sparse_put_element(u_spin_down,&
               hub%down_matrix,hub_proj,hub_proj)

       end do hubprojs_cdft_2

    end do cdft_atoms_2

!====select CDFT-mode and iterate...
    end if CDFT_MODE
!====select CDFT-mode and iterate...


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_build_matrices'
#endif
    end subroutine cdft_update_matrices

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    subroutine cdft_restart_read(hub)

    !==========================================================================!
    ! This subroutine reads the cDFT U-potential from (.cdft) file and update  !
    ! the hub%cdft_u_charge_up/down or hub%cdft_u_spin arrays                  !
    !--------------------------------------------------------------------------!
    ! Written by Gilberto Teobaldi in December 2011                            !
    !==========================================================================!

!gibo-start
      use constants,       only: DP, stdout, VERBOSE, HARTREE_IN_EVS
!gibo-end
      use comms,           only: pub_on_root, pub_root_node_id, comms_bcast, &
                                 comms_abort
!gibo-start
      use hubbard_build,   only: HUBBARD_MODEL
      use hubbard_init,    only: h_species
!gibo-end
!gibo-start
      use rundat, only: pub_rootname, pub_cdft, &
!gibo-start
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge, pub_cdft_group_spin, &
        pub_cdft_group_charge_diff, pub_cdft_group_spin_diff
!gibo-end
      use simulation_cell, only: pub_cell
      use utils,           only: utils_unit, utils_open_unit_check, &
                                 utils_close_unit_check

      implicit none

      ! Arguments
      type(HUBBARD_MODEL), intent(out) :: hub

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: species, hat, g_counter, dummy_int
      character(len=80) :: filename
      !real(KIND=DP), parameter :: zero= 0.0_DP
      real(KIND=DP) :: zero

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering cdft_restart_read'
#endif

      if (pub_on_root) then

         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.cdft'

         open(unit=io_unit,iostat=io_stat,file=filename,&
              access='SEQUENTIAL',form='FORMATTED',position='REWIND', &
              action='READ')


         RESTART_OK: if (io_stat == 0) then

!         write(io_unit,'(a/)') '  site    label       l        Z        Uq(UP) &
!                &    Uq(DOWN)       Us         Nup        Ndown      Nup-Ndown'

           read(io_unit,'(/,/)')  ! skip the 2-line header of .cdft file

           CHARGE_SPIN_CDFT: if (pub_cdft_atom_charge .OR. &
                   pub_cdft_group_charge .OR. pub_cdft_group_charge_diff) then
   
             do species = 1, pub_cell%num_hub_species
               do hat = 1 , pub_cell%nat_hub
   
                if (hub%species_number(hat) == species) then
                    g_counter = 2*hat
                    !write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,     &
                    read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
                      & h_species(species)%hub_species,                 &
                      & h_species(species)%hub_ang_mom,                 &
                      & h_species(species)%hub_charge,                  &
                      & hub%cdft_u_charge_up(hat),                      &
                      & hub%cdft_u_charge_down(hat),                    &
                      !& hub%cdft_u_charge_up(hat)   * HARTREE_IN_EVS,   &
                      !& hub%cdft_u_charge_down(hat) * HARTREE_IN_EVS,   &
                      & zero,                                           &
                      & h_species(species)%cdft_target_up,              &
                      & h_species(species)%cdft_target_down,            &
                      & h_species(species)%cdft_target_spin

                      ! internally convert U-potentials in Hartree
                      hub%cdft_u_charge_up(hat)   = &
                          hub%cdft_u_charge_up(hat)/HARTREE_IN_EVS
                      hub%cdft_u_charge_down(hat) = &
                          hub%cdft_u_charge_down(hat)/HARTREE_IN_EVS
                endif
   
               enddo
             enddo

             !broadcast updated U-potentials
             call comms_bcast(pub_root_node_id, hub%cdft_u_charge_up)
             call comms_bcast(pub_root_node_id, hub%cdft_u_charge_down)
   
           elseif (pub_cdft_atom_spin .OR. pub_cdft_group_spin &
                  .OR. pub_cdft_group_spin_diff) then
   
             do species = 1, pub_cell%num_hub_species
               do hat = 1 , pub_cell%nat_hub
   
                if (hub%species_number(hat) == species) then
                    g_counter = 2*hat
                    !write(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') hat,     &
                    read(io_unit,'(i5,7x,a,5x,i2,7(f12.6))') dummy_int, &
                      & h_species(species)%hub_species,                 &
                      & h_species(species)%hub_ang_mom,                 &
                      & h_species(species)%hub_charge,                  &
                      & zero,                                           &
                      & zero,                                           &
                      & hub%cdft_u_spin(hat),                           &
                      !& hub%cdft_u_spin(hat)        * HARTREE_IN_EVS,   &
                      & h_species(species)%cdft_target_up,              &
                      & h_species(species)%cdft_target_down,            &
                      & h_species(species)%cdft_target_spin

                      ! internally convert U-potentials in Hartree
                      hub%cdft_u_spin(hat)   = &
                          hub%cdft_u_spin(hat)/HARTREE_IN_EVS
                endif
   
               enddo
             enddo

             !broadcast updated U-potentials
             call comms_bcast(pub_root_node_id, hub%cdft_u_spin)
   
           endif CHARGE_SPIN_CDFT

           close(unit=io_unit,iostat=io_stat)
           call utils_close_unit_check('cdft_restart_read',filename,io_stat)

         else
           write(stdout,'(/,3a)') 'WARNING : cdft_restart_read could &
             &not open "', trim(filename),'"'
           write(stdout,'(/a)') '...aborting...'
           call comms_abort
         endif RESTART_OK

      endif

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving cdft_restart_read'
#endif

    end subroutine cdft_restart_read

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


end module cdft_cg
