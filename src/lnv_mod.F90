! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!               Density kernel optimisation module               !
!                                                                !
! This module implements several variants of the LNV method for  !
! optimising the density kernel K for a fixed set of NGWFs. The  !
! density kernel K is expressed in terms of an auxiliary matrix  !
! L by a generalised purifying transformation:                   !
!   K = f(L) [ 3 L.S.L - 2 L.S.L.S.L ]                           !
! where f(L) is a scalar function of L. The variants differ only !
! in the method used to impose the normalisation constraint.     !
!----------------------------------------------------------------!
! Written by Peter Haynes, 17/11/04                              !
! Based on an original version by Chris-Kriton Skylaris in 2000  !
! with subsequent modifications by Chris-Kriton Skylaris,        !
! Arash Mostofi and Peter Haynes.                                !
! Modified by Nicholas Hine to use the NGWF_REP container for    !
! NGWF-related data and matrices, October 2010.                  !
!================================================================!

module lnv

  implicit none

  private

  public :: lnv_denskernel_optimise_cg

contains


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


  subroutine lnv_denskernel_optimise_cg(aux, purkern, ham, lhxc_fine, mu, &
       total_energy, rep, ngwf_basis, hub_proj_basis, hub, &
       grid, localpseudo_fine, core_density_fine, &
       ewald_energy, elements, lnv_threshold, current_maxit_lnv)

    !======================================================================!
    ! This subroutine implements several variants of the LNV method for    !
    ! optimising the density kernel K for a fixed set of NGWFs. The        !
    ! density kernel K is expressed in terms of an auxiliary matrix L by a !
    ! generalised purifying transformation:                                !
    !   K = f(L) [ 3 L.S.L - 2 L.S.L.S.L ]                                 !
    ! where f(L) is a scalar function of L. The variants differ only in    !
    ! the method used to impose the normalisation constraint.              !
    !                                                                      !
    ! Methods implemented so far with references:                          !
    ! * Original Li-Nunes-Vanderbilt scheme: Phys. Rev. B 47, 10891 (1993) !
    ! * Millam-Scuseria variant: J. Chem. Phys. 106, 5569 (1997)           !
    !----------------------------------------------------------------------!
    ! Arguments:                                                           !
    ! aux           (inout): Auxiliary density kernel in SPAM3 format      !
    ! purkern       (inout): Purified density kernel in SPAM3 format       !
    ! ham           (inout): Hamiltonian wrapper (contains SPAM3 matrices) !
    ! lhxc_fine       (out): Local-Hartree-Exhange-Correlation potential   !
    ! mu              (out): Chemical potential type parameter             !
    ! total_energy    (out): Total energy                                  !
    ! rep              (in): NGWF Representation (functions and matrices)  !
    ! ngwf_basis       (in): Function basis type describing the NGWFs.     !
    ! hub_proj_basis   (in): Function basis type describing Hubbard projs  !
    ! grid             (in): The GRID_INFO for all the whole-cell arrays   !
    ! localpseudo_fine (in): Local pseudopotential on fine grid            !
    ! core_density_fine(in): Core Charge on fine grid                      !
    ! ewald_energy     (in): Ewald energy                                  !
    ! lnv_threshold    (in): Convergence threshold for LNV gradient        !
    ! current_maxit_lnv(in): Maximum LNV iterations set in electronic mod  !
    !----------------------------------------------------------------------!
    ! This version written by Peter Haynes, November 2004.                 !
    ! Based on an original version by Chris-Kriton Skylaris in 2000 with   !
    !   subsequent modifications by Chris-Kriton Skylaris, Arash Mostofi   !
    !   and Peter Haynes.                                                  !
    ! Modified by Peter Haynes for parallel SPAM 2, July 2006.             !
    ! Modified to include Nonlinear Core Corrections by Nicholas Hine in   !
    ! January 2009                                                         !
    ! DFT+U added by David O'Regan, April 2009                             !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009            !
    ! Modified by Laura Ratcliff for conduction calculations, Oct 2010     !
    !======================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root, pub_root_node_id, comms_bcast
    use constants, only: DP, max_spins, stdout, VERBOSE, NORMAL, paw_en_size, &
         paw_en_ehart, paw_en_exc, paw_en_exc_dc, paw_en_etxc, paw_en_etxc_dc, &
         paw_en_dijhat, paw_en_dijxc, paw_en_exc_core
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_build_matrix, hamiltonian_energy_components
    use hubbard_build, only: HUBBARD_MODEL, hubbard_energy_total
    use ion, only: ELEMENT
    use kernel, only: pub_ks, kernel_rms_commutator, pub_ksk, &
         kernel_validate_ks, kernel_workspace_invalidate, kernel_normalise, &
         kernel_workspace_allocate, kernel_purify, kernel_fix, &
         kernel_rescale, kernel_validate_ksk, kernel_workspace_deallocate, &
         kernel_rms_err, kernel_occupancy_bounds
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use restart, only : restart_kernel_write
    use rundat, only: pub_any_nl_proj, pub_output_detail, pub_hubbard, &
         old_lnv, lnv_cg_type, write_denskern, pub_usehfx, exact_lnv, &
         pub_lnv_check_trial_steps, lnv_cg_max_step, pub_hub_calculating_u, &
         pub_write_converged_dk_ngwfs, pub_paw, pub_devel_code, &
         pub_cond_calculate
!CW
    use rundat, only : pub_dmft_spoil_kernel,pub_dmft_fully_sc,pub_dmft_fully_sc_h
!END CW
    use services, only: services_flush
    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_rms_element, sparse_trace, &
         sparse_num_element, sparse_create, sparse_product, sparse_copy, &
         sparse_destroy, sparse_axpy, sparse_scale, sparse_transpose
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! ** Arguments
    ! Auxiliary density kernel:
    type(SPAM3), intent(inout) :: aux(pub_cell%num_spins)
    ! Purified density kernel:
    type(SPAM3), intent(inout) :: purkern(pub_cell%num_spins)
    ! Hamiltonian
    type(NGWF_HAM), intent(inout) :: ham
    type(GRID_INFO), intent(inout) :: grid
    ! Local PSP, Hartree and XC potential on grid:
    real(kind=DP), intent(inout) :: lhxc_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(out) :: mu(max_spins)    ! Chemical potential
    real(kind=DP), intent(out) :: total_energy     ! Total energy

    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(in) :: &                 ! Local PSP on grid
         localpseudo_fine(grid%ld1,grid%ld2,grid%max_slabs12)
    real(kind=DP), intent(in) :: &                 ! Core Charge grid
         core_density_fine(grid%ld1,grid%ld2,grid%max_slabs12)
    integer, intent(in) :: current_maxit_lnv       ! Max number of iterations
    real(kind=DP), intent(in) :: ewald_energy      ! Ewald energy
    real(kind=DP), intent(in) :: lnv_threshold     ! LNV threshold
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! ** Local variables
    ! Sparse matrices
    type(SPAM3), allocatable :: trial_aux(:) ! Trial auxiliary matrix

    ! Conjugate gradient variables
    type(SPAM3), allocatable :: co_gradient(:)       ! Covariant gradient
    type(SPAM3), allocatable :: con_gradient(:)      ! Contravariant gradient
    type(SPAM3), allocatable :: con_direction(:)     ! Contra search direction
    type(SPAM3), allocatable :: old_con_direction(:) ! Previous contra direction
    type(SPAM3), allocatable :: old_co_gradient(:)   ! Previous co gradient
    real(kind=DP) :: ne_pure          ! Electron number from purified kernel
    real(kind=DP) :: commutator       ! RMS value of [L,H] commutator
    real(kind=DP) :: rms_gradient     ! RMS value of gradient of E wrt L
    real(kind=DP) :: line_search_step ! Line search optimal step length
    real(kind=DP) :: cg_coeff         ! Conjugate gradients coefficient
    real(kind=DP) :: lhxc_energy      ! Local PSP, Hartree and XC energy
    real(kind=DP) :: mu_term          ! Term associated with chemical potential
    real(kind=DP) :: rms_dkern        ! RMS value of change in density kernel
    real(kind=DP) :: old_g_dot_g      ! Previous scalar product of gradients
    real(kind=DP) :: new_g_dot_g      ! Current scalar product of gradients
    integer :: iteration              ! Iteration number
    integer :: cg_count               ! Number of CG steps since last reset
    integer, parameter :: cg_max=5    ! Maximum number of CG steps before reset
    logical :: converged              ! Flag to indicate convergence

    ! Line search variables
    real(kind=DP) :: Qinitial         ! Initial value of function
    real(kind=DP) :: Qslope           ! Initial slope along search direction
    real(kind=DP) :: &                ! Trial step length
         trial_step=0.5_DP
    real(kind=DP) :: Qtrial           ! Function value at trial step
    real(kind=DP) :: optimal_step     ! Optimal step length
    real(kind=DP) :: Qoptimal         ! Function value at optimal step
    real(kind=DP) :: Qpredict         ! Predicted value at optimal step
    real(kind=DP) :: Qminimum         ! Predicted value of function at minimum
    real(kind=DP) :: total_energy_at_trial_length
    real(kind=DP) :: hubbard_energy   ! Hubbard correction energy

    ! General variables
    real(kind=DP), parameter :: eps=epsilon(1.0_DP)
    real(kind=DP) :: norm_fac(max_spins) ! Normalisation factor
    real(kind=DP) :: spin_fac         ! Spin factor
    integer :: is                     ! Spin counter
    integer :: ierr                   ! Error flag

    ! ndmh: for finite differencing
    integer, parameter :: nfd=10
    real(kind=DP) :: fd_trial_step(nfd)=(/0.00001_DP,0.0001_DP,0.001_DP, &
         0.01_DP,0.1_DP,0.3_DP,0.4_DP,0.5_DP,0.6_DP,1.0_DP/)
    integer :: ifd
    real(kind=DP) :: QFD(nfd)


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: &
         &Entering lnv_denskernel_optimise_cg'
#endif

    ! Flush output
    call services_flush

    ! Start the timer
    call timer_clock('lnv_denskernel_optimise_cg', 1)

    call kernel_workspace_allocate(aux,rep%overlap)

    ! Allocate matrices with the sparsity pattern of the density kernel
    allocate(co_gradient(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('lnv_denskernel_optimise_cg','co_gradient',ierr)
    allocate(con_gradient(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('lnv_denskernel_optimise_cg','con_gradient',ierr)
    allocate(con_direction(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('lnv_denskernel_optimise_cg','con_direction',ierr)
    allocate(old_con_direction(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('lnv_denskernel_optimise_cg', &
         'old_con_direction',ierr)

    if(lnv_cg_type/='LNV_FLETCHER')then
       allocate(old_co_gradient(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('lnv_denskernel_optimise_cg','old_co_gradient',&
            ierr)
    endif

    allocate(trial_aux(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('lnv_denskernel_optimise_cg','trial_aux',ierr)

    do is=1,pub_cell%num_spins
       call sparse_create(co_gradient(is),aux(is))
       call sparse_create(con_gradient(is),aux(is))
       call sparse_create(con_direction(is),aux(is))
       call sparse_create(old_con_direction(is),aux(is))
       if(lnv_cg_type/='LNV_FLETCHER')then
          call sparse_create(old_co_gradient(is),aux(is))
       endif
       call sparse_create(trial_aux(is),aux(is))
    end do

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! Ensure the kernel will be stable before optimisation starts

!CW
! CW COMMENT : HERE MIGHT SPOIL THE DMFT DENS KERNEL FLAG MIGHT BE NEEDED
  write(*,*) '------ NORMALISE KERNEL ------'
  write(*,*) '------ KERNEL FIX ------------'
!END CW

    call kernel_normalise(aux, rep%overlap, rep%inv_overlap, rep%n_occ)
    call kernel_fix(aux, rep%overlap, rep%inv_overlap, miniter=3)

    ! Initialise variables for conjugate gradients
    old_g_dot_g = 1.0e100_DP           ! Set previous product to dummy value
    line_search_step = 0.0_DP          ! Set line search step length to zero
    cg_count = 0                       ! Reset conjugate gradients counter
    converged = .false.                ! Initially no convergence
    spin_fac = 2.0_DP / pub_cell%num_spins

    if (pub_on_root .and. pub_output_detail >= NORMAL) then
       write(stdout,'(a)') '.........................................&
            &.......................................'
       if (exact_lnv) then
          write(stdout,'(a)') '<<<<<<<<<<<<<< LNV (Original version) &
               &density kernel optimisation >>>>>>>>>>>>>>'
       else
          write(stdout,'(a)') '<<<<<<<<<< LNV (Millam-Scuseria variant) &
               &density kernel optimisation >>>>>>>>>>>'
       end if
       write(stdout,'(a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~&
            &~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

       write(stdout,'(a)')'iter|       energy      |  rms DeltaK |&
            &  commutator |LScoef|CGcoef|      Ne    |'
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ** Start of main loop
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    do iteration=1,current_maxit_lnv

       ! Flush output
       call services_flush

       ! cks: ======================= ENERGY AT POINT 0 =======================
!CW
       if(pub_dmft_spoil_kernel.and.iteration==1)then
         total_energy = internal_energy(aux,.true.,spoil_force=.true.)
       else
         total_energy = internal_energy(aux,.true.)
       endif
!END CW
       ! cks: ==================== END ENERGY AT POINT 0 ======================


       ! cks: ******** CALCULATE KOHN-SHAM MATRIX AND MU **********************

       ! ddor: Hamiltonian renewal is not required when calculating DFT+U parameter
       if ((.not. pub_hub_calculating_u) .or. ((pub_hub_calculating_u).and. &
            (iteration .eq. 1))) then
          ! lr408: no need to recalculate ham for a conduction calculation
              if (.not.pub_cond_calculate) then
!CW
! if condition added in from of call ham....
                if(.not.pub_dmft_fully_sc.or.(pub_dmft_fully_sc.and.pub_dmft_fully_sc_h)) call hamiltonian_build_matrix(ham, rep)
!END CW
              end if
          ! Calculate mu
          if (exact_lnv) then
             call internal_lnv_mu_ideal()
          else
             call internal_ms_mu_ideal()
          end if

          do is=1,pub_cell%num_spins
             if (pub_on_root .and. (pub_output_detail == VERBOSE)) &
                  write(stdout,'(a,i1,a,f20.12)') '      ideal mu',is,': ', mu(is)
          end do
       endif

       ! cks: **** END CALCULATE KOHN-SHAM MATRIX AND MU **********************


       ! cks: %%%%%%%%%%%%%%%%%%%%%%%% GRADIENT AT POINT 0 %%%%%%%%%%%%%%%%%%%%
       rms_gradient = internal_gradient_norm()
       ! pdh: hack to make consistent with spin-unpolarised case
       rms_gradient = rms_gradient * pub_cell%num_spins
       ! cks: %%%%%%%%%%%%%%%%%%%%%%% END GRADIENT AT POINT 0 %%%%%%%%%%%%%%%%%


       ! cks: ########################## TEST CONVERGENCE #####################
       converged = internal_test_convergence()
       call comms_bcast(pub_root_node_id,converged)
       if (converged) exit
       ! cks: ######################## END TEST CONVERGENCE ###################

       ! Line search is skipped in the last iteration
       if (iteration /= current_maxit_lnv) then

          ! cks: ===========================LINE SEARCH========================
          ! Initial function value including constraint term
          mu_term = internal_mu_term(aux)
          Qinitial = total_energy + mu_term

          if ( (pub_output_detail == VERBOSE) .and. pub_on_root) &
             write(stdout,'(a,f20.12)') '             Q0: ', Qinitial

          call internal_search_direction

          if ( (pub_output_detail == VERBOSE) .and. pub_on_root) &
             write(stdout,'(a,f20.12)') '   RMS gradient: ', rms_gradient

          ! Sophisticated line search is not suitable initially
          if (iteration > 4 .or. rms_gradient <= 0.05_DP) then

             if (index(pub_devel_code,'LNV_FD')>0) then

                if ( (pub_output_detail < VERBOSE) .and. pub_on_root) &
                     write(stdout,'(a,f20.12)') '             G0: ',Qslope

                ! ndmh: evaluate and print energy at several FD trial steps
                do ifd=1,nfd

                   if (pub_lnv_check_trial_steps) then
                      call internal_kernel_protected_step(trial_aux, &
                           fd_trial_step(ifd),aux,con_direction)
                   else
!
                      ! Take a trial step along the search direction
                      do is=1,pub_cell%num_spins
                         call sparse_copy(trial_aux(is),aux(is))
                         call sparse_axpy(trial_aux(is),con_direction(is), &
                              fd_trial_step(ifd))
                      end do
                      call kernel_workspace_invalidate()

                   end if

                   ! ndmh: %%%%%%%%%%%%%%%%%% ENERGY AT FD LENGTH %%%%%%%%%%%%%%
                   total_energy_at_trial_length = internal_energy(trial_aux, &
                        .false.)

                   mu_term = internal_mu_term(trial_aux)

                   QFD(ifd) = total_energy_at_trial_length + mu_term
!CW
                   if ((pub_output_detail>=VERBOSE).and. &
                        (.not.pub_cond_calculate)) then
!END CW
                      do is=1,pub_cell%num_spins
                         call sparse_scale(purkern(is),norm_fac(is)*spin_fac)
                      end do
                      call hamiltonian_energy_components( &
                           purkern, rep, localpseudo_fine, core_density_fine, &
                           ngwf_basis, hub_proj_basis, hub, ewald_energy, &
                           ham%hfexchange)
                      do is=1,pub_cell%num_spins
                         call sparse_scale(purkern(is), &
                              1.0_DP/(norm_fac(is)*spin_fac))
                      end do
                   end if

                   if (pub_on_root) then
                      write(stdout,'(a,2f20.12)') ' FD step, QFD  : ', &
                           fd_trial_step(ifd),QFD(ifd)
                   end if
                   ! ndmh: %%%%%%%%%%%%%%%%%% END ENERGY AT FD LENGTH %%%%%%%%%%

                end do
             end if

             if (pub_lnv_check_trial_steps) then
                call internal_kernel_protected_step(trial_aux,trial_step, &
                     aux,con_direction)
             else

                ! Take a trial step along the search direction
                do is=1,pub_cell%num_spins
                   call sparse_copy(trial_aux(is),aux(is))
                   call sparse_axpy(trial_aux(is),con_direction(is),trial_step)
                end do
                call kernel_workspace_invalidate()

             end if

             ! cks: %%%%%%%%%%%%%%%%%%%%% ENERGY AT TRIAL LENGTH %%%%%%%%%%%%%%
             total_energy_at_trial_length = internal_energy(trial_aux,.false.)

             mu_term = internal_mu_term(trial_aux)

             Qtrial = total_energy_at_trial_length + mu_term

             if ( (pub_output_detail == VERBOSE) .and. pub_on_root) then
                  write(stdout,'(a,f20.12)') '     Trial step: ', trial_step
                  write(stdout,'(a,f20.12)') '             Q1: ', Qtrial
             end if
             ! cks: %%%%%%%%%%%%%%%%%%% END ENERGY AT TRIAL LENGTH %%%%%%%%%%%%


             ! Calculate the optimal step by fitting a parabola to the two
             ! function evaluations and initial slope
             optimal_step = internal_quadratic_step()


             ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
             ! If the predicted step is negative, the parabolic fit is poor,
             !     so take a second trial step and fit a cubic instead
             ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

             if (optimal_step < 0.0_DP) then

                if ((pub_output_detail == VERBOSE) .and. pub_on_root) &
                     write(stdout,'(a)')'         Quadratic fit unsuccessful: &
                     &proceeding to cubic fit'

                optimal_step = 2.0_DP * trial_step


                ! cks: %%%%%%%%%%%%%%%%%%% ENERGY AT TRIAL LENGTH 2 %%%%%%%%%%%
                if (pub_lnv_check_trial_steps) then
                   call internal_kernel_protected_step(trial_aux,optimal_step, &
                      aux,con_direction)
                else
                   ! Take a second trial step along the search direction
                   do is=1,pub_cell%num_spins
                      call sparse_copy(trial_aux(is),aux(is))
                      call sparse_axpy(trial_aux(is),con_direction(is), &
                           optimal_step)
                   end do
                   call kernel_workspace_invalidate()
                end if

                ! Calculate the energy
                total_energy_at_trial_length = internal_energy(trial_aux,.false.)

                mu_term = internal_mu_term(trial_aux)

                Qoptimal = total_energy_at_trial_length + mu_term

                ! cks: %%%%%%%%%%%%%%%% END ENERGY AT TRIAL LENGTH 2 %%%%%%%%%%

                if ( (pub_output_detail == VERBOSE) .and. pub_on_root) &
                     write(stdout,'(a,f20.12)') '             Q2: ',Qoptimal


                ! Recalculate the optimal step by fitting a cubic to the three
                ! function evaluations and initial slope

                optimal_step = internal_cubic_step()

             else

                ! ndmh: re-ordered next section a bit for improved readability
                ! ndmh: of warnings and supression of warnings at BRIEF level
                if ((pub_output_detail == VERBOSE) .and. pub_on_root) &
                     write(stdout,'(a,f20.12)') '           Qmin: ',Qpredict
                if (optimal_step > lnv_cg_max_step) then

                   ! cks: prevent exceedingly long quadratic steps
                   if (pub_on_root .and. (pub_output_detail >= NORMAL)) then
                      write(stdout,'(a)') 'WARNING in &
                           &lnv_denskernel_optimise_cg: setting quadratic &
                           &optimal_step to safe'
                      write(stdout,'(a,f16.12,a)') 'value (Calculated optimal &
                           &quadratic step=',optimal_step,')'
                   end if
                   cg_coeff = 0.0_DP
                   cg_count = 0
                   optimal_step = 0.15_DP
                elseif ( (pub_output_detail == VERBOSE) .and. pub_on_root) then
                   ! cks: If quadratic step accepted, print verbose info
                   write(stdout,'(a,f16.12,a)') 'Quadratic minimum step length &
                        &(',optimal_step,') accepted'
                end if

             end if


             ! Update trial step for next line minimisation
             if (optimal_step > 0.0_DP) &
                  trial_step = max(sqrt(trial_step*optimal_step),eps)

             line_search_step = optimal_step

          else

             line_search_step = 0.1_DP

          end if

          ! Update the auxiliary matrix
          if (pub_lnv_check_trial_steps) then
              optimal_step = line_search_step
              call internal_kernel_protected_step(trial_aux,line_search_step, &
                   aux,con_direction)
              do is=1,pub_cell%num_spins
                 call sparse_copy(aux(is),trial_aux(is))
              end do
              ! ndmh: reset CG if line step protection truncated the step, since
              ! ndmh: we have not fully minimised in that search direction
              if (line_search_step < optimal_step) then
                 cg_coeff = 0.0_DP
                 cg_count = 0
              end if
          else
             do is=1,pub_cell%num_spins
                call sparse_axpy(aux(is),con_direction(is),line_search_step)
             end do
          end if

          call kernel_workspace_invalidate()

          ! cks: ========================== END OF LINE SEARCH ================

       end if

       call internal_purify_and_print

    end do


    ! Report convergence failure
    if (pub_on_root .and. (.not. converged) .and. pub_output_detail >= NORMAL) &
         write(stdout,'(a,i4,a)') &
         'Finished density kernel iterations (',current_maxit_lnv, ')'


    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ** De-Allocate structures for sparse matrices

    ! Deallocate CG workspace
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(trial_aux(is))
       if(lnv_cg_type/='LNV_FLETCHER')then
          call sparse_destroy(old_co_gradient(is))
       endif
       call sparse_destroy(con_direction(is))
       call sparse_destroy(con_gradient(is))
       call sparse_destroy(co_gradient(is))
       call sparse_destroy(old_con_direction(is))
    end do

    deallocate(trial_aux,stat=ierr)
    call utils_dealloc_check('lnv_denskernel_optimise_cg', &
         'trial_aux',ierr)
    if(lnv_cg_type/='LNV_FLETCHER')then
       deallocate(old_co_gradient,stat=ierr)
       call utils_dealloc_check('lnv_denskernel_optimise_cg', &
            'old_co_gradient',ierr)
    endif
    deallocate(old_con_direction,stat=ierr)
    call utils_dealloc_check('lnv_denskernel_optimise_cg', &
         'old_con_direction',ierr)
    deallocate(con_direction,stat=ierr)
    call utils_dealloc_check('lnv_denskernel_optimise_cg', &
         'con_direction',ierr)
    deallocate(con_gradient,stat=ierr)
    call utils_dealloc_check('lnv_denskernel_optimise_cg','con_gradient', &
         ierr)
    deallocate(co_gradient,stat=ierr)
    call utils_dealloc_check('lnv_denskernel_optimise_cg','co_gradient', &
         ierr)

    ! Deallocate ks and ksk workspace
    call kernel_workspace_deallocate()

    ! Write density kernel to file if requested
    if (write_denskern.and.(.not.pub_write_converged_dk_ngwfs)) &
         call restart_kernel_write(aux,write_cond=pub_cond_calculate)

    ! Stop the timer
    call timer_clock('lnv_denskernel_optimise_cg', 2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving lnv_denskernel_optimise_cg'
#endif

    ! Flush output
    call services_flush

    return

  contains

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

!CW
    real(kind=DP) function internal_energy(kernel,update_ham,spoil_force)
!END CW
      use hamiltonian, only: hamiltonian_dens_dep_matrices

      implicit none

      ! Argument
      type(SPAM3), intent(inout) :: kernel(pub_cell%num_spins)
      logical, intent(in) :: update_ham

      ! Local variables
      integer :: is
      real(kind=DP) :: trks
      real(kind=DP) :: paw_sphere_energies(paw_en_size)
!CW
      logical,optional :: spoil_force
!END CW
      ! Calculate the purified density kernel
      call kernel_purify(purkern,kernel,rep%overlap,rep%inv_overlap,rep%n_occ)
      ! Stored versions of ks and ksk no longer valid (except on first call)
      if (.not.update_ham) call kernel_workspace_invalidate()

      ! Set normalisation factor for revised LNV
      if (exact_lnv) then
         do is=1,pub_cell%num_spins
            trks = sparse_trace(purkern(is),rep%overlap)
            if (rep%n_occ(is) > 0.0_DP) then
               norm_fac(is) = rep%n_occ(is) / trks
            else
               ! ndmh: prevent NaN's when one spin channel has nocc=0
               norm_fac(is) = 1.0_DP
            end if
         enddo
      else
         norm_fac(:) = 1.0_DP
      end if

      ! Convention within this routine is that the auxiliary matrix corresponds
      ! to just one spin state whereas the purified density kernel contains the
      ! spin degeneracy factor.
      do is=1,pub_cell%num_spins
         call sparse_scale(purkern(is),norm_fac(is)*spin_fac)
      end do

      ! ndmh: calculate density dependent energies and matrices
      if (.not.pub_cond_calculate) then
!CW
         call hamiltonian_dens_dep_matrices(ham, lhxc_fine, internal_energy, &
              lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
              ngwf_basis, hub_proj_basis, hub, purkern, ewald_energy, elements,&
              grid, localpseudo_fine, core_density_fine, &
              update_ham, lhxc_fixed=.false.,spoil_force=spoil_force)
!END CW
      else if (pub_cond_calculate) then
         ! lr408: In a conduction calculation, there is no need to recalculate
         ! lr408: any matrices, we just want the trace of the kernel with the
         ! lr408: projected Hamiltonian
         internal_energy = 0.0_dp
         do is=1,pub_cell%num_spins
            internal_energy = internal_energy + &
                 sparse_trace(ham%ham(is),purkern(is))
         end do

      end if

      ! Restore normalisation of purkern
      do is=1,pub_cell%num_spins
         call sparse_scale(purkern(is),1.0_DP/(norm_fac(is)*spin_fac))
      end do

    end function internal_energy

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_lnv_mu_old()

      !===================================================================!
      ! This function calculates mu, the chemical potential-like Lagrange !
      ! multiplier which preserves the electron number to first order in  !
      ! the original version of the LNV.  If there is no                  !
      ! truncation applied to the gradient, then this would be:           !
      !            mu = 6 * tr[LSLS(I-LS)(I-LS)]                          !
      ! which turns out to be proportional to the penalty functional.     !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      !===================================================================!

      implicit none

      ! Local variables
      integer :: is
      real(kind=DP) :: trace
      type(SPAM3) :: lslsc,lsc

      ! Create sparse matrix structures
      call sparse_create(lsc,pub_ks(1))
      call sparse_create(lslsc,pub_ks(1),lsc)

      call kernel_validate_ks(aux,rep%overlap)

      trace = 0.0_DP
      do is=1,pub_cell%num_spins

         ! Calculate (I - L.S)
         call sparse_copy(lsc,pub_ks(is))
         call sparse_scale(lsc,-1.0_DP,1.0_DP)

         ! Calculate L.S.(I - L.S)
         call sparse_product(lslsc,pub_ks(is),lsc)

         ! Calculate tr[L.S.L.S.(I-L.S).(I-L.S)]
         trace = trace + sparse_trace(lslsc,lslsc)

      end do

      ! Function result
      internal_lnv_mu_old = -3.0_DP * spin_fac * trace / ngwf_basis%num

      ! Deallocate workspace
      call sparse_destroy(lslsc)
      call sparse_destroy(lsc)

    end function internal_lnv_mu_old

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_lnv_mu_ideal()

      !===================================================================!
      ! This function calculates mu, the chemical potential-like Lagrange !
      ! multiplier which is used in the revised LNV and is simply given   !
      ! by:                                                               !
      !        mu = tr[LHLS(3I-2LS)] / tr[LSLS(3I-2LS)]                   !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, December 2004                            !
      ! Modified for spin polarisation by Peter Haynes, July 2006         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      implicit none

      ! Local variables
      integer :: is
      real(kind=DP) :: trh,trs
      type(SPAM3) :: lslsc,lh,lsc

      ! Create sparse matrix structures
      call sparse_create(lsc,pub_ks(1))
      call sparse_create(lh,aux(1),ham%ham(1))
      call sparse_create(lslsc,pub_ks(1),lsc)

      call kernel_validate_ks(aux,rep%overlap)

!      trs = 0.0_DP ; trh = 0.0_DP
      do is=1,pub_cell%num_spins

         ! Calculate products L.H, L.S and (3I - 2L.S)
         call sparse_product(lh,aux(is),ham%ham(is))
         call sparse_copy(lsc,pub_ks(is))
         call sparse_scale(lsc,-2.0_DP,3.0_DP)

         ! Calculate L.S.(3I - 2L.S)
         call sparse_product(lslsc,pub_ks(is),lsc)

         ! Calculate traces
         trs    = sparse_trace(lslsc,pub_ks(is))
         trh    = sparse_trace(lslsc,lh)
         if (abs(trs) > tiny(1.0_DP)) then
            mu(is) = trh / trs
         else
            mu(is) = 0.0_DP
         end if

      end do

      ! Deallocate workspace
      call sparse_destroy(lslsc)
      call sparse_destroy(lh)
      call sparse_destroy(lsc)

    end subroutine internal_lnv_mu_ideal

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_ms_mu_ideal()

      !===================================================================!
      ! This function calculates mu, the chemical potential-like Lagrange !
      ! multiplier which makes the contravariant gradient traceless in    !
      ! the Milliam-Scuseria variant of the LNV.  If there is no          !
      ! truncation applied to the gradient, then this would be:           !
      !            mu = 6 * tr[H(L-L.S.L)]/num which                      !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      ! Based on a version by Chris-Kriton Skylaris, 2001                 !
      ! Modified for spin polarisation by Peter Haynes, July 2006         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      implicit none

      ! Local variables
      integer :: is
      real(kind=DP) :: trace
      type(SPAM3) :: lslc

      ! Create sparse matrix structures
      call sparse_create(lslc,pub_ks(1),aux(1))

      call kernel_validate_ks(aux,rep%overlap)
      call kernel_validate_ksk(aux,rep%overlap)

      do is=1,pub_cell%num_spins

         ! Calculate L.S.L - L
         call sparse_copy(lslc,pub_ksk(is))
         call sparse_axpy(lslc,aux(is),-1.0_DP)

         ! Calculate tr[(L.S.L - L)H]
         trace = sparse_trace(lslc,ham%ham(is))
         mu(is)= -3.0_DP * spin_fac * trace / ngwf_basis%num

      end do

      ! Deallocate workspace
      call sparse_destroy(lslc)

    end subroutine internal_ms_mu_ideal

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_gradient_norm()

      implicit none

      ! Local variables
      integer :: is
      real(kind=DP) :: rms_el
      real(kind=DP) :: mu_extra(max_spins)    ! Additional mu resulting from
                                              ! truncation

      ! Start Timer
      call timer_clock('lnv_gradient_norm',1)

      ! Calculate covariant gradient (no tensor correction)
      call internal_co_gradient

      ! Calculate change in mu resulting from truncation for the
      ! covariant gradient
      if (.not. (old_lnv .or. exact_lnv)) then
         call internal_fix_co_direction(co_gradient,mu_extra)
      else
         mu_extra = 0.0_DP
      end if

      ! Apply tensor correction to obtain contravariant gradient
      call internal_con_gradient(mu_extra)

      ! Ensure contravariant gradient preserves electron number constraint
      ! even when matrix multiplication truncation is effective
      if (.not. (old_lnv .or. exact_lnv)) &
           call internal_fix_con_direction(con_gradient)

      ! Calculate inner product between gradients to obtain RMS gradient
      new_g_dot_g = 0.0_DP
      do is=1,pub_cell%num_spins
         new_g_dot_g = new_g_dot_g + &
              sparse_trace(con_gradient(is),co_gradient(is))
      end do

      if (new_g_dot_g < 0.0_DP .and. pub_on_root) write(stdout,'(a)') &
              'WARNING in internal_gradient_norm (lnv_mod.F90): &
              &inverse overlap is not positive definite!'

      internal_gradient_norm = sqrt(abs(new_g_dot_g) / &
           sparse_num_element(con_gradient(1)))

      ! Calculate RMS change in the density kernel from previous iteration
      ! and RMS commutator
      rms_dkern = 0.0_DP
      do is=1,pub_cell%num_spins
         rms_el = sparse_rms_element(old_con_direction(is))
         rms_dkern = rms_dkern + rms_el * rms_el
      end do
      rms_dkern = sqrt(rms_dkern / pub_cell%num_spins) * line_search_step

      commutator = kernel_rms_commutator(aux,ham%ham,rep%overlap)

      ! Stop Timer
      call timer_clock('lnv_gradient_norm',2)

    end function internal_gradient_norm

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_co_gradient

      !======================================================================!
      ! This subroutine calculates the gradient of the LNV functional with   !
      ! respect to the auxiliary matrix.                                     !
      !----------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                               !
      ! Based on an original version by Chris-Kriton Skylaris modified by    !
      ! Arash Mostofi and Peter Haynes                                       !
      ! Spin polarised by Peter Haynes, July 2006                            !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009     !
      ! Improvements to reduce number of matrix products required, by        !
      ! Nicholas Hine in November 2010.                                      !
      !======================================================================!

      implicit none

      ! Local variables
      integer :: is
      type(SPAM3) :: sl, slh, hls, lsc

      ! Note: con_gradient is used as workspace in this routine since it is
      !       not used until this routine has completed

      ! Create sparse matrix structures
      call sparse_create(sl,rep%overlap,aux(1))
      call sparse_create(lsc,aux(1),rep%overlap)
      call sparse_create(slh,sl,ham%ham(1))
      call sparse_create(hls,ham%ham(1),pub_ks(1))

      call kernel_validate_ks(aux,rep%overlap)

      do is=1,pub_cell%num_spins

         ! In revised LNV, replace H by H - mu * S
         if (exact_lnv) call sparse_axpy(ham%ham(is),rep%overlap,-mu(is))

         ! Calculate product S.L
         call sparse_transpose(sl,pub_ks(is))

         ! Calculate products S.L.H and H.L.S
         call sparse_product(slh,sl,ham%ham(is))
         call sparse_transpose(hls,slh)

         ! Calculate H.L.S.(3I - 2L.S) in con_gradient
         call sparse_copy(lsc,pub_ks(is))
         call sparse_scale(lsc,-2.0_DP,3.0_DP)
         call sparse_product(con_gradient(is),hls,lsc)

         ! Calculate (3I - 2S.L).S.L.H in co_gradient
         call sparse_transpose(co_gradient(is),con_gradient(is))

         ! Add up H.L.S.(3I - 2L.S) and (3I - 2S.L).S.L.H in co_gradient
         call sparse_axpy(co_gradient(is),con_gradient(is),1.0_DP)

         ! Calculate product S.L.H.L.S in con_gradient and add
         ! -2 S.L.H.L.S to co_gradient
         call sparse_product(con_gradient(is),slh,pub_ks(is))
         call sparse_axpy(co_gradient(is),con_gradient(is),-2.0_DP)

         ! Add -mu*S (For revised LNV, mu = 0 here).
         ! Total grad: (3I-2S.L).S.L.H + H.L.S.(3I-2L.S) - 2S.L.H.L.S - mu*S
         if (.not. exact_lnv) &
              call sparse_axpy(co_gradient(is),rep%overlap,-mu(is))

         ! Multiply by spin and normalisation factor
         call sparse_scale(co_gradient(is),spin_fac*norm_fac(is))

         ! In revised LNV, return H - mu * S to just H
         if (exact_lnv) call sparse_axpy(ham%ham(is),rep%overlap,mu(is))

      end do

      ! Deallocate workspace
      call sparse_destroy(hls)
      call sparse_destroy(slh)
      call sparse_destroy(lsc)
      call sparse_destroy(sl)

    end subroutine internal_co_gradient

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_con_gradient(mu_extra)

      !===================================================================!
      ! This subroutine converts a covariant gradient to a contravariant  !
      ! gradient by pre- and post-multiplying by the inverse overlap      !
      ! matrix. Note that truncation will introduce errors at this stage. !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      ! Based on an original version by Chris-Kriton Skylaris and         !
      ! modified by Peter Haynes                                          !
      ! Spin polarised by Peter Haynes, July 2006                         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      implicit none

      ! Argument
      integer :: is
      real(kind=DP), intent(in) :: mu_extra(max_spins)  ! additional mu from
                                                        ! truncation

      ! Local variables
      type(SPAM3) :: temp

      ! Allocate temporary sparse matrix
      call sparse_create(temp,co_gradient(1),rep%inv_overlap)

      do is=1,pub_cell%num_spins

         ! Calculate con_grad := co_grad + mu_extra S
         call sparse_copy(con_gradient(is),co_gradient(is))
         call sparse_axpy(con_gradient(is),rep%overlap,mu_extra(is))

         ! Calculate temp := con_grad . S^-1
         call sparse_product(temp,con_gradient(is),rep%inv_overlap)

         ! Calculate S^-1 . (co_grad + mu_extra S) . S^-1 =
         !                        = S^-1 . temp -> con_grad
         call sparse_product(con_gradient(is),rep%inv_overlap,temp)

      end do

      ! Deallocate workspace
      call sparse_destroy(temp)

    end subroutine internal_con_gradient

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_fix_co_direction(co_dir,mu_extra)

      !===================================================================!
      ! This subroutine ensures that a covariant search direction         !
      ! preserves the electron number constraint (to first order).        !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      ! Spin polarised by Peter Haynes, July 2006                         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      implicit none

      ! Argument
      type(SPAM3), intent(inout) :: co_dir(pub_cell%num_spins)  ! co direction
      real(kind=DP), intent(out) :: mu_extra(max_spins)  ! additional mu required

      ! Local variables
      integer :: is
      real(kind=DP) :: trdinvs     ! tr[D.S^-1]
      real(kind=DP) :: trsinvs     ! tr[S.S^-1]

      ! Calculate tr[D.S^-1] and tr[S.S^-1]
      trsinvs = sparse_trace(rep%overlap,rep%inv_overlap)

      trdinvs = 0.0_DP
      do is=1,pub_cell%num_spins

         trdinvs = sparse_trace(co_dir(is),rep%inv_overlap)

         ! Calculate and report correction to mu
         mu_extra(is) = -trdinvs / trsinvs

         if (pub_on_root .and. (pub_output_detail == VERBOSE)) &
              write(stdout,'(a,i1,a,f20.12)') '     fix cov mu', &
              is,': ', mu_extra(is)

      end do


    end subroutine internal_fix_co_direction

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_fix_con_direction(con_dir)

      !===================================================================!
      ! This subroutine ensures that a contravariant search direction     !
      ! preserves the electron number constraint (to first order).        !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      ! Spin polarised by Peter Haynes, July 2006                         !
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!

      implicit none

      ! Argument
      type(SPAM3), intent(inout) :: con_dir(pub_cell%num_spins)

      ! Local variables
      integer :: is
      real(kind=DP) :: trds                ! tr[D.S]
      real(kind=DP) :: trsinvs             ! tr[S^-1.S]
      real(kind=DP) :: mu_extra(max_spins) ! additional mu required

      ! Calculate tr[D.S] and tr[S.S^-1]
      trds = 0.0_DP
      trsinvs = sparse_trace(rep%inv_overlap,rep%overlap)

      do is=1,pub_cell%num_spins

         trds = sparse_trace(con_dir(is),rep%overlap)

         ! Calculate and report correction to mu
         mu_extra(is) = -trds / trsinvs

         if (pub_on_root .and. (pub_output_detail == VERBOSE)) &
              write(stdout,'(a,i1,a,f20.12)') '     fix con mu', &
              is,': ', mu_extra(is)

      end do

      ! Correct gradient
      do is=1,pub_cell%num_spins
         call sparse_axpy(con_dir(is),rep%inv_overlap,mu_extra(is))
      end do

    end subroutine internal_fix_con_direction

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    logical function internal_test_convergence()

      implicit none

      ! Convergence is assumed when either the RMS gradient or RMS commutator
      ! falls below the threshold

      if ((rms_gradient < lnv_threshold .and. iteration > 0) .or. &
           (commutator < lnv_threshold .and. iteration > 1)) then

         if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then

            write(stdout,'(i2,f22.14,2f14.11)') iteration, total_energy, &
                 rms_dkern, commutator
            write(stdout,'(/a)')'                 ......................&
                 &........................'
            write(stdout,'(a,f19.14,a)')'                 | RMS LNV GRADIENT&
                 &= ', rms_gradient, '      |'
            write(stdout,'(a)') '                 | LNV density kernel &
                 &optimisation converged! |'
            write(stdout,'(a)') '                 ~~~~~~~~~~~~~~~~~~~~~~&
                 &~~~~~~~~~~~~~~~~~~~~~~~~'

         end if

         internal_test_convergence = .true.

      else

         internal_test_convergence = .false.

      end if

    end function internal_test_convergence

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_mu_term(trial_aux)
      !===================================================================!
      ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009  !
      !===================================================================!
      implicit none

      ! Argument
      type(SPAM3), intent(in) :: trial_aux(pub_cell%num_spins)

      ! cks: local variable
      real(kind =DP) :: trks

      internal_mu_term = 0.0_DP
      if (exact_lnv) return

      ! Calculate 2 tr(L.S)
      do is=1,pub_cell%num_spins

         trks = sparse_trace(trial_aux(is),rep%overlap)
         trks = trks*spin_fac

         ! Calculate [2 tr(L.S) - Ne]
         trks = trks - spin_fac*real(rep%n_occ(is),kind=DP)

         if ((.not. old_lnv) .and. abs(trks) > 1.0e-6_DP .and. &
              pub_on_root) write(stdout,'(a,f20.12)') &
              'WARNING: electron number wrong by ', trks

         internal_mu_term = internal_mu_term - mu(is)*trks

      end do


      if (pub_on_root .and. (pub_output_detail == VERBOSE)) &
           write(stdout,'(a,f20.12)') '-mu(2tr(KS)-Ne): ', internal_mu_term

    end function internal_mu_term

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_search_direction

      implicit none

      ! Local variables
      integer :: is
      real(kind=DP) :: polak
      real(kind=DP) :: dltne

      ! Check previous scalar product of gradients is not vanishing

      if (abs(old_g_dot_g) < tiny(1.0_DP)) then

         if (pub_on_root) write(stdout,'(a/a)') &
              'WARNING in lnv_denskernel_optimise_cg: &
              &vanishing previous gradient norm -', &
              '  conjugate gradients coefficient reset'

         cg_coeff = 0.0_DP

      else

         ! Calculate correction factor for Polak-Ribiere formula if necessary
         if (lnv_cg_type == 'LNV_FLETCHER') then
            polak = 0.0_DP
         else if (lnv_cg_type == 'LNV_POLAK') then
            polak = 0.0_DP
            do is=1,pub_cell%num_spins
               polak = polak + &
                    sparse_trace(con_gradient(is),old_co_gradient(is))
            end do
         else
            if (pub_on_root) write(stdout,'(a/3a)') &
              'WARNING in lnv_denskernel_optimise_cg: &
              &unknown conjugate gradients method','  "',trim(lnv_cg_type), &
              '" - switching to steepest descents'
            polak = new_g_dot_g
         end if

         ! Calculate conjugate gradients coefficient
         cg_coeff = (new_g_dot_g - polak) / old_g_dot_g

         ! Reset CG coefficient if too large
         if (abs(cg_coeff) > 5.0_DP) then
            if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
                 write(stdout,'(a/a,e10.3,a)') &
                 'WARNING in lnv_denskernel_optimise_cg: ', &
                 '  conjugate gradient coefficient (', &
                 cg_coeff, ') too large - resetting'
            cg_coeff = 0.0_DP
            cg_count = 0
         end if

      end if

      ! Increment CG counter or reset as appropriate
      if (abs(cg_coeff) > eps) then
         cg_count = cg_count + 1
      else
         cg_count = 0
      end if

      ! Reset conjugate gradients after cg_max steps
      if (cg_count > cg_max) then
         cg_coeff = 0.0_DP
         cg_count = 0
      end if

      ! Calculate the line search direction
      do is=1,pub_cell%num_spins
         call sparse_copy(con_direction(is),con_gradient(is))
         ! pdh: hack to keep step lengths the same for spin polarised
         ! pdh: systems: multiply by a factor of two
         call sparse_scale(con_direction(is),-real(pub_cell%num_spins,kind=DP))
         call sparse_axpy(con_direction(is),old_con_direction(is),cg_coeff)
      end do

      ! Ensure that search direction preserves constraint to first order
      ! pdh: (not necessary for current algorithms)
      !      if (.not. (old_lnv .or. exact_lnv)) &
      !           call internal_fix_con_direction(con_direction)
      if ( (pub_output_detail == VERBOSE) .and.(.not. exact_lnv)) then
         dltne = 0.0_DP
         do is=1,pub_cell%num_spins
            dltne = dltne + sparse_trace(con_direction(is),rep%overlap)
         end do
         dltne = dltne * spin_fac
         if (pub_on_root) write(stdout,'(a,f20.12)') 'DltNe tr(Dir*S): ',dltne
      end if

      ! Calculate the function derivative along the search direction
      Qslope = 0.0_DP
      do is=1,pub_cell%num_spins
         Qslope = Qslope + sparse_trace(co_gradient(is),con_direction(is))
      end do

      ! Ensure search direction points 'downhill'
      if (Qslope > 0.0_DP) then

         ! First try steepest descent direction i.e. reset conjugate gradients
         if ( (pub_output_detail >= NORMAL) .and. pub_on_root) &
              write(stdout,'(a/a)') 'WARNING in lnv_denskernel_optimise_cg: &
              &positive line search gradient -', &
              '  resetting conjugate gradients'

         Qslope = 0.0_DP
         do is=1,pub_cell%num_spins
            call sparse_copy(con_direction(is),con_gradient(is))
            call sparse_scale(con_direction(is),-1.0_DP)
            Qslope = Qslope + sparse_trace(co_gradient(is),con_direction(is))
         end do

         ! If this is still positive then things are really pear-shaped i.e.
         ! overlap matrix isn't positive definite or something...
         ! Reverse search direction!
        if (Qslope > 0.0_DP) then
           if ( (pub_output_detail >= NORMAL ) .and. pub_on_root) &
                write(stdout,'(a)') 'WARNING in lnv_denskernel_optimise_cg: &
                &reversing search gradient'
           do is=1,pub_cell%num_spins
              call sparse_scale(con_direction(is),-1.0_DP)
           end do
           Qslope = -Qslope
        end if

     end if

      if ((pub_output_detail == VERBOSE) .and. pub_on_root) &
           write(stdout,'(a,f20.12)') '             G0: ',Qslope

       ! Store quantities for next iteration
      do is=1,pub_cell%num_spins
         call sparse_copy(old_con_direction(is),con_direction(is))
      if (lnv_cg_type == 'LNV_POLAK') &
           call sparse_copy(old_co_gradient(is),co_gradient(is))
      end do
      old_g_dot_g = new_g_dot_g

    end subroutine internal_search_direction

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_kernel_protected_step(trial_kernel,trial_step,kernel, &
         con_dirn)

      !===================================================================!
      ! This subroutine ensures that a chosen line step does not cause    !
      ! kernel occupation numbers to become unstable, by checking the     !
      ! occupancy bounds and RMS occupancy error and reducing the line    !
      ! step until they are within the LNV stability range.               !
      !-------------------------------------------------------------------!
      ! Written by Nicholas Hine, November 2009.                          !
      !===================================================================!

      implicit none

      ! Arguments
      type(SPAM3),intent(inout) :: trial_kernel(pub_cell%num_spins)
      real(kind=DP), intent(inout) :: trial_step
      type(SPAM3),intent(in) :: kernel(pub_cell%num_spins)
      type(SPAM3),intent(in) :: con_dirn(pub_cell%num_spins)

      ! Local Variables
      real(kind=DP), parameter :: minbound = -0.36602540378443864676_DP
      real(kind=DP), parameter :: maxbound = 1.36602540378443864676_DP
      real(kind=DP) :: rms_err, trial_rms_err
      real(kind=DP) :: min_occ(max_spins), max_occ(max_spins)
      real(kind=DP) :: trial_min_occ(max_spins), trial_max_occ(max_spins)
      integer, parameter :: max_trials=8
      integer :: step

      ! ndmh: find rms occupancy error of starting kernel
      rms_err = kernel_rms_err(kernel,rep%overlap)

      ! ndmh: find occupancy bounds of starting kernel
      call kernel_occupancy_bounds(max_occ,min_occ,kernel,rep%overlap)

      ! ndmh: starting at input value of trial_step, keep reducing trial_step
      ! ndmh: until an acceptable value is found
      do step=1,max_trials

         ! Take trial step along the search direction
         do is=1,pub_cell%num_spins
            call sparse_copy(trial_kernel(is),kernel(is))
            call sparse_axpy(trial_kernel(is),con_dirn(is),trial_step)
         end do
         call kernel_workspace_invalidate()

         ! ndmh: find rms occupancy error of trial kernel
         trial_rms_err = kernel_rms_err(trial_kernel,rep%overlap)

         ! ndmh: find occupancy bounds of trial kernel
         call kernel_occupancy_bounds(trial_max_occ,trial_min_occ, &
              trial_kernel,rep%overlap)

         ! ndmh: reduce step length by half if rms err increased by >20% or if
         ! ndmh: the new kernel has occupancies outside stable LNV bounds
         if ((trial_rms_err > rms_err+0.2_DP).or. &
              (any(trial_min_occ(:)<minbound) .or. &
              any(trial_max_occ(:)>maxbound))) then
            if (pub_on_root .and. (pub_output_detail==VERBOSE)) then
               if (step==1) then
                  write(stdout,'(a)') '========================= Kernel &
                       &Line Step Protection =========================='
                  write(stdout,'(a)') '|  Iter |  Trial Step |   RMS &
                       &Occupancy Error |   Occupancy Bounds | Accepted? |'
                  write(stdout,'(i5,f16.6,f24.6,4x,"[",f7.4,":",f8.4,"]")') &
                       0,0.0_DP,rms_err,maxval(max_occ),minval(min_occ)
               end if
               write(stdout,'(i5,f16.6,f24.6,4x,"[",f7.4,":",f8.4,"]",a)') &
                    step,trial_step,trial_rms_err,maxval(trial_max_occ), &
                    minval(trial_min_occ),'        No'
            end if
            trial_step = trial_step * 0.50_DP
         else
            ! ndmh: this trial step was OK, so exit
            if ((step > 1) .and. (pub_output_detail==VERBOSE) .and. &
                 pub_on_root) then
               write(stdout,'(i5,f16.6,f24.6,4x,"[",f7.4,":",f8.4,"]",a)') &
                    step,trial_step,trial_rms_err,maxval(trial_max_occ), &
                    minval(trial_min_occ),'       Yes'
               write(stdout,'(a)') '================================&
                    &================================================'
            end if
            return
         end if

      end do

      ! ndmh: finished loop without finding a stable kernel, so warn about this
      if (pub_on_root) then
         write(stdout,'(a,i2)') 'WARNING in lnv_denskernel_optimise_cg: Trial &
              &step reduced ',max_trials
         write(stdout,'(a)') 'times but kernel may still be unstable'
      end if

    end subroutine internal_kernel_protected_step

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_quadratic_step()

      !===================================================================!
      ! This function fits the quadratic form                             !
      !       2                                                           !
      !    a x  + b x + c = Q(x)                                          !
      !                                                                   !
      ! given Q(0) = Qinitial, Q(trial_step) = Qtrial and Q'(0) = Qslope  !
      !                                                                   !
      ! and then finds the minimum position.                              !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      !===================================================================!

      implicit none

      ! Local variables
      real(kind=DP) :: a,b            ! Coefficients of quadratic equation
      !real(kind=DP) :: qdiff          ! Qtrial - Qinitial qoh: Unused
      real(kind=DP) :: linear_term    ! b * trial_step
      real(kind=DP) :: quadratic_term ! a * trial_step**2

      ! Determine coefficients of quadratic equation

      linear_term = Qslope * trial_step
      quadratic_term = Qtrial - Qinitial - linear_term
      a = quadratic_term / (trial_step * trial_step)
      b = Qslope

      if ((pub_output_detail == VERBOSE) .and. pub_on_root) then
         write(stdout,'(a,f20.12)') ' Quadratic coef: ',quadratic_term
         write(stdout,'(a,f20.12)') ' Linear coef   : ',linear_term
      end if


      ! Find minimum position x = -b / (2 a)
      if (abs(a) > abs(Qinitial)*eps*eps) then

         internal_quadratic_step = -0.5_DP * b / a

         Qpredict = Qinitial + 0.5_DP * b * internal_quadratic_step ! c-b^2/(4a)

      else

         internal_quadratic_step = -1.0_DP
         if ((pub_output_detail == VERBOSE) .and. pub_on_root) &
              write(stdout,'(a)') '  Quadratic term too small: &
              &proceeding to cubic fit'

      end if

    end function internal_quadratic_step

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    real(kind=DP) function internal_cubic_step()

      !===================================================================!
      ! This function fits the cubic form                                 !
      !       3      2                                                    !
      !    a x  + b x  + c x + d = Q(x)                                   !
      !                                                                   !
      ! given Q(0) = Qinitial, Q(trial_step) = Qtrial, Q'(0) = Qslope     !
      ! and now Q(optimal_step) = Qoptimal                                !
      !                                                                   !
      ! and then finds the minimum position.                              !
      !-------------------------------------------------------------------!
      ! Written by Peter Haynes, November 2004                            !
      !===================================================================!

      implicit none

      ! Local variables
      real(kind=DP) :: x1,x2   ! The trial step lengths
      real(kind=DP) :: a,b,c   ! Coefficients of the cubic equation
      real(kind=DP) :: bon3a   ! b / (3 a)
      real(kind=DP) :: con3a   ! c / (3 a)
      real(kind=DP) :: disc    ! Discriminant
      real(kind=DP) :: q

      ! Introduce shorthand working variables:

      !   x1 and x2 are the points at which the energy has been evaluated

      x1 = trial_step
      x2 = optimal_step

      !                                          2
      ! The derivative of the cubic fit is  3 a x  + b x + c

      a = ((x2*x2*Qtrial - x1*x1*Qoptimal)/(x1-x2) &
           + (x1+x2)*Qinitial + x1*x2*Qslope) / (x1*x1*x2*x2)
      b = ((x2*x2*x2*Qtrial - x1*x1*x1*Qoptimal)/(x2-x1) &
           - (x1*x1+x1*x2+x2*x2)*Qinitial - x1*x2*(x1+x2)*Qslope) &
           / (x1*x1*x2*x2)
      c = Qslope

      ! Solve this quadratic equation

      if (abs(a*x2/b) > eps) then                      ! Avoid division by zero

         bon3a = b / (3.0_DP * a)
         con3a = c / (3.0_DP * a)
         disc = bon3a * bon3a - con3a                  ! Discriminant

         if (disc >= 0.0_DP) then

            q = -(bon3a + sign(sqrt(disc),bon3a))
            if (b < 0.0_DP) then
               internal_cubic_step = q
            else
               internal_cubic_step = con3a / q
            end if

            if (pub_on_root .and. (pub_output_detail == VERBOSE)) &
                 write(stdout,'(a,f20.12)') '  Cubic LS coeff: ', &
                 internal_cubic_step

            Qminimum = Qinitial + internal_cubic_step * &
                 (c + internal_cubic_step * (b + internal_cubic_step * a))

            if (pub_on_root .and. (pub_output_detail == VERBOSE)) &
                 write(stdout,'(a,f20.12)') 'Predicted energy: ',Qminimum

            ! Do not allow the step length to become too long

            if ( abs(internal_cubic_step) > lnv_cg_max_step/2.0_DP) then
               if (pub_on_root .and. (pub_output_detail>=NORMAL)) &
                    write(stdout,'(a)') &
                    'WARNING in lnv_denskernel_optimise_cg: &
                    &setting cubic step to safe value'
               internal_cubic_step = 0.15_DP
            else
               if (pub_on_root .and. (pub_output_detail == VERBOSE)) &
                    write(stdout,'(a)') 'Cubic minimum step length accepted'
            end if

         else

            if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
                 write(stdout,'(a/a)') 'WARNING in lnv_denskernel_optimise_cg: &
                 &line search unsuccessful -','  optimal step length set to &
                 &safe value'

            internal_cubic_step = 0.15_DP

         end if

      else

         if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
              write(stdout,'(a/a)') 'WARNING in lnv_denskernel_optimise_cg: &
              &line search unsuccessful -', '  optimal step length set to &
              &safe value'

         internal_cubic_step = 0.15_DP

      end if

    end function internal_cubic_step

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////

    subroutine internal_purify_and_print

      implicit none

      ! Local variables
      character(55) :: fmt
      integer :: is
      integer :: iter,numiter


      if (iteration < current_maxit_lnv) then

         if (exact_lnv) then

               ! Perform some "adaptive" purifications
               if (rms_gradient > 0.00001_DP) then
                  numiter = 3
               else
                  numiter = 1
               end if

               ! cks: First fix Ne of L
               call kernel_rescale(aux,rep%overlap,rep%n_occ, &
                    can_rescale_ks=.true.)
               ! cks: Now fix the idempotency of L
               call kernel_fix(aux, rep%overlap, rep%inv_overlap, numiter)

               ! cks: to ensure L is purified adequately even when the
               ! cks: kernel_fix has exited without doing any iterations
! cks, 16/04/2005: Unfortunately this hinders convergence with denskern truncation.
!               call kernel_purify(aux,aux,overlap,inv_overlap, n_occ)

         else

            if (old_lnv) then

               ! Perform some ad hoc purifications
               if (rms_gradient > 0.00001_DP) then
                  numiter = 3
               else
                  numiter = 1
               end if
               do iter=1,numiter
                  call kernel_purify(aux,aux,rep%overlap,rep%inv_overlap, &
                       rep%n_occ)
                  call kernel_workspace_invalidate()
               end do

            else

               ! Perform some "adaptive" purifications
               if (rms_gradient > 0.00001_DP) then
                  numiter = 3
               else
                  numiter = 1
               end if

               call kernel_fix(aux, rep%overlap, rep%inv_overlap, numiter)

               ! cks: to ensure the kernel is purified adequately even when
               ! cks: the kernel_fix has exited without doing any iterations
               call kernel_purify(aux,aux,rep%overlap,rep%inv_overlap,rep%n_occ)
               ! aux has changed => invalidate pub_ks and pub_ksk
               call kernel_workspace_invalidate()

               ! cks: The "unpurified" denskern must have the correct Ne.
               ! cks: Enforce this in the face of truncation.
               call kernel_rescale(aux,rep%overlap,rep%n_occ)

            end if

         end if

         ! Calculate purified density kernel
         call kernel_purify(purkern,aux,rep%overlap,rep%inv_overlap,rep%n_occ)
         ! pub_ks and pub_ksk are still valid for aux after this

         ! Calculate electron number from purified kernel
         if (exact_lnv) then
            ne_pure = sum(rep%n_occ)
         else
            ne_pure = 0.0_DP
            do is=1,pub_cell%num_spins
               ne_pure = ne_pure + sparse_trace(purkern(is),rep%overlap)
            end do
            ne_pure = ne_pure * spin_fac
         end if

      else

         ! Last iteration
         ne_pure = 0.0_DP
         do is=1,pub_cell%num_spins
            ne_pure = ne_pure + sparse_trace(aux(is),rep%overlap)
         end do
         ne_pure = ne_pure * spin_fac

         cg_coeff = 0.0_DP
         line_search_step = 0.0_DP

      end if

      if (pub_on_root .and. (pub_output_detail >= NORMAL) ) then

         ! ndmh: adapt format string so that it always displays sensible
         ! ndmh: numbers of digits for total energy and Ne
         write(fmt,'(a)')'(i2,'
         if(abs(total_energy)<100000_DP) then
            write(fmt,'(a,a)') trim(fmt),'f22.14'
         else
            write(fmt,'(a,a)') trim(fmt),'f22.12'
         end if
         write(fmt,'(a,a)') trim(fmt),',2f14.11,f7.4,f7.3,'
         if(ne_pure<10000_DP) then
            write(fmt,'(a,a)') trim(fmt),'f14.8)'
         else if (ne_pure<100000_DP) then
            write(fmt,'(a,a)') trim(fmt),'f14.7)'
         else
            write(fmt,'(a,a)') trim(fmt),'f14.5)'
         end if

         write(stdout,fmt) iteration,total_energy,rms_dkern,commutator,&
              line_search_step,cg_coeff,ne_pure
      end if

    end subroutine internal_purify_and_print

    !////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////



  end subroutine lnv_denskernel_optimise_cg

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

end module lnv
