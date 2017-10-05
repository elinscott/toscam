! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                     Penalty functional module                  !
!                                                                !
! This module implements an alternative to the LNV method for    !
! optimising the density kernel K for a fixed set of NGWFs,      !
! observing the idempotency constraint approximately by means of !
! a penalty functional related to McWeeny's purifying            !
! transformation.                                                !
!----------------------------------------------------------------!
! Written by Peter Haynes, 14/02/06                              !
! Based on an earlier version by Peter Haynes and Chris-Kriton   !
!   Skylaris.                                                    !
! Adapted to use SPAM3 matrices by Nicholas Hine, July 2009      !
!================================================================!

module penalty

  implicit none

  private

  public :: penalty_denskernel_optimise_cg

contains

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine penalty_denskernel_optimise_cg(denskern, rep, ham, &
       lhxc_fine, total_energy, localpseudo_fine, core_density_fine, &
       ngwf_basis, hub_proj_basis, hub, ewald_energy, elements, pen_threshold)

    !==========================================================================
    ! This subroutine optimises the density kernel for a fixed set of NGWFs
    ! using the penalty functional to impose idempotency approximately.
    !
    ! See P.D. Haynes and M.C. Payne, Phys. Rev. B 59, 12173 (1999).
    !
    !==========================================================================
    ! Arguments:
    ! denskern (inout)        : Density kernel in row-indexed sparse storage
    ! rep (input)             : NGWF representation (functions and matrices)
    ! lhxc_fine (output)      : Local-Hartree-Exhange-Correlation potential
    ! total_energy (output)   : Total energy
    ! localpseudo_fine (input): Local pseudopotential on fine grid
    ! core_density_fine (input): Core density on fine grid
    ! ngwf_basis (input)      : Function basis of NGWFs
    ! ewald_energy (input)    : Ewald energy
    ! pen_threshold (input)   : Convergence threshold for LNV gradient
    !==========================================================================
    ! Based on lnv_mod.F90 written by Chris-Kriton Skylaris.
    ! Original version written by Peter Haynes.
    ! Parallelised by Peter Haynes and Chris-Kriton Skylaris in December 2003.
    ! Rewritten to use new sparse matrix format by Peter Haynes, February 2006.
    ! Modified to use parallel SPAM 2 throughout by Peter Haynes, July 2006
    ! Modified to include Nonlinear Core Corrections by Nicholas Hine in
    ! January 2009
    ! DFT+U added by David O'Regan, April 2009
    ! Adapted to use SPAM3 matrices by Nicholas Hine, July 2009
    !==========================================================================

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root, pub_root_node_id, comms_bcast
    use constants, only: DP, VERBOSE, stdout, NORMAL, max_spins, paw_en_size, &
         paw_en_ehart, paw_en_exc, paw_en_exc_dc, paw_en_etxc, paw_en_etxc_dc, &
         paw_en_dijhat, paw_en_dijxc, paw_en_exc_core
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix
    use hubbard_build, only: HUBBARD_MODEL, hubbard_energy_total
    use integrals, only: integrals_locpot
    use ion, only: ELEMENT
    use kernel, only: kernel_rms_commutator, pub_ks, kernel_fix, &
         kernel_workspace_invalidate, kernel_workspace_allocate, &
         kernel_validate_ks, kernel_workspace_deallocate, kernel_normalise, &
         kernel_occupancy_bounds
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_ham_create, &
         ngwf_ham_destroy
    use restart, only : restart_kernel_write
    use rundat, only: maxit_pen, pen_param, pub_any_nl_proj, pub_hubbard, &
         pub_output_detail, write_denskern, pub_write_converged_dk_ngwfs, &
         pub_usehfx, pub_aug
    use services, only: services_flush
    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_trace, sparse_num_element, &
         sparse_create, sparse_product, &
         sparse_axpy, sparse_copy, sparse_scale, sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use vdwcorrection, only: pub_dispersion_energy

    implicit none

    ! Arguments: input/output
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
    type(NGWF_HAM), intent(inout) :: ham
    real(kind=DP), intent(out) :: lhxc_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(out) :: total_energy
    type(HUBBARD_MODEL), intent(inout) :: hub

    ! Arguments: input only
    real(kind=DP), intent(in) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(NGWF_REP), intent(in) :: rep
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=DP), intent(in) :: ewald_energy, pen_threshold

    ! Local variables
    real(kind=DP) :: commutator
    real(kind=DP) :: penalty_energy
    real(kind=DP) :: rms_gradient,line_search_coeff,cg_coeff,eps
    real(kind=DP) :: total_energy_at_trial_step,lhxc_energy
    real(kind=DP) :: max_occ(max_spins),min_occ(max_spins)
    real(kind=DP) :: spin_fac
    integer       :: is
    integer       :: ierr
    integer       :: iteration,cg_count,cg_max
    logical       :: converged

    ! Conjugate gradient variables
    type(SPAM3), allocatable :: co_gradient(:)       ! Covariant gradient
    type(SPAM3), allocatable :: con_gradient(:)      ! Contravariant gradient
    type(SPAM3), allocatable :: con_direction(:)     ! Contra search direction
    type(SPAM3), allocatable :: old_con_direction(:) ! Prev contra search dirn
    type(SPAM3), allocatable :: old_co_gradient(:)   ! Prev covariant gradient
    type(SPAM3), allocatable :: trial_denskern(:)    ! Trial density kernel

    !--------------------------------------------------------------------------
    ! pdh: variables introduced purely for line minimisation
    ! pdh: the function to be minimised is generally called 'Q'
    !--------------------------------------------------------------------------
    real(kind=DP) :: Qinitial          ! initial value of function
    real(kind=DP) :: Qslope            ! derivative along search direction
    real(kind=DP) :: trial_step=1.0_DP ! trial step length
    real(kind=DP) :: Qtrial            ! value of function at trial step
    real(kind=DP) :: optimal_step      ! optimal step length
    real(kind=DP) :: Qoptimal          ! value of function at optimal step
    real(kind=DP) :: Qpredict          ! predicted value at optimal step
    real(kind=DP) :: minimum_step      ! minimum step predicted by cubic
    real(kind=DP) :: Qminimum          ! predicted value of function at minimum
    real(kind=DP) :: linear_term
    real(kind=DP) :: quadratic_term
    real(kind=DP) :: xx, yy, aa, bb, aon3b, disc
    logical       :: trial2
    !--------------------------------------------------------------------------
    ! pdh: end variables for line minimisation
    !--------------------------------------------------------------------------

    ! cks : CONVENTIONS:
    ! cks : search directions are ALWAYS CONTRAVARIANT.
    ! cks : gradients can be either covariant or contravariant.
    ! cks : their name will indicate this in each case.

    ! pdh : All density kernels are for spin alpha only, i.e. all spin
    !       degeneracy factors are explicitly included in this module

    ! cks: find time spent in this subroutine
    call timer_clock('penalty_denskernel_optimise_cg',1)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &penalty_denskernel_optimise_cg'
#endif

    call services_flush

    ! ndmh: allocate kernel work arrays
    call kernel_workspace_allocate(denskern,rep%overlap)

    ! Allocate matrices with the sparsity pattern of the density kernel
    allocate(co_gradient(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('penalty_denskernel_optimise_cg','co_gradient',ierr)
    allocate(con_gradient(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('penalty_denskernel_optimise_cg','con_gradient',ierr)
    allocate(con_direction(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('penalty_denskernel_optimise_cg','con_direction', &
         ierr)
    allocate(old_con_direction(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('penalty_denskernel_optimise_cg', &
         'old_con_direction',ierr)
    allocate(old_co_gradient(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('penalty_denskernel_optimise_cg', &
         'old_co_gradient',ierr)

    allocate(trial_denskern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('penalty_denskernel_optimise_cg', &
         'trial_denskern',ierr)

    do is=1,pub_cell%num_spins
       call sparse_create(co_gradient(is),denskern(is))
       call sparse_create(con_gradient(is),denskern(is))
       call sparse_create(con_direction(is),denskern(is))
       call sparse_create(old_con_direction(is),denskern(is))
       call sparse_create(old_co_gradient(is),denskern(is))
       call sparse_create(trial_denskern(is),denskern(is))
    end do


    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ! Make sure the density-matrix corresponds to exactly the right
    ! number of electrons
    spin_fac = 2.0_DP / pub_cell%num_spins
    call kernel_normalise(denskern, rep%overlap, rep%inv_overlap, rep%n_occ)

    ! Initialise conjugate gradients
    rms_gradient = 0.0_DP
    cg_count = 0
    cg_max = 5
    eps = epsilon(1.0_DP)
    converged = .false.

    if (pub_on_root) then
       write(stdout,'(a)') '.........................................&
            &.......................................'
       write(stdout,'(a)') '<<<<<<<<<<<<<<<< Penalty functional &
            &density kernel optimisation >>>>>>>>>>>>>>>>'
       write(stdout,'(a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~&
            &~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

       write(stdout,'(a)')'iter|       energy      |   rms grad  |&
            &  commutator |LScoef|CGcoef|      Ne    |'
    endif

    do iteration=1,maxit_pen

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(2(a,i3))') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: iteration ',iteration,' of ', &
            maxit_pen
#endif

       call services_flush

       do is=1,pub_cell%num_spins
          call sparse_copy(trial_denskern(is),denskern(is))
       end do

       !=========================== ENERGY AT POINT 0 =========================
       total_energy = internal_energy(trial_denskern, .true.)
       penalty_energy = internal_penalty_functional(trial_denskern)
       Qinitial = total_energy + penalty_energy ! for line minimisation

       if (pub_output_detail == VERBOSE .and. pub_on_root) &
            write(stdout,'(3(a,f20.12))') '      Q0=E0+P0: ',Qinitial,'=', &
            total_energy,'+',penalty_energy

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &initial functional evaluation complete'
#endif

       !======================== END ENERGY AT POINT 0 ========================

       !*************** CALCULATE KOHN-SHAM MATRIX ****************************
      call hamiltonian_build_matrix(ham,rep)

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &Kohn-Sham matrix constructed'
#endif

       !*********** END CALCULATE KOHN-SHAM MATRIX ****************************


       ! %%%%%%%%%%%%%%%%%%%%%%%% GRADIENT AT POINT 0 %%%%%%%%%%%%%%%%%%%%%%%%%

       rms_gradient = internal_gradient_norm()
       ! pdh: hack to make consistent with spin-unpolarised case
       rms_gradient = rms_gradient * pub_cell%num_spins


#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &initial gradient calculated'
#endif

       ! %%%%%%%%%%%%%%%%%%%%%%% END GRADIENT AT POINT 0 %%%%%%%%%%%%%%%%%%%%%%


       ! ########################## TEST CONVERGENCE ##########################

       converged = internal_test_convergence()
       call comms_bcast(pub_root_node_id,converged)
       if (converged) exit

       ! ######################## END TEST CONVERGENCE ########################


       !============================= LINE SEARCH =============================

       if (pub_output_detail == VERBOSE .and. pub_on_root) &
            write(stdout,'(a,f20.12)') '    Trial step: ', trial_step

       call internal_search_direction

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &search direction constructed'
#endif

       ! Take a trial step along the search direction
       do is=1,pub_cell%num_spins
          call sparse_copy(trial_denskern(is),denskern(is))
          call sparse_axpy(trial_denskern(is),con_direction(is),trial_step)
       end do
       ! ndmh: mark kernel workspace variables as invalid
       call kernel_workspace_invalidate()

       ! %%%%%%%%%%%%%%%%%%%%% ENERGY AT TRIAL LENGTH %%%%%%%%%%%%%%%%%%%%%%

       total_energy_at_trial_step = internal_energy(trial_denskern,.false.)
       penalty_energy = internal_penalty_functional(trial_denskern)
       Qtrial = total_energy_at_trial_step + penalty_energy ! line minimisation

       if (pub_output_detail == VERBOSE .and. pub_on_root) &
            write(stdout,'(3(a,f20.12))') '      Q1=E1+P1: ',Qtrial,'=', &
            total_energy_at_trial_step,'+',penalty_energy

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &trial functional evaluation complete'
#endif

       ! %%%%%%%%%%%%%%%%%%% END ENERGY AT TRIAL LENGTH %%%%%%%%%%%%%%%%%%%%

       ! Fit a quadratic to find the optimal step

       optimal_step = internal_quadratic_step()

       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !pdh: If the predicted step is negative, the parabolic fit is poor,
       !     so take a second trial step and fit a cubic instead
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

       if (optimal_step < 0.0_DP) then

          if (pub_output_detail == VERBOSE .and. pub_on_root) &
               write(stdout,'(a)')'    !!! Quadratic fit unsuccessful; &
               &proceeding to cubic fit !!!'

          ! Take a second trial step

          trial2 = .true.
          optimal_step = 2.0_DP * trial_step
          do is=1,pub_cell%num_spins
             call sparse_copy(trial_denskern(is),denskern(is))
             call sparse_axpy(trial_denskern(is),con_direction(is),optimal_step)
          end do
          ! ndmh: mark kernel workspace variables as invalid
          call kernel_workspace_invalidate()

       ! ndmh: better protection against dodgy optimal steps
       else if ((optimal_step < trial_step*0.01_DP) .or. &
            (optimal_step > trial_step*20.0_DP)) then

          if (pub_output_detail == VERBOSE .and. pub_on_root) &
               write(stdout,'(a)')'    !!! Quadratic fit may be unreliable; &
               &taking second trial step !!!'

          ! Take a second trial step

          trial2 = .true.
          do is=1,pub_cell%num_spins
             call sparse_copy(trial_denskern(is),denskern(is))
             call sparse_axpy(trial_denskern(is),con_direction(is),optimal_step)
          end do
          ! ndmh: mark kernel workspace variables as invalid
          call kernel_workspace_invalidate()


       else

          trial2 = .false.
          if (pub_output_detail == VERBOSE .and. pub_on_root) &
               write(stdout,'(a,f20.12)') ' Predicted min: ',Qpredict

       end if

       if (trial2) then ! if a second trial step is necessary

          ! %%%%%%%%%%%%%%%%%%% ENERGY AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%%%

          total_energy_at_trial_step = internal_energy(trial_denskern,.false.)
          penalty_energy = internal_penalty_functional(trial_denskern)
          Qoptimal = total_energy_at_trial_step + penalty_energy

          if (pub_output_detail == VERBOSE .and. pub_on_root) &
               write(stdout,'(3(a,f20.12))') '      Q2=E2+P2: ',Qoptimal,'=', &
               total_energy_at_trial_step,'+',penalty_energy

#ifdef DEBUG
          if (pub_on_root) then
             write(stdout,'(a)') 'DEBUG: In penalty_denskernel_optimise_cg: &
                  &second trial functional evaluation complete'
          end if
#endif

          ! %%%%%%%%%%%%%%%% END ENERGY AT TRIAL LENGTH 2 %%%%%%%%%%%%%%%%%%

          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          ! pdh: Now have sufficient information to fit a cubic
          ! ndmh: (as long as we have not already decided to re-do trial step
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

          if (optimal_step == trial_step*2.0_DP) then
            optimal_step = internal_cubic_step()
          else
            ! ndmh: find new optimal step, ignoring first trial step
            Qtrial = Qoptimal
            trial_step = optimal_step
            optimal_step = internal_quadratic_step()
          end if

       else

          ! ndmh: prevent exceedingly long quadratic steps
          if (optimal_step > 3.0_DP) then
             if (pub_on_root .and. (pub_output_detail >= NORMAL)) then
                write(stdout,'(a,f16.12)') &
                    'Calculated optimal quadratic step=',optimal_step
                write(stdout,'(a/a)') &
                    'WARNING in penalty_denskernel_optimise_cg: ', &
                    'setting quadratic optimal_step to safe value'
             end if
             optimal_step = 0.15_DP
             cg_coeff = 0.0_DP
             cg_count = 0
          else if ( (pub_output_detail == VERBOSE) .and. pub_on_root) then
             ! cks: If quadratic step accepted, print verbose info
             write(stdout,'(a)') '    !!! Quadratic minimum step &
                 &length accepted !!!'

          end if

       end if

       ! Update trial step for next line minimisation

       if (optimal_step > 0.0_DP) &
            trial_step = max(sqrt(trial_step*optimal_step),epsilon(1.0_DP))

       line_search_coeff = optimal_step

       ! Set the new density kernel
       do is=1,pub_cell%num_spins
          call sparse_axpy(denskern(is),con_direction(is),line_search_coeff)
       end do
       ! ndmh: mark kernel workspace variables as invalid
       call kernel_workspace_invalidate()

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: In &
            &penalty_denskernel_optimise_cg: &
            &optimal step taken'
#endif

       ! Store gradients for determination of next CG coefficient
       do is=1,pub_cell%num_spins
          call sparse_copy(old_con_direction(is),con_direction(is))
          call sparse_copy(old_co_gradient(is),co_gradient(is))
       end do

       ! ========================== END OF LINE SEARCH ========================

       call internal_print()

    end do

    ! Print out information concerning occupation number errors etc

    total_energy = internal_energy(denskern,.false.)
    penalty_energy = internal_penalty_functional(denskern)
    call kernel_occupancy_bounds(max_occ(:),min_occ(:),denskern,rep%overlap)
    if (pub_on_root) then
       write(stdout,'(/a)') '======================== Penalty functional &
            &minimisation ======================='
       write(stdout,'(a,f20.12)') '   Final total energy     : ',total_energy
       write(stdout,'(a,f20.12)') '   Final penalty energy   : ',penalty_energy
       write(stdout,'(a,f14.6)') '   RMS occupancy error    : ',&
            sqrt(penalty_energy / (pen_param*ngwf_basis%num))

       if (pub_cell%num_spins == 1) then
          write(stdout,'(a,2(f6.3,a))') '   Occupancy bounds       : [', &
               min_occ(1),',',max_occ(1),']'
       else
          do is=1,pub_cell%num_spins
             write(stdout,'(a,i1,a,2(f6.3,a))') '   Occupancy bounds spin ', &
                  is,': [',min_occ(is),',',max_occ(is),']'
          end do
       end if
       write(stdout,'(a/)') '=============================================&
            &==================================='

       if (minval(min_occ(:)) < 0.5_DP*(1.0_DP-sqrt(5.0_DP)) .or. &
            maxval(max_occ(:)) > 0.5_DP*(1.0_DP+sqrt(5.0_DP))) then
          write(stdout,'(/a)') &
               'WARNING: one or more occupancies may be outside the stable range!'
          write(stdout,'(a,i5,a,f10.4)') &
               'WARNING: Recommend maxit_pen >= ',2*maxit_pen, &
               ' and penparam >= ',2.0_DP*pen_param
       end if
    end if

    ! Fix occupancies if necessary
    if (minval(min_occ(:)) < 0.5_DP*(1.0_DP-sqrt(3.0_DP)) .or. &
         maxval(max_occ(:)) > 0.5_DP*(1.0_DP+sqrt(3.0_DP))) then
       call kernel_fix(denskern, rep%overlap, rep%inv_overlap)
       call kernel_occupancy_bounds(max_occ(:),min_occ(:),denskern,rep%overlap)
    end if

    if (pub_on_root) then
       if (minval(min_occ(:)) < 0.5_DP*(1.0_DP-sqrt(3.0_DP)) .or. &
            maxval(max_occ(:)) > 0.5_DP*(1.0_DP+sqrt(3.0_DP))) then
          write(stdout,'(a)') &
               'WARNING: one or more occupancies remains outside the stable range!'
       end if
    end if

    ! Stop if failed to converge
    if (.not. converged .and. pub_on_root) then
        write(stdout,'(a,i4,a)') 'Finished density kernel iterations (', &
        maxit_pen, ')'
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! ** Deallocate structures for sparse matrices

    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(trial_denskern(is))
       call sparse_destroy(old_co_gradient(is))
       call sparse_destroy(con_direction(is))
       call sparse_destroy(con_gradient(is))
       call sparse_destroy(co_gradient(is))
       call sparse_destroy(old_con_direction(is))
    end do

    deallocate(trial_denskern,stat=ierr)
    call utils_dealloc_check('penalty_denskernel_optimise_cg', &
         'trial_denskern',ierr)
    deallocate(old_co_gradient,stat=ierr)
    call utils_dealloc_check('penalty_denskernel_optimise_cg', &
         'old_co_gradient',ierr)
    deallocate(old_con_direction,stat=ierr)
    call utils_dealloc_check('penalty_denskernel_optimise_cg', &
         'old_con_direction',ierr)
    deallocate(con_direction,stat=ierr)
    call utils_dealloc_check('penalty_denskernel_optimise_cg', &
         'con_direction',ierr)
    deallocate(con_gradient,stat=ierr)
    call utils_dealloc_check('penalty_denskernel_optimise_cg','con_gradient', &
         ierr)
    deallocate(co_gradient,stat=ierr)
    call utils_dealloc_check('penalty_denskernel_optimise_cg','co_gradient', &
         ierr)

    ! ndmh: deallocate kernel work arrays
    call kernel_workspace_deallocate()

    ! cks: output density kernel to file if this is requested
    if (write_denskern.and.(.not.pub_write_converged_dk_ngwfs)) &
         call restart_kernel_write(denskern)

    ! Stop the timer
    call timer_clock('penalty_denskernel_optimise_cg',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &penalty_denskernel_optimise_cg'
#endif

    call services_flush

    return

  contains

    !******************************************************
    !******************************************************

    logical function internal_test_convergence()

      ! ########################## TEST CONVERGENCE ##########################

      internal_test_convergence = (rms_gradient < pen_threshold) .and. &
           (iteration > 0)

      if (internal_test_convergence) then
         if (pub_on_root) then
            write(stdout,'(i2,f20.14,2f15.12)') iteration,total_energy,&
                 rms_gradient,commutator
            write(stdout,'(/a)')'                 ......................&
                 &........................'
            write(stdout,'(a,f19.14,a)')'                 | RMS PF GRADIENT   = ',&
                 rms_gradient,'    |'
            write(stdout,'(a)') '                 | PF density kernel &
                 &optimisation converged!  |'
            write(stdout,'(a)') '                 ~~~~~~~~~~~~~~~~~~~~~~&
                 &~~~~~~~~~~~~~~~~~~~~~~~~'
         end if
      endif

      ! ######################## END TEST CONVERGENCE ########################

    end function internal_test_convergence


    !******************************************************
    !******************************************************


    real(kind=DP) function internal_energy(current_denskern,update_ham)

      ! Argument
      type(SPAM3), intent(inout) :: current_denskern(pub_cell%num_spins)
      logical, intent(in) :: update_ham

      ! Local variables
      integer :: is
      real(kind=DP) :: hubbard_energy
      real(kind=DP) :: paw_sphere_energies(paw_en_size)
#ifdef DEBUG
      real(kind=DP) :: ekin, enl
#endif

      ! Put spin degeneracy into density kernel to generate correct
      ! electron density
      do is=1,pub_cell%num_spins
         call sparse_scale(current_denskern(is),spin_fac)
      end do

      ! ndmh: calculate density dependent energies and matrices
      call hamiltonian_dens_dep_matrices(ham, lhxc_fine, internal_energy, &
           lhxc_energy, hubbard_energy, paw_sphere_energies, rep, ngwf_basis, &
           hub_proj_basis, hub, current_denskern, ewald_energy, elements, &
           pub_fine_grid, localpseudo_fine, core_density_fine, update_ham, &
           lhxc_fixed=.false.)

      ! Remove spin degeneracy from density kernel
      if (pub_cell%num_spins == 1) call sparse_scale(current_denskern(1),0.5_DP)

#ifdef DEBUG
      ekin = 0.0_DP ; enl = 0.0_DP
      do is=1,pub_cell%num_spins
         ekin = ekin + sparse_trace(current_denskern(is),rep%kinet)
      end do
      ekin = ekin * spin_fac
      if (pub_any_nl_proj) then
         do is=1,pub_cell%num_spins
            enl = enl + sparse_trace(current_denskern(is),rep%nonlocpot(1))
         end do
         enl = enl * spin_fac
      end if
      if (pub_aug) then
         do is=1,pub_cell%num_spins
            enl = enl + sparse_trace(current_denskern(is),ham%nonlocpot(is))
         end do
         enl = enl * spin_fac
      end if
      if (pub_on_root) then
         write(stdout,'(a,e24.16,a)') &
              'DEBUG: Kinetic energy                  :',ekin,' Ha'
         write(stdout,'(a,e24.16,a)') &
              'DEBUG: Nonlocal pseudopotential energy :',enl,' Ha'
         write(stdout,'(a,e24.16,a)') &
              'DEBUG: Hartree & XC energy             :',lhxc_energy,' Ha'
         if ( pub_hubbard ) write(stdout,'(a,e24.16,a)') &
              'DEBUG: Hubbard DFT+U energy        :',hubbard_energy,' Ha'
      end if
#endif

    end function internal_energy


    !******************************************************
    !******************************************************


    real(kind=DP) function internal_penalty_functional(trial_denskern)

      ! Arguments
      type(SPAM3), intent(in) :: trial_denskern(pub_cell%num_spins)

      ! Local variables
      integer :: is
      !integer :: ierr,i,j,num
      real(kind=DP) :: pen
      type(SPAM3) :: ksksc,ksc

      ! Allocate local arrays
      call sparse_create(ksc,pub_ks(1))
      call sparse_create(ksksc,pub_ks(1),pub_ks(1))

      pen = 0.0_DP

      call kernel_validate_ks(trial_denskern,rep%overlap)

      do is=1,pub_cell%num_spins

         ! Calculate penalty funtional
         ! ndmh: pub_ks is still valid
         call sparse_copy(ksc,pub_ks(is))
         call sparse_scale(ksc,-1.0_DP,1.0_DP)               ! I - K.S
         call sparse_product(ksksc,pub_ks(is),ksc)           ! K.S.(I - K.S)
         pen = pen + sparse_trace(ksksc,ksksc)

      end do

      ! Deallocate local arrays
      call sparse_destroy(ksksc)
      call sparse_destroy(ksc)

      internal_penalty_functional = 0.5_DP * spin_fac * pen_param * pen

    end function internal_penalty_functional



    !******************************************************
    !******************************************************


    real(kind=DP) function internal_gradient_norm()

      implicit none

      ! Local variable
      integer :: is
      real(kind=DP) :: num_els

      ! COVARIANT gradient (not tensor-corrected)
      call internal_co_gradient()

      ! CONTRAVARIANT (tensor-corrected) gradient
      do is=1,pub_cell%num_spins
         call sparse_copy(con_gradient(is),co_gradient(is))
      end do
      call internal_con_gradient(con_gradient)

      ! Project contravariant gradient to conserve electron number
      call internal_correct_ne_gradients(co_gradient,con_gradient, &
           rep%inv_overlap,rep%overlap)

      ! RMS value of gradient is a measure of the error in the current
      ! density kernel
      internal_gradient_norm = 0.0_DP
      num_els = 0.0_DP
      do is=1,pub_cell%num_spins
         internal_gradient_norm = internal_gradient_norm + &
              sparse_trace(con_gradient(is),co_gradient(is))
         num_els = num_els + sparse_num_element(denskern(is))
      end do
      internal_gradient_norm = sqrt(internal_gradient_norm / num_els)

      ! AAM: diagnostics
      commutator = kernel_rms_commutator(denskern,ham%ham,rep%overlap)

    end function internal_gradient_norm


    !******************************************************
    !******************************************************


    subroutine internal_search_direction

      implicit none

      ! Local variable
      integer :: is

      if (iteration > 1 .and. rms_gradient <= 1.01_DP) then
         cg_coeff = internal_polak_cg_coeff(old_con_direction, &
              old_co_gradient,co_gradient,con_gradient)
      else
         cg_coeff = 0.0_DP
      end if

      if (abs(cg_coeff) > eps) then
         cg_count = cg_count + 1
      else
         cg_count = 0
      end if

      ! reset cg_coeff every 'cg_max' steps
      if (cg_count > cg_max) then
         cg_coeff = 0.0_DP
         cg_count = 0
      end if

      ! Determine the line search direction
      do is=1,pub_cell%num_spins
         call sparse_copy(con_direction(is),con_gradient(is))
         ! pdh: hack to keep step lengths the same for spin polarised
         ! pdh: systems: multiply by a factor of two
         call sparse_scale(con_direction(is),-real(pub_cell%num_spins,kind=DP))
         call sparse_axpy(con_direction(is),old_con_direction(is),cg_coeff)
      end do

      ! Find the derivative of the function along the search direction
      Qslope = 0.0_DP
      do is=1,pub_cell%num_spins
         Qslope = Qslope + sparse_trace(co_gradient(is),con_direction(is))
      end do

      ! Check search direction points downhill!!!

      if (Qslope > 0.0_DP) then
         ! First try steepest descent direction i.e. reset conjugate gradients
         if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
              write(stdout,'(a)') '    !!! WARNING: positive line search &
              &gradient: resetting CG !!!'
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
            do is=1,pub_cell%num_spins
               call sparse_scale(con_direction(is),-1.0_DP)
            end do
            Qslope = -Qslope
         end if
      end if

      if (pub_output_detail == VERBOSE .and. pub_on_root) &
           write(stdout,'(a,f20.12)') '            G0: ',Qslope

    end subroutine internal_search_direction


    !******************************************************
    !******************************************************


    subroutine internal_print

      implicit none

      ! Local variable
      character(55) :: fmt
      integer :: is
      real(kind=DP) :: ne_new

      ne_new = 0.0_DP
      do is=1,pub_cell%num_spins
         ne_new = ne_new + sparse_trace(denskern(is),rep%overlap)
      end do
      ne_new = ne_new * spin_fac

      if (pub_on_root) then

         ! ndmh: adapt format string so that it always displays sensible
         ! ndmh: numbers of digits for total energy and Ne
         write(fmt,'(a)')'(i2,'
         if(abs(total_energy)<100000_DP) then
            write(fmt,'(a,a)') trim(fmt),'f22.14'
         else
            write(fmt,'(a,a)') trim(fmt),'f22.12'
         end if
         write(fmt,'(a,a)') trim(fmt),',2f14.11,f7.4,f7.3,'
         if(ne_new<10000_DP) then
            write(fmt,'(a,a)') trim(fmt),'f14.8)'
         else if (ne_new<100000_DP) then
            write(fmt,'(a,a)') trim(fmt),'f14.7)'
         else
            write(fmt,'(a,a)') trim(fmt),'f14.5)'
         end if

         write(stdout,fmt) iteration, total_energy, rms_gradient, commutator, &
             line_search_coeff, cg_coeff, ne_new

      end if

    end subroutine internal_print


    !******************************************************
    !******************************************************


    real(kind=DP) function internal_quadratic_step()

      ! pdh: Three pieces of information: sufficient to fit a quadratic
      ! cks: Q=ax^2+bx+c
      ! cks: c=Qinitial
      linear_term =  Qslope * trial_step  ! cks: bx
      quadratic_term = Qinitial + linear_term - Qtrial  ! cks: ax^2

      if (pub_output_detail == VERBOSE .and. pub_on_root) then
         write(stdout,'(a,f20.12)') 'Quadratic term: ',quadratic_term
         write(stdout,'(a,f20.12)') '   Linear term: ',linear_term
      end if

      if (abs(quadratic_term / Qinitial) > epsilon(1.0_DP)) then
         internal_quadratic_step = &
              0.5_DP * trial_step * linear_term / quadratic_term ! b/2a
         if (pub_output_detail == VERBOSE .and. pub_on_root) &
              write(stdout,'(a,f20.12)') '  Optimal step: ',internal_quadratic_step
         Qpredict = Qinitial + &
              0.25_DP*linear_term**2/quadratic_term ! c+b^2/2a
      else
         internal_quadratic_step = -1.0_DP
         if (pub_output_detail == VERBOSE .and. pub_on_root) &
              write(stdout,'(a)') '     !!! Quadratic term too small;' &
              // 'proceeding to cubic fit !!!'
      end if

    end function internal_quadratic_step


    !******************************************************
    !******************************************************


    real(kind=DP) function internal_cubic_step()

      ! Introduce shorthand working variables:
      !   xx & yy are the points at which the energy has been evaluated
      !   aa and bb are the quadratic and cubic polynomial coefficients
      ! cks: here's what PDH's code means: f= a*x^3+b*x^2+c*x+d
      ! cks: d=Qinitial=f1, c=Qslope=f1', b=aa, a=bb,
      !      f2=Qtrial, f3=Qoptimal

      internal_cubic_step = 0.15_DP
      xx = trial_step ; yy = optimal_step
      aa = ((yy*yy*yy*Qtrial - xx*xx*xx*Qoptimal)/(yy-xx) &
           - (xx*xx+xx*yy+yy*yy)*Qinitial - xx*yy*(xx+yy)*Qslope) &
           / (xx*xx*yy*yy)
      bb = ((yy*yy*Qtrial - xx*xx*Qoptimal)/(xx-yy) &
           + (xx+yy)*Qinitial + xx*yy*Qslope) / (xx*xx*yy*yy)
      disc=-1.0_DP
      if (abs(bb*yy/aa) > epsilon(1.0_DP)) then    ! avoid div by zero
         aon3b = aa / (3.0_DP * bb)
         disc = aon3b * aon3b - Qslope / (3.0_DP * bb)  ! discriminant
         if (disc >= 0.0_DP) then
            minimum_step = -aon3b + sign(sqrt(disc), bb)
            if (pub_output_detail == VERBOSE .and. pub_on_root) &
                 write(stdout,'(a,f20.12)') 'Cubic LS coeff: ',minimum_step
            ! cks: energy according to cubic minimum
            Qminimum = Qinitial + minimum_step * &
                 (Qslope + minimum_step * (aa + minimum_step * bb))
            if (pub_output_detail == VERBOSE .and. pub_on_root) &
                 write(stdout,'(a,f20.12)') ' Predicted min: ',Qminimum
            ! cks: set optimal_step to cubic minimum
            internal_cubic_step = minimum_step
            if (pub_output_detail == VERBOSE .and. pub_on_root) &
                 write(stdout,'(a)') '    !!! Cubic minimum step &
                 &length accepted !!!'
         end if
      end if

      ! Choose safe length if unsuccessful

      if (disc < 0.0_DP ) then
         if (pub_on_root .and. (pub_output_detail >= NORMAL)) then
            write(stdout,'(a)') '    !!! WARNING: Line search was &
                 &unsuccessful: coefficient set manually !!!'
         end if
         !internal_cubic_step = 0.15_DP
      end if

    end function internal_cubic_step


    !******************************************************
    !******************************************************

    ! Written by Arash Mostofi, March 2003
    ! Adapted for penalty functional by Peter Haynes, July 2003
    ! Rearranged to use pre-calculated pub_ks and thus save 1 sparse_product
    ! by Nick Hine, Nov 2007
    subroutine internal_co_gradient

      ! Local variables
      integer :: is
      type(SPAM3) :: sks,ksks

      ! Allocate workspace
      call sparse_create(sks,rep%overlap,pub_ks(1))
      call sparse_create(ksks,pub_ks(1),pub_ks(1))

      call kernel_validate_ks(denskern,rep%overlap)

      do is=1,pub_cell%num_spins

         ! Calculate gradient := H + (alpha/2) S.K.(I - 2 S.K).(I - S.K).S
         call sparse_product(ksks,pub_ks(is),pub_ks(is)) ! KSKS := KS.KS
         call sparse_scale(ksks,2.0_DP,1.0_DP)         ! KSKS := 1+2KS.KS
         call sparse_axpy(ksks,pub_ks(is),-3.0_DP)     ! KSKS := 1-3KS+2KS.KS
         call sparse_product(sks,rep%overlap,pub_ks(is))   ! SKS := S.KS
         call sparse_product(co_gradient(is),sks,ksks) ! G := SKS.(1-3KS+2KSKS)
         call sparse_scale(co_gradient(is),spin_fac * pen_param)
         call sparse_axpy(co_gradient(is),ham%ham(is),spin_fac)

      end do

      ! Deallocate workspace
      call sparse_destroy(ksks)
      call sparse_destroy(sks)

    end subroutine internal_co_gradient


    !******************************************************
    !******************************************************


    subroutine internal_con_gradient(gradient)

      type(SPAM3), intent(inout) :: gradient(pub_cell%num_spins)

      ! Local variables
      integer :: is
      type(SPAM3) :: sinvg

      ! Allocate workspace
      call sparse_create(sinvg,rep%inv_overlap,gradient(1))

      do is=1,pub_cell%num_spins

         ! Calculate contravariant gradient
         call sparse_product(sinvg,rep%inv_overlap,gradient(is))
         call sparse_product(gradient(is),sinvg,rep%inv_overlap)

      end do

      ! Deallocate workspace
      call sparse_destroy(sinvg)

    end subroutine internal_con_gradient


    !******************************************************
    !******************************************************


    subroutine internal_correct_ne_gradients(co_gradient, &
         con_gradient, inv_overlap, overlap)

      type(SPAM3), intent(inout) :: co_gradient(pub_cell%num_spins)
      type(SPAM3), intent(inout) :: con_gradient(pub_cell%num_spins)
      type(SPAM3), intent(in)    :: inv_overlap
      type(SPAM3), intent(in)    :: overlap

      ! Local variables
      integer :: is
      real(kind=DP) :: step
      real(kind=DP) :: rate_of_change_ne ! Rate of change of electron number
      real(kind=DP) :: trunc_rank

      ! Calculate amount of S^-1 to add to correct contravariant gradient
      trunc_rank = sparse_trace(overlap,inv_overlap) * pub_cell%num_spins
      rate_of_change_ne = 0.0_DP
      do is=1,pub_cell%num_spins
         rate_of_change_ne = rate_of_change_ne + &
              sparse_trace(overlap,con_gradient(is))
      end do
      step = -rate_of_change_ne / trunc_rank

      do is=1,pub_cell%num_spins

         ! Correct contravariant gradient
         call sparse_axpy(con_gradient(is),inv_overlap,step)

         ! Correct covariant gradient
         call sparse_axpy(co_gradient(is),overlap,step)

      end do

    end subroutine internal_correct_ne_gradients



    !******************************************************
    !******************************************************

    ! cks: written by cks on 29/6/2001
    real(kind=DP) function internal_polak_cg_coeff(old_con_direction, &
       old_co_gradient,co_gradient,con_gradient)

      ! cks: this subroutine returns the coefficient according to Polak for the
      ! cks: combination of the previous search direction with the negative of
      ! cks: the gradient in the conjugate gradient method.
      ! cks: b_{r+1}=g_{r+1}*(g_{r+1}-g_r)/( p_r* (g_{r+1}-g_r) )
      ! cks: It is written so that covariant vectors are contracted with
      ! cks: contravariant vectors.

      type(SPAM3), intent(in) :: old_con_direction(pub_cell%num_spins)
      type(SPAM3), intent(in) :: old_co_gradient(pub_cell%num_spins)
      type(SPAM3), intent(in) :: co_gradient(pub_cell%num_spins)
      type(SPAM3), intent(in) :: con_gradient(pub_cell%num_spins)

      ! cks: internal declarations
      integer :: is
      real(kind=DP) :: denominator

      denominator = 0.0_DP
      do is=1,pub_cell%num_spins
         denominator = denominator + &
              sparse_trace(old_con_direction(is),co_gradient(is)) - &
              sparse_trace(old_con_direction(is),old_co_gradient(is))
      end do
      ! pdh: hack to maintain consistency with spin polarisation
      denominator = denominator / pub_cell%num_spins

      if (abs(denominator) > 1.0e-6_DP) then ! avoid div by zero

         internal_polak_cg_coeff = 0.0_DP
         do is=1,pub_cell%num_spins
            internal_polak_cg_coeff = internal_polak_cg_coeff + &
                 sparse_trace(con_gradient(is),co_gradient(is)) - &
                 sparse_trace(con_gradient(is),old_co_gradient(is))
         end do
         internal_polak_cg_coeff = internal_polak_cg_coeff / denominator

      else

         if (pub_on_root .and. pub_output_detail >= NORMAL) then
            write(stdout,'(a)') &
                 'WARNING: zero denominator in internal_polak_cg_coef'
            write(stdout,'(a)') &
                 'WARNING: setting internal_polak_cg_coeff to zero'
         end if
         internal_polak_cg_coeff = 0.0_DP

      end if

      if (abs(internal_polak_cg_coeff) > 5.0_DP) then
         if (pub_on_root .and. pub_output_detail >= NORMAL) then
            write(stdout,'(a,f8.5)') 'WARNING: internal_polak_cg_coeff = ', &
                 internal_polak_cg_coeff
            write(stdout,'(a)') &
                 'WARNING: setting internal_polak_cg_coeff to zero'
         end if
         internal_polak_cg_coeff = 0.0_DP
      end if

    end function internal_polak_cg_coeff

  end subroutine penalty_denskernel_optimise_cg

end module penalty
