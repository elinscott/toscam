! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!            Electronic energy optimisation module               !
!                                                                !
! This module contains routines used in the optimisation of the  !
! electronic energy with respect to the NGWFs and the density    !
! kernel.                                                        !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in 2000.           !
! Subsequent modifications by Chris-Kriton Skylaris,             !
! Arash A. Mostofi, Peter D. Haynes and Nicholas D.M. Hine.      !
! Modified by Nicholas Hine to use the NGWF_REP container for    !
! NGWF-related data and matrices, October 2010.                  !
!================================================================!


module electronic

  use constants, only: DP

  implicit none

  private

  public :: electronic_init_denskern
  public :: electronic_energy
  public :: electronic_lagrangian

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  real(kind=DP) function electronic_energy(denskern, pur_denskern, ham, &
       lhxc_fine, mu, rep, ngwf_basis, hub_proj_basis, hub, &
       localpseudo_fine, core_density_fine, ewald_energy, elements, &
       lnv_threshold, current_maxit_lnv, kernel_update)

    !=============================================================!
    ! This function, given a set of NGWFs, returns the total      !
    ! energy after first optimising the density kernel in the     !
    ! LNV scheme (if required).                                   !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2003.              !
    ! Spin polarised by Peter Haynes, July 2006                   !
    ! Modified for NLCC by Nicholas Hine, January 2009            !
    ! DFT+U added by David O'Regan, April 2008                    !
    ! Adapted for SPAM3 and function basis by Nicholas Hine,      !
    ! May-July 2009.                                              !
    ! Adapted for NGWF_REP by Nicholas Hine in October 2010.      !
    ! Modified to calculated the conduction energy by Laura       !
    ! Ratcliff in October 2010.                                   !
    ! Kernel DIIS added by Alvaro Ruiz Serrano, November 2010.    !
    !=============================================================!

    use cell_grid, only: pub_fine_grid, pub_std_grid
    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, max_spins, stdout, paw_en_size
    use function_basis, only: FUNC_BASIS
    use fourier, only: fourier_filter_cell
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix
    use ion, only: ELEMENT
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: kernel_purify, kernel_rescale
    use kernel_diis, only: kernel_diis_calculate, kernel_diis_rescale, &
         kernel_diis_ham_diag
    use lnv, only: lnv_denskernel_optimise_cg
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: exact_lnv, pub_cond_calculate, pub_kernel_diis, &
         pub_nlcc, maxit_lnv_coarse
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_scale, sparse_trace, sparse_copy, SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments (outputs)
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: pur_denskern(pub_cell%num_spins)
    type(NGWF_HAM), intent(inout) :: ham
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(inout) :: mu(max_spins)

    ! Arguments (inputs)
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(in) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(in) :: ewald_energy
    real(kind=DP), intent(in) :: lnv_threshold
    integer, intent(in) :: current_maxit_lnv
    logical, intent(in) :: kernel_update
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Local variables
    real(kind=DP) :: hubbard_energy
    real(kind=DP) :: lhxc_energy
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: homolumo(2,pub_cell%num_spins)
    integer :: is
    logical :: ham_update
    integer :: ierr
    real(kind=DP), allocatable :: localpseudo_std(:,:,:), &
         core_density_std(:,:,:), lhxc_std(:,:,:)


    ! ndmh: optimise density kernel L using an LNV variant
    if ((maxit_lnv_coarse > 0).and.(kernel_update)&
         .and.(.not.pub_kernel_diis)) then

       allocate(localpseudo_std(pub_std_grid%ld1, &
            pub_std_grid%ld2, pub_std_grid%max_slabs12),stat=ierr)
       call utils_alloc_check('lnv_denskernel_opt_cg_coarse', &
            'localpseudo_std',ierr)
       call fourier_filter_cell(localpseudo_fine,localpseudo_std, &
            pub_fine_grid, pub_std_grid)
       if (pub_nlcc) then
          allocate(core_density_std(pub_std_grid%ld1, &
               pub_std_grid%ld2, pub_std_grid%max_slabs12),stat=ierr)
          call utils_alloc_check('lnv_denskernel_opt_cg_coarse', &
               'core_density_std',ierr)
          call fourier_filter_cell(core_density_fine,core_density_std, &
               pub_fine_grid, pub_std_grid)
       else
          allocate(core_density_std(1,1,1),stat=ierr)
          call utils_alloc_check('lnv_denskernel_opt_cg_coarse', &
               'core_density_std',ierr)
       end if
       allocate(lhxc_std(pub_std_grid%ld1, pub_std_grid%ld2, &
            pub_std_grid%max_slabs12),stat=ierr)
       call utils_alloc_check('lnv_denskernel_opt_cg_coarse', &
            'lhxc_std',ierr)

       call lnv_denskernel_optimise_cg( &
            denskern, pur_denskern, ham, lhxc_std, mu, &
            electronic_energy, rep, ngwf_basis, hub_proj_basis, hub, &
            pub_std_grid, localpseudo_std, core_density_std, &
            ewald_energy, elements, lnv_threshold, maxit_lnv_coarse)

       deallocate(lhxc_std,stat=ierr)
       call utils_dealloc_check('lnv_denskernel_opt_cg_coarse', &
            'lhxc_std',ierr)
       deallocate(core_density_std,stat=ierr)
       call utils_dealloc_check('lnv_denskernel_opt_cg_coarse', &
            'core_density_std',ierr)
       deallocate(localpseudo_std,stat=ierr)
       call utils_dealloc_check('lnv_denskernel_opt_cg_coarse', &
            'localpseudo_std',ierr)

    end if

    ! cks: optimise density kernel L using an LNV variant
    if ((current_maxit_lnv > 0).and.(kernel_update)&
         .and.(.not.pub_kernel_diis)) then

       call lnv_denskernel_optimise_cg( &
            denskern, pur_denskern, ham, lhxc_fine, mu, &
            electronic_energy, rep, ngwf_basis, hub_proj_basis, hub, &
            pub_fine_grid, localpseudo_fine, core_density_fine, &
            ewald_energy, elements, lnv_threshold, current_maxit_lnv)

       ! cks: purify denskern and scale before Lagrangian evaluation
       call kernel_purify(pur_denskern, denskern, rep%overlap, &
            rep%inv_overlap, rep%n_occ)

       ! Set normalisation factor for revised LNV
       if (exact_lnv) call kernel_rescale(pur_denskern,rep%overlap, &
            rep%n_occ,silent=.true.)

       ! pdh: changed for spin polarisation
       if (pub_cell%num_spins==1) call sparse_scale(pur_denskern(1),2.0_DP)

    else if (pub_kernel_diis.and.kernel_update) then

       if (.not.pub_cond_calculate) then
          electronic_energy = 0.0_DP
          call kernel_diis_calculate(electronic_energy, mu, pur_denskern, &
               denskern, ham, lhxc_fine, ngwf_basis, hub_proj_basis, hub, &
               rep, localpseudo_fine, core_density_fine, ewald_energy, elements)
       else
          electronic_energy = 0.0_DP
          call kernel_diis_ham_diag(homolumo,pur_denskern, &
               ham%ham,rep%overlap,rep%n_occ)
          if (pub_cell%num_spins==1) call sparse_scale(pur_denskern(1),2.0_DP)
          do is=1,pub_cell%num_spins
             electronic_energy = electronic_energy + &
                  sparse_trace(ham%ham(is),pur_denskern(is))
          end do
       end if

    else ! ndmh: no density kernel optimisation

       ! ars: set flag for ham update
       ham_update = .false.
       if (pub_kernel_diis) ham_update = .true.

       if (.not.pub_kernel_diis) then

          ! cks: purify denskern
          call kernel_purify(pur_denskern, denskern, rep%overlap, &
               rep%inv_overlap, rep%n_occ, fixed_denskern=.true.)

          ! Set normalisation factor for revised LNV
          if (exact_lnv) call kernel_rescale(pur_denskern,rep%overlap, &
               rep%n_occ,silent=.true.)

          ! ndmh: scale purified denskern to include spin-degeneracy
          if (pub_cell%num_spins==1) call sparse_scale(pur_denskern(1),2.0_DP)

       else

          ! ars: rescale density kernel and pur_denskern
          call kernel_diis_rescale(denskern, rep%overlap, rep%n_occ)
          do is = 1, pub_cell%num_spins
             call sparse_copy(pur_denskern(is), denskern(is))
          end do
          if (pub_cell%num_spins==1) call sparse_scale(pur_denskern(1),2.0_DP)

       end if

       ! ndmh: calculate density dependent energies and matrices
       if (.not.pub_cond_calculate) then
          call hamiltonian_dens_dep_matrices(ham, lhxc_fine, &
               electronic_energy, lhxc_energy, hubbard_energy, &
               paw_sphere_energies, rep, ngwf_basis, &
               hub_proj_basis, hub, pur_denskern, ewald_energy, elements, &
               pub_fine_grid, localpseudo_fine, core_density_fine, &
               ham_update, lhxc_fixed=.false.)

          ! ars: update hamiltonian in case of kernel-DIIS
          if(pub_kernel_diis) call hamiltonian_build_matrix(ham,rep)

       else
          ! lr408: No need to recalculate any matrices, we just want the trace
          ! lr408: of the kernel with the projected Hamiltonian
          electronic_energy = 0.0_DP
          do is=1,pub_cell%num_spins
             electronic_energy = electronic_energy + &
                  sparse_trace(ham%ham(is),pur_denskern(is))
          end do
       end if

    endif


  end function electronic_energy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function electronic_lagrangian(Einput, overlap, denskern, &
       ham, mu, n_occ)

    !========================================================!
    ! This function returns the value of the LNV Lagrangian. !
    !--------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2003.         !
    ! Spin polarised by Peter Haynes, July 2006              !
    ! Kernel DIIS by Alvaro Ruiz Serrano, November 2010.     !
    !========================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, max_spins, stdout
    use kernel_diis, only: kernel_diis_lagrangian
    use rundat, only: exact_lnv, pub_kernel_diis
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_trace

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: Einput
    type(SPAM3),   intent(in) :: overlap
    type(SPAM3),   intent(in) :: denskern(pub_cell%num_spins)
    type(SPAM3),   intent(in) :: ham(pub_cell%num_spins)
    real(kind=DP), intent(in) :: mu(max_spins)
    integer,       intent(in) :: n_occ(max_spins)

    ! Local variable
    integer :: is
    real(kind=DP) :: spin_fac


    ! ars: f = spin_fac
    spin_fac = 2.0_DP / pub_cell%num_spins

    ! ars: W = electronic_lagrangian
    electronic_lagrangian = Einput

    if (pub_kernel_diis) then

       ! ars: W = E(K') - f*tr[H(K'SK'-K')], where K' = [Ne/tr(KS)]*K
       call kernel_diis_lagrangian(electronic_lagrangian, overlap, denskern, &
            ham)

    elseif (exact_lnv) then

       ! ars: W = E(K'), where K' = [Ne/tr(KS)]*K
       electronic_lagrangian = electronic_lagrangian

    elseif (.not.exact_lnv) then

       ! ars: W = E(K) - f*mu*(tr(LS)-Nocc)
       do is=1,pub_cell%num_spins
          electronic_lagrangian = electronic_lagrangian - &
               spin_fac * mu(is) * (sparse_trace(denskern(is),overlap) - &
               real(n_occ(is),kind=DP))
       end do

    else

       ! ars: incorrect logic for kernel minimisation - ONETEP stop
       if(pub_on_root) write(stdout,*) &
            "Error in electronic_lagrangian: incorrect kernel optimisation method.&
            & ONETEP stops."
       call comms_abort

    end if

  end function electronic_lagrangian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine electronic_init_denskern(rep, ham, denskern, localpseudo_fine, &
       core_density_fine, lhxc_fine, ngwf_basis, proj_basis, hub_proj_basis, &
       hub, elements, ewald_energy, val_rep, val_ngwf_basis, val_dkn, val_ham, &
       lhxc_fixed_in)

    !=================================================================!
    ! This subroutine initialises the denskern SPAM3, using one of    !
    ! several techniques: diagonalisation of the initial Hamiltonian, !
    ! Palser-Manolopoulos Canonical purification, or minimisation of  !
    ! the penalty functional. It also calls electronic_init_pot which !
    ! initialises the whole-cell arrays.                              !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 21/2/2004.                  !
    ! Modified by Peter D. Haynes on 5/7/2004 for Fourier             !
    ! parallelisation.                                                !
    ! Spin polarised by Peter Haynes, July 2006                       !
    ! Modified for NLCC by Nicholas Hine, January 2009                !
    ! DFT+U added by D. D. O'Regan, April 2009                        !
    ! Adapted for SPAM3 and function basis by Nicholas Hine in        !
    ! May-July 2009.                                                  !
    ! Modified by Jacek Dziedzic on 13/05/2010 to include correction  !
    ! due to smeared ions and implicit solvent.                       !
    ! Reorganised, made to use NGWF_REP and NGWF_HAM, split into two  !
    ! parts and generally tidied-up by Nicholas Hine in October 2010. !
    ! Modified by Laura Ratcliff in October 2010 to allow for         !
    ! use in conduction calculations.                                 !
    ! Modified by Simon Dubois in May 2010 to include various         !
    ! mixing schemes                                                  !
    !=================================================================!

    use augmentation, only: augmentation_screen_dij, aug_projector_denskern, &
         aug_nonlocal_mat
    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, max_spins, VERBOSE, stdout, DN, UP, ANGSTROM, &
         paw_en_size
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_indep_matrices, &
         hamiltonian_proj_cond_matrix, hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix
    use hubbard_build, only: HUBBARD_MODEL
    use integrals, only: integrals_locpot
    use ion, only: ELEMENT
    use kernel, only: kernel_init_core_ham, kernel_christoffel, &
         kernel_basis_transform, kernel_basis_update
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use palser_mano, only: palser_mano_kernel_optimise
    use paw, only: paw_projector_denskern_init
    use penalty, only: penalty_denskernel_optimise_cg
    use restart, only: restart_kernel_read, restart_kernel_write, &
         restart_kernel_compose, retrieve_denskern
    use rundat, only: pub_coulomb_cutoff, pub_output_detail, maxit_pen, &
         lnv_threshold_orig, write_denskern, pub_maxit_palser_mano, &
         read_denskern, cond_read_denskern, pub_any_nl_proj, pub_hubbard, &
         coreham_denskern_guess, pub_hubbard_restart, pub_aug, pub_paw, &
         pub_write_converged_dk_ngwfs, pub_nlcc, pub_cond_calculate, &
         pub_devel_code
    use services, only: services_flush
    use sparse, only: SPAM3, sparse_scale, sparse_copy, &
         sparse_axpy, sparse_create, sparse_destroy, sparse_trace, sparse_show_matrix
    use simulation_cell, only: pub_cell, pub_padded_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    implicit none

    ! Arguments (outputs)
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
    type(NGWF_REP), intent(inout) :: rep
    type(NGWF_HAM), intent(inout) :: ham
    real(kind=DP), intent(inout) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(inout) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)

    ! Arguments (inputs)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=DP), intent(in) :: ewald_energy

    ! lr408: Optional arguments needed for conduction calculation
    type(NGWF_REP), optional, intent(in) :: val_rep
    type(NGWF_HAM), optional, intent(in) :: val_ham
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis
    type(SPAM3), optional, intent(in)      :: val_dkn(pub_cell%num_spins)
    logical, optional, intent(in) :: lhxc_fixed_in

    ! Local Variables
    type(SPAM3), allocatable :: rho_ij(:)
    real(kind=DP) :: filling
    real(kind=DP) :: Ecurr
    real(kind=DP) :: num_elec
    integer :: is
    integer :: ierr
    logical :: init_pot, read_dkn, cond_call, lhxc_fixed

    real(kind=DP) :: electronic_energy, lhxc_energy, hubbard_energy
    real(kind=DP) :: paw_sphere_energies(paw_en_size)

!CW
    real(8) :: ndd
!ENDCW
    ! ------------------------------------------------------------------------------

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &electronic_init_denskern'
#endif

    ! lr408: If optional argument is present, we are initialising a conduction
    ! lr480: kernel, and so the potential does not need initialising.
    ! lr408: Also, read_dkn is set to match read_cond_denskern, rather
    ! lr408: than read_denskern
    if (present(val_rep)) then
       init_pot = .false.
       read_dkn = cond_read_denskern
       cond_call = .true.
       lhxc_fixed = .true.
    else
       init_pot = .true.
       read_dkn = read_denskern
       cond_call = .false.
       lhxc_fixed = .false.
    end if
    if (present(lhxc_fixed_in)) then
       lhxc_fixed = lhxc_fixed_in
       if (lhxc_fixed) init_pot = .false.
    end if

    ! ndmh: Initialise parts of Hamiltonian that do not depend on the density
    ! lr408: If this is a conduction call, we need to pass on the valence rep
    ! lr408: and basis
    call hamiltonian_dens_indep_matrices(rep, ngwf_basis, proj_basis, &
         hub_proj_basis, hub, val_rep, val_ngwf_basis) ! 2 optional arguments

    ! ndmh: Initialise the local pseudopotential, core density, and
    ! ndmh: if required, create the lhxc matrix in SPAM3 format
    if (init_pot) call electronic_init_pot(localpseudo_fine, &
         core_density_fine, lhxc_fine, elements)

    ! cks: initialise density kernel
    if (read_dkn .and. .not. retrieve_denskern) then

       ! ------- READ A DENSITY KERNEL FROM FILE -----------------------
       ! lr408: Pass extra flag to determine which kernel file to read from
       call restart_kernel_read(denskern, read_cond=cond_call)
!CW / CHECKING
#ifdef debug
       do is=1,pub_cell%num_spins
        ndd=sparse_trace(denskern(is),rep%overlap)
        write(*,*) 'CHECKING NUMBER OF ELECTRONS, SPIN / N : ', is,ndd
       enddo
#endif
!END CW
       ! ------- END READ A DENSITY KERNEL FROM FILE--------------------

       ! lr408: Adding safeguard to ensure that the kernel contains the 
       ! lr408: correct number of of 'occupied' states
       if (cond_call) then
          ! Calculate current electron number
          do is=1,pub_cell%num_spins
             num_elec = sparse_trace(denskern(is),rep%overlap)
             if (abs(num_elec/real(rep%n_occ(is),DP)-1.0_DP)>0.05_DP) then
                if (pub_on_root) then
                   write(stdout,'(a,f12.6)') 'Error in electronic_init_denskern: &
                        &number of ''occupied'' orbitals in input kernel = ', &
                        num_elec 
                   write(stdout,'(a,i6,a,i2)') ' which does not match the &
                        &expected value of ',rep%n_occ(is),' for spin=',is
                   write(stdout,'(a)') 'verify that cond_num_extra_states and &
                        &cond_num_extra_its are set correctly.'
                   call comms_abort
                end if
             end if
          end do
       end if


       if (index(pub_devel_code,'DIAGONALISE_ON_RESTART')>0) then
          if (pub_cell%num_spins==1) call sparse_scale(denskern(1),2.0_DP)
          call hamiltonian_dens_dep_matrices(ham, lhxc_fine, &
               electronic_energy, lhxc_energy, hubbard_energy, &
               paw_sphere_energies, rep, ngwf_basis, &
               hub_proj_basis, hub, denskern, ewald_energy, elements, &
               pub_fine_grid, localpseudo_fine, core_density_fine, &
               .true., lhxc_fixed=.false.)

          call hamiltonian_build_matrix(ham, rep)

          do is=1,pub_cell%num_spins
             call kernel_init_core_ham(denskern(is), &
                  ham%ham(is), rep%overlap, rep%n_occ(is))
          end do
       end if

    ! smmd: if required, retrieve density kernel from data in restart module
    else if (retrieve_denskern) then

       if (pub_on_root) write(stdout,'(a)') &
          'electronic_init_denskern : Initialise density kernel from stored &
          &data...  '

       call restart_kernel_compose(denskern, rep, ngwf_basis, &
                  proj_basis,elements)

       if (pub_on_root) write(stdout,'(a)') &
          'electronic_init_denskern : ...density kernel initialised ! '

    else if (coreham_denskern_guess) then

       ! ------- GUESS A DENSITY KERNEL --------------------------------

       ! ndmh: Initialise the PAW nonlocal matrix, if required
       if (pub_aug) then
          ! Allocate storage for projector density kernel
          allocate(rho_ij(pub_cell%num_spins),stat=ierr)
          call utils_alloc_check('electronic_init_pot','rho_ij',ierr)
          do is=1,pub_cell%num_spins
             rho_ij(is)%structure='E'
             call sparse_create(rho_ij(is))
          end do

          ! For conduction calculation, copy in nonlocal energies from val ham
          if (cond_call) then
             do is=1,pub_cell%num_spins
                call sparse_copy(ham%dijhat(is),val_ham%dijhat(is))
             end do
             call aug_projector_denskern(rho_ij,val_dkn,val_rep%sp_overlap)
             do is=1,pub_cell%num_spins
                call sparse_scale(rho_ij(is),pub_cell%spin_fac)
             end do
          end if

          ! If this is a brand-new calculation, guess a kernel
          if (init_pot) then
             if (pub_paw) then
                call paw_projector_denskern_init(rho_ij)
             end if
             call augmentation_screen_dij(ham%dijhat,lhxc_fine,pub_fine_grid)
          end if

          ! Now calculate nonlocal matrix
          call aug_nonlocal_mat(ham%nonlocpot,ham%dijhat,rho_ij, &
               rep%sp_overlap,show_matrices=.true.)

          ! Deallocate storage for projector density kernel
          do is=pub_cell%num_spins,1,-1
             rho_ij(is)%structure='E'
             call sparse_destroy(rho_ij(is))
          end do
          deallocate(rho_ij,stat=ierr)
          call utils_dealloc_check('electronic_init_pot','rho_ij',ierr)

       end if

       do is=1,pub_cell%num_spins

          if (cond_call) then
             if (pub_on_root) then
               write(stdout,'(a,I5,a)') 'Generating conduction density kernel &
                    &for ',rep%n_occ(is),' ''occupied'' orbitals'
             end if
          end if

          if (pub_on_root) then
             if (is == UP) then
                write(stdout,'(a)',advance='no') &
                     'Up spin density kernel initialisation ...'
             elseif (is == DN) then
                write(stdout,'(a)',advance='no') &
                     'Down spin density kernel initialisation ...'
             endif
          endif

          ! ndmh: calculate lhxc matrix
          call integrals_locpot(ham%lhxc(is), &
               rep%ngwfs_on_grid,ngwf_basis,rep%ngwfs_on_grid,ngwf_basis, &
               pub_fine_grid, lhxc_fine(:,:,:,is))

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

          ! lr408: If this is a conduction call, we need to calculate the
          ! lr408: projected conduction Hamiltonian
          ! lr408: Only needed for spin = 1, as both spins are done
          ! lr408: inside routine
          if (cond_call .and. is==1) then
             call hamiltonian_proj_cond_matrix(rep, ham, &
                  val_rep, val_ham, val_dkn)
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

       end do

       ! cks: output density kernel to file if this is requested
       if (write_denskern.and.(.not.pub_write_converged_dk_ngwfs)) &
            call restart_kernel_write(denskern,write_cond=cond_call)

       ! --- END GUESS A DENSITY KERNEL ----------------------------

    else ! last resort: scale S^-1

       ! --- SCALE INVERSE OVERLAP MATRIX --------------------------

       do is=1,pub_cell%num_spins
          filling = real(rep%n_occ(is),DP) / real(ngwf_basis%num,DP)
          call sparse_copy(denskern(is),rep%inv_overlap)
          call sparse_scale(denskern(is),filling)
          if (pub_output_detail == VERBOSE .and. pub_on_root) &
               write(stdout,'(a,i1,a,f6.3)') &
               ' Initialising density kernel for spin ',is, &
               ' to inverse overlap scaled by filling ',filling
       end do

       ! --- END SCALE INVERSE OVERLAP MATRIX ----------------------

    end if

    ! &&&&&&&& USE PENALTY FUNCTIONAL TO IMPROVE DENSITY KERNEL &&&&
    ! ddor: This is unnecessary when restarting with task HUBBARDSCF
    ! lr408: or when restarting with task COND
    if ((maxit_pen > 0) .and. (.not. pub_hubbard_restart ) .and. &
         (.not. pub_cond_calculate)) then

       call penalty_denskernel_optimise_cg(denskern, rep, ham, &
            lhxc_fine, Ecurr, localpseudo_fine, core_density_fine, &
            ngwf_basis, hub_proj_basis, hub, ewald_energy, elements, &
            lnv_threshold_orig)

    end if
    ! &&&&& END USE PENALTY FUNCTIONAL TO IMPROVE DENSITY KERNEL &&&&

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &electronic_init_denskern'
#endif

    call services_flush

  end subroutine electronic_init_denskern


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine electronic_init_pot(localpseudo_fine, core_density_fine, &
       lhxc_fine, elements)

    !===============================================================!
    ! This subroutine initialises the local pseudopotential, the    !
    ! core density, the initial electron density, and uses them to  !
    ! initialise the initial-guess Hamiltonian if required.         !
    !---------------------------------------------------------------!
    ! This subroutine was created by Nicholas Hine in October 2010. !
    ! It was based around parts of the previous routine             !
    ! electronic_init_denskern_pot, which was written by the ODG in !
    ! 2000-2010.                                                    !
    !===============================================================!

    use augmentation, only: aug_nonlocal_mat, augmentation_screen_dij
    use cell_grid, only: pub_fine_grid
    use classical_pot, only: classical_pot_struct_fac
    use comms, only: pub_on_root
    use constants, only: DP, max_spins, VERBOSE, stdout, DN, UP, ANGSTROM
    use cutoff_coulomb, only: cutoff_coulomb_struct_fac, &
         cutoff_coulomb_localpseudo, cutoff_coulomb_hartree, &
         cutoff_coulomb_initial_guess, cutoff_coulomb_core_density, &
         cutoff_coulomb_classical_sfac, pub_padded_grid
    use density, only: density_initial_guess_recip, density_initial_guess_real
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use integrals, only: integrals_locpot
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use paw, only: paw_tcore_hartree_on_grid, paw_tcore_density
    use potential, only: potential_sawtooth_efield
    use pseudopotentials, only: pseudo_make_structure_factor,&
         pseudopotentials_local_on_grid, pseudopotentials_core_density
    use rundat, only: pub_coulomb_cutoff, pub_output_detail,  &
         coreham_denskern_guess, pub_any_nl_proj, pub_constant_efield, &
         pub_hubbard, pub_nlcc, pub_paw, pub_ii_energy_direct, &
         pub_open_localpseudo, pub_is_implicit_solvent, &
         pub_is_smeared_ion_rep, pub_is_include_cavitation, &
         pub_is_dielectric_model, read_denskern, pub_multigrid_hartree, &
         pub_aug, pub_nhat_in_xc, pub_initial_dens_realspace
    use services, only: services_flush
    use sparse, only: SPAM3, sparse_copy, sparse_create, sparse_destroy
    use simulation_cell, only: pub_cell, pub_padded_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    use xc, only: xc_energy_potential

    implicit none

    ! Arguments (outputs)
    real(kind=DP), intent(out) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(out) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC

    real(kind=DP), intent(out) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)

    ! Arguments (inputs)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Local Variables
    integer :: is
    integer :: ierr
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: struct_fac
    complex(kind=DP), allocatable, dimension(:,:,:) :: struct_fac_classical !ars
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: pot_fine
    real(kind=DP) :: xc_energy

    ! Allocate storage for structure factor(s) and calculate it
    if (.not.pub_coulomb_cutoff) then

       allocate(struct_fac(pub_cell%num_pspecies, pub_fine_grid%ld3, &
            pub_fine_grid%ld2, pub_fine_grid%max_slabs23),stat=ierr)
       call utils_alloc_check('electronic_init_pot','struct_fac',ierr)

       ! calculate structure factor for all species
       call pseudo_make_structure_factor(struct_fac,  &    ! output
            elements,pub_fine_grid)                        ! input

       ! cks: allocate structure factor for classical atoms if it is needed
       if (pub_cell%nat_classical.gt.0) then

          allocate(struct_fac_classical(pub_fine_grid%ld3,pub_fine_grid%ld2, &
               pub_fine_grid%max_slabs23),stat=ierr)
          call utils_alloc_check('electronic_init_pot', &
               'struct_fac_classical',ierr)

          ! ars: calculate structure factor for the classical atoms
          call classical_pot_struct_fac(struct_fac_classical,pub_fine_grid)
       else
          allocate(struct_fac_classical(0,0,0),stat=ierr)

          call utils_alloc_check('electronic_init_pot', &
               'struct_fac_classical',ierr)
       end if

       ! calculate the local pseudopotential on the fine grid
       if (.not.pub_paw) then
          call pseudopotentials_local_on_grid(localpseudo_fine, struct_fac, &
               struct_fac_classical,pub_fine_grid,elements)
       else
          call paw_tcore_hartree_on_grid(localpseudo_fine, struct_fac, &
               struct_fac_classical,pub_fine_grid)
       end if

    else ! Cutoff Coulomb - use padding wrappers

       ! allocate padded structure factor
       allocate(struct_fac(pub_padded_cell%num_pspecies, &
            pub_padded_grid%ld3, pub_padded_grid%ld2, &
            pub_padded_grid%max_slabs23),stat=ierr)
       call utils_alloc_check('electronic_init_pot','struct_fac',ierr)

       ! calculate structure factor for all species
       call cutoff_coulomb_struct_fac(struct_fac, &   ! output
            elements)                                 ! input

       ! ndmh: allocate structure factor for classical atoms if it is needed
       if (pub_cell%nat_classical.gt.0) then
          allocate(struct_fac_classical(pub_padded_grid%ld3, &
               pub_padded_grid%ld2, pub_padded_grid%max_slabs23), &
               stat=ierr)
          call utils_alloc_check('electronic_init_pot', &
               'struct_fac_classical',ierr)

          ! ndmh: wrapper to calculate structure factor for the classical atoms
          call cutoff_coulomb_classical_sfac(struct_fac_classical)  ! output
       else
          allocate(struct_fac_classical(0,0,0),stat=ierr)
          call utils_alloc_check('electronic_init_pot', &
               'struct_fac_classical',ierr)
       endif

       ! calculate the local pseudopotential on the fine grid
       call cutoff_coulomb_localpseudo(localpseudo_fine, &
            struct_fac,struct_fac_classical,elements)

    end if

    ! ndmh: calculate the core charge density
    if (pub_nlcc) then

       if (pub_on_root .and. pub_output_detail>=VERBOSE) then
          write(stdout,'(a)',advance='no') 'Calculating core density ...'
       end if

       ! ndmh: use padding wrapper if using coulomb cutoff
       if (.not.pub_coulomb_cutoff) then
          if (.not.pub_paw) then
             call pseudopotentials_core_density(core_density_fine,struct_fac, &
                  pub_fine_grid)
          else
             call paw_tcore_density(core_density_fine,struct_fac,pub_fine_grid)
          end if
       else
          call cutoff_coulomb_core_density(core_density_fine,struct_fac)
       end if

       if (pub_on_root .and. pub_output_detail>=VERBOSE) then
          write(stdout,'(a)') ' ... done'
       end if

    end if

    ! cks: apply constant electric field as sawtooth contribution to local
    ! cks: potential
    if (any(pub_constant_efield(:) /= 0.0_DP)) then
       call potential_sawtooth_efield(localpseudo_fine,pub_fine_grid)
    end if

    if ((coreham_denskern_guess).and.(.not.read_denskern)) then

       ! Allocate temporary arrays for potential and density
       allocate(density_fine(pub_fine_grid%ld1, &
            pub_fine_grid%ld2, pub_fine_grid%max_slabs12, &
            pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('electronic_init_pot','density_fine',ierr)

       ! Calculate initial density guess if required
       if (pub_initial_dens_realspace) then
          call density_initial_guess_real(density_fine, &
               elements, pub_fine_grid, add_aug_den=.true.)
          allocate(pot_fine(pub_fine_grid%ld1, &
               pub_fine_grid%ld2, pub_fine_grid%max_slabs12, &
               pub_cell%num_spins),stat=ierr)
          call utils_alloc_check('electronic_init_pot','pot_fine',ierr)
      else

          if (pub_coulomb_cutoff) then
             call cutoff_coulomb_initial_guess(density_fine, &
                  elements, struct_fac)
          else
             call density_initial_guess_recip(density_fine, &
                  elements, struct_fac, pub_fine_grid)
          end if
       end if

       ! cks: density_fine_init is converted to Hartree potential
       ! jd: Calculate either:
       !     (1) the usual Hartree potential
       !     (2) Hartree potential with cutoff Coulomb
       !     (3) Hartree potential with multigrid, which could be
       !         a) electronic Hartree potential, with open BCs,
       !         b) molecular Hartree potential, with smeared-ions, in vacuo,
       !         c) molecular Hartree potential, in implicit solvent.
       !     All of the above are mutually exclusive.
       !     Iff (3c) and the dielectric is not fixed, but instead depends on
       !     the electronic density, an extra term caused by this dependence
       !     (V_eps(r)+V_cav(r)) is included in the potential to get the correct
       !     functional derivative. The Hartree energy then needs to be
       !     calculated differently, and hartree_via_multigrid() takes care of
       !     that. However, here the actual energy value is not needed.
       if (pub_coulomb_cutoff) then
          call cutoff_coulomb_hartree(lhxc_fine, density_fine)          ! (2)
       else if (pub_multigrid_hartree) then
          call hartree_via_multigrid(lhxc_fine, density_fine)           ! (3)
       else
          call hartree_on_grid(lhxc_fine, density_fine, pub_fine_grid)  ! (1)
       end if

       ! If we have a realistic density guess, we can include the XC potential
       if (pub_initial_dens_realspace) then

          if (pub_paw.and.(.not.pub_nhat_in_xc)) then
             ! Re-generate density but leave out augmentation density this time
             call density_initial_guess_real(density_fine, &
                  elements, pub_fine_grid, add_aug_den=.false.)
             ! Otherwise keep previous density
          end if

          ! ndmh: add on core density before calculating xc potential for NLCC
          if (pub_nlcc) then
             do is=1,pub_cell%num_spins
                density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
                     core_density_fine/real(pub_cell%num_spins,kind=DP)
             end do
          end if

          ! ndmh: no augmentation density to pass in
          call xc_energy_potential(density_fine, xc_energy, pot_fine, &
               pub_fine_grid, 0)

          lhxc_fine(:,:,:,:) = lhxc_fine(:,:,:,:) + pot_fine(:,:,:,:)

          ! Deallocate temporary density array
          deallocate(pot_fine,stat=ierr)
          call utils_dealloc_check('electronic_init_pot','pot_fine',ierr)

       end if

       ! ndmh: Add up local pseudo + Hartree (+ xc if calculated)
       do is=1,pub_cell%num_spins
          lhxc_fine(:,:,:,is) = lhxc_fine(:,:,:,is) + localpseudo_fine
       end do

       deallocate(density_fine,stat=ierr)
       call utils_dealloc_check('electronic_init_pot','density_fine',ierr)

    end if

    ! Deallocate storage for structure factor(s)
    deallocate(struct_fac,stat=ierr)
    call utils_dealloc_check('electronic_init_pot','struct_fac',ierr)
    deallocate(struct_fac_classical,stat=ierr)
    call utils_dealloc_check('electronic_init_pot','struct_fac_classical',ierr)

  end subroutine electronic_init_pot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module electronic
