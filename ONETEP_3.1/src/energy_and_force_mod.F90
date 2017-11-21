! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                      E N E R G Y    A N D    F O R C E                      !
!=============================================================================!
!                                                                             !
! This module calculates the total energy and forces, as well as converged    !
! ngwfs, density kernel and related quantities, for a given ionic             !
! configuration.                                                              !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
!                 The ONETEP code is written and maintained by                !
! Chris-Kriton Skylaris, Arash A. Mostofi, Nicholas Hine and Peter D. Haynes  !
!                                                                             !
!-----------------------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in 2000.                        !
! Improved and parallelised by Chris-Kriton Skylaris in November 2003.        !
! Improvements by Peter D. Haynes in 2004.                                    !
! Turned into a module by Arash A. Mostofi, Version 0.01, 01/11/2004.         !
! Modified by Chris-Kriton Skylaris on 08/11/2004 to work with ppd_strategy.  !
! Modified by Nicholas Hine on 11/08/2008 to remove legacy SPAM denskern      !
! Modified by Nicholas Hine in July 2009 for SPAM3 and function_basis.        !
! Modified by Nicholas Hine in November 2009 to remove use of workspace_mod.  !
! PAW functionality by Nicholas Hine, May-July 2010.                          !
! Conduction NGWF optimisation added by Laura Ratcliff in October 2010.       !
!-----------------------------------------------------------------------------!

module energy_and_force

  implicit none

  private

  public :: energy_and_force_calculate

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine energy_and_force_calculate(total_energy,total_forces,elements, &
       properties_only,return_converged)

    use augmentation, only: augmentation_box_init
    use cdft_cg, only: cdft_u_cg_optimise
    use cell_grid, only: cell_grid_distribute, cell_grid_exit, pub_dbl_grid, &
         pub_fine_grid, pub_std_grid
    use classical_pot, only: classical_elements, classical_pot_ii_energy
    use comms, only: comms_abort, comms_barrier, comms_reduce, pub_on_root, &
         pub_my_node_id, pub_total_num_nodes
    use constants, only: stdout, DP, UP, DN, max_spins, NORMAL, VERBOSE
    use conduction, only: conduction_ngwf_optimise
    use cutoff_coulomb, only: cutoff_coulomb_ii_energy, &
         cutoff_coulomb_init, cutoff_coulomb_exit
    use dense, only: dense_init, dense_exit
    use density, only: density_radial_init, density_radial_exit
    use electronic, only: electronic_init_denskern
    use ewald, only: ewald_calculate_energy
    use fourier, only: fourier_init, fourier_init_cell, fourier_exit, &
         fourier_exit_cell
    use forces, only: forces_calculate
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_init_tight_boxes, function_basis_init_spheres, &
         function_basis_gath_all_tbs, function_basis_init_uni_tb, &
         function_basis_exit_uni_tb, function_basis_estimate_size, &
         function_basis_copy_spheres
    use geometry, only: POINT
    use hf_exchange, only: hf_exchange_set_pub_hfxsw, hf_exchange_init_vmatrix,&
         hf_exchange_fill_vmatrix
    use hubbard_build, only: HUBBARD_MODEL, hubbard_model_init, &
         hubbard_model_exit, hubbard_species_proj, &
         hubbard_build_consist, hubbard_build_consist_exit, &
         hubbard_projector_consistency, hubbard_test_convergence, &
         hubbard_species_exit_proj, &
         hubbard_build_matrices, hubbard_build_matrices_exit
    use hubbard_init, only: hubbard_init_species_exit
    use ion, only: ELEMENT
    use is_solvation, only: implicit_solvent_exit
    use multigrid_methods, only: multigrid_initialise
    use kernel, only: pub_kernel_workspace_allocated
    use ngwf_cg, only: ngwf_cg_optimise
    use ngwf_diis, only: ngwf_diis_optimise
    use ngwf_representation, only: NGWF_REP, ngwf_rep_create, ngwf_rep_destroy, &
         NGWF_HAM, ngwf_ham_create, ngwf_ham_destroy
    use parallel_strategy, only: parallel_strategy_distr_atoms, &
         parallel_strategy_check_atoms, parallel_strategy_exit, &
         parallel_strategy_list_overlaps, pub_elements_on_node, &
         pub_num_atoms_on_node
    use paw, only: paw_read_species, paw_species_exit, paw_species_init_proj, &
         paw_projectors
    use pbc_corrections, only: pbc_corr_check_cell_size, pbc_corr_exit
    use potential, only: potential_add_efield_ion_energy
    use ppd_strategy, only: ppd_strategy_check_and_print, &
         ppd_strategy_fftbox_spec
    use projectors, only: projectors_create_real
    use properties, only: properties_calculate
    use pseudopotentials, only: pseudopotentials_species_exit, &
         pseudopotentials_read_species, pseudo_species_init_proj, &
         nlps_projectors
    use rundat, only: pub_any_nl_proj, pub_charge, pub_coulomb_cutoff, &
         pub_coulomb_cutoff_type, pub_dispersion, pub_do_properties, &
         pub_do_tddft, pub_md_properties, pub_hubbard, &
         pub_hubbard_restart, pub_hubbard_atomsolve, &
         pub_hub_max_iter, pub_cond_calculate, &
         pub_ii_energy_direct, pub_is_implicit_solvent, &
         pub_is_smeared_ion_rep, pub_mt_cutoff, &
         pub_output_detail, pub_fine_grid_scale, pub_dbl_grid_scale, pub_hfxsw,&
         pub_ovlp_for_nonlocal, pub_paw, pub_spin, pub_spin_polarised, &
         pub_write_forces, task, pub_nlcc, pub_devel_code, &
         pub_fine_is_dbl, pub_dbl_is_std, maxit_ngwf_cg, maxit_ngwf_diis, &
         pub_nonsc_forces, pub_aug_den_dim, pub_aug, pub_realspace_projectors, &
         pub_initial_dens_realspace, pub_constant_efield, pub_cdft, &
         pub_ci_cdft, pub_use_aux_ngwfs
!CW
    use rundat, only: pub_dmft_points
!END CW
    use services, only: services_flush, services_print_num_species
    use simulation_cell, only: pub_cell, simulation_cell_fftbox_init, &
         simulation_cell_fftbox_exit
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_copy, &
         sparse_destroy, sparse_exit, sparse_mod_init, &
         sparse_init_blocking_scheme, sparse_scale, &
         BLKS_NGWF, BLKS_PROJ, BLKS_HUB_PROJ, BLKS_COND, BLKS_JOINT, BLKS_AUX
    use ngwfs, only: ngwfs_initialise
    use tddft, only: tddft_calculate
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    use vdwcorrection, only: vdwcorrection_calculate_energy,&
         pub_dispersion_energy
    use xc, only: xc_hfxinit, xc_init, xc_exit, pub_xc_gradient_corrected

    implicit none

    ! Arguments
    type(ELEMENT), intent(inout)   :: elements(pub_cell%nat)
    real(kind=DP), intent(out)  :: total_energy
    real(kind=DP), intent(out)  :: total_forces(1:3,1:pub_cell%nat)
    logical,intent(in),optional :: properties_only
    logical,intent(out),optional :: return_converged

    ! Local Variables
    type(FUNC_BASIS) :: ngwf_basis
    type(FUNC_BASIS) :: aux_ngwf_basis
    type(FUNC_BASIS) :: cond_ngwf_basis
    type(FUNC_BASIS) :: joint_ngwf_basis
    type(FUNC_BASIS) :: proj_basis
    type(FUNC_BASIS) :: hub_proj_basis
    type(NGWF_REP) :: rep
    type(NGWF_HAM) :: ham
    type(HUBBARD_MODEL) :: hub
    type(SPAM3), allocatable, dimension(:) :: denskern
    real(kind=DP), allocatable, dimension(:,:) :: ngwf_nonsc_forces
    real(kind=DP), dimension(:,:,:), allocatable :: localpseudo_fine
    real(kind=DP), dimension(:,:,:), allocatable :: core_density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: lhxc_fine
    real(kind=DP) :: ewald_energy, throwaway
    real(kind=DP) :: proj_sphere_padding
    integer :: n1, n2, n3
    integer :: is
    integer :: ierr
    integer :: local_size, global_size
    logical :: local_properties
    logical :: ngwfs_converged
    logical :: cdft_converged !gibo: logical for cDFT optimisation
    integer :: hub_proj_iteration ! ddor: The DFT+U projector-consistency
                                  ! ddor: iteration we are on
!CW
    logical :: check
!END CW

    ! cks: this initialises the PSEUDO_SPECIES type array
    ! cks: and sets pub_cell%num_projectors
    ! aam: comment: and also sets much of elements
    ! aam: 01/12/2004 - moved here from onetep.F90
    if (.not.pub_paw) then
       call pseudopotentials_read_species(elements) !input/output
    end if

    if (pub_paw) then
       call paw_read_species(elements)
    end if

    ! Print atom counting information
    if (pub_on_root) call services_print_num_species(elements)

    call services_flush

    ! set local variable for properties_only
    local_properties = .false.
    if (present(properties_only)) local_properties = properties_only
    ngwfs_converged = .true. ! set to true if no optimisation requested

    ! Adopt parallel strategy
    call comms_barrier
    if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
         write(stdout,'(a)',advance='no') 'Determining parallel strategy ...'
    call parallel_strategy_distr_atoms(elements)

    ! Distribute the NGWFs
    call function_basis_allocate(ngwf_basis,pub_cell%num_ngwfs,'ngwfs')
    call function_basis_distribute(ngwf_basis,elements)

    ! Exit if there is a node with no NGWFs
    if (minval(ngwf_basis%num_on_node) < 1) then
       if (pub_on_root) then
           write(stdout,*) ' '
           write(stdout,'(a)') 'Too many processors for this system size:'
           write(stdout,'(a)') 'Either reduce the number of processors or &
                &increase the number of atoms'
       end if
       call comms_abort
    end if

    ! Distribute the nonlocal projectors if required
    if (pub_any_nl_proj) then
       call function_basis_allocate(proj_basis,pub_cell%num_projectors, &
            'projs')
       call function_basis_distribute(proj_basis,elements)
    end if

    ! Distribute PAW partial waves if required
    if (pub_paw) then
       call function_basis_allocate(proj_basis,pub_cell%num_pawpws,'pawpws')
       call function_basis_distribute(proj_basis,elements)
    end if

    ! Distribute Hubbard projectors if required
    if (pub_hubbard) then
       call function_basis_allocate(hub_proj_basis,pub_cell%num_hub_proj, &
            'hub_projs')
       call function_basis_distribute(hub_proj_basis,elements)
    end if

    ! Distribute Conduction NGWFs if required
    if (pub_cond_calculate) then
       call function_basis_allocate(cond_ngwf_basis,pub_cell%num_ngwfs_cond, &
            'ngwfs_cond')
       call function_basis_distribute(cond_ngwf_basis,elements)

       call function_basis_allocate(joint_ngwf_basis, pub_cell%num_ngwfs &
            + pub_cell%num_ngwfs_cond,'ngwfs_joint')
       call function_basis_distribute(joint_ngwf_basis,elements)
    end if

    ! Distribute Auxiliary NGWFs if required
    if (pub_use_aux_ngwfs) then
       call function_basis_allocate(aux_ngwf_basis,pub_cell%num_ngwfs_aux, &
            'ngwfs_aux')
       call function_basis_distribute(aux_ngwf_basis,elements)
    end if

    ! Report results of distribution
    if (pub_on_root .and. pub_output_detail == VERBOSE) then
       write(stdout,'(3(a,i6))') '... NGWF load balancing: max ', &
            maxval(ngwf_basis%num_on_node),';  min ', &
            minval(ngwf_basis%num_on_node),';  average ', &
            nint(ngwf_basis%num/real(pub_total_num_nodes,kind=DP))
       if (pub_any_nl_proj.or.pub_paw) then
          write(stdout,'(3(a,i6))') '... Projector load balancing: max ', &
               maxval(proj_basis%num_on_node),';  min ', &
               minval(proj_basis%num_on_node),';  average ', &
               nint(proj_basis%num/real(pub_total_num_nodes,kind=DP))
       end if
       if (pub_hubbard) then
          write(stdout,'(3(a,i6))') '... Hubbard Projector load balancing: max ', &
               maxval(hub_proj_basis%num_on_node),';  min ', &
               minval(hub_proj_basis%num_on_node),';  average ', &
               nint(hub_proj_basis%num/real(pub_total_num_nodes,kind=DP))
       end if
       if (pub_cond_calculate) then
          write(stdout,'(3(a,i6))') '... Conduction NGWF load balancing: max ', &
               maxval(cond_ngwf_basis%num_on_node),';  min ', &
               minval(cond_ngwf_basis%num_on_node),';  average ', &
               nint(cond_ngwf_basis%num/real(pub_total_num_nodes,kind=DP))
          write(stdout,'(3(a,i6))') '... Joint NGWF load balancing: max ', &
               maxval(joint_ngwf_basis%num_on_node),';  min ', &
               minval(joint_ngwf_basis%num_on_node),';  average ', &
               nint(joint_ngwf_basis%num/real(pub_total_num_nodes,kind=DP))
       end if
       if (pub_use_aux_ngwfs) then
          write(stdout,'(3(a,i6))') '... Auxiliary NGWF load balancing: max ', &
               maxval(aux_ngwf_basis%num_on_node),';  min ', &
               minval(aux_ngwf_basis%num_on_node),';  average ', &
               nint(aux_ngwf_basis%num/real(pub_total_num_nodes,kind=DP))
       end if
    end if

    call comms_barrier
    if (pub_on_root .and. (pub_output_detail >= NORMAL)) &
         write(stdout,'(a)') '... done'
    call services_flush

    ! pdh: count number of electrons and determine spin polarisation strategy
    call internal_electrons

!CW
 if(pub_dmft_points>0)then
     INQUIRE(file='store_ham1',EXIST=check)
     if(check) goto 78
 endif
!END CW

    ! Calculate the Ion-Ion energy
    ! jd: Direct summation requested instead of Ewald
    if (pub_ii_energy_direct) then
       if (pub_on_root) then
          write(stdout,'(/a)',advance='no') &
               'Calculating direct Ion-Ion energy ...'
       end if

       if (pub_coulomb_cutoff) then

          ! cks: Different approach depending on whether there are classical atoms
          if (pub_cell%nat_classical > 0) then

             select case (pub_coulomb_cutoff_type)
                ! 0D periodicity
             case('SPHERE','sphere','CYLINDER','cylinder')
                call classical_pot_ii_energy(ewald_energy, & ! output
                     elements)
             case default
                call utils_abort('Illegal cutoff type for cutoff Coulomb with embedding')
             end select

          else
             ! jd: Usual cutoff Coulomb calculation
             call cutoff_coulomb_ii_energy(elements, ewald_energy, &
                  pub_coulomb_cutoff_type)
          endif

       else
          ! jd: Reuse the CC code for direct summation for other calculations
          !     that use direct summation, like MT or smeared-ions
          call cutoff_coulomb_ii_energy(elements, ewald_energy, 'SPHERE')
       end if
    ! jd: Ewald by default
    else
       if (pub_on_root) then
          write(stdout,'(/a)',advance='no') &
               'Calculating Ewald energy ...'
       end if

       if (pub_cell%nat_classical > 0) then
          ! cks: include classical atom charges in Ewald energy
          call ewald_calculate_energy(elements, ewald_energy,.false., &
               classical_elements)
          call ewald_calculate_energy(elements, throwaway, .true.)
       else
          call ewald_calculate_energy(elements, ewald_energy, .false.)
          call ewald_calculate_energy(elements, throwaway, .true.)
       endif

    end if

    ! ndmh: calculate energy of ions in external field
    if (any(abs(pub_constant_efield) > 0.0_DP)) &
         call potential_add_efield_ion_energy(ewald_energy)

    call services_flush

    ! qoh: Dispersion energy
    if (pub_dispersion /=0) then
       call vdwcorrection_calculate_energy(elements)
    else
       pub_dispersion_energy = 0.0_DP
    end if

!CW
 78 continue
!END CW

    ! qoh: Initialise pub_usehfx and pub_hfxfraction this needs to be done before
    ! qoh: pub_fftbox is initialised.
    call xc_hfxinit

    ! cks: FFTBOX dimensions in gridpoints
    call ppd_strategy_fftbox_spec(n1, n2, n3, &     ! output
         elements)                                  ! input

    ! cks: initialise PUB_FFTBOX
    call simulation_cell_fftbox_init(n1, n2, n3)
    call services_flush
    call hf_exchange_set_pub_hfxsw ! Depends on pub_fftbox

    ! PPDS, SPHERES, BOXES, KB_DENOMINATORS
    if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(/a)',advance='no') 'Basis initialisation ...'

    ! CHECK CENTRES
    ! cks: make sure that no atom is outside the simulation cell
    call parallel_strategy_check_atoms(elements)
    if (pub_on_root .and. pub_output_detail == VERBOSE) &
         write(stdout,'(/a)') '... Atom positions lie inside simulation cell'

    ! jd: make sure that the cell is large enough if MT correction is in effect
    if (pub_on_root .and. pub_mt_cutoff /= 0.0_DP) &
         call pbc_corr_check_cell_size(elements)

    ! NGWF SPHERES
    call function_basis_init_spheres(ngwf_basis, pub_elements_on_node)
    if (pub_on_root .and. pub_output_detail == VERBOSE) &
         write(stdout,'(a)') '... NGWF spheres initialised'

    ! PROJECTOR SPHERES
    if (pub_any_nl_proj.or.pub_paw) then
       if (.not.pub_realspace_projectors) then
          call function_basis_init_spheres(proj_basis, pub_elements_on_node)
       else
          proj_sphere_padding = 2.0_DP*maxval(ngwf_basis%spheres(:)%radius)
          call comms_reduce('MAX',proj_sphere_padding)
          call function_basis_init_spheres(proj_basis, pub_elements_on_node, &
               proj_sphere_padding)
       end if
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') '... Projector spheres initialised'
    end if

    ! HUBBARD PROJECTOR SPHERES
    if (pub_hubbard) then
       call function_basis_init_spheres(hub_proj_basis, pub_elements_on_node)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') '... Hubbard Projector spheres initialised'
       ! ndmh: Hubbard Model allocation and setup
       call hubbard_model_init(hub,elements)
    endif

    ! CONDUCTION NGWF SPHERES
    if (pub_cond_calculate) then
       call function_basis_init_spheres(cond_ngwf_basis, pub_elements_on_node)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') '... Conduction NGWF spheres initialised'
       call function_basis_copy_spheres(joint_ngwf_basis,ngwf_basis, &
            cond_ngwf_basis)
    end if

    ! AUXILIARY NGWF SPHERES
    if (pub_use_aux_ngwfs) then
       call function_basis_init_spheres(aux_ngwf_basis, pub_elements_on_node)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') '... Auxiliary NGWF spheres initialised'
    end if

    ! cks: initialisation of tight_boxes
    call function_basis_init_tight_boxes(ngwf_basis)
    if (pub_on_root .and. pub_output_detail == VERBOSE) &
         write(stdout,'(a)') '... NGWF tight boxes initialised'

    ! ndmh: initialisation of tight_boxes for PAW spheres
    if (pub_paw.or.pub_realspace_projectors) then
       call function_basis_init_tight_boxes(proj_basis)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') '... Projector tight boxes initialised'
    end if

    ! ndmh: initialisation of conduction NGWF tight_boxes
    if (pub_cond_calculate) then
       call function_basis_init_tight_boxes(cond_ngwf_basis)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') '... Conduction NGWF tight boxes initialised'
       call function_basis_init_tight_boxes(joint_ngwf_basis)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') '... Joint NGWF tight boxes initialised'
    end if

    ! ndmh: initialisation of auxiliary NGWF tight_boxes
    if (pub_use_aux_ngwfs) then
       call function_basis_init_tight_boxes(aux_ngwf_basis)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') '... Auxiliary NGWF tight boxes initialised'
    end if

    call comms_barrier

    ! cks: gather all tightboxes from nodes and join them in a large array
    call function_basis_gath_all_tbs(ngwf_basis)

    ! lr408: Gather all conduction tightboxes from nodes and join them in a large array
    ! lr408: Pass both NGWF sets to ensure universal tightbox has correct size
    if (pub_cond_calculate) then
       call function_basis_gath_all_tbs(cond_ngwf_basis)
       call function_basis_gath_all_tbs(joint_ngwf_basis)
       call function_basis_init_uni_tb(ngwf_basis,cond_ngwf_basis)
    else
       call function_basis_init_uni_tb(ngwf_basis)
    end if

    if (pub_use_aux_ngwfs) call function_basis_gath_all_tbs(aux_ngwf_basis)

    call function_basis_estimate_size(ngwf_basis,local_size,global_size)
    if (pub_on_root .and. pub_output_detail == VERBOSE) write(stdout,'(a,2i10)') &
         '... NGWF basis size (local,global): ', local_size, global_size

    ! lr408: Estimate conduction basis size as well
    if (pub_cond_calculate) then
       call function_basis_estimate_size(cond_ngwf_basis,local_size,global_size)
       if (pub_on_root .and. pub_output_detail == VERBOSE) write(stdout,'(a,2i10)') &
            '... Conduction NGWF basis size (local,global): ', local_size, global_size
       call function_basis_estimate_size(joint_ngwf_basis,local_size,global_size)
       if (pub_on_root .and. pub_output_detail == VERBOSE) write(stdout,'(a,2i10)') &
            '... Joint NGWF basis size (local,global): ', local_size, global_size
    end if

    call comms_barrier
    if (pub_on_root .and. pub_output_detail == VERBOSE) write(stdout,'(a)') &
         '... Tight-boxes gathered'

    ! cks: initialise nonlocal pseudpotential projectors in
    ! cks: complex fftbox reciprocal representation
    if (pub_any_nl_proj) then
       call pseudo_species_init_proj(nlps_projectors,elements)
    end if
    if (pub_paw) then
       call paw_species_init_proj(paw_projectors,elements)
    end if

    call comms_barrier
    if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)') '... done'
    call services_flush

    ! cks: print psinc grid info and check fftbox and ppd sizes
    call ppd_strategy_check_and_print

    ! ndmh: initialise standard grid
    call cell_grid_distribute(pub_std_grid,pub_cell,1.0_DP,'Coarse grid',.true.)
    call fourier_init_cell(pub_std_grid)

    ! ndmh: see if dbl grid must be initialised separately
    if (pub_dbl_grid_scale > 1.0_DP) then

       ! ndmh: initialise dbl grid
       ! jd: pass .true. as the last param if double is fine,
       !     to (potentially) tweak the distribution for multigrid
       call cell_grid_distribute(pub_dbl_grid,pub_cell,&
            pub_dbl_grid_scale,'Double grid',.true., &
            (pub_fine_grid_scale==pub_dbl_grid_scale))
       call fourier_init_cell(pub_dbl_grid)
       pub_dbl_is_std = .false.

    else

       ! ndmh: copy dbl and fine grids from std grid
       pub_dbl_grid = pub_std_grid
       pub_dbl_is_std = .true.

    end if

    ! ndmh: see if fine grid must be initialised separately
    if (pub_fine_grid_scale/=pub_dbl_grid_scale) then

       ! ndmh: initialise fine grid
       call cell_grid_distribute(pub_fine_grid,pub_cell, &
            pub_fine_grid_scale,'  Fine grid',.true.,.true.)
       call fourier_init_cell(pub_fine_grid)
       pub_fine_is_dbl = .false.

    else

       ! ndmh: copy dbl and fine grids from std grid
       pub_fine_grid = pub_dbl_grid
       pub_fine_is_dbl = .true.

    end if

    ! ndmh: set up parallel strategy and Fourier transforms for padded cell
    if (pub_coulomb_cutoff) call cutoff_coulomb_init(elements)

    if (pub_aug) then
       ! ndmh: Initialise PAW augmentation density box
       call augmentation_box_init(pub_fine_grid)
    end if

    ! ndmh: Initialise FFTbox, tightbox and aug box Fourier transforms
    call fourier_init
    call services_flush

    ! ndmh: Initialise xc module
    call xc_init

    if (pub_aug) then
       ! ndmh: set up num of components of PAW compensation density & gradient
       pub_aug_den_dim = 0
       if (pub_xc_gradient_corrected) pub_aug_den_dim = 3
    end if

    ! Create projectors in real space if required
    if (pub_any_nl_proj.and.pub_realspace_projectors) &
         call projectors_create_real(proj_basis,nlps_projectors)
    if (pub_paw.and.pub_realspace_projectors) &
         call projectors_create_real(proj_basis,paw_projectors)

    if (pub_hubbard) then
       ! ddor: initialise hubbard projectors in
       !      complex fftbox reciprocal representation
       call hubbard_species_proj(hub,elements)

      if (task == 'HUBBARDSCF' .or. pub_hubbard_restart &
            & .or. pub_hubbard_atomsolve) then
          call function_basis_init_tight_boxes(hub_proj_basis,ngwf_basis)
          call function_basis_gath_all_tbs(hub_proj_basis,ngwf_basis)
          call hubbard_build_consist(hub,hub_proj_basis,ngwf_basis)
       endif
    endif

    ! qoh: This subroutine calls sparse_init_blocking_scheme for spherical waves
    if (pub_hfxsw) call hf_exchange_init_vmatrix(elements)

    ! ndmh: initialise dense matrix module
    call dense_init

    ! ndmh: initialise SPAM3 module
    call sparse_init_blocking_scheme(BLKS_NGWF,ngwf_basis%num, &
         ngwf_basis%num_on_node, ngwf_basis%num_on_atom, &
         ngwf_basis%first_on_node, ngwf_basis%first_on_atom, &
         ngwf_basis%atom_of_func, ngwf_basis%node_of_func)
    if (pub_any_nl_proj.or.pub_paw) then
       call sparse_init_blocking_scheme(BLKS_PROJ,proj_basis%num, &
            proj_basis%num_on_node, proj_basis%num_on_atom, &
            proj_basis%first_on_node, proj_basis%first_on_atom, &
            proj_basis%atom_of_func, proj_basis%node_of_func)
    end if
    if (pub_hubbard) then
       call sparse_init_blocking_scheme(BLKS_HUB_PROJ,hub_proj_basis%num, &
            hub_proj_basis%num_on_node, hub_proj_basis%num_on_atom, &
            hub_proj_basis%first_on_node, hub_proj_basis%first_on_atom, &
            hub_proj_basis%atom_of_func, hub_proj_basis%node_of_func)
    end if
    if (pub_cond_calculate) then
       call sparse_init_blocking_scheme(BLKS_COND,cond_ngwf_basis%num, &
            cond_ngwf_basis%num_on_node, cond_ngwf_basis%num_on_atom, &
            cond_ngwf_basis%first_on_node, cond_ngwf_basis%first_on_atom, &
            cond_ngwf_basis%atom_of_func, cond_ngwf_basis%node_of_func)
       call sparse_init_blocking_scheme(BLKS_JOINT,joint_ngwf_basis%num, &
            joint_ngwf_basis%num_on_node, joint_ngwf_basis%num_on_atom, &
            joint_ngwf_basis%first_on_node, joint_ngwf_basis%first_on_atom, &
            joint_ngwf_basis%atom_of_func, joint_ngwf_basis%node_of_func)
    end if
    if (pub_use_aux_ngwfs) then
       call sparse_init_blocking_scheme(BLKS_AUX,aux_ngwf_basis%num, &
            aux_ngwf_basis%num_on_node, aux_ngwf_basis%num_on_atom, &
            aux_ngwf_basis%first_on_node, aux_ngwf_basis%first_on_atom, &
            aux_ngwf_basis%atom_of_func, aux_ngwf_basis%node_of_func)
    end if

    call sparse_mod_init
    call services_flush

    ! ndmh: Allocate SPAM3 structures for overlap matrix, kinetic energy,
    ! ndmh: density kernel, inverse overlap
    call ngwf_rep_create(rep,'',elements)

    allocate(denskern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('energy_and_force_calculate','denskern',ierr)
    do is=1,pub_cell%num_spins
       denskern(is)%structure = 'K'
       call sparse_create(denskern(is))
    end do
    pub_kernel_workspace_allocated = .false.

    ! ddor: initialise SPAM3 matrices for DFT+U if necessary
    if (pub_hubbard) then
       call hubbard_build_matrices(hub,hub_proj_basis)
    endif

    if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)') ' done'

    ! ndmh: allocate storage for local pseudopotential on fine grid
    allocate(localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('energy_and_force_calculate','localpseudo_fine', &
         ierr)

    ! ndmh: allocate storage for core charge if NLCC is present
    if (pub_nlcc) then
       allocate(core_density_fine(pub_fine_grid%ld1, &
            pub_fine_grid%ld2, pub_fine_grid%max_slabs12),stat=ierr)
       call utils_alloc_check('energy_and_force_calculate', &
            'core_density_fine',ierr)
    else
       ! ndmh: no core charge required - allocate dummy storage
       allocate(core_density_fine(1,1,1),stat=ierr)
       call utils_alloc_check('energy_and_force_calculate', &
            'core_density_fine',ierr)
    end if

    ! lr408: Allocate lhxc_fine here so we can avoid unnecessary recalculation
    allocate(lhxc_fine(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('energy_and_force_calculate','lhxc_fine',ierr)

    ! jd: Initialise implicit solvation module, if required
    if (pub_is_implicit_solvent .or. pub_is_smeared_ion_rep) then
       call multigrid_initialise
    end if

    ! ndmh: Allocate and initialise NGWFs
    allocate(rep%ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('energy_and_force_calculate','rep%ngwfs_on_grid', &
         ierr)

    ! ndmh: Moved Hamiltonian creation here so same ham can be re-used for
    ! ndmh: kernel initialisation
    if (.not. pub_cond_calculate) then
       call ngwf_ham_create(ham,rep)
       ! qoh: Initialise vmatrix
       if (pub_hfxsw) call hf_exchange_fill_vmatrix(elements, ham%full_vmatrix,&
            ngwf_basis)
    end if

    ! ndmh: Allocate storage for radial densities
    if (pub_initial_dens_realspace) call density_radial_init

    call ngwfs_initialise(rep%ngwfs_on_grid,ngwf_basis,elements)
    call comms_barrier
    call services_flush

    ! cks: initialise density kernel, local pseudopotential on fine grid
    ! cks: and density-independent matrices
!CW
 if(pub_dmft_points>0)then
     INQUIRE(file='store_ham1',EXIST=check)
     if(check) goto 77
 endif
!END CW

    call electronic_init_denskern(rep, ham, denskern, localpseudo_fine, &
         core_density_fine, lhxc_fine, ngwf_basis, proj_basis, hub_proj_basis, &
         hub, elements, ewald_energy)

!CW
    77 continue
!END CW

    ! ndmh: Deallocate storage for radial densities
    if (pub_initial_dens_realspace) call density_radial_exit

    call comms_barrier

    ! ars: allocate NGWF non self-consistent forces only if required
    if (pub_nonsc_forces) then
       allocate(ngwf_nonsc_forces(1:3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('energy_and_force_calculate','ngwf_nonsc_forces', &
            ierr)
       ngwf_nonsc_forces(:,:) = -999_DP ! ndmh: init to flagged value
    else
       allocate(ngwf_nonsc_forces(0,0),stat=ierr)
       call utils_alloc_check('energy_and_force_calculate','ngwf_nonsc_forces', &
         ierr)
    end if

    ! ddor: Build DFT+U projectors in the case where a restart
    !       file is used with a task such as PROPERTIES
    if ((pub_hubbard_restart .or. pub_hubbard_atomsolve) &
         & .and. (task .ne. 'HUBBARDSCF')) then
       call hubbard_projector_consistency( &
            elements, rep%ngwfs_on_grid, ngwf_basis, &
            proj_basis, hub_proj_basis, hub, rep%overlap, rep%inv_overlap, &
            rep%hub_overlap, rep%hub_overlap_t, &
            rep%sp_overlap, rep%hub_proj_paw_overlap, &
            rep%hub_ngwf_paw_overlap)
    endif

    ! ndmh: standard valence NGWF optimisation (used for most tasks)
    if ((task/='PROPERTIES').and.(task/='COND').and.(task/='HUBBARDSCF').and. &
         (task/='TDDFT').and.(task/='PROPERTIES_COND').and.(.not.local_properties)) then
       ! optimise density kernel and NGWF coefficients
       if ((maxit_ngwf_cg .gt. 0).or.(maxit_ngwf_diis==0).or.pub_nonsc_forces) then
!          ! gibo: unless cDFT is activated proceed as usual,
!          ! gibo: otherwise trigger constraint-optimisation loop
          if ((.not.pub_cdft) .OR. (pub_ci_cdft)) then
             call ngwf_cg_optimise(total_energy, ngwfs_converged, ham, denskern, &
                  rep, ngwf_nonsc_forces, lhxc_fine, ngwf_basis, proj_basis, &
                  hub_proj_basis, hub, elements, ewald_energy, &
                  localpseudo_fine, core_density_fine)
          else
             ! pass the same info passed to subroutine ngwf_cg_optimise()
              call cdft_u_cg_optimise(total_energy, ngwfs_converged, &
                   cdft_converged, ham, denskern, &
                   rep, ngwf_nonsc_forces, lhxc_fine, ngwf_basis, proj_basis, &
                   hub_proj_basis, hub, elements, ewald_energy, &
                   localpseudo_fine, core_density_fine)
          end if

       end if
       if (maxit_ngwf_diis .gt. 0) then
!          ! gibo: unless cDFT is activated proceed as usual,
!          ! gibo: otherwise trigger constraint-optimisation loop
          if ((.not.pub_cdft) .OR. (pub_ci_cdft)) then
              call ngwf_diis_optimise(total_energy, ngwfs_converged, ham, &
                   denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, hub, &
                   elements, ewald_energy, localpseudo_fine, core_density_fine, &
                   lhxc_fine)
          else
             ! pass the same info passed to subroutine ngwf_diis_optimise()
              call cdft_u_cg_optimise(total_energy, ngwfs_converged, &
                   cdft_converged, ham, denskern, &
                   rep, ngwf_nonsc_forces, lhxc_fine, ngwf_basis, proj_basis, &
                   hub_proj_basis, hub, elements, ewald_energy, &
                   localpseudo_fine, core_density_fine)
          endif
       end if
    end if

    if (task == 'HUBBARDSCF') then
       ! ddor: DFT+U energy minimisation with self-consistent NGWF projectors.
       ! Loop over Hubbard projector optimisation counter
       projectoroptimisation : do hub_proj_iteration=1,pub_hub_max_iter
          if (pub_on_root) then
             write(stdout,'(/a/)')'########################################&
                  &########################################'
             write(stdout,'(a,i5)') 'HUBBARD PROJECTOR SCF ITERATION ', &
                  hub_proj_iteration
             write(stdout,'(/a/)')'########################################&
                  &########################################'
          endif
          hub%consistency_iteration = hub_proj_iteration ! public
          ! ddor: Depending on the projector iteration we are on,
          !       carry out some projector and metric mixing,
          !       and generate the new projector-NGWF overlap matrix
          call hubbard_projector_consistency( &
               elements, rep%ngwfs_on_grid, ngwf_basis, &
               proj_basis, hub_proj_basis, hub, rep%overlap, rep%inv_overlap, &
               rep%hub_overlap, rep%hub_overlap_t, rep%sp_overlap, &
               rep%hub_proj_paw_overlap, rep%hub_ngwf_paw_overlap)

          ! gibo: unless CDFT is activated proceed as usual,
          ! gibo: otherwise trigger constraint-optimisation loop
          if ((.not.pub_cdft).OR. (pub_ci_cdft)) then

             ! optimise density kernel and NGWF coefficients
             call ngwf_cg_optimise(total_energy, ngwfs_converged, ham, &
                  denskern, rep, ngwf_nonsc_forces, lhxc_fine, &
                  ngwf_basis, proj_basis, hub_proj_basis, hub, &
                  elements, ewald_energy, localpseudo_fine, core_density_fine)

          else

              call cdft_u_cg_optimise(total_energy, ngwfs_converged, &
                  cdft_converged, ham, &
                  denskern, rep, ngwf_nonsc_forces, lhxc_fine, &
                  ngwf_basis, proj_basis, hub_proj_basis, hub, &
                  elements, ewald_energy, localpseudo_fine, core_density_fine)

          endif

          ! ddor: Call to Hubbard projector optimisation
          !       energy convergence criterion
          if (hubbard_test_convergence(hub,total_energy)) &
               exit projectoroptimisation
       enddo projectoroptimisation
    end if

    if ((task == 'COND').or.(task == 'PROPERTIES_COND')) then
       ! lr408: Call conduction driver routine
       call conduction_ngwf_optimise(total_energy, denskern, rep, ngwf_basis, &
            proj_basis, hub_proj_basis, hub, cond_ngwf_basis, joint_ngwf_basis, &
            elements, ewald_energy, localpseudo_fine, core_density_fine, &
            lhxc_fine,cond_properties_only=(task=='PROPERTIES_COND'))
    end if

    if (((task=='PROPERTIES').or.(local_properties).or.(pub_do_properties) &
         .or.(pub_md_properties)).and.(task/='COND').and. &
         (task/='PROPERTIES_COND')) then
       ! cks: calculate electronic properties
       call properties_calculate(denskern, ham, &
            rep, ngwf_basis, proj_basis, hub_proj_basis, hub, elements, &
            localpseudo_fine, core_density_fine, lhxc_fine, .true., &
            'valence')
    end if

    if ((task == 'TDDFT').or.(pub_do_tddft)) then
       ! ddor: calculate dielectric response from TDDFT
       call tddft_calculate(denskern, rep, ngwf_basis, proj_basis, elements, &
            localpseudo_fine, core_density_fine, ewald_energy)
    end if

    ! cks: calculate forces if required
    if (pub_write_forces .or. task == 'GEOMETRYOPTIMIZATION' .or. &
          task == 'TRANSITIONSTATESEARCH' .or. task == 'MOLECULARDYNAMICS' .or. &
          task == 'PHONON') then
       ! ndmh: Calculate ionic forces
       call forces_calculate(total_forces,denskern,ham,lhxc_fine,rep, &
            ngwf_basis,proj_basis,hub_proj_basis,hub,localpseudo_fine, &
            core_density_fine,elements,ngwf_nonsc_forces)
    end if

    ! Deallocate

    ! ndmh: deallocate storage for NGWFs
    deallocate(rep%ngwfs_on_grid,stat=ierr)
    call utils_dealloc_check('energy_and_force_calculate','rep%ngwfs_on_grid', &
         ierr)

    if (.not. pub_cond_calculate) call ngwf_ham_destroy(ham)

    deallocate(ngwf_nonsc_forces,stat=ierr)
    call utils_dealloc_check('energy_and_force_calculate','ngwf_nonsc_forces', &
         ierr)

    deallocate(lhxc_fine,stat=ierr)
    call utils_dealloc_check('energy_and_force_calculate','lhxc_fine',ierr)

    deallocate(core_density_fine,stat=ierr)
    call utils_dealloc_check('energy_and_force_calculate','core_density_fine', &
         ierr)
    deallocate(localpseudo_fine,stat=ierr)
    call utils_dealloc_check('energy_and_force_calculate','localpseudo_fine', &
         ierr)

    ! cks: memory deallocation for sparse matrices
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(denskern(is))
    end do
    deallocate(denskern,stat=ierr)
    call utils_dealloc_check('energy_and_force_calculate','denskern',ierr)

    call ngwf_rep_destroy(rep)

    ! ddor: DFT+U deallocation if necessary
    if (pub_hubbard) then
       call hubbard_build_matrices_exit(hub)
       ! ddor: deallocate DFT+U projector information
       call hubbard_species_exit_proj(hub)
       ! ddor: deallocate hubbard_atoms type array in DFT+U
       call hubbard_model_exit(hub)
    endif

    ! pdh: clean up sparse module
    call sparse_exit

    ! ndmh: clean up dense module
    call dense_exit

    ! ndmh: clean up cutoff coulomb module
    if (pub_coulomb_cutoff) call cutoff_coulomb_exit

    ! jd: clean up pbc_correction module
    call pbc_corr_exit

    ! jd: clean up is_solvation module
    if (pub_is_implicit_solvent) call implicit_solvent_exit

    ! clean up XC routines
    call xc_exit

    ! clean up FFT routines
    call fourier_exit

    ! clean up whole-cell FFT routines
    call fourier_exit_cell

    ! clean up cell grids
    if (.not.pub_fine_is_dbl) call cell_grid_exit(pub_fine_grid)
    if (.not.pub_dbl_is_std) call cell_grid_exit(pub_dbl_grid)
    call cell_grid_exit(pub_std_grid)

    ! ndmh: Deallocate function basis objects
    if (pub_hubbard) call function_basis_deallocate(hub_proj_basis)
    if (pub_any_nl_proj.or.pub_paw) call function_basis_deallocate(proj_basis)
    call function_basis_deallocate(ngwf_basis)
    if (pub_use_aux_ngwfs) then
       call function_basis_deallocate(aux_ngwf_basis)
    end if
    if (pub_cond_calculate) then
       call function_basis_deallocate(cond_ngwf_basis)
       call function_basis_deallocate(joint_ngwf_basis)
    end if

    call function_basis_exit_uni_tb

    ! clean up parallel strategy
    call parallel_strategy_exit

    ! aam: deallocate pseudo_species type array in pseudopotentials
    if (.not.pub_paw) then
       call pseudopotentials_species_exit
    end if
    ! ndmh: Clean up arrays in the PAW module
    if (pub_paw) then
       call paw_species_exit
    end if

    ! vm: deallocate arrays in cell
    call simulation_cell_fftbox_exit

    ! ndmh: if return value indicating convergence was requested, set it now
    if (present(return_converged)) return_converged = ngwfs_converged

    return

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_electrons

      !======================================================================!
      ! This suboutine counts the number of electrons in the cell and uses   !
      ! this to check the spin polarisation strategy.                        !
      !----------------------------------------------------------------------!
      ! Written by Peter Haynes, July 2006, based on some code from          !
      ! internal_ppd_sph_box_kb_init by Chris-Kriton Skylaris.               !
      ! Modified by Phil Avraam, November 2008 for fractional electronic     !
      ! charges (tidied up by Nicholas Hine)                                 !
      !======================================================================!

      implicit none

      ! Local variables
      integer,parameter :: cfac=4
      real(kind=DP),parameter :: ctol=0.00001_DP
      integer :: atom
      integer :: nelecs
      integer :: nspins
      real(kind=DP) :: cfac_nelecs !, frac_nelecs ! ndmh: removed duplicate

      ! Count number of electrons for charge neutral cell
      cfac_nelecs = 0.0_DP
      do atom=1,pub_cell%nat
         ! pa: number of electrons times common factor
         ! pa: used to round the ionic charge to the nearest 1/cfac
         cfac_nelecs = cfac_nelecs + &
              real(nint(cfac*(elements(atom)%ion_charge)), kind=DP)
      end do
      ! pa: check that the total number of electrons is integer
      if (abs(cfac_nelecs/real(cfac,DP) - pub_charge - &
           nint(cfac_nelecs/real(cfac,DP) - pub_charge)) > ctol) then
         write(stdout,'(a/a)') &
             'WARNING in internal_electrons (energy_and_force_calculate):', &
             ' total electronic charge is non-integer! '
         call comms_abort
      endif

      ! cks: add the total charge of the system (from the input file)
      ! pa: allow for fractional ionic charge but not fractional total charge
      nelecs = nint(cfac_nelecs/cfac - pub_charge)

      ! pdh: check for inconsistencies between charge and spin
      ! pdh: and set number of occupied states accordingly
      if (mod(nelecs,2) /= 0) then
         if (.not. pub_spin_polarised) then
            if (pub_on_root) write(stdout,'(a/a)') &
                 'WARNING in internal_electrons (energy_and_force_calculate):',&
                 '        odd number of electrons in cell - &
                 &spin-polarised calculation will be performed'
            pub_spin_polarised = .true.
         end if
         if (pub_spin == 0) pub_spin = 1
         if (mod(pub_spin,2) == 0) then
            if (pub_on_root) write(stdout,'(a/a)') &
                 'WARNING in internal_electrons (energy_and_force_calculate):',&
                 '        odd number of electrons in cell - &
                 &even spin has been over-ridden to unity'
            pub_spin = 1
         end if
         rep%n_occ(UP) = (nelecs + pub_spin) / 2
         rep%n_occ(DN) = (nelecs - pub_spin) / 2
      else
         if (pub_spin_polarised .and. pub_spin == 0) then
            if (pub_on_root) write(stdout,'(a/a/a)') &
                 'WARNING in internal_electrons (energy_and_force_calculate):',&
                 '        even number of electrons in cell and zero spin', &
                 '        but continuing with spin-polarised calculation'
         end if
         if (mod(pub_spin,2) /= 0) then
            if (pub_on_root) write(stdout,'(a/a)') &
                 'WARNING in internal_electrons (energy_and_force_calculate):',&
                 '        even number of electrons in cell - &
                 &odd spin has been over-ridden to zero'
            pub_spin = 0
         end if
         rep%n_occ(UP) = (nelecs + pub_spin) / 2
         rep%n_occ(DN) = (nelecs - pub_spin) / 2
      end if

      ! pdh: set number of spins
      nspins = 1
      if (pub_spin_polarised) nspins = 2
      pub_cell%num_spins = nspins
      pub_cell%spin_fac = 2.0_DP / real(pub_cell%num_spins,kind=DP)

    end subroutine internal_electrons

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine energy_and_force_calculate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module energy_and_force




