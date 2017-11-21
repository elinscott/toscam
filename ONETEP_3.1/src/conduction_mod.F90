! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                                                                 !
!         Conduction States Calculation Module                    !
!                                                                 !
! This module optimises a set of NGWFs to describe conduction     !
! states and calculates the properties of those states.           !
!-----------------------------------------------------------------!
! Written by Laura Ratcliff on 14/10/2010                         !
! Minor improvements and cleanup by Nicholas Hine, October 2010.  !
!=================================================================!

module conduction

  use constants, only: DP

  implicit none

  private

  public :: conduction_ngwf_optimise

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine conduction_ngwf_optimise(total_energy, val_denskern, val_rep, &
            val_ngwf_basis, proj_basis, hub_proj_basis, hub, cond_ngwf_basis, &
            joint_ngwf_basis, elements, ewald_energy, localpseudo_fine, &
            core_density_fine, lhxc_fine, cond_properties_only)

    !========================================================================!
    ! This subroutine optimises a set of conduction NGWFs by forming a       !
    ! Hamiltonian which projects out valence contributions to the energy.    !
    ! It first initialises the appropriate valence variables and also calls  !
    ! properties_calculate at the end in which the various Hamiltonian forms !
    ! are diagonalised to give the individual eigenvalues.                   !
    !------------------------------------------------------------------------!
    ! total_energy    (out) : The final projected conduction energy          !
    ! val_denskern  (inout) : The ground state density kernel                !
    ! val_rep          (in) : Valence NGWF representation                    !
    ! val_ngwf_basis   (in) : Function basis type for the valence NGWFs      !
    ! proj_basis       (in) : Function basis type for nonlocal               !
    !                              pseudopotential projectors                !
    ! hub_proj_basis   (in) : Function basis type for Hubbard projectors     !
    ! cond_ngwf_basis  (in) : Function basis type for the conduction NGWFs   !
    ! joint_ngwf_basis (in) : Function basis type for the joint NGWFs        !
    ! elements         (in) : all elements in the order they appear in       !
    !                               the input file                           !
    ! ewald_energy     (in) : Ewald energy                                   !
    ! localpseudo_fine (inout) : local pseudopotential on fine grid          !
    ! core_density_fine(inout) : core density on fine grid for NLCC          !
    ! lhxc_fine        (inout) : Local-Hartree-Exhange-Correlation potential !
    !------------------------------------------------------------------------!
    ! Written by Laura Ratcliff in October 2010.                             !
    ! Modified calls to properties to enable the calculation of matrix       !
    ! elements for optical absorption spectra by Laura Ratcliff March 2011.  !
    !========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: comms_barrier, pub_on_root
    use conduction_properties, only: conduction_properties_calculate
    use constants, only: max_spins, stdout, paw_en_size, NORMAL
    use electronic, only: electronic_init_denskern
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_indep_matrices, &
         hamiltonian_dens_dep_matrices, hamiltonian_build_matrix
    use hubbard_build, only: HUBBARD_MODEL
    use ion, only: ELEMENT
    use kernel, only: pub_kernel_workspace_allocated
    use ngwf_cg, only: ngwf_cg_optimise
    use ngwfs, only: ngwfs_initialise
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_rep_create, &
         ngwf_rep_destroy, ngwf_ham_create, ngwf_ham_destroy
    use rundat, only: pub_homo_plot, pub_lumo_plot, &
         pub_homo_dens_plot, pub_lumo_dens_plot, cond_plot_joint_orbitals, &
         cond_plot_vc_orbitals, cond_num_extra_states, cond_num_extra_its, &
         maxit_ngwf_cg, pub_output_detail, pub_aug, pub_hubbard, &
         cond_read_denskern, cond_init_shift
    use services, only: services_flush
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale, &
         sparse_copy
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(FUNC_BASIS), intent(inout) :: cond_ngwf_basis
    type(FUNC_BASIS), intent(in) :: joint_ngwf_basis
    type(ELEMENT), intent(in)    :: elements(pub_cell%nat)
    type(SPAM3), intent(inout)   :: val_denskern(pub_cell%num_spins)
    type(NGWF_REP), intent(in)   :: val_rep
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(inout) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(inout) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(out) :: total_energy
    real(kind=DP), intent(in)  :: ewald_energy
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    logical, intent(in) :: cond_properties_only

    ! Local Variables
    type(NGWF_REP) :: cond_rep
    type(NGWF_HAM) :: cond_ham
    type(NGWF_HAM) :: val_ham
    type(SPAM3), allocatable, dimension(:) :: cond_denskern
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: hubbard_energy
    real(kind=DP) :: lhxc_energy
    integer :: is
    integer :: ierr
    integer :: tmp_maxit_ngwf_cg
    logical :: converged
    ! ars: ngwf_nonsc_forces
    real(kind=DP) :: ngwf_nonsc_forces(1:3,pub_cell%nat)


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &conduction_ngwf_optimise'
#endif

    ! ndmh: ngwf_basis%n_ppds is the number of ppds for the
    ! ndmh: ngwfs belonging to pub_my_node_id
    allocate(cond_rep%ngwfs_on_grid(cond_ngwf_basis%n_ppds*pub_cell%n_pts), &
         stat=ierr)
    call utils_alloc_check('conduction_ngwf_optimise','cond_rep%ngwfs_on_grid',&
         ierr)
    call comms_barrier

    ! ndmh: Allocate SPAM3 structures for overlap matrix, kinetic energy,
    ! ndmh: density kernel, inverse overlap
    call ngwf_rep_create(cond_rep,'c',elements)

    allocate(cond_denskern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('conduction_ngwf_optimise','cond_denskern',ierr)
    do is=1,pub_cell%num_spins
       cond_denskern(is)%structure = 'K'//cond_rep%postfix
       call sparse_create(cond_denskern(is))
    end do
    pub_kernel_workspace_allocated = .false.

    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)') ' done'

    ! lr408: Initialise conduction NGWFs to fireballs or read from file
    call ngwfs_initialise(cond_rep%ngwfs_on_grid,cond_ngwf_basis,elements, &
         read_cond=.true.)
    call comms_barrier
    call services_flush

    ! ndmh: create conduction hamiltonian
    call ngwf_ham_create(cond_ham,cond_rep)

    ! ndmh: create valence hamiltonian
    call ngwf_ham_create(val_ham,val_rep)

    ! lr408: scale for spin degeneracy
    if (pub_cell%num_spins == 1) call sparse_scale(val_denskern(1),2.0_DP)

    ! ndmh: calculate density dependent energies and matrices for valence NGWFS
    call hamiltonian_dens_dep_matrices(val_ham, lhxc_fine, total_energy, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, val_rep, &
         val_ngwf_basis, hub_proj_basis, hub, val_denskern, ewald_energy, &
         elements, pub_fine_grid, localpseudo_fine, core_density_fine, &
         ham_update=.true., lhxc_fixed=.false.)

    call hamiltonian_build_matrix(val_ham, val_rep)

    ! lr408: undo scaling
    if (pub_cell%num_spins == 1) call sparse_scale(val_denskern(1),0.5_DP)

    ! lr408: At some point might be nice to print out the final valence energy
    ! lr408: (components) here, but it's not essential to the calculation so
    ! lr408: for now we'll leave it

    if (.not. cond_properties_only) then

       if (pub_on_root) then
          write(stdout,'(a)') '+-----------------------------------------------&
               &-------------------------------+'
          write(stdout,'(a)') '|                   Starting conduction NGWF &
               &optimisation                      |'
          write(stdout,'(a)') '+-----------------------------------------------&
               &-------------------------------+'
       end if

       ! lr408: If required, do a few NGWF iterations optimising for a higher
       ! lr408: number of conduction states than required as a form of preconditioning
       if (cond_num_extra_states > 0) then
          if (pub_on_root) then
             write(stdout,'(a)') ''
             write(stdout,'(a)') 'Starting conduction NGWF pre-optimisation &
                  &process'
             write(stdout,'(a)') ''
          end if

          tmp_maxit_ngwf_cg = maxit_ngwf_cg
          maxit_ngwf_cg = cond_num_extra_its

          ! lr408: Set up 'occupation' numbers for conduction states
          call internal_cond_occupancies

          ! lr408: Initialise conduction kernel
          call electronic_init_denskern(cond_rep, cond_ham, cond_denskern, &
               localpseudo_fine, core_density_fine, lhxc_fine, &
               cond_ngwf_basis, proj_basis, hub_proj_basis, hub, elements, &
               ewald_energy, &
               val_rep, val_ngwf_basis, val_denskern, val_ham) ! optional args

          call comms_barrier

          ! ndmh: copy in dijhat from valence hamiltonian
          if (pub_aug) then
             do is=1,pub_cell%num_spins
                call sparse_copy(cond_ham%dijhat(is),val_ham%dijhat(is))
             end do
          end if

          ! lr408: Now we can optimise the conduction NGWFs
          call ngwf_cg_optimise(total_energy, converged, cond_ham, &
               cond_denskern, cond_rep, ngwf_nonsc_forces, lhxc_fine, &
               cond_ngwf_basis, proj_basis, hub_proj_basis, hub, &
               elements, ewald_energy, localpseudo_fine, core_density_fine, &
               val_rep, val_ngwf_basis, val_denskern, val_ham)

          maxit_ngwf_cg = tmp_maxit_ngwf_cg
          cond_num_extra_states = 0
          cond_read_denskern = .false.

          if (pub_on_root) then
             write(stdout,'(a)') ''
             write(stdout,'(a)') 'Conduction NGWF pre-optimisation process &
                  &complete.'
             write(stdout,'(a)') 'Proceeding with the main conduction NGWF &
                  &optimisation process.'
             write(stdout,'(a)') ''
          end if

          call comms_barrier

       end if

    end if ! not properties only

    ! lr408: Set up 'occupation' numbers for conduction states
    call internal_cond_occupancies

    ! lr408: Initialise conduction kernel
    ! ndmh: reset cond shift, as we are now re-using same cond_ham as before
    cond_ham%cond_shift = cond_init_shift
    call electronic_init_denskern(cond_rep, cond_ham, cond_denskern, &
         localpseudo_fine, core_density_fine, lhxc_fine, cond_ngwf_basis, &
         proj_basis, hub_proj_basis, hub, elements, ewald_energy, &
         val_rep, val_ngwf_basis, val_denskern, val_ham) ! optional args

    call comms_barrier

    ! ndmh: copy in dijhat from valence hamiltonian
    if (pub_aug) then
       do is=1,pub_cell%num_spins
          call sparse_copy(cond_ham%dijhat(is),val_ham%dijhat(is))
       end do
    end if

    if (.not. cond_properties_only) then

       ! lr408: Now we can optimise the conduction NGWFs
       call ngwf_cg_optimise(total_energy, converged, cond_ham, &
            cond_denskern, cond_rep, ngwf_nonsc_forces, lhxc_fine, &
            cond_ngwf_basis, proj_basis, hub_proj_basis, hub, &
            elements, ewald_energy, localpseudo_fine, core_density_fine, &
            val_rep, val_ngwf_basis, val_denskern, val_ham)


       if (pub_on_root) then
          write(stdout,'(a)') '+-----------------------------------------------&
               &-------------------------------+'
          write(stdout,'(a)') '|                   Conduction NGWF &
               &optimisation complete                      |'
          write(stdout,'(a)') '+-----------------------------------------------&
               &-------------------------------+'
       end if

    end if ! not properties only

    ! call conduction properties module
    call conduction_properties_calculate(val_denskern, cond_denskern, &
         val_rep, cond_rep, val_ham, cond_ham, &
         val_ngwf_basis, proj_basis, hub_proj_basis, hub, cond_ngwf_basis, &
         joint_ngwf_basis, elements, localpseudo_fine, &
         core_density_fine, lhxc_fine)


    ! ndmh: memory deallocation for cond_ham and val_ham
    call ngwf_ham_destroy(val_ham)
    call ngwf_ham_destroy(cond_ham)

    ! ndmh: memory deallocation for cond_rep
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(cond_denskern(is))
    end do
    deallocate(cond_denskern,stat=ierr)
    call utils_dealloc_check('conduction_ngwf_optimise','cond_denskern',ierr)
    call ngwf_rep_destroy(cond_rep)
    deallocate(cond_rep%ngwfs_on_grid,stat=ierr)
    call utils_dealloc_check('conduction_ngwf_optimise', &
         'cond_rep%ngwfs_on_grid',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &conduction_ngwf_optimise'
#endif


    return

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_cond_occupancies

      !======================================================================!
      ! This suboutine determines the number of conduction states of each    !
      ! spin.                                                                !
      !----------------------------------------------------------------------!
      ! Written by Laura Ratcliff, June 2010                                 !
      !======================================================================!

      use constants, only: DN, UP
      use rundat, only: cond_num_states, pub_spin

      implicit none

      ! lr408: Set to arbitrarily equal number of occupied states if not
      ! lr408: specified. Switches which has more if odd number of electrons
      ! lr408: and remains the same if even number of electrons

      if (cond_num_states == 0) then
         cond_rep%n_occ(UP) = val_rep%n_occ(DN)
         cond_rep%n_occ(DN) = val_rep%n_occ(UP)
      else if (cond_num_extra_states /= 0) then
         cond_rep%n_occ(UP) = (cond_num_states + cond_num_extra_states - pub_spin)/2
         cond_rep%n_occ(DN) = (cond_num_states + cond_num_extra_states + pub_spin)/2
      else
         cond_rep%n_occ(UP) = (cond_num_states - pub_spin)/2
         cond_rep%n_occ(DN) = (cond_num_states + pub_spin)/2
      end if

    end subroutine internal_cond_occupancies

  end subroutine conduction_ngwf_optimise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module conduction
