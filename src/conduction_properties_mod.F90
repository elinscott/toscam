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

module conduction_properties

  use constants, only: DP

  implicit none

  private

  public :: conduction_properties_calculate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine conduction_properties_calculate(val_denskern, cond_denskern, &
        val_rep, cond_rep, val_ham, cond_ham, &
        val_ngwf_basis, proj_basis, hub_proj_basis, hub, cond_ngwf_basis, &
        joint_ngwf_basis, elements, localpseudo_fine, &
        core_density_fine, lhxc_fine)

    !========================================================================!
    ! This subroutine calculates properties for a conduction calculation.    !
    !------------------------------------------------------------------------!
    ! Written by Laura Ratcliff in September 2011.                           !
    !========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: comms_barrier, pub_on_root
    use constants, only: stdout, NORMAL, max_spins
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_indep_matrices, &
         hamiltonian_dens_dep_nonsc
    use hubbard_build, only: HUBBARD_MODEL
    use ion, only: ELEMENT
    use kernel, only: pub_kernel_workspace_allocated
    use ngwfs, only: ngwfs_merge_sets
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_rep_create, &
         ngwf_rep_destroy, ngwf_ham_create, ngwf_ham_destroy
    use properties, only: properties_calculate
    use rundat, only: pub_homo_plot, pub_lumo_plot, &
         pub_homo_dens_plot, pub_lumo_dens_plot, cond_plot_joint_orbitals, &
         cond_plot_vc_orbitals, cond_num_extra_states, cond_num_extra_its, &
         pub_output_detail, pub_aug, pub_hubbard
    use services, only: services_flush
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_copy
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: val_ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(FUNC_BASIS), intent(in) :: cond_ngwf_basis
    type(FUNC_BASIS), intent(in) :: joint_ngwf_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(ELEMENT), intent(in)    :: elements(pub_cell%nat)
    type(SPAM3), intent(inout)   :: val_denskern(pub_cell%num_spins)
    type(SPAM3), intent(inout)   :: cond_denskern(pub_cell%num_spins)
    type(NGWF_REP), intent(in)   :: val_rep
    type(NGWF_REP), intent(inout)   :: cond_rep
    type(NGWF_HAM), intent(inout) :: val_ham
    type(NGWF_HAM), intent(inout) :: cond_ham
    real(kind=DP), intent(inout) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(inout) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)

    ! Local Variables
    type(NGWF_REP) :: joint_rep
    type(NGWF_HAM) :: joint_ham
    type(SPAM3), allocatable, dimension(:) :: joint_denskern
    integer :: n_occ_cond(max_spins)
    integer :: is
    integer :: ierr
    integer :: num_tot_eigen
    integer :: tmp_homo_plot, tmp_lumo_plot
    integer :: tmp_homo_dens_plot, tmp_lumo_dens_plot


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &conduction_properties'
#endif

    if (pub_on_root) then
       write(stdout,'(a)') '+--------------------------------------------------&
            &----------------------------+'
       write(stdout,'(a)') '|               Starting conduction properties &
            &calculation                     |'
       write(stdout,'(a)') '+--------------------------------------------------&
            &----------------------------+'
    end if

    ! lr408: Total number of eigenstates to use in spectra calculations
    num_tot_eigen = val_rep%n_occ(1) + cond_rep%n_occ(1)

    ! lr408: Check this isn't greater than total number of states in
    ! lr408: valence NGWF basis
    if (num_tot_eigen > val_ngwf_basis%num) then
       num_tot_eigen = val_ngwf_basis%num
    end if

    ! lr408: Calculate properties for different basis sets/Hamiltonians

    ! lr408: Store original values so they can be reset after temporary overrides
    tmp_homo_plot = pub_homo_plot
    tmp_lumo_plot = pub_lumo_plot
    tmp_homo_dens_plot = pub_homo_dens_plot
    tmp_lumo_dens_plot = pub_lumo_dens_plot

    if (.not. cond_plot_vc_orbitals) then
       pub_homo_plot = -1
       pub_lumo_plot = -1
       pub_homo_dens_plot = -1
       pub_lumo_dens_plot = -1
    end if

    call properties_calculate(val_denskern, val_ham, &
         val_rep, val_ngwf_basis, proj_basis, hub_proj_basis, hub, elements, &
         localpseudo_fine, core_density_fine, lhxc_fine, .false., 'valence', &
         num_tot_eigen)

    ! lr408: Total number of eigenstates to use in spectra calculations
    num_tot_eigen = val_rep%n_occ(1) + cond_rep%n_occ(1)

    ! lr408: Check this isn't greater than total number of states in
    ! lr408: conduction NGWF basis
    if (num_tot_eigen > cond_ngwf_basis%num) then
       num_tot_eigen = cond_ngwf_basis%num
    end if

    pub_homo_plot = -1
    pub_lumo_plot = -1
    pub_homo_dens_plot = -1
    pub_lumo_dens_plot = -1

    ! Unprojected Cond Properties calculation
    n_occ_cond = cond_rep%n_occ
    cond_rep%n_occ = val_rep%n_occ
    call properties_calculate(cond_denskern, cond_ham, &
         cond_rep, cond_ngwf_basis, proj_basis, hub_proj_basis, hub, &
         elements, localpseudo_fine, core_density_fine, lhxc_fine, &
         .false., 'cond')
    cond_rep%n_occ = n_occ_cond

    if (cond_plot_vc_orbitals) then
       pub_lumo_plot = tmp_lumo_plot
       pub_lumo_dens_plot = tmp_lumo_dens_plot
    end if

    ! Projected Cond Properties calculation
    call properties_calculate(cond_denskern, cond_ham, &
         cond_rep, cond_ngwf_basis, proj_basis, hub_proj_basis, hub, &
         elements, localpseudo_fine, core_density_fine, lhxc_fine, &
         .false., 'proj')

    num_tot_eigen = val_rep%n_occ(1) + cond_rep%n_occ(1)

    if (cond_plot_joint_orbitals) then
       pub_homo_plot = tmp_homo_plot
       pub_lumo_plot = tmp_lumo_plot
    else
       pub_homo_plot = -1
       pub_lumo_plot = -1
    end if
    ! lr408: Plotting orbital densities not yet implemented for joint basis
    !pub_homo_dens_plot = tmp_homo_dens_plot
    pub_homo_dens_plot = -1
    pub_lumo_dens_plot = -1

    ! ndmh: Allocate SPAM3 structures for overlap matrix, kinetic energy,
    ! ndmh: density kernel, inverse overlap
    call ngwf_rep_create(joint_rep,'j',elements)

    allocate(joint_denskern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('conduction_properties','joint_denskern',ierr)
    do is=1,pub_cell%num_spins
       joint_denskern(is)%structure = 'K'//joint_rep%postfix
       call sparse_create(joint_denskern(is))
    end do
    pub_kernel_workspace_allocated = .false.

    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a)') ' done'

    ! ndmh: Allocate joint rep NGWFs
    allocate(joint_rep%ngwfs_on_grid(joint_ngwf_basis%n_ppds*pub_cell%n_pts), &
         stat=ierr)
    call utils_alloc_check('conduction_properties', &
        'joint_rep%ngwfs_on_grid',ierr)
    call comms_barrier

    ! lr408: Initialise conduction NGWFs to fireballs or read from file
    if (pub_on_root) write(stdout,'(a)',advance='no') 'Joint Valence + &
         &Conduction NGWF initialisation ...'

    ! ndmh: create a merged set of NGWFs containing the valence and conduction
    ! ndmh: NGWFs previously calculated
    call ngwfs_merge_sets(joint_rep%ngwfs_on_grid,joint_ngwf_basis, &
         val_rep%ngwfs_on_grid,val_ngwf_basis, &
         cond_rep%ngwfs_on_grid, cond_ngwf_basis)
    call comms_barrier
    if (pub_on_root) write(stdout,'(a/)') '... done'
    call services_flush

    ! ndmh: copy in occupation numbers from the valence rep
    joint_rep%n_occ(:) = val_rep%n_occ(:)

    ! ndmh: allocate matrices in the joint hamiltonian
    call ngwf_ham_create(joint_ham,joint_rep)

    ! ndmh: copy in dijhat from valence hamiltonian
    if (pub_aug) then
       do is=1,pub_cell%num_spins
          call sparse_copy(joint_ham%dijhat(is),val_ham%dijhat(is))
       end do
    end if

    ! ndmh: calculate density-kernel-independent parts of the joint hamiltonian
    call hamiltonian_dens_indep_matrices(joint_rep, joint_ngwf_basis, &
         proj_basis, hub_proj_basis, hub, val_rep,val_ngwf_basis)

    ! ndmh: calculate the rest of the joint hamiltonian
    call hamiltonian_dens_dep_nonsc(joint_ham, joint_rep, joint_ngwf_basis, &
         lhxc_fine,hub,val_rep,val_ham,val_denskern)

    ! ndmh: perform a properties calculation using the joint hamiltonian
    call properties_calculate(val_denskern, joint_ham, &
         joint_rep, joint_ngwf_basis, proj_basis, hub_proj_basis, hub, &
         elements, localpseudo_fine, core_density_fine, lhxc_fine, &
         .false., 'joint', num_tot_eigen)

    ! ndmh: memory deallocation for joint rep/ham
    call ngwf_ham_destroy(joint_ham)
    deallocate(joint_rep%ngwfs_on_grid,stat=ierr)
    call utils_dealloc_check('conduction_properties', &
        'joint_rep%ngwfs_on_grid',ierr)
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(joint_denskern(is))
    end do
    deallocate(joint_denskern,stat=ierr)
    call utils_dealloc_check('conduction_properties','joint_denskern',ierr)
    call ngwf_rep_destroy(joint_rep)

    if (pub_on_root) then
       write(stdout,'(a)') '+--------------------------------------------------&
            &----------------------------+'
       write(stdout,'(a)') '|               Conduction properties calculation &
            &complete                     |'
       write(stdout,'(a)') '+--------------------------------------------------&
            &----------------------------+'
    end if

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &conduction_properties'
#endif


    return

  end subroutine conduction_properties_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module conduction_properties
