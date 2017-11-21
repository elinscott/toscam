! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                              Forces Module                                  !
!=============================================================================!
!                                                                             !
! This module calculates the total forces on each ion, after a converged      !
! total energy calculation.                                                   !
!                                                                             !
!-----------------------------------------------------------------------------!
! Module created by Nicholas Hine, 28th June 2010, based largely on existing  !
! routines from energy_and_force_mod, written by Chris-Kriton Skylaris, Arash !
! Mostofi, Peter Haynes and Nicholas Hine over the period 2000-2010.          !
!-----------------------------------------------------------------------------!

module forces

  implicit none

  private

  public :: forces_calculate
  public :: forces_apply_constraints

contains

  subroutine forces_calculate(total_forces,denskern,ham,lhxc_fine, &
       rep,ngwf_basis,projector_basis,hub_proj_basis,hub,localpseudo_fine, &
       core_density_fine,elements,ngwf_nonsc_forces)

    !================================================================!
    ! This subroutine calculates the ionic forces                    !
    !----------------------------------------------------------------!
    ! Written by Nicholas Hine, June 2010, based on the previous     !
    ! routine internal_forces in energy_and_force_calculate, written !
    ! by Arash A. Mostofi, July 2004.                                !
    !================================================================!

    use augmentation, only: augmentation_density_on_grid, &
         aug_nl_calculate_forces, augmentation_density_forces
    use cell_grid, only: pub_fine_grid
    use comms, only: comms_barrier, pub_on_root
    use constants, only: DP, paw_en_size, stdout, VERBOSE
    use classical_pot, only: classical_elements, classical_pot_ii_energy, &
         classical_pot_ii_forces
    use cutoff_coulomb, only: cutoff_coulomb_ii_forces,  &
         cutoff_coulomb_locps_forces, cutoff_coulomb_nlcc_forces
    use density, only: density_on_grid
    use ewald, only: ewald_calculate_forces
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix, hamiltonian_lhxc_calculate
    use hubbard_build, only: HUBBARD_MODEL, hubbard_calculate_forces
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP, NGWF_HAM, ngwf_ham_create, &
         ngwf_ham_destroy
    use paw, only: paw_nlcc_calculate_forces, paw_tcore_hartree_calc_forces
    use pseudopotentials, only: pseudo_local_calculate_forces, &
         pseudo_nlcc_calculate_forces, pseudo_nl_calculate_forces
    use rundat, only: pub_any_nl_proj, pub_coulomb_cutoff, pub_dispersion, &
         pub_hubbard, pub_nlcc, pub_paw, pub_aug, pub_output_detail, &
         pub_write_forces, pub_usehfx, pub_ii_energy_direct, print_qc, &
         pub_nonsc_forces, pub_aug_den_dim, pub_nhat_in_xc
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use vdwcorrection, only: vdwcorrection_calculate_forces

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: total_forces(3,pub_cell%nat)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: projector_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(inout) :: ham
    type(HUBBARD_MODEL), intent(in) :: hub
    real(kind=DP), intent(in) :: localpseudo_fine(:,:,:)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
    type(ELEMENT), intent(in)   :: elements(pub_cell%nat)
    ! ars: ngwf_nonsc_forces
    real(kind=DP), intent(in)  :: ngwf_nonsc_forces(1:3,pub_cell%nat)


    ! Local Variables
    integer :: ierr,atom
    real(kind=DP) :: average_force(1:3)
    real(kind=DP), allocatable, dimension(:,:) :: locps_forces
    real(kind=DP), allocatable, dimension(:,:) :: nlps_forces
    real(kind=DP), allocatable, dimension(:,:) :: nhat_forces
    real(kind=DP), allocatable, dimension(:,:) :: ewald_forces
    real(kind=DP), allocatable, dimension(:,:) :: vdw_forces
    real(kind=DP), allocatable, dimension(:,:) :: nlcc_forces
    real(kind=DP), allocatable, dimension(:,:) :: hub_forces
    real(kind=DP), allocatable, dimension(:,:) :: classical_forces
    ! ndmh: density on fine grid
    real(kind=DP), dimension(:,:,:,:), allocatable :: density_fine
    real(kind=DP), dimension(:,:,:,:,:), allocatable :: nhat_den_grad
    real(kind=DP) :: lhxc_energy   ! dummy for lhxc_calculate

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering forces_calculate'
#endif

    call timer_clock('forces_calculate',1)

    ! Allocate arrays for each force component
    allocate(ewald_forces(3,pub_cell%nat),stat=ierr)
    call utils_alloc_check('forces_calculate','ewald_forces',ierr)
    allocate(locps_forces(3,pub_cell%nat),stat=ierr)
    call utils_alloc_check('forces_calculate','locps_forces',ierr)
    if (pub_any_nl_proj.or.pub_paw) then
       allocate(nlps_forces(3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','nlps_forces',ierr)
    end if
    if (pub_aug) then
       allocate(nhat_forces(3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','nhat_forces',ierr)
    end if
    if (pub_dispersion/=0) then
       allocate(vdw_forces(3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','vdw_forces',ierr)
    end if
    if (pub_nlcc) then
       allocate(nlcc_forces(3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','nlcc_forces',ierr)
    end if
    if (pub_hubbard) then
       allocate(hub_forces(3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','hub_forces',ierr)
    end if


    ! kaw: Allocate array for classical forces
    if (pub_cell%nat_classical > 0) then
       allocate(classical_forces(3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('forces_calculate','classical_forces',ierr)
       classical_forces = 0.0_DP
    end if



    if (pub_coulomb_cutoff) then
       ! Calculate forces due to direct or Ewald contribution
         call cutoff_coulomb_ii_forces(elements,ewald_forces)

         ! kaw: Calculate classical forces
         if (pub_cell%nat_classical > 0) then
            call classical_pot_ii_forces(classical_forces,elements)
         end if

    else
       ! Calculate forces due to Ewald contribution
       call ewald_calculate_forces(elements,ewald_forces,.false.)
       ! Clean up Ewald module allocatables (use locps_forces as temp dummy)
       call ewald_calculate_forces(elements,locps_forces,.true.)
    end if

    ! aam: Allocate charge density slabs
    allocate(density_fine(pub_fine_grid%ld1,pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('forces_calculate','density_fine',ierr)
    if (pub_aug) then
       allocate(nhat_den_grad(pub_fine_grid%ld1,pub_fine_grid%ld2, &
            pub_fine_grid%max_slabs12,pub_cell%num_spins,0:pub_aug_den_dim), &
            stat=ierr)
       call utils_alloc_check('forces_calculate','nhat_den_grad',ierr)
    end if

    ! aam: Scale density kernel to the correct number of electrons
    if (pub_cell%num_spins == 1) call sparse_scale(denskern(1),2.0_DP)

    ! aam: Calculate data-parallelised charge density
    call density_on_grid(density_fine,pub_fine_grid, &
         denskern,rep%ngwf_overlap,rep%ngwfs_on_grid,ngwf_basis)

    ! ndmh: Calculate compensation density
    if (pub_aug) then
       nhat_den_grad = 0.0_DP
       call augmentation_density_on_grid(nhat_den_grad,pub_fine_grid,denskern, &
            rep%sp_overlap)
       density_fine = density_fine + nhat_den_grad(:,:,:,:,0)
    end if

    if (pub_coulomb_cutoff) then

       ! Calculate forces due to local part of ionic pseudopotential
       call cutoff_coulomb_locps_forces(density_fine,elements,locps_forces)

       ! Calculate forces due to the NLCC core charge
       if (pub_nlcc) then
          if (pub_aug) then
             call cutoff_coulomb_nlcc_forces(density_fine, &
                  core_density_fine,elements,nlcc_forces,nhat_den_grad)
          else
             call cutoff_coulomb_nlcc_forces(density_fine, &
                  core_density_fine,elements,nlcc_forces)
          end if
       end if

    else

       ! Calculate forces due to local part of ionic pseudopotential
       if (.not.pub_paw) then
          call pseudo_local_calculate_forces(density_fine,pub_fine_grid, &
               elements,locps_forces)
       else
          call paw_tcore_hartree_calc_forces(density_fine,pub_fine_grid, &
               elements,locps_forces)
       end if

       if (pub_nlcc) then
          ! Calculate forces due to the NLCC core charge
          if (.not.pub_paw) then
             call pseudo_nlcc_calculate_forces(density_fine, &
                  core_density_fine,pub_fine_grid,elements,nlcc_forces)
          else
             call paw_nlcc_calculate_forces(density_fine, &
                  core_density_fine,nhat_den_grad,pub_fine_grid,elements, &
                  nlcc_forces)
          end if
       end if

    end if

    ! ndmh: Calculate nonlocal pseudopotential forces if necessary
    if ((.not.pub_aug).and.pub_any_nl_proj) then

       call pseudo_nl_calculate_forces(nlps_forces, &
            rep%sp_overlap,rep%ngwfs_on_grid,ngwf_basis, &
            projector_basis,denskern)

    ! ndmh: Calculate PAW force terms
    else if (pub_aug) then

       ! ndmh: recalculate LHXC and ham%dijhat
       call hamiltonian_lhxc_calculate(lhxc_fine,lhxc_energy,ham%dijhat, &
            pub_fine_grid,localpseudo_fine,core_density_fine, &
            rep%ngwfs_on_grid,ngwf_basis,denskern,rep%ngwf_overlap, &
            rep%sp_overlap,add_xc_pot=pub_nhat_in_xc)

       ! Calculate PAW nonlocal force terms
       call aug_nl_calculate_forces(nlps_forces, &
            rep%ngwfs_on_grid,ngwf_basis,projector_basis, &
            rep%sp_overlap,rep%inv_overlap,denskern,ham%ham,ham%dijhat)

       ! Calculate PAW compensation density forces
       call augmentation_density_forces(nhat_forces,denskern,rep%sp_overlap, &
            lhxc_fine,pub_fine_grid)

    end if

    ! qoh: Calculate vdw_forces if necessary
    if (pub_dispersion /= 0) then
       call vdwcorrection_calculate_forces(vdw_forces,elements)
    end if

    ! ddor: Calculate DFT+U forces if necessary
    if (pub_hubbard) then
       ! Remove spin-degeneracy in non spin-polarised denskern for DFT+U
       if (pub_cell%num_spins == 1) call sparse_scale(denskern(1), 0.5_DP)
       call hubbard_calculate_forces(hub_forces, &
            rep%ngwfs_on_grid,ngwf_basis,hub_proj_basis,hub, &
            denskern,rep%hub_overlap,rep%hub_overlap_t)
       if (pub_cell%num_spins == 1) call sparse_scale(denskern(1), 2.0_DP)
    endif


    !==============================================================!
    !*** Uncomment in order to convert forces from Eh/a to eV/A ***!
    !     ewald_forces =      ewald_forces*HARTREE_IN_EVS*ANGSTROM !
    !     locps_forces =      locps_forces*HARTREE_IN_EVS*ANGSTROM !
    !      nlps_forces =       nlps_forces*HARTREE_IN_EVS*ANGSTROM !
    !      nhat_forces =       nhat_forces*HARTREE_IN_EVS*ANGSTROM !
    !       vdw_forces =        vdw_forces*HARTREE_IN_EVS*ANGSTROM !
    !      nlcc_forces =       nlcc_forces*HARTREE_IN_EVS*ANGSTROM !
    !      nhat_forces =       nhat_forces*HARTREE_IN_EVS*ANGSTROM !
    !ngwf_nonsc_forces = ngwf_nonsc_forces*HARTREE_IN_EVS*ANGSTROM !
    !==============================================================!

    ! Sum all of the force contributions
    total_forces = 0.0_DP

    !kaw: Add classical forces when relevant
    if ((pub_cell%nat_classical > 0).and.(pub_coulomb_cutoff)) then
       ewald_forces = ewald_forces + classical_forces
    end if
    total_forces = total_forces + ewald_forces + locps_forces

    if (pub_any_nl_proj.or.pub_paw) then
       total_forces = total_forces + nlps_forces
    end if
    if (pub_aug) then
       total_forces = total_forces + nhat_forces
    end if
    if (pub_dispersion/=0) then
       total_forces = total_forces + vdw_forces
    end if
    if (pub_nlcc) then
       total_forces = total_forces + nlcc_forces
    end if
    if (pub_nonsc_forces) then
       if (all(ngwf_nonsc_forces==-999_DP)) then
          if (pub_on_root) write(stdout,'(a)') &
               'WARNING in forces_calculate: NGWF non-self-consistent &
               &forces not present. Omitting from total.'
       else
          total_forces = total_forces + ngwf_nonsc_forces
       end if
    end if
    if (pub_hubbard) then
       total_forces = total_forces + hub_forces
    end if

    ! Calculate the average force
    ! kaw: Account for classical atoms
    average_force = 0.0_dp
    do atom=1,pub_cell%nat - pub_cell%nat_classical
       average_force(:) = average_force(:) + total_forces(:,atom)
    end do
    average_force(:) = average_force(:) / (pub_cell%nat - pub_cell%nat_classical)

    ! Remove the average force
    do atom=1,pub_cell%nat - pub_cell%nat_classical
       total_forces(:,atom) = total_forces(:,atom) - average_force(:)
    end do

    ! Write forces (root node only)
    if (pub_on_root) then
       if (pub_write_forces) then
          if (pub_output_detail == VERBOSE) then
             write(stdout,'(/a)') '********************* Unconstrained &
                  &**********************'
             if (pub_ii_energy_direct) then
                write(stdout,'(a)') '******************** Ion-Ion forces  &
                     &*********************'
             else
                write(stdout,'(a)') '********************* Ewald forces  &
                     &**********************'
             end if
             write(stdout,'(a)') '*                                   &
                  &                     *'
             write(stdout,'(a)') '* Element  Atom         &
                  &Cartesian components (Eh/a)      *'
             write(stdout,'(a)') '* ----------------------------------&
                  &-------------------- *'
             write(stdout,'(a)') '*                       x            &
                  &y            z      *'
             write(stdout,'(a)') '*                                   &
                  &                     *'
             do atom=1,pub_cell%nat
                write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                     elements(atom)%symbol,' ',atom,'   ', &
                     ewald_forces(:,atom),' *'
             end do
             write(stdout,'(a)') '*                                   &
                  &                     *'
             write(stdout,'(a)') '************************************&
                  &**********************'

             write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
                  sum(ewald_forces(1,1:pub_cell%nat)), &
                  sum(ewald_forces(2,1:pub_cell%nat)), &
                  sum(ewald_forces(3,1:pub_cell%nat))

             write(stdout,'(a)') ' '
             write(stdout,'(a)') '********************* Unconstrained &
                  &**********************'
             write(stdout,'(a)') '***************** Local potential forces &
                  &*****************'
             write(stdout,'(a)') '*                                  &
                  &                      *'
             write(stdout,'(a)') '* Element  Atom         &
                  &Cartesian components (Eh/a)      *'
             write(stdout,'(a)') '* ---------------------------------&
                  &--------------------- *'
             write(stdout,'(a)') '*                       x            &
                  &y            z      *'
             write(stdout,'(a)') '*                                  &
                  &                      *'
             do atom=1,pub_cell%nat
                write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                     elements(atom)%symbol,' ',atom,'   ', &
                     locps_forces(:,atom),' *'
             end do
             write(stdout,'(a)') '*                                  &
                  &                      *'
             write(stdout,'(a)') '***********************************&
                  &***********************'

             write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
                  sum(locps_forces(1,1:pub_cell%nat)), &
                  sum(locps_forces(2,1:pub_cell%nat)), &
                  sum(locps_forces(3,1:pub_cell%nat))

             if (pub_aug) then
                write(stdout,'(a)') ' '
                write(stdout,'(a)') '********************* Unconstrained &
                     &**********************'
                write(stdout,'(a)') '************** Compensation density &
                     &forces ***************'
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '* Element  Atom         &
                     &Cartesian components (Eh/a)      *'
                write(stdout,'(a)') '* ----------------------------------&
                     &-------------------- *'
                write(stdout,'(a)') '*                       x            &
                     &y            z      *'
                write(stdout,'(a)') '*                                 &
                     &                       *'
                do atom=1,pub_cell%nat
                   write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                        elements(atom)%symbol,' ',atom,'   ', &
                        nhat_forces(:,atom),' *'
                end do
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '***********************************&
                     &***********************'
                write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
                     sum(nhat_forces(1,1:pub_cell%nat)),  &
                     sum(nhat_forces(2,1:pub_cell%nat)),  &
                     sum(nhat_forces(3,1:pub_cell%nat))

                write(stdout,'(a)') ' '
             end if

             if (pub_paw.or.pub_any_nl_proj) then
                write(stdout,'(a)') '********************* Unconstrained &
                     &**********************'
                write(stdout,'(a)') '*************** Non-local potential &
                     &forces ***************'
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '* Element  Atom         &
                     &Cartesian components (Eh/a)      *'
                write(stdout,'(a)') '* ----------------------------------&
                     &-------------------- *'
                write(stdout,'(a)') '*                       x            &
                     &y            z      *'
                write(stdout,'(a)') '*                                 &
                     &                       *'
                do atom=1,pub_cell%nat
                   write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                        elements(atom)%symbol,' ',atom,'   ', &
                        nlps_forces(:,atom),' *'
                end do
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '***********************************&
                     &***********************'
                write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
                     sum(nlps_forces(1,1:pub_cell%nat)),  &
                     sum(nlps_forces(2,1:pub_cell%nat)),  &
                     sum(nlps_forces(3,1:pub_cell%nat))
             end if

             ! ars: NGWF non self consistent forces
             if(pub_nonsc_forces) then
                write(stdout,'(a)') ' '
                write(stdout,'(a)') '********************* Unconstrained &
                     &**********************'
                write(stdout,'(a)') '************ NGWF non &
                     &self-consistent forces *************'
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '* Element  Atom         &
                     &Cartesian components (Eh/a)      *'
                write(stdout,'(a)') '* ----------------------------------&
                     &-------------------- *'
                write(stdout,'(a)') '*                       x            &
                     &y            z      *'
                write(stdout,'(a)') '*                                 &
                     &                       *'
                do atom=1,pub_cell%nat
                   write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                     elements(atom)%symbol,' ',atom,'   ', &
                         ngwf_nonsc_forces(:,atom),' *'
                end do
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '***********************************&
                     &***********************'
                write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
                     sum(ngwf_nonsc_forces(1,1:pub_cell%nat)),  &
                     sum(ngwf_nonsc_forces(2,1:pub_cell%nat)),  &
                     sum(ngwf_nonsc_forces(3,1:pub_cell%nat))
             end if

             if (pub_dispersion /=0) then
                write(stdout,'(a)') ' '
                write(stdout,'(a)') '********************* Unconstrained &
                     &**********************'
                write(stdout,'(a)') '************************ Dispersion &
                     &forces ***************'
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '* Element  Atom         &
                     &Cartesian components (Eh/a)      *'
                write(stdout,'(a)') '* ----------------------------------&
                     &-------------------- *'
                write(stdout,'(a)') '*                       x            &
                     &y            z      *'
                write(stdout,'(a)') '*                                 &
                     &                       *'
                do atom=1,pub_cell%nat
                   write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                        elements(atom)%symbol,' ',atom,'   ', &
                        vdw_forces(:,atom),' *'
                end do
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '***********************************&
                     &***********************'
                write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
                     sum(vdw_forces(1,1:pub_cell%nat)),  &
                     sum(vdw_forces(2,1:pub_cell%nat)),  &
                     sum(vdw_forces(3,1:pub_cell%nat))
             end if

             if (pub_nlcc) then
                write(stdout,'(a)') ' '
                write(stdout,'(a)') '********************* Unconstrained &
                     &**********************'
                write(stdout,'(a)') '********************** NLCC forces &
                     &**********************'
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '* Element  Atom         &
                     &Cartesian components (Eh/a)      *'
                write(stdout,'(a)') '* ----------------------------------&
                     &-------------------- *'
                write(stdout,'(a)') '*                       x            &
                     &y            z      *'
                write(stdout,'(a)') '*                                 &
                     &                       *'
                do atom=1,pub_cell%nat
                   write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                        elements(atom)%symbol,' ',atom,'   ', &
                        nlcc_forces(:,atom),' *'
                end do
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '***********************************&
                     &***********************'
                write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
                     sum(nlcc_forces(1,1:pub_cell%nat)),  &
                     sum(nlcc_forces(2,1:pub_cell%nat)),  &
                     sum(nlcc_forces(3,1:pub_cell%nat))
             end if

             if (pub_hubbard) then
                write(stdout,'(a)') ' '
                write(stdout,'(a)') '********************* Unconstrained &
                     &**********************'
                write(stdout,'(a)') '********************** DFT+U forces &
                     &**********************'
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '* Element  Atom         &
                     &Cartesian components (Eh/a)      *'
                write(stdout,'(a)') '* ----------------------------------&
                     &-------------------- *'
                write(stdout,'(a)') '*                       x            &
                     &y            z      *'
                write(stdout,'(a)') '*                                 &
                     &                       *'
                do atom=1,pub_cell%nat
                   write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                        elements(atom)%symbol,' ',atom,'   ', &
                        hub_forces(:,atom),' *'
                end do
                write(stdout,'(a)') '*                                  &
                     &                      *'
                write(stdout,'(a)') '***********************************&
                     &***********************'
                write(stdout,'(a,3f13.8)') 'TOTAL:           ',&
                     sum(hub_forces(1,1:pub_cell%nat)),  &
                     sum(hub_forces(2,1:pub_cell%nat)),  &
                     sum(hub_forces(3,1:pub_cell%nat))
             end if
          end if

          write(stdout,'(/a)') '********************* Unconstrained &
               &**********************'
          write(stdout,'(a)') '************************* Forces &
               &*************************'
          write(stdout,'(a)') '*                                &
               &                        *'
          write(stdout,'(a)') '* Element  Atom         Cartesian &
               &components (Eh/a)      *'
          write(stdout,'(a)') '* -------------------------------&
               &----------------------- *'
          write(stdout,'(a)') '*                       x            &
               &y            z      *'
          write(stdout,'(a)') '*                                 &
               &                       *'
          do atom=1,pub_cell%nat
             write(stdout,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
                  elements(atom)%symbol,' ',atom,'   ', &
                  total_forces(:,atom),' *'
          end do
          write(stdout,'(a)') '*                                 &
               &                       *'
          write(stdout,'(a)') '**********************************&
               &************************'
          write(stdout,'(a,3f13.8,a)') '* TOTAL:         ',&
               sum(total_forces(1,1:pub_cell%nat)), &
               sum(total_forces(2,1:pub_cell%nat)), &
               sum(total_forces(3,1:pub_cell%nat)),' *'
          write(stdout,'(a)') '**********************************&
               &************************'
       end if
    end if
    call comms_barrier

    ! aam: Scale density kernel back to half Ne
    ! ddor: This is already done in the case of a DFT+U calculation
    if ((pub_cell%num_spins == 1) .and. (.not. pub_hubbard)) &
         & call sparse_scale(denskern(1), 0.5_DP)

    ! Deallocate
    if (pub_aug) then
       deallocate(nhat_den_grad,stat=ierr)
       call utils_dealloc_check('forces_calculate','nhat_den_grad',ierr)
    end if
    deallocate(density_fine,stat=ierr)
    call utils_dealloc_check('forces_calculate','density_fine',ierr)
    if (pub_hubbard) then
       deallocate(hub_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','hub_forces',ierr)
    end if
    if (pub_nlcc) then
       deallocate(nlcc_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','nlcc_forces',ierr)
    end if
    if (pub_dispersion/=0) then
       deallocate(vdw_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','vdw_forces',ierr)
    end if
    if (pub_aug) then
       deallocate(nhat_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','nhat_forces',ierr)
    end if
    if (pub_paw.or.pub_any_nl_proj) then
       deallocate(nlps_forces,stat=ierr)
       call utils_dealloc_check('forces_calculate','nlps_forces',ierr)
    end if
    deallocate(locps_forces,stat=ierr)
    call utils_dealloc_check('forces_calculate','locps_forces',ierr)
    deallocate(ewald_forces,stat=ierr)
    call utils_dealloc_check('forces_calculate','ewald_forces',ierr)

    if (print_qc) call forces_qc_output(total_forces)

    call timer_clock('forces_calculate',2)

#ifdef DEBUG
  if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving forces_calculate'
#endif

    return
  end subroutine forces_calculate

  subroutine forces_apply_constraints(forces,elements)
    !=========================================================================!
    ! Apply ionic constraints to the ionic forces                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   forces, intent=inout, the ionic forces                                !
    !   elements, intent=in, the element data array                           !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 02/08/2005                                  !
    !=========================================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, stdout
    use ion, only: ELEMENT
    use simulation_cell, only: pub_cell

    implicit none

    real(kind=dp), intent(inout) :: forces(3,1,pub_cell%nat)
    type(element), intent(in)    :: elements(pub_cell%nat)

    ! Local Variables
    integer :: iat
    real(kind=dp) :: proj,norm
    real(kind=dp) :: f(1:3),v(1:3)

    do iat=1,pub_cell%nat

       f(1:3) = forces(1:3,1,iat)
       v(1:3) = elements(iat)%ion_constraint(1:3)

       select case (elements(iat)%ion_constraint_type)
       case ('NONE') ; continue                    ! aam: no constraint
       case ('PLANE')                              ! aam: constrained perp. to v
          proj = f(1)*v(1) + f(2)*v(2) + f(3)*v(3)           ! aam: f.v
          norm = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)           ! aam: v.v
          forces(1:3,1,iat) = f(1:3) - (proj/norm)*v(1:3)    ! aam: f-(f.v/v.v)v
       case ('LINE')                               ! aam: constrained para. to v
          proj = f(1)*v(1) + f(2)*v(2) + f(3)*v(3)           ! aam: f.v
          norm = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)           ! aam: v.v
          forces(1:3,1,iat) = (proj/norm)*v(1:3)             ! aam: (f.v/v.v)v
       case ('FIXED') ; forces(1:3,1,iat) = 0.0_dp ! aam: fixed
       case default
          if (pub_on_root) then
             write(stdout,'(a,i6,a)') 'Error in forces_apply_constraints: &
                  &illegal value for elements(',iat,')%ion_constraint_type'
          endif
       end select

    enddo

  end subroutine forces_apply_constraints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine forces_qc_output(total_forces)

    !=========================================================================!
    ! Prints out quality-control lines <QC> for qc-testing the forces.        !
    !-------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, 19/11/2010.                             !
    !=========================================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use simulation_cell, only: pub_cell


    implicit none

    real(kind=DP), intent(in   ) ::  total_forces(1:3, 1:pub_cell%nat)

    integer :: icomp, iatom

    if(pub_on_root) then

       do icomp = 1, 3
          do iatom = 1, pub_cell%nat
             write(stdout,'(a22, i0, a1, i0, a3,f22.12)') &
                "<QC>    [total_force(", icomp, ",", iatom, ")]:", &
                  total_forces(icomp, iatom)
          end do
       end do

    end if

    call comms_barrier

  end subroutine forces_qc_output

end module forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
