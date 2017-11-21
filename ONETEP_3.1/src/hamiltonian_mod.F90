! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Arash A. Mostofi and Nicholas D.M. Hine
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module hamiltonian

  use constants, only: DP

  private

  public :: hamiltonian_dens_indep_matrices
  public :: hamiltonian_dens_dep_matrices
  public :: hamiltonian_dens_dep_nonsc
  public :: hamiltonian_lhxc_calculate
  public :: hamiltonian_build_matrix
  public :: hamiltonian_proj_cond_matrix
  public :: hamiltonian_energy_components

contains

  subroutine hamiltonian_dens_indep_matrices(rep, &       ! out
       ngwf_basis, projector_basis, hub_proj_basis,hub, & ! in
       val_rep, val_ngwf_basis)                           ! in, optional

    !=============================================================!
    ! This subroutine initialises the matrices which make up the  !
    ! Hamiltonian and are independent of the density kernel.      !
    !-------------------------------------------------------------!
    ! Arguments:                                                  !
    ! rep             (inout) : NGWF Representation (functions    !
    !                           and matrices).                    !
    ! ngwf_basis      (input) : Function basis type for the NGWFs !
    ! projector_basis (input) : Function basis type for nonlocal  !
    !                           pseudopotential projectors        !
    ! hub_proj_basis  (input) : Function basis type for Hubbard   !
    !                           projectors                        !
    ! val_rep         (input) : Valence NGWF Representation       !
    !                           !optional)                        !
    ! val_ngwf_basis  (input) : Function basis type for the       !
    !                           valence NGWFs (optional)          !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 3/7/2001.               !
    ! DFT+U added by David O'Regan, April 2008                    !
    ! Adapted for SPAM3, function basis and new pseudopotential   !
    ! routines by Nicholas Hine, May-July 2009.                   !
    ! Moved to hamiltonian_mod, April 2010.                       !
    ! Calculation of valence-conduction overlap matrix added,     !
    ! with necessary optional arguments by Laura Ratcliff,        !
    ! Oct 2010.                                                   !
    !=============================================================!

    use augmentation, only: augmentation_overlap
    use comms, only: comms_abort, pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use dense, only: DEM,dense_create,dense_destroy,dense_convert,dense_invert
    use function_basis, only: FUNC_BASIS
    use hubbard_build, only: HUBBARD_MODEL, hubbard_projection_mtx
    use integrals, only: integrals_brappd_ketppd, integrals_kinetic
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_projectors
    use projectors, only: projectors_func_ovlp_box
    use pseudopotentials, only: pseudopotentials_nonlocal_mat, nlps_projectors
    use rundat, only: pub_any_nl_proj, pub_hubbard, &
         task, pub_hubbard_restart, pub_hubbard_atomsolve, &
         pub_paw, maxit_hotelling, &
         max_resid_hotelling, pub_output_detail, pub_usp, pub_aug, &
         pub_realspace_projectors
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_scale, sparse_hotelling_init, sparse_hotelling_invert
    use wrappers, only: wrappers_invert_sym_matrix
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: projector_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP), optional, intent(in) :: val_rep
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis

    ! Local Variables
    type(DEM) :: inv_overlap_dens
    real(kind=DP) :: final_max_resid

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering hamiltonian_dens_indep_matrices'
#endif

    ! %%%%%%%%%%%%%%%%%%%%%%%%% OVERLAP MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Starting overlap matrix'
#endif

    call integrals_brappd_ketppd(rep%overlap, &
         rep%ngwfs_on_grid, ngwf_basis, rep%ngwfs_on_grid, ngwf_basis)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Finished overlap matrix'
#endif

    ! %%%%%%%%%%%%%%%%%%%%%%% END OVERLAP MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ========================= KINETIC MATRIX ===============================
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Starting kinetic matrix'
#endif

    call integrals_kinetic(rep%kinet, rep%ngwfs_on_grid, ngwf_basis, &
         rep%ngwfs_on_grid, ngwf_basis)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Finished kinetic matrix'
#endif
    ! ======================= END KINETIC MATRIX =============================

    ! ================= NON-LOCAL PSEUDOPOTENTIAL MATRIX =====================
    ! calculate the ngwf-projector overlap matrix
    if (pub_any_nl_proj) then
#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: Starting nonlocal matrix'
#endif

       if (.not.pub_realspace_projectors) then
          call projectors_func_ovlp_box(rep%sp_overlap, &
               rep%ngwfs_on_grid,ngwf_basis,projector_basis,nlps_projectors)
       else
          call integrals_brappd_ketppd(rep%sp_overlap,rep%ngwfs_on_grid, &
               ngwf_basis,nlps_projectors%projs_on_grid,projector_basis)
       end if

       if (.not.pub_usp) then
          call pseudopotentials_nonlocal_mat(rep%nonlocpot(1),rep%sp_overlap)
       end if

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: Finished nonlocal matrix'
#endif
    else if (pub_paw) then
#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: Starting sp overlap matrix'
#endif
       ! ndmh: Calculate the overlap matrix of the NGWFs and PAW projectors
       call projectors_func_ovlp_box(rep%sp_overlap, &
            rep%ngwfs_on_grid,ngwf_basis,projector_basis,paw_projectors)

       ! ndmh: NB - the nonlocal matrix cannot be calculated at this stage in
       ! ndmh: PAW, since it depends on the density kernel!
#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: Finished sp overlap matrix'
#endif
    end if

    if (pub_aug) then
       ! ndmh: Calculate the augmentation of the overlap matrix due to the
       ! ndmh: augmentation region part of the overlap operator
       call augmentation_overlap(rep%overlap,rep%sp_overlap)
    end if
    ! =============== END NON-LOCAL PSEUDOPOTENTIAL MATRIX ===================


    ! %%%%%%%%%%%%%%%%%%%%%%%%% CROSS OVERLAP MATRIX %%%%%%%%%%%%%%%%%%%%%%%%
    !+ lr408: calculate cross overlap matrix if necessary
    if (present(val_rep)) then

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: Starting cross overlap &
            &matrix'
#endif

       call integrals_brappd_ketppd(rep%cross_overlap, &
            val_rep%ngwfs_on_grid, val_ngwf_basis, rep%ngwfs_on_grid, &
            ngwf_basis)

       if (pub_aug) then
          ! ndmh: Calculate the augmentation of the cross overlap matrix due to
          ! ndmh: the augmentation region part of the overlap operator
          call augmentation_overlap(rep%cross_overlap,val_rep%sp_overlap, &
               rep%sp_overlap)
       end if

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: Finished cross overlap &
            &matrix'
#endif

    end if

    ! %%%%%%%%%%%%%%%%%%%%%%% END CROSS OVERLAP MATRIX %%%%%%%%%%%%%%%%%%%%%%



    ! ====================== DFT+U CORRECTION MATRIX =========================

    ! ddor: DFT+U calculation
    if (pub_hubbard) then

       ! ddor: DFT+U calculation with self-consistency over Hubbard projectors
       consist: if ( ((task == 'HUBBARDSCF') .and.  &
            & (hub%consistency_iteration > 1)) .or. &
            & pub_hubbard_restart .or. pub_hubbard_atomsolve) then

#ifdef DEBUG
          if (pub_on_root) write(stdout,'(a,i4)') 'DEBUG: Hubbard correction &
               &matrix on HUBBARDSCF iteration', hub%consistency_iteration
#endif

          ! ddor: Using the hub_overlap from previous Hubbard
          ! iteration calculate the new hub_overlap from NGWFs.
          call integrals_brappd_ketppd(rep%hub_overlap, &    ! input-output
               rep%ngwfs_on_grid, ngwf_basis, hub%consistency_projs, & !input
               hub_proj_basis)                                     !input

          ! ddor: Augment the <NGWF|proj> matrix for PAW+U or USP calculations
          if (pub_aug) then

             ! ddor: Calculate the overlap matrix of the 
             ! ddor: Hubbard and PAW projectors
             call projectors_func_ovlp_box(rep%hub_proj_paw_overlap, &
                  hub%consistency_projs,hub_proj_basis,&
                  projector_basis,paw_projectors)

             ! ddor: Calculate the augmentation of the 
             ! ddor: overlap matrix due to the
             ! ddor: augmentation region part of the overlap operator
             call augmentation_overlap(rep%hub_overlap,rep%sp_overlap,&
                  rep%hub_proj_paw_overlap)

          end if

          ! calculate the Hubbard projection operator for each site,
          ! using hub_overlap and the tensorial correction matrix
          call hubbard_projection_mtx(.true.,.false., rep%hub_overlap, &
               rep%hub_overlap_t, hub%o_matrix)

       elseif ( ((task .ne. 'HUBBARDSCF') .and. &
            & (.not. pub_hubbard_restart) .and. &
            & (.not. pub_hubbard_atomsolve) ) .or. &
            & ( (task == 'HUBBARDSCF') .and. &
            & ((hub%consistency_iteration .eq. 1) .and. &
            & (.not. pub_hubbard_restart) .and. &
            & (.not. pub_hubbard_atomsolve) ) .and. &
            & hub%dftu_on_first_hubbardscf ) ) then
            ! ddor: Conventional DFT+U calculation

#ifdef DEBUG
          if (pub_on_root) write(stdout,'(a)') 'DEBUG: Starting &
               &hubbard_sp_ovlp_box'
#endif

          ! ddor: Calculate the Hubbard on-site <ngwf|atomic> overlap matrix
          call projectors_func_ovlp_box(rep%hub_overlap, &             ! out
               rep%ngwfs_on_grid,ngwf_basis,hub_proj_basis,hub%projectors) ! in

#ifdef DEBUG
          if (pub_on_root) write(stdout,'(a)') 'DEBUG: Finished &
               &hubbard_sp_ovlp_box'
#endif

          ! ddor: calculate the Hubbard projection operator for each site, with
          ! ddor: no tensorial correction (.false.) since we use atomic
          ! ddor: projectors
          call hubbard_projection_mtx(.false.,.false., rep%hub_overlap, &
               rep%hub_overlap_t, hub%o_matrix)

       endif consist

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: Finished Hubbard &
            &correction matrix'
#endif

    endif ! DFT+U

    ! ====================== DFT+U CORRECTION MATRIX =========================

    ! ======================== INVERSE OVERLAP MATRIX ========================

    ! if this is the first call to this routine...
    if (.not.rep%inv_overlap_init) then
       ! cks: initialise inverse overlap guess if Hotelling recursion
       ! cks: will be used to approximate it.
       ! ndmh: call sparse_mod version of hotelling's algorithm initialisation
       if (maxit_hotelling > 0) call sparse_hotelling_init(rep%inv_overlap, &
            rep%overlap)
       rep%inv_overlap_init = .true.
    end if

    ! cks : approximate inverse overlap by Hotelling's recursion
    if (maxit_hotelling > 0) then

       if (pub_output_detail>=VERBOSE .and. pub_on_root) then
          write(stdout,'(a)')'============ Calculation of NGWF S^-1 using &
               &Hotelling algorithm ================ '
          write(stdout,'(a)')'         Iteration  |  Frobenius norm and Maximum &
               &abs value    (of I-S*S_n^-1)   '
       end if

       ! ndmh: call sparse_mod version of hotelling's algorithm
       call sparse_hotelling_invert(rep%inv_overlap, & !inout
            rep%overlap,show_output=(pub_output_detail>=VERBOSE), &
            max_resid_converged=max_resid_hotelling, &
            num_iter=maxit_hotelling,final_max_resid=final_max_resid)

       ! ndmh: test if Hotelling failed to invert the matrix
       if (final_max_resid > 0.95_DP) then

          if (pub_output_detail>=VERBOSE.and.pub_on_root) then
             write(stdout,'(15x,a)')'Resetting Hotelling algorithm'
          end if

          ! ndmh: re-initialise Hotelling algorithm from scratch
          call sparse_hotelling_init(rep%inv_overlap,rep%overlap)

          ! ndmh: call sparse_mod version of hotelling's algorithm
          call sparse_hotelling_invert(rep%inv_overlap, & !inout
               rep%overlap,show_output=(pub_output_detail>=VERBOSE), &
               max_resid_converged=max_resid_hotelling, &
               num_iter=maxit_hotelling,final_max_resid=final_max_resid)

          ! ndmh: check if it is still broken
          if (final_max_resid > 0.95_DP) then
             call utils_abort('Error in hamiltonian_dens_indep_matrices: &
                  &Inversion of overlap matrix failed')
          end if

       end if

       if (pub_output_detail>=VERBOSE.and.pub_on_root) then
          write(stdout,'(a)')'===================================&
               &============================================='
       end if

    else if (maxit_hotelling == 0) then

       ! ndmh: allocate storage for dense inverse
       call dense_create(inv_overlap_dens,ngwf_basis%num,ngwf_basis%num)

       ! ndmh: copy current overlap matrix into dense inverse overlap
       call dense_convert(inv_overlap_dens, rep%overlap)

       ! ndmh: invert dense matrix
       call dense_invert(inv_overlap_dens)

       ! ndmh: copy dense inverse overlap back to sparse matrix
       call dense_convert(rep%inv_overlap,inv_overlap_dens)

       ! ndmh: deallocate storage for dense inverse
       call dense_destroy(inv_overlap_dens)

    else if (maxit_hotelling < 0) then
       if (pub_on_root) write(stdout,'(a,i6)') &
            'Error in hamiltonian_dens_indep_matrices: &
            &negative maxit_hotelling ',maxit_hotelling
       call comms_abort
    end if

    ! ==================== END INVERSE OVERLAP MATRIX ========================

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hamiltonian_dens_indep_matrices'
#endif

  end subroutine hamiltonian_dens_indep_matrices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CW
  subroutine hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
       lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
       ngwf_basis, hub_proj_basis, hub, denskern, ewald_energy, elements, &
       grid, localpseudo_fine, core_density_fine, ham_update, lhxc_fixed, &
       spoil_force, h_xc_write)
!END CW
    !===================================================================!
    ! This subroutine calculates density-dependent contributions to the !
    ! energy and components of the Hamiltonian matrix.                  !
    ! Updating of the Hamiltonian matrix components is optional, and is !
    ! controlled by the logical flag ham_update                         !
    !-------------------------------------------------------------------!
    ! Written by Nicholas Hine on 22/04/2010, assembled from bits of    !
    ! code formerly in electronic_mod.                                  !
    !===================================================================!
!CW
    use rundat, only: pub_dmft_spoil_kernel
    use sparse, only: sparse_axpy
!END CW
    use augmentation, only: aug_projector_denskern, aug_nonlocal_mat
    use cell_grid, only: GRID_INFO, pub_fine_grid ! ddor-28feb17
    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, paw_en_size, &
         paw_en_ehart, paw_en_exc, paw_en_exc_dc, paw_en_etxc, paw_en_etxc_dc, &
         paw_en_dijhat, paw_en_dijxc, paw_en_exc_core, stdout
    use function_basis, only: FUNC_BASIS
    use hf_exchange, only: hf_exchange_calculate
    use hubbard_build, only: HUBBARD_MODEL, hubbard_ham_matrix, &
         hubbard_energy_total, hubbard_projector_ham
    use integrals, only: integrals_locpot
    use ion, only: ELEMENT
    use is_smeared_ions, only: smeared_ion_E_self, smeared_ion_E_smeared
    use ngwf_representation, only: NGWF_HAM, NGWF_REP
    use rundat, only: pub_usehfx, pub_hubbard, pub_paw, pub_any_nl_proj, &
         pub_is_smeared_ion_rep, pub_is_include_cavitation, pub_aug, pub_usp
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale, &
         sparse_trace,sparse_product
    use utils, only: utils_alloc_check, utils_dealloc_check
    use vdwcorrection, only: pub_dispersion_energy
!CW
!    use IFPORT, only : system
    use sparse,only : sparse_copy
    use rundat , only : pub_dmft_points,pub_dmft_fully_sc,pub_dmft_fully_sc_h, &
        pub_aug_den_dim, pub_nlcc
    use comms, only: pub_my_node_id
    use eigenstates, only: eigenstates_calculate
    use hubbard_build, only: hubbard_dmft_interface
    use restart, only: restart_kernel_write 
    use xc_component, only: xc_energy_potential_comp, xc_init_comp
    use density, only: density_on_grid
    use augmentation, only: augmentation_density_on_grid
!END CW

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(inout) :: ham
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(SPAM3), intent(inout)   :: denskern(pub_cell%num_spins)
    real(kind=DP), intent(out) :: total_energy
    real(kind=DP), intent(out) :: lhxc_energy
    real(kind=DP), intent(out) :: hubbard_energy
    real(kind=DP), intent(out) :: paw_sphere_energies(paw_en_size)
    type(GRID_INFO), intent(inout) :: grid
    real(kind=DP), intent(in) :: localpseudo_fine(grid%ld1,&
         grid%ld2, grid%max_slabs12)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(inout) :: lhxc_fine(grid%ld1,&
         grid%ld2, grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(in) :: ewald_energy
    logical, intent(in) :: ham_update
    logical, intent(in) :: lhxc_fixed
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    logical, optional, intent(in) :: h_xc_write
!CW
    real(8),save        :: dmft_correction
    real(kind=DP)       :: eigen_en(ngwf_basis%num,pub_cell%num_spins)
    integer             :: ressys
!END CW
    ! Local Variables, &
    type(SPAM3), allocatable :: projector_denskern(:)
    integer :: is
    integer :: ierr
!CW
    logical,optional :: spoil_force
    type(SPAM3)      :: dmft_kernel(pub_cell%num_spins),dmft_self(pub_cell%num_spins),dmft_z(pub_cell%num_spins),backupH(pub_cell%num_spins)
    integer          :: godmft
    logical :: partial_dens_present
!END CW

    real(kind=DP), allocatable, dimension(:,:,:) :: h_fine  ! ddor-28feb17
    real(kind=DP), allocatable, dimension(:,:,:,:) :: xc_fine ! ddor-28feb17
    type(SPAM3) :: local_hartree(pub_cell%num_spins) ! ddor-28feb17
    type(SPAM3), allocatable :: local_xc(:) ! ddor-28feb17
    logical :: local_h_xc_write
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: nhat_den_grad
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_fine_partial
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_fine
    real(kind=DP) :: xc_energy_partial
    type(SPAM3) :: denskern_partial(pub_cell%num_spins)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hamiltonian_dens_dep_matrices'
#endif

    local_h_xc_write = .false.
    if (present(h_xc_write)) local_h_xc_write = h_xc_write
    if (local_h_xc_write) then
       allocate(h_fine(pub_fine_grid%ld1, pub_fine_grid%ld2, & ! ddor-28feb17
            pub_fine_grid%max_slabs12),stat=ierr)
       call utils_alloc_check('hamiltonian_dens_dep_matrices','h_fine',ierr)
       allocate(xc_fine(pub_fine_grid%ld1, pub_fine_grid%ld2, & ! ddor-28feb17
            pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_dens_dep_matrices','xc_fine',ierr)
       allocate(local_xc(pub_cell%num_spins),stat=ierr) ! ddor-28feb17
       call utils_alloc_check('hamiltonian_dens_dep_matrices','local_xc',ierr)
       do is=1,pub_cell%num_spins
          call sparse_create(local_xc(is),ham%lhxc(is)) ! ddor-28feb17
          call sparse_create(local_hartree(is),ham%lhxc(is)) ! ddor-28feb17
       enddo
    endif

!CW
  if(pub_dmft_points>0.and.pub_dmft_fully_sc) then
    do is=1,pub_cell%num_spins
       if(.not.ham_update)   call sparse_create(backupH(is), ham%ham(is),iscmplx=.false.)
       if(.not.ham_update)   call   sparse_copy(backupH(is), ham%ham(is))
    enddo
   endif
   if(pub_dmft_points>0.and.pub_dmft_fully_sc) then
    do is=1,pub_cell%num_spins
                             call sparse_create(dmft_self(is)   ,ham%ham(is),iscmplx=.false.)
                             call sparse_create(dmft_z(is)      ,ham%ham(is),iscmplx=.false.)
     if(pub_dmft_fully_sc_h) call sparse_create(dmft_kernel(is),denskern(is),iscmplx=.false.)
    enddo
   endif
   godmft=1   
   10 continue
   dmft_correction=0.d0
!END CW

    ! ndmh: calculate the lhxc potential !
    if (.not.lhxc_fixed) then
!CW
     if(.not.pub_dmft_spoil_kernel.or.present(spoil_force)) then
        if(pub_dmft_spoil_kernel) write(*,*) '... updating LHXC POTENTIALS ...'
!END CW
        if (local_h_xc_write) then
           call hamiltonian_lhxc_calculate(lhxc_fine,lhxc_energy,ham%dijhat, &
             &  grid,localpseudo_fine,core_density_fine,rep%ngwfs_on_grid, &
             &  ngwf_basis,denskern,rep%ngwf_overlap,rep%sp_overlap,    &
             &  add_xc_pot = .true., h_fine = h_fine, xc_fine = xc_fine) ! ddor-28feb17
        else
           call hamiltonian_lhxc_calculate(lhxc_fine,lhxc_energy,ham%dijhat, &
             &  grid,localpseudo_fine,core_density_fine,rep%ngwfs_on_grid, &
             &  ngwf_basis,denskern,rep%ngwf_overlap,rep%sp_overlap)
        endif
!CW
      else
        write(*,*) 'skipping LHXC potentials'        
      endif
!END CW      
    end if

    ! ndmh: only calculate the lhxc matrix if we need to update the Hamiltonian

!CW
   if(ham_update.or.(pub_dmft_points>0.and.pub_dmft_fully_sc))then
!END CW
     ! cks: calculate the lhxc matrix
       do is=1,pub_cell%num_spins
          call integrals_locpot(ham%lhxc(is), &
           &    rep%ngwfs_on_grid,ngwf_basis,rep%ngwfs_on_grid,ngwf_basis, &
           &    grid, lhxc_fine(:,:,:,is))
          if (local_h_xc_write) then ! ddor-28feb17
             call integrals_locpot(local_xc(is), &
                  & rep%ngwfs_on_grid,ngwf_basis,rep%ngwfs_on_grid,ngwf_basis, &
                  & grid, xc_fine(:,:,:,is))
             call restart_kernel_write(local_xc,write_xc=.true.)
          endif
       end do
       if (local_h_xc_write) then ! ddor-28feb17
          call integrals_locpot(local_hartree(1), &
               & rep%ngwfs_on_grid,ngwf_basis,rep%ngwfs_on_grid,ngwf_basis, &
               & grid, h_fine(:,:,:))
          call restart_kernel_write(local_hartree,write_hartree=.true.)
       endif
   end if

!------------------GHB--------------------
!   if 'store_gbooth_partial_kernel1' is present, read this in
    if(pub_my_node_id==1) write(*,*) 'checking if file [x] is present : ','store_gbooth_partial_kernel1' 
    call comms_barrier
    INQUIRE(file='store_gbooth_partial_kernel1',EXIST=partial_dens_present)
    call comms_barrier
    if(.not.partial_dens_present)then
     if(pub_my_node_id==0) write(*,*) 'file is not present. Not calculating subsystem energies.'
    else
     if(pub_my_node_id==0) write(*,*) 'file is present, read from this file'
     ! Allocate denskern_partial
     do is=1,pub_cell%num_spins
        call sparse_create(denskern_partial(is),denskern(is))
     enddo
     call sparse_read_(denskern_partial(1),filename='store_gbooth_partial_kernel1')
     if(pub_my_node_id==0) write(*,*) 'done...'
    endif
    call comms_barrier
    if(pub_cell%num_spins==2) then
        call sparse_copy(denskern_partial(2),denskern_partial(1))
    endif
    !call read_gbooth_partial_kernel
    if(partial_dens_present) then
        ! Calculate partial xc energy.
        ! We have the partial density in denskern_partial
        ! Get it on a grid
        ! %%%%%%%%%%%%%%%PARTIAL DENSITY ON GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        allocate(density_fine_partial(grid%ld1,grid%ld2,grid%max_slabs12, &
             pub_cell%num_spins),stat=ierr)
        call utils_alloc_check('hamiltonian_dens_dep_matrices','density_fine_partial',ierr)
        allocate(density_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
             pub_cell%num_spins),stat=ierr)
        call utils_alloc_check('hamiltonian_dens_dep_matrices','density_fine',ierr)
        call density_on_grid(density_fine_partial, &                      ! output
             grid, denskern_partial, rep%ngwf_overlap, rep%ngwfs_on_grid, ngwf_basis)  ! input
        !Get total density on a grid too!
        call density_on_grid(density_fine, &                      ! output
             grid, denskern, rep%ngwf_overlap, rep%ngwfs_on_grid, ngwf_basis)  ! input
        
        call xc_init_comp()
        xc_energy_partial = 0.0_DP

        ! TODO: Do we want the core_density contribution to the partial
        ! density?!
         ! ndmh: add on core density before calculating xc potential for NLCC
        if (pub_nlcc) then
           do is=1,pub_cell%num_spins
              density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
                   core_density_fine/real(pub_cell%num_spins,kind=DP)
           end do
        end if

        !Now calculate the xc term for this partial density
        ! ndmh: calculate XC energy and potential
        if (pub_aug) then
           allocate(nhat_den_grad(grid%ld1,grid%ld2,grid%max_slabs12, &
                pub_cell%num_spins,0:pub_aug_den_dim),stat=ierr)
           call utils_alloc_check('hamiltonian_lhxc_calculate','nhat_den_grad',ierr)
           nhat_den_grad = 0.0_DP
           call augmentation_density_on_grid(nhat_den_grad,grid,denskern,rep%sp_overlap)
           density_fine(:,:,:,:) = density_fine(:,:,:,:) + nhat_den_grad(:,:,:,:,0)

           call xc_energy_potential_comp(density_fine,density_fine_partial, &
               xc_energy_partial, &
                grid, pub_aug_den_dim, nhat_den_grad)

        else

           ! ndmh: no augmentation density to pass in
           call xc_energy_potential_comp(density_fine,density_fine_partial, &
               xc_energy_partial,grid,0)

        end if
        write(*,*) 'PARTIAL XC ENERGY CONTRIBUTION: ',xc_energy_partial
    endif
!-----------------------------ghb-------------------


    ! qoh: Calculate and include HF exchange if necessary
    if (pub_usehfx) call hf_exchange_calculate(ham%hfexchange,denskern, &
         rep%overlap,rep%ngwfs_on_grid,ngwf_basis,elements,ham%full_vmatrix, &
         .false.)

    ! ndmh: calculate Hubbard energy and Hamiltonian matrix if necessary
    if (pub_hubbard) then
       if (pub_cell%num_spins==1) call sparse_scale(denskern(1),0.5_DP)
       call hubbard_energy_total(hub,hubbard_energy,denskern, &
            hub_proj_basis,rep%hub_overlap,rep%hub_overlap_t)
!CW
       if (ham_update.or.(pub_dmft_points>0.and.pub_dmft_fully_sc)) then
!END CW
          call hubbard_projector_ham(hub,denskern, &
               rep%hub_overlap,rep%hub_overlap_t)
          call hubbard_ham_matrix(hub,ham%hubbard_ham, &
               rep%hub_overlap,rep%hub_overlap_t)
       end if
       if (pub_cell%num_spins==1) call sparse_scale(denskern(1),2.0_DP)
    else
       hubbard_energy = 0.0_DP
    end if

    ! ndmh: calculate nonlocal potential matrix
    if (pub_aug) then

       ! Create temporary matrix for projector denskern
       allocate(projector_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_dens_dep_matrices', &
            'projector_denskern',ierr)
       do is=1,pub_cell%num_spins
          projector_denskern(is)%structure = 'E'
          call sparse_create(projector_denskern(is))
       end do

       ! Create the PAW projector density kernel
       call aug_projector_denskern(projector_denskern,denskern,rep%sp_overlap)

       ! Calculate the PAW nonlocal matrix
       call aug_nonlocal_mat(ham%nonlocpot,ham%dijhat, &
            projector_denskern,rep%sp_overlap,paw_sphere_energies,.false.)

       ! Clean up temporary matrices
       do is=pub_cell%num_spins,1,-1
          call sparse_destroy(projector_denskern(is))
       end do
       deallocate(projector_denskern,stat=ierr)
       call utils_dealloc_check('hamiltonian_dens_dep_matrices', &
            'projector_denskern',ierr)
    end if

    ! Add up local potential, Hartree, exchange-correlation, Ewald, Hubbard and
    ! dispersion energies
    total_energy = lhxc_energy + ewald_energy + hubbard_energy + &
         pub_dispersion_energy

    ! Add up PAW energy terms
    if (pub_paw) then
       total_energy = total_energy &
              - paw_sphere_energies(paw_en_dijhat) &
              - paw_sphere_energies(paw_en_ehart) &
              + paw_sphere_energies(paw_en_exc) &
              - paw_sphere_energies(paw_en_etxc) &
              - paw_sphere_energies(paw_en_dijxc) &
              - paw_sphere_energies(paw_en_exc_core)
    end if

    ! Add kinetic, nonlocal and hfexchange energies
    do is=1,pub_cell%num_spins
       total_energy = total_energy + sparse_trace(denskern(is),rep%kinet)

       ! ndmh: add nonlocal pseudopotential energy if necessary
       if (pub_any_nl_proj.and.(.not.pub_usp)) total_energy = total_energy + &
            sparse_trace(denskern(is),rep%nonlocpot(1))
       if (pub_aug) total_energy = total_energy + &
            sparse_trace(denskern(is),ham%nonlocpot(is))

       ! qoh: include HF exchange if necessary
       if (pub_usehfx) total_energy = total_energy - &
            0.5_DP*sparse_trace(denskern(is),ham%hfexchange(is))

    end do

    ! jd: Include the smeared ion energy corrections, if necessary
    if (pub_is_smeared_ion_rep) then
       total_energy = total_energy + smeared_ion_E_self + smeared_ion_E_smeared
    end if

!CW---------------------------------------------------------------------------!
    if(pub_dmft_points>0.and.pub_dmft_fully_sc.and.godmft==1) then

       call hamiltonian_build_matrix(ham, rep)

       if(pub_my_node_id==0) then 
          ressys=system(" rm store*                   ")
          ressys=system(" rm chem.potential.nmu.iter* ")
       endif
       if(.not.pub_dmft_fully_sc_h)then
            call hubbard_dmft_interface(eigen_en,rep%n_occ,                           &
                                   ham%ham,rep%overlap,rep%inv_overlap,ngwf_basis,    &
                                   hub_proj_basis,hub,rep,elements,denskern,          &
                                   dmft_energy_cor=dmft_correction,dmft_self=dmft_self,dmft_z=dmft_z)
       else
            call hubbard_dmft_interface(eigen_en,rep%n_occ,                           &
                                   ham%ham,rep%overlap,rep%inv_overlap,ngwf_basis,    &
                                   hub_proj_basis,hub,rep,elements,denskern,          &
                                   dmft_energy_cor=dmft_correction,dmft_kernel=dmft_kernel)
       endif
       if(pub_my_node_id==0) then 
          ressys=system(" rm store*                   ")
          ressys=system(" rm chem.potential.nmu.iter* ")
       endif
       if(pub_my_node_id==0) write(*,*) 'DMFT CORRECTION : ', dmft_correction

       ! H -> Z * (H+Sig(oo))

       do is=1,pub_cell%num_spins
          if(ham_update.and..not.pub_dmft_fully_sc_h) then
              call sparse_axpy(ham%ham(is),dmft_self(is),1.d0)
              call sparse_copy(dmft_self(is),ham%ham(is))
              call sparse_product(ham%ham(is),dmft_self(is),dmft_z(is))
          endif
          call sparse_destroy(dmft_z(is))
          call sparse_destroy(dmft_self(is))
       enddo  
       if(pub_dmft_fully_sc_h.and.godmft==1)then
          godmft=2
          goto 10
       endif
    endif

    if(pub_dmft_fully_sc.and.pub_dmft_fully_sc_h)then
      do is=1,pub_cell%num_spins
         call sparse_copy(denskern(is),dmft_kernel(is))
         call sparse_destroy(dmft_kernel(is))
      enddo
    endif

    if(.not.ham_update.and.pub_dmft_points>0.and.pub_dmft_fully_sc) then
      do is=1,pub_cell%num_spins
         call sparse_copy(ham%ham(is),backupH(is))
         call sparse_destroy(backupH(is))
      enddo
    endif

    total_energy = total_energy + dmft_correction
!END CW-----------------------------------------------------------------------!

    if (local_h_xc_write) then
       do is=pub_cell%num_spins,1,-1
          call sparse_destroy(local_xc(is)) ! ddor-28feb17
          call sparse_destroy(local_hartree(is)) ! ddor-28feb17
       enddo
       deallocate(local_xc,stat=ierr) ! ddor-28feb17
       call utils_dealloc_check('hamiltonian_dens_dep_matrices','local_xc',ierr)
       deallocate(xc_fine,stat=ierr) ! ddor-28feb17
       call utils_dealloc_check('hamiltonian_dens_dep_matrices','xc_fine',ierr)
       deallocate(h_fine,stat=ierr) ! ddor-28feb17
       call utils_dealloc_check('hamiltonian_dens_dep_matrices','h_fine',ierr)
    endif
 
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hamiltonian_dens_dep_matrices'
#endif
    contains

     subroutine sparse_read_(mat,filename)
     use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
     implicit none
     type(spam3)                                :: mat
     real(kind=DP), allocatable, dimension(:,:) :: mat_square
     integer                                    :: n1,n2
     character*(*)                              :: filename
        call comms_barrier
        write(6,*) "Get here 1"
        call flush(6)
        n1=NINT(sparse_num_rows(mat))
        write(6,*) "Get here2"
        call flush(6)
        n2=NINT(sparse_num_cols(mat))
        write(6,*) "Get here3"
        call flush(6)
        if(n1==0.and.n2==0)then
          write(*,*) 'SPARSE READ, return, 0-shape matrix'
          return
        endif
        write(6,*) "Get here4"
        call flush(6)
        open(unit=20001,file=filename,form='unformatted')
        write(6,*) "Get here5"
        call flush(6)
        allocate(mat_square(n1,n2),stat=ierr)
        write(6,*) "Get here6"
        call flush(6)
        read(20001) mat_square
        call sparse_convert(mat,mat_square)
        write(6,*) "Get here7"
        call flush(6)
        deallocate(mat_square,stat=ierr)
        close(20001)
     end subroutine

  end subroutine hamiltonian_dens_dep_matrices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_dens_dep_nonsc(ham,rep,ngwf_basis,lhxc_fine, &
       hub,val_rep,val_ham,val_dkn,updated_shift)

    !====================================================================!
    ! Updates the density-matrix dependent parts of the a Hamiltonian    !
    ! when the NGWFs have changed, without recalculating the density     !
    ! (ie non-self-consistently.                                         !
    !--------------------------------------------------------------------!
    ! Written by Laura Ratcliff June 2010 as hamiltonian_cond_ham_update !
    ! Moved to hamiltonian_mod and rearranged by Nicholas Hine in April  !
    ! 2011.                                                              !
    !====================================================================!

    use augmentation, only: aug_projector_denskern, aug_nonlocal_mat
    use cell_grid, only: pub_fine_grid
    use constants, only: paw_en_size
    use function_basis, only: FUNC_BASIS
    use hubbard_build, only: HUBBARD_MODEL, hubbard_ham_matrix
    use integrals, only: integrals_locpot
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_aug, pub_hubbard, pub_usehfx
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout) :: ham
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    type(HUBBARD_MODEL), intent(in) :: hub
    type(NGWF_HAM), intent(in) :: val_ham
    type(NGWF_REP), intent(in) :: val_rep
    type(SPAM3), intent(in) :: val_dkn(pub_cell%num_spins)
    logical, optional, intent(out) :: updated_shift

    ! Local Variables
    type(SPAM3), allocatable :: projector_denskern(:)
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    integer :: ierr
    integer :: is

    ! lr408: calculate the lhxc matrix
    do is=1,pub_cell%num_spins
       call integrals_locpot(ham%lhxc(is), &
            rep%ngwfs_on_grid,ngwf_basis,rep%ngwfs_on_grid,ngwf_basis, &
            pub_fine_grid, lhxc_fine(:,:,:,is))
    end do

    ! ndmh: Calculate and include HF exchange if necessary
    if (pub_usehfx) then
       call utils_abort('Error in hamiltonian_dens_dep_nonsc: Hartree-Fock &
            &not yet supported in conduction calculations')
       !call hf_exchange_calculate(ham%hfexchange,denskern, &
       !     rep%overlap,rep%ngwfs_on_grid,ngwf_basis,.false.)
    end if

    ! ndmh: Calculate and include Hubbard Hamiltonian if necessary
    if (pub_hubbard) then
       call hubbard_ham_matrix(hub,ham%hubbard_ham, &
            rep%hub_overlap,rep%hub_overlap_t)
    end if

    ! ndmh: in PAW, recalculate the nonlocal matrix
    if (pub_aug) then

       ! Create temporary matrix for projector denskern
       allocate(projector_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_dens_dep_nonsc', &
            'projector_denskern',ierr)
       do is=1,pub_cell%num_spins
          projector_denskern(is)%structure = 'E'
          call sparse_create(projector_denskern(is))
       end do

       ! Create the PAW projector density kernel
       call aug_projector_denskern(projector_denskern,val_dkn, &
            val_rep%sp_overlap)
       do is=1,pub_cell%num_spins
          call sparse_scale(projector_denskern(is),pub_cell%spin_fac)
       end do

       ! Calculate the PAW nonlocal matrix
       call aug_nonlocal_mat(ham%nonlocpot,ham%dijhat, &
            projector_denskern,rep%sp_overlap,paw_sphere_energies, &
            show_matrices=.false.)

       ! Clean up temporary matrices
       do is=pub_cell%num_spins,1,-1
          call sparse_destroy(projector_denskern(is))
       end do
       deallocate(projector_denskern,stat=ierr)
       call utils_dealloc_check('hamiltonian_dens_dep_nonsc', &
            'projector_denskern',ierr)
    end if

    ! lr408: Build unprojected conduction hamiltonian
    call hamiltonian_build_matrix(ham, rep)

    ! lr408: Finally build projected matrix
    if (present(updated_shift)) then
       call hamiltonian_proj_cond_matrix(rep, ham, &
            val_rep, val_ham, val_dkn, shift_changed=updated_shift)
    end if

  end subroutine hamiltonian_dens_dep_nonsc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_lhxc_calculate(lhxc_fine, lhxc_energy, dijhat, &
       grid, localpseudo_fine, core_density_fine, ngwfs_on_grid, ngwf_basis, &
       denskern, overlap, sp_overlap, add_xc_pot, h_fine, xc_fine) ! ddor-28feb17

    !=====================================================================!
    ! This subroutine calculates and returns the sum of the               !
    ! local pseudopotential, hartree potential and exchange-correlation   !
    ! potentials and their energy.                                        !
    !---------------------------------------------------------------------!
    !                                                                     !
    !                                                                     !
    !---------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 18/7/2001.                      !
    ! Modified by Chris-Kriton Skylaris on 14/02/2004 so that it          !
    ! is less memory-hungry.                                              !
    ! Modified by Peter Haynes on 1/7/2004 to use fourier parallelisation.!
    ! Modified by Nicholas Hine on 6/2/2009 for NLCC core charges         !
    ! Modified by David D. O'Regan on 6/5/2009 for TDDFT                  !
    ! Modified by Nicholas Hine in July 2009 for function basis type and  !
    ! for SPAM3                                                           !
    ! Modified by Nicholas Hine in July 2010 for PAW                      !
    ! Modified by Nicholas Hine in February 2011 for augmentation_mod     !
    ! Modified by Nicholas Hine in March 2011 to make logic clearer and   !
    ! easier to read.                                                     !
    !=====================================================================!

    use augmentation, only: augmentation_density_on_grid, &
         augmentation_screen_dij
    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
    use cutoff_coulomb, only: cutoff_coulomb_hartree
    use density, only: density_on_grid
    use enthalpy, only: enthalpy_terms
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_coulomb_cutoff, pub_nlcc, pub_paw, &
         pub_multigrid_hartree, pub_aug_den_dim, pub_aug, pub_usp, &
         pub_nhat_in_xc, pub_external_pressure
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc, only: xc_energy_potential

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(SPAM3), intent(in)      :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in)      :: overlap
    type(SPAM3), intent(in)      :: sp_overlap
    type(GRID_INFO), intent(inout)  :: grid
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in) :: localpseudo_fine(grid%ld1,&
         grid%ld2, grid%max_slabs12)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(out) :: lhxc_energy
    real(kind=DP), intent(out) :: lhxc_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12, pub_cell%num_spins)
    type(SPAM3), intent(inout)   :: dijhat(:)
    logical, intent(in), optional :: add_xc_pot
    real(kind=DP), intent(out), optional :: h_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12) ! ddor-28feb17
    real(kind=DP), intent(out), optional :: xc_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12, pub_cell%num_spins) ! ddor-28feb17

    ! Local Variables
    integer :: ierr   ! error flag
    integer :: is
    logical :: loc_add_xc_pot
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: pot_fine
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: nhat_den_grad
    real(kind=DP) :: xc_energy
    real(kind=DP) :: enthalpy_energy
    real(kind=DP) :: locpspot_energy   ! ndmh: energy of density in local pspot
    real(kind=DP) :: hartree_energy    ! jd: Multigrid Hartree energy
    real(kind=DP) :: cavitation_energy ! jd: Implicit solvent cavitation energy

    ! Optional argument to suppress xc potential being added to lhxc_fine
    if (present(add_xc_pot)) then
       loc_add_xc_pot = add_xc_pot
    else
       loc_add_xc_pot = .true.
    end if

    ! Allocate workspace
    allocate(density_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hamiltonian_lhxc_calculate','density_fine',ierr)
    allocate(pot_fine(grid%ld1,grid%ld2,grid%max_slabs12,pub_cell%num_spins), &
         stat=ierr)
    call utils_alloc_check('hamiltonian_lhxc_calculate','pot_fine',ierr)
    if (pub_aug) then
       allocate(nhat_den_grad(grid%ld1,grid%ld2,grid%max_slabs12, &
            pub_cell%num_spins,0:pub_aug_den_dim),stat=ierr)
       call utils_alloc_check('hamiltonian_lhxc_calculate','nhat_den_grad',ierr)
    end if

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY ON GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call density_on_grid(density_fine, &                      ! output
         grid, denskern, overlap, ngwfs_on_grid, ngwf_basis)  ! input

    ! ndmh: for PAW, get augmentation density nhat
    if (pub_aug) then
       nhat_den_grad = 0.0_DP
       call augmentation_density_on_grid(nhat_den_grad,grid,denskern,sp_overlap)
       density_fine(:,:,:,:) = density_fine(:,:,:,:) + nhat_den_grad(:,:,:,:,0)
    end if
    ! %%%%%%%%%%%%%%%%%%%%%%%%%% END DENSITY ON GRID %%%%%%%%%%%%%%%%%%%%%%%%%%



    ! =============== LOCAL PSEUDOPOTENTIAL ENERGY AND POTENTIAL ===============

    locpspot_energy = 0.0_DP
    do is=1,pub_cell%num_spins
       locpspot_energy = locpspot_energy + integrals_product_on_grid( &
            grid,localpseudo_fine(:,:,:),density_fine(:,:,:,is))
       lhxc_fine(:,:,:,is) = localpseudo_fine(:,:,:)
    end do

    ! ============= END LOCAL PSEUDOPOTENTIAL ENERGY AND POTENTIAL =============


    ! ====================== HARTREE ENERGY AND POTENTIAL ======================

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
    !     derivative. The Hartree energy then needs to be
    !     calculated differently, and hartree_via_multigrid() takes care of
    !     that, as well as calculating the cavitation energy.
    hartree_energy = 0.0_DP
    cavitation_energy = 0.0_DP
    if (pub_coulomb_cutoff) then
       call cutoff_coulomb_hartree(pot_fine, density_fine)              ! (2)
    else if (pub_multigrid_hartree) then
       call hartree_via_multigrid(pot_fine, density_fine, &             ! (3)
            hartree_energy, cavitation_energy)
    else
       call hartree_on_grid(pot_fine, density_fine, grid)               ! (1)
    end if
    if (.not.pub_multigrid_hartree) then
       do is=1,pub_cell%num_spins
          hartree_energy = hartree_energy + 0.5_DP*integrals_product_on_grid( &
               grid, pot_fine(:,:,:,is), density_fine(:,:,:,is))
       end do
    end if

    do is=1,pub_cell%num_spins
       lhxc_fine(:,:,:,is) = lhxc_fine(:,:,:,is) + pot_fine(:,:,:,is)
    end do
    if(present(h_fine)) h_fine(:,:,:) = pot_fine(:,:,:,1) ! ddor-28feb17

    ! ==================== END HARTREE ENERGY AND POTENTIAL ====================

    ! ndmh: screen nonlocal projector energies now if the augmentation density
    ! ndmh: is not going to be included in the exchange energy calculation
    if (pub_paw.and.(.not.pub_nhat_in_xc)) then
       call augmentation_screen_dij(dijhat,lhxc_fine,grid)
    end if

    ! ========================= XC ENERGY AND POTENTIAL ========================

    ! ndmh: add on core density before calculating xc potential for NLCC
    if (pub_nlcc) then
       do is=1,pub_cell%num_spins
          density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
               core_density_fine/real(pub_cell%num_spins,kind=DP)
       end do
    end if

    ! ndmh: calculate XC energy and potential
    if (pub_aug) then

       ! ndmh: if we do not want the nhat density in the XC calculation,
       ! ndmh: remove it now and calculate xc potential without it
       if (pub_paw) then
          if (.not.pub_nhat_in_xc) then
             density_fine = density_fine - nhat_den_grad(:,:,:,:,0)
          end if
       end if
       call xc_energy_potential(density_fine, xc_energy, pot_fine, &
            grid, pub_aug_den_dim, nhat_den_grad)

    else

       ! ndmh: no augmentation density to pass in
       call xc_energy_potential(density_fine, xc_energy, pot_fine, grid, 0)

    end if

    if (loc_add_xc_pot) then
       do is=1,pub_cell%num_spins
          lhxc_fine(:,:,:,is) = lhxc_fine(:,:,:,is) + pot_fine(:,:,:,is)
       if(present(h_fine)) xc_fine(:,:,:,is) = pot_fine(:,:,:,is) ! ddor-28feb17
       end do
    end if

    ! ====================== END XC ENERGY AND POTENTIAL =======================

    ! ndmh: screen nonlocal projector energies now if the augmentation density
    ! ndmh: is going to be included in the exchange energy calculation
    if (pub_paw.and.pub_nhat_in_xc) then
       call augmentation_screen_dij(dijhat,lhxc_fine,grid)
    end if

    !==========================ELECTRONIC ENTHALPY==============================

    enthalpy_energy = 0.0_DP
    if (pub_external_pressure>0.0_DP) then
       call enthalpy_terms(pot_fine,enthalpy_energy,density_fine,grid)
       lhxc_fine(:,:,:,:) = lhxc_fine(:,:,:,:) + pot_fine(:,:,:,:)
    end if

    ! ADD UP ENERGIES
    lhxc_energy = xc_energy + locpspot_energy + hartree_energy + &
         cavitation_energy + enthalpy_energy

    ! Deallocate workspace
    if (pub_aug) then
       deallocate(nhat_den_grad,stat=ierr)
       call utils_dealloc_check('hamiltonian_lhxc_calculate','nhat_den_grad', &
            ierr)
    end if
    deallocate(pot_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_lhxc_calculate','pot_fine',ierr)
    deallocate(density_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_lhxc_calculate','density_fine',ierr)

  end subroutine hamiltonian_lhxc_calculate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine hamiltonian_build_matrix(ham,rep)
    !========================================================================!
    ! This subroutine builds the full Hamiltonian matrix in SPAM3 format out !
    ! of its various components, which have already been calculated as SPAM3 !
    ! matrices.                                                              !
    !------------------------------------------------------------------------!
    ! ham        (inout) : Hamiltonian Wrapper (contains lhxc,hubbard_ham,   !
    !                      hfexchange and PAW nonlocpot)                     !
    ! rep           (in) : NGWF Representation (contains kinet and normal    !
    !                      nonlocpot)                                        !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine in April 2010.                                !
    !========================================================================!

    use constants, only: DP
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_any_nl_proj, pub_usehfx, pub_hubbard, pub_paw, &
         pub_aug, pub_usp
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_copy
    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout) :: ham
    type(NGWF_REP), intent(in) :: rep

    ! Local variable
    integer :: is    ! spin loop counter

    do is=1,pub_cell%num_spins

       ! ndmh: copy lhxc matrix to Hamiltonian
       call sparse_copy(ham%ham(is),ham%lhxc(is))

       ! nmdh: add kinetic energy matrix
       call sparse_axpy(ham%ham(is),rep%kinet,1.0_DP)

       ! ndmh: add the nonlocal potential matrix if necessary
       if (pub_any_nl_proj.and.(.not.pub_usp)) &
            call sparse_axpy(ham%ham(is),rep%nonlocpot(1),1.0_DP)
       if (pub_aug) &
            call sparse_axpy(ham%ham(is),ham%nonlocpot(is),1.0_DP)

       ! qoh: add the hfexchange matrix if necessary
       if (pub_usehfx) call sparse_axpy(ham%ham(is),ham%hfexchange(is),-1.0_DP)

       ! ddor: add the hubbard hamiltonian if necessary
       if (pub_hubbard) call sparse_axpy(ham%ham(is),ham%hubbard_ham(is),1.0_DP)

    end do

  end subroutine hamiltonian_build_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_proj_cond_matrix(cond_rep, cond_ham, &
       val_rep, val_ham, val_dkn, cond_dkn, shift_changed)

    !=========================================================================!
    ! This subroutine calculates the projected conduction Hamiltonian matrix  !
    ! for a given shift w such that:                                          !
    ! H_proj = H_c - wS_c - T+K_vH_vK_vT + wT+K_vS_vK_vT                      !
    ! See L.E. Ratcliff, N.D.M. Hine and P.D. Haynes, Phys. Rev. B 84, 165131 !
    ! (2011) for details, specifically Sections II.A to II.D                  !
    !-------------------------------------------------------------------------!
    ! cond_rep       (in) : conduction NGWF representation                    !
    ! val_rep        (in) : valence NGWF representation                       !
    ! cond_ham    (inout) : conduction Hamiltonian matrices                   !
    ! val_ham        (in) : conduction Hamiltonian matrices                   !
    ! val_kernel     (in) : valence density kernel                            !
    ! shift_changed (out) : optional logical output for changing shifts       !
    !-------------------------------------------------------------------------!
    ! Written by Laura Ratcliff in June 2010.                                 !
    ! Modified by Laura Ratcliff in October 2010 to include NGWF_REP and      !
    ! NGWF_HAM types.                                                         !
    ! Modified by Nicholas Hine in January 2012 to add comments and ensure    !
    ! proper working with spin polarisation.                                  !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE, max_spins
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_output_detail, cond_fixed_shift, cond_calc_max_eigen, &
         cond_shift_buffer
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_copy, sparse_create, &
         sparse_destroy, sparse_extremal_eigenvalue, sparse_product, &
         sparse_transpose, sparse_trace, sparse_transpose_structure

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in)     :: cond_rep
    type(NGWF_REP), intent(in)     :: val_rep
    type(NGWF_HAM), intent(inout)  :: cond_ham
    type(NGWF_HAM), intent(in)     :: val_ham
    type(SPAM3), intent(in)        :: val_dkn(pub_cell%num_spins)
    type(SPAM3), optional, intent(in)        :: cond_dkn(pub_cell%num_spins)
    logical, optional, intent(out) :: shift_changed

    ! Local variables
    integer :: is    ! spin loop counter
    real(kind=DP),parameter :: delta_e_thresh =1.0E-12_DP ! energy/atom threshold
    real(kind=DP) :: max_en(max_spins)      ! maximum hamiltonian eigenvalue

    type(SPAM3) :: sinv_ham         ! inverse overlap times hamiltonian
    type(SPAM3) :: trans_cross_olap ! the transpose of the cross overlap
    type(SPAM3) :: uk, kt, ukh, uks, ukhkt, ukskt ! matrix products

#ifdef DEBUG
    real(kind=DP) :: debug_trace    ! matrix trace for debugging purposes
#endif


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: &
         &Entering hamiltonian_proj_cond_matrix'
#endif

    if (present(shift_changed)) shift_changed = .false.

    ! Create transpose of cross_overlap matrix
    call sparse_transpose_structure(trans_cross_olap%structure, &
         cond_rep%cross_overlap)
    call sparse_create(trans_cross_olap)
    call sparse_transpose(trans_cross_olap,cond_rep%cross_overlap)

    do is=1,pub_cell%num_spins

       ! lr408: copy unprojected ham to cond_non_proj_ham
       call sparse_copy(cond_ham%cond_non_proj_ham(is),cond_ham%ham(is))

       ! lr408: If required, calculate maximum eigenvalue of the conduction
       ! lr408: Hamiltonian
       if (cond_calc_max_eigen) then

          call sparse_create(sinv_ham, cond_rep%inv_overlap, &
               cond_ham%cond_non_proj_ham(is))

          ! lr408: S^-1.H
          call sparse_product(sinv_ham, cond_rep%inv_overlap, &
               cond_ham%cond_non_proj_ham(is))

          ! lr408: maximum orbital energy
          call sparse_extremal_eigenvalue(sinv_ham, cond_rep%overlap, &
               max_en(is), delta_e_thresh)

          if ((pub_output_detail == VERBOSE).and.(pub_on_root)) then
             if (pub_cell%num_spins==2) then
                write(stdout,'(a,2f20.14,f20.14)') &
                     'Conduction Hamiltonian max eigenvalues:', max_en(1:2), &
                     cond_ham%cond_shift
             else
                write(stdout,'(a,f20.14,f20.14)') &
                     'Conduction Hamiltonian max eigenvalue:', max_en(1), &
                     cond_ham%cond_shift
             end if
          end if

          call sparse_destroy(sinv_ham)

       end if

    end do

    ! lr408: If the shift is not being kept constant and the maximum
    ! lr408: eigenvalue has gone above the current shift, increase
    ! lr408: the shift accordingly
    if (.not. cond_fixed_shift) then
       if (maxval(max_en) > cond_ham%cond_shift) then
          cond_ham%cond_shift = maxval(max_en) + cond_shift_buffer
          if (pub_output_detail == VERBOSE) then
             if (pub_on_root) write(stdout,'(a,f24.14)') &
                  'Conduction shift updated:', cond_ham%cond_shift
          end if
          if (present(shift_changed)) shift_changed = .true.
       end if
    end if

    ! Allocate storage for matrices required for projection
    call sparse_create(uk,trans_cross_olap,val_dkn(1))
    call sparse_create(kt,val_dkn(1),cond_rep%cross_overlap)
    call sparse_create(ukh,uk,val_ham%ham(1))
    call sparse_create(uks,uk,val_rep%overlap)
    call sparse_create(ukhkt,ukh,kt)
    call sparse_create(ukskt,uks,kt)

    ! Now project out valence Hamiltonian and shift valence states up in energy
    ! (see Eq.11 of PRB 84 165131 (2011))
    do is=1,pub_cell%num_spins

       ! Calculate ukhkt and ukhks
       call sparse_product(uk,trans_cross_olap,val_dkn(is))
       call sparse_product(kt,val_dkn(is),cond_rep%cross_overlap)
       call sparse_product(ukh,uk,val_ham%ham(is))
       call sparse_product(uks,uk,val_rep%overlap)
       call sparse_product(ukhkt,ukh,kt)
       call sparse_product(ukskt,uks,kt)

       ! lr408: add -TKHKT to conduction hamiltonian
       call sparse_axpy(cond_ham%ham(is),ukhkt,-1.0_dp)

       ! lr408: add wTKSKT to conduction hamiltonian
       call sparse_axpy(cond_ham%ham(is),ukskt,cond_ham%cond_shift)

#ifdef DEBUG
       ! lr408: Calculate Tr[T+KSKTM] as a measure of orthogonality with valence
       ! lr408: states
       if (present(cond_dkn)) then
          debug_trace = sparse_trace(ukskt,cond_dkn(is))
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
               &Tr[T+KS_vKTM] = ',debug_trace
          debug_trace = sparse_trace(ukskt,cond_rep%inv_overlap)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: &
               &Tr[T+KS_vKTS^-1] = ',debug_trace
       end if
#endif

    end do

    ! Clean up temporary matrices
    call sparse_destroy(ukskt)
    call sparse_destroy(ukhkt)
    call sparse_destroy(uks)
    call sparse_destroy(ukh)
    call sparse_destroy(kt)
    call sparse_destroy(uk)

    call sparse_destroy(trans_cross_olap)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving hamiltonian_proj_cond_matrix'
#endif


  end subroutine hamiltonian_proj_cond_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hamiltonian_energy_components( &
       pur_denskern, rep, localpseudo_fine, core_density_fine, &   ! input
       ngwf_basis, hub_proj_basis, hub, ewald_energy, hfexchange)  ! input

    !===================================================================!
    ! This subroutine calculates and prints the components that make    !
    ! up the total energy.                                              !
    !-------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 17/4/2001 as                  !
    ! hamiltonian_density_ded_debug and hamiltonian_energy_terms.       !
    ! Modified by Peter Haynes 1/7/2004 for Fourier parallelisation.    !
    ! Modified by Chris-Kriton Skylaris on 7/10/2004 to return the      !
    ! Hartree energy as calculated directly from the charge density.    !
    ! Renamed to hamiltonian_energy_components and rewritten by         !
    ! Chris-Kriton Skylaris on 8/10/2004 in order to print more details !
    ! using less memory.                                                !
    ! Converted to SPAM2 by Peter Haynes, July 2006                     !
    ! Added HF exchange (optional) by Quintin Hill 23/10/2008.          !
    ! Modified by Nicholas Hine on 6/2/2009 for NLCC core charges       !
    ! Modified for DFT+U by D. D. O'Regan in April 2009                 !
    ! Modified by Nicholas Hine in July 2009 for function basis type    !
    !===================================================================!

    use augmentation, only: augmentation_density_on_grid, &
         aug_projector_denskern, augmentation_screen_dij, aug_nonlocal_mat
    use cell_grid, only: pub_fine_grid
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, UP, DN, stdout, max_spins, VERBOSE, paw_en_size, &
         paw_en_dij0, paw_en_ehart, paw_en_exc, paw_en_etxc, paw_en_exc_core
    use cutoff_coulomb, only: cutoff_coulomb_hartree
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use hubbard_build, only: HUBBARD_MODEL, hubbard_energy_total
    use integrals, only: integrals_trace_on_grid, integrals_product_on_grid, &
         integrals_locpot
    use is_smeared_ions, only: smeared_ion_E_self, &
         smeared_ion_E_smeared
    use ngwf_representation, only: NGWF_REP
    use rundat, only: pub_output_detail, pub_dispersion, pub_usehfx, &
         pub_coulomb_cutoff, pub_hubbard, pub_spin, pub_any_nl_proj, &
         pub_ii_energy_direct, pub_nlcc, pub_paw, pub_open_localpseudo, &
         pub_is_implicit_solvent, pub_is_smeared_ion_rep, pub_aug, pub_usp, &
         pub_is_include_cavitation, pub_multigrid_hartree, pub_aug_den_dim, &
         pub_nhat_in_xc, pub_cdft
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_trace, sparse_create, &
         sparse_destroy, sparse_scale
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use vdwcorrection, only: pub_dispersion_energy
    use xc, only: xc_energy_potential

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(SPAM3), intent(inout) :: pur_denskern(pub_cell%num_spins)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(in) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP), intent(in) :: ewald_energy
    type(SPAM3), intent(in) :: hfexchange(:) ! Hartree-Fock exchange matrix
    ! qoh: Warning: hfexchange is only allocated if necessary so may not be
    ! qoh: allocated.

    ! Local Variables
    integer :: is
    integer :: ierr
    integer :: i1,i2,islab12                       ! Grid counters
    real(kind=DP), allocatable, dimension(:,:,:,:) :: density_fine
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: nhat_den_grad
    real(kind=DP), allocatable, dimension(:,:,:,:) :: buffer_fine
    real(kind=DP), allocatable, dimension(:,:,:,:) :: lhxc_fine
    real(kind=DP) :: integrated_density
    real(kind=DP) :: integrated_ps_val_density
    real(kind=DP) :: integrated_spin
    real(kind=DP) :: integrated_field(max_spins)
    real(kind=DP) :: integrated_nhat(max_spins)
    real(kind=DP) :: tr_integrated_density
    real(kind=DP) :: tr_integrated_spin
    real(kind=DP) :: integrated_mod_spin
    real(kind=DP) :: integrated_depleted_spin
    real(kind=DP) :: tr_field(max_spins)
    real(kind=DP) :: hartree_energy
    real(kind=DP) :: nhat_hartree_energy
    real(kind=DP) :: tr_hartree_energy
    real(kind=DP) :: tr_nhat_hartree_energy
    real(kind=DP) :: locpot_energy
    real(kind=DP) :: nhat_locpot_energy
    real(kind=DP) :: tr_locpot_energy
    real(kind=DP) :: tr_nhat_locpot_energy
    real(kind=DP) :: kinetic_energy
    real(kind=DP) :: nonlocpot_energy
    real(kind=DP) :: xc_energy
    real(kind=DP) :: xc_denpot_energy
    real(kind=DP) :: xc_nhat_denpot_energy
    real(kind=DP) :: hfx_energy ! Hartree-Fock exchange energy
    real(kind=DP) :: hubbard_energy ! DFT+U correction energy
    real(kind=DP) :: spin_squared_expectation
    real(kind=DP) :: factor
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: total_energy
    real(kind=DP) :: total_sphere_energy
    type(SPAM3) :: localmat
    type(SPAM3), allocatable :: projector_denskern(:)
    type(SPAM3), allocatable :: local_nonlocpot(:)
    type(SPAM3), allocatable :: local_dijhat(:)
    real(kind=DP) :: smeared_ion_self_energy    ! jd: Extra terms resulting from
    real(kind=DP) :: smeared_ion_nonself_energy ! the smeared ion representation
    real(kind=DP) :: cavitation_energy ! jd: implicit solvent cavitation energy

    ! Start Timer
    call timer_clock('hamiltonian_energy_components',1)

    ! --------------------------------------------------------------------------

    ! Allocate density workspace
    allocate(density_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins), &
         stat=ierr)
    call utils_alloc_check('hamiltonian_energy_components','density_fine',ierr)
    if (pub_aug) then
       allocate(nhat_den_grad(pub_fine_grid%ld1, pub_fine_grid%ld2, &
            pub_fine_grid%max_slabs12, pub_cell%num_spins,0:pub_aug_den_dim), &
            stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'nhat_den_grad',ierr)
    end if

    ! ndmh: allocate SPAM3 matrices for projector density kernel and local
    ! ndmh: workspace with dimensions of PAW projectors, and create projector
    ! ndmh: density kernel
    if (pub_aug) then
       allocate(projector_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'projector_denskern',ierr)
       allocate(local_nonlocpot(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'local_nonlocpot',ierr)
       allocate(local_dijhat(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components', &
            'local_dijhat',ierr)
       do is=1,pub_cell%num_spins
          projector_denskern(is)%structure='E'
          call sparse_create(projector_denskern(is))
          call sparse_create(local_dijhat(is),projector_denskern(is))
       end do
       call aug_projector_denskern(projector_denskern,pur_denskern, &
            rep%sp_overlap)
    end if

    ! ################## CHARGE DENSITY #############################

    ! cks: initialise density_fine
    density_fine = 0.0_DP

    call density_on_grid(density_fine, &                  ! output
         pub_fine_grid, pur_denskern, rep%ngwf_overlap, & ! input
         rep%ngwfs_on_grid, ngwf_basis)                   ! input

    if (pub_aug) then
       nhat_den_grad = 0.0_DP
       call augmentation_density_on_grid(nhat_den_grad,pub_fine_grid, &
            pur_denskern,rep%sp_overlap)
    end if

    ! cks: find number of electrons from valence charge density
    do is=1,pub_cell%num_spins
       integrated_field(is) = &
            integrals_trace_on_grid(density_fine(:,:,:,is),pub_fine_grid)
    end do
    ! ndmh: find number of electrons from compensation density
    if (pub_aug) then
       do is=1,pub_cell%num_spins
          integrated_nhat(is) = &
               integrals_trace_on_grid(nhat_den_grad(:,:,:,is,0), &
               pub_fine_grid)
       end do
    else
       integrated_nhat(:) = 0.0_DP
    end if
    if (pub_cell%num_spins == 1) then
       integrated_ps_val_density = integrated_field(1)
       integrated_density = integrated_field(1) + integrated_nhat(1)
       integrated_spin = 0.0_DP
    else
       integrated_density = integrated_field(UP) + integrated_field(DN)
       integrated_ps_val_density = integrated_density
       integrated_density = integrated_density + &
            integrated_nhat(UP) + integrated_nhat(DN)
       integrated_spin = integrated_field(UP) - integrated_field(DN)
       integrated_spin = integrated_spin + &
            integrated_nhat(UP) - integrated_nhat(DN)
       ! pdh: calculate integrated |spin density|
       integrated_mod_spin = 0.0_DP
       integrated_depleted_spin = 0.0_DP
       do islab12=1,pub_fine_grid%num_my_slabs12
          do i2=1,pub_fine_grid%n2
             do i1=1,pub_fine_grid%n1
                integrated_mod_spin = integrated_mod_spin + &
                     abs(density_fine(i1,i2,islab12,UP) - &
                     density_fine(i1,i2,islab12,DN))
                !ddor: Estimate <S^2>
                if (pub_output_detail == VERBOSE) then
                   integrated_depleted_spin = &
                        integrated_depleted_spin + &
                        DIM(density_fine(i1,i2,islab12,DN),&
                        &density_fine(i1,i2,islab12,UP))
                endif
             end do
          end do
       end do
       call comms_reduce('SUM',integrated_mod_spin)
       integrated_mod_spin = integrated_mod_spin * pub_fine_grid%weight
    end if

    if (pub_output_detail == VERBOSE) then
       do is=1,pub_cell%num_spins
          tr_field(is) = sparse_trace(pur_denskern(is),rep%overlap)
       end do
       if (pub_cell%num_spins == 1) then
          tr_integrated_density = tr_field(1)
          tr_integrated_spin = 0.0_DP
       else
          tr_integrated_density = tr_field(UP) + tr_field(DN)
          tr_integrated_spin = tr_field(UP) - tr_field(DN)
          ! ddor: An estimate of the expectation value of the S^2 operator, used
          !       to correct calculations of splitting between spin-states where
          !       spin contamination has taken place.
          !       Based on Local Density Approximation and derived in
          !       Wang, Becke, Smith. J. Chem. Phys 102 (8) 1995.
          !       The integral of the magnetisation density in the down-spin
          !       excessive areas is subtracted from S(S+1)
          call comms_reduce('SUM',integrated_depleted_spin)
          integrated_depleted_spin = integrated_depleted_spin * &
               pub_fine_grid%weight
          spin_squared_expectation = (0.25_DP * REAL(pub_spin,kind=DP) * &
               &(REAL(pub_spin,kind=DP) + 2.0_DP) ) + integrated_depleted_spin
       end if
    end if

    ! ############### END CHARGE DENSITY ############################


    ! *********** LOCAL PSEUDOPOTENTIAL ENERGIES **********************
    locpot_energy = 0.0_DP
    nhat_locpot_energy = 0.0_DP
    do is=1,pub_cell%num_spins
       locpot_energy = locpot_energy + &
            integrals_product_on_grid(pub_fine_grid,localpseudo_fine, &
            density_fine(:,:,:,is))
       if (pub_aug) then
          nhat_locpot_energy = nhat_locpot_energy + &
               integrals_product_on_grid(pub_fine_grid,localpseudo_fine, &
               nhat_den_grad(:,:,:,is,0))
       end if
    end do
    if (pub_aug) locpot_energy = locpot_energy + nhat_locpot_energy

    if (pub_output_detail == VERBOSE) then
       call sparse_create(localmat,rep%overlap)
       call integrals_locpot(localmat, &                            ! output
            rep%ngwfs_on_grid, ngwf_basis, &                        ! input
            rep%ngwfs_on_grid, ngwf_basis, &                        ! input
            pub_fine_grid, localpseudo_fine)                                       ! input
       tr_locpot_energy = 0.0_DP
       do is=1,pub_cell%num_spins
          tr_locpot_energy = tr_locpot_energy + &
               sparse_trace(pur_denskern(is), localmat)
       end do
       call sparse_destroy(localmat)
       ! ndmh: Add on augmentation density contribution:
       ! ndmh: \int v_ion \hat{n} dr
       if (pub_aug) then
          call augmentation_screen_dij(local_dijhat(1:1),localpseudo_fine, &
               pub_fine_grid)
          tr_nhat_locpot_energy = 0.0_DP
          do is=1,pub_cell%num_spins
             tr_nhat_locpot_energy = tr_nhat_locpot_energy + &
                  sparse_trace(projector_denskern(is),local_dijhat(1))
          end do
          tr_locpot_energy = tr_locpot_energy + tr_nhat_locpot_energy
       end if
    end if
    ! ********END LOCAL PSEUDOPOTENTIAL ENERGIES **********************


    ! cks: allocate simulation cell fine workspace array
    allocate(buffer_fine(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12, pub_cell%num_spins), stat=ierr)
    call utils_alloc_check('hamiltonian_energy_components','buffer_fine',ierr)
    if (pub_aug) then
       allocate(lhxc_fine(pub_fine_grid%ld1, pub_fine_grid%ld2, &
            pub_fine_grid%max_slabs12, pub_cell%num_spins), stat=ierr)
       call utils_alloc_check('hamiltonian_energy_components','lhxc_fine',ierr)
       do is=1,pub_cell%num_spins
          lhxc_fine(:,:,:,is) = localpseudo_fine
       end do
    end if


    ! %%%%%%%%%%%%% HARTREE ENERGIES %%%%%%%%%%%%%%%%%%%%%%%
    ! cks: initialise hartree_fine
    buffer_fine = 0.0_DP

    ! ndmh: Add on compensation density before Hartree calculation
    if (pub_aug) then
       density_fine = density_fine + nhat_den_grad(:,:,:,:,0)
    end if

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
    !     derivative. The Hartree energy then needs to be
    !     calculated differently, and hartree_via_multigrid() takes care of
    !     that, as well as calculating the cavitation energy.
    hartree_energy = 0.0_DP
    cavitation_energy = 0.0_DP
    if (pub_coulomb_cutoff) then
       call cutoff_coulomb_hartree(buffer_fine, density_fine)              ! (2)
    else if (pub_multigrid_hartree) then
       call hartree_via_multigrid(buffer_fine, density_fine, &             ! (3)
            hartree_energy, cavitation_energy)
    else
       call hartree_on_grid(buffer_fine, density_fine, pub_fine_grid)      ! (1)
    end if

    ! jd: If working with the multigrid Hartree, hartree_energy is already
    !     calculated, no need to calculate hartree_energy directly.
    if (.not.pub_multigrid_hartree) then
       do is=1,pub_cell%num_spins
          hartree_energy = hartree_energy + 0.5_DP*integrals_product_on_grid( &
               pub_fine_grid, buffer_fine(:,:,:,is), density_fine(:,:,:,is))
       end do
    end if

    nhat_hartree_energy = 0.0_DP
    if (pub_aug) then

       ! ndmh: add the Hartree potential to the LHXC potential
       lhxc_fine = lhxc_fine + buffer_fine

       ! ndmh: subract off the compensation density to get just the band
       ! ndmh: energy part of the Hartree energy
       density_fine = density_fine - nhat_den_grad(:,:,:,:,0)

       ! ndmh: calculate the energy of the compensation density in the Hartree
       ! ndmh: potential of the total soft density
       nhat_hartree_energy = 0.0_DP
       do is=1,pub_cell%num_spins
          nhat_hartree_energy = nhat_hartree_energy + 0.5_DP* &
               integrals_product_on_grid(pub_fine_grid, buffer_fine(:,:,:,is), &
               density_fine(:,:,:,is))
       end do

    end if

    if (pub_output_detail == VERBOSE) then
       call sparse_create(localmat, rep%overlap)
       tr_hartree_energy = 0.0_DP
       do is=1,pub_cell%num_spins
          call integrals_locpot(localmat, &                            ! output
               rep%ngwfs_on_grid, ngwf_basis, &                        ! input
               rep%ngwfs_on_grid, ngwf_basis, &                        ! input
               pub_fine_grid, buffer_fine(:,:,:,is))                   ! input
          tr_hartree_energy = tr_hartree_energy + &
               0.5_DP * sparse_trace(pur_denskern(is), localmat)
       end do
       call sparse_destroy(localmat)
       ! ndmh: Calculate energy of augmentation density in Hartree potential
       ! ndmh: (1/2) \int v_H[\tilde{n}+\hat{n}] \hat{n} dr
       if (pub_aug) then
          do is=1,pub_cell%num_spins
             call sparse_scale(local_dijhat(is),0.0_DP)
          end do
          call augmentation_screen_dij(local_dijhat,buffer_fine,pub_fine_grid)
          tr_nhat_hartree_energy = 0.0_DP
          do is=1,pub_cell%num_spins
             tr_nhat_hartree_energy = tr_nhat_hartree_energy + &
                  0.5_DP * sparse_trace(projector_denskern(is), &
                  local_dijhat(is))
          end do
          tr_hartree_energy = tr_hartree_energy + tr_nhat_hartree_energy
       end if
    end if
    ! %%%%%%%%%% END HARTREE ENERGIES %%%%%%%%%%%%%%%%%%%%%%%

    ! ============== XC MATRIX & ENERGY =================================

    ! ndmh: in PAW, we need to add on the compensation density before finding
    ! ndmh: the XC potential
    if (pub_aug.and.(pub_nhat_in_xc.or.pub_usp)) then
       density_fine = density_fine + nhat_den_grad(:,:,:,:,0)
    end if

    ! ndmh: add on core density before calculating xc potential for NLCC
    if (pub_nlcc) then
       factor = 1.0_DP/pub_cell%num_spins
       do is=1,pub_cell%num_spins
          density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
               factor*core_density_fine
       end do
    end if

    ! cks: initialise
    buffer_fine = 0.0_DP

    ! ndmh: calculate XC energy and potential
    ! ddor: include optional TDDFT functional switch if necessary
    if (pub_aug) then
       call xc_energy_potential(density_fine, xc_energy, buffer_fine, &
            pub_fine_grid, pub_aug_den_dim, nhat_den_grad)
    else
       call xc_energy_potential(density_fine, xc_energy, buffer_fine, &
            pub_fine_grid, 0)
    end if

    ! ndmh: now subtract off the compensation density and calculate the
    ! ndmh: integral of the pseudo-density \tilde{n} in the XC potential
    if (pub_aug) then

       ! ndmh: zero the dijhat matrices
       do is=1,pub_cell%num_spins
          call sparse_scale(local_dijhat(is),0.0_DP)
       end do

       ! ndmh: calculate screening of dij now if nhat not included in XC
       if (.not.pub_nhat_in_xc) then
          call augmentation_screen_dij(local_dijhat,lhxc_fine,pub_fine_grid)
       end if

       ! Add the XC potential to the LHXC potential
       lhxc_fine = lhxc_fine + buffer_fine

       ! ndmh: calculate screening of dij now if nhat is included in XC
       if (pub_nhat_in_xc) then
          call augmentation_screen_dij(local_dijhat,lhxc_fine,pub_fine_grid)
       end if

       ! ndmh: subract off the compensation density to get just the band
       ! ndmh: energy part of the XC energy
       density_fine = density_fine - nhat_den_grad(:,:,:,:,0)

       ! ndmh: calculate the energy of the PS valence density and of the
       ! ndmh: augmentation density in the XC potential of the total density
       xc_denpot_energy = 0.0_DP
       xc_nhat_denpot_energy = 0.0_DP
       do is=1,pub_cell%num_spins
          xc_denpot_energy = xc_denpot_energy + &
               integrals_product_on_grid(pub_fine_grid,buffer_fine(:,:,:,is), &
               density_fine)
          xc_nhat_denpot_energy = xc_nhat_denpot_energy  + &
               integrals_product_on_grid(pub_fine_grid,buffer_fine(:,:,:,is), &
               nhat_den_grad(:,:,:,is,0))
       end do

    end if

    ! ========== END XC MATRIX & ENERGY =================================

    ! ========== PAW NONLOCAL MATRIX ====================================
    paw_sphere_energies(:) = 0.0_DP
    if (pub_aug) then
       ! ndmh: Create workspace to hold nonlocpot
       do is=1,pub_cell%num_spins
          local_nonlocpot(is)%structure = 'H'//rep%postfix
          call sparse_create(local_nonlocpot(is))
       end do
       call aug_nonlocal_mat(local_nonlocpot,local_dijhat, &
            projector_denskern,rep%sp_overlap,paw_sphere_energies, &
            show_matrices=.true.)
    end if
    ! ========== END PAW NONLOCAL MATRIX ================================

    ! cks: deallocate fine grid workspace
    deallocate(buffer_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components','buffer_fine',ierr)

    ! ndmh: deallocate augmentation workspace
    if (pub_aug) then
       deallocate(lhxc_fine,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components','lhxc_fine',ierr)
       do is=pub_cell%num_spins,1,-1
          call sparse_destroy(local_dijhat(is))
          call sparse_destroy(local_nonlocpot(is))
          call sparse_destroy(projector_denskern(is))
       end do
       deallocate(local_dijhat,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'local_dijhat',ierr)
       deallocate(local_nonlocpot,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'local_nonlocpot',ierr)
       deallocate(projector_denskern,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'projector_denskern',ierr)
    end if

    ! cks: Free up density workspace
    if (pub_aug) then
       deallocate(nhat_den_grad,stat=ierr)
       call utils_dealloc_check('hamiltonian_energy_components', &
            'nhat_den_grad',ierr)
    end if
    deallocate(density_fine,stat=ierr)
    call utils_dealloc_check('hamiltonian_energy_components','density_fine', &
         ierr)

    kinetic_energy = 0.0_DP
    nonlocpot_energy = 0.0_DP
    do is=1,pub_cell%num_spins

       ! cks: calculate kinetic energy
       kinetic_energy = kinetic_energy + &
            sparse_trace(pur_denskern(is), rep%kinet)

       ! cks: calculate non-local potential energy
       ! ndmh: only if there are nonlocal projectors
       if (pub_any_nl_proj) then
          nonlocpot_energy = nonlocpot_energy + &
               sparse_trace(pur_denskern(is), rep%nonlocpot(1))
       end if

    end do

    !qoh: Calculate and include HF exchange if necessary
    hfx_energy = 0.0_DP
    if (pub_usehfx) then
       do is=1,pub_cell%num_spins
          hfx_energy = hfx_energy - 0.5_DP * sparse_trace(pur_denskern(is),&
               hfexchange(is))
       end do
    end if

    ! jd: Include the smeared ion energy corrections, if necessary
    smeared_ion_self_energy = 0.0_DP
    smeared_ion_nonself_energy = 0.0_DP
    if (pub_is_smeared_ion_rep) then
       smeared_ion_self_energy = smeared_ion_E_self
       smeared_ion_nonself_energy = smeared_ion_E_smeared
    end if

    ! =============== DFT+U contribution =======================
    ! ddor: Calculate DFT+U contribution to the total energy
    !       For this we need to remove spin degeneracy factor
    hubbard_energy = 0.0_DP
    if (pub_hubbard) then
       ! ddor: Temporarily remove spin degeneracy factor
       if (pub_cell%num_spins == 1) call sparse_scale(pur_denskern(1),0.5_DP)
       call hubbard_energy_total(hub, hubbard_energy, pur_denskern, &
            hub_proj_basis, rep%hub_overlap, rep%hub_overlap_t)
       if (pub_cell%num_spins == 1) call sparse_scale(pur_denskern(1),2.0_DP)
    endif
    ! ================= end DFT+U contribution ==================

    ! ndmh: add up total energy
    total_energy = kinetic_energy + locpot_energy + nonlocpot_energy &
            + hartree_energy + xc_energy + hfx_energy + hubbard_energy &
            + ewald_energy + pub_dispersion_energy + smeared_ion_self_energy &
            + smeared_ion_nonself_energy + cavitation_energy
    if (pub_paw) then
       total_sphere_energy = paw_sphere_energies(paw_en_dij0) &
            + paw_sphere_energies(paw_en_ehart) &
            + paw_sphere_energies(paw_en_exc) &
            - paw_sphere_energies(paw_en_etxc) &
            - paw_sphere_energies(paw_en_exc_core)
       total_energy = total_energy + total_sphere_energy
    end if

    if (pub_on_root) then
       write(stdout,'(/a)') '===================================================&
            &============================='
       write(stdout,'(11x,a)') '---------------- ENERGY COMPONENTS (Eh) &
            &----------------'
       write(stdout,'(11x,a,f24.14,a)') &
            '| Kinetic                    :',            kinetic_energy,' |'
       if (.not.pub_paw) then

          if (pub_is_smeared_ion_rep) then
             if (pub_open_localpseudo) then
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudo. (local,OBC,corr''d): ', &
                     locpot_energy,' |'
             else
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudopot. (loc.,corrected):', &
                     locpot_energy,' |'
             end if
          else
             if (pub_open_localpseudo) then
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudopot. (local, OBC):    ', &
                     locpot_energy,' |'
             else
                write(stdout,'(a,f24.14,a)') &
                     '           | Pseudopotential (local)    :', &
                     locpot_energy,' |'
             end if
          end if

          write(stdout,'(11x,a,f24.14,a)') &
               '| Pseudopotential (non-local):',       nonlocpot_energy,' |'
          if (pub_is_smeared_ion_rep) then
             write(stdout,'(a,f24.14,a)') &
                  '           | Hartree (molecular)        :', &
                  hartree_energy,' |'
          else
             if (pub_multigrid_hartree) then
                write(stdout,'(a,f24.14,a)') &
                     '           | Hartree (OBC)              :', &
                     hartree_energy,' |'
             else
                write(stdout,'(a,f24.14,a)') &
                     '           | Hartree                    :', &
                     hartree_energy,' |'
             end if
          end if
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exchange-correlation       :',              xc_energy,' |'
       else
          write(stdout,'(11x,a,f24.14,a)') &
               '| Core Hartree               :',          locpot_energy,' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Core-CompCh Hartree        :',     nhat_locpot_energy,' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Hartree                    :',         hartree_energy,' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Hartree (CompCh)           :',    nhat_hartree_energy,' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exch.-cor. (Total)         :',              xc_energy,' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exch.-cor. (vxc.Valence)   :',       xc_denpot_energy,' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exch.-cor. (vxc.CompCh)    :',  xc_nhat_denpot_energy,' |'
       end if
       if (pub_usehfx) write(stdout,'(11x,a,f24.14,a)') &
            '| Hartree-Fock exchange      :',       hfx_energy,' |'
       if ((pub_hubbard).and.(.not.pub_cdft)) then
            write(stdout,'(11x,a,f24.14,a)') &
            '| Hubbard DFT+U correction   :',   hubbard_energy,' |'
       elseif (pub_cdft) then
            write(stdout,'(11x,a,f24.14,a)') &
            '| cDFT [+(DFT+U)] correction :',   hubbard_energy,' |'
       endif
       if (pub_ii_energy_direct) then
          write(stdout,'(11x,a,f24.14,a)') &
               '| Ion-ion (open BC)          :',  ewald_energy,' |'
       else
          write(stdout,'(11x,a,f24.14,a)') &
               '| Ewald                      :',  ewald_energy,' |'
       end if
       if (pub_paw) then
          write(stdout,'(11x,a,f24.14,a)') &
               '| Dijhat                     :',  paw_sphere_energies(1),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Dij0                       :',  paw_sphere_energies(2),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Dijhartree                 :',  paw_sphere_energies(3),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc [AE core]              :',  paw_sphere_energies(4),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc DC [AE core]           :',  paw_sphere_energies(5),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc [PS core]              :',  paw_sphere_energies(6),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc DC [PS core]           :',  paw_sphere_energies(7),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Dijxc [Total]              :',  paw_sphere_energies(8),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Exc [Core]                 :',  paw_sphere_energies(9),' |'
          write(stdout,'(11x,a,f24.14,a)') &
               '| Total [Sphere]             :',     total_sphere_energy,' |'
       end if
       if (pub_dispersion /= 0) write(stdout,'(11x,a,f24.14,a)') &
            '| Dispersion Correction      :', pub_dispersion_energy,&
            ' |'
       if (pub_is_smeared_ion_rep) then
          write(stdout,'(a,f24.14,a)') &
               '           | Smeared ion non-self corr. :', &
               smeared_ion_nonself_energy,' |'
          write(stdout,'(a,f24.14,a)') &
               '           | Smeared ion self corr.     :', &
               smeared_ion_self_energy,' |'
       end if
       if (pub_is_include_cavitation) then
          write(stdout,'(a,f24.14,a)') &
               '           | Solvent cavitation energy  :', &
               cavitation_energy,' |'
       end if

       write(stdout,'(11x,a,f24.14,a)') '| Total                      :',&
            total_energy,' |'
       write(stdout,'(11x,a)') '----------------------------------------&
            &----------------'

       if (pub_output_detail == VERBOSE) then
          write(stdout,'(11x,a)') '------ LOCAL ENERGY COMPONENTS FROM &
               &MATRIX TRACES ------'
          if (.not.pub_paw) then
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Pseudopotential (local)    :', tr_locpot_energy,' |'
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Hartree                    :',tr_hartree_energy,' |'
          else
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Core-Valence Hartree       :', tr_locpot_energy,' |'
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Core-CompCh Hartree        :', tr_nhat_locpot_energy,' |'
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Hartree (Valence)          :',tr_hartree_energy,' |'
             write(stdout,'(11x,a,f24.14,a)') &
                  '| Hartree (CompCh)           :',tr_nhat_hartree_energy,' |'
          end if
          write(stdout,'(11x,a)') '-------------------------------------&
               &-------------------'
       end if

       if (pub_aug) then
          write(stdout,'(11x,a,f24.14)')  &
               'Integrated compensation chrg.:', sum(integrated_nhat( &
               1:pub_cell%num_spins))
       end if

       write(stdout,'(11x,a,f24.14)')  &
            'Integrated density           :', integrated_density
       if (pub_cell%num_spins > 1) then
          write(stdout,'(11x,a,f24.14)') &
               'Integrated spin density      :', integrated_spin
          write(stdout,'(11x,a,f24.14)') &
               'Integrated |spin density|    :', integrated_mod_spin
       end if

       if (pub_output_detail == VERBOSE) then
          if (pub_cell%num_spins > 1) write(stdout,'(11x,a,f24.14)') &
               'Local density approx. <S^2>  :', &
               spin_squared_expectation
          write(stdout,'(11x,a,f24.14)')  &
               'Integrated density tr(KS)    :', tr_integrated_density
          if (pub_cell%num_spins > 1) write(stdout,'(11x,a,f24.14)') &
               'Integrated spin tr(KS)       :',tr_integrated_spin
       end if

       write(stdout,'(a/)') '===================================================&
            &============================='

    end if

    ! Stop Timer
    call timer_clock('hamiltonian_energy_components',2)

  end subroutine hamiltonian_energy_components


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module hamiltonian
