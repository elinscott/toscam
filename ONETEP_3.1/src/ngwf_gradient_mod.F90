!================================================================!
!                                                                !
!                      NGWF gradients module                     !
!                                                                !
! This module calculates gradients of the NGWFs with respect to  !
! their expansion coefficients in the psinc basis.               !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris and                !
! Arash A. Mostofi in 2000 and 2001.                             !
! Subsequent modifications by Chris-Kriton Skylaris,             !
! Arash A. Mostofi, Peter D. Haynes and Nicholas D. M. Hine.     !
!================================================================!


module ngwf_gradient

  use constants, only: DP

  implicit none

  private

  real(kind=DP), allocatable, dimension(:,:,:) :: precond_func_real
  real(kind=DP), allocatable, dimension(:,:,:) :: precond_func_recip
  real(kind=DP), allocatable, dimension(:,:,:) :: smooth_func_real
  real(kind=DP), allocatable, dimension(:,:,:) :: smooth_func_recip

  public :: ngwf_gradient_exit
  public :: ngwf_gradient_lnv

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_lnv(contra_grad, cov_grad,      &        ! output
       denskern, rep, ngwf_basis, proj_basis, hub_proj_basis, &     ! input
       lhxc_fine, ham, hub, mu, &                                   ! input
       elements, val_dkn, val_rep, val_ham, val_ngwf_basis, cond_shift) ! input

    !==========================================================================!
    ! This subroutine calculates the gradient of the LNV total energy          !
    ! function with respect to the expansion coefficients of the NGWFs.        !
    !==========================================================================!
    ! Arguments:                                                               !
    ! contra_grad (output) : Contravariant NGWF gradient in ppd representation !
    !    for the NGWFs of pub_my_node_id.                                      !
    ! cov_grad (output)    : Covariant NGWF gradient in ppd representation     !
    !    for the NGWFs of pub_my_node_id.                                      !
    ! denskern (input) : Density kernel in SPAM3 format.                       !
    ! rep (input) : NGWF Representation (functions and matrices).              !
    ! ngwf_basis (input) : Function basis describing the NGWFs                 !
    ! proj_basis (input) : Function basis describing the nonlocal projectors   !
    ! hub_proj_basis (input) : Function basis describing Hubbard projectors    !
    ! lhxc_fine (input) : Total local potential in the fine grid of the        !
    !    whole simulation cell.                                                !
    ! ham (input) : Hamiltonian matrix in SPAM3 format.                        !
    ! dijhat (input) : Screened nonlocal energies                              !
    ! hub (input) : HUBBARD_MODEL type defining a Hubbard model                !
    ! hfexchange (input) : Hartree-Fock exchange matrix in SPAM3 format.       !
    ! mu (input) : Lagrange multiplier (original version) or fudge             !
    !    parameter (Millam-Scuseria version) of LNV function.                  !
    ! val_dkn (input) : Valence density kernel (optional)                      !
    ! val_rep (input) : Valence NGWF representation (optional)                 !
    ! val_ham (input) : Valence Hamiltonian matrix (optional)                  !
    ! val_ngwf_basis (input) : Function basis describing the valence NGWF      !
    !    basis (optional)                                                      !
    ! cond_shift (input) : Value by which the projected conduction Hamiltonian !
    !    is currently being shifted (optional)                                 !
    !--------------------------------------------------------------------------!
    ! Key internal variables:                                                  !
    !   batch_size: This is set equal to ngwf_grad_batch_size which comes from !
    !     the rundat module. The value of batch size determines how large is   !
    !     the batch of accumulated fftboxes. Increasing this number decreases  !
    !     the communication per processor but increases the allocated memory   !
    !     per processor.                                                       !
    !==========================================================================!
    ! Originally written by Chris-Kriton Skylaris in January 2001,             !
    ! to use a "pair-box".                                                     !
    ! Rewritten by Arash Mostofi in 2003 so that it uses a "triple-box".       !
    ! Rewritten by Chris-Kriton Skylaris on 23/11/2003 so that it runs on      !
    ! parallel computers.                                                      !
    ! Modified by Chris-Kriton Skylaris on 14/02/2005 to implement mixing      !
    ! of occupancy preconditioning.                                            !
    ! Modified by Peter Haynes to use parallel SPAM 2, July 2006               !
    ! DFT+U added by David O'Regan, April 2009                                 !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, July 2009                !
    ! Modified by Laura Ratcliff for conduction calculations, Oct 2010         !
    ! Modified by Alvaro Ruiz Serrano for kernel DIIS, November 2010.          !
    !==========================================================================!

    use cell_grid, only: pub_dbl_grid, pub_fine_grid
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, stdout, max_spins
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT
    use ion, only: ELEMENT
    use hf_exchange, only: hf_exchange_calculate
    use hubbard_build, only: HUBBARD_MODEL
    use potential, only: potential_input_to_workspace
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: ngwf_grad_batch_size, pub_usehfx, pub_hubbard, &
         pub_any_nl_proj, pub_paw, pub_fine_is_dbl, pub_cond_calculate, &
         pub_kernel_diis, precond_real, precond_recip
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_transpose, &
         sparse_transpose_structure
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(in) :: lhxc_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(out) :: contra_grad(ngwf_basis%n_ppds *pub_cell%n_pts)
    real(kind=DP), intent(out) :: cov_grad(ngwf_basis%n_ppds *pub_cell%n_pts)
    type(NGWF_HAM), intent(inout) :: ham !jd: hf_exchange_calculate writes to it
    real(kind=DP), intent(in) :: mu(max_spins)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    ! lr408: Optional conduction arguments
    type(SPAM3), optional, intent(in) :: val_dkn(pub_cell%num_spins)
    type(SPAM3), optional, intent(in) :: val_ham(pub_cell%num_spins)
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis
    type(NGWF_REP), optional, intent(in) :: val_rep
    real(kind=dp), optional, intent(in) :: cond_shift

    ! Local Variables
    type(SPAM3), allocatable :: coeff_mat(:,:)
    type(SPAM3), allocatable :: proj_coeff_mat(:,:)
    type(SPAM3), allocatable :: hub_proj_coeff_mat(:,:) !ddor
    type(SPAM3), allocatable :: cond_coeff_mat(:,:)
    type(SPAM3) :: ps_overlap
    real(kind=DP), allocatable :: fftbox_batch(:,:,:,:,:,:)
    real(kind=DP), allocatable :: lhxc_dbl(:,:,:,:)
    integer, parameter :: nmat=4
    integer :: is
    integer :: ierr
    integer :: imat
    integer :: batch_size
    integer :: batch_count
    integer :: n_batches
    integer :: local_start, local_end, local_len
    integer :: max_current_size ! maximum batch size over all nodes

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering ngwf_gradient_lnv'
#endif

    ! Start timer
    call timer_clock('ngwf_gradient_lnv',1)

    ! ######## INITIALISATIONS #################################
    contra_grad = 0.0_DP
    cov_grad = 0.0_DP
    batch_size = ngwf_grad_batch_size

    ! ndmh: projector-ngwf overlap matrix
    if (pub_any_nl_proj.or.pub_paw) then
       call sparse_transpose_structure(ps_overlap%structure,rep%sp_overlap)
       call sparse_create(ps_overlap)
       call sparse_transpose(ps_overlap,rep%sp_overlap)
    end if

    ! ndmh: allocate storage for coefficient matrices
    allocate(coeff_mat(pub_cell%num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_lnv','coeff_mat',ierr)
    allocate(cond_coeff_mat(pub_cell%num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_lnv','cond_coeff_mat',ierr)
    allocate(proj_coeff_mat(pub_cell%num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_lnv','proj_coeff_mat',ierr)
    allocate(hub_proj_coeff_mat(pub_cell%num_spins,nmat),stat=ierr)
    call utils_alloc_check('ngwf_gradient_lnv','hub_proj_coeff_mat',ierr)

    do is=1,pub_cell%num_spins
       ! ndmh: create elements of coeff_mat array of NGWF coefficient matrices
       call sparse_create(coeff_mat(is,1),denskern(is))
       call sparse_create(coeff_mat(is,2),denskern(is),rep%overlap)
       call sparse_create(coeff_mat(is,3),denskern(is))
       call sparse_create(coeff_mat(is,4),denskern(is),rep%overlap)
       ! ndmh: create elements of proj_coeff_mat array of projector
       ! ndmh: coefficient matrices
       if (pub_paw) then
          call sparse_create(proj_coeff_mat(is,1),ps_overlap,denskern(is))
          call sparse_create(proj_coeff_mat(is,2),ps_overlap,coeff_mat(is,2))
       end if
       if (pub_any_nl_proj.or.pub_paw) then
          call sparse_create(proj_coeff_mat(is,3),ps_overlap,denskern(is))
          call sparse_create(proj_coeff_mat(is,4),ps_overlap,coeff_mat(is,2))
       end if
       ! lr408: Conduction coefficient matrices
       if (pub_cond_calculate) then
          cond_coeff_mat(is,3)%structure = 'Mc'
          call sparse_create(cond_coeff_mat(is,3))
          call sparse_create(cond_coeff_mat(is,4),cond_coeff_mat(is,3), &
               rep%overlap)
       end if
       ! ddor: Hubbard projector coefficient matrices
       if (pub_hubbard) then
          call sparse_create(hub_proj_coeff_mat(is,3), &
               rep%hub_overlap_t,denskern(is))
          call sparse_create(hub_proj_coeff_mat(is,4), &
               rep%hub_overlap_t,coeff_mat(is,2))
       endif
    end do

    ! ndmh: Initialise the NGWF coefficients:
    ! ndmh: coeff_mat(:,1) (qmat) = (3LHL -2LSLHL -2LHLSL)
    ! ndmh: coeff_mat(:,2) (tc_qmat) = (3LHL -2LSLHL -2LHLSL)S
    ! ndmh: coeff_mat(:,3) (pur_denskern) = \tilde{K}
    ! ndmh: coeff_mat(:,4) (tc_pur_denskern) = \tilde{K}S

    ! ars: If kernel DIIS the NGWF coefficients are:
    ! ars: coeff_mat(:,1) (qmat) = Q = -[Ne/tr(KS)]^2 *KHK - 2(mu-[Ne/tr(KS)]^2 * nu] *K
    ! ars: coeff_mat(:,2) (tc_qmat) = QS
    ! ars: coeff_mat(:,3) (pur_denskern) = P = [2Ne/tr(KS)] *K - [Ne/tr(KS)]^2 *KSK
    ! ars: coeff_mat(:,4) (tc_pur_denskern) = PS
    ! ars: mu = tr[KH]/tr[KS]
    ! ars: nu = tr[HKSK]/tr[KS]

    ! ndmh: Initialise (optionally) the nonlocal projector coefficients:
    ! ndmh: proj_coeff_mat(:,1) rQ = O<proj|ngwf>*qmat (PAW only)
    ! ndmh: proj_coeff_mat(:,2) tc_rQ  = O<proj|ngwf>*tc_qmat (PAW only)
    ! ndmh: proj_coeff_mat(:,3) rk = D<proj|ngwf>*pur_denskern (NCPP,PAW)
    ! ndmh: proj_coeff_mat(:,4) tc_rk = D<proj|ngwf>*tc_pur_denskern (NCPP,PAW)
    ! ndmh: where O_ij is the block-diagonal projector overlap matrix and
    ! ndmh: D_ij is the matrix of PAW nonlocal energies (or KB coeffs for NCPPs)

    ! ndmh: Initialise (optionally) the Hubbard projector coefficients:
    ! ndmh: hub_proj_coeff_mat(:,3) hub_wk = G<proj|ngwf>*pur_denskern
    ! ndmh: hub_proj_coeff_mat(:,4) tc_hub_wk = G<proj|ngwf>*tc_pur_denskern
    ! ndmh: where G is the Hubbard Hamiltonian (block-diagonal)

    ! ndmh: Initialise (optionally) the valence NGWF coefficients
    ! ndmh: in a conduction NGWF optimisation,
    ! ndmh: cond_coeff_mat(:,3) cond_grad, cond_coeff_mat(:,4) tc_cond_grad
    ! ndmh: cond_grad = M<val_ngwf|cond_ngwf>*pur_denskern(cond)
    ! ndmh: tc_cond_grad = M<val_ngwf|cond_ngwf>*tc_pur_denskern(cond)
    ! ndmh: where M = w.Kval.Sval.Kval - Kval.Hval.Kval and w is the
    ! ndmh: conduction state shift parameter in the projection, Kval is
    ! ndmh: the valence density kernel, Sval the valence overlap matrix
    ! ndmh: and Hval the valence Hamiltonian.

    call ngwf_gradient_coeffs(coeff_mat(:,1),coeff_mat(:,2),coeff_mat(:,3), &
         coeff_mat(:,4),proj_coeff_mat(:,1),proj_coeff_mat(:,2), &
         proj_coeff_mat(:,3),proj_coeff_mat(:,4), &
         hub_proj_coeff_mat(:,3),hub_proj_coeff_mat(:,4), &
         cond_coeff_mat(:,3),cond_coeff_mat(:,4), &
         denskern,rep%overlap,ham%ham,rep%inv_overlap,ps_overlap,rep%sp_overlap,&
         hub%projector_ham,rep%hub_overlap_t, &
         mu,rep%n_occ,ham%dijhat, &
         val_dkn, val_rep%overlap, val_rep%sp_overlap, val_ham, &
         rep%cross_overlap, cond_shift)


    ! ### END INITIALISATIONS #################################

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Completed initialisation in &
         &ngwf_gradient_lnv'
#endif

    ! ndmh: allocate preconditioner
    if (precond_recip.or.precond_real) then
       call ngwf_grad_init_precond_recip
    end if
    if (precond_real) then
       call ngwf_grad_init_precond_real
    end if

    ! ndmh: in case the fine grid is denser than the double grid, create a 
    ! ndmh: temporary double grid
    allocate(lhxc_dbl(pub_dbl_grid%ld1,pub_dbl_grid%ld2,&
         pub_dbl_grid%max_group_slabs12,pub_cell%num_spins), stat=ierr)
    call utils_alloc_check('ngwf_gradient_lnv','lhxc_dbl',ierr)

    ! ndmh: allocate storage for fftboxes for this batch
    allocate(fftbox_batch(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3,pub_cell%num_spins,nmat,batch_size), stat=ierr)
    call utils_alloc_check('ngwf_gradient_lnv','fftbox_batch',ierr)

    ! ndmh: filter (or just copy) the lhxc potential to the double grid
    do is=1,pub_cell%num_spins
       call potential_input_to_workspace(lhxc_dbl(:,:,:,is), &
            lhxc_fine(:,:,:,is),pub_dbl_grid,pub_fine_grid)
    end do

    ! cks: number of row-steps per row-block
    n_batches = ngwf_basis%max_on_node / batch_size
    if (mod(ngwf_basis%max_on_node,batch_size) > 0) n_batches = n_batches + 1

    ! ndmh: loop over batches of NGWFs
    local_start = 1
    do batch_count=1,n_batches

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a,2(i4,a))') 'DEBUG: Batch loop ',&
            batch_count, ' of ',n_batches,' in ngwf_gradient_lnv'
#endif

       local_end = min(local_start+batch_size-1,ngwf_basis%node_num)
       local_len = local_end - local_start + 1

       ! cks: maximum size of current batch over all nodes
       max_current_size = local_len
       call comms_reduce('MAX', max_current_size)

       ! cks: zero before accumulation
       fftbox_batch(:,:,:,:,:,1:local_len) = 0.0_DP

       ! ndmh: deposit sums of various functions (the NGWFs, the Hamiltonian
       ! ndmh: acting on the NGWFS, nonlocal projectors, Hubbard projectors,
       ! ndmh: valence NGWFs in a conduction optimisation) to the accumulating
       ! ndmh: FFTboxes, then precondition the gradient, extract it from the
       ! ndmh: FFTBoxes to PPD storage, and then shave them according to the
       ! ndmh: NGWFs radii.
       call ngwf_gradient_batch(contra_grad, cov_grad, &  ! output
            fftbox_batch, lhxc_dbl, rep, ngwf_basis, proj_basis, &
            hub_proj_basis, hub, batch_size, local_start, local_end, &
            coeff_mat, proj_coeff_mat, hub_proj_coeff_mat, cond_coeff_mat, &
            ps_overlap, max_current_size, &
            val_ngwf_basis, val_rep)

       local_start = local_start + batch_size
    end do

    ! pdh: deallocate workspace
    deallocate(fftbox_batch,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_lnv','fftbox_batch',ierr)
    deallocate(lhxc_dbl, stat=ierr)
    call utils_dealloc_check('ngwf_gradient_lnv','lhxc_dbl',ierr)

    ! qoh: calculate NGWF gradient for Hartree-Fock exchange
    if (pub_usehfx) call hf_exchange_calculate(ham%hfexchange, &        ! output
         coeff_mat(:,3), rep%overlap, rep%ngwfs_on_grid, ngwf_basis, &  ! input
         elements, ham%full_vmatrix, .true., coeff_mat(:,4), cov_grad, &
         contra_grad, precond_func_recip)                               ! input

    call ngwf_gradient_exit

    ! ndmh: deallocate coefficient matrices
    do is=pub_cell%num_spins,1,-1
       if (pub_hubbard) then !ddor
          call sparse_destroy(hub_proj_coeff_mat(is,4))
          call sparse_destroy(hub_proj_coeff_mat(is,3))
       endif
       if (pub_any_nl_proj.or.pub_paw) then
          call sparse_destroy(proj_coeff_mat(is,4))
          call sparse_destroy(proj_coeff_mat(is,3))
       end if
       if (pub_paw) then
          call sparse_destroy(proj_coeff_mat(is,2))
          call sparse_destroy(proj_coeff_mat(is,1))
       end if
       if (pub_cond_calculate) then
          call sparse_destroy(cond_coeff_mat(is,4))
          call sparse_destroy(cond_coeff_mat(is,3))
       end if
       do imat=4,1,-1
          call sparse_destroy(coeff_mat(is,imat))
       end do
    end do
    deallocate(hub_proj_coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_lnv','hub_proj_coeff_mat',ierr)
    deallocate(proj_coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_lnv','proj_coeff_mat',ierr)
    deallocate(cond_coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_lnv','cond_coeff_mat',ierr)
    deallocate(coeff_mat,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_lnv','coeff_mat',ierr)
    if (pub_any_nl_proj.or.pub_paw) call sparse_destroy(ps_overlap)

    call timer_clock('ngwf_gradient_lnv', 2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving ngwf_gradient_lnv'
#endif

  end subroutine ngwf_gradient_lnv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_batch(contra_grad, cov_grad, &  ! output
       fftbox_batch, lhxc_dbl, rep, &
       ngwf_basis, proj_basis, hub_proj_basis, hub, &
       batch_size, local_start, local_end, &
       coeff_mat, proj_coeff_mat, hub_proj_coeff_mat, cond_coeff_mat, &
       ps_overlap, max_current_size, &
       val_ngwf_basis, val_rep)

    !==========================================================================!
    ! This subroutine returns the contravariant and covariant NGWF             !
    ! gradients for the NGWFs of the current batch. It does this               !
    ! by applying the Hamiltonian operator to the functions accumulated        !
    ! in the fftbox of each batch, then applying kinetic energy                !
    ! preconditioning if required and then by extracting the                   !
    ! relevant ppds from the fftboxes and shaving their values                 !
    ! so that they are non-zero only within their spheres.                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! contra_grad (output) : Contravariant NGWF gradient in ppd representation !
    !    for the NGWFs of pub_my_node_id.                                      !
    ! cov_grad (output)    : Covariant NGWF gradient in ppd representation     !
    !    for the NGWFs of pub_my_node_id.                                      !
    ! fftbox_batch (input) :                                                   !
    ! lhxc_dbl (input) : Total local potential in the double grid for the      !
    !    whole simulation cell.                                                !
    ! rep (input) : NGWF Representation (functions and matrices).              !
    ! ngwf_basis (input) : Function basis describing the NGWFs                 !
    ! proj_basis (input) : Function basis describing the nonlocal projectors   !
    ! hub_proj_basis (input) : Function basis describing Hubbard projectors    !
    ! batch_size (input) : Size of each batch of NGWFs                         !
    ! local_start (input) : First NGWF in this batch on this node              !
    ! local_end (input) : Last NGWF in this batch on this node                 !
    ! coeff_mat (input) : Coefficients of NGWFs to add to gradient             !
    ! proj_coeff_mat (input) : Coefficients of projectors to add to gradient   !
    ! hub_proj_coeff_mat (input) : Coefficients of Hubbard projectors          !
    ! cond_coeff_mat (input) : Coefficients of valence NGWFs for gradient of   !
    !    conduction NGWFs (optional - conduction NGWF optimisations only).     !
    ! val_ngwf_basis (input) : Valence NGWF basis (optional - cond NGWF only). !
    ! val_rep  (input) : NGWF_REP type for valence NGWFs (optional as above.   !
    !--------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 26/12/2003                           !
    ! Modified for various speed enhancements by Nicholas Hine 2007-2009       !
    ! Modified for DFT+U by David O'Regan, 2009.                               !
    ! Modified to not use workspace_mod, Nicholas Hine, November 2009          !
    ! Tidied up, added more comments, reduced memory usage and made flow       !
    ! more logical, Nicholas Hine, October 2010                                !
    !==========================================================================!

    use basis, only: basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox, &
         basis_extract_function_from_box, basis_clean_function
    use cell_grid, only: cell_grid_extract_box, pub_dbl_grid
    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS, function_basis_sum_fftbox_batch
    use simulation_cell, only : pub_cell, pub_fftbox
    use geometry, only: POINT
    use hubbard_build, only: HUBBARD_MODEL
    use ion, only: element
    use ngwf_representation, only: NGWF_REP
    use paw, only: paw_projectors
    use projectors, only: projectors_gradient_batch
    use pseudopotentials, only: nlps_projectors
    use rundat, only: precond_real, precond_recip, task, &
         pub_any_nl_proj, pub_hubbard, pub_hubbard_restart, &
         pub_hubbard_atomsolve, pub_paw, pub_cond_calculate, &
         pub_dbl_grid_scale, smooth_scheme
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: batch_size
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    integer, intent(in) :: local_start, local_end
    integer, intent(in) :: max_current_size
    real(kind=DP), intent(inout) ::contra_grad(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(inout) :: cov_grad(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(inout) :: fftbox_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3, pub_cell%num_spins, &
         4, batch_size)
    real(kind=DP), intent(in) :: lhxc_dbl(pub_dbl_grid%ld1, &
         pub_dbl_grid%ld2, pub_dbl_grid%max_group_slabs12, pub_cell%num_spins)
    type(SPAM3), intent(in) :: coeff_mat(pub_cell%num_spins,4)
    type(SPAM3), intent(in) :: proj_coeff_mat(pub_cell%num_spins,4)
    type(SPAM3), intent(in) :: hub_proj_coeff_mat(pub_cell%num_spins,4) !ddor
    type(SPAM3), intent(in) :: cond_coeff_mat(pub_cell%num_spins,4) !lr408
    type(SPAM3), intent(in) :: ps_overlap
    type(FUNC_BASIS), optional, intent(in) :: val_ngwf_basis
    type(NGWF_REP), optional, intent(in) :: val_rep

    ! Local Variables
    logical :: i_need_potential
    integer :: prev_start1, prev_start2, prev_start3
    integer :: fa_start1, fa_start2, fa_start3
    integer :: fa_cell_start1, fa_cell_start2, fa_cell_start3
    integer :: fftbox_start1_dbl, fftbox_start2_dbl, fftbox_start3_dbl
    integer :: local_fa
    integer :: batch_count
    integer :: is
    integer :: ierr
    integer :: idx_len
    integer, parameter :: fb_box = 1
    integer, parameter :: tc_fb_box = 2
    integer, parameter :: ham_box = 3
    integer, parameter :: tc_ham_box = 4
    real(kind=DP) :: common_fac
    integer, allocatable, dimension(:) :: overlap_idx
    real(kind=DP), allocatable :: func_on_grid_buffer(:)
    real(kind=DP), allocatable, dimension(:,:,:) :: lhxc_fftbox_dbl
    real(kind=DP), allocatable, dimension(:,:,:) :: kin_buffer
    real(kind=DP), allocatable, dimension(:,:,:) :: tc_kin_buffer
    real(kind=DP), allocatable, dimension(:,:,:) :: buffer_dbl
    real(kind=DP), allocatable, dimension(:,:,:,:) :: dwork_box_dbl
    complex(kind=DP), allocatable, dimension(:,:,:) :: zwork_box

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering ngwf_gradient_batch'
#endif

    call timer_clock('ngwf_gradient_batch',1)

    ! pdh: overlap matrix index
    idx_len = sparse_index_length(rep%ngwf_overlap)
    allocate(overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','overlap_idx',ierr)
    allocate(func_on_grid_buffer(ngwf_basis%func_on_grid_buffer_size),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','func_on_grid_buffer', &
         ierr)

    ! ndmh: calculate \sum_\beta \phi_\beta (r) * C^i_\alpha\beta in
    ! ndmh: FFT boxes for the four C^i matrices.
    ! ndmh: Once this is done:
    ! ndmh: Box 1 (fb_box) contains \sum_b Q^ab phi_b
    ! ndmh: Box 2 (tc_fb_box) contains \sum_b (QS)^a_b phi_b
    ! ndmh: Box 3 (ham_box) contains \sum_b K^ab phi_b
    ! ndmh: Box 4 (tc_ham_box) contains \sum_b (KS)^a_b phi_b
    call sparse_generate_index(overlap_idx,rep%ngwf_overlap)
    common_fac = 4.0_DP * pub_cell%weight / pub_cell%num_spins
    call function_basis_sum_fftbox_batch(fftbox_batch, &
         rep%ngwfs_on_grid, ngwf_basis, ngwf_basis, batch_size, &
         local_start, local_end, overlap_idx, idx_len, coeff_mat, &
         fb_box, tc_ham_box, func_on_grid_buffer, common_fac)

    deallocate(func_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','func_on_grid_buffer',ierr)
    deallocate(overlap_idx,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','overlap_idx',ierr)

    ! ndmh: allocate generic complex FFTbox workspace
    allocate(zwork_box(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','zwork_box',ierr)


    ! KKKKKKKKKKK KINETIC KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

    allocate(kin_buffer(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','kin_buffer',ierr)
    allocate(tc_kin_buffer(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','tc_kin_buffer',ierr)

    do is=1,pub_cell%num_spins
       batch_count = 0
       do local_fa=local_start,local_start+max_current_size+1
          batch_count = batch_count + 1

          if (local_fa <= local_end) then
             ! Calculate T on sum of phi_beta
             call ngwf_gradient_kinetic(kin_buffer, tc_kin_buffer, &
                  fftbox_batch(:,:,:,is,ham_box,batch_count),  &
                  fftbox_batch(:,:,:,is,tc_ham_box,batch_count), &
                  zwork_box)

             ! ndmh: add T \sum_b fb to the fb_box
             fftbox_batch(:,:,:,is,fb_box,batch_count) = &
                  fftbox_batch(:,:,:,is,fb_box,batch_count) + &
                  kin_buffer(:,:,:)

             fftbox_batch(:,:,:,is,tc_fb_box,batch_count) = &
                  fftbox_batch(:,:,:,is,tc_fb_box,batch_count) + &
                  tc_kin_buffer(:,:,:)

          end if

       end do
    end do

    deallocate(tc_kin_buffer,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','tc_kin_buffer',ierr)
    deallocate(kin_buffer,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','kin_buffer',ierr)

    ! KKKKKKKK END KINETIC KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK



    ! LLLLLL LHXC POTENTIAL LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    ! cks: location of phi_alpha wrt FFT box --- in the middle
    ! cks: of the box,  or left where it is in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(fa_start1, fa_start2, fa_start3, &
         pub_fftbox%total_pt1, pub_fftbox%total_pt2, pub_fftbox%total_pt3)

    allocate(lhxc_fftbox_dbl(pub_fftbox%total_ld1_dbl,&
         pub_fftbox%total_ld2_dbl,pub_fftbox%total_pt3_dbl),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','lhxc_fftbox_dbl',ierr)
    allocate(buffer_dbl(pub_fftbox%total_ld1_dbl,pub_fftbox%total_ld2_dbl,&
         pub_dbl_grid%max_group_slabs12),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','buffer_dbl',ierr)
    allocate(dwork_box_dbl(pub_fftbox%total_ld1_dbl,pub_fftbox%total_ld2_dbl,&
         pub_fftbox%total_pt3_dbl,2),stat=ierr)
    call utils_alloc_check('ngwf_gradient_batch','dwork_box_dbl',ierr)

    ! cks: initialisation
    buffer_dbl = 0.0_DP

    ! pdh: loop over spins
    do is=1,pub_cell%num_spins

       prev_start1= -1111111
       prev_start2= -2222222
       prev_start3= -3333333
       batch_count = 0
       do local_fa=local_start,local_start+max_current_size+1
          batch_count = batch_count + 1

          if (local_fa <= local_end) then
             call basis_location_func_wrt_cell(fa_cell_start1, &
                  fa_cell_start2,fa_cell_start3,ngwf_basis%tight_boxes(local_fa))
          else
             ! ndmh: prevent passing of unassigned variables in debug mode
             fa_cell_start1 = -1234
             fa_cell_start2 = -1234
             fa_cell_start3 = -1234
          end if

          ! Extract potential from cell to FFT box if position of FFT box
          ! is different from last time.
          i_need_potential = .false.
          if ((fa_cell_start1 /= prev_start1 .or. &
               fa_cell_start2 /= prev_start2 .or. &
               fa_cell_start3 /= prev_start3) .and. local_fa <= local_end) &
               i_need_potential = .true.

          ! ndmh: cell_to_box routine moved to basis_mod
          if (pub_dbl_grid_scale>1.0_DP) then
             fftbox_start1_dbl = 2*(fa_cell_start1 - fa_start1) + 1
             fftbox_start2_dbl = 2*(fa_cell_start2 - fa_start2) + 1
             fftbox_start3_dbl = 2*(fa_cell_start3 - fa_start3) + 1
          else ! ndmh: dbl_grid is scale 1.0
             fftbox_start1_dbl = (fa_cell_start1 - fa_start1) + 1
             fftbox_start2_dbl = (fa_cell_start2 - fa_start2) + 1
             fftbox_start3_dbl = (fa_cell_start3 - fa_start3) + 1
          end if
          call cell_grid_extract_box(lhxc_fftbox_dbl,&
               buffer_dbl, lhxc_dbl(:,:,:,is), pub_dbl_grid, &
               pub_fftbox%total_pt1_dbl, pub_fftbox%total_pt2_dbl, &
               pub_fftbox%total_pt3_dbl, pub_fftbox%total_ld1_dbl, &
               pub_fftbox%total_ld2_dbl, fftbox_start1_dbl, &
               fftbox_start2_dbl, fftbox_start3_dbl, i_need_potential, .true.)

          prev_start1 = fa_cell_start1
          prev_start2 = fa_cell_start2
          prev_start3 = fa_cell_start3

          if (local_fa <= local_end) then

             ! Calculate V_loc on sum of phi_beta
             if (pub_dbl_grid_scale>1.0_DP) then
                ! ndmh: interpolate function sums to fine grid, multiply
                ! ndmh: by fine grid potential, then filter back to coarse grid
                call ngwf_gradient_local( &
                     fftbox_batch(:,:,:,is,ham_box,batch_count), &   ! input-output
                     fftbox_batch(:,:,:,is,tc_ham_box,batch_count),& ! input-output
                     lhxc_fftbox_dbl, lhxc_fftbox_dbl, &             ! input-output
                     dwork_box_dbl)                                  ! workspace
             else
                ! ndmh: dbl_grid_scale:1, so no need to interpolate
                fftbox_batch(:,:,:,is,ham_box,batch_count) = &
                     fftbox_batch(:,:,:,is,ham_box,batch_count) * &
                     lhxc_fftbox_dbl(:,:,:)
                fftbox_batch(:,:,:,is,tc_ham_box,batch_count) = &
                     fftbox_batch(:,:,:,is,tc_ham_box,batch_count) * &
                     lhxc_fftbox_dbl(:,:,:)
             end if

             ! cks: grad = Vloc*fb + fb
             fftbox_batch(:,:,:,is,ham_box,batch_count) = &
                  fftbox_batch(:,:,:,is,ham_box,batch_count) + &
                  fftbox_batch(:,:,:,is,fb_box,batch_count)

             fftbox_batch(:,:,:,is,tc_ham_box,batch_count) = &
                  fftbox_batch(:,:,:,is,tc_ham_box,batch_count) + &
                  fftbox_batch(:,:,:,is,tc_fb_box,batch_count)

          end if

       end do

    end do

    deallocate(dwork_box_dbl,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','dwork_box_dbl',ierr)
    deallocate(buffer_dbl,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','buffer_dbl',ierr)
    deallocate(lhxc_fftbox_dbl,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','lhxc_fftbox_dbl',ierr)

    ! LLLLLL END LHXC POTENTIAL LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL


    ! NNNNNN NON-LOCAL POTENTIAL NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ! ndmh: NCPP version
    if (pub_any_nl_proj) then
       call projectors_gradient_batch(fftbox_batch, ham_box, tc_ham_box, &
            local_start, local_end, batch_size, &
            fa_start1, fa_start2, fa_start3, &
            ngwf_basis, proj_basis, proj_coeff_mat, ps_overlap, nlps_projectors)
    end if
    ! ndmh: PAW version
    if (pub_paw) then
       call projectors_gradient_batch(fftbox_batch, ham_box, tc_ham_box, &
            local_start, local_end, batch_size, &
            fa_start1, fa_start2, fa_start3, &
            ngwf_basis, proj_basis, proj_coeff_mat, ps_overlap, paw_projectors)
    end if
    ! NNNNNN END NON-LOCAL POTENTIAL NNNNNNNNNNNNNNNNNNNNNNNNNNNN


    ! HHHHHH HUBBARD DFT+U POTENTIAL HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if (pub_hubbard) then
       if (((task == 'HUBBARDSCF') .and. (hub%consistency_iteration > 1)) .or.&
            & pub_hubbard_restart .or. pub_hubbard_atomsolve) then

          idx_len = sparse_index_length(rep%hub_overlap_t)
          allocate(overlap_idx(idx_len),stat=ierr)
          call utils_alloc_check('ngwf_gradient_batch','overlap_idx',ierr)
          allocate(func_on_grid_buffer(&
               hub_proj_basis%func_on_grid_buffer_size),stat=ierr)
          call utils_alloc_check('ngwf_gradient_batch',&
               'func_on_grid_buffer',ierr)

          call sparse_generate_index(overlap_idx,rep%hub_overlap_t)
          common_fac = pub_cell%weight
          call function_basis_sum_fftbox_batch(fftbox_batch, &
               hub%consistency_projs, hub_proj_basis, ngwf_basis, &
               batch_size, local_start, &
               local_end, overlap_idx, idx_len, hub_proj_coeff_mat, &
               ham_box, tc_ham_box, func_on_grid_buffer, common_fac)

          deallocate(func_on_grid_buffer,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_batch',&
               'func_on_grid_buffer',ierr)
          deallocate(overlap_idx,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_batch','overlap_idx',ierr)
       else
          call projectors_gradient_batch(fftbox_batch, ham_box, tc_ham_box, &
               local_start, local_end, batch_size, &
               fa_start1, fa_start2, fa_start3, &
               ngwf_basis, hub_proj_basis, hub_proj_coeff_mat, &
               rep%hub_overlap_t, hub%projectors)
       endif
    endif
    ! HHHHHH END HUBBARD DFT+U POTENTIAL HHHHHHHHHHHHHHHHHHHHHHHHH


    ! CCCCCCC CONDUCTION TERM CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    if (pub_cond_calculate) then

       ! lr408: Indexing for cross overlap matrix
       idx_len = sparse_index_length(rep%cross_overlap)
       allocate(overlap_idx(idx_len),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch','overlap_idx',ierr)
       allocate(func_on_grid_buffer( &
            val_ngwf_basis%func_on_grid_buffer_size),stat=ierr)
       call utils_alloc_check('ngwf_gradient_batch', &
            'func_on_grid_buffer',ierr)

       call sparse_generate_index(overlap_idx,rep%cross_overlap)
       common_fac = 4.0_DP * pub_cell%weight / pub_cell%num_spins
       call function_basis_sum_fftbox_batch(fftbox_batch, &
            val_rep%ngwfs_on_grid, val_ngwf_basis, ngwf_basis, batch_size, &
            local_start, local_end, overlap_idx, idx_len, cond_coeff_mat, &
            ham_box, tc_ham_box, func_on_grid_buffer, common_fac)

       deallocate(func_on_grid_buffer,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_batch','func_on_grid_buffer', &
            ierr)
       deallocate(overlap_idx,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_batch','overlap_idx',ierr)

    end if
    ! CCCCCCC END CONDUCTION TERM CCCCCCCCCCCCCCCCCCCCCCCCCCCCC


    ! RRRRRRR RECIPROCAL SPACE PRECONDITIONING RRRRRRRRRRRRRRRR
    if (precond_recip) then
       do is=1,pub_cell%num_spins
          batch_count = 1
          do local_fa=local_start,local_end,2

             ! Copy the covariant gradient into complex workspace array
             if (local_fa < local_end) then
                zwork_box = &
                     cmplx(fftbox_batch(:,:,:,is,tc_ham_box,batch_count), &
                     fftbox_batch(:,:,:,is,tc_ham_box,batch_count+1),kind=DP)
             else
                zwork_box = &
                     cmplx(fftbox_batch(:,:,:,is,tc_ham_box,batch_count), &
                     0.0_DP,kind=DP)
             end if

             ! Forward FFT the covariant gradient to reciprocal space
             call fourier_apply_box('Coarse','Forward',zwork_box)

             ! Apply kinetic energy preconditioning to covariant gradient
             call ngwf_gradient_precond_recip(zwork_box)

             ! Backward FFT the covariant gradient to real space
             call fourier_apply_box('Coarse','Backward',zwork_box)

             ! Copy the preconditioned covariant gradient out of the complex
             ! workspace array
             fftbox_batch(:,:,:,is,tc_ham_box,batch_count) = &
                  real(zwork_box,kind=DP)
             if (local_fa < local_end) &
                  fftbox_batch(:,:,:,is,tc_ham_box,batch_count+1) = &
                  aimag(zwork_box)

             batch_count = batch_count + 2
          end do
       end do
    end if
    ! RRRRRRR RECIPROCAL SPACE PRECONDITIONING RRRRRRRRRRRRRRRR

    ! ndmh: APPLY REAL SPACE PRECONDITIONING, APPLY SPIN-AVERAGING,
    ! ndmh: THEN EXTRACT AND SHAVE GRADIENT
    batch_count = 0
    do local_fa=local_start,local_end
       batch_count = batch_count + 1

       if (precond_real) then
          do is=1,pub_cell%num_spins
             ! cks: precondition the covariant gradient in real space
             call ngwf_grad_precond_rspace(  &
                  fftbox_batch(:,:,:,is,tc_ham_box,batch_count), & ! in/out
                  fa_start1, fa_start2, fa_start3, &
                  ngwf_basis%tight_boxes(local_fa)%tight_pts1, &
                  ngwf_basis%tight_boxes(local_fa)%tight_pts2, &
                  ngwf_basis%tight_boxes(local_fa)%tight_pts3)
          end do
       end if

       ! ndmh: Average gradients over spins for spin-polarised calculations.
       if (pub_cell%num_spins == 2) then
          fftbox_batch(:,:,:,1,ham_box,batch_count) = &
               fftbox_batch(:,:,:,1,ham_box,batch_count) + &
               fftbox_batch(:,:,:,2,ham_box,batch_count)
          fftbox_batch(:,:,:,1,tc_ham_box,batch_count) = &
               fftbox_batch(:,:,:,1,tc_ham_box,batch_count) + &
               fftbox_batch(:,:,:,2,tc_ham_box,batch_count)
       end if

       if (smooth_scheme .ne. 'NONE') then
          ! smmd: smooth the covariant gradient in real space
          call ngwf_grad_smooth_rspace(  &
               fftbox_batch(:,:,:,1,tc_ham_box,batch_count), & ! in/out
               fa_start1, fa_start2, fa_start3, &
               ngwf_basis%tight_boxes(local_fa), &
               ngwf_basis%spheres(local_fa)%centre, &
               ngwf_basis%spheres(local_fa)%radius)
       endif

       ! cks: shaving - stage 1 / extract ppds from the FFT box
       ! extract the contravariant gradient from ham_box to PPDs
       call basis_extract_function_from_box(contra_grad, &
            pub_fftbox%total_ld1,pub_fftbox%total_ld2,pub_fftbox%total_pt3, &
            fftbox_batch(:,:,:,1,ham_box,batch_count), &
            ngwf_basis%spheres(local_fa),&
            ngwf_basis%tight_boxes(local_fa), fa_start1, fa_start2, fa_start3, &
            ngwf_basis%spheres(local_fa)%offset)

       ! extract the covariant gradient from ham_box to PPDs
       call basis_extract_function_from_box(cov_grad, &
            pub_fftbox%total_ld1,pub_fftbox%total_ld2,pub_fftbox%total_pt3, &
            fftbox_batch(:,:,:,1,tc_ham_box,batch_count), &
            ngwf_basis%spheres(local_fa),&
            ngwf_basis%tight_boxes(local_fa), fa_start1, fa_start2, fa_start3, &
            ngwf_basis%spheres(local_fa)%offset)

       ! cks: shaving - stage 2 / zero points outside NGWF sphere in PPD rep.
       ! smmd: if smoothing scheme applied previously, no needs to clean the
       ! smmd: covariant gradient (nor the contravariant gradient which is
       ! smmd: always used within contra.DOT.cov expressions)
       if (pub_cell%n_pts > 1 .and. smooth_scheme .eq. 'NONE') then
          call basis_clean_function(contra_grad, &
               ngwf_basis%spheres(local_fa), ngwf_basis%n_ppds)
          call basis_clean_function(cov_grad, &
               ngwf_basis%spheres(local_fa), ngwf_basis%n_ppds)
       endif

    end do
    ! ndmh: END APPLY REAL SPACE PRECONDITIONING, APPLY SPIN-AVERAGING,
    ! ndmh: THEN EXTRACT AND SHAVE GRADIENT

    ! ndmh: deallocate generic workspace
    deallocate(zwork_box,stat=ierr)
    call utils_dealloc_check('ngwf_gradient_batch','zwork_box',ierr)

    call timer_clock('ngwf_gradient_batch',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving ngwf_gradient_batch'
#endif

  end subroutine ngwf_gradient_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_kinetic(out1,out2,in1,in2,zwork_box)

    !=========================================================!
    ! This subroutine applies the kinetic energy operator to  !
    ! an NGWF in an FFTbox.                                   !
    !---------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris as          !
    ! ngwf_gradient_kinloc_pairbox on 15/6/2001.              !
    ! Rewritten by Arash A. Mostofi and modified to use       !
    ! complex-to-complex FFTs on April 2003                   !
    !=========================================================!

    use fourier,         only: fourier_apply_box_pair
    use kinetic,         only: kinetic_apply_on_box
    use simulation_cell, only: pub_fftbox

    implicit none

    real(kind=DP), intent(out)   :: out1(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out)   :: out2(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(inout) :: in1(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(inout) :: in2(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    complex(kind=DP), intent(out) :: zwork_box(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)

    ! Forward Fourier transform to reciprocal space
    call fourier_apply_box_pair('Coarse','Forward',in1,in2,zwork_box)

    ! Apply kinetic operator in reciprocal space
    call kinetic_apply_on_box(zwork_box)

    ! Backward Fourier transform to real space
    call fourier_apply_box_pair('Coarse','Backward',out1,out2,zwork_box)

    return

  end subroutine ngwf_gradient_kinetic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine ngwf_gradient_local(data_box1, data_box2, pot_dbl1, pot_dbl2, &
       dwork_box_dbl)

    !=========================================================!
    ! This subroutine applies the local potential to a pair   !
    ! of NGWFs in FFTboxes.                                   !
    !---------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris as          !
    ! ngwf_gradient_kinloc_pairbox on 15/6/2001.              !
    ! Rewritten by Arash A. Mostofi and modified to use       !
    ! complex-to-complex FFTs in January 2003.                !
    ! Modified by Chris-Kriton Skylaris on 28/12/2003 to run  !
    ! with the parallel version of ONETEP.                    !
    !=========================================================!

    use fourier,         only: fourier_interpolate,fourier_filter
    use simulation_cell, only: pub_fftbox

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: data_box1(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(inout) :: data_box2(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(in)    :: pot_dbl1(pub_fftbox%total_ld1_dbl, &
         pub_fftbox%total_ld2_dbl,pub_fftbox%total_pt3_dbl)
    real(kind=DP), intent(in)    :: pot_dbl2(pub_fftbox%total_ld1_dbl, &
         pub_fftbox%total_ld2_dbl,pub_fftbox%total_pt3_dbl)
    real(kind=DP), intent(out)   :: dwork_box_dbl(pub_fftbox%total_ld1_dbl, &
         pub_fftbox%total_ld2_dbl,pub_fftbox%total_pt3_dbl,2)


    call fourier_interpolate(data_box1,data_box2,dwork_box_dbl(:,:,:,1), &
         dwork_box_dbl(:,:,:,2))

    dwork_box_dbl(:,:,:,1) = dwork_box_dbl(:,:,:,1) * pot_dbl1
    dwork_box_dbl(:,:,:,2) = dwork_box_dbl(:,:,:,2) * pot_dbl2

    call fourier_filter(dwork_box_dbl(:,:,:,1),dwork_box_dbl(:,:,:,2), &
         data_box1, data_box2)

    return

  end subroutine ngwf_gradient_local


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_smooth_rspace( &
       fftbox_grad, fa_start1, fa_start2, fa_start3, &
       fa_tight_box, sphere_centre, sphere_radius)

    !======================================================================!
    ! Smooth a function inside an FFTbox by convolution in real space      !
    ! with a localized low-pass filter                                     !
    !----------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on January 2011                        !
    !     largely based on ngwf_grad_precond_rspace originally written by  !
    !     Chris-Kriton Skylaris on 29/4/2003.                              !
    !======================================================================!

    use basis, only: basis_copy_tightbox_to_fftbox, function_tight_box, &
           basis_func_centre_wrt_fftbox, basis_location_func_wrt_cell
    use comms, only: pub_on_root
    use constants, only: DP, stdout, two_pi
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, operator(*), operator(.DOT.), operator(+), &
           local_displacement, geometry_distance
    use rundat, only: smooth_scheme, r_smooth
    use simulation_cell, only : pub_fftbox, pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: fa_start1, fa_start2, fa_start3
    type(function_tight_box), intent(in) :: fa_tight_box
    type(POINT), intent(in)   :: sphere_centre
    real(kind=DP), intent(in) :: sphere_radius
    real(kind=DP), intent(inout)  :: fftbox_grad(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)

    ! Local Variables
    integer :: fa_tight_pts1, fa_tight_pts2, fa_tight_pts3
    integer :: row1,row2,row3,step1,step2,step3
    integer :: npt1,npt2,npt3
    integer :: row1_start, row2_start, row3_start
    integer :: period_pt1, period_pt2, period_pt3
    integer :: row1_end, row2_end, row3_end
    integer :: tc1, tc2, tc3
    real(kind=DP), dimension(:,:,:), allocatable  :: smooth_tightbox_grad
    integer :: ierr

    call timer_clock("ngwf_grad_smooth_rspace",1)

    !==========================================================================!
    !===== Preliminaries ======================================================!

    fa_tight_pts1 = fa_tight_box%tight_pts1
    fa_tight_pts2 = fa_tight_box%tight_pts2
    fa_tight_pts3 = fa_tight_box%tight_pts3

    ! smmd: determine number of grid points in smoothing radius
    npt1 = nint( r_smooth/pub_cell%d1 )
    npt2 = nint( r_smooth/pub_cell%d2 )
    npt3 = nint( r_smooth/pub_cell%d3 )

    ! smmd : debug
    if (pub_on_root) then
       write(stdout,*) "Smoothing npts : ", npt1, npt2, npt3
    endif

    ! smmd: allocate tightbox for smooth gradient
    allocate(smooth_tightbox_grad(fa_tight_pts1, fa_tight_pts2, fa_tight_pts3 ),&
         stat=ierr)
    call utils_alloc_check('ngwf_grad_smooth_rspace','smooth_tightbox_grad',ierr)
    smooth_tightbox_grad(1: fa_tight_pts1, 1: fa_tight_pts2, &
         1: fa_tight_pts3) = 0.0_DP

    ! smmd: Initialise smoothing function in real space
    if ( .not. allocated(smooth_func_real) ) then
       call internal_init_smooth()
    endif

    ! smmd: Shave gradient in FFTbox
    call internal_shave_in_fftbox(sphere_radius-r_smooth)

    ! smmd : debug
    if (pub_on_root) then
       write(stdout,*) "Smoothing shave : done ! "
    endif

    !==========================================================================!
    !===== Real space convolution ( gradient & smoothing function ) ===========!

    ! smmd: case without periodic boundary conditions on the fftbox
    if (.not.(pub_fftbox%coin3 .or. &
        pub_fftbox%coin2 .or. pub_fftbox%coin1)) then

       ! cks: loop over all points selected for convolution
       do step3 = -npt3, npt3
          do step2 = -npt2, npt2
             do step1 = -npt1, npt1

                ! cks: convolve only in a spherical region around central point
                ! cks: and only for the points of the (tight) fft_box of fa
                ! cks: which is much smaller than the fft box
                if ( smooth_func_real(step1,step2,step3) .ne. 0.0_dp ) then

                   row3_start =max(1 -step3, fa_start3)
                   row2_start =max(1 -step2, fa_start2)
                   row1_start =max(1 -step1, fa_start1)

                   row3_end = min(pub_fftbox%total_pt3 - step3, &
                        fa_start3 + fa_tight_pts3 -1)
                   row2_end = min(pub_fftbox%total_pt2 - step2, &
                        fa_start2 + fa_tight_pts2 -1)
                   row1_end = min(pub_fftbox%total_pt1 - step1, &
                        fa_start1 + fa_tight_pts1 -1)

                   tc3 = row3_start - fa_start3                      ! smmd
                   do row3= row3_start, row3_end
                      tc3 =tc3+1

                      tc2 = row2_start - fa_start2                   ! smmd
                      do row2= row2_start, row2_end
                         tc2 =tc2 +1

                         tc1 = row1_start - fa_start1                ! smmd
                         do row1= row1_start, row1_end
                            tc1 =tc1 +1

                            smooth_tightbox_grad(tc1, tc2, tc3) = &
                                 smooth_tightbox_grad(tc1, tc2, tc3) &
                                 + fftbox_grad(row1+step1,row2+step2,row3+step3) &
                                 * smooth_func_real(step1,step2,step3)

                         enddo
                      enddo
                   enddo

                endif

             enddo
          enddo
       enddo

    ! smmd: case with periodic boundary conditions on the fftbox
    else

       ! cks: loop over all points selected for convolution
       do step3 = -npt3, npt3
          do step2 = -npt2, npt2
             do step1 = -npt1, npt1

                ! cks: convolve only in a spherical region around central point
                ! cks: and only for the points of the tight fft_box of fa
                ! cks: which is much smaller than the fft box
                if ( smooth_func_real(step1,step2,step3) .ne. 0.0_dp ) then

                   if (pub_fftbox%coin3) then
                      row3_start = fa_start3
                      row3_end = fa_start3 + fa_tight_pts3 -1
                   else
                      row3_start =max(1 -step3, fa_start3)
                      row3_end = min(pub_fftbox%total_pt3 - step3, &
                           fa_start3 + fa_tight_pts3 -1)
                   endif

                   if (pub_fftbox%coin2) then
                      row2_start = fa_start2
                      row2_end = fa_start2 + fa_tight_pts2 -1
                   else
                      row2_start =max(1 -step2, fa_start2)
                      row2_end = min(pub_fftbox%total_pt2 - step2, &
                           fa_start2 + fa_tight_pts2 -1)
                   endif

                   if (pub_fftbox%coin1) then
                      row1_start = fa_start1
                      row1_end = fa_start1 + fa_tight_pts1 -1
                   else
                      row1_start =max(1 -step1, fa_start1)
                      row1_end = min(pub_fftbox%total_pt1 - step1, &
                           fa_start1 + fa_tight_pts1 -1)
                   endif

                   tc3 = 0
                   do row3 = row3_start, row3_end
                      tc3 = tc3 + 1
                      if (pub_fftbox%coin3) then
                         period_pt3 = modulo(row3+step3-1,pub_fftbox%total_pt3)+1
                      else
                         period_pt3 = row3 + step3
                      endif

                      tc2 = 0
                      do row2 = row2_start, row2_end
                         tc2 = tc2 + 1
                         if (pub_fftbox%coin2) then
                            period_pt2 = modulo(row2+step2-1,pub_fftbox%total_pt2)+1
                         else
                            period_pt2 = row2 + step2
                         endif

                         tc1 = 0
                         do row1 = row1_start, row1_end
                            tc1 = tc1 + 1
                            if (pub_fftbox%coin1) then
                               period_pt1 = modulo(row1+step1-1,pub_fftbox%total_pt1)+1
                            else
                               period_pt1 = row1 + step1
                            endif

                            smooth_tightbox_grad(tc1, tc2, tc3) = &
                                 smooth_tightbox_grad(tc1, tc2, tc3) &
                                 + fftbox_grad(period_pt1,period_pt2,period_pt3) &
                                 * smooth_func_real(step1,step2,step3)


                         enddo
                      enddo
                   enddo

                endif

             enddo
          enddo
       enddo

    endif

    ! smmd : debug
    if (pub_on_root) then
       write(stdout,*) "Smoothing convolution : done ! "
    endif

    ! cks: copy the tightbox with the smooth gradient to
    ! cks: the appropriate position in the fftbox that contained the
    ! cks: original gradient.
    call basis_copy_tightbox_to_fftbox(fftbox_grad, &                    !output
         fa_start1, fa_start2, fa_start3, &                               !input
         smooth_tightbox_grad, fa_tight_pts1, fa_tight_pts2, fa_tight_pts3) !input


    ! cks: free up memory
    deallocate(smooth_tightbox_grad,stat=ierr)
    call utils_dealloc_check('ngwf_grad_smooth_rspace','prec_tightbox_grad',&
         ierr)

    ! smmd : debug
    if (pub_on_root) then
       write(stdout,*) "Smoothing complete : done ! "
    endif

    call timer_clock("ngwf_grad_smooth_rspace", 2)


    contains

    !=========================================================================!
    !=========================================================================!

    subroutine internal_shave_in_fftbox(rcut)

      implicit none

      ! Argument
      real(kind=DP), intent(in) :: rcut

      ! Local Variables
      type(POINT) :: fa_centre
      type(POINT) :: current_point
      type(POINT) :: periodic_point
      real(kind=DP) :: distance
      real(kind=DP) :: pt1, pt2, pt3
      integer :: fftbox_start1, fftbox_start2, fftbox_start3
      integer :: period1, period2, period3
      integer :: ip1, ip2, ip3
      logical :: in_sphere

      ! Centre of the current function
      call basis_location_func_wrt_cell(fftbox_start1, fftbox_start2, &
           fftbox_start3, fa_tight_box)
      fa_centre =  basis_func_centre_wrt_fftbox(sphere_centre, &
           fa_start1,fa_start2,fa_start3,fftbox_start1,fftbox_start2, &
           fftbox_start3)

      ! Check for periodic boundary conditions
      period1 = 0
      period2 = 0
      period3 = 0
      if (pub_fftbox%coin1) period1 = 1
      if (pub_fftbox%coin2) period2 = 1
      if (pub_fftbox%coin3) period3 = 1

      ! Loop over the points in tight box
      do tc3 = 1, fa_tight_pts3
         pt3 = real(fa_start3 + tc3 - 1, DP)*pub_cell%d3

         do tc2 = 1, fa_tight_pts2
            pt2 = real(fa_start2 + tc2 - 1, DP)*pub_cell%d2

            do tc1 = 1, fa_tight_pts1
               pt1 = real(fa_start1 + tc1 - 1, DP)*pub_cell%d1

               in_sphere = .false.
               current_point = local_displacement(pub_cell%a1_unit,&
                       pub_cell%a2_unit,pub_cell%a3_unit,pt1,pt2,pt3)

               ! Account for the periodic images of current_point in neighboring cells
               do ip3 = -period3, period3
                  do ip2 = -period2, period2
                     do ip1 = -period1, period1

                        periodic_point = current_point &
                              + real(ip1,DP)*pub_cell%a1 &
                              + real(ip2,DP)*pub_cell%a2 &
                              + real(ip3,DP)*pub_cell%a3
                        distance = geometry_distance(periodic_point,fa_centre)
                        if (distance .lt. rcut) in_sphere = .true.

                     enddo
                  enddo
               enddo

               if (.not. in_sphere) &
                  fftbox_grad(fa_start1+tc1-1,fa_start2+tc2-1,fa_start3+tc3-1) = 0.0d0

            enddo
         enddo
      enddo

    end subroutine internal_shave_in_fftbox

    !=========================================================================!
    !=========================================================================!

    subroutine internal_init_smooth()

      implicit none

      ! Local Variables
      type(POINT) :: distance_vec
      complex(kind=DP), allocatable :: smooth_complex(:,:,:)
      real(kind=DP) :: rr
      integer :: s1, s2, s3

      allocate(smooth_func_real(-npt1:npt1,-npt2:npt2,-npt3:npt3),stat=ierr)
      call utils_alloc_check('ngwf_grad_smooth_rspace','smooth_func_real',&
           ierr)

      if (smooth_scheme == 'BG' .or. &
          smooth_scheme == 'MAURI' .or. &
          smooth_scheme == 'TETER' ) then

         allocate(smooth_complex(pub_fftbox%total_ld1, pub_fftbox%total_ld2, &
              pub_fftbox%total_pt3),stat=ierr)
         call utils_alloc_check('ngwf_grad_smooth_rspace','smooth_complex',&
              ierr)

         if (.not. allocated(smooth_func_recip)) &
              call ngwf_grad_init_smooth_recip

         smooth_complex = cmplx(smooth_func_recip,0.0_DP,kind=DP)

         ! aam: use new fourier routines
         call fourier_apply_box('Coarse','Backward',smooth_complex)

      endif

      do step3 = -npt3, npt3
         s3 = abs(step3)
         do step2 = -npt2, npt2
            s2 = abs(step2)
            do step1 = -npt1, npt1
               s1 = abs(step1)

               distance_vec = real(step1,kind=DP)*pub_cell%d1*pub_cell%a1_unit&
                    + real(step2,kind=DP)*pub_cell%d2*pub_cell%a2_unit &
                    + real(step3,kind=DP)*pub_cell%d3*pub_cell%a3_unit

               ! spherical convolution cut-off radius
               rr = sqrt( distance_vec.DOT.distance_vec  )

               if ( rr.le.r_smooth ) then
                  if (smooth_scheme .eq. 'GAUSSIAN') then
                     smooth_func_real(step1,step2,step3) &
                           = (1/(sqrt(two_pi)*(r_smooth/6)))*exp(-rr**2/(2*(r_smooth/6)**2))
                          != exp(-0.025)*bessel_in(rr,0.025)
                  else
                     smooth_func_real(step1,step2,step3) &
                          = real(smooth_complex(s1+1,s2+1,s3+1),kind=DP)
                  endif
               else
                  smooth_func_real(step1,step2,step3) = 0.0_DP
               endif

            enddo
         enddo
      enddo

      if (smooth_scheme .ne. 'GAUSSIAN') then

         ! free unused memory
         deallocate(smooth_complex,stat=ierr)
         call utils_dealloc_check('ngwf_grad_smooth_rspace','smooth_complex',&
              ierr)
      endif

    end subroutine internal_init_smooth

    !=========================================================================!
    !=========================================================================!

  end subroutine ngwf_grad_smooth_rspace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_init_precond_real

    !======================================================================!
    ! Initialises the kinetic energy preconditioner in real space.         !
    !----------------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 29/4/2003.            !
    ! Modified by Arash A. Mostofi on 26/6/2003, 10/7/2003 and 8/8/2003.   !
    ! Modified for speed by Chris-Kriton Skylaris on 1/3/2004.             !
    ! Modified for accuracy by Simon M.-M. Dubois on 15/04/2011            !
    ! Split to own subroutine by Nicholas Hine on 05/08/2011.              !
    !======================================================================!

    use constants, only: DP
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, operator(*), operator(.DOT.), operator(+)
    use rundat, only: r_precond
    use simulation_cell, only : pub_fftbox, pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Local Variables
    integer :: step1,step2,step3
    integer :: s1,s2,s3,npt1,npt2,npt3
    integer :: ierr
    real(kind=DP) :: rr

    complex(kind=DP), dimension(:,:,:), allocatable :: prec_complex
    type(POINT) :: distance_vec

    ! determine number of grid points in preconditioning radius
    npt1 = nint( r_precond/pub_cell%d1 )
    npt2 = nint( r_precond/pub_cell%d2 )
    npt3 = nint( r_precond/pub_cell%d3 )

    ! ===============================================================
    ! ========= INIT REAL SPACE PRECONDITIONER ======================
    if ( .not. allocated(precond_func_real) ) then

       allocate(precond_func_real(-npt1:npt1,-npt2:npt2,-npt3:npt3),stat=ierr)
       call utils_alloc_check('ngwf_grad_init_precond_real','precond_func_real',&
            ierr)
       allocate(prec_complex(pub_fftbox%total_ld1, pub_fftbox%total_ld2, &
            pub_fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('ngwf_grad_init_precond_real','prec_complex',&
            ierr)

       prec_complex = cmplx(precond_func_recip,0.0_DP,kind=DP)

       ! aam: use new fourier routines
       call fourier_apply_box('Coarse','Backward',prec_complex)

       do step3 = -npt3, npt3
          s3 = abs(step3)
          do step2 = -npt2, npt2
             s2 = abs(step2)
             do step1 = -npt1, npt1
                s1 = abs(step1)

                distance_vec = real(step1,kind=DP)*pub_cell%d1*pub_cell%a1_unit&
                     + real(step2,kind=DP)*pub_cell%d2*pub_cell%a2_unit &
                     + real(step3,kind=DP)*pub_cell%d3*pub_cell%a3_unit

                ! spherical convolution cut-off radius
                rr = sqrt( distance_vec.DOT.distance_vec  )

                if ( rr.le.r_precond ) then
                   precond_func_real(step1,step2,step3) &
                        = real(prec_complex(s1+1,s2+1,s3+1),kind=DP)
                else
                   precond_func_real(step1,step2,step3) = 0.0_DP
                endif

             enddo
          enddo
       enddo


       ! cks: free unused memory
       deallocate(prec_complex,stat=ierr)
       call utils_dealloc_check('ngwf_grad_init_precond_real','prec_complex',&
            ierr)

       ! smmd: no need to scale the preconditioning function
       ! smmd: the following code lines are therefore commented
       !precond_func_real(-npt1:npt1,-npt2:npt2,-npt3:npt3) = &
       !     1.0_DP / precond_func_real(0,0,0) &
       !     * precond_func_real(-npt1:npt1,-npt2:npt2,-npt3:npt3)

    endif
    ! ===== END INIT REAL SPACE PRECONDITIONER ======================
    ! ===============================================================

  end subroutine ngwf_grad_init_precond_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_precond_rspace( &
       fftbox_grad, fa_start1, fa_start2, fa_start3, &
       fa_tight_pts1, fa_tight_pts2, fa_tight_pts3)

    !======================================================================!
    ! Performs kinetic energy preconditioning in real space on a function  !
    ! inside an FFTbox.                                                    !
    !----------------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 29/4/2003.            !
    ! Modified by Arash A. Mostofi on 26/6/2003, 10/7/2003 and 8/8/2003.   !
    ! Modified for speed by Chris-Kriton Skylaris on 1/3/2004.             !
    ! Modified for accuracy by Simon M.-M. Dubois on 15/04/2011            !
    !======================================================================!

    use basis, only: basis_copy_tightbox_to_fftbox
    use constants, only: DP
    use rundat, only: r_precond
    use simulation_cell, only : pub_fftbox, pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: fa_start1, fa_start2, fa_start3
    integer, intent(in) :: fa_tight_pts1
    integer, intent(in) :: fa_tight_pts2
    integer, intent(in) :: fa_tight_pts3
    real(kind=DP), intent(inout)  :: fftbox_grad(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)

    ! Local Variables
    integer :: row1,row2,row3
    integer :: period_pt1,period_pt2,period_pt3
    integer :: step1,step2,step3
    integer :: npt1,npt2,npt3
    integer :: row1_start
    integer :: row2_start
    integer :: row3_start
    integer :: row1_end
    integer :: row2_end
    integer :: row3_end
    integer :: tc1
    integer :: tc2
    integer :: tc3
    integer :: ierr
    real(kind=DP), dimension(:,:,:), allocatable  :: prec_tightbox_grad

    call timer_clock("ngwf_grad_precond_rspace",1)

    allocate(prec_tightbox_grad(fa_tight_pts1, fa_tight_pts2, fa_tight_pts3 ),&
         stat=ierr)
    call utils_alloc_check('ngwf_grad_precond_rspace','prec_tightbox_grad',&
         ierr)

    ! determine number of grid points in preconditioning radius
    npt1 = nint( r_precond/pub_cell%d1 )
    npt2 = nint( r_precond/pub_cell%d2 )
    npt3 = nint( r_precond/pub_cell%d3 )

    prec_tightbox_grad(1: fa_tight_pts1, 1: fa_tight_pts2, &
         1: fa_tight_pts3) = 0.0_DP

    ! smmd: case without periodic boundary conditions on the fftbox
    if (.not.(pub_fftbox%coin3 .or. &
        pub_fftbox%coin2 .or. pub_fftbox%coin1)) then

       ! cks: loop over all points selected for convolution
       do step3 = -npt3, npt3
          do step2 = -npt2, npt2
             do step1 = -npt1, npt1

                ! cks: convolve only in a spherical region around central point
                ! cks: and only for the points of the (tight) fft_box of fa
                ! cks: which is much smaller than the fft box
                ! cks: ( = the whole precond_grad array).
                if ( precond_func_real(step1,step2,step3) .ne. 0.0_dp ) then

                   row3_start =max(1 -step3, fa_start3)
                   row2_start =max(1 -step2, fa_start2)
                   row1_start =max(1 -step1, fa_start1)

                   row3_end = min(pub_fftbox%total_pt3 - step3, &
                        fa_start3 + fa_tight_pts3 -1)
                   row2_end = min(pub_fftbox%total_pt2 - step2, &
                        fa_start2 + fa_tight_pts2 -1)
                   row1_end = min(pub_fftbox%total_pt1 - step1, &
                        fa_start1 + fa_tight_pts1 -1)

                   tc3 = row3_start - fa_start3                      ! smmd
                   do row3= row3_start, row3_end
                      tc3 =tc3+1

                      tc2 = row2_start - fa_start2                   ! smmd
                      do row2= row2_start, row2_end
                         tc2 =tc2 +1

                         tc1 = row1_start - fa_start1                ! smmd
                         do row1= row1_start, row1_end
                            tc1 =tc1 +1

                            prec_tightbox_grad(tc1, tc2, tc3) = &
                                 prec_tightbox_grad(tc1, tc2, tc3) &
                                 + fftbox_grad(row1+step1,row2+step2,row3+step3) &
                                 * precond_func_real(step1,step2,step3)


                         enddo
                      enddo
                   enddo

                endif

             enddo
          enddo
       enddo

    ! smmd: apply the correct periodic boundary conditions on the fftbox
    else

       ! cks: loop over all points selected for convolution
       do step3 = -npt3, npt3
          do step2 = -npt2, npt2
             do step1 = -npt1, npt1

                ! cks: convolve only in a spherical region around central point
                ! cks: and only for the points of the tight fft_box of fa
                ! cks: which is much smaller than the fft box
                if ( precond_func_real(step1,step2,step3) .ne. 0.0_dp ) then

                   if (pub_fftbox%coin3) then
                      row3_start = fa_start3
                      row3_end = fa_start3 + fa_tight_pts3 -1
                   else
                      row3_start =max(1 -step3, fa_start3)
                      row3_end = min(pub_fftbox%total_pt3 - step3, &
                           fa_start3 + fa_tight_pts3 -1)
                   endif

                   if (pub_fftbox%coin2) then
                      row2_start = fa_start2
                      row2_end = fa_start2 + fa_tight_pts2 -1
                   else
                      row2_start =max(1 -step2, fa_start2)
                      row2_end = min(pub_fftbox%total_pt2 - step2, &
                           fa_start2 + fa_tight_pts2 -1)
                   endif

                   if (pub_fftbox%coin1) then
                      row1_start = fa_start1
                      row1_end = fa_start1 + fa_tight_pts1 -1
                   else
                      row1_start =max(1 -step1, fa_start1)
                      row1_end = min(pub_fftbox%total_pt1 - step1, &
                           fa_start1 + fa_tight_pts1 -1)
                   endif

                   tc3 = 0
                   do row3 = row3_start, row3_end
                      tc3 = tc3 + 1
                      if (pub_fftbox%coin3) then
                         period_pt3 = modulo(row3+step3-1,pub_fftbox%total_pt3)+1
                      else
                         period_pt3 = row3 + step3
                      endif

                      tc2 = 0
                      do row2 = row2_start, row2_end
                         tc2 = tc2 + 1
                         if (pub_fftbox%coin2) then
                            period_pt2 = modulo(row2+step2-1,pub_fftbox%total_pt2)+1
                         else
                            period_pt2 = row2 + step2
                         endif

                         tc1 = 0
                         do row1 = row1_start, row1_end
                            tc1 = tc1 + 1
                            if (pub_fftbox%coin1) then
                               period_pt1 = modulo(row1+step1-1,pub_fftbox%total_pt1)+1
                            else
                               period_pt1 = row1 + step1
                            endif

                            prec_tightbox_grad(tc1, tc2, tc3) = &
                                 prec_tightbox_grad(tc1, tc2, tc3) &
                                 + fftbox_grad(period_pt1,period_pt2,period_pt3) &
                                 * precond_func_real(step1,step2,step3)


                         enddo
                      enddo
                   enddo

                endif

             enddo
          enddo
       enddo

    endif


    ! cks: copy the tightbox with the preconditioned gradient to
    ! cks: the appropriate position in the fftbox that contained the
    ! cks: original gradient.
    call basis_copy_tightbox_to_fftbox(fftbox_grad, &                    !output
         fa_start1, fa_start2, fa_start3, &                               !input
         prec_tightbox_grad, fa_tight_pts1, fa_tight_pts2, fa_tight_pts3) !input


    ! cks: free up memory
    deallocate(prec_tightbox_grad,stat=ierr)
    call utils_dealloc_check('ngwf_grad_precond_rspace','prec_tightbox_grad',&
         ierr)

    call timer_clock("ngwf_grad_precond_rspace", 2)


  end subroutine ngwf_grad_precond_rspace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_precond_recip(data_box)

    !======================================================================!
    ! Performs kinetic energy preconditioning in reciprocal space on a     !
    ! function inside an FFTbox.                                           !
    !----------------------------------------------------------------------!
    ! Modified on 08/08/2003 by aam.                                       !
    ! Modified on 28/12/2003 by cks.                                       !
    ! Modified on 15/7/2004 by pdh.                                        !
    !======================================================================!

    use simulation_cell, only: pub_fftbox

    implicit none

    complex(kind=DP), intent(inout) :: data_box(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)

    data_box = data_box * precond_func_recip

  end subroutine ngwf_gradient_precond_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_init_precond_recip

    !======================================================================!
    ! Initialises the preconditioning function in reciprocal space.        !
    !----------------------------------------------------------------------!
    ! Originally written in 2003.                                          !
    !======================================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, stdout
    use rundat, only: precond_scheme, smooth_scheme, k_zero
    use simulation_cell, only: pub_fftbox
    use utils, only: utils_alloc_check

    implicit none

    ! Local variables
    integer :: ierr
    integer :: i1,i2,i3,i1start
    real(kind=DP) :: scale
    real(kind=DP) :: x,temp

    ! cks: test to avoid divisions by zero and other problems
    if (k_zero <= 0.0_DP) then
       if (pub_on_root) write(stdout,'(a,f8.3)') &
            'Error in ngwf_gradient_precond_recip: invalid k_zero =',k_zero
       call comms_abort
    end if

    allocate(precond_func_recip(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('ngwf_grad_init_precond_recip', &
         'precond_func_recip',ierr)

    scale = 1.0_DP / (0.5_DP * k_zero * k_zero)
    precond_func_recip = 0.0_DP

    select case (precond_scheme)
    case ('BG')    ! B-G method
       do i3=1,pub_fftbox%total_pt3
          do i2=1,pub_fftbox%total_pt2
             do i1=1,pub_fftbox%total_pt1
                x = scale * pub_fftbox%recip_grid(5,i1,i2,i3)
                if (smooth_scheme .eq. 'NONE' .or. smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = 1.0_DP / (1.0_DP + x)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(1.0_DP / (1.0_DP + x))
                endif
             end do
          end do
       end do
    case ('MAURI') ! Mauri method
       precond_func_recip(1,1,1) = 1.0_DP   ! G=0
       i1start = 2                          ! miss G=0
       do i3=1,pub_fftbox%total_pt3
          do i2=1,pub_fftbox%total_pt2
             do i1=i1start,pub_fftbox%total_pt1
                x = scale * pub_fftbox%recip_grid(5,i1,i2,i3)
                if (smooth_scheme .eq. 'NONE' .or. smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = min(1.0_DP,1.0_DP/x)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(min(1.0_DP,1.0_DP/x))
                endif
             end do
             i1start = 1
          end do
       end do
    case ('TETER') ! Teter-Allen-Payne method
       do i3=1,pub_fftbox%total_pt3
          do i2=1,pub_fftbox%total_pt2
             do i1=1,pub_fftbox%total_pt1
                x = scale * pub_fftbox%recip_grid(5,i1,i2,i3)
                temp = 27.0_DP + x * (18.0_DP + x * (12.0_DP + 8.0_DP * x))
                if (smooth_scheme .eq. 'NONE' .or. smooth_scheme .eq. 'GAUSSIAN') then
                   precond_func_recip(i1,i2,i3) = temp / (temp + 16.0_DP * x**4)
                else
                   precond_func_recip(i1,i2,i3) = sqrt(temp / (temp + 16.0_DP * x**4))
                endif
             end do
          end do
       end do
    case ('NONE')
       precond_func_recip = 1.0_DP
    case default
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in ngwf_grad_init_precond_recip: &
            &preconditioning scheme "',trim(precond_scheme), &
            '" not recognised'
       call comms_abort
    end select

  end subroutine ngwf_grad_init_precond_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_grad_init_smooth_recip

    !======================================================================!
    ! Initialises the preconditioning function in reciprocal space.        !
    !----------------------------------------------------------------------!
    ! Originally written in 2003.                                          !
    !======================================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, stdout
    use rundat, only: precond_scheme, smooth_scheme, k_smooth
    use simulation_cell, only: pub_fftbox
    use utils, only: utils_alloc_check

    implicit none

    ! Local variables
    integer :: ierr
    integer :: i1,i2,i3,i1start
    real(kind=DP) :: scale
    real(kind=DP) :: x,temp

    if (smooth_scheme .ne. 'SHAVE') then
       allocate(smooth_func_recip(pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
            pub_fftbox%total_pt3),stat=ierr)
       call utils_alloc_check('ngwf_grad_init_smooth_recip', &
            'smooth_func_recip',ierr)

       scale = 1.0_DP / (0.5_DP * k_smooth * k_smooth)
       smooth_func_recip = 0.0_DP

       select case (smooth_scheme)
       case ('BG')    ! B-G method
          do i3=1,pub_fftbox%total_pt3
             do i2=1,pub_fftbox%total_pt2
                do i1=1,pub_fftbox%total_pt1
                   x = scale * pub_fftbox%recip_grid(5,i1,i2,i3)
                   smooth_func_recip(i1,i2,i3) = sqrt(1.0_DP / (1.0_DP + x))
                end do
             end do
          end do
       case ('MAURI') ! Mauri method
          smooth_func_recip(1,1,1) = 1.0_DP   ! G=0
          i1start = 2                          ! miss G=0
          do i3=1,pub_fftbox%total_pt3
             do i2=1,pub_fftbox%total_pt2
                do i1=i1start,pub_fftbox%total_pt1
                   x = scale * pub_fftbox%recip_grid(5,i1,i2,i3)
                   smooth_func_recip(i1,i2,i3) = sqrt(min(1.0_DP,1.0_DP/x))
                end do
                i1start = 1
             end do
          end do
       case ('TETER') ! Teter-Allen-Payne method
          do i3=1,pub_fftbox%total_pt3
             do i2=1,pub_fftbox%total_pt2
                do i1=1,pub_fftbox%total_pt1
                   x = scale * pub_fftbox%recip_grid(5,i1,i2,i3)
                   temp = 27.0_DP + x * (18.0_DP + x * (12.0_DP + 8.0_DP * x))
                   smooth_func_recip(i1,i2,i3) = sqrt(temp / (temp + 16.0_DP * x**4))
                end do
             end do
          end do
       case ('NONE')
          smooth_func_recip = 1.0_DP
       case default
          if (pub_on_root) write(stdout,'(3a)') &
               'Error in ngwf_grad_init_smooth_recip: &
               &smoothing scheme "',trim(smooth_scheme), &
               '" not recognised'
          call comms_abort
       end select
    endif

  end subroutine ngwf_grad_init_smooth_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_exit

    !===========================================!
    ! Free up allocated memory from the module. !
    !-------------------------------------------!
    ! Written by Arash A. Mostofi on 10/7/2003  !
    !===========================================!

    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr

    if (allocated(precond_func_real)) then
       deallocate(precond_func_real,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_exit','precond_func_real',ierr)
    end if

    if (allocated(precond_func_recip)) then
       deallocate(precond_func_recip,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_exit','precond_func_recip',ierr)
    end if

  end subroutine ngwf_gradient_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_gradient_coeffs(qmat,tc_qmat,pur_denskern,tc_pur_denskern, &
       rq,tc_rq,rk,tc_rk,hub_wk,tc_hub_wk,cond_grad,tc_cond_grad,&
       denskern,overlap,ham,inv_overlap,ps_overlap,sp_overlap, &
       hub_projector_ham,hub_overlap_t,mu,n_occ,dijhat, &
       val_dkn, val_olap, val_sp_olap, val_ham, cross_overlap, cond_shift)

    !========================================================================!
    ! This subroutine returns the various matrices which multiply the        !
    ! NGWFs in row-sums required for the NGWF gradient                       !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    !  qmat (inout) : Q matrix in SPAM3 storage.                             !
    !  tc_qmat (inout) : Tensor-corrected Q matrix, QS, in SPAM3 storage.    !
    !  pur_denskern (inout) : Purified density kernel K in SPAM3 storage     !
    !  tc_pur_denskern (inout) : Tensor-corrected purified density kernel    !
    !    KS in SPAM3 storage.                                                !
    !  rk (inout) : Coefficients of nonlocal projectors in gradient in SPAM3 !
    !    storage.                                                            !
    !  rsk (inout) : Tensor-corrected coefficients of nonlocal projectors in !
    !    gradient in SPAM3 storage.                                          !
    !  hub_wk (inout) : Coefficients of DFT+U projectors in gradient in      !
    !  SPAM3 storage.                                                        !
    !  tc_hub_wk (inout) : Coefficients of DFT+U projectors in gradient in   !
    !  SPAM3 storage.                                                        !
    !  ham (input) : Hamiltonian matrix in SPAM3 format.                     !
    !  denskern (input) : Density kernel in SPAM3 format.                    !
    !  overlap (input) : Overlap matrix in SPAM3 format.                     !
    !  inv_overlap (input) : Inverse overlap matrix in SPAM3 format.         !
    !  mu (input) : Lagrange multiplier (original version) or fudge          !
    !    parameter (Millam-Scuseria version) of LNV function.                !
    !  n_occ (input) : number of occupied orbitals for each spin channel.    !
    !  dijhat (input) : Screened part of nonlocal projector energies         !
    !  val_dkn (input) : Valence density kernel (optional)                   !
    !  val_olap (input) : Valence overlap matrix (optional)                  !
    !  val_sp_olap (input) : Valence NGWF-proj overlap matrix (optional)     !
    !  val_ham (input) : Valence Hamiltonian matrix (optional)               !
    !  cross_overlap (input) : Valence-conduction overlap matrix (optional)  !
    !  cond_shift (input) : Value by which the projected conduction          !
    !    Hamiltonian is currently being shifted (optional)                   !
    !------------------------------------------------------------------------!
    ! Written by Nicholas D.M. Hine on 26/07/2009 reusing some code from     !
    !   previous versions of ngwf_gradient_lnv written by Chris-Kriton       !
    !   Skylaris in 2000-2007.                                               !
    ! DFT+U added by David D. O'Regan in September 2009.                     !
    ! Modified for conduction calculation by Laura Ratcliff in October 2010. !
    ! Modified by Alvaro Ruiz Serrano for kernel DIIS in November 2010.      !
    !========================================================================!

    use augmentation, only: aug_projector_denskern
    use cell_grid, only: pub_fine_grid
    use constants, only: DP, max_spins, stdout
    use function_basis, only: FUNC_BASIS
    use hubbard_build, only: hubbard_projector_ham
    use kernel, only: kernel_purify, kernel_rescale
    use kernel_diis, only: kernel_diis_build_pq
    use paw, only: paw_nonlocal_energies, paw_projector_overlap
    use pseudopotentials, only: pseudo_get_dij
    use rundat, only: exact_lnv, occ_mix, pub_hubbard, pub_any_nl_proj, task, &
         pub_paw, pub_cond_calculate,  pub_kernel_diis
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_trace, sparse_index_length, sparse_destroy,&
         sparse_create, sparse_product, sparse_scale, sparse_copy, sparse_axpy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: qmat(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: tc_qmat(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: pur_denskern(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: tc_pur_denskern(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: rq(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: tc_rq(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: rk(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: tc_rk(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: hub_wk(pub_cell%num_spins)  !ddor
    type(SPAM3), intent(inout) :: tc_hub_wk(pub_cell%num_spins) !ddor
    type(SPAM3), intent(inout) :: cond_grad(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: tc_cond_grad(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: overlap
    type(SPAM3), intent(in) :: ham(pub_cell%num_spins)
    type(SPAM3), intent(in) :: inv_overlap
    type(SPAM3), intent(in) :: ps_overlap
    type(SPAM3), intent(in) :: sp_overlap
    type(SPAM3), intent(in) :: hub_overlap_t
    real(kind=DP), intent(in) :: mu(max_spins)
    integer, intent(in) :: n_occ(max_spins)
    type(SPAM3), intent(in) :: dijhat(:)
    type(SPAM3), intent(in) :: hub_projector_ham(:)
    ! lr408: Optional conduction input arguments
    type(SPAM3), optional, intent(in) :: val_dkn(pub_cell%num_spins)
    type(SPAM3), optional, intent(in) :: val_ham(pub_cell%num_spins)
    type(SPAM3), optional, intent(in) :: val_olap
    type(SPAM3), optional, intent(in) :: val_sp_olap
    type(SPAM3), optional, intent(in) :: cross_overlap
    real(kind=dp), optional, intent(in) :: cond_shift

    ! Local Variables
    type(SPAM3) :: tc_buffer
    type(SPAM3) :: oij
    type(SPAM3), allocatable :: rhoij(:)
    type(SPAM3), allocatable :: dij(:)
    type(SPAM3), allocatable :: dij_ps_overlap(:)
    type(SPAM3) :: hub_ps_overlap                !ddor
    ! lr408: Conduction matrices
    type(SPAM3) :: ktm, kt, khv, ksv, ksvkt, khvkt, khvktm, ksvktm
    real(kind=DP) :: spin_fac
    integer :: is
    integer :: ierr


    call timer_clock('ngwf_gradient_coeffs', 1)

    ! cks: Construct purified density kernel K from density kernel L
    if (.not.pub_kernel_diis) then
       call kernel_purify(pur_denskern,denskern,overlap,inv_overlap,n_occ)
    end if

    ! ndmh: Choose means of determining Q matrix according to the procedure
    ! ndmh: used to optimise the density kernel
    if (exact_lnv.and..not.pub_kernel_diis) then

       ! form NGWF gradient matrix
       ! lr408: Extra arguments needed for conduction calculation
       if (pub_cond_calculate) then
          call ngwf_grad_matrix_lnv_robust(qmat, &
               pur_denskern, denskern, ham, overlap, mu, n_occ)
       else
          call ngwf_grad_matrix_lnv_robust(qmat, &
               pur_denskern, denskern, ham, overlap, mu, n_occ)
       end if

       ! ndmh: rescale pur_denskern for exact lnv constraint
       call kernel_rescale(pur_denskern,overlap,n_occ,silent=.true.)

    else if (pub_kernel_diis) then

       call kernel_diis_build_pq(pur_denskern, qmat, denskern, ham, overlap, &
            n_occ, mu)

    else ! neither exact_lnv nor kernel_diis

       call ngwf_grad_matrix_ms_robust(qmat, denskern, ham, overlap, mu)

    end if

    ! lr408: ++++ CONDUCTION GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++++++
    if (pub_cond_calculate) then

       call sparse_create(kt,val_dkn(1),cross_overlap)
       call sparse_create(ktm,kt,pur_denskern(1))
       call sparse_create(khv,val_dkn(1),val_ham(1))
       call sparse_create(ksv,val_dkn(1),val_olap)

       call sparse_create(khvktm,khv,ktm)
       call sparse_create(ksvktm,ksv,ktm)

       do is=1,pub_cell%num_spins

          call sparse_product(kt,val_dkn(is),cross_overlap)
          call sparse_product(ktm,kt,pur_denskern(1))
          call sparse_product(khv,val_dkn(is),val_ham(is))
          call sparse_product(ksv,val_dkn(is),val_olap)

          call sparse_product(khvktm,khv,ktm)
          call sparse_product(ksvktm,ksv,ktm)

          call sparse_scale(khvktm,-1.0_dp)

          ! lr408: -KHVKTM + wKSVKTM
          call sparse_copy(cond_grad(is),khvktm)
          call sparse_axpy(cond_grad(is),ksvktm,cond_shift)

       end do

       call sparse_destroy(ksvktm)
       call sparse_destroy(khvktm)

       call sparse_destroy(khv)
       call sparse_destroy(ksv)
       call sparse_destroy(ktm)
       call sparse_destroy(kt)

    end if
    ! lr408: ++++ END CONDUCTION GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++


    ! Apply tensor correction: the result is not symmetric!
    ! cks: Determine matrices depending on if occ-precond happens
    ! cks: or not.
    do is=1,pub_cell%num_spins
       call sparse_product(tc_pur_denskern(is),pur_denskern(is),overlap)
       call sparse_product(tc_qmat(is),qmat(is),overlap)
       if (pub_cond_calculate) then ! lr408
          call sparse_product(tc_cond_grad(is),cond_grad(is),overlap)
       end if
    end do

    ! cks: ++++ OCC-MIX (SPAM3 MATRICES) ++++++++++++++++++++++++++++++++++++
    ! cks: Mix fraction of occupation preconditioned covariant NGWF gradient
    ! cks: into covariant NGWF gradient.
    if (occ_mix /= 0.0_DP) then

       call sparse_create(tc_buffer,tc_qmat(1))

       do is=1,pub_cell%num_spins
          ! cks: add occ_mix*I to (1.0_DP-occ_mix)*K.S
          call sparse_scale(tc_pur_denskern(is),1.0_DP-occ_mix,occ_mix)

          ! cks: Set tc_buffer to occ_mix*(-1.0_DP)*S^-1.H
          call sparse_product(tc_buffer,inv_overlap,ham(is))
          call sparse_scale(tc_buffer,-occ_mix)

          ! cks: set tc_qmat to occ_mix*(-1.0_DP)*S^-1.H + (1.0_DP-occ_mix)*Q_tc
          call sparse_axpy(tc_buffer,tc_qmat(is),1.0_DP-occ_mix)

          ! lr408: add -occ_mix*wI
          !if (pub_cond_calculate) then
          !   call sparse_scale(tc_buffer,1.0_DP,-occ_mix*cond_shift)
          !end if

          call sparse_copy(tc_qmat(is),tc_buffer)
       end do

       call sparse_destroy(tc_buffer)

       if (pub_cond_calculate) then
          call sparse_create(tc_buffer,tc_cond_grad(1))

          call sparse_create(kt,val_dkn(1),cross_overlap)
          call sparse_create(khv,val_dkn(1),val_ham(1))
          call sparse_create(ksv,val_dkn(1),val_olap)

          call sparse_create(khvkt,khv,kt)
          call sparse_create(ksvkt,ksv,kt)

          do is=1,pub_cell%num_spins

             call sparse_product(kt,val_dkn(is),cross_overlap)
             call sparse_product(khv,val_dkn(is),val_ham(is))
             call sparse_product(ksv,val_dkn(is),val_olap)

             call sparse_product(khvkt,khv,kt)
             call sparse_product(ksvkt,ksv,kt)

             ! lr408: Set tc_buffer to occ_mix*(-KHvKT + wKSvKT)
             call sparse_scale(khvkt,-occ_mix)
             call sparse_copy(tc_buffer,khvkt)
             call sparse_axpy(tc_buffer,ksvkt,occ_mix*cond_shift)

             ! lr408: set tc_qmat to
             ! lr408: occ_mix*(-KHvKT + wKSvKT) + (1.0_DP-occ_mix)*Q_cond_tc
             call sparse_axpy(tc_buffer,tc_cond_grad(is),1.0_DP-occ_mix)
             call sparse_copy(tc_cond_grad(is),tc_buffer)
          end do

          call sparse_destroy(ksvkt)
          call sparse_destroy(khvkt)

          call sparse_destroy(khv)
          call sparse_destroy(ksv)
          call sparse_destroy(kt)

          call sparse_destroy(tc_buffer)
       end if

    end if
    ! cks: ++++ END OCC-MIX (SPAM3 MATRICES) ++++++++++++++++++++++++++++++++

    ! ndmh: ++++ NONLOCAL GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++++++
    if (pub_any_nl_proj.or.pub_paw) then

       ! Create the matrix of nonlocal energies D_ij
       if (pub_paw) then
          allocate(dij_ps_overlap(pub_cell%num_spins),stat=ierr)
          call utils_alloc_check('ngwf_gradient_coeffs','dij_ps_overlap',ierr)
          allocate(dij(pub_cell%num_spins),stat=ierr)
          call utils_alloc_check('ngwf_gradient_coeffs','dij',ierr)
          allocate(rhoij(pub_cell%num_spins),stat=ierr)
          call utils_alloc_check('ngwf_gradient_coeffs','rhoij',ierr)
       else
          allocate(dij_ps_overlap(1),stat=ierr)
          call utils_alloc_check('ngwf_gradient_coeffs','dij_ps_overlap',ierr)
          allocate(dij(1),stat=ierr)
          call utils_alloc_check('ngwf_gradient_coeffs','dij',ierr)
       end if

       ! Norm-Conserving PSP version
       if (pub_any_nl_proj) then

          ! Create matrix structures for D_ij and D_ij.<Proj_j|NGWF_b>
          call sparse_create(dij_ps_overlap(1),ps_overlap)
          dij(1)%structure = 'E'
          call sparse_create(dij(1))

          ! Get Kleinman-Bylander Denominators D_i (nonzero on diagonal only)
          call pseudo_get_dij(dij(1))

          ! Calculate the matrix D_i.<Proj_i|NGWF_b>
          call sparse_product(dij_ps_overlap(1),dij(1),ps_overlap)

          ! Calculate the matrices 4*(D_ij.<Proj_j|NGWF_b>).K^ab
          ! and 4*(D_ij.<Proj_j|NGWF_b>).K^bc.S_ca
          spin_fac = 2.0_DP / pub_cell%num_spins
          do is=1,pub_cell%num_spins
             call sparse_product(rk(is),dij_ps_overlap(1),pur_denskern(is))
             call sparse_scale(rk(is),2.0_DP*spin_fac)
             call sparse_product(tc_rk(is),dij_ps_overlap(1), &
                  tc_pur_denskern(is))
             call sparse_scale(tc_rk(is),2.0_DP*spin_fac)
          end do

          ! Clean up temporary matrices
          call sparse_destroy(dij(1))
          call sparse_destroy(dij_ps_overlap(1))

       end if

       ! PAW version
       if (pub_paw) then

          ! Create matrix structures for rho_ij, D_ij and D_ij.<Proj_j|NGWF_b>
          do is=1,pub_cell%num_spins
             call sparse_create(dij_ps_overlap(is),ps_overlap)
             rhoij(is)%structure = 'E'
             call sparse_create(rhoij(is))
             call sparse_create(dij(is),rhoij(is))
          end do
          call sparse_create(oij,rhoij(1))

          ! Create projector density kernel
          if (pub_cond_calculate) then
             call aug_projector_denskern(rhoij,val_dkn,val_sp_olap)
          else
             call aug_projector_denskern(rhoij,pur_denskern,sp_overlap)
          end if
          do is=1,pub_cell%num_spins
             call sparse_scale(rhoij(is),pub_cell%spin_fac)
          end do

          ! Calculate nonlocal energies (block diagonal matrix)
          ! and add pre-calculated screened part
          call paw_nonlocal_energies(dij,rhoij,show_matrices=.false.)
          do is=1,pub_cell%num_spins
             call sparse_axpy(dij(is),dijhat(is),1.0_DP)
          end do

          ! Calculate projector overlap matrix O_ij
          call paw_projector_overlap(oij)

          ! Calculate the matrix O_ij.<Proj_i|NGWF_b>
          call sparse_product(dij_ps_overlap(1),oij,ps_overlap)
          call sparse_scale(dij_ps_overlap(1),2.0_DP*pub_cell%spin_fac)
          do is=1,pub_cell%num_spins
             ! Calculate the matrices 4*(O_ij.<Proj_j_NGWF_b>).Q^ab
             ! and 4*(O_ij.<Proj_j|NGWF_b>).Q^bc.S_ca
             call sparse_product(rq(is),dij_ps_overlap(1),qmat(is))
             call sparse_product(tc_rq(is),dij_ps_overlap(1),tc_qmat(is))
          end do

          do is=1,pub_cell%num_spins
             ! Calculate the matrix D_ij.<Proj_j|NGWF_b>
             call sparse_product(dij_ps_overlap(is),dij(is),ps_overlap)
             call sparse_scale(dij_ps_overlap(is),2.0_DP*pub_cell%spin_fac)

             ! Calculate the matrices 4*(D_ij.<Proj_j|NGWF_b>).K^ab
             ! and 4*(D_ij.<Proj_j|NGWF_b>).K^bc.S_ca
             call sparse_product(rk(is),dij_ps_overlap(is),pur_denskern(is))
             call sparse_product(tc_rk(is),dij_ps_overlap(is), &
                  tc_pur_denskern(is))
          end do

          ! Both sets of terms will be evaluated together in one fftbox sum, so
          ! add 4*(O_ij.<Proj_j_NGWF_b>).Q^ab to rk(:) and
          ! add 4*(O_ij.<Proj_j|NGWF_b>).Q^bc.S_ca to tc_rk(:)
          do is=1,pub_cell%num_spins
             call sparse_axpy(rk(is),rq(is),1.0_DP)
             call sparse_axpy(tc_rk(is),tc_rq(is),1.0_DP)
          end do

          ! Clean up temporary matrices
          call sparse_destroy(oij)
          do is=pub_cell%num_spins,1,-1
             call sparse_destroy(dij(is))
             call sparse_destroy(rhoij(is))
             call sparse_destroy(dij_ps_overlap(is))
          end do
          deallocate(rhoij,stat=ierr)
          call utils_dealloc_check('ngwf_gradient_coeffs','rhoij',ierr)

       end if

       ! Deallocate SPAM3 arrays
       deallocate(dij,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_coeffs','dij',ierr)
       deallocate(dij_ps_overlap,stat=ierr)
       call utils_dealloc_check('ngwf_gradient_coeffs','dij_ps_overlap',ierr)

    end if
    ! ndmh: ++++ END NONLOCAL GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++

    ! ddor: ++++ DFT+U GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++++++
    if (pub_hubbard) then

       ! Create matrix structure for [H^j_k O^ki + O^jk H_k^i] <Proj_i|NGWF_b>
       call sparse_create(hub_ps_overlap,hub_overlap_t)

       ! Calculate the matrices 4*(G <Proj_i|NGWF_b>).K^ab
       ! and 4*(G <Proj_i|NGWF_b>).S_ac.K^cb
       spin_fac = 2.0_DP / pub_cell%num_spins

       do is=1,pub_cell%num_spins

          ! Calculate the matrix G <Proj_i|NGWF_b>
          call sparse_product(hub_ps_overlap,hub_projector_ham(is), &
               hub_overlap_t)

          ! ddor: Multiply on the right with Kernel and scale
          call sparse_product(hub_wk(is),hub_ps_overlap,pur_denskern(is))
          call sparse_scale(hub_wk(is),2.0_DP*spin_fac)
          call sparse_product(tc_hub_wk(is),hub_ps_overlap,tc_pur_denskern(is))
          call sparse_scale(tc_hub_wk(is),2.0_DP*spin_fac)
       end do

       ! Clean up temporary matrices
       call sparse_destroy(hub_ps_overlap)

    endif
    ! ddor: ++++ END DFT+U GRADIENT COEFFS (SPAM3 MATRICES) ++++++++++++++

    call timer_clock('ngwf_gradient_coeffs', 2)

  end subroutine ngwf_gradient_coeffs


  subroutine ngwf_grad_matrix_ms_robust(gradient, &
       denskern, ham, overlap, mu)

    !=====================================================================!
    ! This subroutine returns the matrix that multiplies the fb-only      !
    ! part of the gradient of the LNV-MS energy as a function             !
    ! of the NGWF expansion coefficients in the psinc basis.              !
    !---------------------------------------------------------------------!
    ! Written by Arash A. Mostofi in March 2003, based on the subroutine  !
    ! ngwf_gradient_plain_matrix_for_fb_robust which was written by       !
    ! Chris-Kriton Skylaris on 21/6/2001 for the ONES code.               !
    ! Rewritten by Peter Haynes for block sparse matrices in spring 2004. !
    ! Modified for conduction calculations by Laura Ratcliff Oct 2010.    !
    !=====================================================================!

    use constants, only: DP, max_spins
    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_transpose, sparse_product, &
         sparse_scale, sparse_axpy, sparse_destroy, sparse_copy

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: gradient(pub_cell%num_spins)
    type(SPAM3), intent(in)  :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in)  :: ham(pub_cell%num_spins)
    type(SPAM3), intent(in)  :: overlap
    real(kind=DP), intent(in)  :: mu(max_spins)

    ! Local Variables
    integer :: is
    type(SPAM3) :: sl,ls,hl,lhl,ltmp

    ! Allocate workspace
    call sparse_create(sl,overlap,denskern(1))
    call sparse_create(ls,denskern(1),overlap)
    call sparse_create(hl,ham(1),denskern(1))
    call sparse_create(lhl,denskern(1),hl)
    call sparse_create(ltmp,denskern(1))

    do is=1,pub_cell%num_spins

       ! Calculate 'L' := 3*L.H.L - 2*(L.H.L.S.L + L.S.L.H.L) - mu*L
       call sparse_product(sl,overlap,denskern(is))
       call sparse_transpose(ls,sl)
       call sparse_product(hl,ham(is),denskern(is))
       call sparse_product(lhl,denskern(is),hl)
       call sparse_copy(gradient(is),denskern(is))
       call sparse_scale(gradient(is),-mu(is))
       call sparse_product(ltmp,lhl,sl)
       call sparse_axpy(gradient(is),ltmp,-2.0_DP)
       call sparse_product(ltmp,ls,lhl)
       call sparse_axpy(gradient(is),ltmp,-2.0_DP)
       call sparse_axpy(gradient(is),lhl,3.0_DP)

       ! lr408: Add extra term => - eM
       !if (present(cond_shift)) then
       !   call sparse_axpy(gradient(is),purkern(is),-cond_shift)
       !end if

    end do

    ! Deallocate workspace
    call sparse_destroy(ltmp)
    call sparse_destroy(lhl)
    call sparse_destroy(hl)
    call sparse_destroy(ls)
    call sparse_destroy(sl)

  end subroutine ngwf_grad_matrix_ms_robust



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine ngwf_grad_matrix_lnv_robust(gradient, &
       purkern, denskern, ham, overlap, mu, n_occ)

    !=====================================================================!
    ! This subroutine returns the matrix that multiplies the fb-only      !
    ! part of the gradient of the original-LNV energy as a function       !
    ! of the NGWF expansion coefficients in the psinc basis.              !
    !---------------------------------------------------------------------!
    ! Written by Arash A. Mostofi in March 2003, based on the subroutine  !
    ! ngwf_gradient_plain_matrix_for_fb_robust which was written by       !
    ! Chris-Kriton Skylaris on 21/6/2001 for the ONES code.               !
    ! Rewritten by Peter Haynes for block sparse matrices in spring 2004. !
    ! Adjusted to allow independent rescaling of different spin channels  !
    ! of LNV gradient to match other parts of the code, by Nicholas Hine  !
    ! in April 2010.                                                      !
    ! Modified for conduction calculations by Laura Ratcliff Oct 2010.    !
    !=====================================================================!

    use comms, only : comms_abort
    use constants, only: DP, max_spins
    use simulation_cell, only  : pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_transpose, sparse_product, &
         sparse_scale, sparse_axpy, sparse_destroy, sparse_trace, sparse_copy

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)   :: gradient(pub_cell%num_spins)
    type(SPAM3), intent(in)    :: purkern(pub_cell%num_spins)
    type(SPAM3), intent(in)    :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in)    :: ham(pub_cell%num_spins)
    type(SPAM3), intent(in)    :: overlap
    real(kind=DP), intent(in)  :: mu(max_spins)
    integer, intent(in) :: n_occ(max_spins)

    ! Local Variables
    integer :: is
    real(kind=DP) :: norm_fac,trks
    type(SPAM3) :: sl,hl,lhl,ktmp,htmp

    ! Allocate workspace
    call sparse_create(sl,overlap,denskern(1))
    call sparse_create(hl,ham(1),denskern(1))
    call sparse_create(lhl,denskern(1),hl)
    call sparse_create(ktmp,purkern(1))
    call sparse_create(htmp,ham(1))

    ! pdh: new LNV requires gradient to be scaled by Nocc / Tr(KS)
    ! ndmh: alteration to allow separate rescaling of each spin channel
    do is=1,pub_cell%num_spins
       trks = sparse_trace(purkern(is),overlap)
       if (trks > tiny(1.0_DP)) then
          norm_fac = n_occ(is) / trks
       else
          norm_fac = 1.0_DP
       end if

       ! Calculate 'L' := 3*L.H'.L - 2*(L.H'.L.S.L + L.S.L.H'.L) - mu*K
       ! ndmh: Removed two sparse_products and replaced them with transposes
       call sparse_copy(htmp,ham(is))
       call sparse_axpy(htmp,overlap,-mu(is))
       call sparse_product(sl,overlap,denskern(is))
       call sparse_product(hl,htmp,denskern(is))
       call sparse_product(lhl,denskern(is),hl)
       call sparse_product(gradient(is),lhl,sl)
       call sparse_scale(gradient(is),-2.0_DP)
       call sparse_transpose(ktmp,gradient(is))
       call sparse_axpy(gradient(is),ktmp,1.0_DP)
       call sparse_axpy(gradient(is),purkern(is),-mu(is))
       call sparse_axpy(gradient(is),lhl,3.0_DP)

       ! pdh: re-scale gradient
       call sparse_scale(gradient(is),norm_fac)

    end do

    ! Deallocate workspace
    call sparse_destroy(htmp)
    call sparse_destroy(ktmp)
    call sparse_destroy(lhl)
    call sparse_destroy(hl)
    call sparse_destroy(sl)

  end subroutine ngwf_grad_matrix_lnv_robust


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module ngwf_gradient
