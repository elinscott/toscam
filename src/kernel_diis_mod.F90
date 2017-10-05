! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

!==============================================================================!
! This module contains the subroutines to perform self-consistent energy       !
! minimisation in the inner loop of ONETEP in the following manner:            !
!                                                                              !
! 1) Diagonalise the Hamiltonian.                                              !
! 2) Build density kernel from the eigenvalues (enforce idempotency).          !
! 3) Perform density kernel DIIS to achieve self-consistency.                  !
!                                                                              !
! Density kernel DIIS is an alternative to LNV in the inner loop of ONETEP.    !
! Hamiltonian diagonalisation only is allowed although is unstable for large   !
! systems.                                                                     !
!------------------------------------------------------------------------------!
! *** Please note that the Hamiltonian diagonalisation and kernel DIIS         !
! functionalities are currently under development.                             !
!------------------------------------------------------------------------------!
! Written by Alvaro Ruiz Serrano in April 2010.                                !
!==============================================================================!


module kernel_diis

!CW
#ifdef GPU_SPEEDUP
 use magma
 use fortran_cuda
#endif
!END CW

  implicit none

  public :: kernel_diis_calculate
  public :: kernel_diis_lagrangian
  public :: kernel_diis_rescale
  public :: kernel_diis_build_pq
  public :: kernel_diis_ham_diag

!CW
! public :: kernel_diis_mix,kernel_diis_shift,kernel_diis_find_ientry,pub_kernel_diis_max,kernel_diis_sparse_destroy
!END CW

contains


  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                           PUBLIC ROUTINES
  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_calculate(total_energy, mu, pur_denskern, &
       denskern, ham, lhxc_fine, ngwf_basis, hub_proj_basis, hub, rep, &
       localpseudo_fine, core_density_fine, ewald_energy, elements)

    !==========================================================================!
    ! Subroutine that performs density kernel mixing and hamiltonian           !
    ! diagonalisation. It minimises the energy with respect to the density     !
    ! kernel. ONETEP inner loop.                                               !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                            !
    ! Updated by Alvaro Ruiz Serrano in October 2010 to combine linear+Pulay   !
    ! mixing schemes.                                                          !
    !==========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root, comms_abort, pub_total_num_nodes
!CW
    use comms, only : pub_my_node_id
!END CW
    use constants, only: stdout, DP, max_spins, VERBOSE
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use hubbard_build, only: HUBBARD_MODEL
    use kernel, only: kernel_normalise
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_max, pub_kernel_diis_maxit, &
         pub_kernel_diis_type, pub_kernel_diis_liter, pub_kernel_diis_c_in, &
         pub_kernel_diis_c_out, print_qc, pub_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_copy, sparse_scale, sparse_max_abs_element, sparse_axpy, sparse_trace
    use utils, only: utils_alloc_check, utils_dealloc_check
    use timer, only: timer_clock


    implicit none

    ! Arguments
    real(kind=DP),    intent(inout) :: total_energy
    real(kind=DP),    intent(  out) :: mu(max_spins)
    type(SPAM3),      intent(inout) :: pur_denskern(1:pub_cell%num_spins)
    type(SPAM3),      intent(inout) :: denskern(1:pub_cell%num_spins)
    type(NGWF_HAM),   intent(inout) :: ham
    real(kind=DP),    intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep
    real(kind=DP),    intent(in   ) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP),    intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP),    intent(in   ) :: ewald_energy
    type(ELEMENT),    intent(in   ) :: elements(pub_cell%nat)


    ! real
    real(kind=DP) :: deltaE
    real(kind=DP), dimension(1:pub_cell%num_spins) :: hks_skh,k_ksk,residual
    real(kind=DP), dimension(2,1:pub_cell%num_spins) :: homolumo

    ! SPAM3
    ! ars: array of input density kernels
    type(SPAM3) :: dkn_in(1:pub_cell%num_spins,1:pub_kernel_diis_max)
    ! ars: array of output density kernels
    type(SPAM3) :: dkn_out(1:pub_cell%num_spins,1:pub_kernel_diis_max)
    ! ars: array of residues
    type(SPAM3) :: residues(1:pub_cell%num_spins,1:pub_kernel_diis_max)
    ! ars: buffer to store the next input density kernel
    type(SPAM3) :: next_dkn_in(1:pub_cell%num_spins)

    ! logical
    logical :: converged, shifted

    ! integer
    integer :: iter, is, ientry


    ! ars: start timer
    call timer_clock('kernel_diis_calculate',1)

    ! ars: iteration header
    call kernel_diis_print_header()

    ! ars: print name of the linear algebra routine
#ifdef SCALAPACK
    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) &
         write(stdout,'(a21, i4, a7)') "Running ScaLAPACK on ", &
         pub_total_num_nodes, " nodes."
#else
    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) &
         write(stdout,'(a21, i4, a45)') "Running LAPACK on ", &
         pub_total_num_nodes, " nodes. Diagonalisation will occur in serial."
#endif


    ! ars: initilise variables
    converged = .false.
    shifted   = .false.
    ientry  = 1
    mu = 0.0_DP

    ! ars: check kernel-DIIS parameters
    call kernel_diis_check_params

    ! ars: rescale to keep the number of electrons constant
    call kernel_diis_rescale(denskern,rep%overlap, rep%n_occ)

    ! ars: sparse_create dkn_in and residues
    call kernel_diis_sparse_create(dkn_in, dkn_out, next_dkn_in, residues)

    ! ars: init dkn_in
    call kernel_diis_init(residues, next_dkn_in, dkn_out, dkn_in, denskern)

!CW
#ifdef GPU_SPEEDUP
    if(pub_my_node_id==0) call choose_gpu(pub_my_node_id)
#endif
!ENDCW

    !----


    ! ars: start kernel_diis loop (ONETEP inner loop)
    kernel_diis_loop:    do iter=1,pub_kernel_diis_maxit



       ! ars: build hamiltonian
       call kernel_diis_build_ham(ham, dkn_in(:,ientry), lhxc_fine, &
            ngwf_basis,hub_proj_basis, hub, rep, localpseudo_fine, &
            core_density_fine, ewald_energy, elements, ham_update=.true., &
            lhxc_fine_fixed=.false.)

       ! ars: 1) diagonalise hamiltonian
       !      2) build new kernel from the eigenvectors
       !      3) calculate mu
       call kernel_diis_ham_diag(homolumo, dkn_out(:,ientry), &
            ham%ham, rep%overlap, rep%n_occ)

       ! ars: calculate mu
       call kernel_diis_mu(mu, dkn_out(:,ientry), ham%ham, rep%overlap)

       ! ars: calculate total energy and check convergence

       call kernel_diis_convergence(total_energy, converged, &
            hks_skh, k_ksk, deltaE, residual, dkn_out, dkn_in,residues,ham,rep,&
            ngwf_basis, hub_proj_basis, hub, lhxc_fine, localpseudo_fine,&
            core_density_fine, ewald_energy, elements, mu, homolumo, iter, &
            ientry)

       ! ars: exit if appropiate
       if(converged.or.(iter.eq.pub_kernel_diis_maxit)) then
          if(pub_on_root.and..not.converged.and.pub_output_detail.ge.VERBOSE) &
               write(stdout,'(/,a40,i3,a12,/)')&
               "==> Kernel-DIIS failed to converge after", iter, " iterations."
          ! ars: copy final kernel to pur_denskern
          do is = 1, pub_cell%num_spins
             call sparse_copy(denskern(is), dkn_out(is, ientry))
             call sparse_copy(pur_denskern(is), denskern(is))
             call sparse_scale(pur_denskern(is), 2.0_DP/pub_cell%num_spins)
          end do
          exit kernel_diis_loop
       end if

       ! ars: perform DIIS
       call kernel_diis_mix(next_dkn_in, dkn_in, dkn_out, residues, &
            rep%overlap, iter, ientry)

       ! ars: shift arrays if we have reached pub_kernel_diis_max iterations
       call kernel_diis_shift(shifted, dkn_in, dkn_out, residues, iter)

       ! ars: find next position in history to work with
       call kernel_diis_find_ientry(ientry, iter, shifted)

       ! ars: copy next_dkn_in to the corresponding ientry of dkn_in
       do is=1, pub_cell%num_spins
          call sparse_copy(dkn_in(is,ientry), next_dkn_in(is))
       end do


    end do kernel_diis_loop


    !----



    ! ars: print <qc> lines
    if(print_qc) call kernel_diis_qc_output(total_energy, &
         deltaE, hks_skh, k_ksk, residual, mu, denskern, pur_denskern, ham%ham)

    ! ars: destroy residues and kernels
    call kernel_diis_sparse_destroy(dkn_in, dkn_out, next_dkn_in, residues)

    ! ars: stop timer
    call timer_clock('kernel_diis_calculate',2)


  end subroutine kernel_diis_calculate



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_lagrangian(electronic_lagrangian, overlap, denskern, &
       ham)

    !==========================================================================!
    ! Evaluates the electronic Lagrangian according to kernel-DIIS.            !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in January 2011.                          !
    !==========================================================================!

    use constants, only: DP, max_spins
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_trace, sparse_create, sparse_destroy, &
         sparse_product, sparse_axpy, sparse_copy


    implicit none


    ! Arguments
    real(kind=DP), intent(inout) :: electronic_lagrangian
    type(SPAM3),   intent(in   ) :: overlap
    type(SPAM3),   intent(in   ) :: denskern(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: ham(pub_cell%num_spins)

    ! Local variable
    integer :: is
    real(kind=DP) :: spin_fac, trhksk
    type(SPAM3) :: ksk, ks


    spin_fac = 2.0_DP / pub_cell%num_spins

    call sparse_create(ks, denskern(1), overlap)
    call sparse_create(ksk, ks, denskern(1))


    ! ars: apply kernel DIIS lagrangian
    !    : at this point, K already contains the rescaling factor
    do is =1, pub_cell%num_spins

       ! ars: calculate tr[H(KSK-K)]
       call sparse_product(ks, denskern(is), overlap)
       call sparse_product(ksk, ks, denskern(is))
       call sparse_axpy(ksk, denskern(is), -1.0_DP)
       trhksk = sparse_trace(ham(is), ksk)

       ! ars: apply Lagrangian
       electronic_lagrangian = electronic_lagrangian - spin_fac * trhksk

    end do

    call sparse_destroy(ksk)
    call sparse_destroy(ks)


  end subroutine kernel_diis_lagrangian



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_rescale(denskern, overlap, n_occ)

    !==================================================================!
    ! This subroutine scales the density kernel to preserve the number !
    ! of electrons.                                                    !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in January 2011.                   !
    !==================================================================!

    use constants, only: DP, max_spins
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_scale, sparse_trace


    implicit none

    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: overlap
    integer,     intent(in   ) :: n_occ(max_spins)

    integer :: is
    real(kind=DP) :: rescale

    do is = 1, pub_cell%num_spins
       rescale = real(n_occ(is),kind=DP) / sparse_trace(denskern(is), overlap)
       call sparse_scale(denskern(is), rescale)
    end do

  end subroutine kernel_diis_rescale



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_build_pq(pmat, qmat, denskern, ham, overlap, n_occ, mu)

    !=====================================================================!
    ! This subroutine returns the P and Q matrices that multiply the H-fb !
    ! and the fb-only parts of the gradient of the kernel DIIS energy as  !
    ! a function of the NGWF expansion coefficients in the psinc basis.   !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in January 2011.                     !
    !=====================================================================!

    use constants, only: DP, max_spins
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_copy, sparse_scale, sparse_axpy, sparse_trace


    implicit none


    ! ars: arguments
    type(SPAM3),   intent(inout) :: pmat(pub_cell%num_spins)
    type(SPAM3),   intent(inout) :: qmat(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: denskern(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: ham(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: overlap
    integer,       intent(in   ) :: n_occ(max_spins)
    real(kind=DP), intent(in   ) :: mu(max_spins)

    ! ars: local variables
    integer :: is
    type(SPAM3) :: kh, ks, ksk
    real(kind=DP) :: nu, rescale


    ! ars: allocate workspace
    call sparse_create(kh,denskern(1), ham(1))
    call sparse_create(ks,denskern(1), overlap)
    call sparse_create(ksk, ks, denskern(1))

    do is = 1, pub_cell%num_spins

       ! ars: build buffers
       call sparse_product(kh, denskern(is), ham(is))
       call sparse_product(ks, denskern(is), overlap)
       call sparse_product(ksk, ks, denskern(is))

       ! ars: build P = [2Ne/tr(KS)] *K - [Ne/tr(KS)]^2 *KSK
       rescale = real(n_occ(is), kind=DP)/sparse_trace(ks)
       call sparse_copy(pmat(is), denskern(is))
       call sparse_scale(pmat(is), 2*rescale)
       call sparse_axpy(pmat(is), ksk, -rescale*rescale)

       ! ars: build Q = -[Ne/tr(KS)]^2 *KHK - 2(mu-[Ne/tr(KS)]^2 * nu] *K
       nu = sparse_trace(ham(is), ksk) / sparse_trace(ks)
       call sparse_product(qmat(is), kh, denskern(is))
       call sparse_scale(qmat(is), rescale*rescale)
       call sparse_axpy(qmat(is), denskern(is), 2.0_DP*(mu(is)-rescale**2*nu))
       call sparse_scale(qmat(is), -1.0_DP)

    end do

    ! ars: deallocate workspace
    call sparse_destroy(kh)
    call sparse_destroy(ks)
    call sparse_destroy(ksk)


  end subroutine kernel_diis_build_pq


  !_____________________________________________________________________________
  !_____________________________________________________________________________
  !
  !                           PRIVATE ROUTINES
  !_____________________________________________________________________________
  !_____________________________________________________________________________




  subroutine kernel_diis_mix(next_dkn_in, dkn_in, dkn_out, residues, overlap, &
       iter, ientry)

    !==========================================================================!
    ! Selects method for kernel DIIS.                                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in November 2010.                         !
    !==========================================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: stdout
    use rundat, only: pub_kernel_diis_type, pub_kernel_diis_max, &
         pub_kernel_diis_liter
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_copy

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: dkn_in(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(in   ) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(in   ) :: residues(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(in   ) :: overlap
    integer,     intent(in   ) :: iter
    integer,     intent(in   ) :: ientry


    integer :: is

    select case(pub_kernel_diis_type)
    case('D')
       ! ars: do diagonalisation only
       do is=1, pub_cell%num_spins
          call sparse_copy(next_dkn_in(is), dkn_out(is, ientry))
       end do

    case('L')
       ! ars: do linear mixing only
       call kernel_diis_linear_mixing(next_dkn_in, dkn_in, &
            dkn_out, ientry)

    case('P')
       !ars: do Pulay mixing. Optional linear mixing for the first iterations.
       if(iter.le.pub_kernel_diis_liter) then
          call kernel_diis_linear_mixing(next_dkn_in, dkn_in, &
               dkn_out, ientry)
       else
          call kernel_diis_pulay_mixing(next_dkn_in, &
               dkn_out, residues, overlap, ientry)
       end if

    case default
       if(pub_on_root) write(stdout,*) "Unknown mixing method. ONETEP stops."
       call comms_abort

    end select

  end subroutine kernel_diis_mix



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~ Linear mixing routines ~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine kernel_diis_linear_mixing(next_dkn_in, dkn_in, dkn_out, ientry)

    !==========================================================================!
    ! Performs linear mixing of density kernels.                               !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in June 2010.                             !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, VERBOSE
    use rundat, only: pub_kernel_diis_max, pub_kernel_diis_c_in, &
         pub_kernel_diis_c_out, pub_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: dkn_in(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(in   ) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    integer,     intent(in   ) :: ientry



    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       write(stdout,'(/,a)') "... PERFORMING LINEAR MIXING ..."
       write(stdout,'(a9,f6.4,a10,f6.4/)') &
            " C_in = ", pub_kernel_diis_c_in, ", C_out = ", pub_kernel_diis_c_out
    end if


    ! ars: mix kernels linearly
    call kernel_diis_linear_mix_kernels(next_dkn_in, dkn_out, &
         dkn_in, ientry)


  end subroutine kernel_diis_linear_mixing



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_linear_mix_kernels(next_dkn_in, dkn_out, dkn_in, ientry)

    !==================================================================!
    ! This subroutine performs a linear mixing of the density kernels. !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in June 2010.                      !
    !==================================================================!

    use constants, only: DP
    use rundat, only: pub_kernel_diis_c_in, pub_kernel_diis_c_out, &
         pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_scale

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(in   ) :: dkn_in(pub_cell%num_spins,pub_kernel_diis_max)
    integer,     intent(in   ) :: ientry

    integer :: is


    do is = 1, pub_cell%num_spins
       ! ars: K_in(n+1) = 0
       call sparse_scale(next_dkn_in(is), 0.0_DP)
       ! ars: K_in(n+1) = c_out*K_out(n)
       call sparse_axpy(next_dkn_in(is), dkn_out(is,ientry), &
            pub_kernel_diis_c_out)
       ! ars: K_in(n+1) = c_out*K_out(n) + c_in*K_in(n)
       call sparse_axpy(next_dkn_in(is), dkn_in(is,ientry), &
            pub_kernel_diis_c_in)
    end do

  end subroutine kernel_diis_linear_mix_kernels



  !_____________________________________________________________________________
  !_____________________________________________________________________________





  !~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~ Pulay mixing routines ~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine kernel_diis_pulay_mixing(next_dkn_in,dkn_out,residues,overlap,ientry)

    !==========================================================================!
    ! Performs Pulay mixing of density kernels.                                !
    ! Based on Pulay, Chem. Phys. Lett. 73, 2, 1980.                           !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in June 2010.                             !
    ! Updated by Alvaro Ruiz Serrano in October 2010 to combine linear+Pulay   !
    ! mixing schemes.                                                          !
    !==========================================================================!


    use comms, only: pub_on_root
    use constants, only: DP, stdout, VERBOSE
    use rundat, only: pub_kernel_diis_max, pub_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check


    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(in   ) :: residues(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(in   ) :: overlap
    integer,     intent(in   ) :: ientry


    ! ars: local variables
    real(kind=DP), allocatable :: Bmat(:,:,:)!ars: linear system to be solved
    real(kind=DP), allocatable :: coeffs(:,:)!ars: Pulay mixing coeffs
    integer :: ierr


    if(pub_on_root.and.pub_output_detail.ge.VERBOSE) then
       write(stdout,'(/,a)') "... PERFORMING PULAY MIXING ..."
       write(stdout,'(a)') "Coefficients and Lagrange multiplier:"
    end if

    ! ars: allocate Bmat and coeffs
    allocate(Bmat(1:pub_cell%num_spins,1:ientry+1,1:ientry+1),stat=ierr)
    call utils_alloc_check('kernel_diis_mix_density', 'Bmat', ierr)
    allocate(coeffs(1:pub_cell%num_spins,1:ientry+1), stat=ierr)
    call utils_alloc_check('kernel_diis_mix_density', 'coeffs', ierr)

    ! ars: calculate B matrix.
    call kernel_diis_pulay_calc_Bmatrix(Bmat, residues, overlap, ientry)

    ! ars: solve linear system Bd=0 (including Lagrange mult)
    call kernel_diis_pulay_find_coeffs(coeffs, Bmat, ientry)
    if(pub_output_detail.ge.VERBOSE) &
         call kernel_diis_pulay_coeffs_summ(coeffs, ientry)

    ! ars: mix kernel
    call kernel_diis_pulay_mix_kernels(next_dkn_in, dkn_out, coeffs, ientry)

    ! ars: deallocate Bmat and coeffs
    deallocate(coeffs, stat=ierr)
    call utils_dealloc_check('kernel_diis_mix_density', 'coeffs', ierr)
    deallocate(Bmat, stat=ierr)
    call utils_dealloc_check('kernel_diis_mix_density', 'Bmat', ierr)


  end subroutine kernel_diis_pulay_mixing



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_pulay_calc_Bmatrix(Bmat,residues, overlap, ientry)

    !==============================================================!
    ! This subroutine calculates matrix B of the system of linear  !
    ! equations to be solved in order to mix the density kernels   !
    ! according to the Pulay method.                               !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in April 2010.                 !
    !==============================================================!

    use constants, only: stdout, DP
    use comms, only: pub_on_root, pub_my_node_id, comms_barrier
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_transpose, sparse_product,&
         sparse_create, sparse_destroy, sparse_trace


    implicit none

    ! ars: arguments
    real(kind=DP), intent(  out) :: Bmat(:,:,:)
    type(SPAM3),   intent(in   ) :: residues(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3),   intent(in   ) :: overlap
    integer,       intent(in   ) :: ientry

    ! ars: local variables
    integer     :: is, ii, jj
    type(SPAM3) :: spam3buffer1, spam3buffer2, spam3buffer3


    ! ars: init spam3buffer blocking scheme
    call sparse_create(spam3buffer1, residues(1,1), overlap)
    call sparse_create(spam3buffer2, residues(1,1))
    call sparse_create(spam3buffer3, spam3buffer1, spam3buffer2)
  

!CW
    if(pub_my_node_id==0) write(*,*) 'COMPUTING PULLAY COEFS'
!END CW
 
    do jj = 1, ientry
       do ii = 1, ientry
          do is = 1, pub_cell%num_spins

             ! ars: spam3buffer1 = R_i x S
!CW
#ifndef GPU_SPEEDUP
             call sparse_product(spam3buffer1, residues(is,ii), overlap)
#else
             call sparse_product_gpu(spam3buffer1, residues(is,ii), overlap)
#endif
!END CW
             ! ars: spam3buffer2 = (R_j)^t
             call sparse_transpose(spam3buffer2,residues(is,jj))

             ! ars: spam3buffer3 = (R_i) x S x (R_j)^t
!CW
#ifndef GPU_SPEEDUP
             call sparse_product(spam3buffer3, spam3buffer1, spam3buffer2)
#else
             call sparse_product_gpu(spam3buffer3, spam3buffer1, spam3buffer2)
#endif
!END CW
             ! ars: Bmat(is,ii,jj) = tr[(R_i) x S x (R_j)^t x S]
             Bmat(is,ii,jj) = sparse_trace(spam3buffer3,overlap)

          end do
       end do
    end do

    ! ars: destroy SPAM3 buffers
    call sparse_destroy(spam3buffer1)
    call sparse_destroy(spam3buffer2)
    call sparse_destroy(spam3buffer3)

    ! ars: add space for Langrange multiplier
    Bmat(:,1:ientry,ientry+1) = -1.0_DP
    Bmat(:,ientry+1,1:ientry) = -1.0_DP
    Bmat(:,ientry+1,ientry+1) =  0.0_DP

  contains

!CW
      subroutine sparse_product_gpu(matc,mata,matb)
#ifdef GPU_SPEEDUP
      use fortran_cuda, only : matmulcuda_r,matmulcuda_r_singleprec_bypass
      use comms, only: comms_bcast,pub_my_node_id
#endif
      use sparse, only: sparse_convert,sparse_num_cols,sparse_num_rows
      implicit none
      integer                                 :: i,nn,n1,n2,n3
      type(SPAM3),intent(in)                  :: mata,matb
      type(SPAM3),intent(inout)               :: matc
      real(8), allocatable, dimension(:,:)    :: matdensa,matdensb,matdensc
      n1 = NINT(sparse_num_rows(mata))  ! row A
      n2 = NINT(sparse_num_cols(mata))  ! sum
      n3 = NINT(sparse_num_cols(matb))  ! col B
      allocate(matdensa(n1,n2),matdensb(n2,n3),matdensc(n1,n3))
      call sparse_convert(matdensa,mata)
      call sparse_convert(matdensb,matb)
#ifdef GPU_SPEEDUP
      if(pub_my_node_id==0) call matmulcuda_r(matdensa,matdensb,matdensc,n1,n2,n3)
      call comms_bcast(0,matdensc)
#endif
      call sparse_convert(matc,matdensc)
      deallocate(matdensa,matdensb,matdensc)
      end subroutine
!END CW

  end subroutine kernel_diis_pulay_calc_Bmatrix



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_pulay_find_coeffs(coeffs,Bmat,ientry)

    !==============================================================!
    ! This subroutine finds the coefficients for the Pulay mixing  !
    ! after solving a system of linear equations Bd=0.             !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in May 2010.                   !
    !==============================================================!

    use comms, only: pub_my_node_id, pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use simulation_cell, only: pub_cell

!CW
    use comms , only : comms_bcast
!END CW

    implicit none

    ! ars: arguments
    real(kind=DP), intent(inout) :: coeffs(:,:)
    real(kind=DP), intent(in   ) :: Bmat(:,:,:)
    integer      , intent(in   ) :: ientry

    ! ars: library wrapper variables
    integer :: INFO, LDA, LDB, M, N, NRHS
    integer :: IPIV(1:ientry+1)
    real(kind=DP) :: A(ientry+1,ientry+1)
    real(kind=DP) :: B(ientry+1)
    character(LEN=1) :: TRANS

    ! ars: local variables
    integer :: is



    ! ars: set up parameters for LAPACK DGETRF
    M = ientry+1
    N = ientry+1
    LDA = ientry+1
    ! ars: set up parameters for LAPACK DGETRS
    TRANS ='N'
    NRHS = 1
    LDB = ientry+1

    ! ars: solve linear system
    do is = 1, pub_cell%num_spins

       ! ars: call LAPACK DGETRF
       A = Bmat(is,:,:)

       call DGETRF(M, N, A, LDA, IPIV, INFO)

       if((INFO.ne.0).and.pub_on_root) write(stdout,*) &
            "Error calling DGETRF. INFO|node = ", INFO, pub_my_node_id

       ! ars: call LAPACK DGETRS
       B(1:ientry) = 0.0_DP
       B(ientry+1) = -1.0_DP

       call DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)

       if((INFO.ne.0).and.pub_on_root) write(stdout,*) &
            "Error calling DGETRS. INFO|node = ", INFO, pub_my_node_id

       ! ars: set coeffs(is) before exit
       coeffs(is,:) = B(:)

    end do

    call comms_barrier



  end subroutine kernel_diis_pulay_find_coeffs



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_pulay_coeffs_summ(coeffs, ientry)

    !==================================================================!
    ! This subroutine prints a summary of the coefficients that are    !
    ! used during the density kernel mixing.                           !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in July 2010.                      !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use simulation_cell, only: pub_cell

    implicit none

    real(kind=DP), intent(in   ) :: coeffs(:,:)
    integer,       intent(in   ) :: ientry

    integer :: is, it

    if(pub_on_root) then
       do it = 1, ientry+1
          do is = 1, pub_cell%num_spins

             if(it.le.ientry) then
                write(stdout,'(a6,i1,a1,i2,a4,f20.12)') &
                     "Coeff(",is,",",it,") = ", coeffs(is,it)
             else
                write(stdout,'(a7,i1,a6,f20.12,/)') &
                     "Lambda(",is,")   = ", coeffs(is,it)
             end if

          end do
       end do
    end if

  end subroutine kernel_diis_pulay_coeffs_summ



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_pulay_mix_kernels(next_dkn_in, dkn_out, coeffs, ientry)

    !==================================================================!
    ! This subroutine mixes the kernels according to the Pulay method. !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in May 2010.                       !
    !==================================================================!

    use constants, only: DP
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_scale

    implicit none

    type(SPAM3),   intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    real(kind=DP), intent(in   ) :: coeffs(:,:)
    integer,       intent(in   ) :: ientry

    integer :: is, counter

    ! ars: K_in(n+1) = 0
    do is = 1, pub_cell%num_spins
       call sparse_scale(next_dkn_in(is),0.0_DP)
    end do

    ! ars: K_in(n+1) = sum_i K_out(i) * d_i
    do counter = 1, ientry
       do is = 1, pub_cell%num_spins
          call sparse_axpy(next_dkn_in(is),dkn_out(is,counter),coeffs(is,counter))
       end do
    end do


  end subroutine kernel_diis_pulay_mix_kernels



  !_____________________________________________________________________________
  !_____________________________________________________________________________





  !~~~~~~~~~~~~~~~~~~
  !~~ Common routines
  !~~~~~~~~~~~~~~~~~~
  subroutine kernel_diis_sparse_create(dkn_in, dkn_out, next_dkn_in, residues)

    !===================================================================!
    ! This subroutine creates the SPAM3 structures for dkn_in,   !
    ! dkn_out, next_dkn_in, residues and hamiltonian_array !
    !-------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                     !
    !===================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_create, SPAM3

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: dkn_in(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: residues(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)

    ! ars: local variables
    integer :: is, iter


    ! ars: create SPAM3 arrays
    do iter = 1, pub_kernel_diis_max
       do is = 1, pub_cell%num_spins
          dkn_in(is,iter)%structure = 'K'
          call sparse_create(dkn_in(is,iter))
          dkn_out(is,iter)%structure = 'K'
          call sparse_create(dkn_out(is,iter))
          residues(is,iter)%structure = 'K'
          call sparse_create(residues(is,iter))
       end do
    end do

    do is=1, pub_cell%num_spins
       next_dkn_in%structure = 'K'
       call sparse_create(next_dkn_in(is))
    end do


  end subroutine kernel_diis_sparse_create



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_sparse_destroy(dkn_in, dkn_out, next_dkn_in, residues)

    !======================================================================!
    ! This subroutine destroys the SPAM3 structures for for dkn_in, !
    ! dkn_out, next_dkn_in, residues and hamiltonian_array.   !
    !----------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in April 2010.                         !
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: stdout
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_destroy, SPAM3


    implicit none


    ! ars: arguments
    type(SPAM3), intent(inout) :: dkn_in(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: residues(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)

    ! ars: local variables
    integer :: is, iter


    ! ars: destroy SPAM3 arrays
    do iter = 1, pub_kernel_diis_max
       do is = 1, pub_cell%num_spins
          call sparse_destroy(dkn_in(is,iter))
          call sparse_destroy(dkn_out(is,iter))
          call sparse_destroy(residues(is,iter))
       end do
    end do

    do is=1, pub_cell%num_spins
       call sparse_destroy(next_dkn_in(is))
    end do


  end subroutine kernel_diis_sparse_destroy




  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_init(residues, next_dkn_in, dkn_out, dkn_in, denskern)

    !====================================================!
    ! This subroutine initialises the first element of   !
    ! dkn_in by copying the latest density kernel !
    ! onto it (coming from the last NGWF iteration).     !
    !----------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in April 2010.       !
    !====================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, DP
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_copy, sparse_scale

    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: residues(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: next_dkn_in(pub_cell%num_spins)
    type(SPAM3), intent(inout) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: dkn_in(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(in   ) :: denskern(pub_cell%num_spins)

    ! ars: local variables
    integer :: is, it


    ! ars: initialise all SPAM3 matrices to zero
    do it = 1, pub_kernel_diis_max
       do is = 1, pub_cell%num_spins
          call sparse_scale(residues(is,it), 0.0_DP)
          call sparse_scale(dkn_in(is,it), 0.0_DP)
          call sparse_scale(dkn_out(is,it), 0.0_DP)
       end do
    end do

    do is = 1, pub_cell%num_spins
       call sparse_scale(next_dkn_in(is), 0.0_DP)
    end do

    ! ars: initialise first ientry of dkn_in
    do is = 1, pub_cell%num_spins
       call sparse_copy(dkn_in(is,1),denskern(is))
    end do

  end subroutine kernel_diis_init



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_build_ham(ham, denskern, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, rep, localpseudo_fine, &
       core_density_fine, ewald_energy, elements, ham_update, lhxc_fine_fixed)

    !==================================================================!
    ! This subroutine calculates the density dependent matrices and    !
    ! builds the Hamiltonian matrix in NGWF representation.            !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in June 2010.                      !
    !==================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: DP, max_spins, paw_en_size
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_dep_matrices,&
         hamiltonian_build_matrix
    use ion, only: ELEMENT
    use hubbard_build, only: HUBBARD_MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_scale, sparse_trace, sparse_max_abs_element

!CW
    use rundat, only : pub_dmft_spoil_kernel
!END CW


    implicit none

    ! ars: arguments
    type(NGWF_HAM),   intent(inout) :: ham
    type(SPAM3),      intent(inout) :: denskern(pub_cell%num_spins)
    real(kind=DP),    intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep
    real(kind=DP),    intent(in   ) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP),    intent(in   ) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP),    intent(in   ) :: ewald_energy
    logical,          intent(in   ) :: ham_update
    logical,          intent(in   ) :: lhxc_fine_fixed
    type(ELEMENT),    intent(in   ) :: elements(pub_cell%nat)

    ! ars: local variables
    real(kind=DP) :: spin_fac
    integer :: is, nspins
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: total_energy
    real(kind=DP) :: lhxc_energy
    real(kind=DP) :: hubbard_energy

!CW
    logical, save :: updated=.false.
!END CW


    ! ars: initialise
    nspins = pub_cell%num_spins
    spin_fac = 2.0_DP/nspins

    ! ars: scale kernel
    do is=1,nspins
       call sparse_scale(denskern(is), spin_fac)
    end do

    ! ars: build hamiltonian:
    !   1: dens_indep_matrices already initialised
    !   2: initiliase density independent matrices
    call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
         ngwf_basis, hub_proj_basis, hub, denskern, ewald_energy, elements, &
         pub_fine_grid, localpseudo_fine, core_density_fine, ham_update, &
         lhxc_fine_fixed)

    !   3: build hamiltonian matrix
    call hamiltonian_build_matrix(ham, rep)

!CW
   if(pub_dmft_spoil_kernel.and..not.updated)then
     updated=.true.
     call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
         ngwf_basis, hub_proj_basis, hub, denskern, ewald_energy, elements, &
         pub_fine_grid, localpseudo_fine, core_density_fine, ham_update, &
         lhxc_fine_fixed,spoil_force=.true.)
     call hamiltonian_build_matrix(ham, rep) 
   endif
!END CW

    ! ars: scale back the kernel
    do is=1,nspins
       call sparse_scale(denskern(is), 1.0_DP/spin_fac)
    end do


  end subroutine kernel_diis_build_ham



  !_____________________________________________________________________________
  !_____________________________________________________________________________




  subroutine kernel_diis_ham_diag(homolumo, denskern, ham, &
       overlap, n_occ)

    !=====================================================================!
    ! This subroutine diagonalises the hamiltonian and builds the output  !
    ! density from the eigenvectors.                                      !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010 as an adapation of the !
    ! existing routine kernel_reset_occupancies.                          !
    !=====================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, max_spins, stdout
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_eigensolve, dense_product, dense_scale
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_num_rows, sparse_trace
    use utils, only : utils_alloc_check, utils_dealloc_check
    use timer, only: timer_clock

    implicit none

    ! ars: arguments
    real(kind=DP), intent(  out) :: homolumo(2,pub_cell%num_spins)
    type(SPAM3),   intent(inout) :: denskern(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: ham(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: overlap
    integer,       intent(in   ) :: n_occ(max_spins)

    ! ars: local variables
    integer   :: ierr
    integer   :: num
    integer   :: is
    integer   :: nspins
    real(kind=DP), allocatable :: eigenvalues(:)
    type(DEM) :: denskern_dens, overlap_dens, eigs_dens, ham_dens



    ! ars: start timer
    call timer_clock('kernel_diis_ham_diag',1)

    ! Establish number of spins
    nspins = pub_cell%num_spins
    num = sparse_num_rows(denskern(1))

    allocate(eigenvalues(num),stat=ierr)
    call utils_alloc_check('kernel_diis_ham_diag','eigenvalues',ierr)


    ! ndmh: create temporary dense matrices
    call dense_create(denskern_dens,num,num)
    call dense_create(ham_dens,num,num) !ars
    call dense_create(overlap_dens,num,num)
    call dense_create(eigs_dens,num,num)

    ! cks: loop over spins
    do is=1,nspins

       eigenvalues(:) = 0.0_DP

       ! ars: expand denskern, ham and overlap to dense matrices
       call dense_convert(denskern_dens,denskern(is))
       call dense_convert(ham_dens,ham(is))
       call dense_convert(overlap_dens,overlap)

       ! ars: solve HC = eSC
       call dense_eigensolve(n_occ(is)+1,eigenvalues,ham_dens,overlap_dens, &
            1,eigs_dens)

       ! ars: calculate homolumo
       homolumo(1,is) = eigenvalues(n_occ(is))
       homolumo(2,is) = eigenvalues(n_occ(is)+1)

       ! ndmh: construct a density kernel from the first n_occ(is) eigenvectors
       call dense_product(denskern_dens,eigs_dens,eigs_dens, &
            transpose_amat=.false.,transpose_bmat=.true., &
            first_k=1,last_k=n_occ(is))

       ! ndmh: convert back to sparse matrix
       call dense_convert(denskern(is),denskern_dens)

    end do

    ! ndmh: clean up temporary matrices
    call dense_destroy(eigs_dens)
    call dense_destroy(overlap_dens)
    call dense_destroy(ham_dens) !ars
    call dense_destroy(denskern_dens)

    deallocate(eigenvalues,stat=ierr)
    call utils_dealloc_check('kernel_diis_ham_diag','eigenvalues',ierr)


    ! ars: stop timer
    call timer_clock('kernel_diis_ham_diag',2)


  end subroutine kernel_diis_ham_diag



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_mu(mu, denskern, ham, overlap)

    !=====================================================================!
    ! This subroutine calculates mu for kernel-rescaling.                 !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in January 2011.                     !
    !=====================================================================!

    use constants, only: max_spins, DP
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_trace

    implicit none

    real(kind=DP), intent(  out) :: mu(max_spins)
    type(SPAM3),   intent(in   ) :: denskern(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: ham(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: overlap

    integer :: is

    ! ars: calculate mu
    do is = 1, pub_cell%num_spins
       mu(is) = sparse_trace(denskern(is), ham(is)) / &
            sparse_trace(denskern(is), overlap)
    end do

  end subroutine kernel_diis_mu


  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_convergence(total_energy, converged, &
       hks_skh, k_ksk, deltaE, residual, dkn_out, dkn_in, residues, ham, rep,&
       ngwf_basis, hub_proj_basis, hub, lhxc_fine, localpseudo_fine,&
       core_density_fine, ewald_energy, elements, mu, homolumo, iter, ientry)

    !=====================================================================!
    ! This subroutine calculates the total energy and the constants of the!
    ! calculation. Then it checks for convergence and prints a summary.   !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in December 2010.                    !
    !=====================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: DP, max_spins
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use hubbard_build, only: HUBBARD_MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only : SPAM3


    implicit none

    ! arguments
    real(kind=DP),    intent(inout) :: total_energy
    logical,          intent(  out) :: converged
    real(kind=DP),    intent(  out) :: hks_skh(pub_cell%num_spins)
    real(kind=DP),    intent(  out) :: k_ksk(pub_cell%num_spins)
    real(kind=DP),    intent(  out) :: deltaE
    real(kind=DP),    intent(  out) :: residual(pub_cell%num_spins)
    type(SPAM3),      intent(inout) :: dkn_out(pub_cell%num_spins,&
         1:pub_kernel_diis_max)
    type(SPAM3),      intent(inout) :: dkn_in(pub_cell%num_spins,&
         1:pub_kernel_diis_max)
    type(SPAM3),      intent(inout) :: residues(pub_cell%num_spins,&
         1:pub_kernel_diis_max)
    type(NGWF_HAM),   intent(inout) :: ham
    type(NGWF_REP),   intent(in   ) :: rep
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP),    intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP),    intent(in   ) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP),    intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP),    intent(in   ) :: ewald_energy
    real(kind=DP),    intent(in   ) :: mu(max_spins)
    real(kind=DP),    intent(in   ) :: homolumo(2,pub_cell%num_spins)
    integer,          intent(in   ) :: iter
    integer,          intent(in   ) :: ientry
    type(ELEMENT),    intent(in   ) :: elements(pub_cell%nat)

    ! real
    real(kind=DP) :: ks_ne(pub_cell%num_spins)
    real(kind=DP) :: total_energy_old

    ! character
    character(len=8) :: okmessage


    ! ars: init
    total_energy_old = total_energy

    ! ars: calculate commutators and number of electrons
    call kernel_diis_constants(hks_skh, k_ksk, ks_ne, &
         dkn_in(:,ientry), ham%ham, rep%overlap)

    ! ars: calculate total energy and print summary
    call kernel_diis_energy(total_energy, ham, dkn_out(:,ientry),&
         lhxc_fine, ngwf_basis, hub_proj_basis, hub, rep, localpseudo_fine, &
         core_density_fine, ewald_energy, elements, ham_update = .false., &
         lhxc_fine_fixed=.false.)

    ! ars: calculate residue
    call kernel_diis_residue_inout(residues(:,ientry), &
         dkn_out(:,ientry), dkn_in(:,ientry))

    ! ars: calculate residual
    call kernel_diis_residual(residual, residues(:,ientry), rep%overlap)

    ! ars: calculate change in energy
    deltaE = total_energy - total_energy_old

    ! ars: check convergence
    call kernel_diis_conv_check(converged, okmessage, &
         deltaE, residual, hks_skh, k_ksk)

    ! ars: print info line
    call kernel_diis_info_line(converged, okmessage, iter, residual, &
         total_energy, homolumo, mu, hks_skh, k_ksk, ks_ne, deltaE)


  end subroutine kernel_diis_convergence



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_constants(hks_skh, k_ksk, ks_ne, denskern, ham, overlap)

    !==================================================================!
    ! This subroutine calculates the constants of the calculation      !
    ! related to the density kernel and prints a summary:              !
    !  1) [HKS,SKH]                                                    !
    !  2) [K,KSK]                                                      !
    !  3) Ne = 2tr[KS]                                                 !
    !------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in July 2010.                      !
    !==================================================================!

    use constants, only: DP, stdout
    use comms, only: pub_on_root
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_rms_element,&
         sparse_trace, sparse_product, sparse_axpy, sparse_transpose

    implicit none

    ! ars: arguments
    real(kind=DP),    intent(  out) :: hks_skh(pub_cell%num_spins)
    real(kind=DP),    intent(  out) :: k_ksk(pub_cell%num_spins)
    real(kind=DP),    intent(  out) :: ks_ne(pub_cell%num_spins)
    type(SPAM3),      intent(in   ) :: denskern(pub_cell%num_spins)
    type(SPAM3),      intent(in   ) :: ham(pub_cell%num_spins)
    type(SPAM3),      intent(in   ) :: overlap

    ! ars: local variables
    integer :: is
    type(SPAM3) :: sk, ksk, skh, hks


    ! ars: create matrices
    call sparse_create(sk, overlap, denskern(1))
    call sparse_create(ksk, denskern(1), sk)
    call sparse_create(skh, sk, ham(1))
    call sparse_create(hks, skh)


    do is = 1, pub_cell%num_spins

       ! ars: calculate matrices
       call sparse_product(sk, overlap, denskern(is))
       call sparse_product(ksk, denskern(is), sk)
       call sparse_product(skh, sk, ham(is))
       call sparse_transpose(hks, skh)

       ! ars: calculate [HKS,SKH]
       call sparse_axpy(hks, skh, -1.0_DP)
       hks_skh(is) = sparse_rms_element(hks)
       hks_skh(is) = abs(hks_skh(is))

       ! ars: calculate [K,KSK]
       call sparse_axpy(ksk, denskern(is), -1.0_DP)
       k_ksk(is) = sparse_rms_element(ksk)
       k_ksk(is) = abs(k_ksk(is))

    end do

    ! ars: destroy matrices
    call sparse_destroy(sk)
    call sparse_destroy(ksk)
    call sparse_destroy(skh)
    call sparse_destroy(hks)

    ! ars: calculate 2tr[KS]
    do is = 1, pub_cell%num_spins
       ks_ne(is) = 2.0_DP*sparse_trace(denskern(is), overlap)/pub_cell%num_spins
    end do


  end subroutine kernel_diis_constants



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_energy(total_energy, ham, denskern, lhxc_fine, &
       ngwf_basis, hub_proj_basis, hub, rep, localpseudo_fine, &
       core_density_fine, ewald_energy, elements, ham_update, lhxc_fine_fixed)

    !===================================================================!
    ! This subroutine calculates the total energy after each DIIS step. !
    !-------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in June 2010.                       !
    !===================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: DP, max_spins, paw_en_size, VERBOSE
    use function_basis, only: FUNC_BASIS
    use hamiltonian, only: hamiltonian_dens_dep_matrices
    use ion, only: ELEMENT
    use hubbard_build, only: HUBBARD_MODEL
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_usehfx, pub_any_nl_proj, pub_paw, pub_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_trace, sparse_create, sparse_destroy, &
         sparse_copy, sparse_scale
    use vdwcorrection, only: pub_dispersion_energy


    ! ars: arguments
    real(kind=DP),    intent(  out) :: total_energy
    type(NGWF_HAM),   intent(inout) :: ham
    type(SPAM3),      intent(inout) :: denskern(pub_cell%num_spins)
    real(kind=DP),    intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    type(FUNC_BASIS), intent(in   ) :: ngwf_basis
    type(FUNC_BASIS), intent(in   ) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(NGWF_REP),   intent(in   ) :: rep
    real(kind=DP),    intent(in   ) :: localpseudo_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP),    intent(in   ) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    real(kind=DP),    intent(in   ) :: ewald_energy
    logical,          intent(in   ) :: ham_update
    logical,          intent(in   ) :: lhxc_fine_fixed
    type(ELEMENT),    intent(in   ) :: elements(pub_cell%nat)

    ! ars: local variables
    integer :: is, nspins
    real(kind=DP) :: spin_fac
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    real(kind=DP) :: lhxc_energy
    real(kind=DP) :: hubbard_energy


    ! ars: initialise
    nspins = pub_cell%num_spins
    spin_fac = 2.0_DP/nspins
    total_energy = 0.0_DP
    lhxc_energy = 0.0_DP
    hubbard_energy = 0.0_DP
    paw_sphere_energies(:) = 0.0_DP


    ! ars: scale kernel to cope with the number of spins
    do is=1,nspins
       call sparse_scale(denskern(is), spin_fac)
    end do

    ! ars: update density dependent matrices
    call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
         lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
         ngwf_basis, hub_proj_basis, hub, denskern, ewald_energy, elements, &
         pub_fine_grid, localpseudo_fine, core_density_fine, ham_update, &
         lhxc_fine_fixed)

    ! ars: scale kernel back
    do is=1,nspins
       call sparse_scale(denskern(is),1.0_DP/spin_fac)
    end do




  end subroutine kernel_diis_energy



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_residue_inout(residue,kern_out,kern_in)

    !====================================================================!
    ! This subroutine calculates the density kernel residue after each   !
    ! kernel DIIS iteration. R_i = K^out_i - K^in_i.                     !
    !--------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                      !
    !====================================================================!

    use comms, only: pub_on_root
    use constants, only: DP
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_copy, sparse_axpy, SPAM3


    implicit none

    ! ars: arguments
    type(SPAM3), intent(inout) :: residue(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: kern_out(pub_cell%num_spins)
    type(SPAM3), intent(in   ) :: kern_in(pub_cell%num_spins)

    ! ars: local variables
    integer :: is, nspins


    ! Establish number of spins
    nspins = pub_cell%num_spins

    do is = 1,nspins
       ! ars: R = K^out
       call sparse_copy(residue(is),kern_out(is))
       ! ars: calculate [R = K^out - K^in]
       call sparse_axpy(residue(is),kern_in(is),-1.0_DP)
    end do


  end subroutine kernel_diis_residue_inout



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_residual(residual,residue, overlap)

    !=====================================================================!
    ! This subroutine checks for convergence of the density kernel.       !
    ! Based on the latest residue matrix only.                            !
    !---------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano in April 2010.                       !
    !=====================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use rundat, only: pub_kernel_diis_thres
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_trace, sparse_create, sparse_destroy, &
         sparse_product, sparse_transpose
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    real(kind=DP), intent(  out) :: residual(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: residue(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: overlap

    integer :: nspins, is
    type(SPAM3) :: spam3buffer1, spam3buffer2, spam3buffer3


    ! Establish number of spins
    nspins = pub_cell%num_spins

    residual(:) = 0.0_DP

    ! ars: init spam3buffer blocking scheme
    call sparse_create(spam3buffer1, residue(1), overlap)
    call sparse_create(spam3buffer2, residue(1))
    call sparse_create(spam3buffer3, spam3buffer1, spam3buffer2)


    do is = 1, pub_cell%num_spins

       ! ars: spam3buffer1 = R x S
       call sparse_product(spam3buffer1, residue(is), overlap)

       ! ars: spam3buffer2 = (R)^t
       call sparse_transpose(spam3buffer2,residue(is))

       ! ars: spam3buffer3 = (R) x S x (R)^t
       call sparse_product(spam3buffer3, spam3buffer1, spam3buffer2)

       ! ars: Res(is) = tr[(R) x S x (R)^t x S]
       residual(is) = sparse_trace(spam3buffer3,overlap)

       ! ars: residual = sqrt(tr[RSRS])
       residual(is) = sqrt(residual(is))

    end do

    ! ars: destroy SPAM3 buffers
    call sparse_destroy(spam3buffer1)
    call sparse_destroy(spam3buffer2)
    call sparse_destroy(spam3buffer3)


  end subroutine kernel_diis_residual



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_conv_check(converged, okmessage, deltaE, residual, &
       hks_skh, k_ksk)

    use constants, only : DP
    use rundat, only: pub_kernel_diis_criteria, pub_kernel_diis_thres
    use simulation_cell, only: pub_cell

    implicit none

    logical,          intent(  out) :: converged
    character(len=8), intent(  out) :: okmessage
    real(kind=DP),    intent(in   ) :: deltaE
    real(kind=DP),    intent(in   ) :: residual(pub_cell%num_spins)
    real(kind=DP),    intent(in   ) :: hks_skh(pub_cell%num_spins)
    real(kind=DP),    intent(in   ) :: k_ksk(pub_cell%num_spins)

    integer :: is

    logical :: conv_criteria(1:4)

    converged = .false.
    conv_criteria(1:4) = .false.
    okmessage(1:8) = '........'


    ! ars: validate convergence criteria:

    ! residual
    if (pub_kernel_diis_criteria(1:1).eq.'1') then
       residual_loop: do is = 1, pub_cell%num_spins
          if (abs(residual(is)).lt.pub_kernel_diis_thres) then
             okmessage(1:2) = 'OK'
             conv_criteria(1) = .true.
          else
             okmessage(1:2) = 'NO'
             conv_criteria(1) = .false.
             exit residual_loop
          end if
       end do residual_loop
    else
       conv_criteria(1) = .true.
    end if

    ! [HKS,SKH]
    if (pub_kernel_diis_criteria(2:2).eq.'1') then
       hks_loop: do is = 1, pub_cell%num_spins
          if (abs(hks_skh(is)).lt.pub_kernel_diis_thres) then
             okmessage(3:4) = 'OK'
             conv_criteria(2) = .true.
          else
             okmessage(3:4) = 'NO'
             conv_criteria(2) = .false.
             exit hks_loop
          end if
       end do hks_loop
    else
       conv_criteria(2) = .true.
    end if

    ! [K,KSK]
    if (pub_kernel_diis_criteria(3:3).eq.'1') then
       ksk_loop: do is = 1, pub_cell%num_spins
          if (abs(k_ksk(is)).lt.pub_kernel_diis_thres) then
             okmessage(5:6) = 'OK'
             conv_criteria(3) = .true.
          else
             conv_criteria(3) = .false.
             okmessage(5:6) = 'NO'
             exit ksk_loop
          end if
       end do ksk_loop
    else
       conv_criteria(3) = .true.
    end if

    ! Energy variations
    if (pub_kernel_diis_criteria(4:4).eq.'1') then
       if(abs(deltaE).lt.pub_kernel_diis_thres) then
          okmessage(7:8) = 'OK'
          conv_criteria(4) = .true.
       else
          okmessage(7:8) = 'NO'
       end if
    else
       conv_criteria(4) = .true.
    end if

    converged = conv_criteria(1).and.conv_criteria(2).and.&
         conv_criteria(3).and.conv_criteria(4)


  end subroutine kernel_diis_conv_check



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_info_line(converged, okmessage, iter, residual, &
       total_energy, homolumo, mu, hks_skh, k_ksk, ks_ne, deltaE)

    !===================================================================!
    ! This subroutine prints convergence information at each DIIS iter. !
    !-------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in December 2010.                   !
    !===================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, NORMAL, VERBOSE, max_spins
    use rundat, only: pub_output_detail, pub_kernel_diis_type, &
         pub_kernel_diis_liter, pub_kernel_diis_criteria
    use simulation_cell, only: pub_cell


    implicit none

    logical,          intent(in   ) :: converged
    character(len=8), intent(in   ) :: okmessage
    integer,          intent(in   ) :: iter
    real(kind=DP),    intent(in   ) :: residual(pub_cell%num_spins)
    real(kind=DP),    intent(in   ) :: total_energy
    real(kind=DP),    intent(in   ) :: homolumo(2,pub_cell%num_spins)
    real(kind=DP),    intent(in   ) :: mu(max_spins)
    real(kind=DP),    intent(in   ) :: hks_skh(pub_cell%num_spins)
    real(kind=DP),    intent(in   ) :: k_ksk(pub_cell%num_spins)
    real(kind=DP),    intent(in   ) :: ks_ne(pub_cell%num_spins)
    real(kind=DP),    intent(in   ) :: deltaE


    integer :: is
    character(len=1) :: diis_type
    character(len=4) :: is_label

    ! ars: print header at each iteration if VERBOSE
    if(pub_output_detail.ge.VERBOSE.and.pub_on_root) write(stdout,'(a80)') &
         "T|Iter|Spin| [HKS,SKH]|   [K,KSK]|  Residual|&
         &         Ne|            Energy(Ha)|"


    ! ars: work out which kind of DIIS has been performed
    if (pub_kernel_diis_type.eq.'D'.or.iter.eq.1) then
       diis_type = 'D'
    elseif ((pub_kernel_diis_type.eq.'L').or.&
         (iter.le.pub_kernel_diis_liter+1)) then
       diis_type = 'L'
    else
       diis_type = 'P'
    end if

    if(pub_cell%num_spins.eq.1) is_label = '   1'
    if(pub_cell%num_spins.gt.1) is_label = ' 1-2'

    ! ars: print type, iter, spin, commutators, residual, Ne and energy.
    !      print HOMO, LUMO and mu if verbose
    if (pub_on_root) then

       if(pub_output_detail.le.NORMAL) then
          write(stdout,'(a1,1x,i4,1x,a4,3(1x,es10.3e2),1x,f11.3,1x,f22.14)') &
               diis_type, iter, is_label, maxval(hks_skh(:)), maxval(k_ksk(:)), &
               maxval(residual(:)), sum(ks_ne(:)), total_energy

       else

          do is = 1 , pub_cell%num_spins

             write(stdout,'(a1,1x,i4,1x,i4,3(1x,es10.3e2),1x,f11.3,1x,f22.14)') &
                  diis_type, iter, is, hks_skh(is), k_ksk(is), residual(is), &
                  ks_ne(is), total_energy

             write(stdout,'(a5,i1,a4,f22.14)') "HOMO(", is, ") = ", homolumo(1,is)
             write(stdout,'(a5,i1,a4,f22.14)') "LUMO(", is, ") = ", homolumo(2,is)
             write(stdout,'(a3,i1,a6,f22.14)') "mu(", is, ")   = ", mu(is)

          end do

       end if

       ! ars: print banner if converged
       if(converged) then
          write(stdout,'(/,a27, i3, a26, f22.14)') &
               "==> Kernel converged after ", iter, &
               " iterations with energy = ", total_energy
          if (pub_output_detail.le.NORMAL) then
             do is = 1, pub_cell%num_spins
                if(pub_kernel_diis_criteria(1:1).eq.'1') &
                     write(stdout,'(a14,i1, a4, f22.14, 3x, a2)') &
                     "==> Residual (",is,") = ", residual(is), okmessage(1:2)
                if(pub_kernel_diis_criteria(2:2).eq.'1') &
                     write(stdout,'(a14,i1, a4, f22.14, 3x, a2)') &
                     "==> [KHS,SKH](",is,") = ", hks_skh(is), okmessage(3:4)
                if(pub_kernel_diis_criteria(3:3).eq.'1') &
                     write(stdout,'(a14,i1, a4, f22.14, 3x, a2)') &
                     "==> [K,KSK]  (",is,") = ", k_ksk(is), okmessage(5:6)
             end do
             if(pub_kernel_diis_criteria(4:4).eq.'1') &
                  write(stdout,'(a19, f22.14, 3x, a2,/)') &
                  "==> deltaE       = ", deltaE, okmessage(7:8)
          end if
       end if

       ! ars: print convergence info if verbose
       if (pub_output_detail.ge.VERBOSE) then
          do is = 1, pub_cell%num_spins
             write(stdout,'(a14,i1, a4, f22.14, 3x, a2)') &
                  "==> Residual (",is,") = ", residual(is), okmessage(1:2)
             write(stdout,'(a14,i1, a4, f22.14, 3x, a2)') &
                  "==> [KHS,SKH](",is,") = ", hks_skh(is), okmessage(3:4)
             write(stdout,'(a14,i1, a4, f22.14, 3x, a2)') &
                  "==> [K,KSK]  (",is,") = ", k_ksk(is), okmessage(5:6)
          end do
          write(stdout,'(a19, f22.14, 3x, a2,/)') &
               "==> deltaE       = ", deltaE, okmessage(7:8)
       end if

    end if


  end subroutine kernel_diis_info_line

  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_shift(shifted, dkn_in, dkn_out, residues, iter)

    !================================================================!
    ! This subroutine shifts the matrices in arrays in case iter has !
    ! reached pub_kernel_diis_max value.                             !
    !----------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in May 2010.                     !
    !================================================================!

    use constants, only: DP, stdout
    use rundat, only: pub_kernel_diis_max
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_copy

    implicit none

    ! ars: arguments
    logical,     intent(  out) :: shifted
    type(SPAM3), intent(inout) :: dkn_in(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: dkn_out(pub_cell%num_spins,pub_kernel_diis_max)
    type(SPAM3), intent(inout) :: residues(pub_cell%num_spins,pub_kernel_diis_max)
    integer,     intent(in   ) :: iter

    ! ars: local variables
    integer :: scntr, is

    shifted = .false.

    ! ars: move old entries one position back
    if (iter.ge.pub_kernel_diis_max) then
       shifted = .true.
       do scntr = 1, pub_kernel_diis_max-1
          do is = 1, pub_cell%num_spins
             call sparse_copy(dkn_in(is,scntr), dkn_in(is, scntr+1))
             call sparse_copy(dkn_out(is,scntr), dkn_out(is, scntr+1))
             call sparse_copy(residues(is,scntr), residues(is, scntr+1))
          end do
       end do
    end if

  end subroutine kernel_diis_shift



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_find_ientry(ientry, iter, shifted)

    !==============================================================!
    ! This subroutine finds the next ientry to work with.           !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in May 2010.                   !
    !==============================================================!

    use comms, only: comms_barrier
    use rundat, only: pub_kernel_diis_max

    implicit none

    ! ars: arguments
    integer, intent(  out) :: ientry
    integer, intent(in   ) :: iter
    logical, intent(in   ) :: shifted


    ! ars: find next position
    if (shifted) then
       ientry = pub_kernel_diis_max
    else
       ientry = iter+1
    end if

    call comms_barrier


  end subroutine kernel_diis_find_ientry



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_print_header()

    !==============================================================!
    ! This subroutine prints the header at each kernel-DIIS iter.  !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in January 2011.               !
    !==============================================================!

    use comms, only: pub_on_root
    use constants, only: stdout, NORMAL
    use rundat, only: pub_output_detail

    implicit none


    ! ars: iteration header
    if (pub_on_root) then
       write(stdout,'(/,a80)') "------------------ &
            &Hamiltonian diagonalisation + kernel DIIS -------------------"
       if(pub_output_detail.le.NORMAL) write(stdout,'(a80)') &
            "T|Iter|Spin| [HKS,SKH]|   [K,KSK]|  Residual|&
            &         Ne|            Energy(Ha)|"
    end if


  end subroutine kernel_diis_print_header



  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_check_params

    !========================================================================!
    ! This subroutine checks the correctness of the kernel-DIIS parameters.  !
    !------------------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in December 2010.                        !
    !========================================================================!

    use comms, only : pub_on_root, comms_barrier, comms_abort
    use constants, only: VERBOSE, stdout, DP
    use rundat, only: pub_kernel_diis_maxit, pub_kernel_diis_max, &
         pub_kernel_diis_liter, pub_kernel_diis_c_in, pub_kernel_diis_c_out, &
         pub_kernel_diis_criteria, pub_kernel_diis_type, pub_kernel_diis_thres,&
         pub_output_detail

    implicit none


    logical, save :: check_success = .false.


    if (.not.check_success) then

       if (pub_on_root.and.pub_output_detail.ge.VERBOSE) &
            write(stdout,'(/,a35)', advance='no') "Checking kernel-DIIS parameters"


       ! ars: check pub_kernel_diis_maxit and correct if negative
       if (pub_kernel_diis_maxit.le.0) then

          if(pub_on_root) then
             write(stdout,'(/,a)') &
                  "*** WARNING: Invalid kernel-DIIS max iterations ****"
             write(stdout,'(a)')    "          kernel_diis_maxit <= 0"
             write(stdout,'(a,/)')  "             Program stops"
          end if

          call comms_abort

       end if

       ! ars: check pub_kernel_diis_max and correct if negative
       if (pub_kernel_diis_max.le.0) then

          if(pub_on_root) then
             write(stdout,'(/,a)') &
                  "*** WARNING: Invalid kernel-DIIS max ***************"
             write(stdout,'(a)')    "             kernel_diis_max <= 0"
             write(stdout,'(a,i3)') &
                  "             kernel_diis_max = ", pub_kernel_diis_max
             write(stdout,'(a,/)') &
                  "             Reseting to safe value kernel_diis_max = 5"
          end if

          pub_kernel_diis_max = 5

       end if

       ! ars: check pub_kernel_diis_max is not greater than pub_kernel_diis_maxit
       if (pub_kernel_diis_max.gt.pub_kernel_diis_maxit) then

          if(pub_on_root) then
             write(stdout,'(/,a)') &
                  "*** WARNING: Invalid kernel-DIIS max ***************"
             write(stdout,'(a)')    "             kernel_diis_max > kernel_diis_maxit"
             write(stdout,'(a,i3)')&
                  "             kernel_diis_max = ", pub_kernel_diis_max
          end if

          if (pub_kernel_diis_maxit.le.5) then
             pub_kernel_diis_max = pub_kernel_diis_maxit
          else
             pub_kernel_diis_max = 5
          end if

          if(pub_on_root) write(stdout,'(a,i1,/)') &
               "             Reseting to safe value kernel_diis_max = ", &
               pub_kernel_diis_max
       end if

       ! ars: check pub_kernel_diis_liter is not greater than pub_kernel_diis_maxit
       if (pub_kernel_diis_liter.gt.pub_kernel_diis_maxit) then

          if(pub_on_root) then
             write(stdout,'(/,a)')  &
                  "*** WARNING: Invalid kernel-DIIS linear iterations ******"
             write(stdout,'(a,i3)') &
                  "             kernel_diis_liter = ", pub_kernel_diis_liter
          end if

          if (pub_kernel_diis_maxit.le.3) then
             pub_kernel_diis_liter = pub_kernel_diis_maxit
          else
             pub_kernel_diis_liter = 3
          end if


          if(pub_on_root) write(stdout,'(a,i1,/)') &
               "             Reseting to safe value kernel_diis_liter = ",&
               pub_kernel_diis_liter
       end if


       ! ars: check C_in and C_out for linear mixing
       if((pub_kernel_diis_c_in.lt.0.0_DP).or.(pub_kernel_diis_c_in.gt.1.0_DP)) then

          if(pub_on_root) then
             write(stdout,'(/,a)')     &
                  "*** WARNING: Invalid kernel-DIIS C_in coefficient ******"
             write(stdout,'(a,f22.12)')&
                  "             kernel_diis_c_in = ", pub_kernel_diis_c_in
             write(stdout,'(a,/)')     &
                  "             Reseting to safe value kernel_diis_c_in=0.8000"
          end if

          pub_kernel_diis_c_in = 0.8_DP
          pub_kernel_diis_c_out = 1.0_DP - pub_kernel_diis_c_in

       end if

       ! ars: check kernel_diis_criteria
       if(  ((pub_kernel_diis_criteria(1:1).ne.'0').and.&
            (pub_kernel_diis_criteria(1:1).ne.'1')).or.&
            ((pub_kernel_diis_criteria(2:2).ne.'0').and.&
            (pub_kernel_diis_criteria(2:2).ne.'1')).or.&
            ((pub_kernel_diis_criteria(3:3).ne.'0').and.&
            (pub_kernel_diis_criteria(3:3).ne.'1')).or.&
            ((pub_kernel_diis_criteria(4:4).ne.'0').and.&
            (pub_kernel_diis_criteria(4:4).ne.'1')).or.&
            (pub_kernel_diis_criteria(1:4).eq.'0000'))then

          if(pub_on_root) then
             write(stdout,'(/,a)') &
                  "*** WARNING: Invalid kernel-DIIS convergence criteria ******"
             write(stdout,'(a,a4)')&
                  "             kernel_diis_criteria = ", pub_kernel_diis_criteria
             write(stdout,'(a,/)') &
                  "             Reseting to safe value kernel_diis_criteria='1000'"

          end if
          pub_kernel_diis_criteria(1:4) = '1000'
       end if

       ! ars: check DIIS type
       if(  (pub_kernel_diis_type(1:1).ne.'P').and.&
            (pub_kernel_diis_type(1:1).ne.'p').and.&
            (pub_kernel_diis_type(1:1).ne.'L').and.&
            (pub_kernel_diis_type(1:1).ne.'l').and.&
            (pub_kernel_diis_type(1:1).ne.'D').and.&
            (pub_kernel_diis_type(1:1).ne.'d')) then

          if(pub_on_root) then
             write(stdout,'(/,a)') &
                  "*** WARNING: Invalid kernel-DIIS method *************"
             write(stdout,'(a,a1)')&
                  "             kernel_diis_type = ", pub_kernel_diis_type
             write(stdout,'(a,/)') &
                  "             Reseting to safe value kernel_diis_type='P'"
          end if

          pub_kernel_diis_type(1:1) = 'P'

       end if

       ! ars: check DIIS threshold
       if(pub_kernel_diis_thres.le.0.0_DP) then
          if(pub_on_root) then
             write(stdout,'(/,a)')     &
                  "*** WARNING: Invalid kernel-DIIS theshold **************"
             write(stdout,'(a,f22.12)')&
                  "             kernel_diis_thres = ", pub_kernel_diis_thres
             write(stdout,'(a,/)')     &
                  "             Reseting to safe value kernel_diis_thres=1E-09"
          end if

          pub_kernel_diis_thres = 0.000000001_DP

       end if

       check_success = .true.
       if (pub_on_root.and.pub_output_detail.ge.VERBOSE) &
            write(stdout,'(a,/)') " ... done"

       call comms_barrier

    end if

  end subroutine kernel_diis_check_params

  !_____________________________________________________________________________
  !_____________________________________________________________________________



  subroutine kernel_diis_qc_output(total_energy, &
       deltaE, hks_skh, k_ksk, residual, mu, denskern, pur_denskern, ham)

    !==============================================================!
    ! This subroutine prints <QC> lines related to kernel-DIIS.    !
    !--------------------------------------------------------------!
    ! Wrtten by Alvaro Ruiz Serrano in December 2010.              !
    !==============================================================!

    use constants, only: stdout, DP, max_spins
    use comms, only: comms_barrier, pub_on_root
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_max_abs_element

    implicit none

    real(kind=DP), intent(in   ) :: total_energy
    real(kind=DP), intent(in   ) :: deltaE
    real(kind=DP), intent(in   ) :: hks_skh(pub_cell%num_spins)
    real(kind=DP), intent(in   ) :: k_ksk(pub_cell%num_spins)
    real(kind=DP), intent(in   ) :: residual(pub_cell%num_spins)
    real(kind=DP), intent(in   ) :: mu(max_spins)
    type(SPAM3),   intent(in   ) :: denskern(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: pur_denskern(pub_cell%num_spins)
    type(SPAM3),   intent(in   ) :: ham(pub_cell%num_spins)

    integer :: is

    real(kind=DP) :: max_denskern(pub_cell%num_spins)
    real(kind=DP) :: max_pur_denskern(pub_cell%num_spins)
    real(kind=DP) :: max_ham(pub_cell%num_spins)

    do is = 1 , pub_cell%num_spins

       max_denskern(is) = sparse_max_abs_element(denskern(is))
       max_pur_denskern(is) = sparse_max_abs_element(pur_denskern(is))
       max_ham(is) = sparse_max_abs_element(ham(is))

    end do

    if(pub_on_root) then

       !write(stdout, '(/,a30, i9)')   '<QC>             [array_size]:', &
       !     ientry

       !write(stdout, '(a30, i9)')     '<QC> [kernel_diis_iterations]:', &
       !     iter

       write(stdout, '(a30, f22.12)') '<QC>           [total_energy]:', &
            total_energy

       write(stdout, '(a30, f22.12)') '<QC>                 [deltaE]:', &
            deltaE


       do is = 1, pub_cell%num_spins

          write(stdout, '(a26, i1, a3, f22.12)') &
               '<QC>             [HKS,SKH(', is, ')]:', hks_skh(is)

          write(stdout, '(a26, i1, a3, f22.12)') &
               '<QC>               [K,KSK(', is, ')]:', k_ksk(is)

          write(stdout, '(a26, i1, a3, f22.12)') &
               '<QC>            [residual(', is, ')]:', residual(is)

          write(stdout, '(a26, i1, a3, f22.12)') &
               '<QC>                  [mu(', is, ')]:', mu(is)

          write(stdout, '(a26, i1, a3, f22.12)') &
               '<QC>        [max(denskern(', is, ')]:', max_denskern(is)

          write(stdout, '(a26, i1, a3, f22.12)') &
               '<QC>    [max(pur_denskern(', is, ')]:', max_pur_denskern(is)

          write(stdout, '(a26, i1, a3, f22.12)') &
               '<QC>             [max(ham(', is, ')]:', max_ham(is)

       end do


    end if

    call comms_barrier

  end subroutine kernel_diis_qc_output



  !_____________________________________________________________________________
  !_____________________________________________________________________________



end module kernel_diis
