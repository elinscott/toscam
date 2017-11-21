! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!         Time-dependent Density Functional Theory module        !
!                                                                !
! This module performs calculations of dielectric                !
! properties using linear-response TDDFT in real time.           !
!----------------------------------------------------------------!
! Written by David D. O'Regan in February 2009.                  !
!================================================================!

module tddft

  implicit none

  private

  public :: tddft_calculate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tddft_calculate( &
       dft_denskern, rep, ngwf_basis, proj_basis, elements, &
       localpseudo_fine, core_density_fine, ewald_energy)

    !======================================================================!
    ! This subroutine propagates the density-matrix in real time subject   !
    ! to an applied perturbation in TDDFT.                                 !
    !----------------------------------------------------------------------!
    ! Written by David D. O'Regan in February 2009.                        !
    !======================================================================!

    use comms, only: comms_abort, pub_on_root
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_locpot
    use ion, only: ELEMENT
    use constants, only: DP, PI, stdout
    use ngwf_representation, only: NGWF_REP
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_product,&
         sparse_destroy, sparse_copy, sparse_scale, &
         sparse_convert, sparse_trace, sparse_rms_element
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
    use rundat, only: pub_tddft_tammdancoff, pub_tddft_hamiltonian_mixing, &
         pub_tddft_maximum_energy, pub_tddft_resolution, &
         pub_tddft_inv_overlap_exact, maxit_hotelling, pub_tddft_damping, &
         pub_tddft_enforced_idempotency, pub_tddft_propagation_method, &
         pub_tddft_dipole_kick_strength, pub_any_nl_proj, pub_paw
    use properties, only: properties_polarisation
    use wrappers, only: wrappers_1d_fft, wrappers_invert_sym_matrix

    implicit none

    ! Unperturbed denskern matrix
    type(SPAM3), intent(inout)   :: dft_denskern(pub_cell%num_spins)
    type(NGWF_REP), intent(inout):: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis  ! function basis for NGWFs
    type(FUNC_BASIS), intent(in) :: proj_basis  ! function basis for Projectors
    ! elements of all atoms in  input file order
    type(ELEMENT), intent(in)    :: elements(:)
    ! local pseudo on this proc slab
    real(kind=DP), intent(in)    :: localpseudo_fine(:,:,:)
    ! core density on this proc slab
    real(kind=DP), intent(in)    :: core_density_fine(:,:,:)
    real(kind=DP), intent(in)    :: ewald_energy

    !<<local variables>>
    type(SPAM3) :: knl           ! kinetic plus nonlocal in SPAM3 format
    type(SPAM3), allocatable :: current_ham(:)! hamiltonian in SPAM3 format
    type(SPAM3), allocatable :: midpoint_ham(:)
    type(SPAM3), allocatable :: current_ham_1(:)
    type(SPAM3), allocatable :: current_ham_2(:)
    type(SPAM3), allocatable :: dft_ham(:)
    type(SPAM3), allocatable :: perturbation(:)
    type(SPAM3), allocatable :: re_denskern(:)
    type(SPAM3), allocatable :: im_denskern(:)
    type(SPAM3), allocatable :: re_del_denskern(:)
    type(SPAM3), allocatable :: im_del_denskern(:)
    type(SPAM3), allocatable :: re_del_denskern_buf(:)
    type(SPAM3), allocatable :: im_del_denskern_buf(:)
    type(SPAM3), allocatable :: re_right(:), im_right(:) ! propagators with H_0
    type(SPAM3), allocatable :: re_left(:), im_left(:)   ! propagators with H_0
    integer :: ierr                    ! memory allocation error flag
    integer :: is                      ! spin loop counter
    ! Cartesian axes along which to apply electric dipole perturbation
    integer :: cartesian
    integer :: timestep                ! The number timestep we are on
    integer :: num_timesteps
    real(kind=DP) :: initial_energy, energy
    real(kind=DP) :: steplength, duration, damping_exponent, phasekick, phasekck
    ! Dipole moment on time axis
    complex(kind=DP), allocatable :: dipole_t(:,:)
    ! Dipole moment on frequency axis
    complex(kind=DP), allocatable :: dipole_w(:,:)
    ! Dipole moment of unperturbed density
    complex(kind=DP), allocatable :: dipole_t_initial(:)
    real(kind=DP), allocatable :: inverse_s_square(:,:)
    real(kind=DP) :: real_output_dipole, momentum
    ! ddor: Internal variables used for self-consistently
    ! enforcing time-reversal symmetry
    integer :: maxit_time_reversal, sc_time_reversal_iter
    !real(kind=DP):: sc_time_reversal_rms, new_sc_time_reversal_rms, &
    real(kind=DP) :: sc_time_reversal_energy, new_sc_time_reversal_energy

    ! ddor: spectrum oversampling
    integer, parameter :: spectrum_finesse_factor = 1

    ! ddor: propagation_scheme_new = .false.
    ! P (t + \Delta t) = exp (- i (H0 + H1) \Delta t ) (P0 + P1)
    ! exp (+ i (H0 + H1) \Delta t )
    ! has [H1 , P1] term in its time derivative

    ! ddor: propagation_scheme_new = .true. given by
    ! P (t + \Delta t) = exp (- i H0 \Delta t ) (P1)  exp (+ i H0 \Delta t ) +
    ! exp (- i H1 \Delta t ) (P0)  exp (+ i H1 \Delta t )
    ! which should respect idempotency to first order so long
    ! as the perturbation does not dirty P1(t = 0+) excessively.
    ! In this case  exp (- i H0 \Delta t ), c.c.,
    ! only needs to be calculated once initially and the propagation subroutines
    ! do not have to deal with an imaginary part of the density kernel

    logical :: propagation_scheme_new = .true.
    if (pub_tddft_tammdancoff) propagation_scheme_new = .false.
    ! ddor: Allow maxit_time_reversal self-consistency iterations
    !       of time-reversal symmetry enforcement scheme
    maxit_time_reversal = 5

    ! cks: start timer
    call timer_clock("tddft_calculate",1)

    if (pub_on_root) write(stdout,'(/a)')'#################################&
         &#######&
         &########################################'
    if (pub_on_root) write(stdout,*) 'WARNING, WARNING: ENTERING EXPERIMENTAL &
         &TDDFT MODULE. THIS FUNCTIONALITY IS NOT YET STABLE!!'
    if (pub_on_root) write(stdout,*) 'WARNING, WARNING: ENTERING EXPERIMENTAL &
         &TDDFT MODULE. THIS FUNCTIONALITY IS NOT YET STABLE!!'
    if (pub_on_root) write(stdout,'(/a)')'#################################&
         &#######&
         &########################################'

    allocate(current_ham(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('tddft_calculate','current_ham',ierr)
    allocate(midpoint_ham(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('tddft_calculate','midpoint_ham',ierr)
    allocate(perturbation(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('tddft_calculate','perturbation',ierr)
    allocate(re_denskern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('tddft_calculate','re_denskern',ierr)
    allocate(im_denskern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('tddft_calculate','im_denskern',ierr)
    if (propagation_scheme_new) then
       allocate(re_del_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','re_del_denskern',ierr)
       allocate(im_del_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','im_del_denskern',ierr)
       allocate(re_right(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','re_right',ierr)
       allocate(im_right(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','im_right',ierr)
       allocate(re_left(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','re_left',ierr)
       allocate(im_left(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','im_left',ierr)
       allocate(dft_ham(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','dft_ham',ierr)
    endif
    if ((pub_tddft_hamiltonian_mixing == 1) .or. &
         &(pub_tddft_hamiltonian_mixing == 2)) then
       allocate(current_ham_1(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','current_ham_1',ierr)
    endif
    if (pub_tddft_hamiltonian_mixing == 2) then
       allocate(current_ham_2(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','current_ham_2',ierr)
    endif
    if (pub_tddft_propagation_method == 'RK_SC_TIME_REVERSAL') then
       allocate(re_del_denskern_buf(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','re_del_denskern_buf',ierr)
       allocate(im_del_denskern_buf(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('tddft_calculate','im_del_denskern_buf',ierr)
    endif

    knl%structure = 'H'
    re_denskern%structure ='K'
    call sparse_create(knl)
    do is=1,pub_cell%num_spins
       call sparse_create(current_ham(is),knl)
       call sparse_create(midpoint_ham(is),knl)
       call sparse_create(re_denskern(is))
       call sparse_create(im_denskern(is),re_denskern(1))
       call sparse_create(perturbation(is),re_denskern(1))
       if (propagation_scheme_new) then
          call sparse_create(re_del_denskern(is),re_denskern(1))
          call sparse_create(im_del_denskern(is),im_denskern(1))
          call sparse_create(dft_ham(is),knl)
       endif
       if ((pub_tddft_hamiltonian_mixing == 1) .or. &
            &(pub_tddft_hamiltonian_mixing == 2)) then
          call sparse_create(current_ham_1(is),knl)
       endif
       if (pub_tddft_hamiltonian_mixing == 2) then
          call sparse_create(current_ham_2(is),knl)
       endif
       if (pub_tddft_propagation_method == 'RK_SC_TIME_REVERSAL') then
          call sparse_create(re_del_denskern_buf(is),re_denskern(1))
          call sparse_create(im_del_denskern_buf(is),im_denskern(1))
       endif
    end do

    call sparse_copy(knl,rep%kinet)
    if (pub_any_nl_proj) then
       call sparse_axpy(knl,rep%nonlocpot(1),1.0_DP)
    end if
    if (pub_paw) then
       call utils_abort('Error in tddft_calculate: PAW TDDFT calculations not &
            &yet supported')
    end if

    num_timesteps = CEILING( pub_tddft_maximum_energy / pub_tddft_resolution )

    ! Allocate work arrays for dipole moment vectors
    allocate(dipole_t(spectrum_finesse_factor*2*num_timesteps,3),stat=ierr)
    call utils_alloc_check('tddft_calculate','dipole_t',ierr)
    allocate(dipole_w(spectrum_finesse_factor*2*num_timesteps,3),stat=ierr)
    call utils_alloc_check('tddft_calculate','dipole_w',ierr)
    allocate(dipole_t_initial(3),stat=ierr)
    call utils_alloc_check('tddft_calculate','dipole_t_initial',ierr)

    dipole_t = (0.0_DP,0.0_DP)
    dipole_w = (0.0_DP,0.0_DP)
    dipole_t_initial = (0.0_DP,0.0_DP)

    duration = 1.0_DP  / pub_tddft_resolution
    steplength = 1.0_DP / pub_tddft_maximum_energy
    ! ddor: pub_tddft_damping is damping factor at end of propagation
    damping_exponent = -1.0_DP * LOG(pub_tddft_damping) / duration
    !       thekick =  (0.5_DP)* &
    !            SQRT( (pub_tddft_dipole_kick_strength(1)**(2.0_DP)) + &
    !                  (pub_tddft_dipole_kick_strength(2)**(2.0_DP)) + &
    !                  (pub_tddft_dipole_kick_strength(3)**(2.0_DP)) )

    ! 2.418 x 10^-2 fs = 1 a.t.u.
    if (pub_on_root) write(stdout,*) 'tddft_calculate: &
         &NUM_TIMESTEPS: ', num_timesteps
    if (pub_on_root) write(stdout,*) 'tddft_calculate: &
         &PROPAGATION DURATION : ', duration* 2.418884326505E-02 ,' fs'
    if (pub_on_root) write(stdout,*) 'tddft_calculate: &
         &PROPAGATION TIME-STEP LENGTH : ', &
         &steplength* 2.418884326505E-02 ,' fs'

!!!! RENEW THE INVERSE OVERLAP PRECISELY !!!!!!
    ! If maxit_hotelling == 0 then this has already been done.
    if ( (pub_tddft_inv_overlap_exact) .and. (maxit_hotelling .gt. 0) )then

       ! cks: calculate exact inverse overlap
       allocate(inverse_s_square(ngwf_basis%num, ngwf_basis%num),stat=ierr)
       call utils_alloc_check('tddft_calculate', &
            'inverse_s_square',ierr)

       ! cks: expand the overlap matrix
       call sparse_convert(inverse_s_square, rep%overlap)

       ! cks: invert overlap matrix
       call wrappers_invert_sym_matrix(inverse_s_square, ngwf_basis%num)

       ! cks: put inv_overlap back to the same
       ! cks: sparsity pattern as the density kernel
       call sparse_convert(rep%inv_overlap, inverse_s_square)

       deallocate(inverse_s_square,stat=ierr)
       call utils_dealloc_check('tddft_calculate', &
            'inverse_s_square',ierr)

    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (pub_on_root) write(stdout,*) 'tddft_idempotency_test'
    call tddft_idempotency_test(0, dft_denskern, rep%overlap)
    ! OUTPUT THE OCCUPANCY AND IDEMPOTENCY ERROR
    !if (pub_on_root) write(stdout,*) 'tddft_idempotency'
    !         call tddft_idempotency(0, dft_denskern, overlap)

    if (pub_on_root) write(stdout,*) '5',&
         sparse_trace(dft_denskern(1),rep%overlap)

    cart:do cartesian=1,3

       ! ddor: Negative kick causes this direction to be omitted
       if (pub_tddft_dipole_kick_strength(cartesian) .lt. 0.0_DP ) cycle cart

       if (pub_on_root) write(stdout,*) 'tddft_calculate: &
            &PERTURBATION IN CARTESIAN DIRECTION: ', cartesian

       ! Apply same perturbation to both spin channels for the time being
       ! This does not alter the energy to first order
       call tddft_perturbation_initiate(perturbation, momentum, &
            dft_denskern, rep%inv_overlap, rep%overlap, cartesian,&
            ngwf_basis, rep%ngwfs_on_grid)

       if (pub_on_root) write(stdout,*) '6',&
            sparse_trace(dft_denskern(1),rep%overlap)

       ! OUTPUT THE OCCUPANCY AND IDEMPOTENCY ERROR
       if (pub_tddft_enforced_idempotency) then
          call tddft_idempotency(0, dft_denskern, rep%overlap, &
               perturbation, dft_denskern,.false.)
       else
          call tddft_idempotency(0, dft_denskern, rep%overlap, perturbation)
       endif

       if (pub_on_root) write(stdout,*) '7',&
            sparse_trace(dft_denskern(1),rep%overlap)

       ! Calculate the Kohn-Sham energy and Hamiltonian
       call tddft_electronic_energy(initial_energy, current_ham,&
            rep%ngwfs_on_grid, dft_denskern,&
            rep%overlap, knl, localpseudo_fine, core_density_fine, ngwf_basis, &
            ewald_energy)

       if (pub_on_root) write(stdout,*) 'tddft_calculated: &
            &ENERGY ON TIMESTEP 0 IS ', initial_energy

       ! Store the Kohn-Sham Hamiltonian
       if (propagation_scheme_new) then
          do is=1,pub_cell%num_spins
             call sparse_copy(dft_ham(is),current_ham(is))
          enddo
       endif

       ! Put initial denskern into working copy
       do is=1,pub_cell%num_spins
          call sparse_copy(re_denskern(is),dft_denskern(is))
          call sparse_copy(im_denskern(is),perturbation(is))
       enddo

       ! POLARISATION
       if (pub_on_root) write(stdout,*) 'tddft_calculate: &
            &DIPOLE MOMENT ON TIMESTEP 0'
       call properties_polarisation(rep,ngwf_basis,proj_basis, &
            elements,dft_denskern,cartesian,real_output_dipole)
       ! END POLARISATION
       dipole_t_initial(cartesian) = &
            CMPLX(real_output_dipole,0.0_DP,kind=DP)

       phasekick = 0.0_DP
       do is=1,pub_cell%num_spins
          phasekick = phasekick + sparse_trace(im_denskern(is),rep%overlap)
       end do
       if (pub_on_root) write(stdout,*) 'tddft_calculate: &
            &Tr[Im[\rho]]', phasekick
       phasekick = ATAN(phasekick /  REAL(SUM(rep%n_occ(:)))  )
       if (pub_on_root) write(stdout,*) 'tddft_calculate: &
            &ARCTAN ( Tr[Im[\rho]] / Tr[Re[\rho]])', phasekick
       phasekck = 0.0_DP
       phasekick = 0.0_DP
       do is=1,pub_cell%num_spins
          phasekck = phasekck + sparse_rms_element(re_denskern(is))
       enddo
       do is=1,pub_cell%num_spins
          phasekick = phasekick + sparse_rms_element(im_denskern(is))
       enddo
       phasekick = ATAN(phasekick / phasekck)
       if (pub_on_root) write(stdout,*) 'tddft_calculate: &
            &ARCTAN ( RMS[Im[\rho]] / RMS[Re[\rho]])', phasekick

       ! Calculate  exp (+/- i H0 \Delta t ) once and for all
       if (propagation_scheme_new) then
          if (pub_on_root) write(stdout,*) 'Calculating propagator with &
               &Kohn-Sham Hamiltonian'
          cn4: if (pub_tddft_propagation_method == 'CRANKNICHOLSON') then
             call tddft_crank_nicholson_prop_init(&
                  rep%inv_overlap, dft_ham, steplength, &
                  re_right, im_right, re_left, im_left)
          else cn4
             call tddft_runge_kutta_prop_init(&
                  rep%inv_overlap, dft_ham, steplength, &
                  re_right, im_right, re_left, im_left)
          endif cn4
          if (pub_on_root) write(stdout,*) 'sparse_trace(x) where &
               & x = exp (+/- i H0 S^-1 \Delta t ), &
               &exp (+/- i S^-1 H0 \Delta t ) for spin alpha',&
               sparse_trace(re_right(1)),sparse_trace(im_right(1)), &
               sparse_trace(re_left(1)),sparse_trace(im_left(1))
       endif

       interval:do timestep=1,num_timesteps

          ! Mix in some of the previous Hamiltonian to extrapolate to
          ! U(T + 0.5 \Delta T, T) = exp (- i H(T + 0.5 \Delta T) \Delta T )
          ! in order to approximately preserve time-reversal symmetry.
          ! For TDDFT_HAMILTONIAN MIXING == 2 we have
          ! H(T + 0.5 \Delta T) =
          ! 0.125 [3H(T-2\Delta T)-10H(T-\Delta T )+15H(T)]
          ! For TDDFT_HAMILTONIAN MIXING == 1 we have
          ! H(T + 0.5 \Delta T) = 0.5 [-H(T-\Delta T )+3H(T)]

          ! If using enforced time-reversal symmetry both H(T)
          ! and an approximation to H(T + \Delta T) are needed
          ! Mix in some of the previous Hamiltonian to extrapolate to
          ! U(T + \Delta T, T) = exp (- i H(T + \Delta T) \Delta T / 2 )
          !                      exp (- i H(T) \Delta T / 2 )
          ! in order to approximately preserve time-reversal symmetry.
          ! For TDDFT_HAMILTONIAN MIXING == 2 we have
          ! H(T + \Delta T) =
          ! [ H(T-2\Delta T)-3H(T-\Delta T )+3H(T)]
          ! For TDDFT_HAMILTONIAN MIXING == 1 we have
          ! H(T +\Delta T) = [-H(T-\Delta T )+2H(T)]

          mix: if (pub_tddft_hamiltonian_mixing == 2) then
             if ((pub_tddft_propagation_method == 'RK_TIME_REVERSAL') .or. &
                  &(pub_tddft_propagation_method == 'RK_SC_TIME_REVERSAL')) then
                if (pub_on_root) write(stdout,*) &
                     &'pub_tddft_hamiltonian_mixing == 2 and so&
                     & H(T + \Delta T) = [ H(T-2\Delta T)&
                     &-3 H(T-\Delta T )+ 3 H(T)]'
             else
                if (pub_on_root) write(stdout,*) &
                     &'pub_tddft_hamiltonian_mixing == 2 and so&
                     & H(T + 0.5 \Delta T) = 0.125 [3H(T-2\Delta T)&
                     &-10H(T-\Delta T )+15H(T)]'
             endif
             st2:if (timestep == 1) then
                do is=1,pub_cell%num_spins
                   call sparse_copy(midpoint_ham(is),current_ham(is))
                enddo
             elseif (timestep == 2) then
                do is=1,pub_cell%num_spins
                   call sparse_copy(midpoint_ham(is),current_ham(is))
                   if ((pub_tddft_propagation_method == &
                        &'RK_TIME_REVERSAL') .or. &
                        &(pub_tddft_propagation_method == &
                        &'RK_SC_TIME_REVERSAL')) then
                      call sparse_scale(midpoint_ham(is),2.0_DP)
                      call sparse_axpy(midpoint_ham(is),&
                           current_ham_1(is),-1.0_DP)
                   else
                      call sparse_scale(midpoint_ham(is),1.5_DP)
                      call sparse_axpy(midpoint_ham(is),&
                           current_ham_1(is),-0.5_DP)
                   endif
                enddo
             elseif (timestep .gt. 2) then
                do is=1,pub_cell%num_spins
                   call sparse_copy(midpoint_ham(is),current_ham(is))
                   if ((pub_tddft_propagation_method == &
                        &'RK_TIME_REVERSAL') .or. &
                        &(pub_tddft_propagation_method == &
                        'RK_SC_TIME_REVERSAL')) then
                      call sparse_scale(midpoint_ham(is),3.0_DP)
                      call sparse_axpy(midpoint_ham(is),&
                           current_ham_1(is),-3.0_DP)
                      call sparse_axpy(midpoint_ham(is),&
                           current_ham_2(is),1.0_DP)
                   else
                      call sparse_scale(midpoint_ham(is),1.875_DP)
                      call sparse_axpy(midpoint_ham(is),&
                           current_ham_1(is),-1.25_DP)
                      call sparse_axpy(midpoint_ham(is),&
                           current_ham_2(is),0.375_DP)
                   endif
                enddo
             endif st2
          elseif (pub_tddft_hamiltonian_mixing == 1) then
             if ((pub_tddft_propagation_method == &
                  &'RK_TIME_REVERSAL') .or. &
                  &(pub_tddft_propagation_method == &
                  &'RK_SC_TIME_REVERSAL')) then
                if (pub_on_root) write(stdout,*) &
                     &'pub_tddft_hamiltonian_mixing == 1 and so&
                     & H(T + \Delta T) = &
                     &0.5 [-H(T-\Delta T )+3H(T)]'
             else
                if (pub_on_root) write(stdout,*) &
                     &'pub_tddft_hamiltonian_mixing == 1 and so&
                     & H(T + 0.5 \Delta T) = &
                     &[-H(T-\Delta T )+2H(T)]'
             endif
             st1:if (timestep == 1) then
                do is=1,pub_cell%num_spins
                   call sparse_copy(midpoint_ham(is),current_ham(is))
                enddo
             elseif (timestep .ge. 2) then
                do is=1,pub_cell%num_spins
                   call sparse_copy(midpoint_ham(is),current_ham(is))
                   if ((pub_tddft_propagation_method == &
                        &'RK_TIME_REVERSAL') .or. &
                        &(pub_tddft_propagation_method == &
                        &'RK_SC_TIME_REVERSAL')) then
                      call sparse_scale(midpoint_ham(is),2.0_DP)
                      call sparse_axpy(midpoint_ham(is),&
                           current_ham_1(is),-1.0_DP)
                   else
                      call sparse_scale(midpoint_ham(is),1.5_DP)
                      call sparse_axpy(midpoint_ham(is),&
                           current_ham_1(is),-0.5_DP)
                   endif
                enddo
             endif st1
          elseif (pub_tddft_hamiltonian_mixing == 0) then
             ! Propagate with the retarded Hamiltonian
             ! from the previous density
             if (pub_on_root) write(stdout,*) &
                  &'pub_tddft_hamiltonian_mixing == 0 and so&
                  & H(T + 0.5 \Delta T) = H(T)'
             do is=1,pub_cell%num_spins
                call sparse_copy(midpoint_ham(is),current_ham(is))
             enddo
          endif mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PROPAGATION IN FULL LINEAR-RESPONSE APPROXIMATION !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  P (t + \Delta t) = exp (- i H0 \Delta t ) (P1)  exp (+ i H0 \Delta t ) + !
!!!!                     exp (- i H1 \Delta t ) (P0)  exp (+ i H1 \Delta t )   !
!!!!   in the new propagation scheme                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  P (t + \Delta t) = exp (- i (H0 + H1) \Delta t ) (P0 + P1)  !!!!!!!!!!!!!!
!!!!                     exp (+ i (H0 + H1) \Delta t )            !!!!!!!!!!!!!!
!!!!  in the old propagation scheme                               !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          propscheme: if (propagation_scheme_new) then

             ! Find the first-order change in the density matrix
             ! P1 = P - P_DFT.
             ! When on the first iteration,this is zero.
             do is=1,pub_cell%num_spins
                call sparse_copy(re_del_denskern(is),re_denskern(is))
                call sparse_axpy(re_del_denskern(is),&
                     dft_denskern(is),-1.0_DP)
                call sparse_copy(im_del_denskern(is),im_denskern(is))
                ! The perturbation counts as part of P1
                !   call sparse_axpy(im_del_denskern(is),&
                !        perturbation(is),-1.0_DP)
             enddo

             if (pub_on_root) write(stdout,*) &
                  &'Computing  exp (- i H0 \Delta t )&
                  & (P1)  exp (+ i H0 \Delta t )'

             ! The [H0,P1] component exp (- i H0 \Delta t ) (P1)
             ! exp (+ i H0 \Delta t ) is computed
             call tddft_generic_propagation(re_del_denskern, &
                  im_del_denskern,&
                  re_right, im_right, re_left, im_left, .true.)

             ! Copy this into the active density-matrix P_new
             do is=1,pub_cell%num_spins
                call sparse_copy(re_denskern(is),re_del_denskern(is))
                !im_del_denskern redundant here
                call sparse_copy(im_denskern(is),im_del_denskern(is))
             enddo

             if (pub_on_root) write(stdout,*) &
                  &'Re[Tr[K S]] and Im[Tr[K S]] for spin alpha &
                  & where K comes from the [H0, P1] component',&
                  sparse_trace(re_denskern(1),rep%overlap),&
                  sparse_trace(im_denskern(1),rep%overlap)

             ! This propagation then computes the [H1,P0] component,
             ! where there is no imaginary part to the Hamiltonian
             ! exp (- i H1 \Delta t ) (P0)  exp (+ i H1 \Delta t )

             ! Just use unperturbed (real) DFT kernel P0
             ! to generate the propagator
             do is=1,pub_cell%num_spins
                call sparse_copy(re_del_denskern(is),dft_denskern(is))
                !call sparse_copy(im_del_denskern(is),perturbation(is))
                call sparse_scale(im_del_denskern(is),0.0_DP)
             enddo
             ! Just use the first-order change in Hamiltonian H1 = H - H0
             do is=1,pub_cell%num_spins
                call sparse_axpy(midpoint_ham(is),dft_ham(is),-1.0_DP)
             enddo

             if (pub_on_root) write(stdout,*) &
                  &'Computing  exp (- i H1 \Delta t )&
                  & (P0)  exp (+ i H1 \Delta t )'

             ! ddor: Explicit Crank-Nicholson
             cn2: if (pub_tddft_propagation_method == 'CRANKNICHOLSON') then
                call tddft_crank_nicholson_prop(re_del_denskern, &
                     im_del_denskern,&
                     rep%inv_overlap, midpoint_ham, steplength, .true.)
             elseif (pub_tddft_propagation_method == 'RK_SC_TIME_REVERSAL') then
                ! ddor: Enforced time-reversal symmetry with Runge-Kutta 4
                do is=1,pub_cell%num_spins
                   call sparse_axpy(current_ham(is),dft_ham(is),-1.0_DP)
                enddo
                call tddft_runge_kutta_prop(re_del_denskern, &
                     im_del_denskern,&
                     rep%inv_overlap, current_ham, steplength*0.5_DP, .true.)
                do is=1,pub_cell%num_spins
                   call sparse_axpy(current_ham(is),dft_ham(is),1.0_DP)
                enddo
                !sc_time_reversal_rms = 0.0_DP
                !do is=1,pub_cell%num_spins
                !   sc_time_reversal_rms = sc_time_reversal_rms + &
                !        sparse_rms_element(midpoint_ham(is))
                !enddo
                self_consistency1: do sc_time_reversal_iter = 1, &
                     maxit_time_reversal

                   sc_time_reversal_energy = 0.0_DP

                   do is=1,pub_cell%num_spins
                      call sparse_copy(re_del_denskern_buf(is),&
                           re_del_denskern(is))
                      call sparse_copy(im_del_denskern_buf(is),&
                           im_del_denskern(is))
                   enddo
                   call tddft_runge_kutta_prop(re_del_denskern_buf, &
                        im_del_denskern_buf,&
                        rep%inv_overlap, midpoint_ham, steplength*0.5_DP, .true.)
                   do is=1,pub_cell%num_spins
                      call sparse_axpy(re_del_denskern_buf(is),&
                           re_denskern(is),1.0_DP)
                   enddo
                   ! Calculate the new energy and Hamiltonian
                   call tddft_electronic_energy(new_sc_time_reversal_energy,&
                        midpoint_ham,&
                        rep%ngwfs_on_grid, re_del_denskern_buf,&
                        rep%overlap, knl, localpseudo_fine, &
                        core_density_fine, ngwf_basis, ewald_energy)
                   if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                        & DEVIATION FROM INITIAL ENERGY ON &
                        &SELF-CONSISTENCY STEP ', &
                        sc_time_reversal_iter ,' IS ', &
                        new_sc_time_reversal_energy - initial_energy,'Hartree'
                   !new_sc_time_reversal_rms = 0.0_DP
                   !do is=1,pub_cell%num_spins
                   !   new_sc_time_reversal_rms = sc_time_reversal_rms + &
                   !        sparse_rms_element(midpoint_ham(is))
                   !enddo
                   if (ABS(sc_time_reversal_energy-&
                        &new_sc_time_reversal_energy) &
                        & .lt. 1.0e-15 ) exit self_consistency1
                   !sc_time_reversal_rms = new_sc_time_reversal_rms
                   sc_time_reversal_energy = new_sc_time_reversal_energy
                   ! Just use the first-order change in Hamiltonian H1 = H - H0
                   do is=1,pub_cell%num_spins
                      call sparse_axpy(midpoint_ham(is),dft_ham(is),-1.0_DP)
                   enddo
                enddo self_consistency1
                do is=1,pub_cell%num_spins
                   call sparse_axpy(re_del_denskern_buf(is),&
                        re_denskern(is),-1.0_DP)
                   call sparse_copy(re_del_denskern(is),re_del_denskern_buf(is))
                   call sparse_copy(im_del_denskern(is),im_del_denskern_buf(is))
                enddo
             elseif (pub_tddft_propagation_method == 'RK_TIME_REVERSAL') then
                ! ddor: Approximate enforced time-reversal symmetry
                !       fourth-order Runge-Kutta
                do is=1,pub_cell%num_spins
                   call sparse_axpy(current_ham(is),dft_ham(is),-1.0_DP)
                enddo
                call tddft_runge_kutta_prop(re_del_denskern, &
                     im_del_denskern,&
                     rep%inv_overlap, current_ham, steplength*0.5_DP, .true.)
                do is=1,pub_cell%num_spins
                   call sparse_axpy(current_ham(is),dft_ham(is),1.0_DP)
                enddo
                call tddft_runge_kutta_prop(re_del_denskern, &
                     im_del_denskern,&
                     rep%inv_overlap, midpoint_ham, steplength*0.5_DP, .true.)
             else ! ddor: Exponential-midpoint-rule with Runge-Kutta 4
                call tddft_runge_kutta_prop(re_del_denskern, &
                     im_del_denskern,&
                     rep%inv_overlap, midpoint_ham, steplength, .true.)
             endif cn2

             if (pub_on_root) write(stdout,*) &
                  &'Re[Tr[K S]] and Im[Tr[K S]] for spin alpha &
                  & where \delta K comes from &
                  &the [H1, P0] component',&
                  sparse_trace(re_del_denskern(1),rep%overlap),&
                  sparse_trace(im_del_denskern(1),rep%overlap)

             ! P (t + \Delta t) = exp (- i H0 \Delta t ) (P1)
             !                    exp (+ i H0 \Delta t ) +
             !                    exp (- i H1 \Delta t ) (P0)
             !                    exp (+ i H1 \Delta t )
             ! is then obtained
             do is=1,pub_cell%num_spins
                call sparse_axpy(re_denskern(is),re_del_denskern(is),1.0_DP)
                call sparse_axpy(im_denskern(is),im_del_denskern(is),1.0_DP)
             enddo

          else propscheme

             ! P (t + \Delta t) = exp (- i (H0 + H1) \Delta t ) (P0 + P1)
             !                    exp (+ i (H0 + H1) \Delta t )
             ! ddor: Explicit Crank-Nicholson
             cn3: if (pub_tddft_propagation_method == 'CRANKNICHOLSON') then
                call tddft_crank_nicholson_prop(re_denskern, im_denskern,&
                     rep%inv_overlap, midpoint_ham, steplength, .true.)
             elseif (pub_tddft_propagation_method == 'RK_SC_TIME_REVERSAL') then
                ! ddor: Enforced time-reversal symmetry with Runge-Kutta 4
                call tddft_runge_kutta_prop(re_denskern, &
                     im_denskern,&
                     rep%inv_overlap, current_ham, steplength*0.5_DP, .true.)
                !sc_time_reversal_rms = 0.0_DP
                !do is=1,pub_cell%num_spins
                !   sc_time_reversal_rms = sc_time_reversal_rms + &
                !        sparse_rms_element(midpoint_ham(is))
                !enddo
                self_consistency2: do sc_time_reversal_iter = 1, &
                     maxit_time_reversal

                   sc_time_reversal_energy = 0.0_DP

                   do is=1,pub_cell%num_spins
                      call sparse_copy(re_del_denskern_buf(is),re_denskern(is))
                      call sparse_copy(im_del_denskern_buf(is),im_denskern(is))
                   enddo
                   call tddft_runge_kutta_prop(re_del_denskern_buf, &
                        im_del_denskern_buf,&
                        rep%inv_overlap, midpoint_ham, steplength*0.5_DP, .true.)
                   ! Calculate the new energy and Hamiltonian
                   call tddft_electronic_energy(new_sc_time_reversal_energy,&
                        midpoint_ham,&
                        rep%ngwfs_on_grid, re_del_denskern_buf,&
                        rep%overlap, knl, localpseudo_fine, &
                        core_density_fine, ngwf_basis, ewald_energy)
                   if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                        & DEVIATION FROM INITIAL ENERGY ON &
                        &SELF-CONSISTENCY STEP ', &
                        sc_time_reversal_iter ,' IS ', &
                        new_sc_time_reversal_energy - initial_energy,'Hartree'
                   !new_sc_time_reversal_rms = 0.0_DP
                   !do is=1,pub_cell%num_spins
                   !   new_sc_time_reversal_rms = sc_time_reversal_rms + &
                   !        sparse_rms_element(midpoint_ham(is))
                   !enddo
                   if (ABS(sc_time_reversal_energy-new_sc_time_reversal_energy) &
                        & .lt. 1.0e-15 ) exit self_consistency2
                   !sc_time_reversal_rms = new_sc_time_reversal_rms
                   sc_time_reversal_energy = new_sc_time_reversal_energy
                enddo self_consistency2
                do is=1,pub_cell%num_spins
                   call sparse_copy(re_denskern(is),re_del_denskern_buf(is))
                   call sparse_copy(im_denskern(is),im_del_denskern_buf(is))
                enddo
             elseif (pub_tddft_propagation_method == 'RK_TIME_REVERSAL') then
                ! ddor: Enforced time-reversal symmetry with Runge-Kutta 4
                call tddft_runge_kutta_prop(re_denskern, im_denskern,&
                     rep%inv_overlap, current_ham, steplength*0.5_DP, .true.)
                call tddft_runge_kutta_prop(re_denskern, im_denskern,&
                     rep%inv_overlap, midpoint_ham, steplength*0.5_DP, .true.)
             else ! Exponential-midpoint-rule with Runge-Kutta 4
                call tddft_runge_kutta_prop(re_denskern, im_denskern,&
                     rep%inv_overlap, midpoint_ham, steplength, .true.)
             endif cn3

          endif propscheme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!  END PROPAGATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! OUTPUT THE OCCUPANCY AND IDEMPOTENCY ERROR
          if (pub_tddft_enforced_idempotency) then
             call tddft_idempotency(timestep, re_denskern, rep%overlap, &
                  im_denskern, dft_denskern, .true.)
          else
             call tddft_idempotency(timestep, re_denskern, rep%overlap, &
                  im_denskern)
          endif

          ! IF MIXING HAMILTONIAN, STORE BEFORE EVALUATING NEW ONE
          if ((pub_tddft_hamiltonian_mixing == 2) .and.(timestep .ge. 2)) then
             do is=1,pub_cell%num_spins
                call sparse_copy(current_ham_2(is),current_ham_1(is))
             enddo
          endif
          if ((pub_tddft_hamiltonian_mixing == 2) .or.&
               &(pub_tddft_hamiltonian_mixing == 1)) then
             do is=1,pub_cell%num_spins
                call sparse_copy(current_ham_1(is),current_ham(is))
             enddo
          endif

          ! If carrying out self-consistent propagation then
          ! the last calculation of the Hamiltonian is good.
          if  (pub_tddft_propagation_method == 'RK_SC_TIME_REVERSAL') then
             do is=1,pub_cell%num_spins
                call sparse_copy(current_ham(is),midpoint_ham(is))
             enddo
             if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                  & DEVIATION FROM INITIAL ENERGY ON TIMESTEP ', &
                  timestep ,' IS ', new_sc_time_reversal_energy - &
                  &initial_energy,'Hartree'
          else
             ! Calculate the new energy and Hamiltonian
             call tddft_electronic_energy(energy, current_ham,&
                  rep%ngwfs_on_grid, re_denskern,&
                  rep%overlap, knl, localpseudo_fine, &
                  core_density_fine, ngwf_basis, ewald_energy)
             ! OUTPUT THE ENERGY
             if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                  & DEVIATION FROM INITIAL ENERGY ON TIMESTEP ', &
                  timestep ,' IS ', energy - initial_energy,'Hartree'
          endif

          ! Calculate the change in dipole moment
          if (pub_on_root) write(stdout,*) 'tddft_calculate: &
               &DIPOLE MOMENT ON TIMESTEP ',timestep
          call properties_polarisation(rep,ngwf_basis,proj_basis, &
               elements,re_denskern,cartesian,real_output_dipole)
          ! END POLARISATION

          ! Time-dependent change in dipole moment
          dipole_t(timestep+1,cartesian) = &
               CMPLX(real_output_dipole,0.0_DP,kind=DP)
          dipole_t(timestep+1,cartesian) = dipole_t(timestep+1,cartesian) - &
               dipole_t_initial(cartesian)

          if (pub_on_root) write(stdout,*) &
               &'Change in time-dependent dipole moment (',timestep,&
               ',',cartesian,') =', dipole_t(timestep+1,cartesian)

       enddo interval

       !if (pub_on_root) then
       !   write(stdout,*) 'Time-dependent polarisation(',&
       !        cartesian,')'
       !   do timestep = 1, (spectrum_finesse_factor*num_timesteps)+1
       !      write(stdout,*) timestep,REAL(dipole_t(timestep,cartesian))
       !   enddo
       !endif


       ! Calculating polarisability and
       ! damping so that transitions with
       ! infinite lifetime may be transformed
       do timestep = 1, num_timesteps
          !if (pub_on_root) write(stdout,*) 'Damping factor on step',&
          !     timestep,' is ',&
          !     EXP( -1.0_DP * damping_exponent * steplength * &
          !     &REAL(timestep,kind=DP) / duration )

          !dipole_t(timestep+1,cartesian) =  dipole_t(timestep+1,cartesian) * &
          !     CMPLX((-0.5_DP/momentum),0.0_DP,kind=DP) * &
          !     EXP( CMPLX( -1.0_DP * damping_exponent * steplength * &
          !     &REAL(timestep,kind=DP) / duration , 0.0_DP ,kind=DP) )

          ! ddor: Experimentally using cubic filter instead of
          !       exponential damping
          !if (pub_on_root) write(stdout,*) 'Damping factor on step',&
          !     timestep,' is ',&
          !     &1.0_DP - (3.0_DP*((REAL(timestep,kind=DP)/&
          !     &REAL(num_timesteps,kind=DP))**2.0_DP)) + &
          !     &(2.0_DP*((REAL(timestep,kind=DP)/&
          !     &REAL(num_timesteps,kind=DP))**3.0_DP))

          dipole_t(timestep+1,cartesian) =  dipole_t(timestep+1,cartesian) * &
               CMPLX((-1.0_DP/momentum),0.0_DP,kind=DP) * &
               & ( 1.0_DP - (3.0_DP*((REAL(timestep,kind=DP)/&
               &REAL(num_timesteps,kind=DP))**2.0_DP)) + &
               &(2.0_DP*((REAL(timestep,kind=DP)/&
               &REAL(num_timesteps,kind=DP))**3.0_DP)) )

          ! ddor: Symmetry about midpoint as well as T=0
          ! dipole_t((2*num_timesteps*spectrum_finesse_factor)-&
          ! timestep+1,cartesian) = &
          !     dipole_t(timestep+1,cartesian) * CMPLX(-1.0_DP,0.0_DP,kind=DP)

          !if (pub_on_root) write(stdout,*) &
          !     &'Damped change in time-dependent polarisability (',&
          !     timestep,&
          !     ',',cartesian,') =', dipole_t(timestep+1,cartesian)
       enddo

       !do spectrum_finesse_counter = 1,(spectrum_finesse_factor-1)
       !   do timestep = 1, (2*num_timesteps)
       !      dipole_t(timestep+(spectrum_finesse_counter*2*num_timesteps),&
       !            &cartesian) = dipole_t(timestep,cartesian)
       !   enddo
       !enddo

       ! Frequency-dependent polarisability
       call wrappers_1d_fft('F', spectrum_finesse_factor*2*num_timesteps, &
            dipole_t(:,cartesian), dipole_w(:,cartesian))
       ! ddor: Make into sine transform
       !dipole_w(:,cartesian) = CMPLX(AIMAG(dipole_w(:,cartesian)),&
       !REAL(dipole_w(:,cartesian)),kind=DP)
       !if (pub_on_root) then
       !   write(stdout,*) 'Real part of damped change in polarisability_w(',&
       !        cartesian,')'
       !   do timestep = 1, (spectrum_finesse_factor*num_timesteps)+1
       !      write(stdout,*) timestep,REAL(dipole_w(timestep,cartesian))
       !   enddo
       !   write(stdout,*) 'Imaginary damped change in &
       !        &polarisability_w(',cartesian,')'
       !   do timestep = 1, (spectrum_finesse_factor*num_timesteps)+1
       !      write(stdout,*) timestep,AIMAG(dipole_w(timestep,cartesian))
       !   enddo
       !endif

       ! ddor: Sanity check
       !call slow_sine_transform(2*num_timesteps, steplength, &
       !     dipole_t(:,cartesian), dipole_w(:,cartesian))
       !if (pub_on_root) then
       !   write(stdout,*) 'Real part of damped change in polarisability_w(',&
       !        cartesian,')'
       !   do timestep = 1, num_timesteps+1
       !      write(stdout,*) timestep,REAL(dipole_w(timestep,cartesian))
       !   enddo
       !   write(stdout,*) 'Imaginary part of damped change in &
       !        &polarisability_w(',cartesian,')'
       !   do timestep = 1, num_timesteps+1
       !      write(stdout,*) timestep,AIMAG(dipole_w(timestep,cartesian))
       !   enddo
       !endif


    end do cart

    call tddft_observables(dipole_w, dipole_t,num_timesteps, &
         spectrum_finesse_factor, duration)!, rep%n_occ)

    ! Deallocate work arrays for dipole moment vectors
    deallocate(dipole_t_initial,stat=ierr)
    call utils_dealloc_check('tddft_calculate','dipole_t_initial',ierr)
    deallocate(dipole_w,stat=ierr)
    call utils_dealloc_check('tddft_calculate','dipole_w',ierr)
    deallocate(dipole_t,stat=ierr)
    call utils_dealloc_check('tddft_calculate','dipole_t',ierr)

    do is=pub_cell%num_spins,1,-1
       if (propagation_scheme_new) then
          call sparse_destroy(im_left(is))
          call sparse_destroy(re_left(is))
          call sparse_destroy(im_right(is))
          call sparse_destroy(re_right(is))
       endif
       if (pub_tddft_propagation_method == 'RK_SC_TIME_REVERSAL') then
          call sparse_destroy(im_del_denskern_buf(is))
          call sparse_destroy(re_del_denskern_buf(is))
       endif
       if (pub_tddft_hamiltonian_mixing == 2) then
          call sparse_destroy(current_ham_2(is))
       endif
       if ((pub_tddft_hamiltonian_mixing == 1) .or. &
            &(pub_tddft_hamiltonian_mixing == 2)) then
          call sparse_destroy(current_ham_1(is))
       endif
       if (propagation_scheme_new) then
          call sparse_destroy(dft_ham(is))
          call sparse_destroy(im_del_denskern(is))
          call sparse_destroy(re_del_denskern(is))
       endif
       call sparse_destroy(perturbation(is))
       call sparse_destroy(im_denskern(is))
       call sparse_destroy(re_denskern(is))
       call sparse_destroy(midpoint_ham(is))
       call sparse_destroy(current_ham(is))
    end do
    call sparse_destroy(knl)

    if (pub_tddft_propagation_method == 'RK_SC_TIME_REVERSAL') then
       deallocate(im_del_denskern_buf,stat=ierr)
       call utils_dealloc_check('tddft_calculate','im_del_denskern_buf',ierr)
       deallocate(re_del_denskern_buf,stat=ierr)
       call utils_dealloc_check('tddft_calculate','re_del_denskern_buf',ierr)
    endif
    if (pub_tddft_hamiltonian_mixing == 2) then
       deallocate(current_ham_2,stat=ierr)
       call utils_dealloc_check('tddft_calculate','current_ham_1',ierr)
    endif
    if ((pub_tddft_hamiltonian_mixing == 1) .or. &
         &(pub_tddft_hamiltonian_mixing == 2)) then
       deallocate(current_ham_1,stat=ierr)
       call utils_dealloc_check('tddft_calculate','current_ham_2',ierr)
    endif
    if (propagation_scheme_new) then
       deallocate(dft_ham,stat=ierr)
       call utils_dealloc_check('tddft_calculate','dft_ham',ierr)
       deallocate(im_left,stat=ierr)
       call utils_dealloc_check('tddft_calculate','im_left',ierr)
       deallocate(re_left,stat=ierr)
       call utils_dealloc_check('tddft_calculate','re_left',ierr)
       deallocate(im_right,stat=ierr)
       call utils_dealloc_check('tddft_calculate','im_right',ierr)
       deallocate(re_right,stat=ierr)
       call utils_dealloc_check('tddft_calculate','re_right',ierr)
       deallocate(im_del_denskern,stat=ierr)
       call utils_dealloc_check('tddft_calculate','im_del_denskern',ierr)
       deallocate(re_del_denskern,stat=ierr)
       call utils_dealloc_check('tddft_calculate','re_del_denskern',ierr)
    endif
    deallocate(im_denskern,stat=ierr)
    call utils_dealloc_check('tddft_calculate','im_denskern',ierr)
    deallocate(re_denskern,stat=ierr)
    call utils_dealloc_check('tddft_calculate','re_denskern',ierr)
    deallocate(midpoint_ham,stat=ierr)
    call utils_dealloc_check('tddft_calculate','current_ham',ierr)
    deallocate(current_ham,stat=ierr)
    call utils_dealloc_check('tddft_calculate','midpoint_ham',ierr)

    ! cks: stop timer
    call timer_clock("tddft_calculate",2)

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine tddft_perturbation_initiate(im_delta_denskern,  momentum, &
         initial_denskern, inverse_overlap,overlap, kick_axis,&
         ngwf_basis, ngwfs_on_grid)

      !========================================================================!
      ! This subroutine takes in the initial (ground-state) density matrix and !
      ! applies to it a phase shift corresponding to giving all electrons      !
      ! a small momentum boost of magnitude k = pub_tddft_dipole_kick_strength !
      ! in the cartesian direction kick_axis. The output is (purely imaginary) !
      ! first-order change in the density matrix  K0(t=0+) = -i[k.r , K0(t=0)].!
      !                                                                        !
      ! The density matrix is not repurified following perturbation in this    !
      ! subroutine. The phase shift is transformed into NGWF basis prior to    !
      ! application. The momentum boost is given in the input file in terms of !
      ! its imaginary energy in Hartree for each cartesian direction, though it!
      ! may be useful to consider the corresponding phase angle induced        !
      ! in order to ensure the linear-response regime is retained.             !
      !========================================================================!
      ! Arguments:                                                             !
      ! im_delta_denskern (output): Imaginary part of first-order change in the!
      !                             density matrix just after perturbation.    !
      !========================================================================!
      ! Written by David D. O'Regan in February 2009.                          !
      ! Modified 22/01/2011 by Nicholas Hine to use cell_grid_real_pt routine  !
      ! to save on storage.                                                    !
      !========================================================================!

      use cell_grid, only: cell_grid_real_pt, pub_fine_grid
      use integrals, only: integrals_locpot
      use rundat, only: pub_tddft_dipole_kick_strength
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create, sparse_axpy, &
           sparse_destroy, sparse_copy, sparse_product
      use timer, only: timer_clock
      use geometry, only: geometry_magnitude

      implicit none

      type(SPAM3), intent(inout) :: im_delta_denskern(pub_cell%num_spins)
      real(kind=DP), intent(out) :: momentum
      ! Unperturbed density matrix
      type(SPAM3), intent(in) :: initial_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: inverse_overlap ! inverse overlap matrix
      type(SPAM3), intent(in) :: overlap ! overlap matrix
      ! First order change in density matrix immediately
      ! following perturbation is imaginary K1(t=0) = -i[k.r , K0(t=0)]
      integer, intent(in) :: kick_axis
      type(FUNC_BASIS), intent(in) :: ngwf_basis ! function basis for NGWFs
      real(kind=DP), intent(in) :: ngwfs_on_grid(:) ! ngwfs on this proc


      ! Local variables
      integer :: i1,i2,islab12,is      ! Loop counters
      real(kind=DP), allocatable, dimension(:,:,:,:) :: perturbation_fine
      ! matrix elements of the dipole pulse
      type(SPAM3), allocatable :: dipole_kick(:)
      type(SPAM3) :: initial_polarisation, initial_polarisation_contra
      type(SPAM3) :: polarisation_initial, contra_polarisation_initial
      real(kind=DP) :: imaginary_energy
      real(kind=DP) :: rpt(3)

      ! start timer
      call timer_clock("tddft_perturbation_initiate",1)

      if (pub_on_root) write(stdout,*) 'Entering tddft_perturbation_initiate'

      ! Make some space in memory
      allocate(perturbation_fine(pub_fine_grid%ld1, pub_fine_grid%ld2, &
           pub_fine_grid%max_slabs12, pub_cell%num_spins),stat=ierr)
      call utils_alloc_check('tddft_perturbation_initiate',&
           'perturbation_fine',ierr)
      allocate(dipole_kick(pub_cell%num_spins),stat=ierr)
      call utils_alloc_check('tddft_perturbation_initiate','dipole_kick',ierr)
      dipole_kick%structure = 'S'
      do is=1,pub_cell%num_spins
         call sparse_create(dipole_kick(is))
      end do
      call sparse_create(initial_polarisation,&
           initial_denskern(1),dipole_kick(1))
      call sparse_create(initial_polarisation_contra,&
           initial_polarisation,inverse_overlap)
      call sparse_create(polarisation_initial,&
           dipole_kick(1),initial_denskern(1))
      call sparse_create(contra_polarisation_initial,&
           inverse_overlap,polarisation_initial)

      perturbation_fine(:,:,:,:) = 0.0_DP

      select case (kick_axis)
      case (1)
         momentum = PI * pub_tddft_dipole_kick_strength(kick_axis) / &
              geometry_magnitude(pub_cell%a1)
      case (2)
         momentum = PI * pub_tddft_dipole_kick_strength(kick_axis) / &
              geometry_magnitude(pub_cell%a2)
      case (3)
         momentum = PI * pub_tddft_dipole_kick_strength(kick_axis) / &
              geometry_magnitude(pub_cell%a3)
      end select
      if (pub_on_root) write(stdout,*) 'Momentum kick in direction (',&
           kick_axis,'):',momentum,'a_0^-1'
      imaginary_energy = (momentum**2.0_DP)/2.0_DP
      if (pub_on_root) write(stdout,*) 'Imaginary energy kick in direction (',&
           kick_axis,'):',imaginary_energy,'Ha'

      if (pub_tddft_dipole_kick_strength(kick_axis) .ne. 0.0_DP) then
         ! Loop over real-space grid on this node
         do islab12=1,pub_fine_grid%num_my_slabs12
            do i2=1,pub_fine_grid%n2
               do i1=1,pub_fine_grid%n1
                  call cell_grid_real_pt(rpt,i1,i2,islab12,pub_fine_grid)

                  perturbation_fine(i1,i2,islab12,1) = &
                       perturbation_fine(i1,i2,islab12,1) + &
                       momentum*rpt(kick_axis)

                  if (pub_cell%num_spins == 2) then
                     if (pub_tddft_dipole_kick_strength(kick_axis) &
                          &.gt. 0.0_DP) then
                        perturbation_fine(i1,i2,islab12,2) = &
                             perturbation_fine(i1,i2,islab12,1)
                     elseif (pub_tddft_dipole_kick_strength(kick_axis) &
                          &.lt. 0.0_DP) then
                        perturbation_fine(i1,i2,islab12,2) = &
                             -1.0_DP * perturbation_fine(i1,i2,islab12,1)
                     endif
                  endif

               end do
            end do
         end do

         ! Calculate the dipole perturbation matrix
         do is=1,pub_cell%num_spins
            call integrals_locpot(dipole_kick(is), &
                 ngwfs_on_grid,ngwf_basis,ngwfs_on_grid,ngwf_basis, &
                 pub_fine_grid,perturbation_fine(:,:,:,is))

            call sparse_product(initial_polarisation,&
                 initial_denskern(is),dipole_kick(is))
            call sparse_product(initial_polarisation_contra,&
                 initial_polarisation,&
                 inverse_overlap)
            call sparse_copy(im_delta_denskern(is),&
                 initial_polarisation_contra)

            call sparse_product(polarisation_initial,&
                 dipole_kick(is),initial_denskern(is))
            call sparse_product(contra_polarisation_initial,&
                 inverse_overlap,&
                 polarisation_initial)
            call sparse_axpy(im_delta_denskern(is),&
                 contra_polarisation_initial,-1.0_DP)
         end do
      endif

      do is=1,pub_cell%num_spins
         if (pub_on_root) write(stdout,*) &
              &' sparse_trace(overlap,im_delta_denskern(',is,')',&
              sparse_trace(overlap,im_delta_denskern(is))
      enddo

      ! Clear memory
      call sparse_destroy(contra_polarisation_initial)
      call sparse_destroy(polarisation_initial)
      call sparse_destroy(initial_polarisation_contra)
      call sparse_destroy(initial_polarisation)
      do is=pub_cell%num_spins,1,-1
         call sparse_destroy(dipole_kick(is))
      end do
      deallocate(dipole_kick,stat=ierr)
      call utils_dealloc_check('tddft_perturbation_initiate','dipole_kick',ierr)
      deallocate(perturbation_fine,stat=ierr)
      call utils_dealloc_check('tddft_perturbation_initiate',&
           &'perturbation_fine',ierr)

      if (pub_on_root) write(stdout,*) 'Leaving tddft_perturbation_initiate'

      ! end timer
      call timer_clock("tddft_perturbation_initiate",2)

    end subroutine tddft_perturbation_initiate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine tddft_electronic_energy(energy, ham, ngwfs_on_grid, re_denskern,&
         overlap, kinnonloc, localpseudo_fine, core_density_fine, ngwf_basis, &
         ewald_energy)

      !==================================================================!
      ! This function returns the total electronic energy for a given    !
      ! set of NGWFs with a given real part of the density matrix.       !
      ! Hartree, Hartree-Fock, XC are independent of Im[K] to 1st order  !
      ! Based on electronic_energy_no_kern_updt subroutine.              !
      !------------------------------------------------------------------!
      ! Written by David D. O'Regan in February 2009.                    !
      !==================================================================!

      use cell_grid, only: pub_fine_grid
      use constants, only: DP, stdout
      use hamiltonian, only: hamiltonian_lhxc_calculate
      use rundat, only: pub_usehfx, pub_aug
      use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale, &
           sparse_trace
      use simulation_cell, only: pub_cell
      use timer, only: timer_clock
      use vdwcorrection, only: pub_dispersion_energy
      use xc, only: xc_exit, xc_init

      implicit none

      ! Arguments
      real(kind=DP), intent(out) :: energy ! The energy
      type(SPAM3), intent(inout) :: ham(pub_cell%num_spins) ! The Hamiltonian
      type(FUNC_BASIS), intent(in) :: ngwf_basis
      real(kind=DP), intent(in) :: &
           ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
      ! Real part of the density kernel
      type(SPAM3), intent(inout) :: re_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: overlap
      ! Kinetic and non-local pseudo contributions to Hamiltonian
      type(SPAM3), intent(in) :: kinnonloc
      real(kind=DP), intent(in) :: localpseudo_fine(pub_fine_grid%ld1, &
           pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
      real(kind=DP), intent(in) :: core_density_fine(:,:,:)
      ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
      real(kind=DP), intent(in) :: ewald_energy

      ! Local variables
      ! lhxc on this proc slab
      real(kind=DP), allocatable, dimension(:,:,:,:) :: lhxc_fine
      real(kind=DP) :: lhxc_energy       ! pseudo+hartree+xc energy
      integer :: is                      ! spin loop counter
      integer :: ierr                    ! memory allocation error flag
      type(SPAM3), allocatable :: hfexchange(:)     ! HF exchange matrix
      ! local pseudo + Hartree + XC matrix in SPAM3 format
      type(SPAM3) :: lhxc
      type(SPAM3), allocatable :: dijhat(:)

      ! start timer
      call timer_clock("tddft_electronic_energy",1)

      if (pub_on_root) write(stdout,*) 'Entering tddft_electronic_energy'

      ! Make some space in memory
      allocate(lhxc_fine(pub_fine_grid%ld1, pub_fine_grid%ld2, &
           pub_fine_grid%max_slabs12, pub_cell%num_spins),stat=ierr)
      call utils_alloc_check('tddft_electronic_energy','lhxc_fine',ierr)
      call sparse_create(lhxc,overlap)
      if (pub_usehfx) then
         allocate(hfexchange(pub_cell%num_spins),stat=ierr)
         call utils_alloc_check('tddft_electronic_energy','hfexchange',ierr)
         do is=1, pub_cell%num_spins
            call sparse_create(hfexchange(is), re_denskern(is))
         end do
      endif
      if (pub_aug) then
         allocate(dijhat(pub_cell%num_spins),stat=ierr)
         call utils_alloc_check('tddft_electronic_energy','dijhat',ierr)
         do is=1, pub_cell%num_spins
            dijhat(is)%structure = 'E'
            call sparse_create(dijhat(is))
         end do
      end if

      ! scale density matrix for unpolarised calculation
      if (pub_cell%num_spins == 1) call sparse_scale(re_denskern(1),2.0_DP)

      ! cks: calculate the lhxc potential on the fine grid
      call xc_exit
      call xc_init(.true.)
      call hamiltonian_lhxc_calculate(lhxc_fine,lhxc_energy,dijhat, &
           pub_fine_grid, &
           localpseudo_fine,core_density_fine,ngwfs_on_grid,ngwf_basis, &
           re_denskern,rep%overlap,rep%sp_overlap)
      call xc_exit
      call xc_init(.false.)

      ! cks: calculate the density-dependent part of the Hamiltonian matrix
      do is=1,pub_cell%num_spins
         call integrals_locpot(lhxc, &
              ngwfs_on_grid,ngwf_basis,ngwfs_on_grid,ngwf_basis, &
              pub_fine_grid, lhxc_fine(:,:,:,is))
         call sparse_copy(ham(is),lhxc)
      end do

      ! kinetic + nonlocal potential energy
      energy = 0.0_DP
      do is=1,pub_cell%num_spins
         energy = energy + sparse_trace(re_denskern(is),kinnonloc)
         call sparse_axpy(ham(is),kinnonloc,1.0_DP)
      end do

      !Calculate and include HF exchange if necessary
      if (pub_usehfx) then
         !       call hf_exchange_calculate(hfexchange,re_denskern,overlap,&
         !       ngwfs_on_grid,ngwf_spheres,n_ngwf_ppds,.false.)
         !         ! Add HF contribution to energy and Hamiltonian
         !         do is=1,pub_cell%num_spins
         !            energy = energy - 0.5_DP*sparse_trace(re_denskern(is),&
         !                 hfexchange(is))
         !            call sparse_axpy(ham(is),hfexchange(is),-1.0_DP)
         !         end do
         if (pub_on_root) write(stdout,*) 'Error in tddft_electronic_energy: &
              &Time-Dependent Hartree-Fock not yet implemented'
         call comms_abort
      end if

      ! kinetic + nonlocal + local potential, Hartree,
      ! exchange-correlation and Ewald energies
      energy = energy + lhxc_energy +  ewald_energy + pub_dispersion_energy

      ! undo scaling
      if (pub_cell%num_spins == 1) call sparse_scale(re_denskern(1),0.5_DP)

      ! Clear memory
      if (pub_aug) then
         do is=pub_cell%num_spins,1,-1
            call sparse_destroy(dijhat(is))
         end do
         deallocate(dijhat,stat=ierr)
         call utils_dealloc_check('tddft_electronic_energy','dijhat',ierr)
      end if
      if (pub_usehfx) then
         do is=pub_cell%num_spins,1,-1
            call sparse_destroy(hfexchange(is))
         end do
         deallocate(hfexchange,stat=ierr)
         call utils_dealloc_check('tddft_electronic_energy','hfexchange',&
              ierr)
      endif
      call sparse_destroy(lhxc)
      deallocate(lhxc_fine,stat=ierr)
      call utils_dealloc_check('tddft_electronic_energy','lhxc_fine',ierr)

      if (pub_on_root) write(stdout,*) 'Leaving tddft_electronic_energy'

      ! End timer
      call timer_clock("tddft_electronic_energy",2)

    end subroutine tddft_electronic_energy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine tddft_runge_kutta_prop(re_density, im_density,&
         inverse_overlap, hamiltonian, tddft_timestep, im_in)

      !========================================================================!
      ! This subroutine takes in a Hamiltonian and inverse overlap matrix and  !
      ! computes the 4th order Runge-Kutta expansion of the right-propagator   !
      ! exp(+ i H inv_s dt) giving the real and imaginary parts separately.    !
      !                                                                        !
      ! The input density matrix is then updated with this propagator          !
      !                                                                        !
      ! Matrix multiplication is in the NGWF basis with sparsity patterns      !
      ! determined by the user parameter TDDFT_SPARSITY_LEVEL.                 !
      ! The sparsity of re_prop and im_prop is determined elsewhere.           !
      ! 0 : full square matrices used to calculate propagator.                 !
      ! 1 : Sparsity patterns up to order (H inv_s)^4 allowed                  !
      ! 2 : Sparsity patterns up to order (H inv_s)^2 allowed                  !
      ! 3 : Sparsity patterns up to order (H inv_s) allowed                    !
      !                                                                        !
      ! Timestep dt determined by user parameter TDDFT_MAXIMUM_ENERGY          !
      !========================================================================!
      ! Local variables:                                                       !
      ! im_prop:          Imaginary part of 4th order of right-propagator      !
      !                   Im [exp(+ i H inv_s dt) ] =                          !
      !                   + (H inv_s dt) - (1/3!) (H inv_s dt)^3               !
      ! re_prop:          Re part of 4th order right-hand propagator           !
      !                   Re [exp(+ i H inv_s dt) ] =                          !
      !                   1 - (1/2!) (H inv_s dt)^2 + (1/4!) (H inv_s dt)^4    !
      !========================================================================!
      ! Written by David D. O'Regan in February 2009.                          !
      !========================================================================!

      use rundat, only: pub_tddft_sparsity_level
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_transpose,&
           sparse_destroy, sparse_copy, sparse_product, sparse_scale
      use timer, only: timer_clock

      implicit none

      ! Real part of density matrix
      type(SPAM3), intent(inout) :: re_density(pub_cell%num_spins)
      ! Imaginary part
      type(SPAM3), intent(inout) :: im_density(pub_cell%num_spins)
      ! inverse overlap matrix
      type(SPAM3), intent(in) :: inverse_overlap
      ! Hamiltonian matrix
      type(SPAM3), intent(in) :: hamiltonian(pub_cell%num_spins)
      ! Length of time interval \Delta t
      real(kind=DP), intent(in) :: tddft_timestep
      ! work with an imaginary part of input denskern
      logical, intent(in) :: im_in

      ! Local variables
      integer :: is      ! spin counter
      real(kind=DP) :: onefac, twofac, threefac, fourfac

      ! Powers of (H inv_s)
      type(SPAM3) :: hs, hshs, hshshs, hshshshs
      ! Real part of propagator from right:
      ! 1 - (1/2) (H inv_s dt)^2 + (1/4!) (H inv_s dt)^4
      type(SPAM3) :: re_prop
      ! Imaginary part of propagator from right:
      ! + (H inv_s dt) - (1/3!) (H inv_s dt)^3
      type(SPAM3) :: im_prop
      ! Real part of propagator from left:
      ! 1 - (1/2!) (inv_s H dt)^2 + (1/4!) (inv_s H dt)^4
      type(SPAM3) :: re_prop_adjoint
      ! Imaginary part of propagator from left:
      ! - (inv_s H dt) + (1/3!) (inv_s H dt)^3
      type(SPAM3) :: im_prop_adjoint
      ! General purpose matrices of sparsity K HSHSHSHS
      type(SPAM3) :: re_dens_re_prop, re_dens_im_prop
      type(SPAM3) :: im_dens_re_prop, im_dens_im_prop
      ! and HSHSHSHS K HSHSHSHS
      type(SPAM3) :: aux_prop_dens_prop

      ! start timer
      call timer_clock("tddft_runge_kutta_prop",1)

      if (pub_on_root) write(stdout,*) 'Entering tddft_runge_kutta_prop'

      onefac = tddft_timestep
      twofac = -0.5_DP * (tddft_timestep**2.0_DP)
      threefac = (tddft_timestep**3.0_DP) / (-6.0_DP)
      fourfac = (tddft_timestep**4.0_DP) / (24.0_DP)

      ! Make some space in memory
      select case (pub_tddft_sparsity_level)
      case (0) ! FULL SQUARE MATRICES
         hs%structure = 'D'
         call sparse_create(hs)
         call sparse_create(hshs,hs)
         call sparse_create(hshshs,hs)
         call sparse_create(hshshshs,hs)
      case (1) ! ORDER (H inv_S)^4 ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(hshs,hs,hs)
         call sparse_create(hshshs,hshs,hs)
         call sparse_create(hshshshs,hshshs,hs)
      case (2) ! ORDER (H inv_S)^2 ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(hshs,hs,hs)
         call sparse_create(hshshs,hshs)
         call sparse_create(hshshshs,hshs)
      case (3) ! ORDER (H inv_S) ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(hshs,hs)
         call sparse_create(hshshs,hs)
         call sparse_create(hshshshs,hs)
      case default
         if (pub_on_root) write(stdout,*) 'Values of TDDFT_SPARSITY_LEVEL', &
              ' may lie in the range 0 to 3 inclusive.',&
              ' Setting to most conservative default 0.'
         hs%structure = 'D'
         call sparse_create(hs)
         call sparse_create(hshs,hs)
         call sparse_create(hshshs,hs)
         call sparse_create(hshshshs,hs)
      end select

      ! In all cases use biggest created structure for propagator
      ! and its adjoint
      call sparse_create(re_prop,hshshshs)
      call sparse_create(im_prop,hshshshs)
      call sparse_create(re_prop_adjoint,hshshshs)
      call sparse_create(im_prop_adjoint,hshshshs)
      call sparse_create(re_dens_re_prop,re_density(1),re_prop)
      call sparse_create(re_dens_im_prop,re_density(1),re_prop)
      if (im_in) then
         call sparse_create(im_dens_re_prop,re_density(1),re_prop)
         call sparse_create(im_dens_im_prop,re_density(1),re_prop)
      endif
      call sparse_create(aux_prop_dens_prop,re_prop_adjoint,re_dens_re_prop)

      ! SPARSITY PATTERNS FOR PROP_DENS_PROP
      ! CASE (0) : FKF
      ! CASE (1) : HSHSHSHS K HSHSHSHS
      ! CASE (2) : HSHS K HSHS
      ! CASE (3) : HS K HS
      ! May wish to consider explicilty creating (SH)^n patterns

      ! Calculate some powers of H(is) inv_s
      do is=1,pub_cell%num_spins

         call sparse_product(hs,hamiltonian(is),inverse_overlap)
         call sparse_product(hshs,hs,hs)
         call sparse_product(hshshs,hshs,hs)
         call sparse_product(hshshshs,hshshs,hs)

         ! REAL PART
         !  (H inv_s)^2
         call sparse_copy(re_prop,hshs)
         !  1 - (1/2!) (H inv_s dt)^2
         call sparse_scale(re_prop,twofac,1.0_DP)
         ! 1 - (1/2!) (H inv_s dt)^2 + (1/4!) (H inv_s dt)^4
         call sparse_axpy(re_prop,hshshshs,fourfac)
         !  1 - (1/2!) (inv_s H dt)^2 + (1/4!) (inv_s H dt)^4
         call sparse_transpose(re_prop_adjoint,re_prop)

         ! IMAGINARY PART
         ! H inv_s
         call sparse_copy(im_prop,hs)
         ! (H inv_s dt)
         call sparse_scale(im_prop,onefac,0.0_DP)
         ! (H inv_s dt) - (1/3!) (H inv_s dt)^3
         call sparse_axpy(im_prop,hshshs,threefac)
         ! (inv_s Hdt) - (1/3!) (inv_s H dt)^3
         call sparse_transpose(im_prop_adjoint,im_prop)
         ! -(inv_s H dt) + (1/3!) (inv_s H dt)^3
         call sparse_scale(im_prop_adjoint,-1.0_DP,0.0_DP)

         ! Re[P] Re[EXP(H)] and  Re[P] Im[EXP(H)]
         call sparse_product(re_dens_re_prop,re_density(is),re_prop)
         call sparse_product(re_dens_im_prop,re_density(is),im_prop)

         if (im_in) then
            ! Im[P] Re[EXP(H)] and  Im[P] Im[EXP(H)]
            call sparse_product(im_dens_re_prop,im_density(is),re_prop)
            call sparse_product(im_dens_im_prop,im_density(is),im_prop)
         endif

         ! Re[P_out] from Re[P_in]
         ! Re[EXP(-H)] Re[P] Re[EXP(H)]
         call sparse_product(re_density(is),re_prop_adjoint,re_dens_re_prop)
         call sparse_product(aux_prop_dens_prop,im_prop_adjoint,re_dens_im_prop)
         ! -Im[EXP(-H)] Re[P] Im[EXP(H)]
         call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)

         ! Im[P_out] from Re[P_in]
         ! i Im[EXP(-H)] Re[P] Re[EXP(H)]
         call sparse_product(im_density(is),im_prop_adjoint,re_dens_re_prop)
         call sparse_product(aux_prop_dens_prop,re_prop_adjoint,re_dens_im_prop)
         ! i Re[EXP(-H)] Re[P] Im[EXP(H)]
         call sparse_axpy(im_density(is),aux_prop_dens_prop,1.0_DP)

         if (im_in) then
            ! Re[P_out] from Im[P_in]
            call sparse_product(aux_prop_dens_prop,&
                 re_prop_adjoint,im_dens_im_prop)
            ! -Re[EXP(-H)] Im[P] Im[EXP(H)]
            call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)
            call sparse_product(aux_prop_dens_prop,&
                 im_prop_adjoint,im_dens_re_prop)
            ! -Im[EXP(-H)] Im[P] Re[EXP(H)]
            call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)

            ! Im[P_out] from Im[P_in]
            call sparse_product(aux_prop_dens_prop,&
                 re_prop_adjoint,im_dens_re_prop)
            ! i Re[EXP(-H)] Im[P] Re[EXP(H)]
            call sparse_axpy(im_density(is),aux_prop_dens_prop,1.0_DP)
            call sparse_product(aux_prop_dens_prop,&
                 im_prop_adjoint,im_dens_im_prop)
            ! -i Im[EXP(-H)] Im[P] Im[EXP(H)]
            call sparse_axpy(im_density(is),aux_prop_dens_prop,-1.0_DP)
         endif

      end do

      ! Clear memory
      ! EXP(-H) P EXP(H)
      call sparse_destroy(aux_prop_dens_prop)
      ! P EXP(H)
      if (im_in) then
         call sparse_destroy(im_dens_im_prop)
         call sparse_destroy(im_dens_re_prop)
      endif
      call sparse_destroy(re_dens_im_prop)
      call sparse_destroy(re_dens_re_prop)
      ! EXP(-H)
      call sparse_destroy(im_prop_adjoint)
      call sparse_destroy(re_prop_adjoint)
      ! EXP(H)
      call sparse_destroy(im_prop)
      call sparse_destroy(re_prop)
      ! (HS^-1)^n
      call sparse_destroy(hshshshs)
      call sparse_destroy(hshshs)
      call sparse_destroy(hshs)
      call sparse_destroy(hs)

      if (pub_on_root) write(stdout,*) 'Leaving tddft_runge_kutta_prop'

      ! end timer
      call timer_clock("tddft_runge_kutta_prop",2)

    end subroutine tddft_runge_kutta_prop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine tddft_runge_kutta_prop_init(&
         inverse_overlap, hamiltonian, tddft_timestep, &
         re_prop, im_prop, re_prop_adjoint, im_prop_adjoint)

      !========================================================================!
      ! This subroutine takes in a Hamiltonian and inverse overlap matrix and  !
      ! computes the 4th order Runge-Kutta expansion of the right-propagator   !
      ! exp(+ i H inv_s dt) giving the real and imaginary parts separately.    !
      !                                                                        !
      ! Matrix multiplication is in the NGWF basis with sparsity patterns      !
      ! determined by the user parameter TDDFT_SPARSITY_LEVEL.                 !
      ! The sparsity of re_prop and im_prop is determined elsewhere.           !
      ! 0 : full square matrices used to calculate propagator.                 !
      ! 1 : Sparsity patterns up to order (H inv_s)^4 allowed                  !
      ! 2 : Sparsity patterns up to order (H inv_s)^2 allowed                  !
      ! 3 : Sparsity patterns up to order (H inv_s) allowed                    !
      !                                                                        !
      ! Timestep dt determined by user parameter TDDFT_MAXIMUM_ENERGY          !
      !========================================================================!
      ! Local variables:                                                       !
      ! im_prop:          Imaginary part of 4th order of right-propagator      !
      !                   Im [exp(+ i H inv_s dt) ] =                          !
      !                   + (H inv_s dt) - (1/3!) (H inv_s dt)^3               !
      ! re_prop:          Re part of 4th right-hand propagator                 !
      !                   Re [exp(+ i H inv_s dt) ] =                          !
      !                   1 - (1/2!) (H inv_s dt)^2 + (1/4!) (H inv_s dt)^4    !
      !========================================================================!
      ! Written by David D. O'Regan in March 2009.                             !
      !========================================================================!

      use rundat, only: pub_tddft_sparsity_level
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_transpose,&
           sparse_destroy, sparse_copy, sparse_product, sparse_scale
      use timer, only: timer_clock

      implicit none

      type(SPAM3), intent(in) :: inverse_overlap ! inverse overlap matrix
      ! Hamiltonian matrix
      type(SPAM3), intent(in) :: hamiltonian(pub_cell%num_spins)
      real(kind=DP), intent(in) :: tddft_timestep

      ! Local variables
      integer :: is      ! spin counter
      real(kind=DP) :: onefac, twofac, threefac, fourfac

      ! Powers of (H inv_s)
      type(SPAM3) :: hs, hshs, hshshs, hshshshs
      ! Real part of propagator from right:
      ! 1 - (1/2!) (H inv_s dt)^2 + (1/4!) (H inv_s dt)^4
      type(SPAM3), intent(inout) :: re_prop(pub_cell%num_spins)
      ! Imaginary part of propagator from right:
      ! + (H inv_s dt) - (1/3!) (H inv_s dt)^3
      type(SPAM3), intent(inout) :: im_prop(pub_cell%num_spins)
      ! Real part of propagator from left:
      ! 1 - (1/2!) (inv_s H dt)^2 + (1/4!) (inv_s H dt)^4
      type(SPAM3), intent(inout) :: re_prop_adjoint(pub_cell%num_spins)
      ! Imaginary part of propagator from left:
      ! - (inv_s H dt) + (1/3!) (inv_s H dt)^3
      type(SPAM3), intent(inout) :: im_prop_adjoint(pub_cell%num_spins)

      ! start timer
      call timer_clock("tddft_runge_kutta_prop_init",1)

      if (pub_on_root) write(stdout,*) 'Entering tddft_runge_kutta_prop_init'

      onefac = tddft_timestep
      twofac = -0.5_DP * (tddft_timestep**2.0_DP)
      threefac = (tddft_timestep**3.0_DP) / (-6.0_DP)
      fourfac = (tddft_timestep**4.0_DP) / (24.0_DP)

      ! Make some space in memory
      select case (pub_tddft_sparsity_level)
      case (0) ! FULL SQUARE MATRICES
         hs%structure = 'D'
         call sparse_create(hs)
         call sparse_create(hshs,hs)
         call sparse_create(hshshs,hs)
         call sparse_create(hshshshs,hs)
      case (1) ! ORDER (H inv_S)^4 ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(hshs,hs,hs)
         call sparse_create(hshshs,hshs,hs)
         call sparse_create(hshshshs,hshshs,hs)
      case (2) ! ORDER (H inv_S)^2 ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(hshs,hs,hs)
         call sparse_create(hshshs,hshs)
         call sparse_create(hshshshs,hshs)
      case (3) ! ORDER (H inv_S) ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(hshs,hs)
         call sparse_create(hshshs,hs)
         call sparse_create(hshshshs,hs)
      case default
         if (pub_on_root) write(stdout,*) 'Values of TDDFT_SPARSITY_LEVEL', &
              ' may lie in the range 0 to 3 inclusive.',&
              ' Setting to most conservative default 0.'
         hs%structure = 'D'
         call sparse_create(hs)
         call sparse_create(hshs,hs)
         call sparse_create(hshshs,hs)
         call sparse_create(hshshshs,hs)
      end select

      ! In all cases use biggest created structure for propagator
      ! and its adjoint
      do is=1,pub_cell%num_spins
         call sparse_create(re_prop(is),hshshshs)
         call sparse_create(im_prop(is),hshshshs)
         call sparse_create(re_prop_adjoint(is),hshshshs)
         call sparse_create(im_prop_adjoint(is),hshshshs)
      enddo

      ! May wish to consider explicilty creating (SH)^n patterns

      ! Calculate some powers of H(is) inv_s
      do is=1,pub_cell%num_spins

         call sparse_product(hs,hamiltonian(is),inverse_overlap)
         call sparse_product(hshs,hs,hs)
         call sparse_product(hshshs,hshs,hs)
         call sparse_product(hshshshs,hshshs,hs)

         ! REAL PART
         !  (H inv_s)^2
         call sparse_copy(re_prop(is),hshs)
         !  1 - (1/2!) (H inv_s dt)^2
         call sparse_scale(re_prop(is),twofac,1.0_DP)
         ! 1 - (1/2!) (H inv_s dt)^2 + (1/4!) (H inv_s dt)^4
         call sparse_axpy(re_prop(is),hshshshs,fourfac)
         !  1 - (1/2!) (inv_s H dt)^2 + (1/4!) (inv_s H dt)^4
         call sparse_transpose(re_prop_adjoint(is),re_prop(is))

         ! IMAGINARY PART
         ! H inv_s
         call sparse_copy(im_prop(is),hs)
         ! (H inv_s dt)
         call sparse_scale(im_prop(is),onefac,0.0_DP)
         ! (H inv_s dt) - (1/3!) (H inv_s dt)^3
         call sparse_axpy(im_prop(is),hshshs,threefac)
         ! (inv_s Hdt) - (1/3!) (inv_s H dt)^3
         call sparse_transpose(im_prop_adjoint(is),im_prop(is))
         ! -(inv_s H dt) + (1/3!) (inv_s H dt)^3
         call sparse_scale(im_prop_adjoint(is),-1.0_DP,0.0_DP)

      end do

      ! Clear memory
      ! (HS^-1)^n
      call sparse_destroy(hshshshs)
      call sparse_destroy(hshshs)
      call sparse_destroy(hshs)
      call sparse_destroy(hs)

      if (pub_on_root) write(stdout,*) 'Leaving tddft_runge_kutta_prop_init'

      ! end timer
      call timer_clock("tddft_runge_kutta_prop_init",2)

    end subroutine tddft_runge_kutta_prop_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine tddft_crank_nicholson_prop(re_density, im_density,&
         inverse_overlap, hamiltonian, tddft_timestep, im_in)

      !========================================================================!
      ! This subroutine takes in a Hamiltonian and inverse overlap matrix and  !
      ! computes the 1st order Pade approximant to the right-propagator        !
      ! exp(+ i H inv_s dt) giving the real and imaginary parts separately.    !
      !                                                                        !
      ! The result is the unitary Crank-Nicholson propagator given by          !
      !                                                                        !
      !         1 + i 0.5 \Delta T H(T + 0.5 \Delta T) S^-1                    !
      ! U_CN =  ___________________________________________  + O(\Delta T^2)   !
      !                                                                        !
      !         1 - i 0.5 \Delta T H(T + 0.5 \Delta T) S^-1                    !
      !                                                                        !
      ! where the matrix inverse is carried out using Hotelling's algorithm.   !
      !                                                                        !
      ! The input density matrix is then updated with this propagator          !
      !                                                                        !
      ! Matrix multiplication is in the NGWF basis with sparsity patterns      !
      ! determined by the user parameter TDDFT_SPARSITY_LEVEL.                 !
      ! The sparsity of re_prop and im_prop is determined elsewhere.           !
      ! 0   : full square matrices used to calculate propagator.               !
      ! 1   : Sparsity patterns up to order (H inv_s)^2 and (inv_s H)^2 allowed!
      ! 2   : Sparsity patterns up to order (H inv_s)^2 allowed                !
      ! 3   : Sparsity patterns up to order (H inv_s) allowed                  !
      !                                                                        !
      ! In cases 2,3 the complex conjugate of the right-propagators is used    !
      ! for the left-propagator, with loss of information due to sparsity.     !
      !                                                                        !
      ! Timestep dt determined by user parameter TDDFT_MAXIMUM_ENERGY          !
      !========================================================================!
      !========================================================================!
      ! Written by David D. O'Regan in March 2009.                             !
      !========================================================================!

      use rundat, only: pub_tddft_sparsity_level, pub_tddft_maxit_hotelling
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_transpose,&
           sparse_destroy, sparse_copy, sparse_product, sparse_scale
      use timer, only: timer_clock
      use wrappers, only: wrappers_invert_sym_cmatrix

      implicit none

      ! Real part of density matrix
      type(SPAM3), intent(inout) :: re_density(pub_cell%num_spins)
      ! Imaginary part
      type(SPAM3), intent(inout) :: im_density(pub_cell%num_spins)
      ! inverse overlap matrix
      type(SPAM3), intent(in) :: inverse_overlap
      ! Hamiltonian matrix
      type(SPAM3), intent(in) :: hamiltonian(pub_cell%num_spins)
      real(kind=DP), intent(in) :: tddft_timestep
      ! work with an imaginary part of input denskern
      logical, intent(in) :: im_in

      ! For computing exact inverse matrices
      type(SPAM3) :: spam_squarebuffer
      real(kind=DP), allocatable :: squarebuffer(:,:)
      complex(kind=DP), allocatable :: c_squarebuffer(:,:)

      ! Local variables
      integer :: is      ! spin counter
      real(kind=DP) :: onefac

      ! Powers of (H inv_s) and denominator
      type(SPAM3) :: hs, im_hs, hs_inv, im_hs_inv
      ! Matrices for the denominator ! Not necessary if using left = right*
      type(SPAM3) :: sh, im_sh, sh_inv, im_sh_inv
      ! Real part of propagator from right:
      type(SPAM3) :: re_prop
      ! Imaginary part of propagator from right:
      type(SPAM3) :: im_prop
      ! Real part of propagator from left:
      type(SPAM3) :: re_prop_adjoint
      ! Imaginary part of propagator from left:
      type(SPAM3) :: im_prop_adjoint

      ! General purpose matrices of sparsity K HSHSHSHS
      type(SPAM3) :: re_dens_re_prop, re_dens_im_prop
      type(SPAM3) :: im_dens_re_prop, im_dens_im_prop
      type(SPAM3) :: aux_prop_dens_prop

      ! start timer
      call timer_clock("tddft_crank_nicholson_prop",1)

      if (pub_on_root) write(stdout,*) 'Entering tddft_crank_nicholson_prop'

      onefac = 0.5_DP * tddft_timestep

      if (pub_tddft_maxit_hotelling .eq. 0) then
         allocate(squarebuffer(ngwf_basis%num, ngwf_basis%num),stat=ierr)
         call utils_alloc_check('tddft_crank_nicholson_prop', &
              'squarebuffer',ierr)
         allocate(c_squarebuffer(ngwf_basis%num, ngwf_basis%num),stat=ierr)
         call utils_alloc_check('tddft_crank_nicholson_prop', &
              'c_squarebuffer',ierr)
         spam_squarebuffer%structure = 'D'
         call sparse_create(spam_squarebuffer)
      endif

      ! Make some space in memory
      select case (pub_tddft_sparsity_level)
      case (0) ! FULL SQUARE MATRICES
         hs%structure = 'D'
         call sparse_create(hs)
         call sparse_create(im_hs,hs)
         call sparse_create(hs_inv,hs)
         call sparse_create(im_hs_inv,hs)
         call sparse_create(re_prop,hs)
         call sparse_create(im_prop,hs)

         call sparse_create(sh,hs)
         call sparse_create(im_sh,hs)
         call sparse_create(sh_inv,hs)
         call sparse_create(im_sh_inv,hs)

         call sparse_create(re_prop_adjoint,hs)
         call sparse_create(im_prop_adjoint,hs)

      case (1) ! ORDER (H inv_S)^2 AND (inv_s H)^2 ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(im_hs,hs)
         call sparse_create(hs_inv,hs)
         call sparse_create(im_hs_inv,hs)
         call sparse_create(re_prop,hs,hs)
         call sparse_create(im_prop,hs,hs)

         call sparse_create(sh,inverse_overlap,hamiltonian(1))
         call sparse_create(im_sh,sh)
         call sparse_create(sh_inv,hs)
         call sparse_create(im_sh_inv,hs)

         call sparse_create(re_prop_adjoint,sh,sh)
         call sparse_create(im_prop_adjoint,sh,sh)

      case (2) ! ORDER (H inv_S)^2 ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(im_hs,hs)
         call sparse_create(hs_inv,hs)
         call sparse_create(im_hs_inv,hs)
         call sparse_create(re_prop,hs,hs)
         call sparse_create(im_prop,hs,hs)

         call sparse_create(re_prop_adjoint,hs,hs)
         call sparse_create(im_prop_adjoint,hs,hs)
      case (3) ! ORDER (H inv_S) ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(im_hs,hs)
         call sparse_create(hs_inv,hs)
         call sparse_create(im_hs_inv,hs)
         call sparse_create(re_prop,hs)
         call sparse_create(im_prop,hs)

         call sparse_create(re_prop_adjoint,hs)
         call sparse_create(im_prop_adjoint,hs)

      case default
         if (pub_on_root) write(stdout,*) 'Values of TDDFT_SPARSITY_LEVEL', &
              ' may lie in the range 0 to 3 inclusive.',&
              ' Setting to most conservative default 0.'
         hs%structure = 'D'
         call sparse_create(hs)
         call sparse_create(im_hs,hs)
         call sparse_create(re_prop,hs)
         call sparse_create(im_prop,hs)

         call sparse_create(sh,hs)
         call sparse_create(im_sh,hs)
         call sparse_create(sh_inv,hs)
         call sparse_create(im_sh_inv,hs)

         call sparse_create(re_prop_adjoint,hs)
         call sparse_create(im_prop_adjoint,hs)
      end select

      ! In all cases use biggest created structure for propagator
      ! and its adjoint
      call sparse_create(re_dens_re_prop,re_density(1),re_prop)
      call sparse_create(re_dens_im_prop,re_density(1),re_prop)
      if (im_in) then
         call sparse_create(im_dens_re_prop,re_density(1),re_prop)
         call sparse_create(im_dens_im_prop,re_density(1),re_prop)
      endif
      call sparse_create(aux_prop_dens_prop,re_prop_adjoint,re_dens_re_prop)

      ! SPARSITY PATTERNS FOR PROP_DENS_PROP
      ! CASE (0) : FKF
      ! CASE (1) : SHSH K HSHS
      ! CASE (2) : HSHS K HSHS
      ! CASE (3) : HS K HS

      do is=1,pub_cell%num_spins ! The main spin loop

         call sparse_product(hs,hamiltonian(is),inverse_overlap)
         call sparse_copy(im_hs,hs)
         call sparse_scale(hs,0.0_DP,1.0_DP) ! 1
         call sparse_scale(im_hs,-1.0_DP * onefac) ! - i 0.5 \Delta T H inv_s

         ! Compute (1 - i 0.5 \Delta T H inv_s)^-1
         if (pub_tddft_maxit_hotelling .eq. 0) then

            ! Calculate exact inverse by inverting
            ! the symmetric (inv_s - i 0.5 \Delta T inv_s H inv_s)
            call sparse_product(spam_squarebuffer,inverse_overlap,hs)
            call sparse_convert(squarebuffer,spam_squarebuffer)
            c_squarebuffer = CMPLX(squarebuffer, 0.0_DP, kind=DP)
            call sparse_product(spam_squarebuffer,inverse_overlap,im_hs)
            call sparse_convert(squarebuffer,spam_squarebuffer)
            c_squarebuffer = c_squarebuffer + &
                 CMPLX(0.0_DP,squarebuffer, kind=DP)

            ! Invert complex symmetric indefinite matrix
            call wrappers_invert_sym_cmatrix(c_squarebuffer, ngwf_basis%num)

            ! postmultiply by inv_s to remove the extra factor of S
            ! And put back into original sparsity pattern of H inv_s
            squarebuffer = REAL(c_squarebuffer)
            call sparse_convert(spam_squarebuffer, squarebuffer)
            call sparse_product(hs,spam_squarebuffer,inverse_overlap)
            squarebuffer = AIMAG(c_squarebuffer)
            call sparse_convert(spam_squarebuffer, squarebuffer)
            call sparse_product(im_hs,spam_squarebuffer,inverse_overlap)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         elseif  (pub_tddft_maxit_hotelling .gt. 0) then

            call tddft_hotelling_initial(hs_inv, hs, im_hs_inv, im_hs)
            call tddft_apply_hotelling(hs_inv, hs, im_hs_inv, &
                 im_hs, num_iter=pub_tddft_maxit_hotelling)

         elseif  (pub_tddft_maxit_hotelling .lt. 0) then

            write(stdout,*)'TDDFT_MAXIT HOTELLING LESS THAN 0. ONETEP stops'
            call comms_abort

         endif

         ! + i 0.5 \Delta T H inv_s for right numerator
         call sparse_scale(im_hs,-1.0_DP)
         ! REAL PART
         call sparse_product(re_prop,im_hs,im_hs_inv)
         call sparse_scale(re_prop,-1.0_DP)
         !Re[ (1 + i 0.5 \Delta T H inv_s) (1 - i 0.5 \Delta T H inv_s)^-1 ]
         call sparse_axpy(re_prop,hs_inv,1.0_DP)
         ! IMAGINARY PART
         call sparse_product(im_prop,im_hs,hs_inv)
         !Im[ (1 + i 0.5 \Delta T H inv_s) (1 - i 0.5 \Delta T H inv_s)^-1 ]
         call sparse_axpy(im_prop,im_hs_inv,1.0_DP)

         select case (pub_tddft_sparsity_level)
         case (0,1)
            call sparse_product(sh,inverse_overlap,hamiltonian(is))
            call sparse_copy(im_sh,sh)
            call sparse_scale(sh,0.0_DP,1.0_DP) ! 1
            call sparse_scale(im_sh,onefac) ! + i 0.5 \Delta T inv_s H

            ! Compute (1 + i 0.5 \Delta T inv_s H)^-1
            if  (pub_tddft_maxit_hotelling .eq. 0) then

               ! Calculate exact inverse by inverting
               ! the symmetric (inv_s + i 0.5 \Delta T inv_s H inv_s)
               call sparse_product(spam_squarebuffer,sh,inverse_overlap)
               call sparse_convert(squarebuffer,spam_squarebuffer)
               c_squarebuffer = CMPLX(squarebuffer, 0.0_DP, kind=DP)
               call sparse_product(spam_squarebuffer,im_sh,inverse_overlap)
               call sparse_convert(squarebuffer,spam_squarebuffer)
               c_squarebuffer = c_squarebuffer + &
                    CMPLX(0.0_DP,squarebuffer, kind=DP)

               ! Invert complex symmetric indefinite matrix
               call wrappers_invert_sym_cmatrix(c_squarebuffer, ngwf_basis%num)

               ! premultiply by inv_s to remove the extra factor of S
               ! And put back into original sparsity pattern of H inv_s
               squarebuffer = REAL(c_squarebuffer)
               call sparse_convert(spam_squarebuffer, squarebuffer)
               call sparse_product(sh,inverse_overlap,spam_squarebuffer)
               squarebuffer = AIMAG(c_squarebuffer)
               call sparse_convert(spam_squarebuffer, squarebuffer)
               call sparse_product(im_sh,inverse_overlap,spam_squarebuffer)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif  (pub_tddft_maxit_hotelling .gt. 0) then

               call tddft_hotelling_initial(sh_inv, sh, im_sh_inv, im_sh)
               call tddft_apply_hotelling(sh_inv, sh, &
                    im_sh_inv, im_sh, num_iter=pub_tddft_maxit_hotelling)

            elseif  (pub_tddft_maxit_hotelling .lt. 0) then

               write(stdout,*)'TDDFT_MAXIT HOTELLING LESS THAN 0. ONETEP stops'
               call comms_abort

            endif

            ! - i 0.5 \Delta T inv_s H for left numerator
            call sparse_scale(im_sh,-1.0_DP)
            ! REAL PART
            call sparse_product(re_prop_adjoint,im_sh_inv,im_sh)
            call sparse_scale(re_prop_adjoint,-1.0_DP)
            !Re[ (1 + i 0.5 \Delta T inv_s H)^-1  (1 - i 0.5 \Delta T inv_s H) ]
            call sparse_axpy(re_prop_adjoint,sh_inv,1.0_DP)
            ! IMAGINARY PART
            call sparse_product(im_prop_adjoint,sh_inv,im_sh)
            !Im[ (1 + i 0.5 \Delta T inv_s H)^-1  (1 - i 0.5 \Delta T inv_s H) ]
            call sparse_axpy(im_prop_adjoint,im_sh_inv,1.0_DP)

         case (2,3)
            call sparse_transpose(re_prop_adjoint,re_prop)
            call sparse_transpose(im_prop_adjoint,im_prop)
            call sparse_scale(im_prop_adjoint,-1.0_DP)

         case default
            if (pub_on_root) write(stdout,*) 'Values of TDDFT_SPARSITY_LEVEL', &
                 ' may lie in the range 0 to 3 inclusive.',&
                 ' Setting to most conservative default 0.'
            call sparse_product(sh,inverse_overlap,hamiltonian(is))
            call sparse_copy(im_sh,sh)
            call sparse_scale(sh,0.0_DP,1.0_DP) ! 1
            call sparse_scale(im_sh,onefac) ! + i 0.5 \Delta T inv_s H

            ! Compute (1 + i 0.5 \Delta T inv_s H)^-1
            if (pub_tddft_maxit_hotelling .eq. 0) then

               ! Calculate exact inverse by inverting
               ! the symmetric (inv_s + i 0.5 \Delta T inv_s H inv_s)
               call sparse_product(spam_squarebuffer,sh,inverse_overlap)
               call sparse_convert(squarebuffer,spam_squarebuffer)
               c_squarebuffer = CMPLX(squarebuffer, 0.0_DP, kind=DP)
               call sparse_product(spam_squarebuffer,im_sh,inverse_overlap)
               call sparse_convert(squarebuffer,spam_squarebuffer)
               c_squarebuffer = c_squarebuffer + &
                    CMPLX(0.0_DP,squarebuffer, kind=DP)

               ! Invert complex symmetric indefinite matrix
               call wrappers_invert_sym_cmatrix(c_squarebuffer, ngwf_basis%num)

               ! premultiply by inv_s to remove the extra factor of S
               ! And put back into original sparsity pattern of H inv_s
               squarebuffer = REAL(c_squarebuffer)
               call sparse_convert(spam_squarebuffer, squarebuffer)
               call sparse_product(sh,inverse_overlap,spam_squarebuffer)
               squarebuffer = AIMAG(c_squarebuffer)
               call sparse_convert(spam_squarebuffer, squarebuffer)
               call sparse_product(im_sh,inverse_overlap,spam_squarebuffer)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif  (pub_tddft_maxit_hotelling .gt. 0) then

               call tddft_hotelling_initial(sh_inv, sh, im_sh_inv, im_sh)
               call tddft_apply_hotelling(sh_inv, sh, im_sh_inv, &
                    im_sh,num_iter=maxit_hotelling)

            elseif  (pub_tddft_maxit_hotelling .lt. 0) then

               write(stdout,*)'TDDFT_MAXIT HOTELLING LESS THAN 0. ONETEP stops'
               call comms_abort

            endif

            ! - i 0.5 \Delta T inv_s H for left numerator
            call sparse_scale(im_sh,-1.0_DP)
            ! REAL PART
            call sparse_product(re_prop_adjoint,im_sh,im_sh_inv)
            call sparse_scale(re_prop_adjoint,-1.0_DP)
            !Re[  (1 - i 0.5 \Delta T inv_s H) (1 + i 0.5 \Delta T inv_s H)^-1 ]
            call sparse_axpy(re_prop_adjoint,sh_inv,1.0_DP)
            ! IMAGINARY PART
            call sparse_product(im_prop_adjoint,im_sh,sh_inv)
            !Im[  (1 - i 0.5 \Delta T inv_s H) (1 + i 0.5 \Delta T inv_s H)^-1 ]
            call sparse_axpy(im_prop_adjoint,im_sh_inv,1.0_DP)

         end select

         ! Re[P] Re[EXP(H)] and  Re[P] Im[EXP(H)]
         call sparse_product(re_dens_re_prop,re_density(is),re_prop)
         call sparse_product(re_dens_im_prop,re_density(is),im_prop)

         if (im_in) then
            ! Im[P] Re[EXP(H)] and  Im[P] Im[EXP(H)]
            call sparse_product(im_dens_re_prop,im_density(is),re_prop)
            call sparse_product(im_dens_im_prop,im_density(is),im_prop)
         endif

         ! Re[P_out] from Re[P_in]
         ! Re[EXP(-H)] Re[P] Re[EXP(H)]
         call sparse_product(re_density(is),re_prop_adjoint,re_dens_re_prop)
         call sparse_product(aux_prop_dens_prop,im_prop_adjoint,re_dens_im_prop)
         ! -Im[EXP(-H)] Re[P] Im[EXP(H)]
         call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)

         ! Im[P_out] from Re[P_in]
         ! i Im[EXP(-H)] Re[P] Re[EXP(H)]
         call sparse_product(im_density(is),im_prop_adjoint,re_dens_re_prop)
         call sparse_product(aux_prop_dens_prop,re_prop_adjoint,re_dens_im_prop)
         ! i Re[EXP(-H)] Re[P] Im[EXP(H)]
         call sparse_axpy(im_density(is),aux_prop_dens_prop,1.0_DP)

         if (im_in) then
            ! Re[P_out] from Im[P_in]
            call sparse_product(aux_prop_dens_prop,&
                 re_prop_adjoint,im_dens_im_prop)
            ! - Re[EXP(-H)] Im[P] Im[EXP(H)]
            call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)
            call sparse_product(aux_prop_dens_prop,&
                 im_prop_adjoint,im_dens_re_prop)
            ! - Im[EXP(-H)] Im[P] Re[EXP(H)]
            call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)

            ! Im[P_out] from Im[P_in]
            call sparse_product(aux_prop_dens_prop,&
                 re_prop_adjoint,im_dens_re_prop)
            ! i Re[EXP(-H)] Im[P] Re[EXP(H)]
            call sparse_axpy(im_density(is),aux_prop_dens_prop,1.0_DP)
            call sparse_product(aux_prop_dens_prop,&
                 im_prop_adjoint,im_dens_im_prop)
            ! -i Im[EXP(-H)] Im[P] Im[EXP(H)]
            call sparse_axpy(im_density(is),aux_prop_dens_prop,-1.0_DP)
         endif

      end do

      ! Clear memory
      ! EXP(-H) P EXP(H)
      call sparse_destroy(aux_prop_dens_prop)
      ! P EXP(H)
      if (im_in) then
         call sparse_destroy(im_dens_im_prop)
         call sparse_destroy(im_dens_re_prop)
      endif
      call sparse_destroy(re_dens_im_prop)
      call sparse_destroy(re_dens_re_prop)
      ! EXP(-H)
      call sparse_destroy(im_prop_adjoint)
      call sparse_destroy(re_prop_adjoint)

      select case (pub_tddft_sparsity_level)
      case (0,1)
         call sparse_destroy(im_sh_inv)
         call sparse_destroy(sh_inv)
         call sparse_destroy(im_sh)
         call sparse_destroy(sh)
      case (2,3)
      case default
         call sparse_destroy(im_sh_inv)
         call sparse_destroy(sh_inv)
         call sparse_destroy(im_sh)
         call sparse_destroy(sh)
      end select
      ! EXP(H)
      call sparse_destroy(im_prop)
      call sparse_destroy(re_prop)
      call sparse_destroy(im_hs_inv)
      call sparse_destroy(hs_inv)
      call sparse_destroy(im_hs)
      call sparse_destroy(hs)

      if (pub_tddft_maxit_hotelling .eq. 0) then
         call sparse_destroy(spam_squarebuffer)
         deallocate(c_squarebuffer,stat=ierr)
         call utils_dealloc_check('tddft_crank_nicholson_prop', &
              'c_squarebuffer',ierr)
         deallocate(squarebuffer,stat=ierr)
         call utils_dealloc_check('tddft_crank_nicholson_prop', &
              'squarebuffer',ierr)
      endif

      if (pub_on_root) write(stdout,*) 'Leaving tddft_crank_nicholson_prop'

      ! end timer
      call timer_clock("tddft_crank_nicholson_prop",2)

    end subroutine tddft_crank_nicholson_prop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine tddft_crank_nicholson_prop_init(&
         inverse_overlap, hamiltonian, tddft_timestep, &
         re_prop, im_prop, re_prop_adjoint, im_prop_adjoint)

      !========================================================================!
      ! This subroutine takes in a Hamiltonian and inverse overlap matrix and  !
      ! computes the 1st order Pade approximant to the right-propagator        !
      ! exp(+ i H inv_s dt) giving the real and imaginary parts separately.    !
      !                                                                        !
      ! The result is the unitary Crank-Nicholson propagator given by          !
      !                                                                        !
      !         1 + i 0.5 \Delta T H(T + 0.5 \Delta T) S^-1                    !
      ! U_CN =  ___________________________________________  + O(\Delta T^2)   !
      !                                                                        !
      !         1 - i 0.5 \Delta T H(T + 0.5 \Delta T) S^-1                    !
      !                                                                        !
      ! where the matrix inverse is carried out using Hotelling's algorithm.   !
      !                                                                        !
      ! Matrix multiplication is in the NGWF basis with sparsity patterns      !
      ! determined by the user parameter TDDFT_SPARSITY_LEVEL.                 !
      ! The sparsity of re_prop and im_prop is determined elsewhere.           !
      ! 0   : full square matrices used to calculate propagator.               !
      ! 1   : Sparsity patterns up to order (H inv_s)^2 and (inv_s H)^2 allowed!
      ! 2   : Sparsity patterns up to order (H inv_s)^2 allowed                !
      ! 3   : Sparsity patterns up to order (H inv_s) allowed                  !
      !                                                                        !
      ! In cases 2,3 the complex conjugate of the right-propagators is used    !
      ! for the left-propagator, with loss of information due to sparsity.     !
      !                                                                        !
      ! Timestep dt determined by user parameter TDDFT_MAXIMUM_ENERGY          !
      !========================================================================!
      !========================================================================!
      ! Written by David D. O'Regan in March 2009.                             !
      !========================================================================!

      use rundat, only: pub_tddft_sparsity_level, pub_tddft_maxit_hotelling
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_transpose,&
           sparse_destroy, sparse_copy, sparse_product, sparse_scale
      use timer, only: timer_clock
      use wrappers, only : wrappers_invert_sym_cmatrix

      implicit none

      type(SPAM3), intent(in) :: inverse_overlap ! inverse overlap matrix
      ! Hamiltonian matrix
      type(SPAM3), intent(in) :: hamiltonian(pub_cell%num_spins)
      real(kind=DP), intent(in) :: tddft_timestep

      ! Local variables
      integer :: is      ! spin counter
      real(kind=DP) :: onefac

      ! For computing exact inverse matrices
      type(SPAM3) :: spam_squarebuffer
      real(kind=DP), allocatable :: squarebuffer(:,:)
      complex(kind=DP), allocatable :: c_squarebuffer(:,:)

      ! Powers of (H inv_s) and denominator
      type(SPAM3) :: hs, im_hs, hs_inv, im_hs_inv
      ! Matrices for the denominator ! Not necessary if using left = right*
      type(SPAM3) :: sh, im_sh, sh_inv, im_sh_inv
      ! Real part of propagator from right:
      type(SPAM3), intent(inout) :: re_prop(pub_cell%num_spins)
      ! Imaginary part of propagator from right:
      type(SPAM3), intent(inout) :: im_prop(pub_cell%num_spins)
      ! Real part of propagator from left:
      type(SPAM3), intent(inout) :: re_prop_adjoint(pub_cell%num_spins)
      ! Imaginary part of propagator from left:
      type(SPAM3), intent(inout) :: im_prop_adjoint(pub_cell%num_spins)

      ! start timer
      call timer_clock("tddft_crank_nicholson_prop_init",1)

      if (pub_on_root) write(stdout,*) &
           &'Entering tddft_crank_nicholson_prop_init'

      onefac = 0.5_DP * tddft_timestep

      if (pub_tddft_maxit_hotelling .eq. 0) then
         allocate(squarebuffer(ngwf_basis%num, ngwf_basis%num),stat=ierr)
         call utils_alloc_check('tddft_crank_nicholson_prop', &
              'squarebuffer',ierr)
         allocate(c_squarebuffer(ngwf_basis%num, ngwf_basis%num),stat=ierr)
         call utils_alloc_check('tddft_crank_nicholson_prop', &
              'c_squarebuffer',ierr)
         spam_squarebuffer%structure = 'D'
         call sparse_create(spam_squarebuffer)
      endif

      ! Make some space in memory
      select case (pub_tddft_sparsity_level)
      case (0) ! FULL SQUARE MATRICES
         hs%structure = 'D'
         call sparse_create(hs)
         call sparse_create(im_hs,hs)
         call sparse_create(hs_inv,hs)
         call sparse_create(im_hs_inv,hs)
         do is=1,pub_cell%num_spins
            call sparse_create(re_prop(is),hs)
            call sparse_create(im_prop(is),hs)
         enddo

         call sparse_create(sh,hs)
         call sparse_create(im_sh,hs)
         call sparse_create(sh_inv,hs)
         call sparse_create(im_sh_inv,hs)

         do is=1,pub_cell%num_spins
            call sparse_create(re_prop_adjoint(is),hs)
            call sparse_create(im_prop_adjoint(is),hs)
         enddo

      case (1) ! ORDER (H inv_S)^2 AND (inv_s H)^2 ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(im_hs,hs)
         call sparse_create(hs_inv,hs)
         call sparse_create(im_hs_inv,hs)
         do is=1,pub_cell%num_spins
            call sparse_create(re_prop(is),hs,hs)
            call sparse_create(im_prop(is),hs,hs)
         enddo

         call sparse_create(sh,inverse_overlap,hamiltonian(1))
         call sparse_create(im_sh,sh)
         call sparse_create(sh_inv,hs)
         call sparse_create(im_sh_inv,hs)

         do is=1,pub_cell%num_spins
            call sparse_create(re_prop_adjoint(is),sh,sh)
            call sparse_create(im_prop_adjoint(is),sh,sh)
         enddo

      case (2) ! ORDER (H inv_S)^2 ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(im_hs,hs)
         call sparse_create(hs_inv,hs)
         call sparse_create(im_hs_inv,hs)
         do is=1,pub_cell%num_spins
            call sparse_create(re_prop(is),hs,hs)
            call sparse_create(im_prop(is),hs,hs)

            call sparse_create(re_prop_adjoint(is),hs,hs)
            call sparse_create(im_prop_adjoint(is),hs,hs)
         enddo

      case (3) ! ORDER (H inv_S) ALLOWED
         call sparse_create(hs,hamiltonian(1),inverse_overlap)
         call sparse_create(im_hs,hs)
         call sparse_create(hs_inv,hs)
         call sparse_create(im_hs_inv,hs)
         do is=1,pub_cell%num_spins
            call sparse_create(re_prop(is),hs)
            call sparse_create(im_prop(is),hs)

            call sparse_create(re_prop_adjoint(is),hs)
            call sparse_create(im_prop_adjoint(is),hs)
         enddo

      case default
         if (pub_on_root) write(stdout,*) 'Values of TDDFT_SPARSITY_LEVEL', &
              ' may lie in the range 0 to 3 inclusive.',&
              ' Setting to most conservative default 0.'
         hs%structure = 'D'
         call sparse_create(hs)
         call sparse_create(im_hs,hs)
         do is=1,pub_cell%num_spins
            call sparse_create(re_prop(is),hs)
            call sparse_create(im_prop(is),hs)
         enddo

         call sparse_create(sh,hs)
         call sparse_create(im_sh,hs)
         call sparse_create(sh_inv,hs)
         call sparse_create(im_sh_inv,hs)

         do is=1,pub_cell%num_spins
            call sparse_create(re_prop_adjoint(is),hs)
            call sparse_create(im_prop_adjoint(is),hs)
         enddo

      end select

      do is=1,pub_cell%num_spins ! The main spin loop

         call sparse_product(hs,hamiltonian(is),inverse_overlap)
         call sparse_copy(im_hs,hs)
         call sparse_scale(hs,0.0_DP,1.0_DP) ! 1
         call sparse_scale(im_hs,-1.0_DP * onefac) ! - i 0.5 \Delta T H inv_s

         ! Compute (1 - i 0.5 \Delta T H inv_s)^-1
         if (pub_tddft_maxit_hotelling .eq. 0) then

            ! Calculate exact inverse by inverting
            ! the symmetric (inv_s - i 0.5 \Delta T inv_s H inv_s)
            call sparse_product(spam_squarebuffer,inverse_overlap,hs)
            call sparse_convert(squarebuffer,spam_squarebuffer)
            c_squarebuffer = CMPLX(squarebuffer, 0.0_DP, kind=DP)
            call sparse_product(spam_squarebuffer,inverse_overlap,im_hs)
            call sparse_convert(squarebuffer,spam_squarebuffer)
            c_squarebuffer = c_squarebuffer + &
                 CMPLX(0.0_DP,squarebuffer, kind=DP)

            ! Invert complex symmetric indefinite matrix
            call wrappers_invert_sym_cmatrix(c_squarebuffer, ngwf_basis%num)

            ! postmultiply by inv_s to remove the extra factor of S
            ! And put back into original sparsity pattern of H inv_s
            squarebuffer = REAL(c_squarebuffer)
            call sparse_convert(spam_squarebuffer, squarebuffer)
            call sparse_product(hs,spam_squarebuffer,inverse_overlap)
            squarebuffer = AIMAG(c_squarebuffer)
            call sparse_convert(spam_squarebuffer, squarebuffer)
            call sparse_product(im_hs,spam_squarebuffer,inverse_overlap)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         elseif  (pub_tddft_maxit_hotelling .gt. 0) then

            call tddft_hotelling_initial(hs_inv, hs, im_hs_inv, im_hs)
            call tddft_apply_hotelling(hs_inv, hs, im_hs_inv, &
                 im_hs, num_iter=pub_tddft_maxit_hotelling)

         elseif  (pub_tddft_maxit_hotelling .lt. 0) then

            write(stdout,*)'TDDFT_MAXIT HOTELLING LESS THAN 0. ONETEP stops'
            call comms_abort

         endif

         call sparse_scale(im_hs,-1.0_DP) ! + i 0.5 \Delta T H inv_s
         ! REAL PART
         call sparse_product(re_prop(is),im_hs,im_hs_inv)
         call sparse_scale(re_prop(is),-1.0_DP)
         !Re[ (1 + i 0.5 \Delta T H inv_s) (1 - i 0.5 \Delta T H inv_s)^-1 ]
         call sparse_axpy(re_prop(is),hs_inv,1.0_DP)
         ! IMAGINARY PART
         call sparse_product(im_prop(is),im_hs,hs_inv)
         !Im[ (1 + i 0.5 \Delta T H inv_s) (1 - i 0.5 \Delta T H inv_s)^-1 ]
         call sparse_axpy(im_prop(is),im_hs_inv,1.0_DP)

         select case (pub_tddft_sparsity_level)
         case (0,1)
            call sparse_product(sh,inverse_overlap,hamiltonian(is))
            call sparse_copy(im_sh,sh)
            call sparse_scale(sh,0.0_DP,1.0_DP) ! 1
            call sparse_scale(im_sh,onefac) ! + i 0.5 \Delta T inv_s H

            ! Compute (1 + i 0.5 \Delta T inv_s H)^-1
            if  (pub_tddft_maxit_hotelling .eq. 0) then

               ! Calculate exact inverse by inverting
               ! the symmetric (inv_s + i 0.5 \Delta T inv_s H inv_s)
               call sparse_product(spam_squarebuffer,sh,inverse_overlap)
               call sparse_convert(squarebuffer,spam_squarebuffer)
               c_squarebuffer = CMPLX(squarebuffer, 0.0_DP, kind=DP)
               call sparse_product(spam_squarebuffer,im_sh,inverse_overlap)
               call sparse_convert(squarebuffer,spam_squarebuffer)
               c_squarebuffer = c_squarebuffer + &
                    CMPLX(0.0_DP,squarebuffer, kind=DP)

               ! Invert complex symmetric indefinite matrix
               call wrappers_invert_sym_cmatrix(c_squarebuffer, ngwf_basis%num)

               ! premultiply by inv_s to remove the extra factor of S
               ! And put back into original sparsity pattern of H inv_s
               squarebuffer = REAL(c_squarebuffer)
               call sparse_convert(spam_squarebuffer, squarebuffer)
               call sparse_product(sh,inverse_overlap,spam_squarebuffer)
               squarebuffer = AIMAG(c_squarebuffer)
               call sparse_convert(spam_squarebuffer, squarebuffer)
               call sparse_product(im_sh,inverse_overlap,spam_squarebuffer)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif  (pub_tddft_maxit_hotelling .gt. 0) then

               call tddft_hotelling_initial(sh_inv, sh, im_sh_inv, im_sh)
               call tddft_apply_hotelling(sh_inv, sh, &
                    im_sh_inv, im_sh, num_iter=pub_tddft_maxit_hotelling)

            elseif  (pub_tddft_maxit_hotelling .lt. 0) then

               write(stdout,*)'TDDFT_MAXIT HOTELLING LESS THAN 0. ONETEP stops'
               call comms_abort

            endif

            call sparse_scale(im_sh,-1.0_DP) ! - i 0.5 \Delta T inv_s H

            ! REAL PART
            call sparse_product(re_prop_adjoint(is),im_sh_inv,im_sh)
            call sparse_scale(re_prop_adjoint(is),-1.0_DP)
            !Re[ (1 + i 0.5 \Delta T inv_s H)^-1  (1 - i 0.5 \Delta T inv_s H) ]
            call sparse_axpy(re_prop_adjoint(is),sh_inv,1.0_DP)
            ! IMAGINARY PART
            call sparse_product(im_prop_adjoint(is),sh_inv,im_sh)
            !Im[ (1 + i 0.5 \Delta T inv_s H)^-1  (1 - i 0.5 \Delta T inv_s H) ]
            call sparse_axpy(im_prop_adjoint(is),im_sh_inv,1.0_DP)

         case (2,3)
            call sparse_transpose(re_prop_adjoint(is),re_prop(is))
            call sparse_transpose(im_prop_adjoint(is),im_prop(is))
            call sparse_scale(im_prop_adjoint(is),-1.0_DP)

         case default
            if (pub_on_root) write(stdout,*) 'Values of TDDFT_SPARSITY_LEVEL', &
                 ' may lie in the range 0 to 3 inclusive.',&
                 ' Setting to most conservative default 0.'
            call sparse_product(sh,inverse_overlap,hamiltonian(is))
            call sparse_copy(im_sh,sh)
            call sparse_scale(sh,0.0_DP,1.0_DP) ! 1
            call sparse_scale(im_sh,onefac) ! + i 0.5 \Delta T inv_s H

            ! Compute (1 + i 0.5 \Delta T inv_s H)^-1
            if (pub_tddft_maxit_hotelling .eq. 0) then

               ! Calculate exact inverse by inverting
               ! the symmetric (inv_s + i 0.5 \Delta T inv_s H inv_s)
               call sparse_product(spam_squarebuffer,sh,inverse_overlap)
               call sparse_convert(squarebuffer,spam_squarebuffer)
               c_squarebuffer = CMPLX(squarebuffer, 0.0_DP, kind=DP)
               call sparse_product(spam_squarebuffer,im_sh,inverse_overlap)
               call sparse_convert(squarebuffer,spam_squarebuffer)
               c_squarebuffer = c_squarebuffer + &
                    CMPLX(0.0_DP,squarebuffer, kind=DP)

               ! Invert complex symmetric indefinite matrix
               call wrappers_invert_sym_cmatrix(c_squarebuffer, ngwf_basis%num)

               ! premultiply by inv_s to remove the extra factor of S
               ! And put back into original sparsity pattern of H inv_s
               squarebuffer = REAL(c_squarebuffer)
               call sparse_convert(spam_squarebuffer, squarebuffer)
               call sparse_product(sh,inverse_overlap,spam_squarebuffer)
               squarebuffer = AIMAG(c_squarebuffer)
               call sparse_convert(spam_squarebuffer, squarebuffer)
               call sparse_product(im_sh,inverse_overlap,spam_squarebuffer)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif  (pub_tddft_maxit_hotelling .gt. 0) then

               call tddft_hotelling_initial(sh_inv, sh, im_sh_inv, im_sh)
               call tddft_apply_hotelling(sh_inv, sh, im_sh_inv, &
                    im_sh,num_iter=maxit_hotelling)

            elseif  (pub_tddft_maxit_hotelling .lt. 0) then

               write(stdout,*)'TDDFT_MAXIT HOTELLING LESS THAN 0. ONETEP stops'
               call comms_abort

            endif

            call sparse_scale(im_sh,-1.0_DP) ! - i 0.5 \Delta T inv_s H
            ! REAL PART
            call sparse_product(re_prop_adjoint(is),im_sh,im_sh_inv)
            call sparse_scale(re_prop_adjoint(is),-1.0_DP)
            !Re[  (1 - i 0.5 \Delta T inv_s H) (1 + i 0.5 \Delta T inv_s H)^-1 ]
            call sparse_axpy(re_prop_adjoint(is),sh_inv,1.0_DP)
            ! IMAGINARY PART
            call sparse_product(im_prop_adjoint(is),im_sh,sh_inv)
            !Im[  (1 - i 0.5 \Delta T inv_s H) (1 + i 0.5 \Delta T inv_s H)^-1 ]
            call sparse_axpy(im_prop_adjoint(is),im_sh_inv,1.0_DP)

         end select


      end do

      ! Clear memory
      select case (pub_tddft_sparsity_level)
      case (0,1)
         call sparse_destroy(im_sh_inv)
         call sparse_destroy(sh_inv)
         call sparse_destroy(im_sh)
         call sparse_destroy(sh)
      case (2,3)
      case default
         call sparse_destroy(im_sh_inv)
         call sparse_destroy(sh_inv)
         call sparse_destroy(im_sh)
         call sparse_destroy(sh)
      end select
      ! EXP(H)
      call sparse_destroy(im_hs_inv)
      call sparse_destroy(hs_inv)
      call sparse_destroy(im_hs)
      call sparse_destroy(hs)

      if (pub_tddft_maxit_hotelling .eq. 0) then
         call sparse_destroy(spam_squarebuffer)
         deallocate(c_squarebuffer,stat=ierr)
         call utils_dealloc_check('tddft_crank_nicholson_prop', &
              'c_squarebuffer',ierr)
         deallocate(squarebuffer,stat=ierr)
         call utils_dealloc_check('tddft_crank_nicholson_prop', &
              'squarebuffer',ierr)
      endif

      if (pub_on_root) write(stdout,*) 'Leaving tddft_crank_nicholson_prop_init'

      ! end timer
      call timer_clock("tddft_crank_nicholson_prop_init",2)

    end subroutine tddft_crank_nicholson_prop_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine tddft_generic_propagation(re_density, im_density,&
         re_prop, im_prop, re_prop_adjoint, im_prop_adjoint, im_in)

      !========================================================================!
      ! In this subroutine an  input density matrix is updated with a propagato!
      ! Matrix multiplication is carried out without truncation.               !
      !========================================================================!
      !========================================================================!
      ! Written by David D. O'Regan in March 2009.                             !
      !========================================================================!

      use rundat, only: pub_tddft_sparsity_level
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_transpose,&
           sparse_destroy, sparse_copy, sparse_product, sparse_scale
      use timer, only: timer_clock

      implicit none

      ! Real part of density matrix
      type(SPAM3), intent(inout) :: re_density(pub_cell%num_spins)
      ! Imaginary part
      type(SPAM3), intent(inout) :: im_density(pub_cell%num_spins)
      ! work with an imaginary part of input denskern
      logical, intent(in) :: im_in
      ! Real part of propagator from right:
      type(SPAM3), intent(in) :: re_prop(pub_cell%num_spins)
      ! Imaginary part of propagator from right:
      type(SPAM3), intent(in) :: im_prop(pub_cell%num_spins)
      ! Real part of propagator from left:
      type(SPAM3), intent(in) :: re_prop_adjoint(pub_cell%num_spins)
      ! Imaginary part of propagator from left:
      type(SPAM3), intent(in) :: im_prop_adjoint(pub_cell%num_spins)

      ! Local variables
      integer :: is      ! spin counter

      ! General purpose matrices of sparsity prop K prop
      type(SPAM3) :: re_dens_re_prop, re_dens_im_prop
      type(SPAM3) :: im_dens_re_prop, im_dens_im_prop
      type(SPAM3) :: aux_prop_dens_prop

      ! start timer
      call timer_clock("tddft_generic_propagation",1)

      if (pub_on_root) write(stdout,*) 'Entering tddft_generic_propagation'

      ! In all cases use biggest created structure for propagator
      ! and its adjoint
      call sparse_create(re_dens_re_prop,re_density(1),re_prop(1))
      call sparse_create(re_dens_im_prop,re_density(1),re_prop(1))
      if (im_in) then
         call sparse_create(im_dens_re_prop,re_density(1),re_prop(1))
         call sparse_create(im_dens_im_prop,re_density(1),re_prop(1))
      endif
      call sparse_create(aux_prop_dens_prop,re_prop_adjoint(1),re_dens_re_prop)

      ! SPARSITY PATTERNS FOR PROP_DENS_PROP when using Crank-Nicholson
      ! CASE (0) : FKF
      ! CASE (1) : SHSH K HSHS
      ! CASE (2) : HSHS K HSHS
      ! CASE (3) : HS K HS
      ! SPARSITY PATTERNS FOR PROP_DENS_PROP when using Runge-Kutta
      ! CASE (0) : FKF
      ! CASE (1) : HSHSHSHS K HSHSHSHS
      ! CASE (2) : HSHS K HSHS
      ! CASE (3) : HS K HS
      ! May wish to consider explicilty creating (SH)^n patterns in latter case

      do is=1,pub_cell%num_spins ! The main spin loop

         ! Re[P] Re[EXP(H)] and  Re[P] Im[EXP(H)]
         call sparse_product(re_dens_re_prop,re_density(is),re_prop(is))
         call sparse_product(re_dens_im_prop,re_density(is),im_prop(is))

         if (im_in) then
            ! Im[P] Re[EXP(H)] and  Im[P] Im[EXP(H)]
            call sparse_product(im_dens_re_prop,im_density(is),re_prop(is))
            call sparse_product(im_dens_im_prop,im_density(is),im_prop(is))
         endif

         ! Re[P_out] from Re[P_in]
         ! Re[EXP(-H)] Re[P] Re[EXP(H)]
         call sparse_product(re_density(is),re_prop_adjoint(is),re_dens_re_prop)
         call sparse_product(aux_prop_dens_prop,&
              im_prop_adjoint(is),re_dens_im_prop)
         ! -Im[EXP(-H)] Re[P] Im[EXP(H)]
         call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)

         ! Im[P_out] from Re[P_in]
         ! i Im[EXP(-H)] Re[P] Re[EXP(H)]
         call sparse_product(im_density(is),im_prop_adjoint(is),re_dens_re_prop)
         call sparse_product(aux_prop_dens_prop,&
              re_prop_adjoint(is),re_dens_im_prop)
         ! i Re[EXP(-H)] Re[P] Im[EXP(H)]
         call sparse_axpy(im_density(is),aux_prop_dens_prop,1.0_DP)

         if (im_in) then
            ! Re[P_out] from Im[P_in]
            call sparse_product(aux_prop_dens_prop,&
                 re_prop_adjoint(is),im_dens_im_prop)
            ! - Re[EXP(-H)] Im[P] Im[EXP(H)]
            call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)
            call sparse_product(aux_prop_dens_prop,&
                 im_prop_adjoint(is),im_dens_re_prop)
            ! - Im[EXP(-H)] Im[P] Re[EXP(H)]
            call sparse_axpy(re_density(is),aux_prop_dens_prop,-1.0_DP)

            ! Im[P_out] from Im[P_in]
            call sparse_product(aux_prop_dens_prop,&
                 re_prop_adjoint(is),im_dens_re_prop)
            ! i Re[EXP(-H)] Im[P] Re[EXP(H)]
            call sparse_axpy(im_density(is),aux_prop_dens_prop,1.0_DP)
            call sparse_product(aux_prop_dens_prop,&
                 im_prop_adjoint(is),im_dens_im_prop)
            ! -i Im[EXP(-H)] Im[P] Im[EXP(H)]
            call sparse_axpy(im_density(is),aux_prop_dens_prop,-1.0_DP)
         endif

      end do

      ! Clear memory
      ! EXP(-H) P EXP(H)
      call sparse_destroy(aux_prop_dens_prop)
      ! P EXP(H)
      if (im_in) then
         call sparse_destroy(im_dens_im_prop)
         call sparse_destroy(im_dens_re_prop)
      endif
      call sparse_destroy(re_dens_im_prop)
      call sparse_destroy(re_dens_re_prop)

      if (pub_on_root) write(stdout,*) 'Leaving tddft_generic_propagation'

      ! end timer
      call timer_clock("tddft_generic_propagation",2)

    end subroutine tddft_generic_propagation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine slow_sine_transform(steps, steplength_t, &
         time_domain, frequency_domain)

      use constants, only: DP, PI

      implicit none

      integer, intent(in) :: steps
      real(kind=DP), intent(in) :: steplength_t
      complex(kind=DP), intent(in) :: time_domain(steps)
      complex(kind=DP), intent(out) :: frequency_domain(steps) !(2*steps)+1)

      ! Internal variables
      integer :: count_w, count_t
      real(kind=DP) :: time
      real(kind=DP) :: frequency
      real(kind=DP) :: steplength_w

      steplength_w = PI / steplength_t
      frequency_domain = (0.0_DP,0.0_DP)

      do count_w = 1, steps ! should strictly be (2*steps)+1
         frequency = steplength_w * (count_w-1)
         do count_t = 1, steps
            time = steplength_t * (count_t-1)
            frequency_domain = time_domain * &
                 CMPLX(SIN( frequency * time ),0.0_DP,kind=DP)
         enddo
      enddo

    end subroutine slow_sine_transform

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine tddft_observables(polarisability, polarisability_t, &
         timestep_number, finesse_fac, proptime)!, electrons)

      !========================================================================!
      ! This subroutine accepts a frequency-dependent polarisability alpha(w)  !
      ! discretised on a grid of n points separated by deltaw = (2 pi / T).    !
      !                                                                        !
      ! With this we calculate:  the Dipole Strength Function                  !
      ! S(w) = (2 m / pi e^2 hbar) w Im [alpha (w)]                            !
      !                                                                        !
      ! The photoabsorption cross-section sigma(w) = (4 pi / c) w Im[alpha (w)]!
      !                                                                        !
      ! The static polarisability alpha(0) = (2/pi) \int dw (Im[alpha(w)]/w)   !
      !                                                                        !
      ! The Thomas-Reiche-Kuhn f-sum \int dE S(E) = \sum_i f_i = N             !
      !                                                                        !
      ! The linear dielectric function epsilon(w) = 1 + 4 pi N alpha (w)       !
      !                                                                        !
      ! The conductivity sigma(w) = w Im[ epsilon (w) ]                        !
      !                                                                        !
      ! The refraction index given by the Lorentz-Lorenz relation              !
      ! (4 pi / 3) N alpha(w) = (n^2(w) - 1) / (n^2(w) + 2)                    !
      !                                                                        !
      !========================================================================!
      !========================================================================!
      ! Written by David D. O'Regan in March 2009.                             !
      !========================================================================!

      use timer, only: timer_clock
      use constants, only: DP, PI, stdout
      use utils, only: utils_unit
      use comms, only: comms_abort, comms_barrier, pub_on_root
      use rundat, only: pub_rootname

      implicit none

      ! Polarisability alpha(w)
      integer, intent(in) :: timestep_number
      ! Oversampling factor of frequency_dependent_spectum
      integer, intent(in) :: finesse_fac
      complex(kind=DP), intent(in) :: polarisability(&
           &finesse_fac*2*timestep_number,3)
      ! Polarisability in time Domain
      complex(kind=DP), intent(in) :: polarisability_t(&
           &finesse_fac*2*timestep_number,3)
      real(kind=DP), intent(in) :: proptime
      !integer, intent(in) :: electrons(pub_cell%num_spins)

      real(kind=DP) :: deltaw
      real(kind=DP), parameter :: twouponpi = 0.5_DP / PI
      real(kind=DP), allocatable :: dipolestrength(:,:)
      real(kind=DP) :: frequency, fsum, density, staticpolarisability(4)
      integer :: step, dirn, point
      integer :: n ! Frequency steps to sample
      ! Refraction index and dielectric function including static values
      complex(kind=DP) :: fac, refraction((finesse_fac*timestep_number)+1)
      complex(kind=DP) :: dielectric((finesse_fac*timestep_number)+1)
      complex(kind=DP) :: invdielectric((finesse_fac*timestep_number)+1)
      complex(kind=DP) :: staticrefraction, staticdielectric
      complex(kind=DP) :: staticinvdielectric
      !complex(kind=DP) :: fourpidensity

      ! For writing output
      character(len=256) :: output_file  ! file names
      integer :: output_unit ! fortran output unit number

      ! start timer
      call timer_clock("tddft_observables",1)

      if (pub_on_root) write(stdout,*) 'Entering tddft_observables'

      n = finesse_fac * timestep_number

      allocate(dipolestrength(n+1,4),stat=ierr)
      call utils_alloc_check('tddft_observables','dipolestrength',ierr)

      deltaw = 2.0_DP * PI / (proptime * REAL(finesse_fac,kind=DP) )
      !? PI/proptime?
      fsum = 0.0_DP
      staticpolarisability(:) = 0.0_DP
      density = 1.0_DP ! molecules per unit volume, needs addressing
      ! Extra factors of 1/3 normalise the trace
      !fourpidensity = CMPLX(4.0_DP * PI * density / 3.0_DP,0.0_DP,kind=DP)
      fac = CMPLX(4.0 * PI * density / 9.0_DP ,0.0_DP,kind=DP)

      ! ddor: Note dipole_w(1) = \alpha(w=0)
      !            dipole_w(n+1) = \alpha(w_max)

      ! Loop over sampled frequencies
      do step=1,(n+1)
         frequency = deltaw * REAL(step-1,kind=DP)
         ! Loop over cartesian directions
         do dirn=1,3

            ! Dipole strength function
            dipolestrength(step,dirn) = AIMAG(polarisability(step,dirn)) &
                 &* frequency
            ! Static polarisability

            if (step .ne. 1) then
               staticpolarisability(dirn) =  staticpolarisability(dirn) + &
                    (  AIMAG(polarisability(step,dirn)) / REAL(step-1,kind=DP) )
            else
               staticpolarisability(dirn) = 0.0_DP
            endif

         enddo

         ! Put trace of Dipole strength function into 4th component
         dipolestrength(step,4) = ( dipolestrength(step,1) + &
              dipolestrength(step,2) + &
              dipolestrength(step,3) ) / (3.0_DP,0.0_DP)

         fsum = fsum + dipolestrength(step,4)

         ! Dielectric stores the trace of the polarisability
         dielectric(step) =  polarisability(step,1) + &
              polarisability(step,2) + polarisability(step,3)

         ! N^2 = (1 + 2 fac alpha) / (1 - fac alpha) by Lorentz-Lorenz
         refraction(step) = ( CMPLX(1.0_DP, 0.0_DP,kind=DP) + &
              ( CMPLX(2.0_DP, 0.0_DP,kind=DP) * &
              fac * dielectric(step)  ) ) / &
              ( CMPLX(1.0_DP, 0.0_DP,kind=DP) - ( fac * dielectric(step)  ) )

         refraction(step) = SQRT( refraction(step) )

         dielectric(step) = CMPLX(1.0_DP,0.0_DP,kind=DP) + &
              (fac * dielectric(step) ) ! Linear dielectric function
         invdielectric(step) = &
              dielectric(step) ** (-1.0_DP) ! Linear response function

      enddo ! Loop over frequency invervals

      ! Put trace of static polarisability into 4th component
      staticpolarisability(4) = &
           ( staticpolarisability(1) + staticpolarisability(2) + &
           staticpolarisability(3) ) / 3.0_DP

      ! Correctly normalise the results
      staticpolarisability = staticpolarisability * twouponpi
      dipolestrength =  dipolestrength * twouponpi ! Hartree^(-1) units
      !* 100.0_DP / REAL(SUM(electrons(:)))
      fsum = fsum * deltaw * twouponpi

      ! Use static polarisability to approximate static dielectric function and
      ! refractive index
      staticrefraction =  ( CMPLX(1.0_DP, 0.0_DP,kind=DP) + &
           ( CMPLX(2.0_DP, 0.0_DP,kind=DP) &
           * fac * staticpolarisability(4)  ) ) / &
           ( CMPLX(1.0_DP, 0.0_DP,kind=DP) - &
           ( fac * staticpolarisability(4)  ) )
      staticrefraction = SQRT( staticrefraction )
      ! Static Linear dielectric function
      staticdielectric = CMPLX(1.0_DP,0.0_DP,kind=DP) + &
           ( fac * staticpolarisability(4) )
      staticinvdielectric = staticdielectric ** (-1.0_DP)

      ! Having calculated the f-sum, multiply dipolestrength by
      ! ( 2 pi^2 e^2 hbar ) / ( c m ) in order to obtain the
      ! photoabsorption cross-section sigma(w)
      dipolestrength = dipolestrength * 2.0_DP * PI * PI / 137.03599907_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Print out some observables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (pub_on_root) then
         !write(stdout,*) 'frequency, Real part of &
         !     &frequency-dependent dielectric function'
         !do point = 1,(n+1)
         !   write(stdout,*) REAL(point-1,kind=DP)*deltaw,REAL(dielectric(point))
         !enddo
         !write(stdout,*) 'frequency, Imaginary part of &
         !     &frequency-dependent dielectric function'
         !do point = 1,(n+1)
         !   write(stdout,*) REAL(point-1,kind=DP)*deltaw,&
         !        &AIMAG(dielectric(point))
         !enddo
         !write(stdout,*) 'frequency, Real part of inverse &
         !     &frequency-dependent dielectric function'
         !do point = 1,(n+1)
         !   write(stdout,*) REAL(point-1,kind=DP)*deltaw,&
         !        &REAL(invdielectric(point))
         !enddo
         !write(stdout,*) 'frequency, Imaginary part of inverse &
         !     &frequency-dependent dielectric function'
         !do point = 1,(n+1)
         !   write(stdout,*) REAL(point-1,kind=DP)*deltaw,&
         !        &AIMAG(invdielectric(point))
         !enddo
         !write(stdout,*) 'frequency, Real part of &
         !     &frequency-dependent refractive index'
         !do point = 1,(n+1)
         !   write(stdout,*) REAL(point-1,kind=DP)*deltaw,&
         !        &REAL(refraction(point))
         !enddo
         !write(stdout,*) 'frequency, Imaginary part of &
         !     &frequency-dependent refractive index'
         !do point = 1,(n+1)
         !   write(stdout,*) REAL(point-1,kind=DP)*deltaw,&
         !        &AIMAG(refraction(point))
         !enddo
         !do dirn = 1,3
         !   write(stdout,*) 'frequency, Photoabsorption &
         !        &cross-section in direction ',dirn
         !   do point = 1,(n+1)
         !      write(stdout,*) REAL(point-1,kind=DP)*deltaw,&
         !           &dipolestrength(point,dirn)
         !   enddo
         !enddo
         !write(stdout,*) 'frequency, Trace of photoabsorption cross-section'
         !do point = 1,(n+1)
         !   write(stdout,*) REAL(point-1,kind=DP)*deltaw,dipolestrength(point,4)
         !enddo
         write(stdout,*) 'Static polarisability in each direction and its trace'
         do point = 1,4
            write(stdout,*) staticpolarisability(point)
         enddo
         write(stdout,*) 'Static dielectric function', staticdielectric
         write(stdout,*) 'Static inverse dielectric function', &
              &staticinvdielectric
         write(stdout,*) 'Static refractive index', staticrefraction
      endif
      if (pub_on_root) write(stdout,*) 'The charge of the Thomas-Reiche-Kuhn &
           & f-sum rule which has been recovered is: ',fsum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   The remainder of this subroutine based on  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   properties_dos_gp by C.-K. Skylaris        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (pub_on_root) then

         if (pub_on_root) write(stdout,'(a)')' '

         write(stdout,'(a)')&
              '===================== Photoabsorption cross-section calculation &
              &====================='

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !write(output_file,*)trim(pub_rootname)//'_RDF.txt'

         !output_file = adjustl(output_file)

         !cks: print output warning
         !write(stdout,'(3a)',advance ='no') &
         !     'Writing Real Dielectric Function "', trim(output_file),'" ...'

         !cks: get a unit number that is free
         !output_unit = utils_unit()

         !open(unit=output_unit, form="formatted" ,file=trim(output_file), &
         !     action="write", err=90 )

         !cks: write first line
         !write(output_unit,'(a)',err =100)'#  Energy (Ha) |  &
         !     &Real part of linear Dielectric function (a.u.)'

         !do point = 1,(n+1)

         !   write(output_unit,'(f20.6,2x,f30.26)',err =100) &
         !        &deltaw * REAL(point-1,kind=DP), REAL(dielectric(point))

         !enddo

         !close(unit =output_unit, err=110)

         ! cks: notify of end of output
         !write(stdout,*)' done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(output_file,*)trim(pub_rootname)//'_IDF.txt'

         output_file = adjustl(output_file)

         !cks: print output warning
         write(stdout,'(3a)',advance ='no') &
              'Writing Imaginary Dielectric Function "', &
              &trim(output_file),'" ...'

         !cks: get a unit number that is free
         output_unit = utils_unit()

         open(unit=output_unit, form="formatted" ,file=trim(output_file), &
              action="write", err=90 )

         !cks: write first line
         write(output_unit,'(a)',err =100)'#  Energy (Ha) |  &
              &Imaginary part of linear Dielectric function (a.u.)'

         do point = 1,(n+1)

            write(output_unit,'(f20.6,2x,f30.26)',err =100) &
                 &deltaw * REAL(point-1,kind=DP), AIMAG(dielectric(point))

         enddo

         close(unit =output_unit, err=110)

         ! cks: notify of end of output
         write(stdout,*)' done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(output_file,*)trim(pub_rootname)//'_RRI.txt'

         output_file = adjustl(output_file)

         !cks: print output warning
         write(stdout,'(3a)',advance ='no') &
              'Writing Real Refractive Index "', trim(output_file),'" ...'

         !cks: get a unit number that is free
         output_unit = utils_unit()

         open(unit=output_unit, form="formatted" ,file=trim(output_file), &
              action="write", err=90 )

         !cks: write first line
         write(output_unit,'(a)',err =100)'#  Energy (Ha) |  &
              &Real part of Refractive Index (a.u.)'

         do point = 1,(n+1)

            write(output_unit,'(f20.6,2x,f30.26)',err =100) &
                 &deltaw * REAL(point-1,kind=DP), REAL(refraction(point))

         enddo

         close(unit =output_unit, err=110)

         ! cks: notify of end of output
         write(stdout,*)' done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(output_file,*)trim(pub_rootname)//'_IRI.txt'

         output_file = adjustl(output_file)

         !cks: print output warning
         write(stdout,'(3a)',advance ='no') &
              'Writing Imaginary Refractive Index "', trim(output_file),'" ...'

         !cks: get a unit number that is free
         output_unit = utils_unit()

         open(unit=output_unit, form="formatted" ,file=trim(output_file), &
              action="write", err=90 )

         !cks: write first line
         write(output_unit,'(a)',err =100)'#  Energy (Ha) |  &
              &Imaginary part of Refractive Index (a.u.)'

         do point = 1,(n+1)

            write(output_unit,'(f20.6,2x,f30.26)',err =100) &
                 &deltaw * REAL(point-1,kind=DP), AIMAG(refraction(point))

         enddo

         close(unit =output_unit, err=110)

         ! cks: notify of end of output
         write(stdout,*)' done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(output_file,*)trim(pub_rootname)//'_polarisability.txt'

         output_file = adjustl(output_file)

         !cks: print output warning
         write(stdout,'(3a)',advance ='no') &
              'Writing Trace of polarisability in time "', &
              &trim(output_file),'" ...'

         !cks: get a unit number that is free
         output_unit = utils_unit()

         open(unit=output_unit, form="formatted" ,file=trim(output_file), &
              action="write", err=90 )

         !cks: write first line
         write(output_unit,'(a)',err =100)'#  Energy (Ha) |  &
              &Isotropic Time-dependent polarisability (a.u.)'

         do point = 1,(2*n)

            write(output_unit,'(f20.6,2x,f30.26)',err =100) &
                 &deltaw * REAL(point-1,kind=DP), &
                 &(( REAL(polarisability_t(point,1)) + &
                 & REAL(polarisability_t(point,2)) + &
                 & REAL(polarisability_t(point,3)) ) / 3.0_DP )

         enddo

         close(unit =output_unit, err=110)

         ! cks: notify of end of output
         write(stdout,*)' done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(output_file,*)trim(pub_rootname)//'_PCS.txt'

         output_file = adjustl(output_file)

         !cks: print output warning
         write(stdout,'(3a)',advance ='no') &
              'Writing "', trim(output_file),'" ...'

         !cks: get a unit number that is free
         output_unit = utils_unit()

         open(unit=output_unit, form="formatted" ,file=trim(output_file), &
              action="write", err=90 )

         !cks: write first line
         write(output_unit,'(a)',err =100)'#  Energy (Ha) |  &
              &Photoabsorption cross-section (a.u.)'

         do point = 1,(n+1)

            write(output_unit,'(f20.6,2x,f30.26)',err =100) &
                 &deltaw * REAL(point-1,kind=DP), dipolestrength(point,4)

         enddo

         close(unit =output_unit, err=110)

         ! cks: notify of end of output
         write(stdout,*)' done'

         write(stdout,'(a)')'================================&
              &================================================'

      endif
      call comms_barrier

      return


90    write(stdout,*)'Problem openning file in tddft_observables. ONETEP stops'
      call comms_abort
100   write(stdout,*)'Problem writing to file in tddft_observables. &
           &ONETEP stops'
      call comms_abort
110   write(stdout,*)'Problem closing file in tddft_observables. ONETEP stops'
      call comms_abort

      deallocate(dipolestrength,stat=ierr)
      call utils_dealloc_check('tddft_observables','dipolestrength',ierr)

      if (pub_on_root) write(stdout,*) 'Leaving tddft_observables'

      ! end timer
      call timer_clock("tddft_observables",2)

    end subroutine tddft_observables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine tddft_idempotency(step, realpart, overlap, &
         imaginarypart, original, kickswitch)

      !========================================================================!
      ! Calculates the idempotency and occupancy error for a given density     !
      ! kernel and overlap matrix.                                             !
      ! When also passed original (real only) density matrix it calulates the  !
      ! difference and projects out the part of the difference which violates  !
      ! preservation of idempotency to 1st order as                            !
      ! Pnew =Pold P0 + P0 Pold - 2 P0 Pold P0                                 !
      ! Thus, only occupied-virtual and virtual-occupied transitions allowed.  !
      !========================================================================!
      ! Written by David D. O'Regan in March 2009.                             !
      !========================================================================!

      use timer, only: timer_clock
      use constants, only: DP, stdout
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_product,&
           sparse_destroy, sparse_copy, sparse_scale, &
           sparse_trace, sparse_rms_element

      implicit none

      integer, intent(in) :: step ! timestep we are on
      ! Re[Denskern matrix]
      type(SPAM3), intent(inout) :: realpart(pub_cell%num_spins)
      ! Overlap matrix
      type(SPAM3), intent(in) :: overlap
      ! Im[Denskern matrix]
      type(SPAM3), optional, intent(inout) :: imaginarypart(pub_cell%num_spins)
      ! DFT Denskern matrix
      type(SPAM3), optional, intent(in) :: original(pub_cell%num_spins)
      ! Turns off attempt to purify DFT Denskern
      logical, optional, intent(in) :: kickswitch

      type(SPAM3) :: re_denskern_overlap, im_denskern_overlap
      real(kind=DP) :: idempotency_error, im_idempotency_error, nocc
      real(kind=DP) :: trace_re_s_re_s, trace_im_s_im_s
      ! For enforcing idempotency
      ! The DFT density kernel is real
      type(SPAM3) :: delta_denskern, KSK, KSKS, KSKSK

      ! start timer
      call timer_clock("tddft_idempotency",1)

      call sparse_create(re_denskern_overlap,realpart(1),overlap)
      if (present(imaginarypart)) then
         call sparse_create(im_denskern_overlap,imaginarypart(1),overlap)
      endif
      if (present(original)) then
         call sparse_create(delta_denskern,realpart(1))
         call sparse_create(KSK,re_denskern_overlap,realpart(1))
         call sparse_create(KSKS,KSK,overlap)
         call sparse_create(KSKSK,KSKS,realpart(1))
      endif

      nocc = 0.0_DP
      idempotency_error = 0.0_DP
      im_idempotency_error = 0.0_DP
      trace_re_s_re_s = 0.0_DP
      trace_im_s_im_s = 0.0_DP

      ! OUTPUT THE OCCUPANCY AND IDEMPOTENCY ERROR

      im1:if (present(imaginarypart)) then
         do is=1,pub_cell%num_spins
            call sparse_product(re_denskern_overlap,realpart(is),overlap)
            call sparse_product(im_denskern_overlap,imaginarypart(is),overlap)
            nocc = nocc + sparse_trace(re_denskern_overlap)
            trace_re_s_re_s = trace_re_s_re_s + &
                 sparse_trace(re_denskern_overlap,re_denskern_overlap)
            trace_im_s_im_s = trace_im_s_im_s + &
                 sparse_trace(im_denskern_overlap,im_denskern_overlap)
            idempotency_error = idempotency_error + &
                 trace_re_s_re_s - trace_im_s_im_s
            im_idempotency_error =  im_idempotency_error + &
                 sparse_trace(im_denskern_overlap,re_denskern_overlap) + &
                 sparse_trace(re_denskern_overlap,im_denskern_overlap) - &
                 sparse_trace(im_denskern_overlap)
         end do
         idempotency_error =  (idempotency_error - nocc)/nocc
         im_idempotency_error =  im_idempotency_error / nocc
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &Nocc ON TIMESTEP ',step,' IS ',nocc
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &Tr[ReK S ReK S] ON TIMESTEP ',step,' IS ',  &
              & trace_re_s_re_s
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &Tr[ImK S ImK S] ON TIMESTEP ',step,' IS ',  &
              & trace_im_s_im_s
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &(Re[Tr[\rho^2 - \rho]])/(Tr[\rho]) ON TIMESTEP ',&
              &step,' IS ',idempotency_error
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &(Im[Tr[\rho^2 - \rho]])/(Tr[\rho]) ON TIMESTEP ',&
              &step,' IS ',im_idempotency_error
         !idempotency_error = SQRT( (idempotency_error**2.0_DP) +&
         !     (im_idempotency_error**2.0_DP))
         !if (pub_on_root) write(stdout,*) 'tddft_calculate: &
         !     &(|Tr[\rho^2 - \rho]|)/(Tr[\rho]) ON TIMESTEP ',&
         !     &step,' IS ',idempotency_error
      else im1
         do is=1,pub_cell%num_spins
            call sparse_product(re_denskern_overlap,realpart(is),overlap)
            if (pub_on_root) write(stdout,*) &
                 &'sparse_trace(re_denskern_overlap)',&
                 sparse_trace(re_denskern_overlap)
            nocc = nocc + sparse_trace(re_denskern_overlap)
            if (pub_on_root) write(stdout,*) 'nocc',nocc
            idempotency_error =  idempotency_error + &
                 sparse_trace(re_denskern_overlap,re_denskern_overlap)
            if (pub_on_root) write(stdout,*) 'idempotency_error',&
                 &idempotency_error
         end do
         idempotency_error =  (idempotency_error - nocc)/nocc
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &Nocc ON TIMESTEP ',step,' IS ',nocc
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &(Tr[\rho^2 - \rho])/(Tr[\rho]) ON TIMESTEP ',&
              &step,' IS ',idempotency_error
      endif im1

      ! If handed the DFT density kernel then use it to
      ! project out the part of the
      ! first-order change to the density matrix which
      ! does not respect idempotency
      pure:if (present(original)) then

         ! No point in trying to purify the real part of the DFT density kernel
         if (  (.not.present(kickswitch)) .or. &
              &( (present(kickswitch)) .and. (kickswitch) ) ) then

            if (pub_on_root) write(stdout,*) &
                 &'Purifying real first-order change in density matrix'

            do is=1,pub_cell%num_spins

               !See how much the density changed since time t=0
               call sparse_copy(delta_denskern,realpart(is))
               call sparse_axpy(delta_denskern,original(is),-1.0_DP)

               ! K0 S
               call sparse_product(re_denskern_overlap,original(is),overlap)
               ! K0 S K1
               call sparse_product(KSK,re_denskern_overlap,delta_denskern)
               ! K0 S K1 S
               call sparse_product(KSKS,KSK,overlap)
               ! K0 S K1 S K0
               call sparse_product(KSKSK,KSKS,original(is))
               ! - 2 K0 S K1 S K0
               call sparse_scale(KSKSK,-2.0_DP)

               ! Add K0 S K1 into the buffer
               call sparse_axpy(KSKSK,KSK,1.0_DP)
               ! K1 S
               call sparse_product(re_denskern_overlap,delta_denskern, overlap)
               ! K1 S K0
               call sparse_product(KSK,re_denskern_overlap,original(is))
               ! Add K1 S K0 into the buffer
               call sparse_axpy(KSKSK,KSK,1.0_DP)
               ! Add K0 into the buffer
               call sparse_axpy(KSKSK,original(is),1.0_DP)
               ! The purified density matrix
               call sparse_copy(realpart(is),KSKSK)

            enddo

         endif

         im2: if (present(imaginarypart)) then

            if (pub_on_root) write(stdout,*) &
                 &'Purifying imaginary first-order change in density matrix'

            do is=1,pub_cell%num_spins

               !See how much the density changed since time t=0
               call sparse_copy(delta_denskern,imaginarypart(is))

               ! K0 S
               call sparse_product(im_denskern_overlap,original(is),overlap)
               ! K0 S K1
               call sparse_product(KSK,im_denskern_overlap,delta_denskern)
               ! K0 S K1 S
               call sparse_product(KSKS,KSK,overlap)
               ! K0 S K1 S K0
               call sparse_product(KSKSK,KSKS,original(is))
               ! - 2 K0 S K1 S K0
               call sparse_scale(KSKSK,-2.0_DP)

               ! Add K0 S K1 into the buffer
               call sparse_axpy(KSKSK,KSK,1.0_DP)
               ! K1 S
               call sparse_product(im_denskern_overlap,delta_denskern, overlap)
               ! K1 S K0
               call sparse_product(KSK,im_denskern_overlap,original(is))
               ! Add K1 S K0 into the buffer
               call sparse_axpy(KSKSK,KSK,1.0_DP)
               ! The purified density matrix
               call sparse_copy(imaginarypart(is),KSKSK)

            enddo

         endif im2

         nocc = 0.0_DP
         idempotency_error = 0.0_DP
         im_idempotency_error = 0.0_DP
         trace_re_s_re_s = 0.0_DP
         trace_im_s_im_s = 0.0_DP

         ! OUTPUT THE NEW OCCUPANCY AND IDEMPOTENCY ERROR

         im3:if (present(imaginarypart)) then
            do is=1,pub_cell%num_spins
               call sparse_product(re_denskern_overlap,realpart(is),overlap)
               call sparse_product(im_denskern_overlap,&
                    imaginarypart(is),overlap)
               nocc = nocc + sparse_trace(re_denskern_overlap)
               trace_re_s_re_s = trace_re_s_re_s + &
                    sparse_trace(re_denskern_overlap,re_denskern_overlap)
               trace_im_s_im_s = trace_im_s_im_s + &
                    sparse_trace(im_denskern_overlap,im_denskern_overlap)
               idempotency_error =  idempotency_error + &
                    trace_re_s_re_s - trace_im_s_im_s
               im_idempotency_error =  im_idempotency_error + &
                    sparse_trace(im_denskern_overlap,re_denskern_overlap) + &
                    sparse_trace(re_denskern_overlap,im_denskern_overlap) - &
                    sparse_trace(im_denskern_overlap)
            end do
            idempotency_error =  (idempotency_error - nocc)/nocc
            im_idempotency_error =  im_idempotency_error / nocc
            if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                 &purified Nocc ON TIMESTEP ',step,' IS ',nocc
            if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                 &purified (Re[Tr[\rho^2 - \rho]])/(Tr[\rho]) ON TIMESTEP ',&
                 step,' IS ',idempotency_error
            if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                 &purified (Im[Tr[\rho^2 - \rho]])/(Tr[\rho]) ON TIMESTEP ',&
                 step,' IS ',im_idempotency_error
            !idempotency_error = SQRT( (idempotency_error**2.0_DP) +&
            !     (im_idempotency_error**2.0_DP))
            !if (pub_on_root) write(stdout,*) 'tddft_calculate: &
            !     &purified (|Tr[\rho^2 - \rho]|)/(Tr[\rho]) ON TIMESTEP ',&
            !     step,' IS ',idempotency_error
         else im3
            do is=1,pub_cell%num_spins
               call sparse_product(re_denskern_overlap,realpart(is),overlap)
               nocc = nocc + sparse_trace(re_denskern_overlap)
               idempotency_error =  idempotency_error + &
                    sparse_trace(re_denskern_overlap,re_denskern_overlap)
            end do
            idempotency_error =  (idempotency_error - nocc)/nocc
            if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                 &purified Nocc ON TIMESTEP ',step,' IS ',nocc
            if (pub_on_root) write(stdout,*) 'tddft_calculate: &
                 &purified (Tr[\rho^2 - \rho])/(Tr[\rho]) ON TIMESTEP ',&
                 step,' IS ',idempotency_error
         endif im3


      endif pure

      ! 3rd order purification of density kernel
      !    call kernel_purify(pur_denskern, denskern, &
      !    overlap, inv_overlap, n_occ)
      !  RESCALE K
      !    call kernel_rescale(aux,overlap,n_occ,can_rescale_ks=.true.)

      if (present(original)) then
         call sparse_destroy(KSKSK)
         call sparse_destroy(KSKS)
         call sparse_destroy(KSK)
         call sparse_destroy(delta_denskern)
      endif
      if (present(imaginarypart)) then
         call sparse_destroy(im_denskern_overlap)
      endif
      call sparse_destroy(re_denskern_overlap)

      ! end timer
      call timer_clock("tddft_idempotency",2)

    end subroutine tddft_idempotency

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine tddft_idempotency_test(step, realpart, overlap, &
         imaginarypart)

      !========================================================================!
      ! Calculates the idempotency and occupancy error for a given density     !
      ! kernel and overlap matrix.                                             !
      !========================================================================!
      ! Written by David D. O'Regan in March 2009.                             !
      !========================================================================!

      use timer, only: timer_clock
      use constants, only: DP, stdout
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create, sparse_axpy, sparse_product,&
           sparse_destroy, sparse_copy, sparse_scale, &
           sparse_trace, sparse_rms_element

      implicit none

      integer, intent(in) :: step ! timestep we are on
      ! Re[Denskern matrix]
      type(SPAM3), intent(in) :: realpart(pub_cell%num_spins)
      ! Overlap matrix
      type(SPAM3), intent(in) :: overlap
      ! Re[Denskern matrix]
      type(SPAM3), optional, intent(in) :: imaginarypart(pub_cell%num_spins)

      type(SPAM3) :: re_denskern_overlap, im_denskern_overlap
      real(kind=DP) :: idempotency_error, im_idempotency_error, nocc
      real(kind=DP) :: trace_re_s_re_s, trace_im_s_im_s

      ! start timer
      call timer_clock("tddft_idempotency_test",1)

      call sparse_create(re_denskern_overlap,realpart(1),overlap)
      if (present(imaginarypart)) then
         call sparse_create(im_denskern_overlap,imaginarypart(1),overlap)
      endif

      nocc = 0.0_DP
      idempotency_error = 0.0_DP
      im_idempotency_error = 0.0_DP
      trace_re_s_re_s = 0.0_DP
      trace_im_s_im_s = 0.0_DP

      ! OUTPUT THE OCCUPANCY AND IDEMPOTENCY ERROR

      im1:if (present(imaginarypart)) then
         do is=1,pub_cell%num_spins
            call sparse_product(re_denskern_overlap,realpart(is),overlap)
            call sparse_product(im_denskern_overlap,imaginarypart(is),overlap)
            nocc = nocc + sparse_trace(re_denskern_overlap)
            trace_re_s_re_s = trace_re_s_re_s + &
                 sparse_trace(re_denskern_overlap,re_denskern_overlap)
            trace_im_s_im_s = trace_im_s_im_s + &
                 sparse_trace(im_denskern_overlap,im_denskern_overlap)
            idempotency_error =  idempotency_error + &
                 trace_re_s_re_s - trace_im_s_im_s
            im_idempotency_error =  im_idempotency_error + &
                 sparse_trace(im_denskern_overlap,re_denskern_overlap) + &
                 sparse_trace(re_denskern_overlap,im_denskern_overlap) - &
                 sparse_trace(im_denskern_overlap)
         end do
         idempotency_error =  (idempotency_error - nocc)/nocc
         im_idempotency_error =  im_idempotency_error / nocc
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &Nocc ON TIMESTEP ',step,' IS ',nocc
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &Tr[ReK S ReK S] ON TIMESTEP ',step,' IS ',  &
              trace_re_s_re_s
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &Tr[ImK S ImK S] ON TIMESTEP ',step,' IS ',  &
              trace_im_s_im_s
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &(Re[Tr[\rho^2 - \rho]])/(Tr[\rho]) ON TIMESTEP ',&
              step,' IS ',idempotency_error
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &(Im[Tr[\rho^2 - \rho]])/(Tr[\rho]) ON TIMESTEP ',&
              step,' IS ',im_idempotency_error
         !idempotency_error = SQRT( (idempotency_error**2.0_DP) +&
         !     (im_idempotency_error**2.0_DP))
         !if (pub_on_root) write(stdout,*) 'tddft_calculate: &
         !     &(|Tr[\rho^2 - \rho]|)/(Tr[\rho]) ON TIMESTEP ',&
         !     step,' IS ',idempotency_error
      else im1
         do is=1,pub_cell%num_spins
            call sparse_product(re_denskern_overlap,realpart(is),overlap)
            if (pub_on_root) write(stdout,*) 'sparse_trace(realpart,overlap)',&
                 sparse_trace(realpart(is),overlap)
            nocc = nocc + sparse_trace(re_denskern_overlap)
            if (pub_on_root) write(stdout,*) 'nocc',nocc
            idempotency_error =  idempotency_error + &
                 sparse_trace(re_denskern_overlap,re_denskern_overlap)
            if (pub_on_root) write(stdout,*) 'idempotency_error',&
                 idempotency_error
         end do
         idempotency_error =  (idempotency_error - nocc)/nocc
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &Nocc ON TIMESTEP ',step,' IS ',nocc
         if (pub_on_root) write(stdout,*) 'tddft_calculate: &
              &(Tr[\rho^2 - \rho])/(Tr[\rho]) ON TIMESTEP ',&
              step,' IS ',idempotency_error
      endif im1

      if (present(imaginarypart)) then
         call sparse_destroy(im_denskern_overlap)
      endif
      call sparse_destroy(re_denskern_overlap)

      ! end timer
      call timer_clock("tddft_idempotency_test",2)

    end subroutine tddft_idempotency_test

  end subroutine tddft_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tddft_hotelling_initial(inv_s, overlap, im_inv_s, im_overlap)

    !============================================================!
    !  This subroutine creates an initial guess for the          !
    !  calculation of the inverse overlap matrix of NGWFs        !
    !  by Hotelling's method, based on the theory in the         !
    !  paper by T. Ozaki, Phys. Rev. B., vol 64, page 195110.    !
    !------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 31/7/2003   !
    ! for the ONETEP program.                                    !
    ! Modified to support SPAM 2 by Peter Haynes on 21/7/2004.   !
    ! Modified to accept an optional second matrix, intended for !
    ! imaginary components by David D. O'Regan on 7/5/2009.      !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009  !
    !============================================================!


    use comms, only: pub_on_root, comms_reduce, pub_my_node_id
    use constants, only: DP, stdout
    use sparse, only: SPAM3, sparse_scale, sparse_get_col, sparse_clr_col, &
         sparse_copy, sparse_first_elem_on_node, sparse_num_rows
    use simulation_cell, only : pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(SPAM3), intent(inout) :: inv_s
    type(SPAM3), intent(in)  :: overlap
    ! ddor: Optional 2nd matrices to invert
    type(SPAM3), optional, intent(inout) :: im_inv_s  ! output
    type(SPAM3), optional, intent(in) :: im_overlap ! input

    ! cks: local variables
    real(kind=DP) :: sigma, sum_col
    real(kind=DP), allocatable, dimension(:) :: olap_row
    ! ddor: Used when im_overlap is present
    real(kind=DP), allocatable, dimension(:) :: im_olap_row
    integer :: nn, row, col
    integer :: ierr

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: &
         &Entering tddft_hotelling_initial'
#endif

    ! Allocate workspace
    nn = sparse_num_rows(overlap)
    allocate(olap_row(nn),stat=ierr)
    call utils_alloc_check('tddft_hotelling_initial','olap_row',ierr)
    ! ddor: Allocate im_olap_row when im_overlap is present
    if(present(im_overlap)) then
       allocate(im_olap_row(nn),stat=ierr)
       call utils_alloc_check('tddft_hotelling_initial','im_olap_row',ierr)
    endif

    sigma = 0.0_DP
    olap_row = 0.0_DP
    ! ndmh: initialise imaginary column if necessary
    if (present(im_overlap)) then
       im_olap_row = 0.0_DP
    end if

    ! Loop over all cols of S on this processor
    do col=sparse_first_elem_on_node(pub_my_node_id,overlap,'C'), &
         sparse_first_elem_on_node(pub_my_node_id+1,overlap,'C') - 1

       ! Get column of S
       call sparse_get_col(olap_row,overlap,col)
       ! ddor: Get column of im_overlap if necessary
       if(present(im_overlap)) then
          call sparse_get_col(im_olap_row,im_overlap,col)
       endif

       sum_col = 0.0_DP
       ! ddor: include imaginary component of matrix if necessary
       if(present(im_overlap)) then
          do row=1,nn
             sum_col = sum_col + &
                  SQRT( (olap_row(row))**2.0_DP + (im_olap_row(row))**2.0_DP )
          end do
       else
         do row=1,nn
            sum_col = sum_col + abs(olap_row(row))
         end do
       end if
       sigma = max(sigma,sum_col)

       call sparse_clr_col(olap_row,overlap,col)
       ! ddor: Clear columns of im_olap_row if required
       if(present(im_overlap)) then
          call sparse_clr_col(im_olap_row,overlap,col)
       endif

    end do

    call comms_reduce('MAX',sigma)

    ! Deallocate workspace (ddor: deallocation for im_olap_row if required)
    if(present(im_overlap)) then
       deallocate(im_olap_row,stat=ierr)
       call utils_dealloc_check('tddft_hotelling_initial','im_olap_row',ierr)
    endif
    deallocate(olap_row,stat=ierr)
    call utils_dealloc_check('tddft_hotelling_initial','olap_row',ierr)

    ! Copy overlap into inv_s
    call sparse_copy(inv_s,overlap)
    ! ddor: Copy im_overlap into im_inv_s if required
    if(present(im_overlap)) then
       call sparse_copy(im_inv_s,im_overlap)
    endif

    if (sigma > epsilon(1.0_DP)) then
       sigma = 1.0_DP / (sigma*sigma)
    else
       sigma = 0.001_DP
       if (pub_on_root) write(stdout,'(a)') &
            'WARNING in tddft_hotelling_initial: zero overlap matrix'
    end if

    ! cks: scale by sigma
    call sparse_scale(inv_s,sigma)
    ! ddor: scale im_overlap by sigma if required
    if(present(im_overlap)) then
       call sparse_scale(im_inv_s,sigma)
    endif

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving tddft_hotelling_initial'
#endif

  end subroutine tddft_hotelling_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tddft_apply_hotelling(inv_s,overlap,im_inv_s,im_overlap,num_iter)

    !=============================================================!
    !  This subroutine applies successive quadratically           !
    !  convergent Hotelling iterations to improve an approximate  !
    !  inverse overlap matrix, based on the theory in the         !
    !  paper by T. Ozaki, Phys. Rev. B., vol 64, page 195110.     !
    !  The iterations continue (unless they exceed the            !
    !  maxit_hotelling input parameter) until convergence         !
    !  to machine precision is reached - but taking into          !
    !  account limitations arising from the truncation of the     !
    !  inverse overlap matrix.                                    !
    !-------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 31/7/2003    !
    ! for the ONETEP program.                                     !
    ! Modified to support SPAM 2 by Peter Haynes on 21/7/2004.    !
    ! Modified to test for convergence by Chris-Kriton Skylaris   !
    ! on 4/10/2004.                                               !
    ! Minor modifications for parallel SPAM 2 by Peter Haynes     !
    ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009   !
    !=============================================================!

    use comms, only: pub_on_root,comms_abort
    use constants, only: DP, VERBOSE, stdout
    use rundat, only: pub_output_detail, max_resid_hotelling, &
         pub_tddft_max_resid_hotelling
    use sparse, only: SPAM3, sparse_scale, sparse_product, sparse_create, &
         sparse_destroy, sparse_num_element, sparse_max_abs_element, &
         sparse_rms_element, sparse_copy, sparse_axpy
    use timer, only: timer_clock

    implicit none

    type(SPAM3), intent(inout) :: inv_s
    type(SPAM3), intent(in) :: overlap
    integer, optional, intent(in) :: num_iter
    ! ddor: Optional 2nd matrices to invert
    type(SPAM3), optional, intent(inout) :: im_inv_s
    type(SPAM3), optional, intent(in) :: im_overlap

    ! cks: local variables
    integer :: loc_num_iter
    integer :: iter
    logical :: quitearly  ! pdh: flag
    real(kind=DP) :: max_resid  ! maximum value of residual
    real(kind=DP) :: frob_norm
    real(kind=DP) :: previous_frob_norm
    type(SPAM3) :: sk,ktmp
    ! ddor: Additional matrices needed when imaginary part is present
    type(SPAM3) :: im_sk,im_ktmp,sk_buffer
    real(kind=DP) :: im_frob_norm ! ddor

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering tddft_apply_hotelling'
#endif

    ! cks: start timer
    call timer_clock('tddft_apply_hotelling',1)

    ! cks: initialisations
    previous_frob_norm = huge(1.0_DP)

    if (present(num_iter)) then
       loc_num_iter = num_iter
    else
       loc_num_iter = 1
    end if

    quitearly = .false.

    if (pub_output_detail == VERBOSE .and. pub_on_root) then
       write(stdout,'(a)')'============ Calculation of NGWF S^-1 using &
            &Hotelling algorithm ================ '
       write(stdout,'(a)')'         Iteration  |  Frobenius norm and Maximum &
            &abs value    (of I-S*S_n^-1)   '
    end if

    ! Allocate workspace
    call sparse_create(ktmp,inv_s)
    call sparse_create(sk,overlap,inv_s)
    ! ddor: Allocate workspace needed for imaginary part
    if(present(im_overlap)) then
       call sparse_create(im_ktmp,inv_s)
       call sparse_create(im_sk,overlap,inv_s)
       call sparse_create(sk_buffer,overlap,inv_s)
    endif

    ! pdh: Start of Hotelling iteration loop
    hotel_sparse_loop: do iter=1,loc_num_iter

       ! SK := S*M
       call sparse_product(sk,overlap,inv_s)
       ! ddor: if imaginary part is present
       if(present(im_overlap)) then
          call sparse_product(sk_buffer,im_overlap,im_inv_s)
          call sparse_axpy(sk,sk_buffer,-1.0_DP) !Re[SK] = ReS ReK - ImS ImK
          call sparse_product(im_sk,im_overlap,inv_s)
          call sparse_product(sk_buffer,overlap,im_inv_s)
          call sparse_axpy(im_sk,sk_buffer,1.0_DP) !Im[SK] = ImS ReK + ReS ImK
       endif

       ! cks: SK := I-S*M
       call sparse_scale(sk,-1.0_DP,1.0_DP)
       ! ddor: if imaginary part is present
       if(present(im_overlap)) then
          call sparse_scale(im_sk,-1.0_DP,0.0_DP)
          ! Re[SK]-> 1-Re[SK], Im[SK]-> -Im[SK]
       endif

       ! cks: Calculate Frobenius norm of SK -->0
       frob_norm = sparse_rms_element(sk) * sqrt(sparse_num_element(sk))
       ! ddor: include the imaginary contribution to the norm if required
       if(present(im_overlap)) then !ddor
          im_frob_norm = sparse_rms_element(im_sk) * &
               sqrt(sparse_num_element(im_sk))
          frob_norm = SQRT(frob_norm**2.0_DP + im_frob_norm**2.0_DP)
       endif

       ! cks: Maximum element of residual I-S*S_n^-1
       ! ddor: Use maximum element of either real or imaginary part if required
       if(present(im_overlap)) then
          max_resid = MAX(sparse_max_abs_element(sk),&
               sparse_max_abs_element(im_sk))
       else
          max_resid = sparse_max_abs_element(sk)
       endif


       ! cks: Test for convergence due to machine precision
       ! cks: or inverse overlap truncation
       if (pub_output_detail == VERBOSE .and. pub_on_root) then
          write(stdout,'(t12,i5,tr5,e16.8,tr3,e16.8)') iter, frob_norm, &
               max_resid
       end if

       ! ddor: If the imaginary part is present then use
       !        pub_tddft_max_resid_hotelling
       !       as a more stringent criterion on the residual
       if(present(im_overlap)) then
          if (frob_norm >= previous_frob_norm .or. ( &
               (max_resid <= max_resid_hotelling) .and. &
               (max_resid <= pub_tddft_max_resid_hotelling) ) ) then
             quitearly = .true.
             exit
          else
             previous_frob_norm = frob_norm
          end if
       else
          if (frob_norm >= previous_frob_norm .or. &
               max_resid <= max_resid_hotelling) then
             quitearly = .true.
             exit
          else
             previous_frob_norm = frob_norm
          endif
       end if

       ! cks: SK := I +SK := 2I -S*M
       call sparse_scale(sk,1.0_DP,1.0_DP)

       ! M := M*(2*I-S*M)
       call sparse_copy(ktmp,inv_s)
       call sparse_product(inv_s,ktmp,sk)
       ! ddor: Carry out matrix products for imaginary part
       if(present(im_overlap)) then !ddor
          call sparse_copy(im_ktmp,im_inv_s)
          call sparse_product(im_inv_s,im_ktmp,sk)
          call sparse_product(sk_buffer,im_ktmp,im_sk)
          call sparse_axpy(inv_s,sk_buffer,-1.0_DP)
          call sparse_product(sk_buffer,ktmp,im_sk)
          call sparse_axpy(im_inv_s,sk_buffer,1.0_DP)
       endif

    end do hotel_sparse_loop

    if (pub_on_root .and. .not. quitearly) write(stdout,*) &
         'WARNING: max Hotelling iterations exceeded. Last norm=',frob_norm

    ! Deallocate workspace
    ! ddor: Carry out deallocation for optional matrices
    if(present(im_overlap)) then
       call sparse_destroy(sk_buffer)
       call sparse_destroy(im_sk)
       call sparse_destroy(im_ktmp)
    endif
    call sparse_destroy(sk)
    call sparse_destroy(ktmp)

    if (pub_output_detail == VERBOSE .and. pub_on_root) then
       write(stdout,'(a)')'===================================&
            &============================================='
    endif

    ! cks: stop timer
    call timer_clock('tddft_apply_hotelling',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving tddft_apply_hotelling'
#endif

  end subroutine tddft_apply_hotelling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tddft
