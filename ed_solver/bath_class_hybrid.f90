MODULE bath_class_hybrid

   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   !$$ MAKE HYBRIDIZATION FROM BATH $$
   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   ! USE bath_class_vec
   ! USE minimization_wrapping
   ! use stringmanip
   use correl_class, only: correl_type
   use bath_class,   only: bath_type

   IMPLICIT NONE

   TYPE(correl_type),        PRIVATE, SAVE     :: hybrid2fit
   TYPE(bath_type),          PRIVATE, SAVE     :: batht
   integer,                  PRIVATE          :: icount_ = 0

   public :: bath2hybrid
   public :: hybrid2bath

contains

   function fill_mat_from_vector(n, vec)

      implicit none

      integer :: n, i, j, k
      real(8) :: vec(n*(n + 1)/2), fill_mat_from_vector(n, n)

      k = 0
      do i = 1, n
         do j = i, n
            k = k + 1
            fill_mat_from_vector(i, j) = vec(k)
            if(i /= j) fill_mat_from_vector(j, i) = fill_mat_from_vector(i, j)
         enddo
      enddo
      if(k /= size(vec))then
         write(*, *) 'ERROR fill_mat_from_vector'
         stop
      endif
   end function

   subroutine bath2hybrid(BATH, FREQTYPE, WRITE_HYBRID, short, cptvec, &
        cpt_build_matrix, Vweight)

      use correl_class,        only: glimpse_correl
      use bath_class,          only: bath_type, nambu_eb, nambu_vbc
      use genvar,              only: dbl, fermionic, iproc, retarded
      use globalvar_ed_solver, only: bath_nearest_hop, cpt_upper_bound, &
           cutoff_rvb, diag_bath, diag_v, epsilon_cpt, &
           fit_green_func_and_not_hybrid, fast_fit, fit_nw, fmos, &
           force_no_pairing, force_nupdn_basis, freeze_pole_lambda, &
           freeze_poles_delta, freeze_poles_delta_iter, &
           freeze_poles_delta_niter, frozen_poles, ncpt_approx, &
           ncpt_approx_tot, ncpt_chain_coup, ncpt_para, ncpt_tot, &
           supersc_state, t_cpt
      use impurity_class,      only: Eccc
      use masked_matrix_class, only: delete_masked_matrix, masked_matrix_type
      use matrix,              only: bande_mat, eigenvector_matrix, invmat, &
                                     matmul_x

      implicit none

      CHARACTER(LEN = 9), INTENT(IN) :: FREQTYPE
      LOGICAL, OPTIONAL, INTENT(IN)  :: WRITE_HYBRID, short
      TYPE(bath_type)           :: BATH
      TYPE(masked_matrix_type)  :: EbNambu, VbcNambu
      COMPLEX(DBL)              :: Eb((BATH%Nb + ncpt_approx_tot)*2, (BATH%Nb &
                                   + ncpt_approx_tot)*2)
      COMPLEX(DBL)              :: Ebm1((BATH%Nb + ncpt_approx_tot)*2, &
                                   (BATH%Nb + ncpt_approx_tot)*2)
      COMPLEX(DBL)              :: Vbc((BATH%Nb + ncpt_approx_tot)*2, &
                                   BATH%Nc*2), Vcb(BATH%Nc*2, (BATH%Nb + &
                                   ncpt_approx_tot)*2)
      REAL(DBL)                 :: ratio, eigenvalues(2*(BATH%Nb + &
                                   ncpt_approx_tot))
      COMPLEX(DBL), POINTER     :: hybrid(:,:,:) => NULL(), freq(:) => NULL()
      INTEGER                   :: iw, Nw, i, j, k, Nb, Ntot, Nc, bl, mm, ll, &
                                   mmstep
      LOGICAL                   :: write_, block, diagmat
      REAL(8), TARGET, OPTIONAL :: cptvec(ncpt_para*ncpt_tot)
      REAL(8), POINTER          :: cptup(:) => NULL(), cptdn(:) => NULL(), &
                                   cptbcs(:) => NULL()
      LOGICAL, OPTIONAL         :: cpt_build_matrix
      REAL(8), OPTIONAL         :: Vweight


      Nb = BATH%Nb
      Ntot = Nb + ncpt_approx_tot
      Nc = BATH%Nc

      write_ = .false.
      IF(PRESENT(WRITE_HYBRID)) write_ = WRITE_HYBRID

      CALL Nambu_Eb (EbNambu, BATH%Eb, BATH%Pb)
      CALL Nambu_Vbc(VbcNambu, BATH%Vbc, BATH%PVbc)

      if(ncpt_approx_tot == 0)then

         Eb  = -EbNambu%rc%mat
         Vbc =  VbcNambu%rc%mat
         Vcb =  TRANSPOSE(conjg(Vbc))

      else

         Vbc = 0.d0
         Vcb = 0.d0
         Eb = 0.d0


         Eb( 1: Nb, 1 :Nb) = EbNambu%rc%mat( 1: Nb, 1: Nb)
         Eb(Ntot + 1:Ntot + Nb, Ntot + 1:Ntot + Nb) = EbNambu%rc%mat(Nb + &
              1:2*Nb, Nb + 1:2*Nb)
         Vbc( 1 :Nb, 1:2*Nc) = VbcNambu%rc%mat( 1: Nb, 1:2*Nc)
         Vbc(Ntot + 1:Ntot + Nb, 1:2*Nc) = VbcNambu%rc%mat(Nb + 1:2*Nb, &
              1:2*Nc)
         Vcb                                 =   TRANSPOSE(conjg(Vbc))

         if(present(cptvec))then
            cptup => cptvec(           1:  ncpt_tot)
            if(supersc_state .and. .not.force_nupdn_basis .and. &
                 .not.force_no_pairing)then
               if(ncpt_para == 3)then
                  cptdn => cptvec(  ncpt_tot + 1:2*ncpt_tot)
                  cptbcs => cptvec(2*ncpt_tot + 1:3*ncpt_tot)
               else
                  cptdn => cptvec(           1:  ncpt_tot)
                  cptbcs => cptvec(  ncpt_tot + 1:2*ncpt_tot)
               endif
            else
               if(ncpt_para == 2)then
                  cptdn => cptvec(ncpt_tot + 1:2*ncpt_tot)
               else
                  cptdn => cptvec(         1:  ncpt_tot)
               endif
            endif

            do i = 1, ncpt_chain_coup*Nb
               if(abs(cptup(i)) > cpt_upper_bound) cptup(i) = &
                    cptup(i)/abs(cptup(i))*cpt_upper_bound
               if(abs(cptdn(i)) > cpt_upper_bound) cptdn(i) = &
                    cptdn(i)/abs(cptdn(i))*cpt_upper_bound
               if(supersc_state .and. .not.force_nupdn_basis .and. &
                    .not.force_no_pairing)then
                  if(abs(cptbcs(i)) > cpt_upper_bound) cptbcs(i) = &
                       cptbcs(i)/abs(cptbcs(i))*cpt_upper_bound
               endif
            enddo

            if(present(Vweight)) then
               Vweight = sum(abs(cptup( 1:ncpt_chain_coup*Nb ))**2) + &
                    sum(abs(cptdn( 1:ncpt_chain_coup*Nb ))**2)
               if(supersc_state .and. .not.force_nupdn_basis .and. &
                    .not.force_no_pairing)then
                  Vweight = Vweight + sum(abs(cptbcs( 1:ncpt_chain_coup*Nb &
                       ))**2)
               endif
            endif

            if(present(cpt_build_matrix))then
               if(.not.allocated(epsilon_cpt)) allocate(epsilon_cpt(2*(Ntot + &
                    Nc), 2*(Ntot + Nc)))
               if(.not.allocated(T_cpt)) allocate( T_cpt(2*(Ntot + Nc), &
                    2*(Ntot + Nc)))
               epsilon_cpt = 0.d0

               do k = 0, ncpt_chain_coup-1
                  do i = 1, Nb
                     j = Nb + k + (i-1)*ncpt_approx + 1
                     epsilon_cpt(i, j) = cptup(i + k*Nb)
                     epsilon_cpt(j, i) = conjg(epsilon_cpt(i, j))
                  enddo
                  do i = 1, Nb
                     j = Nb + k + (i-1)*ncpt_approx + 1
                     epsilon_cpt(Ntot + i, Ntot + j) = -cptdn(i + k*Nb)
                     epsilon_cpt(Ntot + j, Ntot + i) = conjg(epsilon_cpt(Ntot &
                          + i, Ntot + j))
                  enddo
               enddo
               if(supersc_state .and. .not.force_nupdn_basis .and. &
                    .not.force_no_pairing)then
                  do k = 0, ncpt_chain_coup-1
                     do i = 1, Nb
                        j = Nb + k + (i-1)*ncpt_approx + 1
                        epsilon_cpt(i, Ntot + j) = cptbcs(i + k*Nb)
                        epsilon_cpt(j, Ntot + i) = cptbcs(i + k*Nb)
                        epsilon_cpt(Ntot + j, i) = cptbcs(i + k*Nb)
                        epsilon_cpt(Ntot + i, j) = cptbcs(i + k*Nb)
                     enddo
                  enddo
               endif
               T_cpt = epsilon_cpt
               epsilon_cpt = 0.d0
               epsilon_cpt(2*Nc + 1:2*Nc + 2*Ntot, 2*Nc + 1:2*Nc + 2*Ntot) = &
                    T_cpt(1:2*Ntot, 1:2*Ntot)
               T_cpt = 0.d0
            endif
            do k = 0, ncpt_chain_coup-1
               do i = 1, Nb
                  j = Nb + k + (i-1)*ncpt_approx + 1
                  Eb(i, j) = cptup(i + k*Nb)
                  Eb(j, i) = conjg(Eb(i, j))
               enddo
               do i = 1, Nb
                  j = Nb + k + (i-1)*ncpt_approx + 1
                  Eb(Ntot + i, Ntot + j) = -cptdn(i + k*Nb)
                  Eb(Ntot + j, Ntot + i) = conjg(Eb(Ntot + i, Ntot + j))
               enddo
            enddo
            if(supersc_state .and. .not.force_nupdn_basis .and. &
                 .not.force_no_pairing)then
               do k = 0, ncpt_chain_coup-1
                  do i = 1, Nb
                     j = Nb + k + (i-1)*ncpt_approx + 1
                     Eb(i, Ntot + j) = cptbcs(i + k*Nb)
                     Eb(j, Ntot + i) = cptbcs(i + k*Nb)
                     Eb(Ntot + j, i) = cptbcs(i + k*Nb)
                     Eb(Ntot + i, j) = cptbcs(i + k*Nb)
                  enddo
               enddo
            endif
            do i = 1, Nb
               bl = ncpt_approx*(ncpt_approx + 1)/2
               j = ncpt_chain_coup*Nb + (i-1)*ncpt_approx
               Eb(j + 1:j + ncpt_approx, j + 1:j + ncpt_approx) = &
                    fill_mat_from_vector(ncpt_approx, cptup(ncpt_chain_coup*Nb &
                    + (i-1)*bl + 1:ncpt_chain_coup*Nb + i*bl))
            enddo
            do i = 1, Nb
               bl = ncpt_approx*(ncpt_approx + 1)/2
               j = ncpt_chain_coup*Nb + (i-1)*ncpt_approx
               Eb(Ntot + j + 1:Ntot + j + ncpt_approx, Ntot + j + 1:Ntot + j + &
                    ncpt_approx) = &
                    -transpose(fill_mat_from_vector(ncpt_approx, &
                    cptdn(ncpt_chain_coup*Nb + (i-1)*bl + 1:ncpt_chain_coup*Nb &
                    + i*bl)))
            enddo
            if(supersc_state .and. .not.force_nupdn_basis .and. &
                 .not.force_no_pairing)then
               do i = 1, Nb
                  bl = ncpt_approx*(ncpt_approx + 1)/2
                  j = ncpt_chain_coup*Nb + (i-1)*ncpt_approx
                  Eb(j + 1:j + ncpt_approx, Ntot + j + 1:Ntot + j + &
                       ncpt_approx) = fill_mat_from_vector(ncpt_approx, &
                       cptbcs(ncpt_chain_coup*Nb + (i-1)*bl + &
                       1:ncpt_chain_coup*Nb + i*bl))
                  Eb(Ntot + j + 1:Ntot + j + ncpt_approx, j + 1:j + &
                       ncpt_approx) = Eb(j + 1:j + ncpt_approx, Ntot + j + &
                       1:Ntot + j + ncpt_approx)
               enddo
            endif
            if(present(cpt_build_matrix))then
               T_cpt                                             = 0.d0
               T_cpt ( 1: Nc, 2*Nc + 1:2*Nc + Nb ) = Vcb ( 1 : Nc, 1: Nb )
               T_cpt ( Nc + 1:2*Nc, 2*Nc + Ntot + 1:2*Nc + Ntot + Nb ) = Vcb ( &
                    Nc + 1 : 2*Nc, Ntot + 1:Ntot + Nb )
               T_cpt = T_cpt + transpose(conjg(T_cpt))
               T_cpt (    1:2*Nc,           1:2*Nc             ) = Eccc%rc%mat
               T_cpt (2*Nc + 1:2*Nc + 2*Ntot, 2*Nc + 1:2*Nc + 2*Ntot)     = Eb
               if(maxval(abs( T_cpt - transpose(conjg( T_cpt )) )) > 1.d-4)then
                  write(*, *) 'error T_cpt non hermitian'
                  stop
               endif
            endif

         endif

         Eb = -Eb

      endif

      !-------------------------------------!
      SELECT CASE (FREQTYPE)
      CASE(FERMIONIC)
         hybrid => BATH%hybrid%fctn
         freq   => BATH%hybrid%freq%vec
      CASE(RETARDED)
         hybrid => BATH%hybridret%fctn
         freq   => BATH%hybridret%freq%vec
         CASE DEFAULT
         write(*, *) "ERROR IN bath2hybrid: ALLOWED FREQUENCY TYPES ARE " // &
              FERMIONIC // " AND " // RETARDED
         stop
      END SELECT
      !-------------------------------------!

      Nw = SIZE(freq)
      if(present(short))then
         if(fit_nw > 0)Nw = min(Nw, fit_nw)
      endif

      !-----------------------------------------------!
      ! Delta(iw) = Vcb * ( iw 1 - Ebath )^(-1) * Vbc !
      !-----------------------------------------------!

      hybrid = 0.0_DBL
      block = .not.BATH%SUPER .or. (maxval(abs(BATH%Pb%rc%mat)) < cutoff_rvb &
           .and. maxval(abs(BATH%PVbc(1)%rc%mat)) < cutoff_rvb .and. &
           maxval(abs(BATH%PVbc(2)%rc%mat)) < cutoff_rvb)
      diagmat =  diag_bath .and. .not.bath_nearest_hop

      !---------------------------------------------------------------------------!
      !---------------------------------------------------------------------------!
      if(freeze_poles_delta)then
         if(ncpt_tot > 0)then
            write(*, *) 'ERROR freezing pole with CPT approximation, not &
                 &implemented...'
            stop
         endif
         if(.not.allocated(FROZEN_POLES)) allocate(FROZEN_POLES(2*BATH%Nb))
         Ebm1 = -Eb
         CALL eigenvector_matrix(lsize = size(Ebm1, 1), mat = Ebm1, vaps = &
              eigenvalues, eigenvec = Ebm1)

         Vcb = MATMUL(Vcb, Ebm1)
         Vbc = MATMUL(conjg(transpose(Ebm1)), Vbc)
         ratio = dble(freeze_poles_delta_iter - 1) / &
              dble(freeze_poles_delta_niter) * freeze_pole_lambda
         if(freeze_poles_delta_iter == 1)then
            FROZEN_POLES = eigenvalues
         else
            eigenvalues = FROZEN_POLES
            do i = 1, BATH%Nb-1
               eigenvalues(i) = FROZEN_POLES(i) + ratio*(FROZEN_POLES(i + &
                    1)-FROZEN_POLES(i))
            enddo
            do i = BATH%Nb + 1, 2*BATH%Nb-1
               eigenvalues(i) = FROZEN_POLES(i) + ratio*(FROZEN_POLES(i + &
                    1)-FROZEN_POLES(i))
            enddo
         endif
!#ifndef OPENMP_MPI_SAFE
!$OMP PARALLEL PRIVATE(iw)
!$OMP DO
!#endif
         do iw = 1, Nw
            hybrid(:, :, iw) = MATMUL(Vcb, MATMUL_x(aa = bande_mat(2*BATH%Nb, &
                 1.d0/(freq(iw)-eigenvalues)), bb = Vbc, IdL = .true.))
         enddo
!#ifndef OPENMP_MPI_SAFE
!$OMP END DO
!$OMP END PARALLEL
!#endif

      else
         if(.not.fast_fit .or. (fit_green_func_and_not_hybrid .and. &
              present(short)))then

            if(fit_green_func_and_not_hybrid .and. present(short))then
               !#ifndef OPENMP_MPI_SAFE
               !$OMP PARALLEL PRIVATE(iw, Ebm1, i) SHARED(Eb, freq)
               !$OMP DO
               !#endif
               do iw = 1, Nw
                  Ebm1 = Eb
                  do i = 1, size(Eb, 1)
                     Ebm1(i, i) = Ebm1(i, i) + freq(iw)
                  enddo
                  CALL invmat(n = size(Ebm1, 1), mat = Ebm1, block_matrix = &
                       block, diagmat = diagmat)
                  if(.not.diagmat)then
                     hybrid(:, :, iw) = MATMUL(Vcb, MATMUL_x(aa = Ebm1, bb = &
                          Vbc, a_by_block = block))
                  else
                     hybrid(:, :, iw) = MATMUL(Vcb, MATMUL_x(aa = Ebm1, bb = &
                          Vbc, IdL = .true.))
                  endif
                  hybrid(:, :, iw) = -hybrid(:, :, iw)
                  do i = 1, size(hybrid, 1)
                     hybrid(i, i, iw) = hybrid(i, i, iw) + freq(iw)
                  enddo
                  hybrid(:, :, iw) = hybrid(:, :, iw) - Eccc%rc%mat
                  CALL invmat(n = size(hybrid, 1), mat = hybrid(:, :, iw), &
                       block_matrix = block, diagmat = diagmat)
               enddo
               !#ifndef OPENMP_MPI_SAFE
               !$OMP END DO
               !$OMP END PARALLEL
               !#endif
            else
               mmstep = size(Eb, 1)/size(hybrid, 1)
               !#ifndef OPENMP_MPI_SAFE
               !$OMP PARALLEL PRIVATE(iw, Ebm1, i, mm, ll) SHARED(Eb, freq)
               !$OMP DO
               !#endif
               do iw = 1, Nw
                  if(.not.(diag_V .and. diagmat))then
                     Ebm1 = Eb
                     do i = 1, size(Eb, 1)
                        Ebm1(i, i) = Ebm1(i, i) + freq(iw)
                     enddo
                     CALL invmat(n = size(Ebm1, 1), mat = Ebm1, block_matrix = &
                          block, diagmat = diagmat)
                     if(.not.diagmat)then
                        hybrid(:, :, iw) = MATMUL(Vcb, MATMUL_x(aa = Ebm1, bb &
                             = Vbc, a_by_block = block))
                     else
                        hybrid(:, :, iw) = MATMUL(Vcb, MATMUL_x(aa = Ebm1, bb &
                             = Vbc, IdL = .true.))
                     endif
                  else
                     hybrid(:, :, iw) = 0.d0
                     do mm = 1, size(hybrid, 1)
                        hybrid(mm, mm, iw) = 0.d0
                        if(fmos)then
                           do ll = (mm-1)*mmstep + 1, mm*mmstep
                              hybrid(mm, mm, iw) = hybrid(mm, mm, iw) + &
                                   conjg(Vcb(mm, ll))*Vcb(mm, ll)/(freq(iw) + &
                                   Eb(ll, ll))
                           enddo
                        else
                           do ll = 1, size(Eb, 1)
                              hybrid(mm, mm, iw) = hybrid(mm, mm, iw) + &
                                   conjg(Vcb(mm, ll))*Vcb(mm, ll)/Ebm1(ll, ll)
                           enddo
                        endif
                     enddo
                  endif
               enddo
               !#ifndef OPENMP_MPI_SAFE
               !$OMP END DO
               !$OMP END PARALLEL
               !#endif
            endif
         else
            Ebm1 = -Eb
            CALL eigenvector_matrix(lsize = size(Ebm1, 1), mat = Ebm1, vaps = &
                 eigenvalues, eigenvec = Ebm1)

            Vcb = MATMUL(Vcb, Ebm1)
            Vbc = MATMUL(conjg(transpose(Ebm1)), Vbc)
!#ifndef OPENMP_MPI_SAFE
!$OMP PARALLEL PRIVATE(iw)
!$OMP DO
!#endif
            do iw = 1, Nw
               hybrid(:, :, iw) = MATMUL(Vcb, MATMUL_x(aa = &
                    bande_mat(size(eigenvalues), 1.d0/(freq(iw)-eigenvalues)), &
                    bb = Vbc, IdL = .true.))
            enddo
!#ifndef OPENMP_MPI_SAFE
!$OMP END DO
!$OMP END PARALLEL
!#endif
         endif
      endif

      IF(write_ .AND. iproc == 1) CALL glimpse_correl(BATH%hybrid)

      CALL delete_masked_matrix(EbNambu)
      CALL delete_masked_matrix(VbcNambu)

      return
   end subroutine

   subroutine hybrid2bath(bath)

      use bath_class_vec,        only: bath2vec, vec2bath
      use common_def,            only: dump_message, reset_timer, timer_fortran
      use correl_class,          only: average_correlations, copy_correl, &
           new_correl
      use genvar,                only: dbl, fermionic, log_unit, no_mpi, &
           rank, retarded, size2
      use bath_class,            only: bath_type, copy_bath, new_bath
      use globalvar_ed_solver,   only: average_g, dist_max, &
           fit_green_func_and_not_hybrid, fit_meth, iterdmft, mask_average, &
           ncpt_flag_two_step_fit, ncpt_para, ncpt_tot, niter_search_max, niter_search_max_0, &
           param_input, param_output, search_step, skip_fit_, tolerance, &
           use_specific_set_parameters, verbose_graph
      use impurity_class,        only: Eccc
      use mpirout,               only: mpibarrier, mpibcast
      use matrix,                only: invmat, write_array
      use minimization_wrapping, only: minimize_func_wrapper
      use random,                only: dran_tab

      implicit none

      TYPE(bath_type), INTENT(INOUT) :: bath
      INTEGER                :: spin, iparam, ii, j, Nw, i, iw
      INTEGER                :: start_hybrid2bath, istep
      REAL(DBL)              :: dist_min, dist_test
      REAL(DBL), ALLOCATABLE :: test(:)

      IF(.NOT.ASSOCIATED(bath%Eb)) STOP "ERROR IN hybrid2bath : INPUT ISNT &
           &ASSOCIATED!"

      CALL reset_timer(start_hybrid2bath)

      ! CREATE THE REFERENCE HYBRIDIZATION WE WANT TO FIT
      IF(.NOT.ASSOCIATED(hybrid2fit%fctn)) CALL new_correl(hybrid2fit, &
           bath%hybrid)
      CALL copy_correl(hybrid2fit, bath%hybrid)
      hybrid2fit%title = 'hybrid2fit'

      if(abs(average_G) >= 4)then
         call average_correlations(bath%hybrid, hybrid2fit%fctn, average_G >= &
              0, MASK_AVERAGE)
         bath%hybrid%fctn = hybrid2fit%fctn
         call bath2vec(bath)
      endif

      if(fit_green_func_and_not_hybrid) then
         write(*, *) ' ... FITTING GREEN FUNCTION ... '
         call write_array( Eccc%rc%mat(:, :), 'Ec for Fit', unit = 6, short = &
              .true.)
         Nw = size(hybrid2fit%fctn, 3)
         write(*, *) 'Nw = ', Nw
         do iw = 1, Nw
            hybrid2fit%fctn(:, :, iw) = - hybrid2fit%fctn(:, :, iw)
            do i = 1, size(hybrid2fit%fctn, 1)
               hybrid2fit%fctn(i, i, iw) = hybrid2fit%fctn(i, i, iw) + &
                    hybrid2fit%freq%vec(iw)
            enddo
            hybrid2fit%fctn(:, :, iw) = hybrid2fit%fctn(:, :, iw) - &
                 Eccc%rc%mat(:, :)
            CALL invmat(n = size(hybrid2fit%fctn, 1), mat = hybrid2fit%fctn(:, &
                 :, iw))
         enddo
         write(*, *) 'Weiss Field Function to fit defined'
         call write_array( Eccc%rc%mat(:, :), 'Ec for Fit', unit = 6, short = &
              .true.)
         call write_array( real(hybrid2fit%fctn(:, :, 1)), 'Re Weiss to fit &
              &(iw = 1)', unit = 6, short = .true.)
         call write_array( aimag(hybrid2fit%fctn(:, :, 1)), 'Re Weiss to fit &
              &(iw = 1)', unit = 6, short = .true.)
      endif

      ! CREATE THE RUNNING BATH

      IF(.NOT.ASSOCIATED(batht%Eb)) CALL new_bath(batht, bath)
      CALL copy_bath(batht, bath)
      CALL bath2vec(batht)
      CALL bath2vec(bath)

      CALL dump_message(TEXT = "############################")
      CALL dump_message(TEXT = "### BEGIN HYBRID => BATH ###")
      CALL dump_message(TEXT = "############################")

      if(size2 > 1 .and. .not.no_mpi)then
         write(*, *) 'BEGIN HYBRIDIZATION FIT, RANK = ', rank
         call mpibarrier
      endif

      write(log_unit, *) "parameters for minimization = "
      write(log_unit, *) 'nparam, search_step, dist_max, Niter_search_max : '
      write(log_unit, *) bath%nparam, bath%search_step, bath%dist_max, &
           bath%Niter_search_max

      if(allocated(test)) deallocate(test)
      allocate(test(bath%nparam + ncpt_para*ncpt_tot))
      test                           =   0.d0
      test(1:bath%nparam)            =   bath%vec(1:bath%nparam)
      if(ncpt_tot > 0)then
         do i = 1, ncpt_tot
            test(bath%nparam + i) = (-1.d0 + 2.d0*dran_tab(i))/2.d0/1000.d0
            if(ncpt_para == 2)then
               test(bath%nparam + i + ncpt_tot)   = test(bath%nparam + i)
            endif
         enddo
      endif

      if(use_specific_set_parameters .and. (ncpt_tot == 0 .or. &
           .not.ncpt_flag_two_step_fit)) test(1:bath%nparam + &
           ncpt_para*ncpt_tot) = param_input(1:bath%nparam + &
           ncpt_para*ncpt_tot)

      if(.not.skip_fit_)then
         if(iterdmft == 1)then
            istep = Niter_search_max_0
         else
            istep = bath%Niter_search_max
         endif

         write(*, *) 'minimizing with [x] total parameters - rank : ', &
              bath%nparam + ncpt_para*ncpt_tot, rank
         write(*, *) 'size of array test (parameters) - rank : ', size(test), &
              rank

         if(size(test) == 0)then
            write(*, *) 'ERROR - array test has zero dimension in bath_hybrid'
            write(*, *) '        rank = ', rank
            stop
         endif

         icount_ = 0

         if(size2 > 1 .and. .not.no_mpi)then
            write(*, *) 'CALLING WRAPPER, RANK = ', rank
            call mpibarrier
         endif

#ifdef OPENMP_MPI_SAFE___
         if(rank == 0 .or. no_mpi)then
            write(*, *) 'COMPUTING DISTANCE_FUNC'
            call distance_func(dist_test, size(test), test)
            write(*, *) 'INITIAL DISTANCE : ', dist_test
            write(*, *) 'CALLING minimize_func .....'
            call minimize_func_wrapper(distance_func, test, bath%nparam + &
                 ncpt_para*ncpt_tot, FIT_METH, istep, dist_min, bath%dist_max, &
                 bath%search_step, use_mpi = .false.)
         endif

         write(*, *) 'RANK WAITING FOR MASTER : ', rank
         write(*, *) 'no_mpi flag             : ', no_mpi

         if(.not.no_mpi)then
            call mpibarrier
            call mpibcast(test)
         endif
#else
         write(*, *) 'COMPUTING DISTANCE_FUNC'
         call distance_func(dist_test, size(test), test)
         write(*, *) 'INITIAL DISTANCE : ', dist_test
         write(*, *) 'CALLING minimize_func .....'
         call minimize_func_wrapper(distance_func, test, bath%nparam + &
              ncpt_para*ncpt_tot, FIT_METH, istep, dist_min, bath%dist_max, &
              bath%search_step, use_mpi = .true.)
#endif
      else
         write(*, *) 'SKIPPING FIT'
         dist_min = 0.d0
      endif

      if(size2 > 1 .and. .not.no_mpi)then
         write(*, *) 'COLLECTING OPTIMAL PARAMETERS, RANK = ', rank
         call mpibarrier
      endif

      bath%vec(1:bath%nparam)                        = test(1:bath%nparam)
      param_output = 0.d0
      param_output(1:bath%nparam + ncpt_para*ncpt_tot) = test(1:bath%nparam + &
           ncpt_para*ncpt_tot)

      CALL vec2bath(bath)

      if(ncpt_tot == 0)then
         CALL bath2hybrid(bath, FERMIONIC)
         CALL plot_some_results
         CALL bath2hybrid(bath, RETARDED)
      else
         !to store epsilon_cpt and T_cpt
         CALL bath2hybrid(bath, FERMIONIC, cptvec = test(bath%nparam + &
              1:size(test)), cpt_build_matrix = .true.)
         CALL plot_some_results
         !we solve and obtain the GF of the disconnected problem
         CALL bath2hybrid(bath, FERMIONIC)
         CALL bath2hybrid(bath, RETARDED)
      endif

      deallocate(test)

      if(verbose_graph)then
         ! call plotarray(real(bath%hybridret%freq%vec), &
         !      real(bath%hybridret%fctn(1,1,:)), 'ED hybrid Bath re bath')
         ! call plotarray(real(bath%hybridret%freq%vec), &
         ! real(bath%hybridret%fctn(2, 2, :)), 'ED hybrid Bath re bath2')
         ! call plotarray(real(bath%hybridret%freq%vec), &
         ! aimag(bath%hybridret%fctn(1, 1, :)), 'ED hybrid Bath im bath')
         ! call plotarray(real(bath%hybridret%freq%vec), &
         ! aimag(bath%hybridret%fctn(2, 2, :)), 'ED hybrid Bath im bath2')
      endif

      write(*, *) 'min bath parameters are : ', param_output
      WRITE(log_unit, '(2(a, f0.12), a)') "# END of conjugate gradient: &
           &dist_min = ", dist_min, " (tolerance = ", bath%dist_max, ")"
      CALL timer_fortran(start_hybrid2bath, "### HYBRID => BATH TOOK ")

   contains

      subroutine plot_some_results()

         use stringmanip,         only: tostring
         use globalvar_ed_solver, only: fit_all_elements_show_graphs, iterdmft

         implicit none

         character(12) :: lab
         integer       :: kkk

         if(rank /= 0) return
         !call PGSUBP(4, 4)
         do j = 1, size(bath%hybrid%fctn(:, 1, 1))
            do kkk = 1, size(bath%hybrid%fctn(1, :, 1))

               if(j == kkk .or. fit_all_elements_show_graphs)then
                  lab = TRIM(ADJUSTL(toString(j))) // "_" // &
                       TRIM(ADJUSTL(toString(kkk))) // "_iter" // &
                       TRIM(ADJUSTL(toString(iterdmft))) // "_"
                  if(.not.fit_green_func_and_not_hybrid)then
                     !!call plotarray( aimag(hybrid2fit%freq%vec),
                     ! real(bath%hybrid%fctn(j, kkk, :)),
                     ! aimag(hybrid2fit%freq%vec), real(hybrid2fit%fctn(j, kkk,
                     ! :)), 'FIT_ED_re_bath' // lab, inset = .true., nn = 100)
                     !!call plotarray( aimag(hybrid2fit%freq%vec),
                     ! aimag(bath%hybrid%fctn(j, kkk, :)),
                     ! aimag(hybrid2fit%freq%vec), aimag(hybrid2fit%fctn(j, kkk,
                     ! :)), 'FIT_ED_im_bath' // lab, inset = .true., nn = 100)
                  else
                     !!call plotarray( aimag(hybrid2fit%freq%vec),
                     ! real(batht%hybrid%fctn(j, kkk, :)),
                     ! aimag(hybrid2fit%freq%vec), real(hybrid2fit%fctn(j, kkk,
                     ! :)), 'FIT_ED_re_bath_t' // lab, inset = .true., nn = 100)
                     !!call plotarray( aimag(hybrid2fit%freq%vec),
                     ! aimag(batht%hybrid%fctn(j, kkk, :)),
                     ! aimag(hybrid2fit%freq%vec), aimag(hybrid2fit%fctn(j, kkk,
                     ! :)), 'FIT_ED_im_bath_t' // lab, inset = .true., nn = 100)
                  endif
               endif
            enddo
         enddo
         !call PGSUBP(1, 1)

      end subroutine

   end subroutine

   subroutine distance_func(dist, n, vec)

      ! EXTRACT THE NEW BATH PARAMETERS OUT OF VECTOR 'vec'

      use genvar,              only: dbl, fermionic, log_unit, rank, size2
      use globalvar_ed_solver, only: cpt_lagrange, ncpt_para, ncpt_tot, &
           verbose_graph
      use bath_class_vec,      only: vec2bath
      use correl_class,        only: diff_correl

      implicit none

      REAL(DBL), INTENT(OUT) :: dist
      INTEGER, INTENT(IN)    :: n
      REAL(DBL), INTENT(IN)  :: vec(n)
      REAL(DBL) :: Vweight

      if(n < batht%nparam)then
         write(*, *) 'error in routine sent to conjg gradient'
         write(*, *) 'dimension given as input : ', n
         write(*, *) 'batht%nparam             : ', batht%nparam
         stop 'termine'
      endif

      batht%vec(1:batht%nparam) = vec(1:batht%nparam)
      CALL vec2bath(batht)

      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      !$$ COMPUTE THE NEW HYBRIZATION FUNCTION $$
      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      if(ncpt_tot == 0)then
         CALL bath2hybrid(batht, FERMIONIC, WRITE_HYBRID = .false., short = .true.)
      else
         if(n == batht%nparam)then
            write(*, *) 'ERROR using CPT approximation, total parameters &
                 &mismatch1'
            stop
         endif
         if(ncpt_para*ncpt_tot + batht%nparam /= n)then
            write(*, *) 'ERROR using CPT approximation, total parameters &
                 &mismatch2'
            stop
         endif
         CALL bath2hybrid(batht, FERMIONIC, WRITE_HYBRID = .false., short = .true., &
              cptvec = vec(batht%nparam + 1:n), Vweight = Vweight)
      endif

      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      !$$ COMPUTE THE DISTANCE TO THE DESIRED HYBRIZATION FUNCTION $$
      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      CALL diff_correl(dist, batht%hybrid, hybrid2fit)

      if(ncpt_tot /= 0)then
         dist = dist + Vweight*cpt_lagrange
      endif

      icount_ = icount_ + 1
      if(mod(icount_, 200) == 0) write(log_unit, *) ' dist = ', dist, icount_, &
           rank
      if(mod(icount_, 20) == 0) write(*, *) ' dist = ', dist, icount_, rank

      if(mod(icount_, 100) == 0 .and. verbose_graph .and. size2 == 1) call &
           plot_it_

   end subroutine

   subroutine plot_it_()

      implicit none

      ! call simplehist(dble(batht%vec), 'vec', display = 5)
      ! call plotarray( aimag(hybrid2fit%freq%vec), real(batht%hybrid%fctn(1,
      ! 1, :)), aimag(hybrid2fit%freq%vec), real(hybrid2fit%fctn(1, 1, :)),
      ! 're bath' , inset = .true., nn = 100, display = 1)
      ! call plotarray( aimag(hybrid2fit%freq%vec), real(batht%hybrid%fctn(2,
      ! 2, :)), aimag(hybrid2fit%freq%vec), real(hybrid2fit%fctn(2, 2, :)),
      ! 're bath2', inset = .true., nn = 100, display = 2)
      ! call plotarray( aimag(hybrid2fit%freq%vec), aimag(batht%hybrid%fctn(1,
      ! 1, :)), aimag(hybrid2fit%freq%vec), aimag(hybrid2fit%fctn(1, 1, :)),
      ! 'im bath', inset = .true., nn = 100, display = 3)
      ! call plotarray( aimag(hybrid2fit%freq%vec), aimag(batht%hybrid%fctn(2,
      ! 2, :)), aimag(hybrid2fit%freq%vec), aimag(hybrid2fit%fctn(2, 2, :)),
      ! 'im bath2', inset = .true., nn = 100, display = 4)

   end subroutine

end module
