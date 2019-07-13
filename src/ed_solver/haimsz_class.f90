module haimsz_class

   use fermion_hilbert_class, only: fermion_sector_type
   use genvar, only: DP

   !-------------------------------------------------------!
   ! FULL HAMILTONIAN IN FERMION SECTOR WITH TOTAL SPIN Sz !
   !-------------------------------------------------------!

   IMPLICIT NONE

   private

   INTERFACE HAIMsz_mult
      MODULE PROCEDURE HAIMsz_multc, HAIMsz_multr
   END INTERFACE

   INTERFACE HAIMsz_mult_fly
      MODULE PROCEDURE HAIMsz_mult_fly_r, HAIMsz_mult_fly_c
   END INTERFACE

   REAL(DP), ALLOCATABLE                     :: QUART_INT_SZ(:)
   TYPE(fermion_sector_type), POINTER, PRIVATE :: sector_sz => NULL()
   INTEGER, PRIVATE :: istatemin, istatemax, dchunk, &
                       dimen
   LOGICAL, PRIVATE :: offdiag_coulomb

   public :: delete_HAIMsz
   public :: HAIMsz_mult
   public :: HAIMsz_mult_fly
   public :: new_HAIMsz

contains

   subroutine new_HAIMsz(AIM, sector_in)

      use fermion_hilbert_class, only: fermion_sector_type
      use masked_matrix_class_mod, only: masked_real_matrix_type, &
         new_masked_real_matrix
      use genvar, only: iproc
      use matrix, only: new_diag
      use aim_class, only: aim_type
      use haim2_class, only: new_haim2

      implicit none

      TYPE(AIM_type), INTENT(IN)                    :: AIM
      TYPE(fermion_sector_type), INTENT(IN), TARGET :: sector_in
      TYPE(masked_real_matrix_type) :: U
      LOGICAL                       :: is_diag(AIM%Nc, AIM%Nc)
      INTEGER                       :: spin, start_tabH

      sector_sz => sector_in
      CALL new_HAIM2(AIM, sector_sz)

      CALL new_masked_real_matrix(U, AIM%impurity%U)

      U%MASK%mat = .false.
      WHERE (U%mat /= 0.0_DP)
      U%MASK%mat = .true.
      END WHERE

      offdiag_coulomb = ANY(U%MASK%mat .AND. (.NOT. is_diag))

      IF (ALLOCATED(QUART_INT_SZ)) DEALLOCATE (QUART_INT_SZ)
      ALLOCATE (QUART_INT_SZ(sector_sz%chunk(iproc)))
      CALL tab_HAIMsz(U)

      dimen = sector_sz%dimen
      istatemin = sector_sz%istatemin(iproc)
      istatemax = sector_sz%istatemax(iproc)
      dchunk = sector_sz%chunk(iproc)

      CALL new_diag(is_diag, AIM%Nc)

   end subroutine

   subroutine delete_HAIMsz()

      use haim2_class, only: delete_haim2

      implicit none

      NULLIFY (sector_sz)
      CALL delete_HAIM2()
      IF (ALLOCATED(QUART_INT_SZ)) DEALLOCATE (QUART_INT_SZ)
   end subroutine

   subroutine tab_HAIMsz(U)

      ! TABULATE THE FULL HAMILTONIAN IN THE SECTOR OF TOTAL SPIN Sz

      use masked_matrix_class_mod, only: masked_real_matrix_type
      use fermion_ket_class, only: fermion_ket_type
      use globalvar_ed_solver, only: jhund, jhund_slater_type, jjmatrix, &
         on_fly
      use genvar, only: iproc
      use haim2_class, only: tab_haim2, AIM2sz

      implicit none

      TYPE(masked_real_matrix_type), INTENT(IN) :: U
      INTEGER, POINTER       :: IMPiorb(:) => NULL()
      INTEGER                :: site1, site2, n1, n2, Nc
      INTEGER                :: istate, istateloc, istatemin, istatemax, state
      LOGICAL                :: bn(2), bn2(2), bn1(2)
      TYPE(fermion_ket_type) :: ket_in

      istatemin = sector_sz%istatemin(iproc)
      istatemax = sector_sz%istatemax(iproc)
      IMPiorb => AIM2sz%IMPiorb
      Nc = AIM2sz%Nc

      !------------------!
      ! QUADRATIC PART   !
      !------------------!

      if (.not. ON_FLY) CALL tab_HAIM2(AIM2sz, sector_sz) ! TABULATE QUADRATIC
      ! PART

      !---------------------------------------------------------------!
      ! QUARTIC PART (IMPURITY ONLY)                                  !
      ! WARNING: QUADRATIC PART WAS WRITTEN IN NAMBU BASIS I.E.       !
      ! SPIN DOWN PART WAS PARTICLE-HOLE TRANSFORMED                  !
      ! SO WE TAKE CARE OF THAT ALSO WHEN TABULATING THE QUARTIC PART !
      !---------------------------------------------------------------!

      QUART_INT_SZ = 0.0_DP

      DO istate = istatemin, istatemax ! span the local chunk of the sector

         istateloc = istate - istatemin + 1
         state = sector_sz%state(istate)

         ! ONSITE COULOMB
         DO site1 = 1, Nc
            bn1 = (/BTEST(state, IMPiorb(site1) - 1),.NOT. BTEST(state, &
                                                                 IMPiorb(site1 + Nc) - 1)/)
            n1 = COUNT(bn1)
            IF (n1 == 2 .AND. U%MASK%mat(site1, site1)) QUART_INT_SZ(istateloc) &
               = QUART_INT_SZ(istateloc) + U%mat(site1, site1)

            ! INTER-SITE COULOMB
            IF (offdiag_coulomb .AND. n1 /= 0) THEN
               DO site2 = site1 + 1, Nc
                  IF (U%MASK%mat(site1, site2)) THEN
                     bn = (/BTEST(state, IMPiorb(site2) - 1),.NOT. BTEST(state, &
                                                                         IMPiorb(site2 + Nc) - 1)/)
                     n2 = COUNT(bn)
                     QUART_INT_SZ(istateloc) = QUART_INT_SZ(istateloc) + &
                                               U%mat(site1, site2)*n1*n2
                  ENDIF
               ENDDO
            ENDIF

            IF (abs(Jhund) > 1.d-5 .and. Jhund_Slater_type) THEN
               DO site2 = site1 + 1, Nc
                  bn2 = (/BTEST(state, IMPiorb(site2) - 1),.NOT. BTEST(state, &
                                                                       IMPiorb(site2 + Nc) - 1)/)
                  if (bn1(1) .and. bn2(1)) then
                     QUART_INT_SZ(istateloc) = QUART_INT_SZ(istateloc) - &
                                               JJmatrix(site1, site2)
                  endif
                  if (bn1(2) .and. bn2(2)) then
                     QUART_INT_SZ(istateloc) = QUART_INT_SZ(istateloc) - &
                                               JJmatrix(site1, site2)
                  endif
               ENDDO
            ENDIF

         ENDDO

      ENDDO

   end subroutine

   subroutine HAIMsz_multc(vec_out, vec_in)

      ! WE COMPUTE THE RELEVANT CHUNK OF vec_out = H * vec_in

      use fermion_Hilbert_class, only: HILBERT_SPACE_SPLITED_AMONG_NODES
      use genvar, only: dp, ierr, iproc, size2
      use openmpmod, only: omp_get_num_threads, omp_get_thread_num, &
         openmp_split_array
      use globalvar_ed_solver, only: flag_mpi_greens, jhund, JHund_Slater_type, &
         jjmatrix, open_mp, use_cc
      use haim2_class, only: AIM2sz, diagsz, noffsz, offdiagsz, rankoffsz
      use linalg, only: long_sum
      use lockmod, only: MAXT
      use mpi
      use mpi_mod, only: mpibcast, split

      implicit none

      COMPLEX(DP), INTENT(INOUT) :: vec_out(:)
      COMPLEX(DP), INTENT(IN)    :: vec_in(:)
      COMPLEX(DP), ALLOCATABLE :: vec_tot_out(:)
      COMPLEX(DP)              :: hoffdiag
      INTEGER                   :: istateloc, istate, jstate, noff, irank, &
                                   istatemin0, istatemax0, imin_(MAXT), &
                                   imax_(MAXT), TID, istatemin_, j
      INTEGER(8)                :: nhoffdiag
      INTEGER                   :: Nc
      integer                   :: site1, site2, n1, n2, i1, i2, ket_in_, &
                                   ket_out_, jj, fermion_sign_
      logical                   :: bn1(2), bn2(2), go_for_omp
      INTEGER, POINTER          :: IMPiorb(:) => NULL()
      logical                   :: emptychunks

      if (dimen == 1) then
         vec_out(1) = 1
         return
      endif

      IMPiorb => AIM2sz%IMPiorb

      Nc = AIM2sz%Nc

      emptychunks = any(sector_sz%chunk == 0)

      if (.not. USE_CC .or. emptychunks) vec_out = 0.d0

      if (sector_sz%chunk(iproc) == 0) then
         write (*, *) 'EMPTY SECTOR - SIZE VEC_OUT : ', size(vec_out)
         goto 35
      endif

#ifndef OPENMP_MPI
      if (OPEN_MP) call openmp_split_array(dimen, imin_, imax_)
#else
      if (OPEN_MP) then
         call openmp_split_array(istatemax - istatemin + 1, imin_, imax_)
         where (imin_ /= 0) imin_ = imin_ + istatemin - 1
         where (imax_ /= 0) imax_ = imax_ + istatemin - 1
      endif
#endif

      go_for_omp = .false.
      if (OPEN_MP) go_for_omp = minval(imax_) > 0 .and. minval(sector_sz%chunk) &
                                > 0 .and. minval(imin_) > 0

!$OMP       PARALLEL IF(go_for_omp) PRIVATE(site1, site2, n1, n2, i1, i2, &
!$OMP          ket_in_, ket_out_, jj, fermion_sign_, bn1, bn2, hoffdiag, TID, &
!$OMP          istatemin_, istatemin0, istatemax0, istate, istateloc, noff, irank, &
!$OMP          jstate, nhoffdiag) SHARED(vec_out, vec_in, noffsz, rankoffsz, &
!$OMP          offdiagsz)

      if (OPEN_MP .and. go_for_omp) then
         TID = OMP_GET_THREAD_NUM() + 1
         istatemin0 = imin_(TID)
         istatemax0 = imax_(TID)
         istatemin_ = 1
#ifdef OPENMP_MPI
         istatemin_ = istatemin
#endif
#ifndef OPENMP_MPI
         if (istatemin0 > 1) then
#else
            if (istatemin0 > istatemin_) then
#endif
#ifdef OPENMP_MPI
               nhoffdiag = long_sum(noffsz(1:istatemin0 - istatemin_))
#else
               nhoffdiag = long_sum(noffsz(1:istatemin0 - 1))
#endif
            else
               nhoffdiag = 0
            endif
         else
            istatemin0 = istatemin
            istatemax0 = istatemax
            istatemin_ = istatemin
            nhoffdiag = 0
         endif

         if (OPEN_MP .and. go_for_omp) then
            if (OMP_GET_NUM_THREADS() /= MAXT) Then
               write (*, *) 'ERROR obtained number of threads is different from &
                    &MAXT'
               write (*, *) 'MAXT = ', MAXT
               write (*, *) 'NUMBER OF THREADS = ', OMP_GET_NUM_THREADS()
               stop
            endif
         endif

         if (istatemax0 == 0) goto 1192

         DO istate = istatemin0, istatemax0
#ifdef OPENMP_MPI
            if (istatemax0 == 0) cycle
#endif
            istateloc = istate - istatemin_ + 1
            if (.not. USE_CC) then
               vec_out(istate) = vec_out(istate) + QUART_INT_SZ(istateloc)* &
                                 vec_in(istate)
            else
               vec_out(istate) = QUART_INT_SZ(istateloc)*vec_in(istate)
            endif
            vec_out(istate) = vec_out(istate) + diagsz(istateloc)* &
                              vec_in(istate)
            noff = noffsz(istateloc)

            DO irank = 1, noff
               jstate = rankoffsz(nhoffdiag + irank)
               hoffdiag = offdiagsz(nhoffdiag + irank)
               if (.not. USE_CC) vec_out(jstate) = vec_out(jstate) + hoffdiag* &
                                                   vec_in(istate)
               vec_out(istate) = vec_out(istate) + conjg(hoffdiag)* &
                                 vec_in(jstate)
            ENDDO

            nhoffdiag = nhoffdiag + noff
         ENDDO

         IF (abs(Jhund) > 1.d-5 .and. offdiag_coulomb) THEN

            !this term is -2 x Jhund x Si Sj (spin operators)

            DO istate = istatemin0, istatemax0

               istateloc = istate - istatemin_ + 1

               DO site1 = 1, Nc

                  bn1 = (/BTEST(sector_sz%state(istate), IMPiorb(site1) - 1), &
                          BTEST(sector_sz%state(istate), IMPiorb(site1 + Nc) - 1)/)
                  n1 = COUNT(bn1)

                  DO site2 = site1 + 1, Nc

                     bn2 = (/BTEST(sector_sz%state(istate), IMPiorb(site2) - 1), &
                             BTEST(sector_sz%state(istate), IMPiorb(site2 + &
                                                                    Nc) - 1)/)
                     n2 = COUNT(bn2)

                     IF (n1 == 2 .and. n2 == 0) THEN
                        ket_in_ = sector_sz%state(istate)
                        ket_out_ = IBCLR(ket_in_, IMPiorb(site1) - 1)
                        ket_out_ = IBCLR(ket_out_, IMPiorb(site1 + Nc) - 1)
                        ket_out_ = IBSET(ket_out_, IMPiorb(site2 + Nc) - 1)
                        ket_out_ = IBSET(ket_out_, IMPiorb(site2) - 1)
                        i1 = min(IMPiorb(site1), IMPiorb(site2))
                        i2 = max(IMPiorb(site1), IMPiorb(site2))
                        fermion_sign_ = 1
                        DO jj = i1 + 1, i2 - 1
                           IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                              -fermion_sign_
                        ENDDO
                        i1 = min(IMPiorb(site2 + Nc), IMPiorb(site1 + Nc))
                        i2 = max(IMPiorb(site2 + Nc), IMPiorb(site1 + Nc))
                        DO jj = i1 + 1, i2 - 1
                           IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                              -fermion_sign_
                        ENDDO
                        jstate = sector_sz%rank(ket_out_)
                        vec_out(istate) = vec_out(istate) - &
                                          vec_in(jstate)*JJmatrix(site1, &
                                                                  site2)*dble(fermion_sign_)
                     ENDIF

                     IF (n1 == 0 .and. n2 == 2) THEN
                        ket_in_ = sector_sz%state(istate)
                        ket_out_ = IBCLR(ket_in_, IMPiorb(site2) - 1)
                        ket_out_ = IBCLR(ket_out_, IMPiorb(site2 + Nc) - 1)
                        ket_out_ = IBSET(ket_out_, IMPiorb(site1 + Nc) - 1)
                        ket_out_ = IBSET(ket_out_, IMPiorb(site1) - 1)
                        i1 = min(IMPiorb(site1), IMPiorb(site2))
                        i2 = max(IMPiorb(site1), IMPiorb(site2))
                        fermion_sign_ = 1
                        DO jj = i1 + 1, i2 - 1
                           IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                              -fermion_sign_
                        ENDDO
                        i1 = min(IMPiorb(site2 + Nc), IMPiorb(site1 + Nc))
                        i2 = max(IMPiorb(site2 + Nc), IMPiorb(site1 + Nc))
                        DO jj = i1 + 1, i2 - 1
                           IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                              -fermion_sign_
                        ENDDO
                        jstate = sector_sz%rank(ket_out_)
                        vec_out(istate) = vec_out(istate) - &
                                          vec_in(jstate)*JJmatrix(site1, &
                                                                  site2)*dble(fermion_sign_)
                     ENDIF

                     if (.not. Jhund_Slater_type) then
                        ! -2J Si^z Sj^z : up up -- > -J/2 (ni-1) (nj-1)
                        vec_out(istate) = vec_out(istate) - 0.50d0* &
                                          JJmatrix(site1, site2)*dble(n1 - 1)*dble(n2 - 1 &
                                                                                   )*vec_in(istate)
                     endif

                  ENDDO
               ENDDO

            ENDDO
         ENDIF

         IF (abs(Jhund) > 1.d-5 .and. offdiag_coulomb) THEN

            if (Jhund < 0.) then
               write (*, *) 'Weird, hund coupling is negative'
               stop
            endif

            !this term is Jhund x double_hopping(up + dn) from site1 to site2

            DO istate = istatemin0, istatemax0
               istateloc = istate - istatemin_ + 1

               DO site1 = 1, Nc
                  bn1 = (/BTEST(sector_sz%state(istate), IMPiorb(site1) - 1), &
                          BTEST(sector_sz%state(istate), IMPiorb(site1 + Nc) - 1)/)
                  IF (bn1(1) .and. .not. bn1(2)) THEN
                     DO site2 = 1, Nc
                        if (site1 /= site2) then
                           bn2 = (/BTEST(sector_sz%state(istate), &
                                         IMPiorb(site2) - 1), &
                                   BTEST(sector_sz%state(istate), IMPiorb(site2 + &
                                                                          Nc) - 1)/)
                           IF (bn2(2) .and. .not. bn2(1)) THEN
                              ket_in_ = sector_sz%state(istate)
                              ket_out_ = IBCLR(ket_in_, IMPiorb(site1) - 1) !up
                              ket_out_ = IBCLR(ket_out_, IMPiorb(site2 + &
                                                                 Nc) - 1) !do
                              ket_out_ = IBSET(ket_out_, IMPiorb(site1 + &
                                                                 Nc) - 1) !do
                              ket_out_ = IBSET(ket_out_, IMPiorb(site2) - 1) !up

                              i1 = min(IMPiorb(site1), IMPiorb(site2))
                              i2 = max(IMPiorb(site1), IMPiorb(site2))
                              fermion_sign_ = 1
                              DO jj = i1 + 1, i2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              i1 = min(IMPiorb(site1 + Nc), IMPiorb(site2 + Nc))
                              i2 = max(IMPiorb(site1 + Nc), IMPiorb(site2 + Nc))
                              DO jj = i1 + 1, i2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              vec_out(istate) = vec_out(istate) - &
                                                vec_in(jstate)*JJmatrix(site1, &
                                                                        site2)*dble(fermion_sign_)
                           ENDIF
                        endif
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

1192     continue

!$OMP       END PARALLEL

35       continue

         if (HILBERT_SPACE_SPLITED_AMONG_NODES) then
            if (FLAG_MPI_GREENS > 0) then
               write (*, *) 'you try to split the hilbert space amongst nodes, &
                    &but also'
               write (*, *) 'at the same time you are parallelizing the &
                    &computation of green functions'
               stop 'critical'
            endif
            if (USE_CC .and. .not. emptychunks) then
               call mpibcast(vec_out, sector_sz%chunk(:), &
                             [(sum(sector_sz%chunk(1:j - 1)), j=1, size2)])
            else
               if (allocated(vec_tot_out)) deallocate (vec_tot_out)
               allocate (vec_tot_out(dimen))
               vec_tot_out = 0.
               call MPI_ALLREDUCE(vec_out, vec_tot_out, dimen, &
                                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
               vec_out = vec_tot_out
               if (allocated(vec_tot_out)) deallocate (vec_tot_out)
            endif
         endif

      end subroutine

      subroutine HAIMsz_multr(vec_out, vec_in)

         ! WE COMPUTE THE RELEVANT CHUNK OF vec_out = H * vec_in

         use fermion_Hilbert_class, only: HILBERT_SPACE_SPLITED_AMONG_NODES
         use genvar, only: dp, ierr, iproc, size2
         use globalvar_ed_solver, only: flag_mpi_greens, jhund, open_mp, use_cc
         use haim2_class, only: diagsz, noffsz, offdiagsz, rankoffsz
         use linalg, only: long_sum
         use lockmod, only: MAXT
         use mpi
         use mpi_mod, only: mpibcast, split
         use openmpmod, only: omp_get_num_threads, omp_get_thread_num, &
            openmp_split_array

         implicit none

         REAL(DP), INTENT(INOUT) :: vec_out(:)
         REAL(DP), INTENT(IN)    :: vec_in(:)
         REAL(DP), ALLOCATABLE :: vec_tot_out(:)
         REAL(DP)              :: hoffdiag
         INTEGER                :: istateloc, istate, jstate, noff, irank, &
                                   istatemin0, istatemax0, imin_(MAXT), &
                                   imax_(MAXT), TID, istatemin_, j
         INTEGER(8)             :: nhoffdiag
         LOGICAL                :: go_for_omp
         logical                :: emptychunks

         if (dimen == 1) then
            vec_out(1) = 1
            return
         endif

         IF (abs(Jhund) > 1.d-5) THEN
            write (*, *) 'Jhund and Sz basis : not yet implemented'
            stop
         ENDIF

         emptychunks = any(sector_sz%chunk == 0)

         if (.not. USE_CC .or. emptychunks) vec_out = 0.d0

         if (sector_sz%chunk(iproc) == 0) then
            write (*, *) 'EMPTY SECTOR - SIZE VEC_OUT : ', size(vec_out)
            goto 35
         endif

#ifndef OPENMP_MPI
         if (OPEN_MP) call openmp_split_array(dimen, imin_, imax_)
#else
         if (OPEN_MP) then
            call openmp_split_array(istatemax - istatemin + 1, imin_, imax_)
            where (imin_ /= 0) imin_ = imin_ + istatemin - 1
            where (imax_ /= 0) imax_ = imax_ + istatemin - 1
         endif
#endif

         go_for_omp = .false.
         if (OPEN_MP) go_for_omp = minval(imax_) > 0 .and. minval(sector_sz%chunk) &
                                   > 0 .and. minval(imin_) > 0

!$OMP       PARALLEL IF(go_for_omp) &
!$OMP       PRIVATE(hoffdiag, TID, istatemin_, istatemin0, istatemax0, istate, &
!$OMP            istateloc, noff, irank, jstate, nhoffdiag) &
!$OMP       SHARED(vec_out, vec_in, noffsz, rankoffsz, offdiagsz)

         if (OPEN_MP .and. go_for_omp) then
            TID = OMP_GET_THREAD_NUM() + 1
            istatemin0 = imin_(TID)
            istatemax0 = imax_(TID)
            istatemin_ = 1
#ifdef OPENMP_MPI
            istatemin_ = istatemin
#endif
#ifndef OPENMP_MPI
            if (istatemin0 > 1) then
#else
               if (istatemin0 > istatemin_) then
#endif
#ifdef OPENMP_MPI
                  nhoffdiag = long_sum(noffsz(1:istatemin0 - istatemin_))
#else
                  nhoffdiag = long_sum(noffsz(1:istatemin0 - 1))
#endif
               else
                  nhoffdiag = 0
               endif
            else
               istatemin0 = istatemin
               istatemax0 = istatemax
               istatemin_ = istatemin
               nhoffdiag = 0
            endif

            if (OPEN_MP .and. go_for_omp) then
               if (OMP_GET_NUM_THREADS() /= MAXT) Then
                  write (*, *) 'ERROR obtained number of threads is different from &
                       &MAXT'
                  write (*, *) 'MAXT = ', MAXT
                  write (*, *) 'NUMBER OF THREADS = ', OMP_GET_NUM_THREADS()
                  stop
               endif
            endif

            if (istatemax0 == 0) goto 1192

            DO istate = istatemin0, istatemax0
#ifdef OPENMP_MPI
               if (istatemax0 == 0) cycle
#endif
               istateloc = istate - istatemin_ + 1
               if (.not. USE_CC) then
                  vec_out(istate) = vec_out(istate) + QUART_INT_SZ(istateloc)* &
                                    vec_in(istate)
               else
                  vec_out(istate) = QUART_INT_SZ(istateloc)*vec_in(istate)
               endif
               vec_out(istate) = vec_out(istate) + diagsz(istateloc)* &
                                 vec_in(istate)
               noff = noffsz(istateloc)

               DO irank = 1, noff
                  jstate = rankoffsz(nhoffdiag + irank)
                  hoffdiag = offdiagsz(nhoffdiag + irank)
                  if (.not. USE_CC) vec_out(jstate) = vec_out(jstate) + hoffdiag* &
                                                      vec_in(istate)
                  vec_out(istate) = vec_out(istate) + hoffdiag*vec_in(jstate)
               ENDDO

               nhoffdiag = nhoffdiag + noff

            ENDDO

1192        continue

!$OMP       END PARALLEL

35          continue

            if (HILBERT_SPACE_SPLITED_AMONG_NODES) then
               if (FLAG_MPI_GREENS > 0) then
                  write (*, *) 'you try to split the hilbert space amongst nodes, &
                       &but also'
                  write (*, *) 'at the same time you are parallelizing the &
                       &computation of green functions'
                  stop 'critical'
               endif
               if (USE_CC .and. .not. emptychunks) then
                  call mpibcast(vec_out, sector_sz%chunk(:), &
                                [(sum(sector_sz%chunk(1:j - 1)), j=1, size2)])
               else
                  if (allocated(vec_tot_out)) deallocate (vec_tot_out)
                  allocate (vec_tot_out(dimen))
                  vec_tot_out = 0.
                  call MPI_ALLREDUCE(vec_out, vec_tot_out, dimen, &
                                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                  vec_out = vec_tot_out
                  if (allocated(vec_tot_out)) deallocate (vec_tot_out)
               endif
            endif

         end subroutine

         subroutine HAIMsz_mult_fly_r(vec_out, vec_in)

            ! WE COMPUTE THE RELEVANT CHUNK OF vec_out = H * vec_in

            use fermion_Hilbert_class, only: HILBERT_SPACE_SPLITED_AMONG_NODES
            use genvar, only: dp, ierr, iproc, size2
            use haim2_class, only: AIM2sz
            use lockmod, only: MAXT
            use openmpmod, only: omp_get_num_threads, omp_get_thread_num, &
               openmp_split_array
            use mpi
            use mpi_mod, only: mpibcast, split
            use globalvar_ed_solver, only: flag_mpi_greens, jhund, open_mp, use_cc

            implicit none

            REAL(DP), INTENT(INOUT) :: vec_out(:)
            REAL(DP), INTENT(IN)    :: vec_in(:)
            REAL(DP), ALLOCATABLE :: vec_tot_out(:)
            REAL(DP)              :: hoffdiag
            INTEGER                :: istateloc, istate, jstate, irank, istatemin0, &
                                      istatemax0, imin_(MAXT), imax_(MAXT), TID, &
                                      istatemin_, j, iorb, jorb
            INTEGER                :: IMPnorbs, BATHnorbs, i, ket_in_, ket_out_, &
                                      norbs_, jj, n1, n2
            INTEGER(2)             :: fermion_sign_
            INTEGER, POINTER       :: IMPiorb(:) => NULL(), BATHiorb(:) => NULL()
            LOGICAL                :: go_for_omp
            logical                :: emptychunks

            if (dimen == 1) then
               vec_out(1) = 1
               return
            endif

            if (.not. USE_CC) stop 'error on fly needs USE_CC'

            IF (abs(Jhund) > 1.d-5) THEN
               write (*, *) 'Jhund and Sz basis : not yet implemented'
               stop
            ENDIF

            IMPiorb => AIM2sz%IMPiorb
            IMPnorbs = AIM2sz%IMPnorbs
            BATHiorb => AIM2sz%BATHiorb
            BATHnorbs = AIM2sz%BATHnorbs

            emptychunks = any(sector_sz%chunk == 0)

            if (emptychunks) vec_out = 0.d0

            if (sector_sz%chunk(iproc) == 0) then
               write (*, *) 'EMPTY SECTOR - SIZE VEC_OUT : ', size(vec_out)
               goto 35
            endif

#ifndef OPENMP_MPI
            if (OPEN_MP) call openmp_split_array(dimen, imin_, imax_)
#else
            if (OPEN_MP) then
               call openmp_split_array(istatemax - istatemin + 1, imin_, imax_)
               where (imin_ /= 0) imin_ = imin_ + istatemin - 1
               where (imax_ /= 0) imax_ = imax_ + istatemin - 1
            endif
#endif

            go_for_omp = .false.
            if (OPEN_MP) go_for_omp = minval(imax_) > 0 .and. minval(sector_sz%chunk) &
                                      > 0 .and. minval(imin_) > 0

            !$OMP PARALLEL IF(go_for_omp) PRIVATE(hoffdiag, TID, istatemin_, &
            !$OMP    istatemin0, istatemax0, istate, istateloc, irank, jstate, &
            !$OMP    iorb, jorb, i, jj, j, ket_in_, ket_out_, fermion_sign_, &
            !$OMP    norbs_) SHARED(vec_out, vec_in, QUART_INT_SZ)

            if (OPEN_MP .and. go_for_omp) then
               TID = OMP_GET_THREAD_NUM() + 1
               istatemin0 = imin_(TID)
               istatemax0 = imax_(TID)
               istatemin_ = 1
#ifdef OPENMP_MPI
               istatemin_ = istatemin
#endif
            else
               istatemin0 = istatemin
               istatemax0 = istatemax
               istatemin_ = istatemin
            endif

            if (OPEN_MP .and. go_for_omp) then
               if (OMP_GET_NUM_THREADS() /= MAXT) Then
                  write (*, *) 'ERROR obtained number of threads is different from &
                       &MAXT'
                  write (*, *) 'MAXT = ', MAXT
                  write (*, *) 'NUMBER OF THREADS = ', OMP_GET_NUM_THREADS()
                  stop
               endif
            endif

            if (istatemax0 == 0) goto 1192

            DO istate = istatemin0, istatemax0
#ifdef OPENMP_MPI
               if (istatemax0 == 0) cycle
#endif
               istateloc = istate - istatemin_ + 1

               vec_out(istate) = QUART_INT_SZ(istateloc)*vec_in(istate)

               ket_in_ = sector_sz%state(istate)
               norbs_ = sector_sz%norbs

               DO iorb = 1, IMPnorbs
                  IF (BTEST(ket_in_, IMPiorb(iorb) - 1)) vec_out(istate) = &
                     vec_out(istate) + AIM2sz%Ec%rc%mat(iorb, iorb)*vec_in(istate)
               ENDDO
               DO iorb = 1, BATHnorbs
                  IF (BTEST(ket_in_, BATHiorb(iorb) - 1)) vec_out(istate) = &
                     vec_out(istate) + AIM2sz%Eb%rc%mat(iorb, iorb)*vec_in(istate)
               ENDDO

               DO jorb = 1, IMPnorbs
                  IF (BTEST(ket_in_, IMPiorb(jorb) - 1)) then
                     DO iorb = 1, IMPnorbs
                        IF (iorb /= jorb .AND. AIM2sz%Ec%rc%MASK%mat(iorb, jorb)) THEN
                           ket_out_ = IBCLR(ket_in_, IMPiorb(jorb) - 1)
                           IF (.NOT. BTEST(ket_out_, IMPiorb(iorb) - 1)) THEN
                              ket_out_ = IBSET(ket_out_, IMPiorb(iorb) - 1)
                              n1 = min(IMPiorb(jorb), IMPiorb(iorb))
                              n2 = max(IMPiorb(jorb), IMPiorb(iorb))
                              fermion_sign_ = 1
                              DO jj = n1 + 1, n2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              hoffdiag = AIM2sz%Ec%rc%mat(iorb, jorb)* &
                                         fermion_sign_
                              vec_out(istate) = vec_out(istate) + hoffdiag* &
                                                vec_in(jstate)
                           ENDIF
                        ENDIF
                     ENDDO

                  ENDIF
               ENDDO

               DO jorb = 1, BATHnorbs

                  if (BTEST(ket_in_, BATHiorb(jorb) - 1)) then
                     DO iorb = 1, BATHnorbs
                        IF (iorb /= jorb .AND. AIM2sz%Eb%rc%MASK%mat(iorb, jorb)) THEN
                           ket_out_ = IBCLR(ket_in_, BATHiorb(jorb) - 1)
                           IF (.NOT. BTEST(ket_out_, BATHiorb(iorb) - 1)) THEN
                              ket_out_ = IBSET(ket_out_, BATHiorb(iorb) - 1)
                              n1 = min(BATHiorb(jorb), BATHiorb(iorb))
                              n2 = max(BATHiorb(jorb), BATHiorb(iorb))
                              fermion_sign_ = 1
                              DO jj = n1 + 1, n2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              hoffdiag = AIM2sz%Eb%rc%mat(iorb, jorb)*fermion_sign_
                              vec_out(istate) = vec_out(istate) + hoffdiag* &
                                                vec_in(jstate)
                           ENDIF
                        ENDIF
                     ENDDO

                  ENDIF
               ENDDO

               DO jorb = 1, IMPnorbs

                  if (BTEST(ket_in_, IMPiorb(jorb) - 1)) then
                     DO iorb = 1, BATHnorbs
                        IF (AIM2sz%Vbc%rc%MASK%mat(iorb, jorb)) THEN
                           ket_out_ = IBCLR(ket_in_, IMPiorb(jorb) - 1)
                           IF (.NOT. BTEST(ket_out_, BATHiorb(iorb) - 1)) THEN
                              ket_out_ = IBSET(ket_out_, BATHiorb(iorb) - 1)
                              n1 = min(IMPiorb(jorb), BATHiorb(iorb))
                              n2 = max(IMPiorb(jorb), BATHiorb(iorb))
                              fermion_sign_ = 1
                              DO jj = n1 + 1, n2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              hoffdiag = AIM2sz%Vbc%rc%mat(iorb, jorb)*fermion_sign_
                              vec_out(istate) = vec_out(istate) + hoffdiag* &
                                                vec_in(jstate)
                           ENDIF
                        ENDIF
                     ENDDO

                  ENDIF
               ENDDO

               DO iorb = 1, BATHnorbs

                  if (BTEST(ket_in_, BATHiorb(iorb) - 1)) then
                     DO jorb = 1, IMPnorbs
                        IF (AIM2sz%Vbc%rc%MASK%mat(iorb, jorb)) THEN
                           ket_out_ = IBCLR(ket_in_, BATHiorb(iorb) - 1)
                           IF (.NOT. BTEST(ket_out_, IMPiorb(jorb) - 1)) THEN
                              ket_out_ = IBSET(ket_out_, IMPiorb(jorb) - 1)
                              n1 = min(BATHiorb(iorb), IMPiorb(jorb))
                              n2 = max(BATHiorb(iorb), IMPiorb(jorb))
                              fermion_sign_ = 1
                              DO jj = n1 + 1, n2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              hoffdiag = AIM2sz%Vbc%rc%mat(iorb, jorb)*fermion_sign_
                              vec_out(istate) = vec_out(istate) + hoffdiag* &
                                                vec_in(jstate)
                           ENDIF
                        ENDIF
                     ENDDO

                  ENDIF
               ENDDO

            ENDDO

1192        continue

!$OMP       END PARALLEL

35          continue

            if (HILBERT_SPACE_SPLITED_AMONG_NODES) then
               if (FLAG_MPI_GREENS > 0) then
                  write (*, *) 'you try to split the hilbert space amongst nodes, but &
                       &also'
                  write (*, *) 'at the same time you are parallelizing the &
                       &computation of green functions'
                  stop 'critical'
               endif

               if (USE_CC .and. .not. emptychunks) then
                  call mpibcast(vec_out, sector_sz%chunk(:), &
                                [(sum(sector_sz%chunk(1:j - 1)), j=1, size2)])
               else
                  if (allocated(vec_tot_out)) deallocate (vec_tot_out)
                  allocate (vec_tot_out(dimen))
                  vec_tot_out = 0.
                  call MPI_ALLREDUCE(vec_out, vec_tot_out, dimen, &
                                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                  vec_out = vec_tot_out
                  if (allocated(vec_tot_out)) deallocate (vec_tot_out)
               endif

            endif

         end subroutine

         subroutine HAIMsz_mult_fly_c(vec_out, vec_in)

            ! WE COMPUTE THE RELEVANT CHUNK OF vec_out = H * vec_in

            use fermion_Hilbert_class, only: HILBERT_SPACE_SPLITED_AMONG_NODES
            use genvar, only: dp, ierr, iproc, size2
            use haim2_class, only: AIM2sz
            use lockmod, only: MAXT
            use openmpmod, only: omp_get_num_threads, omp_get_thread_num, &
               openmp_split_array
            use mpi
            use mpi_mod, only: mpibcast, split
            use globalvar_ed_solver, only: flag_mpi_greens, jhund, open_mp, use_cc

            implicit none

            COMPLEX(DP), INTENT(INOUT) :: vec_out(:)
            COMPLEX(DP), INTENT(IN)    :: vec_in(:)
            COMPLEX(DP), ALLOCATABLE :: vec_tot_out(:)
            COMPLEX(DP)              :: hoffdiag
            INTEGER                   :: istateloc, istate, jstate, irank, &
                                         istatemin0, istatemax0, imin_(MAXT), &
                                         imax_(MAXT), TID, istatemin_, j, iorb, jorb
            INTEGER                   :: IMPnorbs, BATHnorbs, i, ket_in_, ket_out_, &
                                         norbs_, jj, n1, n2
            INTEGER(2)                :: fermion_sign_
            INTEGER, POINTER          :: IMPiorb(:) => NULL(), BATHiorb(:) => NULL()
            LOGICAL                   :: go_for_omp
            logical                   :: emptychunks

            if (dimen == 1) then
               vec_out(1) = 1
               return
            endif

            IF (abs(Jhund) > 1.d-5) THEN
               write (*, *) 'Jhund and Sz basis : not yet implemented'
               stop
            ENDIF

            if (.not. USE_CC) stop 'error on fly needs USE_CC'

            IMPiorb => AIM2sz%IMPiorb
            IMPnorbs = AIM2sz%IMPnorbs
            BATHiorb => AIM2sz%BATHiorb
            BATHnorbs = AIM2sz%BATHnorbs

            emptychunks = any(sector_sz%chunk == 0)

            if (emptychunks) vec_out = 0.d0

            if (sector_sz%chunk(iproc) == 0) then
               write (*, *) 'EMPTY SECTOR - SIZE VEC_OUT : ', size(vec_out)
               goto 35
            endif

#ifndef OPENMP_MPI
            if (OPEN_MP) call openmp_split_array(dimen, imin_, imax_)
#else
            if (OPEN_MP) then
               call openmp_split_array(istatemax - istatemin + 1, imin_, imax_)
               where (imin_ /= 0) imin_ = imin_ + istatemin - 1
               where (imax_ /= 0) imax_ = imax_ + istatemin - 1
            endif
#endif

            go_for_omp = .false.
            if (OPEN_MP) go_for_omp = minval(imax_) > 0 .and. minval(sector_sz%chunk) &
                                      > 0 .and. minval(imin_) > 0

!$OMP       PARALLEL IF(go_for_omp) PRIVATE(hoffdiag, TID, istatemin_, &
!$OMP            istatemin0, istatemax0, istate, istateloc, irank, jstate, &
!$OMP            iorb, jorb, i, jj, j, ket_in_, ket_out_, fermion_sign_, &
!$OMP            norbs_) SHARED(vec_out, vec_in, QUART_INT_SZ)

            if (OPEN_MP .and. go_for_omp) then
               TID = OMP_GET_THREAD_NUM() + 1
               istatemin0 = imin_(TID)
               istatemax0 = imax_(TID)
               istatemin_ = 1
#ifdef OPENMP_MPI
               istatemin_ = istatemin
#endif
            else
               istatemin0 = istatemin
               istatemax0 = istatemax
               istatemin_ = istatemin
            endif

            if (OPEN_MP .and. go_for_omp) then
               if (OMP_GET_NUM_THREADS() /= MAXT) Then
                  write (*, *) 'ERROR obtained number of threads is different from &
                       &MAXT'
                  write (*, *) 'MAXT = ', MAXT
                  write (*, *) 'NUMBER OF THREADS = ', OMP_GET_NUM_THREADS()
                  stop
               endif
            endif

            if (istatemax0 == 0) goto 1192

            DO istate = istatemin0, istatemax0
#ifdef OPENMP_MPI
               if (istatemax0 == 0) cycle
#endif
               istateloc = istate - istatemin_ + 1

               vec_out(istate) = QUART_INT_SZ(istateloc)*vec_in(istate)

               ket_in_ = sector_sz%state(istate)
               norbs_ = sector_sz%norbs

               DO iorb = 1, IMPnorbs
                  IF (BTEST(ket_in_, IMPiorb(iorb) - 1)) vec_out(istate) = &
                     vec_out(istate) + AIM2sz%Ec%rc%mat(iorb, iorb)*vec_in(istate)
               ENDDO
               DO iorb = 1, BATHnorbs
                  IF (BTEST(ket_in_, BATHiorb(iorb) - 1)) vec_out(istate) = &
                     vec_out(istate) + AIM2sz%Eb%rc%mat(iorb, iorb)*vec_in(istate)
               ENDDO

               DO jorb = 1, IMPnorbs
                  IF (BTEST(ket_in_, IMPiorb(jorb) - 1)) then
                     DO iorb = 1, IMPnorbs
                        IF (iorb /= jorb .AND. AIM2sz%Ec%rc%MASK%mat(iorb, jorb)) THEN
                           ket_out_ = IBCLR(ket_in_, IMPiorb(jorb) - 1)
                           IF (.NOT. BTEST(ket_out_, IMPiorb(iorb) - 1)) THEN
                              ket_out_ = IBSET(ket_out_, IMPiorb(iorb) - 1)
                              n1 = min(IMPiorb(jorb), IMPiorb(iorb))
                              n2 = max(IMPiorb(jorb), IMPiorb(iorb))
                              fermion_sign_ = 1
                              DO jj = n1 + 1, n2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              hoffdiag = AIM2sz%Ec%rc%mat(iorb, jorb)* &
                                         fermion_sign_
                              vec_out(istate) = vec_out(istate) + conjg(hoffdiag)* &
                                                vec_in(jstate)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO

               DO jorb = 1, BATHnorbs

                  if (BTEST(ket_in_, BATHiorb(jorb) - 1)) then
                     DO iorb = 1, BATHnorbs
                        IF (iorb /= jorb .AND. AIM2sz%Eb%rc%MASK%mat(iorb, jorb)) THEN
                           ket_out_ = IBCLR(ket_in_, BATHiorb(jorb) - 1)
                           IF (.NOT. BTEST(ket_out_, BATHiorb(iorb) - 1)) THEN
                              ket_out_ = IBSET(ket_out_, BATHiorb(iorb) - 1)
                              n1 = min(BATHiorb(jorb), BATHiorb(iorb))
                              n2 = max(BATHiorb(jorb), BATHiorb(iorb))
                              fermion_sign_ = 1
                              DO jj = n1 + 1, n2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              hoffdiag = AIM2sz%Eb%rc%mat(iorb, jorb)* &
                                         fermion_sign_
                              vec_out(istate) = vec_out(istate) + conjg(hoffdiag)* &
                                                vec_in(jstate)
                           ENDIF
                        ENDIF
                     ENDDO

                  ENDIF
               ENDDO

               DO jorb = 1, IMPnorbs

                  if (BTEST(ket_in_, IMPiorb(jorb) - 1)) then
                     DO iorb = 1, BATHnorbs
                        IF (AIM2sz%Vbc%rc%MASK%mat(iorb, jorb)) THEN
                           ket_out_ = IBCLR(ket_in_, IMPiorb(jorb) - 1)
                           IF (.NOT. BTEST(ket_out_, BATHiorb(iorb) - 1)) THEN
                              ket_out_ = IBSET(ket_out_, BATHiorb(iorb) - 1)
                              n1 = min(IMPiorb(jorb), BATHiorb(iorb))
                              n2 = max(IMPiorb(jorb), BATHiorb(iorb))
                              fermion_sign_ = 1
                              DO jj = n1 + 1, n2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              hoffdiag = AIM2sz%Vbc%rc%mat(iorb, jorb)* &
                                         fermion_sign_
                              vec_out(istate) = vec_out(istate) + conjg(hoffdiag)* &
                                                vec_in(jstate)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO

               DO iorb = 1, BATHnorbs

                  if (BTEST(ket_in_, BATHiorb(iorb) - 1)) then
                     DO jorb = 1, IMPnorbs
                        IF (AIM2sz%Vbc%rc%MASK%mat(iorb, jorb)) THEN
                           ket_out_ = IBCLR(ket_in_, BATHiorb(iorb) - 1)
                           IF (.NOT. BTEST(ket_out_, IMPiorb(jorb) - 1)) THEN
                              ket_out_ = IBSET(ket_out_, IMPiorb(jorb) - 1)
                              n1 = min(BATHiorb(iorb), IMPiorb(jorb))
                              n2 = max(BATHiorb(iorb), IMPiorb(jorb))
                              fermion_sign_ = 1
                              DO jj = n1 + 1, n2 - 1
                                 IF (BTEST(ket_out_, jj - 1)) fermion_sign_ = &
                                    -fermion_sign_
                              ENDDO
                              jstate = sector_sz%rank(ket_out_)
                              hoffdiag = AIM2sz%Vbc%rc%mat(iorb, jorb)* &
                                         fermion_sign_
                              vec_out(istate) = vec_out(istate) + hoffdiag* &
                                                vec_in(jstate)
                           ENDIF
                        ENDIF
                     ENDDO

                  ENDIF
               ENDDO

            ENDDO

1192        continue

!$OMP       END PARALLEL

35          continue

            if (HILBERT_SPACE_SPLITED_AMONG_NODES) then
               if (FLAG_MPI_GREENS > 0) then
                  write (*, *) 'you try to split the hilbert space amongst nodes, but &
                       &also'
                  write (*, *) 'at the same time you are parallelizing the &
                       &computation of green functions'
                  stop 'critical'
               endif

               if (USE_CC .and. .not. emptychunks) then
                  call mpibcast(vec_out, sector_sz%chunk(:), &
                                [(sum(sector_sz%chunk(1:j - 1)), j=1, size2)])
               else
                  if (allocated(vec_tot_out)) deallocate (vec_tot_out)
                  allocate (vec_tot_out(dimen))
                  vec_tot_out = 0.
                  call MPI_ALLREDUCE(vec_out, vec_tot_out, dimen, &
                                     MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
                  vec_out = vec_tot_out
                  if (allocated(vec_tot_out)) deallocate (vec_tot_out)
               endif
            endif

         end subroutine

      end module
