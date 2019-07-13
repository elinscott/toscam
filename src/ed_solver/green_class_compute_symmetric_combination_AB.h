
   SUBROUTINE reshuffle_vecAB(greenAB, greenBA, greenAA, greenBB, &
                              reshuffle_dyn, BA)

      use common_def, only: find_rank
      use correl_class, only: correl_type
      use genvar, only: DP, log_unit
      use green_class, only: green_type
      use mask_class, only: mask_type
      use masked_matrix_class, only: masked_matrix_type

      TYPE(green_type), INTENT(INOUT)  :: greenAB
      TYPE(green_type), INTENT(INOUT)  :: greenBA
      TYPE(green_type), INTENT(IN)     :: greenAA
      TYPE(green_type), INTENT(IN)     :: greenBB
      INTEGER                          :: nnn, ijj, jii
      LOGICAL, INTENT(IN)     :: reshuffle_dyn, BA
      COMPLEX(DP)                     :: swapstat
      TYPE(masked_matrix_type)         :: statAA(2, 2), statAB(2, 2), &
                                          statBB(2, 2), statBA(2, 2)
      COMPLEX(DP), POINTER            :: vec_ii(:) => NULL(), vec_jj(:) => NULL(), vec_ij(:) => NULL(), &
                                          vec_ji(:) => NULL()
      COMPLEX(DP), POINTER            :: tvec_ji(:) => NULL(), tvec_ij(:) => NULL()
      LOGICAL                          :: already_diag, already_diagstat
#ifdef _complex
      COMPLEX(DP), POINTER            :: vecstat_ii => NULL(), vecstat_jj => NULL(), vecstat_ij => NULL()
      COMPLEX(DP), POINTER            :: vecstat_ji => NULL(), tvecstat_ji => NULL(), tvecstat_ij => NULL()
#else
      REAL(DP), POINTER            :: vecstat_ii => NULL(), vecstat_jj => NULL(), vecstat_ij => NULL()
      REAL(DP), POINTER            :: vecstat_ji => NULL(), tvecstat_ji => NULL(), tvecstat_ij => NULL()
#endif
      COMPLEX(DP)                     :: swap_stat, swapvec(greenAB%Nw), tswapvec(greenAB%Nw)
      TYPE(correl_type)                :: correlAA(2, 2), correlBB(2, 2), &
                                          correlAB(2, 2), correlBA(2, 2)
      TYPE(mask_type)                  :: MASKAB(2, 2), MASKBA(2, 2), &
                                          MASKAA(2, 2), MASKBB(2, 2)
      INTEGER                          :: iorb, jorb, ij, ii, jj, ji, tij, &
                                          tji, ipm, jpm, iw, miw

      !-------------------------------------------------------------------------!
      ! We calculated : < 0| (A^+ - B) (A - B^+ ) |0 > instead of <0|A^B|0>
      !           and   < 0| (A^+ - A) (A - A^+ ) |0 > !
      !-------------------------------------------------------------------------!

      call new_masks_()
      call fill_the_blanks_()

      !========================================================================!
      !=================================== NORMAL =============================!
      !========================================================================!

      DO ipm = 1, 2
         DO iorb = 1, greenAB%N
            DO jorb = 1, greenAB%N
               IF (MASKAB(ipm, 3 - ipm)%mat(iorb, jorb)) THEN
                  ! THIS INCLUDES THE DIAGONAL CASE !
                  ij = find_rank(MASKAB(ipm, 3 - ipm)%imat(iorb, jorb), &
                                 MASKAB(ipm, 3 - ipm)%ivec)
                  ii = find_rank(MASKAA(ipm, 3 - ipm)%imat(iorb, iorb), &
                                 MASKAA(ipm, 3 - ipm)%ivec)
                  jj = find_rank(MASKBB(ipm, 3 - ipm)%imat(jorb, jorb), &
                                 MASKBB(ipm, 3 - ipm)%ivec)
                  ! EQUAL-TIME
                  vecstat_ij => greenAB%correlstat(ipm, 3 - ipm)%rc%vec(ij)
                  vecstat_ii => greenAA%correlstat(ipm, 3 - ipm)%rc%vec(ii)
                  vecstat_jj => greenBB%correlstat(ipm, 3 - ipm)%rc%vec(jj)
                  vecstat_ij = 0.5_DP*(vecstat_ij - vecstat_ii - vecstat_jj)
                  ! DYNAMIC
                  IF (reshuffle_dyn .AND. greenAB%compute(ipm, 3 - ipm)) THEN
                     vec_ij => greenAB%correl(ipm, 3 - ipm)%vec(ij, :)
                     vec_ii => greenAA%correl(ipm, 3 - ipm)%vec(ii, :)
                     vec_jj => greenBB%correl(ipm, 3 - ipm)%vec(jj, :)
                     vec_ij = 0.5_DP*(vec_ij - vec_ii - vec_jj)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

         IF (BA) THEN ! THIS INCLUDES THE DIAGONAL CASE !
            DO iorb = 1, greenBA%N
               DO jorb = 1, greenBA%N
                  IF (MASKBA(ipm, 3 - ipm)%mat(iorb, jorb)) THEN
                     ij = find_rank(MASKBA(ipm, 3 - ipm)%imat(iorb, jorb), &
                                    MASKBA(ipm, 3 - ipm)%ivec)
                     ii = find_rank(MASKBB(ipm, 3 - ipm)%imat(iorb, iorb), &
                                    MASKBB(ipm, 3 - ipm)%ivec)
                     jj = find_rank(MASKAA(ipm, 3 - ipm)%imat(jorb, jorb), &
                                    MASKAA(ipm, 3 - ipm)%ivec)
                     ! EQUAL-TIME
                     vecstat_ij => greenBA%correlstat(ipm, 3 - ipm)%rc%vec(ij)
                     vecstat_ii => greenBB%correlstat(ipm, 3 - ipm)%rc%vec(ii)
                     vecstat_jj => greenAA%correlstat(ipm, 3 - ipm)%rc%vec(jj)
                     vecstat_ij = 0.5_DP*(vecstat_ij - vecstat_ii - &
                                           vecstat_jj)
                     ! DYNAMIC
                     IF (reshuffle_dyn .AND. greenBA%compute(ipm, 3 - ipm)) THEN
                        vec_ij => greenBA%correl(ipm, 3 - ipm)%vec(ij, :)
                        vec_ii => greenBB%correl(ipm, 3 - ipm)%vec(ii, :)
                        vec_jj => greenAA%correl(ipm, 3 - ipm)%vec(jj, :)
                        vec_ij = 0.5_DP*(vec_ij - vec_ii - vec_jj)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDDO

#ifdef _complex

      ! THEN DISENTANGLE G(i, j) AND G(j, i) IN THE COMPLEX CASE !

      IF (.NOT. BA) then
         write (*, *) "ERROR IN reshuffle GAB : < BA > AND < AB > ARE &
              &INDEPENDANT IF H IS COMPLEX"
         STOP
      ENDIF
      !----------------------------------------------!
      !     G[Ai,Bj]=   (G[AiBj]+G[BjAi])/2 if i<=j  !
      ! but G[Ai,Bj]= I*(G[AiBj]-G[BjAi])/2 if i> j  !
      !                                              !
      ! and G[Bi,Aj]=   (G[BiAj]+G[AjBi])/2 if i< j  !
      ! but G[Bi,Aj]= I*(G[BiAj]-G[AjBi])/2 if i>=j  !
      !----------------------------------------------!

      !AB-BA

      DO ipm = 1, 2

         ! THAT INCLUDES DIAGONAL CASE !

         DO iorb = 1, greenAB%N
            DO jorb = 1, greenAB%N

               IF (ff(iorb, jorb, diag=.true.) .and. MASKAB(ipm, &
                                                            3 - ipm)%mat(iorb, jorb)) THEN

                  ij = find_rank(MASKAB(ipm, 3 - ipm)%imat(iorb, jorb), &
                                 MASKAB(ipm, 3 - ipm)%ivec)
                  ji = find_rank(MASKBA(ipm, 3 - ipm)%imat(jorb, iorb), &
                                 MASKBA(ipm, 3 - ipm)%ivec)

                  vecstat_ij => greenAB%correlstat(ipm, 3 - ipm)%rc%vec(ij)
                  vecstat_ji => greenBA%correlstat(ipm, 3 - ipm)%rc%vec(ji)

                  swap_stat = vecstat_ij + imi*vecstat_ji
                  vecstat_ij = vecstat_ij - imi*vecstat_ji
                  vecstat_ji = swap_stat

                  IF (reshuffle_dyn .AND. greenAB%compute(ipm, 3 - ipm)) THEN
                     vec_ij => greenAB%correl(ipm, 3 - ipm)%vec(ij, :)
                     vec_ji => greenBA%correl(ipm, 3 - ipm)%vec(ji, :)
                     swapvec = vec_ij + imi*vec_ji
                     vec_ij = vec_ij - imi*vec_ji
                     vec_ji = swapvec
                  ENDIF

                  ij = find_rank(MASKBA(ipm, 3 - ipm)%imat(iorb, jorb), &
                                 MASKBA(ipm, 3 - ipm)%ivec)
                  ji = find_rank(MASKAB(ipm, 3 - ipm)%imat(jorb, iorb), &
                                 MASKAB(ipm, 3 - ipm)%ivec)

                  ! EQUAL TIME
                  vecstat_ij => greenBA%correlstat(ipm, 3 - ipm)%rc%vec(ij)
                  vecstat_ji => greenAB%correlstat(ipm, 3 - ipm)%rc%vec(ji)
                  swap_stat = vecstat_ij + imi*vecstat_ji
                  vecstat_ij = vecstat_ij - imi*vecstat_ji
                  vecstat_ji = swap_stat

                  IF (reshuffle_dyn .AND. greenAB%compute(ipm, 3 - ipm)) THEN
                     vec_ij => greenBA%correl(ipm, 3 - ipm)%vec(ij, :)
                     vec_ji => greenAB%correl(ipm, 3 - ipm)%vec(ji, :)
                     swapvec = vec_ij + imi*vec_ji
                     vec_ij = vec_ij - imi*vec_ji
                     vec_ji = swapvec
                  ENDIF

               ENDIF

            ENDDO
         ENDDO
      ENDDO

#endif

      !========================================================================!
      !================================= ANOMALOUS ============================!
      !========================================================================!

      DO ipm = 1, 2
         DO iorb = 1, greenAB%N
            DO jorb = 1, greenAB%N

               IF (MASKAB(ipm, ipm)%mat(iorb, jorb)) THEN

                  !----------------------------------------------------!
                  ! RESHUFFLE ALSO DIAGONAL PART                       !
                  ! G[Ai,Bj] <- ( G[Ai,Bj] - G[Ai,ai] - G[bj,Bj] ) / 2 !
                  ! G[ai,bj] <- ( G[ai,bj] - G[ai,Ai] - G[Bj,bj] ) / 2 !
                  !----------------------------------------------------!

                  ij = find_rank(MASKAB(ipm, ipm)%imat(iorb, jorb), MASKAB( &
                                 ipm, ipm)%ivec)
                  ii = find_rank(MASKAA(ipm, 3 - ipm)%imat(iorb, iorb), MASKAA( &
                                 ipm, 3 - ipm)%ivec)
                  jj = find_rank(MASKBB(3 - ipm, ipm)%imat(jorb, jorb), &
                                 MASKBB(3 - ipm, ipm)%ivec)

                  ! EQUAL TIME
                  vecstat_ij => greenAB%correlstat(ipm, ipm)%rc%vec(ij)
                  vecstat_ii => greenAA%correlstat(ipm, 3 - ipm)%rc%vec(ii)
                  vecstat_jj => greenBB%correlstat(3 - ipm, ipm)%rc%vec(jj)
                  vecstat_ij = 0.5_DP*vecstat_ij - 0.5_DP*coef1*vecstat_ii - &
                               0.5_DP*coef2*vecstat_jj

                  ! DYNAMIC
                  IF (reshuffle_dyn .AND. greenAB%compute(ipm, ipm)) THEN
                     vec_ij => greenAB%correl(ipm, ipm)%vec(ij, :)
                     vec_ii => greenAA%correl(ipm, 3 - ipm)%vec(ii, :)
                     vec_jj => greenBB%correl(3 - ipm, ipm)%vec(jj, :)
                     vec_ij = 0.5_DP*vec_ij - 0.5_DP*coef1*vec_ii - &
                              0.5_DP*coef2*vec_jj
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

#ifdef _complex

      DO ipm = 1, 2
         IF (BA) THEN

            DO iorb = 1, greenBA%N
               DO jorb = 1, greenBA%N
                  IF (MASKBA(ipm, ipm)%mat(iorb, jorb)) THEN

                     !-----------------------------------------------------!
                     ! G[Bi,Aj] <- ( G[Bi,Aj] - G[Bi,bi] - G[aj,Aj] ) / 2  !
                     ! G[bi,aj] <- ( G[bi,aj] - G[bi,Bi] - G[Aj,aj] ) / 2  !
                     !-----------------------------------------------------!

                     ij = find_rank(MASKBA(ipm, ipm)%imat(iorb, jorb), &
                                    MASKBA(ipm, ipm)%ivec)
                     ii = find_rank(MASKBB(ipm, 3 - ipm)%imat(iorb, iorb), &
                                    MASKBB(ipm, 3 - ipm)%ivec)
                     jj = find_rank(MASKAA(3 - ipm, ipm)%imat(jorb, jorb), &
                                    MASKAA(3 - ipm, ipm)%ivec)

                     ! EQUAL TIME
                     vecstat_ij => greenBA%correlstat(ipm, ipm)%rc%vec(ij)
                     vecstat_ii => greenBB%correlstat(ipm, 3 - ipm)%rc%vec(ii)
                     vecstat_jj => greenAA%correlstat(3 - ipm, ipm)%rc%vec(jj)
                     vecstat_ij = 0.5_DP*vecstat_ij - 0.5_DP*coef1*vecstat_ii - &
                                  0.5_DP*coef2*vecstat_jj

                     ! DYNAMIC
                     IF (reshuffle_dyn .AND. greenBA%compute(ipm, ipm)) THEN
                        vec_ij => greenBA%correl(ipm, ipm)%vec(ij, :)
                        vec_ii => greenBB%correl(ipm, 3 - ipm)%vec(ii, :)
                        vec_jj => greenAA%correl(3 - ipm, ipm)%vec(jj, :)
                        vec_ij = 0.5_DP*vec_ij - 0.5_DP*coef1*vec_ii - &
                                 0.5_DP*coef2*vec_jj
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO

         ENDIF
      ENDDO

#endif

      CALL PREPARE_ANOMALOUS()
      !   CALL DUMP_OUTPUT_TEST(greenAB, greenBA, BA)

      !========================================================================!
      !=============================== ANOMALOUS OFFDIAG ======================!
      !========================================================================!

#ifdef _complex

      nnn = greenAB%N

      DO ipm = 1, 2
         DO iorb = 1, greenAB%N
            DO jorb = 1, greenAB%N
               IF (ff(iorb, jorb, diag=.true.) .and. MASKAB(ipm, &
                                                            ipm)%mat(iorb, jorb)) THEN

                  ! THAT INCLUDES THE DIAGONAL CASE !

                  write (log_unit, *) 'ipm iorb, jorb : ', ipm, iorb, jorb, &
                     greenAB%compute(ipm, ipm), greenBA%compute(ipm, ipm)

                  ij = find_rank(MASKAB(ipm, ipm)%imat(iorb, jorb), MASKAB( &
                                 ipm, ipm)%ivec)
                  ji = find_rank(MASKAB(ipm, ipm)%imat(jorb, iorb), MASKAB( &
                                 ipm, ipm)%ivec)
                  tij = find_rank(MASKBA(3 - ipm, 3 - ipm)%imat(iorb, jorb), &
                                  MASKBA(3 - ipm, 3 - ipm)%ivec)
                  tji = find_rank(MASKBA(3 - ipm, 3 - ipm)%imat(jorb, iorb), &
                                  MASKBA(3 - ipm, 3 - ipm)%ivec)

                  tvecstat_ji => statBA(3 - ipm, 3 - ipm)%rc%vec(tji)
                  vecstat_ij => greenAB%correlstat(ipm, ipm)%rc%vec(ij)
                  ! -> AB_11_ij
                  vecstat_ij = vecstat_ij - imi*tvecstat_ji

                  IF (reshuffle_dyn .AND. greenAB%compute(ipm, ipm)) THEN
                     !-------------------------------------------!
                     ! F[Ai, Bj] - imi * TF[Bj, Ai] = G[Ai, Bj]  !
                     !-------------------------------------------!
                     tvec_ji => correlBA(3 - ipm, 3 - ipm)%vec(tji, :)
                     vec_ij => greenAB%correl(ipm, ipm)%vec(ij, :) ! -> AB_11_ij
                     vec_ij = vec_ij - imi*tvec_ji
                  ENDIF

                  if (iorb /= jorb) then
                     tvecstat_ji => greenAB%correlstat(ipm, ipm)%rc%vec(ji)
                     ! -> AB_11_ji
                     vecstat_ij => statBA(3 - ipm, 3 - ipm)%rc%vec(tij)
                     tvecstat_ji = vecstat_ij + imi*tvecstat_ji

                     IF (reshuffle_dyn .AND. greenAB%compute(ipm, ipm)) THEN
                        !-------------------------------------------!
                        ! TF[Bi, Aj] + imi *  F[Aj, Bi] = G[Aj, Bi] !
                        !-------------------------------------------!
                        tvec_ji => greenAB%correl(ipm, ipm)%vec(ji, :)
                        ! -> AB_11_ji
                        vec_ij => correlBA(3 - ipm, 3 - ipm)%vec(tij, :)
                        tvec_ji = vec_ij + imi*tvec_ji
                     endif
                  endif

               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO ipm = 1, 2

         DO iorb = 1, greenBA%N
            DO jorb = 1, greenBA%N
               IF (ff(iorb, jorb, diag=.true.) .and. MASKBA(ipm, &
                                                            ipm)%mat(iorb, jorb)) THEN

                  ! THAT INCLUDES THE DIAGONAL CASE !

                  write (log_unit, *) 'ipm iorb, jorb : ', ipm, iorb, jorb, &
                     greenAB%compute(ipm, ipm), greenBA%compute(ipm, ipm)

                  ij = find_rank(MASKBA(ipm, ipm)%imat(iorb, jorb), MASKBA( &
                                 ipm, ipm)%ivec)
                  ji = find_rank(MASKBA(ipm, ipm)%imat(jorb, iorb), MASKBA( &
                                 ipm, ipm)%ivec)

                  tij = find_rank(MASKAB(3 - ipm, 3 - ipm)%imat(iorb, jorb), &
                                  MASKAB(3 - ipm, 3 - ipm)%ivec)
                  tji = find_rank(MASKAB(3 - ipm, 3 - ipm)%imat(jorb, iorb), &
                                  MASKAB(3 - ipm, 3 - ipm)%ivec)

                  if (iorb /= jorb) then
                     tvecstat_ji => statAB(3 - ipm, 3 - ipm)%rc%vec(tji)
                     vecstat_ij => greenBA%correlstat(ipm, ipm)%rc%vec(ij)
                     ! -> BA _11_ij
                     vecstat_ij = vecstat_ij - imi*tvecstat_ji
                     IF (reshuffle_dyn .AND. greenBA%compute(ipm, ipm)) THEN
                        !-----------------------------------------------!
                        ! F[Bi, Aj] - imi * TF[Aj, Bi] =   G[Bi, Aj]    !
                        !-----------------------------------------------!
                        tvec_ji => correlAB(3 - ipm, 3 - ipm)%vec(tji, :)
                        vec_ij => greenBA%correl(ipm, ipm)%vec(ij, :)
                        ! -> BA _11_ij
                        vec_ij = vec_ij - imi*tvec_ji
                     endif
                  endif

                  tvecstat_ji => greenBA%correlstat(ipm, ipm)%rc%vec(ji)
                  ! -> BA_11_ji
                  vecstat_ij => statAB(3 - ipm, 3 - ipm)%rc%vec(tij)
                  tvecstat_ji = vecstat_ij + imi*tvecstat_ji
                  IF (reshuffle_dyn .AND. greenBA%compute(ipm, ipm)) THEN
                     !-----------------------------------------------!
                     ! TF[Ai, Bj] + imi *  F[Bj, Ai] =   G[Bj, Ai]   !
                     !-----------------------------------------------!
                     tvec_ji => greenBA%correl(ipm, ipm)%vec(ji, :)
                     ! -> BA_11_ji
                     vec_ij => correlAB(3 - ipm, 3 - ipm)%vec(tij, :)
                     tvec_ji = vec_ij + imi*tvec_ji
                  endif
               ENDIF
            ENDDO
         ENDDO
      ENDDO

#endif

      !===========================================================================!
      !================================== THE END  ===============================!
      !===========================================================================!

#ifndef _complex
      IF (reshuffle_dyn) call fill_the_blanks()
#endif
      call cleanup()

   contains

#include "green_class_compute_symmetric_combination_AB_tools.h"

   end subroutine

