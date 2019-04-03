   SUBROUTINE reshuffle_vecAA(green, reshuffle_dyn)

      use green_class,  only: green_type
      use common_def,   only: find_rank
      use correl_class, only: correl_type, delete_correl, transform_correl
      use mask_class,   only: delete_mask, mask_type, new_mask

      TYPE(green_type)      :: green
      LOGICAL               :: reshuffle_dyn
      COMPLEX(DBL), POINTER :: vec_ii(:) => NULL(), vec_jj(:) => NULL(), &
                               vec_ij(:) => NULL()
      COMPLEX(DBL), POINTER :: vec_ji(:) => NULL(), tvec_ji(:) => NULL(), &
                               tvec_ij(:) => NULL()
#ifdef _complex
      COMPLEX(DBL), POINTER :: vecstat_ii => NULL(), vecstat_jj => NULL(), &
                               vecstat_ij => NULL()
      COMPLEX(DBL), POINTER :: vecstat_ji => NULL(), tvecstat_ji => NULL(), &
                               tvecstat_ij => NULL()
#else
      REAL(DBL),    POINTER :: vecstat_ii => NULL(), vecstat_jj => NULL(), &
                               vecstat_ij => NULL()
      REAL(DBL),    POINTER :: vecstat_ji => NULL(), tvecstat_ji => NULL(), &
                               tvecstat_ij => NULL()
#endif
      LOGICAL,      POINTER :: MASKpp(:, :) => NULL(), MASKmm(:, :) => NULL(), &
                               MASKpm(:, :) => NULL(), MASKmp(:, :) => NULL()
      INTEGER,      POINTER :: IMASKpp(:, :) => NULL(), IMASKmm(:, :) => &
                               NULL(), IMASKpm(:, :) => NULL(), IMASKmp(:, :) &
                               => NULL()
      INTEGER,      POINTER :: IVECpp(:) => NULL(), IVECmm(:) => NULL(), &
                               IVECpm(:) => NULL(), IVECmp(:) => NULL()
      COMPLEX(DBL)          :: swap_stat, swapvec(green%Nw)
      TYPE(correl_type)     :: correl(2, 2)
      TYPE(mask_type)       :: MASK(2, 2)
      INTEGER               :: iorb, jorb, ij, ii, jj, ji, tij, tji, ipm, jpm, &
                               iw, miw

      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL new_mask(MASK(ipm, jpm), green%correlstat(ipm, jpm)%rc%MASK)
         ENDDO
      ENDDO

      !========================================================================!
      !========================    NORMAL    ==================================!
      !========================================================================!

      DO ipm = 1, 2
         DO iorb = 1, green%N
            DO jorb = 1, green%N
               IF(iorb /= jorb .AND. MASK(ipm, 3-ipm)%mat(iorb, jorb))THEN

                  ii = find_rank(MASK(ipm, 3-ipm)%imat(iorb, iorb), MASK(ipm, &
                       3-ipm)%ivec)
                  jj = find_rank(MASK(ipm, 3-ipm)%imat(jorb, jorb), MASK(ipm, &
                       3-ipm)%ivec)
                  ij = find_rank(MASK(ipm, 3-ipm)%imat(iorb, jorb), MASK(ipm, &
                       3-ipm)%ivec)

                  vecstat_ii => green%correlstat(ipm, 3-ipm)%rc%vec(ii)
                  vecstat_jj => green%correlstat(ipm, 3-ipm)%rc%vec(jj)
                  vecstat_ij => green%correlstat(ipm, 3-ipm)%rc%vec(ij)
                  vecstat_ij =  0.5_DBL * ( vecstat_ij - vecstat_ii - vecstat_jj )

                  IF(reshuffle_dyn .AND. green%compute(ipm, 3-ipm))THEN
                     vec_ii   => green%correl(ipm, 3-ipm)%vec(ii, :)
                     vec_jj   => green%correl(ipm, 3-ipm)%vec(jj, :)
                     vec_ij   => green%correl(ipm, 3-ipm)%vec(ij, :)
                     vec_ij   =  0.5_DBL * ( vec_ij - vec_ii - vec_jj )
                  ENDIF

               ENDIF
            ENDDO
         ENDDO
      ENDDO

#ifdef _complex

      ! THEN DISENTANGLE G(i, j) AND G(j, i) IN THE COMPLEX CASE
      !     G(i, j) =   (G(i, j) + G(j, i))/2 if i < j
      ! but G(i, j) = I*(G(i, j) - G(j, i))/2 if i > j

      DO ipm = 1, 2
         MASKpm  => green%correlstat(ipm, 3-ipm)%rc%MASK%mat
         IMASKpm => green%correlstat(ipm, 3-ipm)%rc%MASK%imat
         IVECpm  => green%correlstat(ipm, 3-ipm)%rc%MASK%ivec
         DO iorb = 1, green%N
            DO jorb = iorb + 1, green%N
               IF(MASK(ipm, 3-ipm)%mat(iorb, jorb))THEN
                  ij = find_rank(MASK(ipm, 3-ipm)%imat(iorb, jorb), MASK(ipm, &
                       3-ipm)%ivec)
                  ji = find_rank(MASK(ipm, 3-ipm)%imat(jorb, iorb), MASK(ipm, &
                       3-ipm)%ivec)
                  ! EQUAL TIME
                  vecstat_ij => green%correlstat(ipm, 3-ipm)%rc%vec(ij)
                  vecstat_ji => green%correlstat(ipm, 3-ipm)%rc%vec(ji)
                  swap_stat       = vecstat_ij + imi * vecstat_ji
                  vecstat_ij = vecstat_ij - imi * vecstat_ji
                  vecstat_ji = swap_stat
                  ! DYNAMIC
                  IF(reshuffle_dyn .AND. green%compute(ipm, 3-ipm))THEN
                     vec_ij   => green%correl(ipm, 3-ipm)%vec(ij, :)
                     vec_ji   => green%correl(ipm, 3-ipm)%vec(ji, :)
                     swapvec  =  vec_ij + imi * vec_ji
                     vec_ij   =  vec_ij - imi * vec_ji
                     vec_ji   = swapvec
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
#endif

      DO ipm = 1, 2
         IF(reshuffle_dyn .AND. green%compute(ipm, 3-ipm))THEN
            if(.not.green%compute(3-ipm, ipm)) CALL &
                 transform_correl(green%correl(3-ipm, ipm), green%correl(ipm, &
                 3-ipm))
         ENDIF
      ENDDO

      !========================================================================!
      !========================== ANOMALOUS ===================================!
      !========================================================================!

      DO ipm = 1, 2
         DO iorb = 1, green%N
            DO jorb = 1, green%N
               IF(MASK(ipm, ipm)%mat(iorb, jorb))THEN

                  ! RESHUFFLE ALSO DIAGONAL PART
                  ! G[Ai, Aj] <- ( G[Ai, Aj] - G[Ai, ai] - G[aj, Aj] ) / 2
                  ! G[ai, aj] <- ( G[ai, aj] - G[ai, Ai] - G[Aj, aj] ) / 2

                  ij = find_rank(MASK( ipm, ipm)%imat(iorb, jorb), MASK( ipm, &
                       ipm)%ivec)
                  ii = find_rank(MASK( ipm, 3-ipm)%imat(iorb, iorb), MASK( &
                       ipm, 3-ipm)%ivec)
                  jj = find_rank(MASK(3-ipm, ipm)%imat(jorb, jorb), &
                       MASK(3-ipm, ipm)%ivec)

                  ! EQUAL TIME
                  vecstat_ij => green%correlstat(  ipm,  ipm)%rc%vec(ij)
                  vecstat_ii => green%correlstat(  ipm, 3-ipm)%rc%vec(ii)
                  vecstat_jj => green%correlstat(3-ipm,  ipm)%rc%vec(jj)
                  vecstat_ij =  0.5_DBL * ( vecstat_ij - vecstat_ii - vecstat_jj )

                  ! DYNAMIC
                  IF(reshuffle_dyn .AND. green%compute(ipm, ipm))THEN
                     vec_ij => green%correl(  ipm,  ipm)%vec(ij, :)
                     vec_ii => green%correl(  ipm, 3-ipm)%vec(ii, :)
                     vec_jj => green%correl(3-ipm,  ipm)%vec(jj, :)
                     vec_ij =  0.5_DBL * ( vec_ij - vec_ii - vec_jj )
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

#ifdef _complex

      ! THEN DISENTANGLE G[Ai, Aj] AND G[Aj, Ai] IN THE COMPLEX CASE FIRST WE
      ! GENERATE THE MISSING CORRELATIONS

      DO ipm = 1, 2
         IF(reshuffle_dyn .AND. green%compute(ipm, ipm))THEN
            if(.not.green%compute(3-ipm, 3-ipm))then
               CALL new_correl(correl(3-ipm, 3-ipm), green%correl(ipm, ipm))
               CALL transform_correl(correl(3-ipm, 3-ipm), green%correl(ipm, &
                    ipm))
            else
               CALL new_correl(correl(3-ipm, 3-ipm), green%correl(3-ipm, &
                    3-ipm))
            endif
         ENDIF
      ENDDO

      write(log_unit, *) ' ... reshuffle vecAA ... '

      DO ipm = 1, 2
         DO iorb = 1, green%N
            DO jorb = iorb, green%N
               IF(MASK(ipm, ipm)%mat(iorb, jorb))THEN

                  ij = find_rank(MASK( ipm, ipm)%imat(iorb, jorb), MASK( ipm, &
                       ipm)%ivec)
                  tij = find_rank(MASK(3-ipm, 3-ipm)%imat(iorb, jorb), &
                       MASK(3-ipm, 3-ipm)%ivec)
                  ji = find_rank(MASK( ipm, ipm)%imat(jorb, iorb), MASK( ipm, &
                       ipm)%ivec)
                  tji = find_rank(MASK(3-ipm, 3-ipm)%imat(jorb, iorb), &
                       MASK(3-ipm, 3-ipm)%ivec)

                  vecstat_ij  => green%correlstat(  ipm,  ipm)%rc%vec( ij)
                  tvecstat_ji => green%correlstat(3-ipm, 3-ipm)%rc%vec(tji)
                  vecstat_ji  => green%correlstat(  ipm,  ipm)%rc%vec( ji)
                  tvecstat_ij => green%correlstat(3-ipm, 3-ipm)%rc%vec(tij)

                  swap_stat   = vecstat_ij + imi * tvecstat_ji
                  vecstat_ij  = vecstat_ij - imi * tvecstat_ji
                  tvecstat_ji = swap_stat

                  IF(reshuffle_dyn .AND. green%compute(ipm, ipm))THEN
                     vec_ij  => green%correl(  ipm,  ipm)%vec( ij, :)
                     tvec_ji =>       correl(3-ipm, 3-ipm)%vec(tji, :)
                     swapvec =  vec_ij + imi * tvec_ji
                     vec_ij  =  vec_ij - imi * tvec_ji
                     tvec_ji =  swapvec
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

#endif

      !============================================================================!
      !============================ THE END =======================================!
      !============================================================================!

      DO ipm = 1, 2
         IF(reshuffle_dyn .AND. green%compute(ipm, ipm))THEN
            if(.not.green%compute(3-ipm, 3-ipm)) CALL &
                 transform_correl(green%correl(3-ipm, 3-ipm), &
                 green%correl(ipm, ipm))
         ENDIF
      ENDDO

      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL delete_correl(correl(ipm, jpm))
            CALL delete_mask(MASK(ipm, jpm))
         ENDDO
      ENDDO

   END SUBROUTINE

