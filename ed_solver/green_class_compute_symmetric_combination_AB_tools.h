
   subroutine printintermid(unit)

      implicit none

      integer :: unit

      write(unit,*) '==============================='
      write(unit, *) 'ipm jpm : ', ipm, ipm
      write(unit, *) ' iorb, jorb : ', iorb, jorb
      write(unit, *) 'AA, 1 2 ii : ', greenAA%correl(  ipm, 3-ipm)%vec(ii, 4)
      write(unit, *) 'AA, 2 1 ii : ', greenAA%correl(  3-ipm, ipm)%vec(ii, 4)
      write(unit, *) 'BB, 2 1 jj : ', greenBB%correl(3-ipm,  ipm)%vec(jj, 4)
      write(unit, *) 'BB  1 2 jj : ', greenBB%correl(ipm, 3-  ipm)%vec(jj, 4)
      write(unit, *) 'AA, 1 2 jj : ', greenAA%correl(  ipm, 3-ipm)%vec(jj, 4)
      write(unit, *) 'AA, 2 1 jj : ', greenAA%correl(  3-ipm, ipm)%vec(jj, 4)
      write(unit, *) 'BB, 2 1 ii : ', greenBB%correl(3-ipm,  ipm)%vec(ii, 4)
      write(unit, *) 'BB  1 2 ii : ', greenBB%correl(ipm, 3-  ipm)%vec(ii, 4)
      write(unit, *) 'AB  1 1 ij : ', greenAB%correl(ipm, ipm)%vec(ij, 4)
      write(unit, *) 'BA  1 1 ij : ', greenBA%correl(ipm, ipm)%vec(ij, 4)
      write(unit, *) 'AB  2 2 ij : ', greenAB%correl(3-ipm, 3-ipm)%vec(ij, 4)
      write(unit, *) 'BA  2 2 ij : ', greenBA%correl(3-ipm, 3-ipm)%vec(ij, 4)
      write(unit,*) '==============================='
   end subroutine

   subroutine prepare_anomalous
   
      implicit none

      integer :: ipm

      IF(reshuffle_dyn)THEN
         DO ipm = 1, 2
            if(greenAB%compute(ipm, ipm))then
               CALL new_masked_matrix(statAB(ipm, ipm), &
                    greenAB%correlstat(ipm, ipm))
            endif
            if(BA)then
               if(greenBA%compute(ipm, ipm))then
                  CALL new_masked_matrix(statBA(ipm, ipm), &
                       greenBA%correlstat(ipm, ipm))
               endif
            endif
            IF(greenAB%compute(ipm, ipm))THEN
               CALL new_correl(correlAB(ipm, ipm), greenAB%correl(ipm, ipm))
            ENDIF
            if(BA)then
               IF(greenBA%compute(ipm, ipm))THEN
                  CALL new_correl(correlBA(ipm, ipm), greenBA%correl(ipm, ipm))
               ENDIF
            endif
         ENDDO
      ENDIF

   end subroutine

   subroutine dump_Ai_Aj
      DO ipm = 1, 2
         DO iorb = 1, greenAB%N
            DO jorb = 1, greenAB%N
               IF(MASKAB(ipm, ipm)%mat(iorb, jorb))THEN
                  ij = find_rank(MASKAB( ipm, ipm)%imat(iorb, jorb), MASKAB( &
                       ipm, ipm)%ivec)
                  ji = find_rank(MASKAB( ipm, ipm)%imat(jorb, iorb), MASKAB( &
                       ipm, ipm)%ivec)
                  write(log_unit, *) '----------iorb, jorb :---------------- &
                       &', iorb, jorb
                  write(log_unit, *) 'i- > j 1 1 AB ', greenAB%correlstat(1, &
                       1)%rc%vec(ij)
                  write(log_unit, *) 'j- > i 2 2 AB ', greenAB%correlstat(2, &
                       2)%rc%vec(ji)
                  if(BA)then
                     write(log_unit, *) 'j- > i 1 1 BA ', &
                          greenBA%correlstat(1, 1)%rc%vec(ji)
                     write(log_unit, *) 'j- > i 2 2 BA ', &
                          greenBA%correlstat(2, 2)%rc%vec(ji)
                  endif
               endif
            enddo
         enddo
      enddo

      write(log_unit,*) '======================='
      DO ipm = 1, 2
         DO iorb = 1, greenAB%N
            DO jorb = 1, greenAB%N
               IF(MASKAB(ipm, ipm)%mat(iorb, jorb))THEN
                  ij = find_rank(MASKAB( ipm, ipm)%imat(iorb, jorb), MASKAB( &
                       ipm, ipm)%ivec)
                  ji = find_rank(MASKAB( ipm, ipm)%imat(jorb, iorb), MASKAB( &
                       ipm, ipm)%ivec)
                  write(log_unit, *) '----------iorb, jorb :---------------- &
                       &', iorb, jorb
                  write(log_unit, *) 'i- > j 1 1 AB ', greenAB%correl(1, &
                       1)%vec(ij, 4)
                  write(log_unit, *) 'j- > i 2 2 AB ', greenAB%correl(2, &
                       2)%vec(ji, 4)
                  if(BA)then
                     write(log_unit, *) 'j- > i 1 1 BA ', greenBA%correl(1, &
                          1)%vec(ji, 4)
                     write(log_unit, *) 'j- > i 2 2 BA ', greenBA%correl(2, &
                          2)%vec(ji, 4)
                  endif
               endif
            enddo
         enddo
      enddo

   end subroutine

   logical function ff(iorb, jorb, diag)

      implicit none

      integer :: iorb, jorb
      logical, optional :: diag

      ff = (jorb > iorb)
      if(present(diag)) ff =  jorb >= iorb

   end function

   subroutine new_masks_()

      implicit none

      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL new_mask(MASKAA(ipm, jpm), greenAA%correlstat(ipm, &
                 jpm)%rc%MASK)
            CALL new_mask(MASKBB(ipm, jpm), greenBB%correlstat(ipm, &
                 jpm)%rc%MASK)
            CALL new_mask(MASKAB(ipm, jpm), greenAB%correlstat(ipm, &
                 jpm)%rc%MASK)
            IF(BA) CALL new_mask(MASKBA(ipm, jpm), greenBA%correlstat(ipm, &
                 jpm)%rc%MASK)
         ENDDO
      ENDDO
   end subroutine

   subroutine cleanup

      implicit none

      integer :: ipm, jmp
      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL delete_mask(MASKAA(ipm, jpm))
            CALL delete_mask(MASKBB(ipm, jpm))
            CALL delete_mask(MASKAB(ipm, jpm))
            CALL delete_mask(MASKBA(ipm, jpm))
            CALL delete_correl(correlAA(ipm, jpm))
            CALL delete_correl(correlBB(ipm, jpm))
            CALL delete_correl(correlAB(ipm, jpm))
            CALL delete_correl(correlBA(ipm, jpm))
         ENDDO
      ENDDO

   end subroutine

   subroutine clean_it()

      implicit none

      integer :: ipm, jpm
      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL delete_correl(correlAA(ipm, jpm))
            CALL delete_masked_matrix(statAA(ipm, jpm))
            CALL delete_correl(correlBB(ipm, jpm))
            CALL delete_masked_matrix(statBB(ipm, jpm))
            CALL delete_correl(correlAB(ipm, jpm))
            CALL delete_masked_matrix(statAB(ipm, jpm))
            CALL delete_correl(correlBA(ipm, jpm))
            CALL delete_masked_matrix(statBA(ipm, jpm))
         ENDDO
      ENDDO
   end subroutine

   subroutine fill_the_blanks()

      implicit none

      integer :: ipm, jmp

      DO ipm = 1, 2
         IF(.not.greenAA%compute(3-ipm, 3-ipm) .and. greenAA%compute(ipm, &
              ipm)) CALL transform_correl(greenAA%correl(3-ipm, 3-ipm), &
              greenAA%correl(ipm, ipm))
         IF(.not.greenBB%compute(3-ipm, 3-ipm) .and. greenBB%compute(ipm, &
              ipm)) CALL transform_correl(greenBB%correl(3-ipm, 3-ipm), &
              greenBB%correl(ipm, ipm))
         IF(.not.greenAB%compute(3-ipm, 3-ipm) .and. greenAB%compute(ipm, &
              ipm)) CALL transform_correl(greenAB%correl(3-ipm, 3-ipm), &
              greenAB%correl(ipm, ipm))
         if(BA)then
            IF(.not.greenAA%compute(3-ipm, 3-ipm) .and. greenBA%compute(ipm, &
                 ipm)) CALL transform_correl(greenBA%correl(3-ipm, 3-ipm), &
                 greenBA%correl(ipm, ipm))
         endif
      ENDDO

   end subroutine

   subroutine fill_the_blanks_()

      implicit none

      integer :: ipm

      DO ipm = 1, 2
         if(greenAA%compute(ipm, 3-ipm))then
            IF(.not.greenAA%compute(3-ipm, ipm)) CALL &
                 transform_correl(greenAA%correl(3-ipm, ipm), &
                 greenAA%correl(ipm, 3-ipm))
         endif
         if(greenBB%compute(ipm, 3-ipm))then
            IF(.not.greenBB%compute(3-ipm, ipm)) CALL &
                 transform_correl(greenBB%correl(3-ipm, ipm), &
                 greenBB%correl(ipm, 3-ipm))
         endif
      ENDDO

   end subroutine
