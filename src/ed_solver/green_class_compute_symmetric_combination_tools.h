   subroutine dump_output_test(greenAB, greenBA, BA)

      use common_def, only: find_rank
      use correl_class, only: transform_correl
      use green_class, only: green_type
      use mask_class, only: mask_type, new_mask

      TYPE(green_type), INTENT(INOUT) :: greenAB
      TYPE(green_type), INTENT(INOUT) :: greenBA
      integer                            :: iorb, jorb, ij, ji, ipm, si, ssi, &
                                            sssi, i, j, k, ssi_, si_, sssi_
      TYPE(mask_type)                    :: MASK(2, 2)
      complex(8)                         :: xxx(8), aaa
      LOGICAL                            :: BA

      CALL new_mask(MASK(1, 1), greenAB%correlstat(1, 1)%rc%MASK)
      CALL new_mask(MASK(1, 1), greenAB%correlstat(1, 2)%rc%MASK)
      if (BA) then
         write (96, *) 'MASKS BA --:', greenBA%correlstat(1, 1)%rc%MASK%mat
         write (96, *) 'MASKS BA + + :', greenBA%correlstat(2, 2)%rc%MASK%mat
      endif
      write (96, *) 'MASKS BA --:', greenAB%correlstat(1, 1)%rc%MASK%mat
      write (96, *) 'MASKS BA + + :', greenAB%correlstat(2, 2)%rc%MASK%mat
      if (BA) then
         write (96, *) 'MASKS BA - + :', greenBA%correlstat(1, 2)%rc%MASK%mat
         write (96, *) 'MASKS BA + -:', greenBA%correlstat(2, 1)%rc%MASK%mat
      endif
      write (96, *) 'MASKS AB - + :', greenAB%correlstat(1, 2)%rc%MASK%mat
      write (96, *) 'MASKS AB + -:', greenAB%correlstat(2, 1)%rc%MASK%mat

      if (.false.) then
         write (96, *) ' ============ SECTOR particle/particle ============= '
         DO iorb = 1, greenAB%N
            DO jorb = 1, greenAB%N
               ij = find_rank(MASK(1, 1)%imat(iorb, jorb), MASK(1, 1)%ivec)
               ji = find_rank(MASK(1, 1)%imat(jorb, iorb), MASK(1, 1)%ivec)
               write (96, *) 'link  : ', iorb, jorb
               write (96, *) ' ---------------AB----------------'
               write (96, *) 'ij -- : ', greenAB%correlstat(1, 1)%rc%vec(ij)
               write (96, *) 'ji -- : ', greenAB%correlstat(1, 1)%rc%vec(ji)
               write (96, *) 'ij + + : ', greenAB%correlstat(2, 2)%rc%vec(ij)
               write (96, *) 'ji + + : ', greenAB%correlstat(2, 2)%rc%vec(ji)
               if (BA) then
                  write (96, *) ' ---------------BA----------------'
                  write (96, *) 'ij -- : ', greenBA%correlstat(1, 2)%rc%vec(ij)
                  write (96, *) 'ji -- : ', greenBA%correlstat(1, 2)%rc%vec(ji)
                  write (96, *) 'ij + + : ', greenBA%correlstat(2, 1)%rc%vec(ij)
                  write (96, *) 'ji + + : ', greenBA%correlstat(2, 1)%rc%vec(ji)
               endif
               write (96, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
            ENDDO
         ENDDO
      endif

      write (96, *)
      write (96, *) ' ---------------- > NOW DYNAMIC CORRELATIONS < &
           &---------------------------'
      write (96, *) ' ---------------- > NOW DYNAMIC CORRELATIONS < &
           &---------------------------'
      write (96, *)

      DO iorb = 1, greenAB%N
         DO jorb = 1, greenAB%N
            ij = find_rank(MASK(1, 1)%imat(iorb, jorb), MASK(1, 1)%ivec)
            ji = find_rank(MASK(1, 1)%imat(jorb, iorb), MASK(1, 1)%ivec)
            write (96, *) 'link  : ', iorb, jorb
            write (96, *) ' ---------------AB----------------'
            write (96, *) 'ij -- : ', greenAB%correl(1, 1)%vec(ij, 4)
            write (96, *) 'ji -- : ', greenAB%correl(1, 1)%vec(ji, 4)
            write (96, *) 'ij + + : ', greenAB%correl(2, 2)%vec(ij, 4)
            write (96, *) 'ji + + : ', greenAB%correl(2, 2)%vec(ji, 4)
            if (BA) then
               write (96, *) ' ---------------BA----------------'
               write (96, *) 'ij -- : ', greenBA%correl(1, 1)%vec(ij, 4)
               write (96, *) 'ji -- : ', greenBA%correl(1, 1)%vec(ji, 4)
               write (96, *) 'ij + + : ', greenBA%correl(2, 2)%vec(ij, 4)
               write (96, *) 'ji + + : ', greenBA%correl(2, 2)%vec(ji, 4)
            endif
            write (96, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

            write (97, *) 'link  : ', iorb, jorb
            write (97, *) ' ---------------AB----------------'
            write (97, *) 'ij -- : ', greenAB%correl(1, 1)%vec(ij, 600)
            write (97, *) 'ji -- : ', greenAB%correl(1, 1)%vec(ji, 600)
            write (97, *) 'ij + + : ', greenAB%correl(2, 2)%vec(ij, 600)
            write (97, *) 'ji + + : ', greenAB%correl(2, 2)%vec(ji, 600)
            if (BA) then
               write (97, *) ' ---------------BA----------------'
               write (97, *) 'ij -- : ', greenBA%correl(1, 1)%vec(ij, 600)
               write (97, *) 'ji -- : ', greenBA%correl(1, 1)%vec(ji, 600)
               write (97, *) 'ij + + : ', greenBA%correl(2, 2)%vec(ij, 600)
               write (97, *) 'ji + + : ', greenBA%correl(2, 2)%vec(ji, 600)
            endif
            write (97, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

            if (.false.) then
               xxx(1) = (greenAB%correl(1, 1)%vec(ij, 4))
               xxx(2) = (greenAB%correl(1, 1)%vec(ji, 4))
               xxx(3) = (greenAB%correl(2, 2)%vec(ij, 4))
               xxx(4) = (greenAB%correl(2, 2)%vec(ji, 4))
               if (BA) then
                  xxx(5) = (greenBA%correl(1, 1)%vec(ij, 4))
                  xxx(6) = (greenBA%correl(1, 1)%vec(ji, 4))
                  xxx(7) = (greenBA%correl(2, 2)%vec(ij, 4))
                  xxx(8) = (greenBA%correl(2, 2)%vec(ji, 4))
               endif
               do i = 1, 8
                  do j = 1, 8
                     do k = 1, 8
                        if (i /= j .and. j /= k) then
                           if (abs(abs(xxx(i)) - abs(xxx(j))) > 1.d-3 .and. &
                               abs(abs(xxx(j)) - abs(xxx(k))) > 1.d-3 .and. &
                               abs(abs(xxx(i)) - abs(xxx(k))) > 1.d-3) then
                              do si = -2, 2
                                 do ssi = -2, 2
                                    do sssi = -2, 2
                                       do si_ = -2, 2
                                          do ssi_ = -2, 2
                                             do sssi_ = -2, 2
                                                aaa = abs(cmplx(dble(si), &
                                                                dble(si_))*xxx(i) + &
                                                          cmplx(dble(ssi), &
                                                                dble(ssi_))*xxx(j) - &
                                                          cmplx(dble(sssi), &
                                                                dble(sssi_))*xxx(k))
                                                if (abs(aaa) < 1.d-5 .and. &
                                                    abs(si) + abs(ssi) + &
                                                    abs(sssi) > 0) then
                                                   write (96, *) '-------------------------------------'
                                                   write (96, *) 'elements are &
                                                        &equal : ', i, j, k
                                                   write (96, *)
                                                   write (96, *) 'x1 : ', &
                                                      xxx(i), si, si_
                                                   write (96, *) ' + x2 : ', &
                                                      xxx(j), ssi, ssi_
                                                   write (96, *) ' = x3 : ', &
                                                      xxx(k), sssi, sssi_
                                                   write (96, *) 'difference : &
                                                        &', aaa
                                                endif
                                             enddo
                                          enddo
                                       enddo
                                    enddo
                                 enddo
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               enddo
            endif

            if (.true.) then
               do ipm = 1, 2
                  if (BA) CALL transform_correl(greenBA%correl(ipm, ipm), &
                                                greenBA%correl(ipm, ipm))
                  CALL transform_correl(greenAB%correl(ipm, ipm), &
                                        greenAB%correl(ipm, ipm))
               enddo
               write (96, *) 'link  : ', iorb, jorb
               write (96, *) ' ---------------AB----------------'
               write (96, *) 'ij -- : ', greenAB%correl(1, 1)%vec(ij, 4)
               write (96, *) 'ji -- : ', greenAB%correl(1, 1)%vec(ji, 4)
               write (96, *) 'ij + + : ', greenAB%correl(2, 2)%vec(ij, 4)
               write (96, *) 'ji + + : ', greenAB%correl(2, 2)%vec(ji, 4)
               if (BA) then
                  write (96, *) ' ---------------BA----------------'
                  write (96, *) 'ij -- : ', greenBA%correl(1, 1)%vec(ij, 4)
                  write (96, *) 'ji -- : ', greenBA%correl(1, 1)%vec(ji, 4)
                  write (96, *) 'ij + + : ', greenBA%correl(2, 2)%vec(ij, 4)
                  write (96, *) 'ji + + : ', greenBA%correl(2, 2)%vec(ji, 4)
               endif
               write (96, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

               write (97, *) 'link  : ', iorb, jorb
               write (97, *) ' ---------------AB----------------'
               write (97, *) 'ij -- : ', greenAB%correl(1, 1)%vec(ij, 600)
               write (97, *) 'ji -- : ', greenAB%correl(1, 1)%vec(ji, 600)
               write (97, *) 'ij + + : ', greenAB%correl(2, 2)%vec(ij, 600)
               write (97, *) 'ji + + : ', greenAB%correl(2, 2)%vec(ji, 600)
               if (BA) then
                  write (97, *) ' ---------------BA----------------'
                  write (97, *) 'ij -- : ', greenBA%correl(1, 1)%vec(ij, 600)
                  write (97, *) 'ji -- : ', greenBA%correl(1, 1)%vec(ji, 600)
                  write (97, *) 'ij + + : ', greenBA%correl(2, 2)%vec(ij, 600)
                  write (97, *) 'ji + + : ', greenBA%correl(2, 2)%vec(ji, 600)
               endif
               write (97, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

               do ipm = 1, 2
                  if (BA) CALL transform_correl(greenBA%correl(ipm, ipm), &
                                                greenBA%correl(ipm, ipm))
                  CALL transform_correl(greenAB%correl(ipm, ipm), &
                                        greenAB%correl(ipm, ipm))
               enddo
            endif

         ENDDO
      ENDDO

   END SUBROUTINE
