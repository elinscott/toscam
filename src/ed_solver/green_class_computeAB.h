   subroutine scan_indices()

      use globalvar_ed_solver, only: Neigen

      implicit none

      if (allocated(indices_state_sector)) deallocate (indices_state_sector)
      allocate (indices_state_sector(GS%nsector*2*2*2*Neigen*greenAB%N*greenAB%N, &
                                     7))
      ktot = 0
      indices_state_sector = 0
      DO isector = 1, GS%nsector
         DO ipm = 1, 2
            DO jpm = 1, 2
               IF (greenAB%compute(ipm, jpm)) then
                  IF (associated(greenAB%correl(ipm, jpm)%MM%MASK%mat)) then
                     IF (ANY(greenAB%correl(ipm, jpm)%MM%MASK%mat)) then
                        DO iph = 1, 2
                           DO ieigen = 1, GS%es(isector)%lowest%neigen
                              iorb_f = 0
                              jorb_f = 0
                              do iiorb = 1, greenAB%N
                                 do jjorb = 1, greenAB%N
                                    IF (MASKAB(ipm, jpm)%mat(iiorb, jjorb)) THEN
                                       if (iorb_f == 0 .and. jorb_f == 0) then
                                          iorb_f = iiorb
                                          jorb_f = jjorb
                                       endif
                                       ktot = ktot + 1
                                       if (messages3 .and. rank == 0) write (*, &
                                                                             *) 'SECTOR TO CONSIDER IN MPI : ', &
                                          isector, ipm, jpm, iph, ieigen, &
                                          iiorb, jjorb
                                       call save_indices
                                       if (ktot + 1 > &
                                           size(indices_state_sector, 1)) &
                                          then
                                          write (*, *) 'MPI_GREENS : exceed &
                                               &size of indices'
                                          stop 'critical'
                                       endif
                                    endif
                                 enddo
                              enddo
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   end subroutine

   subroutine save_indices()

      implicit none

      indices_state_sector(ktot, 1) = isector
      indices_state_sector(ktot, 2) = ipm
      indices_state_sector(ktot, 3) = jpm
      indices_state_sector(ktot, 4) = iph
      indices_state_sector(ktot, 5) = ieigen
      indices_state_sector(ktot, 6) = iiorb
      indices_state_sector(ktot, 7) = jjorb
   end subroutine

   subroutine fix_indices()

      implicit none

      isector = indices_state_sector(kk, 1)
      ipm = indices_state_sector(kk, 2)
      jpm = indices_state_sector(kk, 3)
      iph = indices_state_sector(kk, 4)
      iieigen = indices_state_sector(kk, 5)
      iiorb = indices_state_sector(kk, 6)
      jjorb = indices_state_sector(kk, 7)
   end subroutine

   subroutine do_ba()

      use globalvar_ed_solver, only: dEmax0, FLAG_FULL_ED_GREEN
      use H_class, only: delete_H
      use linalg, only: dexpc

      implicit none

      if (BA) then
         call build_H_ba
         DO ieigen = 1, es%lowest%neigen
            if (FLAG_MPI_GREENS /= 2 .or. ieigen == iieigen) then
               if (FLAG_FULL_ED_GREEN) then
                  if (abs(beta*(es%lowest%eigen(ieigen)%val - E0)) > dEmax0) goto &
                     38
               endif
               eigen => es%lowest%eigen(ieigen)
               boltz = DEXPc(-beta*(eigen%val - E0))/Zpart
               call apply_creatba
               call greenba_
38             continue
            endif
         ENDDO
         IF (compute_dyn_correl) CALL delete_H()
      endif
   end subroutine

   subroutine do_ab()

      use globalvar_ed_solver, only: dEmax0, FLAG_FULL_ED_GREEN
      use H_class, only: delete_H
      use linalg, only: dexpc

      implicit none

      call build_H_ab()
      DO ieigen = 1, es%lowest%neigen
         if (FLAG_MPI_GREENS /= 2 .or. ieigen == iieigen) then
            if (FLAG_FULL_ED_GREEN) then
               if (abs(beta*(es%lowest%eigen(ieigen)%val - E0)) > dEmax0) goto 38
            endif
            eigen => es%lowest%eigen(ieigen)
            boltz = DEXPc(-beta*(eigen%val - E0))/Zpart
            call apply_creatab()
            call greenab_()
38          continue
         endif
      ENDDO
      IF (compute_dyn_correl) CALL delete_H()
   end subroutine

   subroutine build_H_ba()

      use common_def, only: dump_message
      use H_class, only: new_H

      implicit none

      IF (compute_dyn_correl) THEN
         SELECT CASE (iph)
         CASE (1)
            CALL dump_message(TEXT="### PARTICLE PART")
            CALL new_H(AIM, Bpm_es(3 - ipm)%sector)
         CASE (2)
            CALL dump_message(TEXT="### HOLE PART")
            CALL new_H(AIM, Bpm_es(ipm)%sector)
         END SELECT
      ENDIF
      SELECT CASE (iph)
      CASE (1)
         if (associated(Bpm_es(3 - ipm)%sector%sz)) isec_back = &
            Bpm_es(3 - ipm)%sector%sz%dimen
      CASE (2)
         if (associated(Bpm_es(ipm)%sector%sz)) isec_back = &
            Bpm_es(ipm)%sector%sz%dimen
      END SELECT
      write (*, *) 'build H BA, sector size : ', isec_back
   end subroutine

   subroutine build_H_ab()

      use common_def, only: dump_message
      use H_class, only: new_H

      implicit none

      IF (compute_dyn_correl) THEN
         SELECT CASE (iph)
         CASE (1)
            CALL dump_message(TEXT="### PARTICLE PART")
            CALL new_H(AIM, Apm_es(3 - ipm)%sector)
         CASE (2)
            CALL dump_message(TEXT="### HOLE PART")
            CALL new_H(AIM, Apm_es(ipm)%sector)
         END SELECT
      ENDIF
      SELECT CASE (iph)
      CASE (1)
         if (associated(Apm_es(3 - ipm)%sector%sz)) isec_back = &
            Apm_es(3 - ipm)%sector%sz%dimen
      CASE (2)
         if (associated(Apm_es(ipm)%sector%sz)) isec_back = &
            Apm_es(ipm)%sector%sz%dimen
      END SELECT
      write (*, *) 'build H AB, sector size for A : ', isec_back
   end subroutine

   subroutine print_label()

      use genvar, only: log_unit
      use matrix, only: write_array

      implicit none

      write (log_unit, *) ' &
           &%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write (log_unit, *) ' ------ > PART/HOLE CONTRIB   : ', iph
      write (log_unit, *) ' ------ > C C                 : ', ipm, jpm
      write (log_unit, *) ' ------ >  up, dn, tot sites : ', uup, ddn, itot
      write (log_unit, *) ' ------ >  Sz                 : ', ssz
      write (log_unit, *) ' ------ > neigen              : ', es%lowest%neigen
      write (log_unit, *) ' ------ > dim space of the GS : ', &
         es%lowest%eigen(:)%dim_space
      write (log_unit, *) ' ------ > Energies : ', es%lowest%eigen(:)%val
      write (log_unit, *) ' ------ > MASK greenAB        : '
      call write_array(greenAB%correl(ipm, jpm)%MM%MASK%mat, ' green AB ')
      write (log_unit, *) ' ------ > MASKAB : ', MASKAB(ipm, jpm)%ivec
      if (BA) then
         write (log_unit, *) ' ------ > MASK greenBA        : '
         call write_array(greenBA%correl(ipm, jpm)%MM%MASK%mat, ' green BA ')
      endif
      write (log_unit, *) ' &
           &%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   end subroutine

   subroutine init_data()

      use common_def, only: dump_message, find_rank, reset_timer
      use genvar, only: log_unit
      use mask_class, only: new_mask

      implicit none

      write (*, *) 'init_data'
      write (*, *) 'part1'
      IF (BA) THEN
         IF (greenBA%Nw /= greenAB%Nw .OR. greenAA%Nw /= greenAB%Nw .OR. &
             greenBB%Nw /= greenAB%Nw) then
            write (*, *) "ERROR IN compute_greenAB: INCONSISTENT Nw!"
            STOP
         ENDIF
         IF (greenBA%N /= greenAB%N .OR. greenAA%N /= greenAB%N .OR. greenBB%N &
             /= greenAB%N) then
            write (*, *) "ERROR IN compute_greenAB: INCONSISTENT N!"
            STOP
         ENDIF
      ENDIF
      CALL reset_timer(compute_green_timer)
      IF (BA) THEN
         CALL dump_message(TEXT="### START COMPUTING AB "// &
              TRIM(greenAB%title)//" AND "//TRIM(greenBA%title)//" GREEN &
              &S FUNCTIONS")
      ELSE
         CALL dump_message(TEXT="### START COMPUTING AB "// &
                           TRIM(greenAB%title)//" GREEN S FUNCTION")
      ENDIF

      ORBMASKvec = .false.
      write (*, *) 'part2'
      DO ipm = 1, 2
         DO jpm = 1, 2
            ! EQUAL-TIME
            write (*, *) 'building equal time', ipm, jpm
            CALL new_mask(MASKAB(ipm, jpm), greenAB%correlstat(ipm, &
                                                               jpm)%rc%MASK)
            write (*, *) 'new mask built'
            if (BA) then
               write (*, *) 'building BA mask'
               CALL new_mask(MASKBA(ipm, jpm), greenBA%correlstat(ipm, &
                                                                  jpm)%rc%MASK)
            endif
            ! VECTOR MASK OF ORBITALS: MASK(i) = T IF NEED TO APPLY C(i)
            write (*, *) 'building vector mask'
            DO iorb = 1, greenAB%N
               IF (ANY(MASKAB(ipm, jpm)%mat(iorb, :)) .OR. ANY(MASKAB(ipm, &
                                                                      jpm)%mat(:, iorb))) ORBMASKvec(ipm, jpm, iorb) = .true.
            ENDDO
         ENDDO
      ENDDO

      !-----------------------------!
      ! SET MATRIX ELEMENTS TO ZERO !
      !-----------------------------!
      DO ipm = 1, 2
         DO jpm = 1, 2
            DO iorb = 1, greenAB%N
               DO jorb = 1, greenAB%N

                  IF (MASKAB(ipm, jpm)%mat(iorb, jorb)) THEN
                     iind = find_rank(MASKAB(ipm, jpm)%imat(iorb, jorb), &
                                      MASKAB(ipm, jpm)%ivec)
                     if (greenAB%compute(ipm, jpm)) then
                        greenAB%correlstat(ipm, jpm)%rc%vec(iind) = 0.0_DBL
                        greenAB%correl(ipm, jpm)%vec(iind, :) = 0.0_DBL
                     endif
                  ENDIF
                  if (BA) then
                     IF (MASKBA(ipm, jpm)%mat(iorb, jorb)) THEN
                        iind = find_rank(MASKBA(ipm, jpm)%imat(iorb, jorb), &
                                         MASKBA(ipm, jpm)%ivec)
                        if (greenBA%compute(ipm, jpm)) then
                           greenBA%correlstat(ipm, jpm)%rc%vec(iind) = 0.0_DBL
                           greenBA%correl(ipm, jpm)%vec(iind, :) = 0.0_DBL
                        endif
                     ENDIF

                  endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO iorb = 1, greenAB%N

         IF (ANY(ORBMASKvec(:, :, iorb))) THEN
            if (associated(greenAB%Amean)) then
               greenAB%Amean(iorb, :) = 0.0_DBL
               greenAB%Bmean(iorb, :) = 0.0_DBL
            endif
            if (associated(greenBA%Amean)) then
               greenBA%Amean(iorb, :) = 0.0_DBL
               greenBA%Bmean(iorb, :) = 0.0_DBL
            endif
         endif
      ENDDO

   end subroutine

   subroutine init_sector()

      use eigen_sector_class, only: new_eigensector, &
         not_commensurate_sector_
      use genvar, only: log_unit, pm
      use sector_class, only: equal_sector

      implicit none

      es => GS%es(isector)

      if (AIM%BATH%SUPER) then
         uup = -1
         ddn = -1
         ssz = es%sector%sz%npart
      else
         ssz = 0
         uup = es%sector%updo%up%npart
         ddn = es%sector%updo%down%npart
      endif

      write (log_unit, *) '-----------------------------------------------------------------------------'
      write (log_unit, *) '------ > PARSE EIGEN SECTORS : ', isector, GS%nsector
      write (log_unit, *) ' ...... eigenvalues ...... : ', &
         es%lowest%eigen(:)%val
      write (log_unit, *) ' ...... lowest neig. ...... : ', es%lowest%neigen
      write (log_unit, *) ' ...... Sz           ...... : ', ssz
      write (log_unit, *) ' ...... eigenvector ...... : ', (/(maxval(abs( &
                                                                     es%lowest%eigen(ii)%vec%rc)), ii=1, es%lowest%neigen)/)
      write (log_unit, *) '-----------------------------------------------------------------------------'

      if (es%lowest%neigen == 0) then
         write (*, *) 'ERROR in green_computeAB : es%lowest%neigen = 0'
         stop
      endif
      if (associated(GS%es(isector)%lowest%eigen)) then
         if (GS%es(isector)%lowest%neigen == 0) then
            write (*, *) 'ERROR in green_computeAB : GS%es%lowest%neigen = 0'
            stop
         endif
      endif

      NOT_COMMENSURATE = .false.
      ! FIRST WE CREATE THE TARGET SECTORS
      DO iph = 1, 2
         CALL Asector(Asec, pm(iph), es%sector)
         CALL new_eigensector(Apm_es(iph), Asec)
         if (not_commensurate_sector_(Asec)) NOT_COMMENSURATE = .true.
         CALL Bsector(Bsec, pm(iph), es%sector)
         CALL new_eigensector(Bpm_es(iph), Bsec)
         if (not_commensurate_sector_(Bsec)) NOT_COMMENSURATE = .true.
      ENDDO

      ! CONSISTENCY CHECK
      DO ipm = 1, 2
         DO jpm = 1, 2
            IF (greenAB%compute(ipm, jpm)) THEN
               IF (.NOT. equal_sector(Apm_es(3 - ipm)%sector, Bpm_es( &
                                      jpm)%sector) .OR. .NOT. equal_sector(Apm_es(ipm)%sector, &
                                                                           Bpm_es(3 - jpm)%sector)) then
                  write (*, *) "ERROR IN green_class_computeAB: INCONSISTENT &
                       &SECTORS!"
                  stop
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   end subroutine

   subroutine clean_everything()

      use common_def, only: find_rank, timer_fortran
      use correl_class, only: vec2correl
      use eigen_class, only: delete_eigen
      use eigen_sector_class, only: delete_eigensector
      use globalvar_ed_solver, only: donot_compute_holepart
      use mask_class, only: delete_mask
      use masked_matrix_class, only: vec2masked_matrix
      use sector_class, only: delete_sector

      implicit none

      if (donot_compute_holepart) then
         DO iorb = 1, greenAB%N
            DO jorb = 1, greenAB%N
               IF (MASKAB(1, 2)%mat(iorb, jorb)) THEN
                  iind = find_rank(MASKAB(1, 2)%imat(iorb, jorb), MASKAB(1, &
                                                                         2)%ivec)
                  if (maxval(abs(greenAB%correlstat(1, 2)%rc%vec(:))) < 1.d-15) &
                     greenAB%correlstat(1, 2)%rc%vec(iind) = &
                     1.d0 - greenAB%correlstat(2, 1)%rc%vec(iind)
                  if (maxval(abs(greenAB%correlstat(2, 1)%rc%vec(:))) < 1.d-15) &
                     greenAB%correlstat(2, 1)%rc%vec(iind) = &
                     1.d0 - greenAB%correlstat(1, 2)%rc%vec(iind)
               ENDIF
               if (BA) then
                  IF (MASKBA(1, 2)%mat(iorb, jorb)) THEN
                     iind = find_rank(MASKBA(1, 2)%imat(iorb, jorb), MASKBA(1, &
                                                                            2)%ivec)
                     if (maxval(abs(greenBA%correlstat(1, 2)%rc%vec(:))) < &
                         1.d-15) greenBA%correlstat(1, 2)%rc%vec(iind) = &
                        1.d0 - greenBA%correlstat(2, 1)%rc%vec(iind)
                     if (maxval(abs(greenBA%correlstat(2, 1)%rc%vec(:))) < &
                         1.d-15) greenBA%correlstat(2, 1)%rc%vec(iind) = &
                        1.d0 - greenBA%correlstat(1, 2)%rc%vec(iind)
                  ENDIF
               endif
            ENDDO
         ENDDO
      endif

      !-----------------!
      ! EXPAND MATRICES !
      !-----------------!

      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL vec2masked_matrix(greenAB%correlstat(ipm, jpm))
            IF (BA) CALL vec2masked_matrix(greenBA%correlstat(ipm, jpm))
            IF (compute_dyn_correl .AND. (greenAB%compute(ipm, jpm) .OR. &
                                          greenAB%compute(3 - ipm, 3 - jpm))) THEN
               CALL vec2correl(greenAB%correl(ipm, jpm))
               IF (BA) CALL vec2correl(greenBA%correl(ipm, jpm))
            ENDIF
         ENDDO
      ENDDO

      !----------!
      ! CLEAN UP !
      !----------!

      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL delete_mask(MASKAB(ipm, jpm))
            CALL delete_mask(MASKBA(ipm, jpm))
         ENDDO
      ENDDO
      CALL delete_sector(Asec)
      CALL delete_sector(Bsec)
      DO iph = 1, 2
         CALL delete_eigensector(Apm_es(iph))
         CALL delete_eigensector(Bpm_es(iph))
      ENDDO
      CALL delete_eigen(ABpmsym)

      IF (BA) THEN
         CALL timer_fortran(compute_green_timer, "### COMPUTING "// &
              TRIM(greenAB%title)//" AND "//TRIM(greenBA%title)//" &
              &GREEN'S FUNCTIONS TOOK")
      ELSE
         CALL timer_fortran(compute_green_timer, "### COMPUTING "// &
                            TRIM(greenAB%title)//" GREEN'S FUNCTION TOOK")
      ENDIF

   end subroutine

   subroutine greenba_()

      use common_def, only: c2s, find_rank, i2c, &
         reset_timer, timer_fortran
      use genvar, only: log_unit
      use linalg, only: k_to_ij
      use green_class_compute_symmetric, only: symmetric_combineAB
      use green_class_compute_dynamic, only: compute_dynamic
      use rcvector_class, only: norm_rcvector

      implicit none

      real(8) :: normvec
      integer :: iisector, i_, v(2), start, step

      write (log_unit, *) '---------------------------------'
      write (log_unit, *) '......calculating BA terms.......'
      write (log_unit, *) ' BOLTZMAN WEIGHT      : ', boltz
      write (log_unit, *) '---------------------------------'
      write (log_unit, *) '---------------------------------'

      if (FLAG_MPI_GREENS == 1) then
         start = rank
         step = size2
      else
         start = 0
         step = 1
      endif

      DO i_ = start + 1, greenBA%N**2, step
         v = k_to_ij(greenBA%N, i_)
         iorb = v(1)
         jorb = v(2)

         if (.not. (FLAG_MPI_GREENS == 2) .or. (iorb == iiorb .and. jorb == &
                                                jjorb)) then
            IF (MASKBA(ipm, jpm)%mat(iorb, jorb)) THEN

               iind = find_rank(MASKBA(ipm, jpm)%imat(iorb, jorb), MASKBA(ipm, &
                                                                          jpm)%ivec)
               if (iind > 0) then

                  write (log_unit, *) '....orbitals....', iorb, jorb, iind
                  write (log_unit, *) 'how many elements in MASK ---- > ', &
                     size(MASKBA(ipm, jpm)%ivec)

                  CALL symmetric_combineAB(ABpmsym, iorb, jorb, ipm, jpm, iph, &
                                           Bpm_es, Apm_es, BA=.true., isec_back=isec_back, iisector &
                                           =iisector, GS=GS, ssz=issz)

                  normvec = 0.d0

                  ! EQUAL-TIME CORRELATIONS
                  IF (iph == 1) THEN
                     normvec = norm_rcvector(ABpmsym%vec)
                     greenBA%correlstat(ipm, jpm)%rc%vec(iind) = &
                        greenBA%correlstat(ipm, jpm)%rc%vec(iind) + &
                        boltz*normvec**2
                  ENDIF

                  ! DYNAMIC  CORRELATIONS
                  IF (compute_dyn_correl .AND. greenBA%compute(ipm, jpm)) THEN
                     CALL reset_timer(compute_dyn_timer)
                     CALL compute_dynamic(iph, dyn, greenBA%correl(ipm, &
                                                                   jpm)%freq, ABpmsym, greenBA%correl(ipm, jpm)%stat, &
                                          greenBA%title, normvec, iisector, GS)
                     greenBA%correl(ipm, jpm)%vec(iind, :) = &
                        greenBA%correl(ipm, jpm)%vec(iind, :) + boltz*dyn
                     CALL timer_fortran(compute_dyn_timer, "# COMPUTING MATRIX &
                          &ELEMENT BA "//c2s(i2c(iind))//"/"// &
                          c2s(i2c(MASKBA(ipm, jpm)%nind))//" OF "// &
                          TRIM(greenBA%title)//" TOOK ")
                  ENDIF

               endif
            endif
         endif
      ENDDO

   end subroutine

   subroutine greenab_()

      use common_def, only: c2s, find_rank, i2c, &
         reset_timer, timer_fortran
      use genvar, only: log_unit, pm
      use green_class_compute_dynamic, only: compute_dynamic
      use green_class_compute_symmetric, only: symmetric_combineAB
      use linalg, only: k_to_ij
      use rcvector_class, only: norm_rcvector

      implicit none

      real(8) :: normvec
      integer :: iisector, i_, v(2), start, step

      ! THEN WE FORM THEIR SYMMETRIC COMBINATIONS & COMPUTE CORRELATIONS FIRST
      ! THE AB CASE !

      write (log_unit, *) '---------------------------------'
      write (log_unit, *) '......calculating AB terms.......'
      write (log_unit, *) ' BOLTZMAN WEIGHT      : ', boltz
      write (log_unit, *) '---------------------------------'
      write (log_unit, *) '---------------------------------'

      if (FLAG_MPI_GREENS == 1) then
         start = rank
         step = size2
      else
         start = 0
         step = 1
      endif

      DO i_ = start + 1, greenAB%N**2, step
         v = k_to_ij(greenAB%N, i_)
         iorb = v(1)
         jorb = v(2)

         if (.not. (FLAG_MPI_GREENS == 2) .or. (iorb == iiorb .and. jorb == &
                                                jjorb)) then
            IF (MASKAB(ipm, jpm)%mat(iorb, jorb)) THEN

               iind = find_rank(MASKAB(ipm, jpm)%imat(iorb, jorb), MASKAB(ipm, &
                                                                          jpm)%ivec)
               if (iind > 0) then

                  CALL symmetric_combineAB(ABpmsym, iorb, jorb, ipm, jpm, iph, &
                                           Apm_es, Bpm_es, BA=.false., isec_back=isec_back, iisector &
                                           =iisector, GS=GS, ssz=issz)

                  write (log_unit, *) 'iorb, jorb, max AB symmetrized :', &
                     maxval(abs(ABpmsym%vec%rc))

                  normvec = 0.d0

                  ! EQUAL-TIME CORRELATIONS
                  IF (iph == 1) THEN
                     normvec = norm_rcvector(ABpmsym%vec)
                     greenAB%correlstat(ipm, jpm)%rc%vec(iind) = &
                        greenAB%correlstat(ipm, jpm)%rc%vec(iind) + &
                        boltz*normvec**2
                  ENDIF

                  ! DYNAMIC  CORRELATIONS
                  IF (compute_dyn_correl .AND. greenAB%compute(ipm, jpm)) THEN
                     CALL reset_timer(compute_dyn_timer)
                     CALL compute_dynamic(iph, dyn, greenAB%correl(ipm, &
                                                                   jpm)%freq, ABpmsym, greenAB%correl(ipm, jpm)%stat, &
                                          greenAB%title, normvec, iisector, GS)
                     greenAB%correl(ipm, jpm)%vec(iind, :) = &
                        greenAB%correl(ipm, jpm)%vec(iind, :) + boltz*dyn
                     CALL timer_fortran(compute_dyn_timer, "# COMPUTING MATRIX &
                          &ELEMENT "//pm(iph)//c2s(i2c(iind))//"/"// &
                          c2s(i2c(MASKAB(ipm, jpm)%nind))//" OF "// &
                          TRIM(greenAB%title)//" TOOK ")
                  ENDIF

               endif
            endif
         endif
      ENDDO

      return
   end subroutine

   subroutine apply_creatab()

      use common_def, only: reset_timer, timer_fortran
      use eigen_class, only: delete_eigenlist, rank_eigen_in_list
      use genvar, only: log_unit, logfile, pm
      use mpi_mod, only: MPI_DOT_PRODUCT
      use sector_class, only: equal_sector

      implicit none

      !------------------------------------------------------!
      ! APPLY CREATION/ANNIHILATION OPERATORS ON THE STATES  !
      !------------------------------------------------------!

      write (log_unit, *) 'apply operator on the state in this sector, on &
           &orbitals : ', ORBMASKvec(ipm, jpm, :)

      DO kpm = 1, 2
         CALL reset_timer(apply_timer)
         CALL delete_eigenlist(Apm_es(kpm)%lowest)
         CALL applyA(Apm_es(kpm), pm(kpm), ORBMASKvec(ipm, jpm, :), es, ieigen)
         CALL timer_fortran(apply_timer, "# APPLYING A"//pm(kpm)//" ON &
              &EIGENSTATE TOOK")
         write (logfile, *) 'applying A gives xx states : ', &
            size(Apm_es(kpm)%lowest%eigen)

         CALL reset_timer(apply_timer)
         CALL delete_eigenlist(Bpm_es(kpm)%lowest)
         CALL applyB(Bpm_es(kpm), pm(kpm), ORBMASKvec(ipm, jpm, :), es, ieigen)
         CALL timer_fortran(apply_timer, "# APPLYING B"//pm(kpm)//" ON &
              &EIGENSTATE TOOK")
         write (logfile, *) 'applying B gives xx states : ', &
            size(Bpm_es(kpm)%lowest%eigen)

         if (.not. (FLAG_MPI_GREENS == 2) .or. (iiorb == iorb_f .and. jjorb == &
                                                jorb_f)) then
            IF (equal_sector(es%sector, Apm_es(kpm)%sector)) THEN
               write (log_unit, *) '.... SECTOR IS STABLE UNDER OPERATOR A ....'
               DO iorb = 1, greenAB%N
                  IF (ORBMASKvec(ipm, jpm, iorb)) THEN
                     OPeigen => &
                        Apm_es(kpm)%lowest%eigen(rank_eigen_in_list(iorb, &
                                                                    Apm_es(kpm)%lowest))
                     greenAB%Amean(iorb, kpm) = greenAB%Amean(iorb, kpm) + &
                                                boltz*MPI_DOT_PRODUCT(eigen%vec%rc, &
                                                                      OPeigen%vec%rc, split=USE_TRANSPOSE_TRICK_MPI)
                  ENDIF
               ENDDO
            ELSE
               write (log_unit, *) '.... SECTOR IS NOT STABLE UNDER OPERATOR A &
                    &....'
            ENDIF
            IF (equal_sector(es%sector, Bpm_es(kpm)%sector)) THEN
               write (log_unit, *) '.... SECTOR IS STABLE UNDER OPERATOR B ....'
               DO iorb = 1, greenAB%N
                  IF (ORBMASKvec(ipm, jpm, iorb)) THEN
                     OPeigen => &
                        Bpm_es(kpm)%lowest%eigen(rank_eigen_in_list(iorb, &
                                                                    Bpm_es(kpm)%lowest))
                     greenAB%Bmean(iorb, kpm) = greenAB%Bmean(iorb, kpm) + &
                                                boltz*MPI_DOT_PRODUCT(eigen%vec%rc, &
                                                                      OPeigen%vec%rc, split=USE_TRANSPOSE_TRICK_MPI)
                  ENDIF
               ENDDO
            ELSE
               write (log_unit, *) '.... SECTOR IS NOT STABLE UNDER OPERATOR B &
                    &....'
            ENDIF
         endif

      ENDDO

   end subroutine

   subroutine apply_creatba()

      use common_def, only: reset_timer, timer_fortran
      use eigen_class, only: delete_eigenlist, rank_eigen_in_list
      use genvar, only: log_unit, logfile, pm
      use mpi_mod, only: MPI_DOT_PRODUCT
      use sector_class, only: equal_sector

      implicit none

      !------------------------------------------------------!
      ! APPLY CREATION/ANNIHILATION OPERATORS ON THE STATES  !
      !------------------------------------------------------!

      write (log_unit, *) 'apply operator on the state in this sector, on &
           &orbitals : ', ORBMASKvec(ipm, jpm, :)

      DO kpm = 1, 2
         CALL reset_timer(apply_timer)
         CALL delete_eigenlist(Apm_es(kpm)%lowest)
         CALL applyA(Apm_es(kpm), pm(kpm), ORBMASKvec(ipm, jpm, :), es, ieigen)
         CALL timer_fortran(apply_timer, "# APPLYING A"//pm(kpm)//" ON &
              &EIGENSTATE TOOK")
         write (logfile, *) 'applying A gives xx states : ', &
            size(Apm_es(kpm)%lowest%eigen)

         CALL reset_timer(apply_timer)
         CALL delete_eigenlist(Bpm_es(kpm)%lowest)
         CALL applyB(Bpm_es(kpm), pm(kpm), ORBMASKvec(ipm, jpm, :), es, ieigen)
         CALL timer_fortran(apply_timer, "# APPLYING B"//pm(kpm)//" ON &
              &EIGENSTATE TOOK")
         write (logfile, *) 'applying B gives xx states : ', &
            size(Bpm_es(kpm)%lowest%eigen)

         if (.not. (FLAG_MPI_GREENS == 2) .or. (iiorb == iorb_f .and. jjorb == &
                                                jorb_f)) then
            IF (equal_sector(es%sector, Apm_es(kpm)%sector)) THEN
               write (log_unit, *) '.... SECTOR IS STABLE UNDER OPERATOR A ....'
               DO iorb = 1, greenAB%N
                  IF (ORBMASKvec(ipm, jpm, iorb)) THEN
                     OPeigen => &
                        Apm_es(kpm)%lowest%eigen(rank_eigen_in_list(iorb, &
                                                                    Apm_es(kpm)%lowest))
                     greenBA%Amean(iorb, kpm) = greenBA%Amean(iorb, kpm) + &
                                                boltz*MPI_DOT_PRODUCT(eigen%vec%rc, &
                                                                      OPeigen%vec%rc, split=USE_TRANSPOSE_TRICK_MPI)
                  ENDIF
               ENDDO
            ELSE
               write (log_unit, *) '.... SECTOR IS NOT STABLE UNDER OPERATOR A &
                    &....'
            ENDIF
            IF (equal_sector(es%sector, Bpm_es(kpm)%sector)) THEN
               write (log_unit, *) '.... SECTOR IS STABLE UNDER OPERATOR B ....'
               DO iorb = 1, greenAB%N
                  IF (ORBMASKvec(ipm, jpm, iorb)) THEN
                     OPeigen => &
                        Bpm_es(kpm)%lowest%eigen(rank_eigen_in_list(iorb, &
                                                                    Bpm_es(kpm)%lowest))
                     greenBA%Bmean(iorb, kpm) = greenBA%Bmean(iorb, kpm) + &
                                                boltz*MPI_DOT_PRODUCT(eigen%vec%rc, &
                                                                      OPeigen%vec%rc, split=USE_TRANSPOSE_TRICK_MPI)
                  ENDIF
               ENDDO
            ELSE
               write (log_unit, *) '.... SECTOR IS NOT STABLE UNDER OPERATOR B &
                    &....'
            ENDIF
         endif

      ENDDO

   end subroutine
