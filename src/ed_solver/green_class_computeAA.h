   subroutine scan_indices()

      use globalvar_ed_solver, only: Neigen

      implicit none

      if (allocated(indices_state_sector)) deallocate (indices_state_sector)
      allocate (indices_state_sector(GS%nsector*2*2*2*Neigen*green%N*green%N, 7))
      ktot = 0
      indices_state_sector = 0
      DO isector = 1, GS%nsector
         DO ipm = 1, 2
            DO jpm = 1, 2
               IF (GREEN%compute(ipm, jpm)) then
                  IF (associated(green%correl(ipm, jpm)%MM%MASK%mat)) then
                     IF (ANY(green%correl(ipm, jpm)%MM%MASK%mat)) then
                        DO iph = 1, iph_max
                           DO ieigen = 1, GS%es(isector)%lowest%neigen
                              iorb_f = 0
                              jorb_f = 0
                              do iiorb = 1, green%N
                                 do jjorb = 1, green%N
                                    IF (MASK(ipm, jpm)%mat(iiorb, jjorb)) THEN
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
      ieigen = indices_state_sector(kk, 5)
      iiorb = indices_state_sector(kk, 6)
      jjorb = indices_state_sector(kk, 7)
   end subroutine

   subroutine loop_over_states()

      use globalvar_ed_solver, only: dEmax0, FLAG_FULL_ED_GREEN
      use linalg, only: dexpc

      implicit none

      if (FLAG_FULL_ED_GREEN) then
         if (abs(beta*(es%lowest%eigen(ieigen)%val - E0)) > dEmax0) goto 38
      endif
      eigen => es%lowest%eigen(ieigen)
      boltz = DEXPc(-beta*(eigen%val - E0))/Zpart
      if (messages3) write (*, *) 'APPLY OPERATOR ON STATE', rank
      call apply_creat
      if (messages3) write (*, *) 'CARRY ON LANCZOS', rank
      call lanczos_part
38    continue
   end subroutine

   subroutine clean_everything()

      use common_def, only: timer_fortran
      use correl_class, only: vec2correl
      use eigen_class, only: delete_eigen
      use eigen_sector_class, only: delete_eigensector
      use mask_class, only: delete_mask
      use masked_matrix_class, only: vec2masked_matrix
      use sector_class, only: delete_sector

      implicit none

      integer :: i

      !--------------------!
      ! EXPAND TO MATRICES !
      !--------------------!

      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL vec2masked_matrix(green%correlstat(ipm, jpm))
            IF (compute_dyn_correl .AND. (green%compute(ipm, jpm) .OR. &
                                          green%compute(3 - ipm, 3 - jpm))) CALL &
               vec2correl(green%correl(ipm, jpm))
         ENDDO
      ENDDO

      !----------!
      ! CLEAN UP !
      !----------!

      DO ipm = 1, 2
         DO jpm = 1, 2
            CALL delete_mask(MASK(ipm, jpm))
         ENDDO
      ENDDO
      CALL delete_sector(Asec)
      DO iph = 1, 2
         CALL delete_eigensector(Apm_es(iph))
      ENDDO
      CALL delete_eigen(Apmsym)
      CALL timer_fortran(compute_green_timer, "### COMPUTING "// &
                         TRIM(green%title)//" GREEN S FUNCTION TOOK")

   end subroutine

   subroutine lanczos_part()

      use common_def, only: c2s, find_rank, i2c, &
         reset_timer, timer_fortran
      use eigen_class, only: rank_eigen_in_list
      use genvar, only: pm
      use green_class_compute_symmetric, only: symmetric_combineAA
      use green_class_compute_dynamic, only: compute_dynamic
      use linalg, only: k_to_ij
      use rcvector_class, only: norm_rcvector

      implicit none

      real(8) :: normvec
      integer :: iisector, i_, v(2), iii, start, step
      logical :: do_orb

      !------------------------------------------------------------------!
      ! THEN WE FORM THEIR SYMMETRIC COMBINATIONS & COMPUTE CORRELATIONS !
      !------------------------------------------------------------------!

      if (FLAG_MPI_GREENS == 1) then
         start = rank
         step = size2
      else
         start = 0
         step = 1
      endif

      DO i_ = start + 1, green%N**2, step
         v = k_to_ij(green%N, i_)
         iorb = v(1)
         jorb = v(2)

         if (messages3) write (*, *) 'computing iorb, jorb = ', iorb, jorb

         do_orb = MASK(ipm, jpm)%mat(iorb, jorb)
         if (present(keldysh_level)) then
            do_orb = keldysh_level == iorb .and. keldysh_level == jorb
         endif

         IF (do_orb) THEN
            if (.not. (FLAG_MPI_GREENS == 2) .or. (iorb == iiorb .and. jorb == &
                                                   jjorb)) then
               if (messages3) write (*, *) 'FIND RANK IN MASK', rank, iorb, jorb
               iind = find_rank(MASK(ipm, jpm)%imat(iorb, jorb), MASK(ipm, &
                                                                      jpm)%ivec)
               if (messages3) write (*, *) 'SYM COMBINATION', rank, iind
               CALL symmetric_combineAA(Apmsym, iorb, jorb, ipm, jpm, iph, &
                                        Apm_es, isec_back, iisector, GS, issz)
               if (present(keldysh_level)) then
                  Apmsym%vec%rc = &
                     Apm_es(1)%lowest%eigen(rank_eigen_in_list(jorb, &
                                                               Apm_es(1)%lowest))%vec%rc
               endif
               normvec = 0.d0
               !-------------------------!
               ! EQUAL-TIME CORRELATIONS !
               !-------------------------!
               if (messages3) write (*, *) 'COMPUTE EQUAL TIMES', rank
               IF (iph == 1) THEN
                  normvec = norm_rcvector(Apmsym%vec)
                  green%correlstat(ipm, jpm)%rc%vec(iind) = &
                     green%correlstat(ipm, jpm)%rc%vec(iind) + &
                     boltz*normvec**2
               ENDIF
               !-----------------------!
               ! DYNAMIC  CORRELATIONS !
               !-----------------------!
               if (messages3) write (*, *) 'COMPUTE DYN COR', rank
               if (messages3) write (*, *) 'compute_dyn_correl', &
                  compute_dyn_correl
               if (messages3) write (*, *) 'compute_ipm_jpm', &
                  green%compute(ipm, jpm)

               IF (compute_dyn_correl .AND. green%compute(ipm, jpm)) THEN
                  CALL reset_timer(compute_dyn_timer)
                  CALL compute_dynamic(iph, dyn, green%correl(ipm, jpm)%freq, &
                                       Apmsym, green%correl(ipm, jpm)%stat, green%title, &
                                       normvec, iisector, GS, keldysh_level)
                  if (present(GS_out)) then
                     if (.not. present(keldysh_level)) then
                        write (*, *) 'ERROR, some inputs for keldysh are missing'
                        stop
                     endif
                     if (size(GS_out%es(isector)%lowest%eigen(ieigen)%vec%rc) &
                         /= size(Apmsym%vec%rc)) then
                        write (*, *) 'ERROR Keldysh, problem in copying the &
                             &eigenvalue to new time frame'
                        write (*, *) 'Size of sectors do not match'
                        stop
                     endif
                     if (iorb /= jorb) then
                        write (*, *) 'ERROR KELDYSH, should only be diagional &
                             &orbitals'
                        stop
                     endif
                     GS_out%es(isector)%lowest%eigen(ieigen)%vec%rc = &
                        Apmsym%vec%rc
                  endif
                  green%correl(ipm, jpm)%vec(iind, :) = green%correl(ipm, &
                                                                     jpm)%vec(iind, :) + boltz*dyn
                  CALL timer_fortran(compute_dyn_timer, "# COMPUTING MATRIX &
                       &ELEMENT "//pm(iph)//c2s(i2c(iind))//"/"// &
                       c2s(i2c(MASK(ipm, jpm)%nind))//" TOOK ")
               ENDIF

            ENDIF
         ENDIF

      ENDDO

      return
   end subroutine

   subroutine apply_creat()

      use common_def, only: reset_timer, timer_fortran
      use eigen_class, only: delete_eigenlist, rank_eigen_in_list
      use genvar, only: pm, rank
      use mpi_mod, only: MPI_DOT_PRODUCT
      use sector_class, only: equal_sector

      implicit none

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: entering green_class_computeAA_apply_creat"
#endif

      !------------------------------------------------------!
      ! APPLY CREATION/ANNIHILATION OPERATORS ON THE STATES  !
      !------------------------------------------------------!

      DO kpm = 1, 2

         CALL reset_timer(applyA_timer)
         CALL delete_eigenlist(Apm_es(kpm)%lowest)
         CALL applyA(Apm_es(kpm), pm(kpm), ORBMASKvec(ipm, jpm, :), es, ieigen)
         if (Apm_es(kpm)%lowest%neigen == 0) then
            write (*, '(a,i3,a)') "Error on rank ", rank, " in computeAA: Apm_es neigen = 0"
            stop
         endif
         CALL timer_fortran(applyA_timer, "# APPLYING A"//pm(kpm)//" ON &
              &EIGENSTATE TOOK")

         !--------------------------------------------------------------------------------!
         ! IF SECTOR IS STABLE UNDER OPERATOR, WE ALSO COMPUTE THE MEAN VALUE
         ! OF OPERATOR !
         !--------------------------------------------------------------------------------!

         if (.not. (FLAG_MPI_GREENS == 2) .or. (iiorb == iorb_f .and. jjorb == &
                                                jorb_f)) then
            IF (equal_sector(es%sector, Apm_es(kpm)%sector)) THEN
               write (log_unit, *) '.... SECTOR IS STABLE ....'
               DO iorb = 1, green%N

                  IF (ORBMASKvec(ipm, jpm, iorb)) THEN
                     OPeigen => &
                        Apm_es(kpm)%lowest%eigen(rank_eigen_in_list(iorb, &
                                                                    Apm_es(kpm)%lowest))
                     green%Amean(iorb, kpm) = green%Amean(iorb, kpm) + boltz* &
                                              MPI_DOT_PRODUCT(eigen%vec%rc, OPeigen%vec%rc, split &
                                                              =USE_TRANSPOSE_TRICK_MPI)
                  ENDIF
               ENDDO
            ENDIF
         endif

      ENDDO

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: leaving green_class_computeAA_apply_creat"
#endif

      return

   end subroutine

   subroutine init_iph()

      use common_def, only: dump_message
      use H_class, only: new_H

      implicit none

      write (log_unit, *) ' &
           &%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write (log_unit, *) ' ------ > IPH IPM JPM         : ', iph, ipm, jpm
      write (log_unit, *) ' ------ >  up, dn, tot sites : ', uup, ddn, itot
      write (log_unit, *) ' ------ > neigen              : ', es%lowest%neigen
      write (log_unit, *) ' &
           &%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

      IF (compute_dyn_correl .AND. green%compute(ipm, jpm)) THEN
         SELECT CASE (iph)
         CASE (1)
            CALL dump_message(TEXT="### PARTICLE PART")
            CALL new_H(AIM, Apm_es(jpm)%sector)
         CASE (2)
            CALL dump_message(TEXT="### HOLE PART")
            CALL new_H(AIM, Apm_es(3 - jpm)%sector)
         END SELECT
      ENDIF

      IF (green%compute(ipm, jpm)) THEN
         SELECT CASE (iph)
         CASE (1)
            if (associated(Apm_es(jpm)%sector%sz)) isec_back = &
               Apm_es(jpm)%sector%sz%dimen
         CASE (2)
            if (associated(Apm_es(3 - jpm)%sector%sz)) isec_back = &
               Apm_es(3 - jpm)%sector%sz%dimen
         END SELECT
      ENDIF

      if (es%lowest%neigen == 0) stop 'ERROR in green_computeAA : &
           &es%lowest%neigen = 0'

   end subroutine

   subroutine init_sector()

      use eigen_sector_class, only: new_eigensector, &
         not_commensurate_sector_
      use genvar, only: pm

      implicit none

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: entering green_class_computeAA_init_sector"
#endif

      es => GS%es(isector)

      if (AIM%BATH%SUPER) then
         uup = -1
         ddn = -1
      else
         uup = GS%es(isector)%sector%updo%up%npart
         ddn = GS%es(isector)%sector%updo%down%npart
      endif

      write (log_unit, *) '------ > PARSE EIGEN SECTORS : ', isector, GS%nsector
      write (log_unit, *) ' ...... eigenvalues ...... : ', &
         es%lowest%eigen(:)%val
      write (log_unit, *) ' ...... lowest neig. ...... : ', es%lowest%neigen
      write (log_unit, *) ' ...... eigenvector ...... : ', (/(maxval(abs( &
                                                                     es%lowest%eigen(ii)%vec%rc)), ii=1, es%lowest%neigen)/)

      if (es%lowest%neigen == 0) then
         write (*, *) 'ERROR in green_computeAA : es%lowest%neigen = 0'
         stop
      end if
      if (associated(GS%es(isector)%lowest%eigen)) then
         if (GS%es(isector)%lowest%neigen == 0) then
            write (*, *) 'ERROR in green_computeAA : GS%es%lowest%neigen = 0'
            stop
         end if
      end if

      !------------------------------!
      ! WE CREATE THE TARGET SECTORS !
      !------------------------------!

      NOT_COMMENSURATE = .false.
      DO iph = 1, 2
         CALL Asector(Asec, pm(iph), es%sector)
         CALL new_eigensector(Apm_es(iph), Asec)
         if (associated(es%sector%updo)) then
            write (log_unit, *) &
               '-------------------------------------------------------'
            write (log_unit, *) 'chunk up', iph, es%sector%updo%up%chunk, &
               es%sector%updo%up%title
            write (log_unit, *) 'chunk dn', iph, es%sector%updo%down%chunk, &
               es%sector%updo%down%title
            if (associated(Apm_es(iph)%sector%updo)) then
               write (log_unit, *) 'chunk up A', iph, &
                  Apm_es(iph)%sector%updo%up%chunk, &
                  Apm_es(iph)%sector%updo%up%title
               write (log_unit, *) 'chunk dn A', iph, &
                  Apm_es(iph)%sector%updo%down%chunk, &
                  Apm_es(iph)%sector%updo%down%title
            endif
            write (log_unit, *) &
               '-------------------------------------------------------'
         endif
         if (not_commensurate_sector_(Asec)) NOT_COMMENSURATE = .true.
      ENDDO

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: leaving green_class_computeAA_init_sector"
#endif

   end subroutine

   subroutine init_data()

      use common_def, only: dump_message, find_rank, reset_timer
      use genvar, only: log_unit
      use globalvar_ed_solver, only: force_para_state, para_state
      use eigen_sector_class, only: gsenergy, partition
      use mask_class, only: new_mask

      implicit none

      integer :: ika1, ikb1

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: entering green_class_computeAA_init_data"
#endif

      iph_max = 2
      if (present(keldysh_level)) iph_max = 1

      itot = AIM%bath%Nb + AIM%impurity%Nc

      compute_dyn_correl = .true.
      IF (PRESENT(COMPUTE_DYN)) compute_dyn_correl = COMPUTE_DYN

      CALL reset_timer(compute_green_timer)
      CALL dump_message(TEXT="### START COMPUTING AA "//TRIM(green%title) &
                        //" GREEN S FUNCTION")

      ! INITIALIZATION

      ! BUILD MASKS
      ORBMASKvec = .false.
      DO ipm = 1, 2
         DO jpm = 1, 2
            ! EQUAL-TIME
            CALL new_mask(MASK(ipm, jpm), green%correlstat(ipm, jpm)%rc%MASK)
            ! VECTOR MASK OF ORBITALS: MASK(i) = T IF NEED TO APPLY C(i)
            DO iorb = 1, green%N
               IF (ANY(MASK(ipm, jpm)%mat(iorb, :)) .OR. ANY(MASK(ipm, &
                                                                  jpm)%mat(:, iorb))) ORBMASKvec(ipm, jpm, iorb) = .true.
            ENDDO
         ENDDO
      ENDDO

      ! SET ONLY THOSE MATRIX ELEMENTS THAT WE ARE GOING TO COMPUTE TO ZERO
      DO ipm = 1, 2
         DO jpm = 1, 2
            DO iorb = 1, green%N
               DO jorb = 1, green%N
                  IF (MASK(ipm, jpm)%mat(iorb, jorb)) THEN
                     iind = find_rank(MASK(ipm, jpm)%imat(iorb, jorb), &
                                      MASK(ipm, jpm)%ivec)
                     if (iind == 0) then
                        write (*, *) 'ERROR index is zero! not supposed to &
                             &happen'
                        write (*, *) 'MASK OF THE GF (in green_moduleAA) : '
                        do ika1 = 1, green%N
                           do ikb1 = 1, green%N
                              write (*, *) 'MASK % IMAT (orb1, orb2) : ', ika1, &
                                 ikb1, MASK(ipm, jpm)%imat(ika1, ikb1)
                           enddo
                        enddo
                        write (*, *) 'FOR IPM JPM : ', ipm, jpm
                        write (*, *) 'MASK vec is : ', MASK(ipm, jpm)%ivec
                        write (*, *) 'GF vec is : ', green%correlstat(ipm, &
                                                                      jpm)%rc%MASK%ivec
                        write (*, *) 'FOR IPH_MAX : ', iph_max
                        write (*, *) 'MY RANK     : ', rank
                        write (*, *) 'KELDYSH?    : ', present(keldysh_level)
                        write (*, *) 'force para state / para state flags :', &
                           force_para_state, para_state
                        write (*, *) 'GREEN FUNCTION - CORREL - STATIC PART &
                             &MASK : '
                        do ika1 = 1, green%N
                           do ikb1 = 1, green%N
                              write (*, *) ' GREEN % IMAT (orb1, orb2) ', ika1, &
                                 ikb1, green%correlstat(ipm, &
                                                        jpm)%rc%MASK%imat(ika1, ikb1)
                           enddo
                        enddo

                        stop
                     endif
                     ! EQUAL-TIME
                     green%correlstat(ipm, jpm)%rc%vec(iind) = 0.0_DP
                     ! DYNAMIC
                     IF (compute_dyn_correl .AND. green%compute(ipm, jpm)) &
                        green%correl(ipm, jpm)%vec(iind, :) = 0.0_DP
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      ! MEAN
      DO iorb = 1, green%N
         IF (ANY(ORBMASKvec(:, :, iorb))) THEN
            green%Amean(iorb, :) = 0.0_DP
         ENDIF
      ENDDO

      ! COMPUTE Zpart
      Zpart = partition(beta, GS)
      E0 = GSenergy(GS)

      write (log_unit, *) '-------------------'
      write (log_unit, *) 'COMPUTING GREEN FUNCTION MASK (IPM, JPM) &
           &(hole-particle channels)'
      write (log_unit, *) GREEN%compute(1, :)
      write (log_unit, *) GREEN%compute(2, :)
      write (log_unit, *) '-------------------'

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: leaving green_class_computeAA_init_data"
#endif

   end subroutine
