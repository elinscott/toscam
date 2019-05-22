   subroutine init_correlations(impurity, rbeta, FILEIN, rdelta_width, wmin, &
                                wmax, Niw, Nww)

      use correl_class, only: new_correl
      use common_def, only: close_safe, open_safe, skip_line
      use impurity_class, only: impurity_type, symmetry
      use readable_vec_class, only: new_readable_veclist_from_file
      use genvar, only: bosonic, dbl, fermionic, log_unit
      use globalvar_ed_solver, only: average_g, dens, do_keldysh, &
         donot_compute_holepart, donot_compute_holepart_spm, &
         flag_all_green_func_computed, flag_build_correl_low_part, &
         force_nupdn_basis, force_para_state, force_singlet_state, &
         mask_average, mask_average_sigma, min_all_bath_param, nvecout, &
         para_state, superconducting_state, supersc_state
      use green_class, only: new_green
      use mask_class, only: mask_type, new_mask
      use matrix, only: write_array
      use string5, only: get_unit
      use frequency_class, only: freq_type

      implicit none

      TYPE(impurity_type), INTENT(IN) :: impurity
      CHARACTER(LEN=*), INTENT(IN)  :: FILEIN
      INTEGER              :: UNIT
      TYPE(mask_type)      :: MASK
      INTEGER, ALLOCATABLE :: IMASK_tmp(:, :), IMASK(:, :), IMASKup(:, :), &
                              IMASKdo(:, :), IMASKpair(:, :)
      INTEGER, ALLOCATABLE :: IMASKNS(:, :), IMASKP3(:, :), IMASKP4(:, :)
      REAL(DBL)            :: wmin, wmax, wSzmax, wSpmmax, wNmax, wP3max, &
                              wP4max, rbeta
      INTEGER              :: Nc, ipm, jpm, site, spin, iP3, iP3_, nP3, iP4, &
                              iP4_, nP4, ff, Niw, Nww, jj, i, j, ii
      LOGICAL              :: IS_HERM, to_compute(2, 2), to_compute_in(2, 2)
      TYPE(freq_type)      :: dummyfreq
      real(DBL)            :: rdelta_width
      INTEGER, POINTER     :: mat_temp(:, :)
      INTEGER              :: i_temp

      Nc = impurity%Nc
      width = abs(rdelta_width)
      beta = rbeta

      if (allocated(MASK_AVERAGE)) deallocate (MASK_AVERAGE)
      allocate (MASK_AVERAGE(2*Nc, 2*Nc))
      if (allocated(MASK_AVERAGE_SIGMA)) deallocate (MASK_AVERAGE_SIGMA)
      allocate (MASK_AVERAGE_SIGMA(2*Nc, 2*Nc))

      write (log_unit, *) ' ======== BETA ========== '
      write (log_unit, *) beta
      write (log_unit, *) ' ====== width idelta ==== '
      write (log_unit, *) width
      write (log_unit, *) ' ==== impurity sites ==== '
      write (log_unit, *) Nc
      write (log_unit, *) ' ======================== '

      if (beta < 1.d-15) stop 'error wrong temperature'

      CALL open_safe(UNIT, FILEIN, 'UNKNOWN', 'READ', get_unit=.true.)
      CALL skip_line(UNIT, 3)

#ifdef _complex
      IS_HERM = .false. ! G /= TRANSPOSE(G) IF H COMPLEX
#else
      IS_HERM = .true. ! G  = TRANSPOSE(G) IF H REAL
#endif

      write (log_unit, *) ' READ NAMBU-LIKE IMASKS '

      CALL skip_line(UNIT, 5)
      if (allocated(IMASK)) deallocate (IMASK)
      ALLOCATE (IMASK(Nc*2, Nc*2))

      READ (UNIT, *) (IMASK(site, :), site=1, Nc*2)

      MASK_AVERAGE_SIGMA = IMASK

      if (FLAG_ALL_GREEN_FUNC_COMPUTED) then
         jj = 0
         do ii = 1, 2*Nc
            do j = ii, 2*NC
               jj = jj + 1
               IMASK(ii, j) = jj
            enddo
         enddo
      endif

      if (force_nupdn_basis .or. .not. supersc_state .or. min_all_bath_param > &
          0) then
         IMASK(Nc + 1:2*Nc, 1:Nc) = 0
         IMASK(1:Nc, Nc + 1:2*Nc) = 0
         MASK_AVERAGE_SIGMA(Nc + 1:2*Nc, 1:Nc) = 0
         MASK_AVERAGE_SIGMA(1:Nc, Nc + 1:2*Nc) = 0
      endif

      SUPER = (ANY(IMASK(1:Nc, Nc + 1:Nc*2) /= 0) .OR. ANY(IMASK(Nc + 1:Nc*2, &
                                                                 1:Nc) /= 0))
      superconducting_state = SUPER

      !--------------------------------------------------------!
      if ((force_para_state .or. para_state) .and. (.not. SUPER .or. &
                                                    force_singlet_state)) then
         IMASK(Nc + 1:2*Nc, Nc + 1:2*Nc) = IMASK(1:Nc, 1:Nc)
         MASK_AVERAGE_SIGMA(Nc + 1:2*Nc, Nc + 1:2*Nc) = &
            MASK_AVERAGE_SIGMA(1:Nc, 1:Nc)
      endif
      if (force_singlet_state) then
         do i = 1, Nc - 1
            do j = i + 1, Nc
               IMASK(j, Nc + i) = IMASK(i, Nc + j)
               MASK_AVERAGE_SIGMA(j, Nc + i) = MASK_AVERAGE_SIGMA(i, Nc + j)
            enddo
         enddo
      endif
      !--------------------------------------------------------!

      call impose_imask_sym_singlet_para
      call build_lower_part_(FLAG_BUILD_CORREL_LOW_PART, IMASK)
      call build_lower_part_(FLAG_BUILD_CORREL_LOW_PART, MASK_AVERAGE_SIGMA)

      MASK_AVERAGE = IMASK

      if (abs(average_G) /= 0) then
         if (average_G > 0) then
            MASK_AVERAGE = IMASK
            IMASK = 0
            jj = 0

            do ii = 1, Nc
               do j = ii, NC
                  if (MASK_AVERAGE(ii, j) /= 0) then
                     jj = jj + 1
                     IMASK(ii, j) = jj
                  endif
               enddo
            enddo
            do ii = Nc + 1, 2*Nc
               do j = ii, 2*NC
                  if (MASK_AVERAGE(ii, j) /= 0) then
                     jj = jj + 1
                     IMASK(ii, j) = jj
                  endif
               enddo
            enddo
            do ii = 1, Nc
               do j = Nc + ii, 2*NC
                  if (MASK_AVERAGE(ii, j) /= 0) then
                     jj = jj + 1
                     IMASK(ii, j) = jj
                  endif
               enddo
            enddo
            call impose_imask_sym_singlet_para
            call build_lower_part_(.true., IMASK)
         else
            MASK_AVERAGE = IMASK
            jj = 0
            do ii = 1, 2*Nc
               if (MASK_AVERAGE(ii, ii) /= 0) then
                  if (jj == 0) then
                     jj = maxval(abs(IMASK))
                  else
                     jj = jj + 1
                     IMASK(ii, ii) = jj
                  endif
               endif
            enddo
         endif
      endif

      call write_array(IMASK, ' MASK NAMBU GREEN FUNC CORRELATIONS ', unit= &
                       log_unit)
      call write_array(MASK_AVERAGE, ' MASK AVERAGE ', unit=log_unit)
      call write_array(MASK_AVERAGE_SIGMA, ' MASK AVERAGE SIGMA', unit= &
                       log_unit)

      write (log_unit, *) ' MAKE NAMBU MASK '
      CALL new_mask(MASK, Nc*2, Nc*2, IMASK=IMASK)
      call write_array(MASK%mat, ' MASK NAMBU ', unit=log_unit)
      call write_array(MASK%imat, ' MASK NAMBU ELEMENTS ', unit=log_unit)
      if (Niw <= 0) stop 'init correlations Niw = 0'
      if (Nww <= 0) stop 'init correlations Nww = 0'
      if (beta < 1.d-6) stop 'error SNAMBU beta is zero'
      if (width < 1.d-6) stop 'error SNAMBU width real axis is zero'
      write (log_unit, *) ' COMPUTE THE ANOMALOUS Green Func : ', SUPER
      write (log_unit, *) ' Force nup ndn basis ? : ', force_nupdn_basis

      !-----------------------!
      ! FLAG FOR NORMAL PART  !
      !-----------------------!

      write (log_unit, *) ' BUILDING NAMBU SELF-ENERGY '
      CALL new_correl(SNAMBU, "S", Nc*2, Niw, beta, STAT=FERMIONIC, IMASK= &
                      IMASK)
      CALL new_correl(SNAMBUret, "Sret", Nc*2, Nww, wmin, wmax, width, STAT= &
                      FERMIONIC, IMASK=IMASK)

      write (log_unit, *) ' NAMBU GREEN S FUNCTION '
      to_compute(1, :) = (/.false., .false./)
      to_compute(2, :) = (/.true., .false./)

      CALL new_green(GNAMBU, to_compute, "G", Nc*2, Niw, beta, FERMIONIC, &
                     IMASK=IMASK)
      CALL new_green(GNAMBUret, to_compute, "Gret", Nc*2, Nww, wmin, wmax, &
                     width, FERMIONIC, IMASK=IMASK)
      CALL new_green(GNAMBUN, to_compute, "GN", Nc*2, Niw, beta, FERMIONIC, &
                     IMASK=IMASK)

      write (log_unit, *) ' NORMAL GREEN FCTN '
      write (log_unit, *) ' SPIN UP '

      to_compute(1, :) = (/.false., .true./)
      to_compute(2, :) = (/.true., .false./)

      if (.not. donot_compute_holepart) then
         to_compute_in = to_compute
      else
         to_compute_in = to_compute
         to_compute_in(1, :) = (/.false., .false./)
         write (log_unit, *) 'DO NOT COMPUTE HOLE PART SPIN UP, ASSUME SYMMETRY &
              &(IPM, JPM) (1, 0) and (0, 1)'
         write (log_unit, *) 'to compute is : '
         write (log_unit, *) to_compute_in(1, :)
         write (log_unit, *) to_compute_in(2, :)
      endif

      if (allocated(IMASKup)) deallocate (IMASKup)
      ALLOCATE (IMASKup(Nc, Nc))
      IMASKup = 0
      IMASKup = IMASK(1:Nc, 1:Nc)
      CALL new_green(G(1), to_compute, "G(sz = 1)", Nc, Niw, beta, FERMIONIC, &
                     IMASK=IMASKup, force_compute=to_compute_in)
      CALL new_green(Gret(1), to_compute, "Gret(sz = 1)", Nc, Nww, wmin, wmax, &
                     width, FERMIONIC, IMASK=IMASKup, force_compute=to_compute_in)
      CALL new_green(GN(1), to_compute, "GN(sz = 1)", Nc, Niw, beta, &
                     FERMIONIC, IMASK=IMASKup, force_compute=to_compute_in)

      if (do_keldysh) then
         to_compute(1, :) = (/.false., .true./)
         to_compute(2, :) = (/.false., .false./)
         if (allocated(IMASK_tmp)) deallocate (IMASK_tmp)
         allocate (IMASK_tmp(Nc, Nc))
         IMASK_tmp = 0
         do ii = 1, Nc
            IMASK_tmp(ii, ii) = ii
         enddo
         CALL new_green(GKret(1), to_compute, "GKret(sz = 1)", Nc, Nww, wmin, &
                        wmax, width, FERMIONIC, IMASK=IMASK_tmp, force_compute= &
                        to_compute_in)
      endif

      to_compute_in(1, :) = (/.false., .true./)
      to_compute_in(2, :) = (/.true., .false./)
      CALL new_green(GNF(1), to_compute, "GNF1", Nc, Niw, beta, FERMIONIC, &
                     IMASK=IMASKup, AB=.true., force_compute=to_compute_in)
      CALL new_green(GNF(2), to_compute, "GNF2", Nc, Niw, beta, FERMIONIC, &
                     IMASK=IMASKup, AB=.true., force_compute=to_compute_in)

      write (log_unit, *) ' SPIN DOWN '
      to_compute(1, :) = (/.false., .true./)
      to_compute(2, :) = (/.true., .false./)

      if (.not. donot_compute_holepart) then
         to_compute_in = to_compute
      else
         to_compute_in(1, :) = (/.false., .false./)
         to_compute_in = to_compute
      endif
      if (para_state .and. .not. SUPER) then
         write (*, *) 'WARNING : paramagnetic state or supra - flags para and &
              &super are : ', para_state, SUPER
         write (*, *) '          turning off spin down GF '
         to_compute_in = .false.
      endif

      if (allocated(IMASKdo)) deallocate (IMASKdo)
      ALLOCATE (IMASKdo(Nc, Nc))
      IMASKdo = 0
      IMASKdo = IMASK(Nc + 1:Nc*2, Nc + 1:Nc*2)

      CALL new_green(G(2), to_compute, "G(sz = -1)", Nc, Niw, beta, FERMIONIC, &
                     IMASK=IMASKdo, force_compute=to_compute_in)
      CALL new_green(Gret(2), to_compute, "Gret(sz = -1)", Nc, Nww, wmin, &
                     wmax, width, FERMIONIC, IMASK=IMASKdo, force_compute= &
                     to_compute_in)
      CALL new_green(GN(2), to_compute, "GN(sz = -1)", Nc, Niw, beta, &
                     FERMIONIC, IMASK=IMASKdo, force_compute=to_compute_in)

      if (do_keldysh) then
         to_compute(1, :) = (/.false., .true./)
         to_compute(2, :) = (/.false., .false./)
         if (allocated(IMASK_tmp)) deallocate (IMASK_tmp)
         allocate (IMASK_tmp(Nc, Nc))
         IMASK_tmp = 0
         do ii = 1, Nc
            IMASK_tmp(ii, ii) = ii
         enddo
         CALL new_green(GKret(2), to_compute, "GKret(sz = -1)", Nc, Nww, wmin, &
                        wmax, width, FERMIONIC, IMASK=IMASK_tmp, force_compute= &
                        to_compute_in)
      endif

      write (log_unit, *) ' cleans redundancies '
      DO ipm = 1, 2
         DO jpm = 1, 2
            IF (G(2)%compute(ipm, jpm)) THEN
               G(2)%correl(ipm, jpm)%MM%MASK%mat = MASK%mat(Nc + 1:Nc*2, Nc + &
                                                            1:Nc*2)
               Gret(2)%correl(ipm, jpm)%MM%MASK%mat = MASK%mat(Nc + 1:Nc*2, Nc &
                                                               + 1:Nc*2)
            ENDIF
            G(2)%correlstat(ipm, jpm)%rc%MASK%mat = MASK%mat(Nc + 1:Nc*2, Nc &
                                                             + 1:Nc*2)
            Gret(2)%correlstat(ipm, jpm)%rc%MASK%mat = MASK%mat(Nc + 1:Nc*2, &
                                                                Nc + 1:Nc*2)
         ENDDO
      ENDDO

      write (log_unit, *) 'MASK IMAT G(1) IPM = JPM = 1'
      ipm = 1
      jpm = 1

      IF (G(1)%compute(ipm, jpm)) THEN
         do i = 1, Nc
            do j = 1, Nc
               write (log_unit, *) i, j, G(1)%correl(ipm, jpm)%MM%MASK%imat(i, j)
            enddo
         enddo
         write (log_unit, *)
      endif
      write (log_unit, *) 'MASK IMAT G(1) IPM = JPM = 2'
      ipm = 2
      jpm = 2

      IF (G(1)%compute(ipm, jpm)) THEN
         do i = 1, Nc
            do j = 1, Nc
               write (log_unit, *) i, j, G(1)%correl(ipm, jpm)%MM%MASK%imat(i, j)
            enddo
         enddo
         write (log_unit, *)
      endif
      write (log_unit, *) 'MASK IMAT G(2) IPM = JPM = 1'
      ipm = 1
      jpm = 1

      IF (G(2)%compute(ipm, jpm)) THEN
         do i = 1, Nc
            do j = 1, Nc
               write (log_unit, *) i, j, G(2)%correl(ipm, jpm)%MM%MASK%imat(i, j)
            enddo
         enddo
         write (log_unit, *)
      endif
      write (log_unit, *) 'MASK IMAT G(2) IPM = JPM = 2'
      ipm = 2
      jpm = 2

      IF (G(2)%compute(ipm, jpm)) THEN
         do i = 1, Nc
            do j = 1, Nc
               write (log_unit, *) i, j, G(2)%correl(ipm, jpm)%MM%MASK%imat(i, j)
            enddo
         enddo
         write (log_unit, *)
      endif

      do ipm = 1, 2
         do jpm = 1, 2
            IF (G(1)%compute(ipm, jpm)) THEN
               mat_temp => G(1)%correl(ipm, jpm)%MM%MASK%imat
               i_temp = maxval(G(1)%correl(ipm, jpm)%MM%MASK%ivec)
               if (any(mat_temp > i_temp)) then
                  write (*, *) 'ERROR : element outside range, correl init, &
                       &spin up'
                  write (*, *) 'vec = ', G(1)%correl(ipm, jpm)%MM%MASK%ivec
                  write (*, *) 'force para state / para state flags :', &
                     force_para_state, para_state
                  stop
               endif
            endif
         enddo
      enddo

      do ipm = 1, 2
         do jpm = 1, 2
            IF (G(2)%compute(ipm, jpm)) THEN
               mat_temp => G(2)%correl(ipm, jpm)%MM%MASK%imat
               i_temp = maxval(G(2)%correl(ipm, jpm)%MM%MASK%ivec)
               if (any(mat_temp > i_temp)) then
                  write (*, *) 'ERROR : element outside range, correl init, &
                       &spin down'
                  write (*, *) 'vec = ', G(2)%correl(ipm, jpm)%MM%MASK%ivec
                  write (*, *) 'force para state / para state flags :', &
                     force_para_state, para_state
                  stop
               endif
            endif
         enddo
      enddo

      !--------------------------!
      ! FLAG FOR ANOMALOUS PART  !
      !--------------------------!

      IF (SUPER) THEN

         write (log_unit, *) ' ==== will compute anomalous Green functions &
              &GF(1) ==== '

         to_compute(1, :) = (/.true., .false./)
         to_compute(2, :) = (/.false., .true./)

#ifdef _complex
         ! < ci[down](z) * cj[up] > IS INDEPENDANT !
         to_compute_in(1, :) = [.true., .false.]
         to_compute_in(2, :) = [.false., .true.]
#else
         ! COMPUTE ONLY : < ci[up](z) * cj[down] >
         if (donot_compute_holepart) then
            to_compute_in(1, :) = [.true., .false.]
            to_compute_in(2, :) = [.false., .false.]
         else
            to_compute_in(1, :) = [.true., .false.]
            to_compute_in(2, :) = [.false., .true.]
         endif
#endif

         if (allocated(IMASKpair)) deallocate (IMASKpair)
         ALLOCATE (IMASKpair(Nc, Nc))
         IMASKpair = 0
         IMASKpair = IMASK(1:Nc, Nc + 1:Nc*2)
         CALL new_green(GF(1), to_compute, "Fupdo", Nc, Niw, beta, FERMIONIC, &
                        IMASK=IMASKpair, AB=.true., force_compute=to_compute_in)
         CALL new_green(GFret(1), to_compute, "Fupdoret", Nc, Nww, wmin, wmax, &
                        width, FERMIONIC, IMASK=IMASKpair, AB=.true., force_compute= &
                        to_compute_in)

#ifdef _complex
         write (log_unit, *) ' ==== will compute anomalous Green functions &
              &GF(2) (complex parameters) ==== '
         CALL new_green(GF(2), to_compute, "Fdoup", Nc, Niw, beta, FERMIONIC, &
                        IMASK=IMASKpair, AB=T, force_compute=to_compute_in)
         CALL new_green(GFret(2), to_compute, "Fdoupret", Nc, Nww, wmin, wmax, &
                        width, FERMIONIC, IMASK=IMASKpair, AB=.true., force_compute= &
                        to_compute_in)
         if (.not. GF(1)%compute(1, 1) .or. .not. GF(1)%compute(2, 2)) stop 'RVB &
              &in ED SOLVER but not computing anomalous part GF1, some error?'
         if (.not. GF(2)%compute(1, 1) .or. .not. GF(2)%compute(2, 2)) stop 'RVB &
              &in ED SOLVER but not computing anomalous part GF2, some error?'
#endif

      ENDIF

      !------------------------!
      ! BOSONIC GREEN S FCTNS  !
      !------------------------!

      ! N-SS CORRELATIONS

      ! FLAGS
      CALL skip_line(UNIT, 3)
      READ (UNIT, *) compute_Sz
      READ (UNIT, *) compute_Spm
      READ (UNIT, *) compute_N

      ! MASK
      if (allocated(IMASKNS)) deallocate (IMASKNS)
      ALLOCATE (IMASKNS(Nc, Nc))
      IMASKNS = 0
      READ (UNIT, *) (IMASKNS(site, :), site=1, Nc)

      ! REAL (RETARDED) FREQ.
      READ (UNIT, *) wSzmax
      READ (UNIT, *) wSpmmax
      READ (UNIT, *) wNmax

      write (log_unit, *) '----------- bosonic green functions ---------------'
      write (log_unit, *) '  wSzmax, sSpmmax, wNmax : ', wSzmax, wSpmmax, wNmax
      write (log_unit, *) '  beta                 : ', beta
      write (log_unit, *) '  width                : ', width
      write (log_unit, *) '---------------------------------------------------'

      if (Niw <= 0) stop 'init correlations Niw = 0'
      if (Nww <= 0) stop 'init correlations Nww = 0'

      to_compute(1, :) = (/.false., .false./)  !    COMPU.true.E   ONLY
      to_compute(2, :) = (/.true., .false./)  ! < Sz[i](i*w) * Sz[j] >

      CALL new_green(Sz, to_compute, "Sz", Nc, Niw, beta, BOSONIC, IMASK= &
                     IMASKNS)
      CALL new_green(Szret, to_compute, "Szret", Nc, Nww, -wSzmax, wSzmax, &
                     width, BOSONIC, IMASK=IMASKNS)

      to_compute(1, :) = (/.false., .true./)  ! COMPUTE    < S-[i](i*w) * S + [j] >
      to_compute(2, :) = (/.true., .false./)  ! AND DEDUCE < S + [i](i*w) * S-[j] >

      if (.not. donot_compute_holepart_spm) then
         to_compute_in = to_compute
      else
         to_compute_in = to_compute
         to_compute_in(1, :) = (/.false., .false./)
      endif
      CALL new_green(Spm, to_compute, "Spm", Nc, Niw, beta, BOSONIC, IMASK= &
                     IMASKNS, force_compute=to_compute_in)
      CALL new_green(Spmret, to_compute, "Spmret", Nc, Nww, -wSpmmax, wSpmmax, &
                     width, BOSONIC, IMASK=IMASKNS, force_compute=to_compute_in)

      to_compute(1, :) = (/.false., .false./)  !    COMPUTE   ONLY
      to_compute(2, :) = (/.true., .false./)  ! < N[i](i*w) * N[j] >

      CALL new_green(N, to_compute, "N", Nc, Niw, beta, BOSONIC, IMASK= &
                     IMASKNS)
      CALL new_green(Nret, to_compute, "Nret", Nc, Nww, -wNmax, wNmax, width, &
                     BOSONIC, IMASK=IMASKNS)

      ! P3 AND CHI CORRELATIONS

      write (log_unit, *) 'P3 and CHI correlations, FILEIN, unit : ', FILEIN, &
         UNIT
      CALL skip_line(UNIT, 3)
      READ (UNIT, *) compute_P3
      write (log_unit, *) 'compute_P3 : ', compute_P3

      ! LIST OF TRIPLETS
      READ (UNIT, *) nP3
      write (log_unit, *) 'how many triplets : ', nP3

      IF (nP3 > 0) then
         if (allocated(triplets)) deallocate (triplets)
         ALLOCATE (triplets(nP3, 3))
         CALL skip_line(UNIT, 1)
         READ (UNIT, *) (iP3, triplets(iP3, 1:3), iP3_=1, nP3)
         CALL skip_line(UNIT, 1)
         if (allocated(IMASKP3)) deallocate (IMASKP3)
         ALLOCATE (IMASKP3(nP3, nP3))
         IMASKP3 = 0
         READ (UNIT, *) (IMASKP3(iP3, :), iP3=1, nP3)
      else
         call skip_line(UNIT, 2)
      ENDIF

      ! FREQ.
      READ (UNIT, *) wP3max
      write (log_unit, *) '... max frequency ...', wP3max

      ! P3 CORRELATIONS
      to_compute(1, :) = (/.true., .true./) ! COMPUTE < P3-[i](i*w) * P3-[j] > AND <
      ! P3-[i](i*w) * P3 + [j] >
      to_compute(2, :) = (/.true., .true./) ! AND DEDUCE < P3 + [i](i*w) * P3 + [j] >
      ! AND < P3 + [i](i*w) * P3-[j] >
      if (nP3 > 0) then
         CALL new_green(P3, to_compute, "P3", Nc, Niw, beta, BOSONIC, IMASK= &
                        IMASKP3)
         CALL new_green(P3ret, to_compute, "P3ret", Nc, Nww, -wP3max, wP3max, &
                        width, BOSONIC, IMASK=IMASKP3)
      endif

      ! .... THEN DEDUCE CHI CORRELATIONS ....
      if (nP3 > 0) then
         CALL new_correl(CHI, "CHI", Nc, Niw, beta, BOSONIC, IMASK=IMASKP3)
         CALL new_correl(CHIret, "CHIret", Nc, Nww, -wP3max, wP3max, width, &
                         BOSONIC, IMASK=IMASKP3)
      endif

      ! P4 CORRELATIONS
      CALL skip_line(UNIT, 3)
      READ (UNIT, *) compute_P4
      write (log_unit, *) "compute P4 : ", compute_P4

      ! LIST OF QUADRUPLETS
      READ (UNIT, *) nP4
      write (log_unit, *) 'how many quadruplets : ', nP4

      IF (nP4 > 0) then
         if (allocated(quadruplets)) deallocate (quadruplets)
         ALLOCATE (quadruplets(nP4, 4))
         CALL skip_line(UNIT, 1)
         READ (UNIT, *) (iP4, quadruplets(iP4, 1:4), iP4_=1, nP4)
         if (allocated(IMASKP4)) deallocate (IMASKP4)
         ALLOCATE (IMASKP4(nP4, nP4))
         IMASKP4 = 0
         CALL skip_line(UNIT, 1)
         READ (UNIT, *) (IMASKP4(iP4, :), iP4=1, nP4)
      else
         call skip_line(UNIT, 2)
      endif

      ! ... FREQ ... !

      READ (UNIT, *) wP4max
      write (log_unit, *) '... max frequency ...', wP4max

      to_compute(1, :) = (/.true., .true./) ! COMPUTE < P4 + [i](i*w) * P4 + [j] >, < P4
      ! + [i](i*w) * P4-[j] >
      to_compute(2, :) = (/.true., .true./) ! AND < P4-[i](i*w) * P4 + [j] >, <
      ! P4-[i](i*w) * P4-[j] >
      if (nP4 > 0) then
         CALL new_green(P4, to_compute, "P4", Nc, Niw, beta, BOSONIC, IMASK= &
                        IMASKP4)
         CALL new_green(P4ret, to_compute, "P4ret", Nc, Nww, -wP4max, wP4max, &
                        width, BOSONIC, IMASK=IMASKP4)
      endif

      ! INITIALIZE EQUAL-TIME QUANTITIES
      if (allocated(dens)) deallocate (dens)
      ALLOCATE (dens(Nc))

      ! READ IMPURITY STATES WHOSE WEIGHT WE WANT TO COMPUTE COPY NECESSARY
      ! PARAMETERS TO COMPUTE WEIGHTS
      CALL new_readable_veclist_from_file(vec_list, UNIT, impurity%iorb, NAMBU &
                                          =SUPER, nvecout=nvecout)
      CALL close_safe(UNIT)

   contains

      subroutine impose_imask_sym_singlet_para()

         use globalvar_ed_solver, only: force_para_state, force_singlet_state

         implicit none

         if (force_para_state .and. (.not. SUPER .or. force_singlet_state)) then
            IMASK(Nc + 1:2*Nc, Nc + 1:2*Nc) = IMASK(1:Nc, 1:Nc)
         endif
         if (force_singlet_state) then
            do i = 1, Nc - 1
               do j = i + 1, Nc
                  IMASK(j, Nc + i) = IMASK(i, Nc + j)
                  IMASK(Nc + j, i) = IMASK(Nc + i, j)
               enddo
            enddo
         endif
         !--------------------------------------------------------!

      end subroutine

      subroutine build_lower_part_(do_it_, MASK)

         use genvar, only: logfile

         implicit none

         logical :: do_it_
         integer :: MASK(:, :)

#ifdef _complex
         if (do_it_) then
            jj = maxval(abs([(MASK(ii, ii:2*Nc), ii=1, 2*Nc)]))
            write (logfile, *) 'max in upper triangle in mask : ', jj
            do ii = 2, 2*Nc
               do j = 1, ii - 1
                  if (MASK(j, ii) /= 0) then
                     jj = jj + 1
                     MASK(ii, j) = jj
                  endif
               enddo
            enddo
         endif
#else
         if (do_it_) then
            do ii = 2, 2*Nc
               do j = 1, ii - 1
                  MASK(ii, j) = MASK(j, ii)
               enddo
            enddo
         endif
#endif

      end subroutine

   end subroutine

