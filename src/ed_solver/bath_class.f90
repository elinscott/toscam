MODULE bath_class

   use correl_class, only: correl_type
   use genvar, only: DP
   use masked_matrix_class, only: masked_matrix_type

   implicit none

   TYPE bath_type
      !$$$$$$$$$$$$$$$
      !$$BATH TYPE$$
      !$$$$$$$$$$$$$$$
      ! FILEOUT
      CHARACTER(LEN=100) :: fileout = '\0'
      ! NUMBER OF BATH SITES
      INTEGER :: Nb = 0
      ! NUMBER OF 1-PARTICLE ORBITALS
      INTEGER :: norbs = 0
      ! SIZE OF REDUCED HILBERT SPACE OF THE BATH
      INTEGER :: nstates = 0
      ! RANK OF ORBITALS
      INTEGER, POINTER :: iorb(:, :) => NULL()
      ! NUMBER OF IMPURITY SITES (FOR HYBRIDIZATION)
      INTEGER :: Nc = 0
      LOGICAL :: SUPER = .false. ! TRUE IF SUPERCONDUCTING
      ! CONDUCTION BATH ENERGY
      TYPE(masked_matrix_type), POINTER :: Eb(:) => NULL() ! in (site, site) basis for a given spin
      ! HYBRIDIZATION
      TYPE(masked_matrix_type), POINTER :: Vbc(:) => NULL()! in (site[impurity], site[bath]) basis
      TYPE(masked_matrix_type), POINTER :: PVbc(:) => NULL()! in (site[impurity], site[bath]) basis
      ! PAIRING ENERGY
      TYPE(masked_matrix_type)          :: Pb       ! in (site, site) basis
      ! TOTAL NUMBER OF (REAL) PARAMETERS
      INTEGER   :: nparam = 0
      REAL(DP), POINTER :: vec(:) => NULL()
      ! HYBRID => BATH PARAMETERS
      INTEGER   :: Niter_search_max = 0 ! max # of iterations in search routine (conj. grad.)
      REAL(DP) :: search_step = 0.0_DP ! small step in initial search direction
      REAL(DP) :: dist_max = 0.0_DP ! max. error on hybridization functions
      ! HYBRIDIZATION FUNCTION
      TYPE(correl_type) :: hybrid, hybridret
   END TYPE

   INTERFACE new_bath
      MODULE PROCEDURE new_bath_from_scratch
      MODULE PROCEDURE new_bath_from_old
   END INTERFACE

   public :: bath_type
   public :: copy_bath
   public :: nambu_eb
   public :: nambu_vbc
   public :: nbathparam_
   public :: new_bath
   public :: read_bath

contains

   subroutine which_basis_to_use(BATH)

      implicit none

      TYPE(bath_type), INTENT(INOUT) :: BATH

      BATH%SUPER = which_basis_to_use__(BATH%Pb%rc%MASK%mat)
   end subroutine

   logical function which_basis_to_use__(mat)

      use genvar, only: rank
      use globalvar_ed_solver, only: force_nupdn_basis, force_sz_basis, para_state, supersc_state

      implicit none

      logical :: mat(:, :)

      which_basis_to_use__ = .false.
      if (.not. which_basis_to_use__ .and. para_state) which_basis_to_use__ = &
         .false.
      if (.not. which_basis_to_use__ .and. ANY(mat)) which_basis_to_use__ = &
         .true.
      if (.not. which_basis_to_use__ .and. supersc_state) which_basis_to_use__ &
         = .true.
      if (.not. which_basis_to_use__ .and. force_sz_basis) &
         which_basis_to_use__ = .true.
      if (which_basis_to_use__ .and. force_nupdn_basis) which_basis_to_use__ &
         = .false.

      write (145 + rank, *) ' ---- > USING SUPERCONDUCTIVITY ? : ', &
         which_basis_to_use__
      write (145 + rank, *) 'any          : ', ANY(mat)
      write (145 + rank, *) 'para         : ', para_state
      write (145 + rank, *) 'supersc      : ', supersc_state
      write (145 + rank, *) 'force_sz     : ', force_sz_basis
      write (145 + rank, *) 'force_nupndn : ', force_nupdn_basis

   end function

   subroutine new_bath_from_scratch(BATH, Nb, Nc, IMASKE, IMASKP, IMASKV, &
                                    IMASKPV)

      ! CREATE BATH FROM SCRATCH

      use genvar, only: cspin
      use globalvar_ed_solver, only: istati
      use masked_matrix_class, only: clean_redundant_imask, new_masked_matrix
      use linalg, only: ramp

      implicit none

      TYPE(bath_type), INTENT(INOUT) :: BATH
      INTEGER, INTENT(IN)            :: Nb
      INTEGER, INTENT(IN)            :: Nc
      INTEGER, OPTIONAL, INTENT(IN)  :: IMASKE(Nb, Nb, 2), IMASKP(Nb, Nb), &
                                        IMASKV(Nb, Nc, 2), IMASKPV(Nb, Nc, 2)
      INTEGER :: spin

      CALL delete_bath(BATH)
      BATH%Nb = Nb
      BATH%norbs = Nb*2
      BATH%nstates = 2**BATH%norbs

      ! ORDERING OF ORBITALS WITH INCREASING RANK |(site, up) > |(site, down) >

      if (ASSOCIATED(BATH%iorb)) DEALLOCATE (BATH%iorb, STAT=istati)
      ALLOCATE (BATH%iorb(Nb, 2))

      CALL ramp(BATH%iorb(:, 1), Nb)
      BATH%iorb(:, 2) = BATH%iorb(:, 1) + Nb
      BATH%Nc = Nc

      BATH%fileout = ''
      BATH%nparam = 0

      ! CONDUCTION BATH ENERGY
      IF (PRESENT(IMASKE)) THEN
         IF (ASSOCIATED(BATH%Eb)) DEALLOCATE (BATH%Eb, STAT=istati)
         ALLOCATE (BATH%Eb(SIZE(IMASKE, 3)))
         DO spin = 1, SIZE(IMASKE, 3)
            CALL new_masked_matrix(BATH%Eb(spin), "Eb(sz = "// &
                                   TRIM(cspin(spin))//")", Nb, Nb, IMASK=IMASKE(:, :, spin), &
                                   IS_HERM=.true.)
         ENDDO
         CALL clean_redundant_imask(BATH%Eb)
      ENDIF

      ! HYBRIDIZATION
      IF (PRESENT(IMASKV)) THEN
         IF (ASSOCIATED(BATH%Vbc)) DEALLOCATE (BATH%Vbc, STAT=istati)
         ALLOCATE (BATH%Vbc(SIZE(IMASKV, 3)))
         DO spin = 1, SIZE(IMASKV, 3)
            CALL new_masked_matrix(BATH%Vbc(spin), "Vbc(sz = "// &
                                   TRIM(cspin(spin))//")", Nb, Nc, IMASK=IMASKV(:, :, spin))
         ENDDO
         CALL clean_redundant_imask(BATH%Vbc)
      ENDIF

      IF (PRESENT(IMASKPV)) THEN
         IF (ASSOCIATED(BATH%PVbc)) DEALLOCATE (BATH%PVbc, STAT=istati)
         ALLOCATE (BATH%PVbc(SIZE(IMASKPV, 3)))
         DO spin = 1, SIZE(IMASKPV, 3)
            CALL new_masked_matrix(BATH%PVbc(spin), "PVbc(sz = "// &
                                   TRIM(cspin(spin))//")", Nb, Nc, IMASK=IMASKPV(:, :, &
                                                                                 spin))
         ENDDO
         CALL clean_redundant_imask(BATH%PVbc)
      ENDIF

      ! PAIRING ENERGY
      IF (PRESENT(IMASKP)) THEN
         CALL new_masked_matrix(BATH%Pb, "Pb", Nb, Nb, IMASK=IMASKP)
         call which_basis_to_use(BATH)
      ENDIF

   end subroutine

   subroutine new_bath_from_old(BATHOUT, BATHIN)

      use globalvar_ed_solver, only: istati

      implicit none

      TYPE(bath_type), INTENT(INOUT) :: BATHOUT
      TYPE(bath_type), INTENT(IN)    :: BATHIN
      INTEGER, ALLOCATABLE :: IMASKE(:, :, :), IMASKV(:, :, :), IMASKP(:, :), &
                              IMASKPV(:, :, :)
      INTEGER              :: spin

      IF (.NOT. ASSOCIATED(BATHIN%Eb)) STOP "ERROR IN new_bath_from_old: INPUT &
           &ISNT ALLOCATED!"

      IF (ALLOCATED(IMASKE)) DEALLOCATE (IMASKE, STAT=istati)
      ALLOCATE (IMASKE(BATHIN%Nb, BATHIN%Nb, 2))
      IF (ALLOCATED(IMASKV)) DEALLOCATE (IMASKV, STAT=istati)
      ALLOCATE (IMASKV(BATHIN%Nb, BATHIN%Nc, 2))
      IF (ALLOCATED(IMASKPV)) DEALLOCATE (IMASKPV, STAT=istati)
      ALLOCATE (IMASKPV(BATHIN%Nb, BATHIN%Nc, 2))
      IF (ALLOCATED(IMASKP)) DEALLOCATE (IMASKP, STAT=istati)
      ALLOCATE (IMASKP(BATHIN%Nb, BATHIN%Nb))

      DO spin = 1, 2
         IMASKE(:, :, spin) = BATHIN%Eb(spin)%rc%MASK%imat
         IMASKV(:, :, spin) = BATHIN%Vbc(spin)%rc%MASK%imat
         IMASKPV(:, :, spin) = BATHIN%PVbc(spin)%rc%MASK%imat
      ENDDO

      if (ASSOCIATED(BATHIN%Pb%rc%MASK%imat)) then
         IMASKP = BATHIN%Pb%rc%MASK%imat
      endif

      CALL new_bath_from_scratch(BATHOUT, BATHIN%Nb, BATHIN%Nc, IMASKE= &
                                 IMASKE, IMASKP=IMASKP, IMASKV=IMASKV, IMASKPV=IMASKPV)
      CALL copy_bath(BATHOUT, BATHIN)

      IF (ALLOCATED(IMASKE)) DEALLOCATE (IMASKE, STAT=istati)
      IF (ALLOCATED(IMASKV)) DEALLOCATE (IMASKV, STAT=istati)
      IF (ALLOCATED(IMASKPV)) DEALLOCATE (IMASKPV, STAT=istati)
      IF (ALLOCATED(IMASKP)) DEALLOCATE (IMASKP, STAT=istati)

   end subroutine

   subroutine copy_bath(BATHOUT, BATHIN)

      ! CREATE BATHOUT BY COPYING EXISTING BATHIN

      use correl_class, only: copy_correl, new_correl
      use masked_matrix_class, only: copy_masked_matrix

      implicit none

      TYPE(bath_type), INTENT(INOUT) :: BATHOUT
      TYPE(bath_type), INTENT(IN)    :: BATHIN
      INTEGER :: spin

      BATHOUT%Nb = BATHIN%Nb
      BATHOUT%norbs = BATHIN%norbs
      BATHOUT%nstates = BATHIN%nstates
      BATHOUT%iorb = BATHIN%iorb
      BATHOUT%Nc = BATHIN%Nc
      BATHOUT%fileout = BATHIN%fileout
      BATHOUT%nparam = BATHIN%nparam

      IF (ASSOCIATED(BATHIN%vec)) THEN
         IF (.NOT. ASSOCIATED(BATHOUT%vec)) ALLOCATE (BATHOUT%vec(BATHOUT%nparam))
         BATHOUT%vec = BATHIN%vec
      ENDIF

      ! CONDUCTION BATH ENERGY
      DO spin = 1, SIZE(BATHIN%Eb)
         CALL copy_masked_matrix(BATHOUT%Eb(spin), BATHIN%Eb(spin))
      ENDDO

      ! HYBRIDIZATION
      DO spin = 1, SIZE(BATHIN%Vbc)
         CALL copy_masked_matrix(BATHOUT%Vbc(spin), BATHIN%Vbc(spin))
      ENDDO

      DO spin = 1, SIZE(BATHIN%PVbc)
         CALL copy_masked_matrix(BATHOUT%PVbc(spin), BATHIN%PVbc(spin))
      ENDDO

      BATHOUT%SUPER = BATHIN%SUPER

      ! PAIRING ENERGY
      IF (BATHOUT%SUPER) CALL copy_masked_matrix(BATHOUT%Pb, BATHIN%Pb)
      IF (ASSOCIATED(BATHIN%hybrid%fctn)) THEN
         IF (.NOT. ASSOCIATED(BATHOUT%hybrid%fctn)) THEN
            CALL new_correl(BATHOUT%hybrid, BATHIN%hybrid)
         ELSE
            CALL copy_correl(BATHOUT%hybrid, BATHIN%hybrid)
         ENDIF
      ENDIF
      IF (ASSOCIATED(BATHIN%hybridret%fctn)) THEN
         IF (.NOT. ASSOCIATED(BATHOUT%hybridret%fctn)) THEN
            CALL new_correl(BATHOUT%hybridret, BATHIN%hybridret)
         ELSE
            CALL copy_correl(BATHOUT%hybridret, BATHIN%hybridret)
         ENDIF
      ENDIF
      BATHOUT%Niter_search_max = BATHIN%Niter_search_max
      BATHOUT%search_step = BATHIN%search_step
      BATHOUT%dist_max = BATHIN%dist_max

   end subroutine

   subroutine delete_bath(BATH)

      use correl_class, only: delete_correl
      use globalvar_ed_solver, only: istati
      use masked_matrix_class, only: delete_masked_matrix

      implicit none

      TYPE(bath_type), INTENT(INOUT) :: BATH
      INTEGER :: spin

      IF (ASSOCIATED(BATH%iorb)) DEALLOCATE (BATH%iorb, STAT=istati)
      IF (ASSOCIATED(BATH%vec)) DEALLOCATE (BATH%vec, STAT=istati)

      IF (ASSOCIATED(BATH%Eb)) THEN
         DO spin = 1, SIZE(BATH%Eb)
            CALL delete_masked_matrix(bath%Eb(spin))
         ENDDO
         DEALLOCATE (BATH%Eb, STAT=istati)
      ENDIF

      IF (ASSOCIATED(BATH%Vbc)) THEN
         DO spin = 1, SIZE(BATH%Vbc)
            CALL delete_masked_matrix(bath%Vbc(spin))
         ENDDO
         DEALLOCATE (BATH%Vbc, STAT=istati)
      ENDIF

      IF (ASSOCIATED(BATH%PVbc)) THEN
         DO spin = 1, SIZE(BATH%PVbc)
            CALL delete_masked_matrix(bath%PVbc(spin))
         ENDDO
         DEALLOCATE (BATH%PVbc, STAT=istati)
      ENDIF

      IF (BATH%SUPER) CALL delete_masked_matrix(BATH%Pb)
      CALL delete_correl(bath%hybrid)
      CALL delete_correl(bath%hybridret)

   end subroutine

   subroutine read_raw_bath(BATH, FILEIN)

      ! $$ READ RAW BATH PARAMETERS

      use common_def, only: close_safe, open_safe
      use genvar, only: log_unit
      use globalvar_ed_solver, only: min_all_bath_param
      use masked_matrix_class, only: clean_redundant_imask, &
         read_raw_masked_matrix
      use string5, only: get_unit

      implicit none

      TYPE(bath_type), INTENT(INOUT) :: BATH
      CHARACTER(LEN=*), INTENT(IN) :: FILEIN
      INTEGER :: UNIT
      INTEGER :: Nb, Nc, spin, nspin
      LOGICAL :: read_hybrid

      CALL open_safe(UNIT, FILEIN, 'UNKNOWN', 'READ', get_unit=.true.)

      ! READ BATH SCALARS/ALLOCATE BATH
      READ (UNIT, *) Nb
      READ (UNIT, *) Nc

      CALL new_bath_from_scratch(BATH, Nb, Nc)

      ! READ BATH
      READ (UNIT, *) BATH%fileout
      READ (UNIT, *) BATH%nparam
      READ (UNIT, *) nspin

      ALLOCATE (BATH%Eb(nspin))
      DO spin = 1, nspin
         CALL read_raw_masked_matrix(BATH%Eb(spin), UNIT)
      ENDDO
      CALL clean_redundant_imask(BATH%Eb)
      READ (UNIT, *) nspin

      ALLOCATE (BATH%Vbc(nspin))
      DO spin = 1, nspin
         CALL read_raw_masked_matrix(BATH%Vbc(spin), UNIT)
      ENDDO
      CALL clean_redundant_imask(BATH%Vbc)

      ALLOCATE (BATH%PVbc(nspin))
      DO spin = 1, nspin
         CALL read_raw_masked_matrix(BATH%PVbc(spin), UNIT)
      ENDDO
      CALL clean_redundant_imask(BATH%PVbc)

      CALL read_raw_masked_matrix(BATH%Pb, UNIT)

      call which_basis_to_use(BATH)

      write (log_unit, *) ' GO WITH SUPERCONDUCTIVITY ? ', BATH%SUPER
      READ (UNIT, *) BATH%Niter_search_max
      READ (UNIT, *) BATH%search_step

      READ (UNIT, *) BATH%dist_max

      CALL close_safe(UNIT)

   end subroutine

   subroutine write_raw_bath(BATH)

      ! $$ WRITE RAW BATH PARAMETERS

      use masked_matrix_class, only: write_raw_masked_matrix
      use common_def, only: close_safe, open_safe
      use string5, only: get_unit

      implicit none

      TYPE(BATH_type), INTENT(IN) :: BATH
      INTEGER :: UNIT
      INTEGER :: spin

      CALL open_safe(UNIT, BATH%fileout, 'UNKNOWN', 'WRITE', get_unit=.true.)

      ! WRITE BATH

      WRITE (UNIT, *) BATH%Nb ! redundancy to allow preallocation in
      ! read_from_raw_file
      WRITE (UNIT, *) BATH%Nc ! redundancy to allow preallocation in
      ! read_from_raw_file
      WRITE (UNIT, *) BATH%fileout
      WRITE (UNIT, *) BATH%nparam
      WRITE (UNIT, *) SIZE(BATH%Eb)

      DO spin = 1, SIZE(BATH%Eb)
         CALL write_raw_masked_matrix(BATH%Eb(spin), UNIT)
      ENDDO
      WRITE (UNIT, *) SIZE(BATH%Vbc)
      DO spin = 1, SIZE(BATH%Vbc)
         CALL write_raw_masked_matrix(BATH%Vbc(spin), UNIT)
      ENDDO
      WRITE (UNIT, *) SIZE(BATH%PVbc)
      DO spin = 1, SIZE(BATH%PVbc)
         CALL write_raw_masked_matrix(BATH%PVbc(spin), UNIT)
      ENDDO
      CALL write_raw_masked_matrix(BATH%Pb, UNIT)
      WRITE (UNIT, *) BATH%Niter_search_max
      WRITE (UNIT, *) BATH%search_step
      WRITE (UNIT, *) BATH%dist_max
      CALL close_safe(UNIT)

   end subroutine

   subroutine read_bath(BATH, FILEIN, Nc)

      use genvar, only: cspin, log_unit
      use globalvar_ed_solver, only: bath_nearest_hop, diag_bath, dist_max, &
         fmos, force_no_bcs_pairing, force_pairing_from_mask, &
         force_para_state, force_singlet_state, istati, min_all_bath_param, &
         niter_search_max, para_state, pairing_imp_to_bath, search_step, &
         start_para
      use masked_matrix_class, only: clean_redundant_imask, &
         fill_masked_matrix, new_masked_matrix, test_masked_matrix_hermitic
      use common_def, only: close_safe, open_safe, skip_line
      use matrix, only: write_array
      use random, only: random_float_from_interval, random_complex_from_interval
      use string5, only: get_unit

      implicit none

      TYPE(bath_type), INTENT(INOUT) :: BATH
      CHARACTER(LEN=*), INTENT(IN) :: FILEIN
      INTEGER, INTENT(IN)            :: Nc
      INTEGER, ALLOCATABLE :: IMASKE(:, :, :), IMASKV(:, :, :), IMASKP(:, :), &
                              IMASKPV(:, :, :)
      INTEGER              :: Nb, index_in, ib, jb, spin, mu, iind_, iind, &
                              fac_cplx, nind
      INTEGER              :: Nitermax, ff
      REAL(DP)            :: step, dmax
      CHARACTER(LEN=5)   :: sym_bath
      INTEGER              :: UNIT, UNIT2, k
#ifdef _complex
      COMPLEX(DP)         :: val
#else
      REAL(DP)            :: val
#endif

      if (min_all_bath_param == 0) then
         CALL open_safe(UNIT, FILEIN, 'UNKNOWN', 'READ', get_unit=.true.)
         CALL skip_line(UNIT, 1)
         READ (UNIT, *) Nb
      else
         Nb = abs(min_all_bath_param)
      endif

      if (Nb < 1) stop 'error number of site in the bath'

      BATH%Niter_search_max = Niter_search_max ! max # of iterations in search routine (conj. grad.)
      BATH%search_step = search_step ! small step in initial search direction
      BATH%dist_max = dist_max ! max. error on hybridization functions

      CALL new_bath_from_scratch(BATH, Nb, Nc)

      IF (ALLOCATED(IMASKE)) DEALLOCATE (IMASKE, STAT=istati)
      ALLOCATE (IMASKE(Nb, Nb, 2))
      IMASKE = 0

      if (min_all_bath_param == 0) then
         CALL skip_line(UNIT, 3)
         DO spin = 1, 2
            CALL skip_line(UNIT, 1)
            READ (UNIT, *) ((IMASKE(ib, jb, spin), jb=1, Nb), ib=1, Nb)
         ENDDO
      else
         do spin = 1, 2
            if (para_state .or. spin == 1 .or. force_para_state) k = 1
            do ib = 1, Nb
               do jb = ib, Nb
                  if ((.not. diag_bath .or. ib == jb) .or. (bath_nearest_hop &
                                                            .and. abs(ib - jb) <= 1)) then
                     IMASKE(ib, jb, spin) = k
                     k = k + 1
                  endif
               enddo
            enddo
         enddo
      endif
      call write_array(IMASKE, ' MASK BATH PARAMETERS EPSILON', unit= &
                       log_unit)

      IF (ASSOCIATED(BATH%Eb)) DEALLOCATE (BATH%Eb, STAT=istati)
      ALLOCATE (BATH%Eb(2))

      DO spin = 1, 2
         CALL new_masked_matrix(BATH%Eb(spin), "Eb(sz = "//TRIM(cspin(spin)) &
                                //")", Nb, Nb, IMASK=IMASKE(:, :, spin), IS_HERM=.true.)
      ENDDO
      CALL clean_redundant_imask(BATH%Eb)

      if (min_all_bath_param == 0) then
         CALL skip_line(UNIT, 1)
         DO iind_ = 1, BATH%Eb(1)%rc%MASK%nind + BATH%Eb(2)%rc%MASK%nind
            READ (UNIT, *) iind, val
            DO spin = 1, 2
               CALL fill_masked_matrix(BATH%Eb(spin), iind, val)
            ENDDO
         ENDDO
      else
         if (start_para) then
            nind = BATH%Eb(1)%rc%MASK%nind
         else
            nind = BATH%Eb(1)%rc%MASK%nind + BATH%Eb(2)%rc%MASK%nind
         end if
         DO iind_ = 1, nind
            iind = iind_
#ifdef _complex
            val = random_complex_from_interval(-0.5d0, 0.0d0, 0.5d0, 1.0d0, same_across_tasks=.true.)
#else
            val = random_float_from_interval(-0.5d0, 0.5d0, same_across_tasks=.true.)
#endif
            DO spin = 1, 2
               CALL fill_masked_matrix(BATH%Eb(spin), iind, val)
               if (start_para) call fill_masked_matrix(BATH%Eb(spin), iind + BATH%Eb(1)%rc%MASK%nind, val)
            ENDDO
         ENDDO
      endif
      call write_array(BATH%Eb(1)%rc%mat, ' BATH PARAMETERS EPSILON UP spin', &
                       unit=log_unit, short=.true.)
      call write_array(BATH%Eb(2)%rc%mat, ' BATH PARAMETERS EPSILON DN spin', &
                       unit=log_unit, short=.true.)

      DO spin = 1, 2
         write (*, *) 'test hermiticity bath spin :', spin
         CALL test_masked_matrix_hermitic(BATH%Eb(spin))
      ENDDO

      write (*, *) '====================================================='
      write (*, *) 'read baths . There are ... sites in the bath : ', Nb
      write (*, *) '====================================================='

      IF (ALLOCATED(IMASKP)) DEALLOCATE (IMASKP, STAT=istati)
      ALLOCATE (IMASKP(Nb, Nb))
      IMASKP = 0

      if (min_all_bath_param == 0 .or. force_pairing_from_mask) then
         if (min_all_bath_param /= 0 .and. force_pairing_from_mask) then
            CALL open_safe(UNIT2, './ED/pairing', 'UNKNOWN', 'READ', get_unit &
                           =.true.)
            READ (UNIT2, *) ((IMASKP(ib, jb), jb=1, Nb), ib=1, Nb)
            call close_safe(UNIT2)
         else
            CALL skip_line(UNIT, 4)
            READ (UNIT, *) ((IMASKP(ib, jb), jb=1, Nb), ib=1, Nb)
         endif
      elseif (min_all_bath_param < 0) then
         if (.not. force_no_bcs_pairing) then

            if (.not. force_singlet_state) then
               k = 1
               do ib = 1, Nb
                  do jb = 1, Nb
                     IMASKP(ib, jb) = k
                     k = k + 1
                  enddo
               enddo
            endif
            if (force_singlet_state) then
               k = 1
               DO ib = 1, Nb
                  IMASKP(ib, ib) = 0
                  do jb = 1, ib - 1
                     IMASKP(jb, ib) = k
                     IMASKP(ib, jb) = IMASKP(jb, ib)
                     k = k + 1
                  enddo
               enddo
            endif
         endif
      endif

      CALL new_masked_matrix(BATH%Pb, "Pb", Nb, Nb, IMASK=IMASKP)
      call which_basis_to_use(BATH)

      if (.not. force_no_bcs_pairing) then
         IF (BATH%SUPER) THEN
            IF (BATH%Pb%rc%MASK%nind > 0) then
               if (min_all_bath_param == 0) then
                  write (*, *) 'read bath parameters ', BATH%Pb%rc%MASK%nind
                  CALL skip_line(UNIT, 1)
                  DO iind_ = 1, BATH%Pb%rc%MASK%nind
                     READ (UNIT, *) iind, val
                     CALL fill_masked_matrix(BATH%Pb, iind, val)
                  ENDDO
               elseif (min_all_bath_param < 0) then
                  DO iind_ = 1, BATH%Pb%rc%MASK%nind
                     iind = iind_
#ifdef _complex
                     val = random_complex_from_interval(-0.5d0, 0.0d0, 0.5d0, 1.0d0, same_across_tasks=.true.)
#else
                     val = random_float_from_interval(-0.5d0, 0.5d0, same_across_tasks=.true.)
#endif
                     CALL fill_masked_matrix(BATH%Pb, iind, val)
                  ENDDO
               endif
            ENDIF
         ENDIF
      endif

      write (log_unit, *) 'MIN ALL BATH PARAM = ', min_all_bath_param
      call write_array(BATH%Pb%rc%MASK%imat, ' RVB PARAMETERS ', unit= &
                       log_unit)
      call write_array(BATH%Pb%rc%mat, ' RVB PARAMETERS ', unit=log_unit)

      IF (ALLOCATED(IMASKV)) DEALLOCATE (IMASKV, STAT=istati)
      ALLOCATE (IMASKV(Nb, Nc, 2))
      IMASKV = 0
      IF (ALLOCATED(IMASKPV)) DEALLOCATE (IMASKPV, STAT=istati)
      ALLOCATE (IMASKPV(Nb, Nc, 2))
      IMASKPV = 0

      if (min_all_bath_param == 0) then
         CALL skip_line(UNIT, 3)
         DO spin = 1, 2
            CALL skip_line(UNIT, 1)
            READ (UNIT, *) ((IMASKV(ib, mu, spin), ib=1, Nb), mu=1, Nc)
         ENDDO
      else
         do spin = 1, 2
            if (para_state .or. spin == 1 .or. force_para_state) k = 1
            do ib = 1, Nb
               do mu = 1, Nc
                  if (fmos) Then
                     if (.not. (ib > (mu - 1)*Nb/Nc .and. ib <= mu*Nb/Nc)) &
                        cycle
                  endif
                  IMASKV(ib, mu, spin) = k
                  k = k + 1
               enddo
            enddo
         enddo
      endif
      call write_array(IMASKV, ' MASK BATH PARAMETERS Vbc ', unit=log_unit)

      if (BATH%SUPER) then
         if (PAIRING_IMP_TO_BATH) then
            if (min_all_bath_param == 0) then
               DO spin = 1, 2
                  CALL skip_line(UNIT, 1)
                  READ (UNIT, *) ((IMASKPV(ib, mu, spin), ib=1, Nb), mu=1, &
                                  Nc)
               ENDDO
            else
               if (.not. force_no_bcs_pairing) then
                  k = 1
                  do spin = 1, 2
                     do ib = 1, Nb
                        do mu = 1, Nc
                           IMASKPV(ib, mu, spin) = k
                           k = k + 1
                        enddo
                     enddo
                  enddo
                  if (force_singlet_state) then
                     IMASKPV(:, :, 2) = IMASKPV(:, :, 1)
                  endif
               endif
               call write_array(IMASKPV, ' MASK BATH PARAMETERS PVbc ', unit &
                                =log_unit)
            endif
         endif
      endif

      IF (ASSOCIATED(BATH%Vbc)) DEALLOCATE (BATH%Vbc, STAT=istati)
      ALLOCATE (BATH%Vbc(2))
      DO spin = 1, 2
         CALL new_masked_matrix(BATH%Vbc(spin), "Vbc(sz = "// &
                                TRIM(cspin(spin))//")", Nb, Nc, IMASK=IMASKV(:, :, spin))
      ENDDO
      CALL clean_redundant_imask(BATH%Vbc)

      IF (ASSOCIATED(BATH%PVbc)) DEALLOCATE (BATH%PVbc, STAT=istati)
      ALLOCATE (BATH%PVbc(2))
      DO spin = 1, 2
         CALL new_masked_matrix(BATH%PVbc(spin), "PVbc(sz = "// &
                                TRIM(cspin(spin))//")", Nb, Nc, IMASK=IMASKPV(:, :, spin))
      ENDDO
      CALL clean_redundant_imask(BATH%PVbc)

      if (BATH%Vbc(1)%rc%MASK%nind + BATH%Vbc(2)%rc%MASK%nind > 0) then
         if (min_all_bath_param == 0) then
            CALL skip_line(UNIT, 1)
            DO iind_ = 1, BATH%Vbc(1)%rc%MASK%nind + BATH%Vbc(2)%rc%MASK%nind
               READ (UNIT, *) iind, val
               DO spin = 1, 2
                  CALL fill_masked_matrix(BATH%Vbc(spin), iind, val)
               ENDDO
            ENDDO
         else
            if (start_para) then
               nind = BATH%Eb(1)%rc%MASK%nind
            else
               nind = BATH%Eb(1)%rc%MASK%nind + BATH%Eb(2)%rc%MASK%nind
            end if
            DO iind_ = 1, nind
               iind = iind_
#ifdef _complex
               val = random_complex_from_interval(-0.5d0, 0.0d0, 0.5d0, 1.0d0, same_across_tasks=.true.)
#else
               val = random_float_from_interval(-0.5d0, 0.5d0, same_across_tasks=.true.)
#endif
               DO spin = 1, 2
                  CALL fill_masked_matrix(BATH%Vbc(spin), iind, val)
                  if (start_para) call fill_masked_matrix(BATH%Vbc(spin), iind + BATH%Vbc(1)%rc%MASK%nind, val)
               ENDDO
            ENDDO
         endif
      endif
      call write_array(BATH%Vbc(1)%rc%mat, ' BATH PARAMETERS Vbc UP spin', &
                       unit=log_unit, short=.true.)
      call write_array(BATH%Vbc(2)%rc%mat, ' BATH PARAMETERS Vbc DN spin', &
                       unit=log_unit, short=.true.)

      if (.not. force_no_bcs_pairing) then
         if (BATH%SUPER) then
            if (PAIRING_IMP_TO_BATH) then
               if (BATH%PVbc(1)%rc%MASK%nind + BATH%PVbc(2)%rc%MASK%nind > 0) then
                  if (min_all_bath_param == 0) then
                     CALL skip_line(UNIT, 1)
                     DO iind_ = 1, BATH%PVbc(1)%rc%MASK%nind + &
                        BATH%PVbc(2)%rc%MASK%nind
                        READ (UNIT, *) iind, val
                        DO spin = 1, 2
                           CALL fill_masked_matrix(BATH%PVbc(spin), iind, val)
                        ENDDO
                     ENDDO
                     CALL close_safe(UNIT)
                  else
                     DO iind_ = 1, BATH%PVbc(1)%rc%MASK%nind + &
                        BATH%PVbc(2)%rc%MASK%nind
                        iind = iind_
#ifdef _complex
                        val = random_complex_from_interval(-0.5d0, 0.0d0, 0.5d0, 1.0d0, same_across_tasks=.true.)
#else
                        val = random_float_from_interval(-0.5d0, 0.5d0, same_across_tasks=.true.)
#endif
                        DO spin = 1, 2
                           CALL fill_masked_matrix(BATH%PVbc(spin), iind, val)
                        ENDDO
                     ENDDO
                  endif
               endif
               call write_array(BATH%PVbc(1)%rc%mat, ' BATH PARAMETERS PVbc &
                    &UP spin', unit=log_unit, short=.true.)
               call write_array(BATH%PVbc(2)%rc%mat, ' BATH PARAMETERS PVbc &
                    &DN spin', unit=log_unit, short=.true.)
            endif
         endif
      endif

      index_in = INDEX(FILEIN, '.in')
      BATH%fileout = FILEIN(1:index_in - 1)//'.out'
      BATH%nparam = nbathparam(BATH)

      IF (ALLOCATED(IMASKE)) DEALLOCATE (IMASKE, STAT=istati)
      IF (ALLOCATED(IMASKV)) DEALLOCATE (IMASKV, STAT=istati)
      IF (ALLOCATED(IMASKP)) DEALLOCATE (IMASKP, STAT=istati)
      IF (ALLOCATED(IMASKPV)) DEALLOCATE (IMASKPV, STAT=istati)

   end subroutine

   function nbathparam(bath)

      INTEGER                     :: nbathparam
      TYPE(bath_type), INTENT(IN) :: BATH
      INTEGER                     :: spin

      write (*, *) '------------------------------------------------------------'
      write (*, *) ' nbathparam : Pb, Eb, Vbc '
      write (*, *) BATH%Pb%rc%MASK%nind, BATH%Eb%rc%MASK%nind, BATH%Vbc%rc%MASK%nind
      write (*, *) '------------------------------------------------------------'

#ifdef _complex
      nbathparam = offset__(BATH%Pb)
      do spin = 1, 2
         nbathparam = nbathparam + offset__(BATH%Eb(spin)) + offset__(BATH%Vbc(spin))
      enddo
      do spin = 1, 2
         nbathparam = nbathparam + offset__(BATH%PVbc(spin))
      enddo
#else
      nbathparam = BATH%Pb%rc%MASK%nind
      DO spin = 1, 2
         nbathparam = nbathparam + BATH%Eb(spin)%rc%MASK%nind + &
                      BATH%Vbc(spin)%rc%MASK%nind
      ENDDO
      do spin = 1, 2
         nbathparam = nbathparam + BATH%PVbc(spin)%rc%MASK%nind
      enddo
#endif

   contains
      integer function offset__(MM)

         use masked_matrix_class, only: masked_matrix_type

         implicit none

         integer                  :: ninddiag, nindoffdiag, nind
         type(masked_matrix_type) :: MM

         offset__ = 0

#ifdef _complex
         ninddiag = MM%rc%MASKdiag%nind
         nindoffdiag = MM%rc%MASKoffdiag%nind
#endif
         nind = MM%rc%MASK%nind

         if (nind /= 0) then
            IF (.NOT. MM%rc%IS_HERM) THEN
               offset__ = offset__ + 2*nind
            ELSE
               offset__ = offset__ + ninddiag
               offset__ = offset__ + 2*nindoffdiag
            ENDIF
         endif
      end function

   end function

   function nbathparam_(BATH)

      ! COMPUTE NUMBER OF BATH PARAMETERS

      implicit none

      TYPE(bath_type), INTENT(IN) :: BATH
      INTEGER :: nbathparam_
      INTEGER :: spin

      write (*, *) '... compute number of bath parameters : '

      nbathparam_ = BATH%Pb%rc%MASK%nind
      DO spin = 1, 2
         nbathparam_ = nbathparam_ + BATH%Eb(spin)%rc%MASK%nind + &
                       BATH%Vbc(spin)%rc%MASK%nind
      ENDDO
      DO spin = 1, 2
         nbathparam_ = nbathparam_ + BATH%PVbc(spin)%rc%MASK%nind
      ENDDO

#ifdef _complex
      nbathparam_ = nbathparam_ + BATH%Pb%rc%MASK%nind
      DO spin = 1, 2
         nbathparam_ = nbathparam_ + BATH%Eb(spin)%rc%MASKoffdiag%nind + &
                       BATH%Vbc(spin)%rc%MASKoffdiag%nind
      ENDDO
      DO spin = 1, 2
         nbathparam_ = nbathparam_ + BATH%PVbc(spin)%rc%MASKoffdiag%nind
      ENDDO
#endif

   end function

   subroutine write_bath(BATH, UNIT, FILEIN)

      use common_def, only: dump_message
      use genvar, only: log_unit
      use masked_matrix_class, only: write_masked_matrix

      implicit none

      TYPE(bath_type), INTENT(IN)              :: BATH
      INTEGER, OPTIONAL, INTENT(IN)            :: UNIT
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FILEIN
      INTEGER :: unit_, spin

      IF (.NOT. ASSOCIATED(BATH%Eb)) STOP "ERROR IN write_bath: INPUT ISNT &
           &ASSOCIATED!"

      CALL dump_message(UNIT=UNIT, TEXT="#######################")
      CALL dump_message(UNIT=UNIT, TEXT="### BATH PARAMETERS ###")
      CALL dump_message(UNIT=UNIT, TEXT="#######################")

      IF (PRESENT(FILEIN)) CALL dump_message(UNIT=UNIT, TEXT="### READ FROM &
           &FILE "//TRIM(ADJUSTL(FILEIN)))

      unit_ = log_unit
      IF (PRESENT(UNIT)) unit_ = UNIT

      WRITE (unit_, '(a, I0)') "# NUMBER OF SITES IN THE IMPURITY = ", BATH%Nc
      WRITE (unit_, '(a, I0)') "# TOTAL NUMBER OF BATH PARAMETERS = ", &
         BATH%nparam

      ! BATH ENERGY
      DO spin = 1, 2
         CALL write_masked_matrix(BATH%Eb(spin), UNIT=UNIT)
      ENDDO

      ! PAIRING ENERGY IF SUPERCONDUCTING
      IF (BATH%SUPER) CALL write_masked_matrix(BATH%Pb, UNIT=UNIT)

      ! HYBRIDIZATION
      DO spin = 1, 2
         CALL write_masked_matrix(BATH%Vbc(spin), UNIT=UNIT)
      ENDDO

      DO spin = 1, 2
         CALL write_masked_matrix(BATH%PVbc(spin), UNIT=UNIT)
      ENDDO

      ! HYBRID => BATH PARAMETERS
      CALL dump_message(UNIT=UNIT, TEXT="# HYBRID => BATH FIT:")
      WRITE (unit_, '(a, I0, a, E10.3)') "# Nparam = ", BATH%nparam, " &
           &parameters, TOLERANCE = ", BATH%dist_max
      WRITE (unit_, '(a, I0, a, E10.3)') "# Niter = ", BATH%Niter_search_max, " &
           &iterations with STEP = ", BATH%search_step

      CALL flush(unit_)

   end subroutine

   subroutine Nambu_Eb(EbNambu, Eb, Pb)

      ! MAKE NAMBU ENERGY MATRIX

      use globalvar_ed_solver, only: energy_global_shift2
      use masked_matrix_class, only: masked_matrix_type, new_masked_matrix
      use matrix, only: diag

      implicit none

      TYPE(masked_matrix_type), INTENT(INOUT) :: EbNambu
      TYPE(masked_matrix_type), INTENT(IN)    :: Eb(:)
      TYPE(masked_matrix_type), INTENT(IN)    :: Pb
      INTEGER :: Nb ! to clarify only

      IF (SIZE(Eb) == 0) STOP "ERROR IN Nambu_Eb: INPUT Eb ISNT ALLOCATED!"
      IF (Pb%rc%n1 == 0) STOP "ERROR IN Nambu_Eb: INPUT Pb ISNT ALLOCATED!"

      Nb = Eb(1)%rc%n1
      CALL new_masked_matrix(EbNambu, "EbNambu", Nb*2, Nb*2, IS_HERM=.true.)
      EbNambu%rc%mat = 0.0_DP

      ! DIAGONAL BLOCKS: CONDUCTION ENERGY
      EbNambu%rc%mat(1:Nb, 1:Nb) = Eb(1)%rc%mat
      EbNambu%rc%mat(Nb + 1:Nb*2, Nb + 1:Nb*2) = -TRANSPOSE(Eb(2)%rc%mat)

      energy_global_shift2 = -sum(diag(EbNambu%rc%mat(Nb + 1:Nb*2, Nb + &
                                                      1:Nb*2)))

      ! OFF-DIAGONAL BLOCKS: PAIRING ENERGY
      EbNambu%rc%mat(1:Nb, Nb + 1:Nb*2) = Pb%rc%mat

#ifdef _complex
      EbNambu%rc%mat(Nb + 1:Nb*2, 1:Nb) = TRANSPOSE(CONJG(Pb%rc%mat))
#else
      EbNambu%rc%mat(Nb + 1:Nb*2, 1:Nb) = TRANSPOSE(Pb%rc%mat)
#endif

   end subroutine

   subroutine Nambu_Vbc(VbcNambu, Vbc, PVbc)

      ! MAKE NAMBU HYBRIDIZATION MATRIX

      use globalvar_ed_solver, only: pairing_imp_to_bath
      use masked_matrix_class, only: masked_matrix_type, new_masked_matrix

      implicit none

      TYPE(masked_matrix_type), INTENT(INOUT) :: VbcNambu
      TYPE(masked_matrix_type), INTENT(IN)    :: Vbc(:), PVbc(:)
      INTEGER :: Nc, Nb ! to clarify only

      IF (SIZE(Vbc) == 0) STOP "ERROR IN Nambu_Vbc: INPUT Vbc ISNT ALLOCATED!"

      Nb = Vbc(1)%rc%n1
      Nc = Vbc(1)%rc%n2

      CALL new_masked_matrix(VbcNambu, "VbcNambu", Nb*2, Nc*2)
      VbcNambu%rc%mat = 0.0_DP
      VbcNambu%rc%mat(1:Nb, 1:Nc) = Vbc(1)%rc%mat
#ifdef _complex
      VbcNambu%rc%mat(Nb + 1:Nb*2, Nc + 1:Nc*2) = -CONJG(Vbc(2)%rc%mat)
#else
      VbcNambu%rc%mat(Nb + 1:Nb*2, Nc + 1:Nc*2) = -Vbc(2)%rc%mat
#endif

      if (PAIRING_IMP_TO_BATH) then
         VbcNambu%rc%mat(Nb + 1:2*Nb, 1:Nc) = PVbc(1)%rc%mat
         VbcNambu%rc%mat(1:Nb, Nc + 1:2*Nc) = PVbc(2)%rc%mat
      endif

   end subroutine

end module
