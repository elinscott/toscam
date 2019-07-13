MODULE masked_matrix_class_mod

   use mask_class, only: mask_type
   use genvar, only: DP, DP, log_unit, FERMIONIC, BOSONIC

   IMPLICIT NONE

   private

   INTEGER                                :: istati

   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   !$$REAL MASKED MATRIX CLASS$$
   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   TYPE masked_real_matrix_type
      ! NAME OF MATRIX
      CHARACTER(LEN=100)       :: title = '\0'

      ! FULL 2D MATRIX
      INTEGER                  :: n1 = 0, n2 = 0
      ! REAL MATRIX
      REAL(DP), POINTER :: mat(:, :) => NULL()
      ! MASKS
      TYPE(mask_type)          :: MASK

      ! VECTOR OF INDEPENDANT MATRIX ELEMENTS
      REAL(DP), POINTER :: vec(:) => NULL()
      ! FLAGS
      LOGICAL                  :: IS_HERM = .false.    ! TRUE IF HERMITIC MATRIX
   END TYPE

   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   !$$COMPLEX MASKED MATRIX CLASS$$
   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   TYPE masked_cplx_matrix_type
      ! NAME OF MATRIX
      CHARACTER(LEN=100)       :: title = '\0'

      ! FULL 2D MATRIX
      INTEGER                  :: n1 = 0, n2 = 0
      ! COMPLEX MATRIX
      COMPLEX(DP), POINTER :: mat(:, :) => NULL()
      ! MASKS
      TYPE(mask_type)          :: MASK

      ! VECTOR OF INDEPENDANT MATRIX ELEMENTS
      COMPLEX(DP), POINTER :: vec(:) => NULL()

      ! ADDITIONAL MASKS/VECTORS IF HERMITIC
      LOGICAL                  :: IS_HERM = .false.    ! TRUE IF HERMITIC MATRIX
      TYPE(mask_type)          :: MASKdiag
      REAL(DP), POINTER ::  vecdiag(:) => NULL()
      TYPE(mask_type)          :: MASKoffdiag
      COMPLEX(DP), POINTER ::  vecoffdiag(:) => NULL()
   END TYPE

   INTERFACE new_masked_cplx_matrix
      MODULE PROCEDURE new_masked_cplx_matrix_from_scratch
      MODULE PROCEDURE new_masked_cplx_matrix_from_old
   END INTERFACE

   INTERFACE new_masked_real_matrix
      MODULE PROCEDURE new_masked_real_matrix_from_scratch
      MODULE PROCEDURE new_masked_real_matrix_from_old
   END INTERFACE

   INTERFACE copy_masked_matrix_
      MODULE PROCEDURE copy_masked_cplx_matrix, copy_masked_real_matrix
   END INTERFACE

   INTERFACE delete_masked_matrix_
      MODULE PROCEDURE delete_masked_cplx_matrix, delete_masked_real_matrix
   END INTERFACE

   INTERFACE write_raw_masked_matrix_
      MODULE PROCEDURE write_raw_masked_cplx_matrix
      MODULE PROCEDURE write_raw_masked_real_matrix
   END INTERFACE

   INTERFACE read_raw_masked_matrix_
      MODULE PROCEDURE read_raw_masked_cplx_matrix, read_raw_masked_real_matrix
   END INTERFACE

   INTERFACE write_masked_matrix_
      MODULE PROCEDURE write_masked_cplx_matrix, write_masked_real_matrix
   END INTERFACE

   INTERFACE masked_matrix2vec_
      MODULE PROCEDURE masked_cplx_matrix2vec, masked_real_matrix2vec
   END INTERFACE

   INTERFACE vec2masked_matrix_
      MODULE PROCEDURE vec2masked_cplx_matrix, vec2masked_real_matrix
   END INTERFACE

   INTERFACE fill_masked_matrix_
      MODULE PROCEDURE fill_masked_cplx_matrix, fill_masked_real_matrix
   END INTERFACE

   INTERFACE build_mask_
      MODULE PROCEDURE build_mask_cplx, build_mask_real
   END INTERFACE

   INTERFACE test_masked_matrix_hermitic_
      MODULE PROCEDURE test_masked_cplx_matrix_hermitic
   END INTERFACE

   INTERFACE pad_masked_matrix_
      MODULE PROCEDURE pad_masked_cplx_matrix, pad_masked_real_matrix
   END INTERFACE

   INTERFACE filter_masked_matrix_
      MODULE PROCEDURE filter_masked_cplx_matrix, filter_masked_real_matrix
   END INTERFACE

   INTERFACE slice_masked_matrix_
      MODULE PROCEDURE slice_masked_cplx_matrix, slice_masked_real_matrix
   END INTERFACE

   INTERFACE clean_redundant_imask_
      MODULE PROCEDURE clean_redundant_imask_cplx, clean_redundant_imask_real
   END INTERFACE

   INTERFACE transform_masked_matrix_
      MODULE PROCEDURE transform_masked_cplx_matrix
      MODULE PROCEDURE transform_masked_real_matrix
   END INTERFACE

   public :: build_mask_
   public :: clean_redundant_imask_
   public :: copy_masked_matrix_
   public :: copy_masked_real_matrix
   public :: delete_masked_cplx_matrix
   public :: delete_masked_matrix_
   public :: delete_masked_real_matrix
   public :: fill_masked_matrix_
   public :: fill_masked_real_matrix
   public :: filter_masked_matrix_
   public :: gather_diag_offdiag_vec
   public :: masked_cplx_matrix2vec
   public :: masked_cplx_matrix_type
   public :: masked_matrix2vec_
   public :: masked_real_matrix_type
   public :: new_masked_cplx_matrix
   public :: new_masked_cplx_matrix_from_old
   public :: new_masked_cplx_matrix_from_scratch
   public :: new_masked_real_matrix
   public :: new_masked_real_matrix_from_old
   public :: new_masked_real_matrix_from_scratch
   public :: pad_masked_matrix_
   public :: read_raw_masked_cplx_matrix
   public :: read_raw_masked_matrix_
   public :: slice_masked_cplx_matrix
   public :: slice_masked_matrix_
   public :: test_masked_cplx_matrix_hermitic
   public :: test_masked_real_matrix_symmetric
   public :: transform_masked_matrix_
   public :: vec2masked_cplx_matrix
   public :: vec2masked_matrix_
   public :: write_masked_matrix_
   public :: write_masked_real_matrix
   public :: write_raw_masked_cplx_matrix
   public :: write_raw_masked_matrix_

contains

   subroutine new_masked_cplx_matrix_from_scratch(MM, title, n1, n2, IMASK, &
                                                  IS_HERM)

      ! CREATE MASKED MATRIX MM FROM SCRATCH

      use mask_class, only: new_mask, HERM

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MM
      CHARACTER(LEN=*), INTENT(IN)               :: title
      INTEGER, INTENT(IN)                          :: n1, n2
      LOGICAL, OPTIONAL, INTENT(IN)                :: IS_HERM
      INTEGER, OPTIONAL, INTENT(IN)                :: IMASK(n1, n2)
      CHARACTER(LEN=9) :: SYMMETRY

      CALL delete_masked_cplx_matrix(MM)

      MM%title = title(1:MIN(100, LEN_TRIM(title)))

      MM%IS_HERM = .false.
      IF (PRESENT(IS_HERM)) THEN
         IF (IS_HERM .AND. n1 /= n2) then
            write (*, *) "ERROR in new_masked_cplx_matrix: "// &
               TRIM(ADJUSTL(title))//" ISNT SQUARE! IS_HERM IRRELEVANT!"
            stop
         endif
         MM%IS_HERM = IS_HERM
      ENDIF

      MM%n1 = n1
      MM%n2 = n2
      ALLOCATE (MM%mat(n1, n2))
      MM%mat = 0.0_DP

      SYMMETRY = ''
      IF (MM%IS_HERM) SYMMETRY = HERM
      CALL new_mask(MM%MASK, n1, n2, IMASK=IMASK, SYMMETRY=SYMMETRY)

      IF (MM%IS_HERM) THEN
         ! DIAGONAL
         CALL new_mask(MM%MASKdiag, n1, n2, SYMMETRY=SYMMETRY)
         ! OFF-DIAGONAL
         CALL new_mask(MM%MASKoffdiag, n1, n2, SYMMETRY=SYMMETRY)
      ENDIF

      CALL build_mask_cplx(MM)

   end subroutine

   subroutine new_masked_cplx_matrix_from_old(MMOUT, MMIN)

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MMOUT
      TYPE(masked_cplx_matrix_type), INTENT(IN)    :: MMIN

      CALL delete_masked_cplx_Matrix(MMOUT)
      IF (.NOT. ASSOCIATED(MMIN%mat)) then
         write (*, *) "ERROR IN new_masked_cplx_matrix_from_old: INPUT ISNT &
              &ALLOCATED!"
         stop
      endif
      CALL new_masked_cplx_matrix_from_scratch(MMOUT, MMIN%title, MMIN%n1, &
                                               MMIN%n2, IMASK=MMIN%MASK%imat, IS_HERM=MMIN%IS_HERM)
      CALL copy_masked_cplx_matrix(MMOUT, MMIN)
   end subroutine

   subroutine copy_masked_cplx_matrix(MMOUT, MMIN)

      ! CREATE MMOUT BY COPYING EXISTING MMIN

      use mask_class, only: copy_mask

      implicit none

      TYPE(masked_cplx_matrix_type) :: MMOUT
      TYPE(masked_cplx_matrix_type) :: MMIN
      integer                       :: i, j

      IF (.NOT. ASSOCIATED(MMIN%mat)) STOP "ERROR IN copy_masked_cplx_matrix: &
           &INPUT ISNT ALLOCATED!"

      MMOUT%IS_HERM = MMIN%IS_HERM

      ! FULL MATRIX
      ! SIZE
      MMOUT%n1 = MMIN%n1
      MMOUT%n2 = MMIN%n2

      ! MASKS
      CALL copy_mask(MMOUT%MASK, MMIN%MASK)

      ! MATRIX
      IF (.NOT. ASSOCIATED(MMIN%mat)) STOP "ERROR IN copy_masked_cplx_matrix: &
           &INPUT ISNT ALLOCATED!"
      IF (.NOT. ASSOCIATED(MMOUT%mat)) then
         ALLOCATE (MMOUT%mat(MMOUT%n1, MMOUT%n2))
      else
         if (ANY(shape(MMOUT%mat) .ne. shape(MMIN%mat))) then
            DEALLOCATE (MMOUT%mat, STAT=istati)
            ALLOCATE (MMOUT%mat(MMOUT%n1, MMOUT%n2))
         endif
      endif

      MMOUT%mat = MMIN%mat

      ! VECTOR
      IF (ASSOCIATED(MMIN%vec)) THEN
         IF (ASSOCIATED(MMOUT%vec) .AND. SIZE(MMOUT%vec) /= MMIN%MASK%nind) &
            DEALLOCATE (MMOUT%vec, STAT=istati)
         IF (.NOT. ASSOCIATED(MMOUT%vec)) ALLOCATE (MMOUT%vec(MMIN%MASK%nind))
         MMOUT%vec = MMIN%vec
      ENDIF

      ! ADDITIONAL MASKS IF IS_HERM = T
      IF (MMOUT%IS_HERM) THEN
         CALL copy_mask(MMOUT%MASKdiag, MMIN%MASKdiag)
         CALL copy_mask(MMOUT%MASKoffdiag, MMIN%MASKoffdiag)
         IF (ASSOCIATED(MMIN%vecdiag)) THEN
            IF (ASSOCIATED(MMOUT%vecdiag) .AND. SIZE(MMOUT%vecdiag) /= &
                SIZE(MMIN%vecdiag)) DEALLOCATE (MMOUT%vecdiag, STAT=istati)
            IF (.NOT. ASSOCIATED(MMOUT%vecdiag)) &
               ALLOCATE (MMOUT%vecdiag(MMIN%MASKdiag%nind))
            MMOUT%vecdiag = MMIN%vecdiag
         ENDIF
         IF (ASSOCIATED(MMIN%vecoffdiag)) THEN
            IF (ASSOCIATED(MMOUT%vecoffdiag) .AND. SIZE(MMOUT%vecoffdiag) /= &
                SIZE(MMIN%vecoffdiag)) DEALLOCATE (MMOUT%vecoffdiag, STAT= &
                                                   istati)
            IF (.NOT. ASSOCIATED(MMOUT%vecoffdiag)) &
               ALLOCATE (MMOUT%vecoffdiag(MMIN%MASKoffdiag%nind))
            MMOUT%vecoffdiag = MMIN%vecoffdiag
         ENDIF
      ENDIF

      return
   end subroutine

   subroutine delete_masked_cplx_matrix(MM)

      use mask_class, only: delete_mask

      implicit none

      TYPE(masked_cplx_matrix_type) :: MM

      IF (ASSOCIATED(MM%mat)) DEALLOCATE (MM%mat, STAT=istati)
      CALL delete_mask(MM%MASK)
      IF (ASSOCIATED(MM%vec)) DEALLOCATE (MM%vec, STAT=istati)
      CALL delete_mask(MM%MASKdiag)
      IF (ASSOCIATED(MM%vecdiag)) DEALLOCATE (MM%vecdiag, STAT=istati)
      CALL delete_mask(MM%MASKoffdiag)
      IF (ASSOCIATED(MM%vecoffdiag)) DEALLOCATE (MM%vecoffdiag, STAT=istati)
      IF (ASSOCIATED(MM%mat)) NULLIFY (MM%mat)
      IF (ASSOCIATED(MM%vec)) NULLIFY (MM%vec)
      IF (ASSOCIATED(MM%vecdiag)) NULLIFY (MM%vecdiag)
      IF (ASSOCIATED(MM%vecoffdiag)) NULLIFY (MM%vecoffdiag)
   end subroutine

   subroutine masked_cplx_matrix2vec(MM)

      ! PACK MATRIX INTO VECTOR (CREATED IF NECESSARY) ACCORDING TO MASKS

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MM

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN masked_cplx_matrix2vec: INPUT &
           &ISNT ALLOCATED!"

      ! INTIALIZE VECTORS IF NECESSARY
      IF (.NOT. ASSOCIATED(MM%vec)) THEN ! VECTOR OF INDPDT ELEMTS ISNT YET
         ! INITIALIZED
         ! nind IS SUPPOSED ALREADY COMPUTED
         IF (MM%MASK%nind /= 0) ALLOCATE (MM%vec(MM%MASK%nind))
         IF (MM%IS_HERM) THEN
            ! DIAGONAL PART ONLY
            IF (MM%MASKdiag%nind /= 0) ALLOCATE (MM%vecdiag(MM%MASKdiag%nind))
            ! OFF-DIAGONAL PART ONLY
            IF (MM%MASKoffdiag%nind /= 0) &
               ALLOCATE (MM%vecoffdiag(MM%MASKoffdiag%nind))
         ENDIF
      ENDIF

      ! PACK MATRIX ELEMENTS
      if (ANY(MM%MASK%mat)) MM%vec = PACK(MM%mat, MM%MASK%mat)
      IF (MM%IS_HERM) THEN
         if (MM%MASKdiag%nind /= 0) MM%vecdiag = PACK(MM%mat, &
                                                      MM%MASKdiag%mat)
         if (MM%MASKoffdiag%nind /= 0) MM%vecoffdiag = PACK(MM%mat, &
                                                            MM%MASKoffdiag%mat)
      ENDIF

   end subroutine

   subroutine vec2masked_cplx_matrix(MM, MMEXT)

      ! UNPACK MATRIX INTO VECTOR (ASSUMED TO EXIST) ACCORDING TO MASKS

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT)        :: MM
      TYPE(masked_cplx_matrix_type), INTENT(IN), OPTIONAL :: MMEXT
      INTEGER :: iind

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN vec2masked_cplx_matrix: INPUT &
           &ISNT ALLOCATED!"

      IF (PRESENT(MMEXT)) THEN
         DO iind = 1, MMEXT%MASK%nind
            CALL fill_masked_cplx_matrix(MM, MMEXT%MASK%ivec(iind), &
                                         MMEXT%vec(iind))
         ENDDO
      ELSE
         if (size(MM%vec) < MM%MASK%nind .or. size(MM%MASK%ivec) < &
             MM%MASK%nind) then
            write (*, *) 'vec2masked_cplx error : '
            write (*, *) 'size mm%vec and mm%ivec : ', size(MM%vec), &
               size(MM%MASK%ivec)
            write (*, *) 'nind : ', MM%MASK%nind
            stop
         endif
         DO iind = 1, MM%MASK%nind ! NOTHING HAPPENS IF nind = 0
            CALL fill_masked_cplx_matrix(MM, MM%MASK%ivec(iind), MM%vec(iind))
         ENDDO
      ENDIF

      IF (MM%IS_HERM) then
         CALL test_masked_cplx_matrix_hermitic(MM)
      ENDIF

   end subroutine

   subroutine fill_masked_cplx_matrix(MM, iind, val)

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MM
      INTEGER, INTENT(IN)                          :: iind
      COMPLEX(DP) :: val
      INTEGER      :: i1, i2

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN fill_masked_cplx_matrix: &
           &INPUT ISNT ALLOCATED!"
      DO i1 = 1, MM%n1
         DO i2 = 1, MM%n2
            IF (MM%MASK%imat(i1, i2) == iind) then
               if (i1 == i2 .and. MM%IS_HERM) val = real(val)
               MM%mat(i1, i2) = val
            endif
            IF (MM%MASK%imat(i1, i2) == -iind) then
               if (i1 == i2 .and. MM%IS_HERM) val = real(val)
               MM%mat(i1, i2) = -val
            endif
         ENDDO
      ENDDO
      ! SYMMETRIZE
      IF (MM%IS_HERM) THEN
         DO i1 = 1, MM%n1
            DO i2 = 1, MM%n2
               IF (MM%MASK%imat(i1, i2) == 0 .AND. MM%MASK%imat(i2, i1) /= 0) &
                  MM%mat(i1, i2) = CONJG(MM%mat(i2, i1))
            ENDDO

         ENDDO
      ENDIF
   end subroutine

   subroutine gather_diag_offdiag_vec(MM)

      ! GATHER VECTORS OF INDPDT DIAGONAL AND OFF-DIAGONAL ELEMTS IN A SINGLE
      ! ONE

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MM
      INTEGER :: iind

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN gather_diag_offdiag_vec: &
           &INPUT ISNT ALLOCATED!"

      ! DIAGONAL PART
      DO iind = 1, MM%MASKdiag%nind
         WHERE (MM%MASK%ivec == MM%MASKdiag%ivec(iind))
         MM%vec = MM%vecdiag(iind)
         END WHERE
      ENDDO

      ! OFF-DIAGONAL PART
      DO iind = 1, MM%MASKoffdiag%nind
         WHERE (MM%MASK%ivec == MM%MASKoffdiag%ivec(iind))
         MM%vec = MM%vecoffdiag(iind)
         END WHERE
      ENDDO

   end subroutine

   subroutine build_mask_cplx(MM, IMASK)

      ! BUILD LOGICAL MASK FROM INTEGER MASK

      use mask_class, only: build_logical_mask, filter_mask
      use matrix, only: new_diag

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MM
      INTEGER, OPTIONAL, INTENT(IN)                :: IMASK(:, :)
      LOGICAL :: is_diag(MM%n1, MM%n1)

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN build_mask_cplx: INPUT ISNT &
           &ALLOCATED!"
      CALL build_logical_mask(MM%MASK, IMASK=IMASK)
      IF (MM%IS_HERM) THEN
         ! WE ALREADY HAVE A iELL BEHAVED MASK IN MM%MASK
         CALL new_diag(is_diag, MM%n1)
         ! BUILD DIAGONAL MASK
         CALL filter_mask(MM%MASKdiag, MM%MASK, is_diag)
         ! BUILD OFF-DIAGONAL MASK
         CALL filter_mask(MM%MASKoffdiag, MM%MASK,.NOT. is_diag)
      ENDIF

   end subroutine

   subroutine test_masked_cplx_matrix_hermitic(MM)

      ! TEST HERMITICITY OF MATRIX

      use common_def, only: c2s, dump_message, i2c
      use matrix, only: write_array

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(IN) :: MM
      LOGICAL :: is_cplx_hermitic

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN &
           &test_masked_cplx_matrix_hermitic: INPUT ISNT ALLOCATED!"

      IF (MM%n1 == MM%n2) THEN
         is_cplx_hermitic = ALL(abs(MM%mat - TRANSPOSE(CONJG(MM%mat))) < 1.d-8)
      ELSE
         CALL dump_message(TEXT="ERROR IN test_masked_cplx_matrix_hermitic: &
              &MATRIX "//TRIM(ADJUSTL(MM%title))//" ISNT SQUARE!")
         CALL dump_message(TEXT="n1 = "//c2s(i2c(MM%n1))//" n2 = "// &
                           c2s(i2c(MM%n2)))
      ENDIF

      IF (.NOT. is_cplx_hermitic) THEN
         CALL dump_message(TEXT="ERROR IN test_masked_cplx_matrix_hermitic: &
              &MATRIX "//TRIM(ADJUSTL(MM%title))//" ISNT HERMITIC!")
         CALL write_masked_cplx_matrix(MM, SHOW_MASK=.true.)
         call write_array(MM%mat - TRANSPOSE(CONJG(MM%mat)), 'M-M\dag', UNIT=6)
         CALL dump_message(TEXT=" n, m "//c2s(i2c(MM%n1))//" "// &
                           c2s(i2c(MM%n2)))
         write (*, *) 'mat :'
         write (*, *) MM%mat
         write (*, *) 'hermit conj :'
         write (*, *) TRANSPOSE(CONJG(MM%mat))
         write (*, *) 'MATRIX IS NOT HERMITIC'
         STOP 'not hermitic stop'
      ENDIF

   end subroutine

   subroutine write_masked_cplx_matrix(MM, SHOW_MASK, UNIT, SHORT)

      ! WRITE COMPLEX MASKED MATRIX MM(n1, n2)

      use common_def, only: c2s, dump_message, i2c
      use matrix, only: write_array

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(IN) :: MM
      LOGICAL, OPTIONAL, INTENT(IN)             :: SHOW_MASK, SHORT
      INTEGER, OPTIONAL, INTENT(IN)             :: UNIT
      INTEGER              :: i1, i2, unit_
      LOGICAL              :: show_mask_, short_
      CHARACTER(LEN=400) :: fmt_MM

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN write_masked_cplx_matrix: &
           &INPUT ISNT ALLOCATED!"

      show_mask_ = .true. ! DEFAULT: SHOW INTEGER MASK
      IF (PRESENT(SHOW_MASK)) show_mask_ = SHOW_MASK

      ! PRINT INTEGER MASK
      IF (show_mask_) THEN
         IF (MM%MASK%n1 <= MM%MASK%n2) THEN
            CALL write_array(MM%MASK%imat, "# "//TRIM(ADJUSTL(MM%title)), &
                             UNIT=UNIT)
         ELSE
            CALL write_array(TRANSPOSE(MM%MASK%imat), "# "// &
                             TRIM(ADJUSTL(MM%title))//" (transposed)", UNIT=UNIT)
         ENDIF
      ELSE
         CALL dump_message(TEXT="# "//TRIM(ADJUSTL(MM%title)), UNIT=UNIT)
      ENDIF

      ! PRINT MATRIX
      IF (ANY(MM%MASK%mat)) THEN
         IF (show_mask_) CALL dump_message(TEXT="# with", UNIT=UNIT)

         short_ = .false. ! DEFAULT: LONG FORMAT
         IF (PRESENT(SHORT)) short_ = SHORT

         IF (short_) THEN
            WRITE (fmt_MM, *) '(2x, a, 2(a, f0.6), a)'
         ELSE
            WRITE (fmt_MM, *) '(2x, a, 2(a, f0.16), a)'
         ENDIF

         unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
         IF (PRESENT(UNIT)) unit_ = UNIT

         DO i2 = 1, MM%n2
            DO i1 = 1, MM%n1
               IF (MM%MASK%mat(i1, i2)) WRITE (unit_, fmt_MM) '('// &
                  c2s(i2c(i1))//', '//c2s(i2c(i2))//') = '// &
                  c2s(i2c(MM%MASK%imat(i1, i2)))//' = ', '(', &
                  real(MM%mat(i1, i2), kind=DP), ', ', AIMAG(MM%mat(i1, i2)), ')'
            ENDDO
         ENDDO

      ENDIF

      CALL flush(unit_)
   end subroutine

   subroutine write_raw_masked_cplx_matrix(MM, UNIT)

      ! WRITE COMPLEX MASKED MATRIX MM(n1, n2)

      use mask_class, only: write_raw_mask

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(IN) :: MM
      INTEGER, INTENT(IN)                       :: UNIT

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN write_raw_masked_cplx_matrix: &
           &INPUT ISNT ALLOCATED!"
      WRITE (UNIT, *) MM%title
      WRITE (UNIT, *) MM%n1, MM%n2
      WRITE (UNIT, *) MM%IS_HERM
      CALL write_raw_mask(MM%MASK, UNIT)
      IF (MM%IS_HERM) THEN
         CALL write_raw_mask(MM%MASKdiag, UNIT)
         CALL write_raw_mask(MM%MASKoffdiag, UNIT)
      ENDIF
      WRITE (UNIT, *) MM%mat
   end subroutine

   subroutine read_raw_masked_cplx_matrix(MM, UNIT)

      use mask_class, only: read_raw_mask

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MM
      INTEGER, INTENT(IN)                          :: UNIT
      CHARACTER(LEN=100) :: title
      INTEGER              :: n1, n2
      LOGICAL              :: IS_HERM

      CALL delete_masked_cplx_matrix(MM)
      READ (UNIT, *) title
      READ (UNIT, *) n1, n2
      READ (UNIT, *) IS_HERM
      CALL new_masked_cplx_matrix(MM, TRIM(ADJUSTL(title)), n1, n2, IS_HERM= &
                                  IS_HERM)
      CALL read_raw_mask(MM%MASK, UNIT)
      IF (MM%IS_HERM) THEN
         CALL read_raw_mask(MM%MASKdiag, UNIT)
         CALL read_raw_mask(MM%MASKoffdiag, UNIT)
      ENDIF
      READ (UNIT, *) MM%mat
   end subroutine

   subroutine pad_masked_cplx_matrix(MMOUT, MMIN)

      use common_def, only: find_rank

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT)        :: MMOUT
      TYPE(masked_cplx_matrix_type), INTENT(IN), OPTIONAL :: MMIN
      INTEGER :: i1, i2, rankin, rankout

      IF (PRESENT(MMIN)) THEN
         ! FIRST COPY VECTOR OF INDEPENDANT ELEMENTS
         DO rankin = 1, MMIN%MASK%nind
            rankout = find_rank(MMIN%MASK%ivec(rankin), MMOUT%MASK%ivec)
            IF (rankout /= 0) MMOUT%vec(rankout) = MMIN%vec(rankin)
         ENDDO
         ! THEN EXPAND TO FULL MATRIX
         CALL vec2masked_cplx_matrix(MMOUT)
      ENDIF
      DO i1 = 1, MMOUT%n1
         DO i2 = 1, MMOUT%n2
            IF (MMOUT%MASK%mat(i1, i2)) THEN
               CALL fill_masked_cplx_matrix(MMOUT, MMOUT%MASK%imat(i1, i2), &
                                            MMOUT%mat(i1, i2))
            ENDIF
         ENDDO
      ENDDO
   end subroutine

   subroutine filter_masked_cplx_matrix(diag, MM, FILTER)

      use mask_class, only: filter_mask

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: diag
      TYPE(masked_cplx_matrix_type), INTENT(IN)    :: MM
      LOGICAL, INTENT(IN)                          :: FILTER(:, :)

      CALL delete_masked_cplx_matrix(diag)
      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN filter_masked_cplx_matrix: &
           &INPUT ISNT ALLOCATED!"
      IF (SIZE(FILTER, 1) /= MM%n1 .OR. SIZE(FILTER, 2) /= MM%n2) STOP "ERROR IN &
           &filter_masked_cplx_matrix: INCONSISTENT DIMENSIONS!"
      CALL new_masked_cplx_matrix(diag, MM)
      CALL filter_mask(diag%MASK, MM%MASK, FILTER)
      diag%mat = MERGE(MM%mat, CMPLX(0.0_DP, 0.0_DP, 8), FILTER)
   end subroutine

   subroutine slice_masked_cplx_matrix(MMOUT, MMIN, rbounds, cbounds)

      use mask_class, only: slice_mask
      use common_def, only: c2s, i2c

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MMOUT
      TYPE(masked_cplx_matrix_type), INTENT(IN)    :: MMIN
      INTEGER, INTENT(IN)                          :: rbounds(2), cbounds(2)
      CHARACTER(LEN=100) :: cslice
      INTEGER              :: n1slice, n2slice
      LOGICAL              :: IS_HERMslice

      IF (.NOT. ASSOCIATED(MMIN%mat)) STOP "ERROR IN slice_masked_cplx_matrix: &
           &INPUT ISNT ALLOCATED!"
      IF (rbounds(1) < LBOUND(MMIN%mat, 1) .OR. rbounds(2) > UBOUND(MMIN%mat, 1)) &
         STOP "ERROR IN slice_masked_cplx_matrix: INCONSISTENT ROW BOUNDS!"
      IF (cbounds(1) < LBOUND(MMIN%mat, 2) .OR. cbounds(2) > UBOUND(MMIN%mat, 2)) &
           STOP "ERROR IN slice_masked_cplx_matrix: INCONSISTENT COLUMN &
           &BOUNDS!"
      n1slice = rbounds(2) - rbounds(1) + 1
      n2slice = cbounds(2) - cbounds(1) + 1
      IF (n1slice <= 0) STOP "ERROR IN slice_masked_cplx_matrix: NON-NULL ROW &
           &DIMENSIONS REQUIRED!"
      IF (n2slice <= 0) STOP "ERROR IN slice_masked_cplx_matrix: NON-NULL &
           &COLUMN DIMENSIONS REQUIRED!"
      IS_HERMslice = .false.
      IF (ALL(rbounds == cbounds)) IS_HERMslice = MMIN%IS_HERM ! SLICE =
      ! DIAGONAL BLOCK
      cslice = c2s(i2c(rbounds(1)))//"_"//c2s(i2c(rbounds(2)))//"_"// &
               c2s(i2c(cbounds(1)))//"_"//c2s(i2c(cbounds(2)))
      CALL new_masked_cplx_matrix(MMOUT, TRIM(ADJUSTL(MMIN%title))//"_"// &
                                  TRIM(ADJUSTL(cslice)), n1slice, n2slice, IS_HERM=IS_HERMslice)
      CALL slice_mask(MMOUT%MASK, MMIN%MASK, rbounds, cbounds)

      IF (MMOUT%IS_HERM) THEN
         CALL slice_mask(MMOUT%MASKdiag, MMIN%MASKdiag, rbounds, cbounds)
         CALL slice_mask(MMOUT%MASKoffdiag, MMIN%MASKoffdiag, rbounds, cbounds)
      ENDIF

      IF (ASSOCIATED(MMOUT%mat)) DEALLOCATE (MMOUT%mat, STAT=istati)
      IF (ASSOCIATED(MMOUT%vec)) DEALLOCATE (MMOUT%vec, STAT=istati)

      MMOUT%mat => MMIN%mat(rbounds(1):rbounds(2), cbounds(1):cbounds(2))
      MMOUT%vec => MMIN%vec

   end subroutine

   subroutine clean_redundant_imask_cplx(MM)

      use mask_class, only: mask2vec

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MM(:)
      INTEGER :: iMM, jMM, i1, i2, j1, j2, iind
      INTEGER :: nMM, n1, n2

      nMM = SIZE(MM)
      n1 = MM(1)%n1
      n2 = MM(1)%n2
      DO iMM = nMM, 2, -1
         DO i1 = 1, n1
            DO i2 = 1, n2
               IF (MM(iMM)%MASK%mat(i1, i2)) THEN
                  iind = MM(iMM)%MASK%imat(i1, i2)
                  DO jMM = 1, iMM - 1
                     DO j1 = 1, n1
                        DO j2 = 1, n2
                           IF (MM(jMM)%MASK%mat(j1, &
                                                j2) .AND. MM(jMM)%MASK%imat(j1, j2) == iind) THEN
                              MM(iMM)%MASK%mat(i1, i2) = .false.
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CALL mask2vec(MM(iMM)%MASK) ! update vector of independant elements
      ENDDO
   end subroutine

   subroutine transform_masked_cplx_matrix(MMOUT, MMIN, STAT)

      ! MMOUT = (-1)^F * ( TRANSPOSE(MMIN) - Id )
      ! I.E.  ENFORCE (ANTI-)COMMUTATOR
      ! IF    MMIN [ij] = < a[i]*A[j] >
      ! THEN  MMOUT[ij] = < A[i]*a[j] >

      use matrix, only: new_id
      use mask_class, only: mask2vec

      implicit none

      TYPE(masked_cplx_matrix_type), INTENT(INOUT) :: MMOUT
      TYPE(masked_cplx_matrix_type), INTENT(IN)    :: MMIN
      CHARACTER(LEN=*), INTENT(IN)               :: STAT
      REAL(8) :: Id(MMIN%n1, MMIN%n1)

      IF (MMIN%n1 /= MMOUT%n1 .OR. MMIN%n2 /= MMOUT%n2) STOP "ERROR IN &
           &transform_masked_cplx_matrix: INCONSISTENT DIMENSIONS!"
      CALL new_Id(Id)
      SELECT CASE (STAT)
      CASE (FERMIONIC)
         MMOUT%mat = -(TRANSPOSE(MMIN%mat) - Id)
      CASE (BOSONIC)
         MMOUT%mat = (TRANSPOSE(MMIN%mat) - Id)
      END SELECT
      CALL masked_cplx_matrix2vec(MMOUT)
      MMOUT%MASK%imat = TRANSPOSE(MMOUT%MASK%imat)
      CALL mask2vec(MMOUT%MASK)
   end subroutine

   subroutine new_masked_real_matrix_from_scratch(MM, title, n1, n2, IMASK, &
                                                  IS_HERM)

      ! CREATE MASKED MATRIX MM FROM SCRATCH

      use mask_class, only: new_mask, SYM

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MM
      CHARACTER(LEN=*), INTENT(IN)               :: title
      INTEGER, INTENT(IN)                          :: n1, n2
      INTEGER, OPTIONAL, INTENT(IN)                :: IMASK(n1, n2)
      LOGICAL, OPTIONAL, INTENT(IN)                :: IS_HERM
      CHARACTER(LEN=9) :: SYMMETRY

      CALL delete_masked_real_matrix(MM)

      MM%title = title(1:MIN(100, LEN_TRIM(title)))

      MM%IS_HERM = .false.
      IF (PRESENT(IS_HERM)) THEN
         IF (IS_HERM .AND. n1 /= n2) then
            write (*, *) "ERROR in new_masked_real_matrix: "// &
               TRIM(ADJUSTL(title))//" ISNT SQUARE! IS_HERM IRRELEVANT!"
            stop
         endif
         MM%IS_HERM = IS_HERM
      ENDIF

      ! FULL MATRIX
      ! SIZE
      MM%n1 = n1
      MM%n2 = n2
      ! MATRIX
      ALLOCATE (MM%mat(n1, n2))
      MM%mat = 0.0_DP

      ! MASKS
      SYMMETRY = ''
      IF (MM%IS_HERM) SYMMETRY = SYM
      CALL new_mask(MM%MASK, n1, n2, IMASK=IMASK, SYMMETRY=SYMMETRY)

      ! REDUNDANT
      CALL build_mask_real(MM)
   end subroutine

   subroutine new_masked_real_matrix_from_old(MMOUT, MMIN)

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MMOUT
      TYPE(masked_real_matrix_type), INTENT(IN)    :: MMIN

      CALL delete_masked_real_matrix(MMOUT)
      IF (.NOT. ASSOCIATED(MMIN%mat)) STOP "ERROR IN &
           &new_masked_real_matrix_from_old: INPUT ISNT ALLOCATED!"
      CALL new_masked_real_matrix_from_scratch(MMOUT, MMIN%title, MMIN%n1, &
                                               MMIN%n2, IMASK=MMIN%MASK%imat, IS_HERM=MMIN%IS_HERM)
      CALL copy_masked_real_matrix(MMOUT, MMIN)
   end subroutine

   subroutine copy_masked_real_matrix(MMOUT, MMIN)

      ! CREATE MMOUT BY COPYING EXISTING MMIN

      use mask_class, only: copy_mask

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MMOUT
      TYPE(masked_real_matrix_type), INTENT(IN)    :: MMIN

      IF (.NOT. ASSOCIATED(MMIN%mat)) STOP "ERROR IN copy_masked_real_matrix: &
           &INPUT ISNT ALLOCATED!"
      IF (.NOT. ASSOCIATED(MMOUT%mat)) STOP "ERROR IN copy_masked_real_matrix: &
           &OUTPUT ISNT ALLOCATED!"

      MMOUT%IS_HERM = MMIN%IS_HERM

      ! FULL MATRIX

      ! SIZE
      MMOUT%n1 = MMIN%n1
      MMOUT%n2 = MMIN%n2
      ! MATRIX
      IF (.NOT. ASSOCIATED(MMOUT%mat)) ALLOCATE (MMOUT%mat(MMOUT%n1, MMOUT%n2))
      MMOUT%mat = MMIN%mat
      ! MASKS
      CALL copy_mask(MMOUT%MASK, MMIN%MASK)

      ! VECTOR
      IF (ASSOCIATED(MMIN%vec)) THEN
         IF (ASSOCIATED(MMOUT%vec) .AND. SIZE(MMOUT%vec) /= SIZE(MMIN%vec)) &
            DEALLOCATE (MMOUT%vec, STAT=istati)
         IF (.NOT. ASSOCIATED(MMOUT%vec)) ALLOCATE (MMOUT%vec(MMOUT%MASK%nind))
         MMOUT%vec = MMIN%vec
      ENDIF
   end subroutine

   subroutine delete_masked_real_matrix(MM)

      use mask_class, only: delete_mask

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MM

      IF (ASSOCIATED(MM%mat)) DEALLOCATE (MM%mat, STAT=istati)
      IF (ASSOCIATED(MM%vec)) DEALLOCATE (MM%vec, STAT=istati)
      !CEDRIC
      !NULLIFY(MM%mat)
      !NULLIFY(MM%vec)

      CALL delete_mask(MM%MASK)

   end subroutine

   subroutine masked_real_matrix2vec(MM)

      ! PACK MATRIX INTO VECTOR (CREATED IF NECESSARY) ACCORDING TO MASKS

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MM

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN masked_real_matrix2vec: INPUT &
           &ISNT ALLOCATED!"

      ! INTIALIZE VECTORS IF NECESSARY (nind IS SUPPOSED KNOWN)

      IF (.NOT. ASSOCIATED(MM%vec) .AND. MM%MASK%nind /= 0) &
         ALLOCATE (MM%vec(MM%MASK%nind))

      ! PACK MATRIX ELEMENTS

      IF (MM%MASK%nind /= 0) MM%vec = PACK(MM%mat, MM%MASK%mat)
   end subroutine

   subroutine vec2masked_real_matrix(MM, MMEXT)

      ! UNPACK MATRIX INTO VECTOR (ASSUMED TO EXIST) ACCORDING TO MASKS

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT)        :: MM
      TYPE(masked_real_matrix_type), INTENT(IN), OPTIONAL :: MMEXT
      INTEGER :: iind

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN vec2masked_real_matrix: INPUT &
           &ISNT ALLOCATED!"
      IF (PRESENT(MMEXT)) THEN
         DO iind = 1, MMEXT%MASK%nind
            CALL fill_masked_real_matrix(MM, MMEXT%MASK%ivec(iind), &
                                         MMEXT%vec(iind))
         ENDDO
      ELSE
         !MM%mat = 0.0_DP
         DO iind = 1, MM%MASK%nind ! NOTHING HAPPENS IF nind = 0
            CALL fill_masked_real_matrix(MM, MM%MASK%ivec(iind), MM%vec(iind))
         ENDDO
      ENDIF
      IF (MM%IS_HERM) CALL test_masked_real_matrix_symmetric(MM)
   end subroutine

   subroutine fill_masked_real_matrix(MM, iind, val)

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MM
      INTEGER, INTENT(IN)                          :: iind
      REAL(DP) :: val
      INTEGER   :: i1, i2

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN fill_masked_real_matrix: &
           &INPUT ISNT ALLOCATED!"
      DO i1 = 1, MM%n1
         DO i2 = 1, MM%n2
            IF (MM%MASK%imat(i1, i2) == iind) MM%mat(i1, i2) = val
            IF (MM%MASK%imat(i1, i2) == -iind) MM%mat(i1, i2) = -val
         ENDDO
      ENDDO
      ! SYMMETRIZE
      IF (MM%IS_HERM) THEN
         DO i1 = 1, MM%n1
            DO i2 = 1, MM%n2
               IF (MM%MASK%imat(i1, i2) == 0 .AND. MM%MASK%imat(i2, i1) /= 0) &
                  MM%mat(i1, i2) = MM%mat(i2, i1)
            ENDDO
         ENDDO
      ENDIF
   end subroutine

   subroutine build_mask_real(MM, IMASK)

      ! BUILD LOGICAL MASK FROM INTEGER MASK

      use mask_class, only: build_logical_mask

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MM
      INTEGER, OPTIONAL, INTENT(IN)                :: IMASK(:, :)

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN build_mask_real: INPUT ISNT &
           &ALLOCATED!"

      CALL build_logical_mask(MM%MASK, IMASK=IMASK)

   end subroutine

   subroutine test_masked_real_matrix_symmetric(MM)

      ! TEST HERMITICITY OF MATRIX

      use common_def, only: c2s, dump_message, i2c
      use matrix, only: write_array

      implicit none

      TYPE(masked_real_matrix_type), INTENT(IN) :: MM
      LOGICAL :: is_real_symmetric

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN &
           &test_masked_real_matrix_symmetric: INPUT ISNT ALLOCATED!"

      IF (MM%n1 == MM%n2) THEN
         is_real_symmetric = ALL(abs(MM%mat - TRANSPOSE(MM%mat)) < 1.d-8)
      ELSE
         CALL dump_message(TEXT="ERROR IN test_masked_real_matrix_symmetric: &
              &MATRIX "//TRIM(ADJUSTL(MM%title))//" ISNT SQUARE!")
         CALL dump_message(TEXT="n1 = "//c2s(i2c(MM%n1))//" n2 = "// &
                           c2s(i2c(MM%n2)))
      ENDIF

      IF (.NOT. is_real_symmetric) THEN
         CALL dump_message(TEXT="ERROR IN test_masked_real_matrix_symmetric: &
              &MATRIX "//TRIM(ADJUSTL(MM%title))//" ISNT SYMMETRIC!")
         CALL write_masked_real_matrix(MM, SHOW_MASK=.true.)
         call write_array(MM%mat - TRANSPOSE(MM%mat), 'M-M\dag', UNIT=6)
         STOP 'critical error'
      ENDIF

   end subroutine

   subroutine write_masked_real_matrix(MM, SHOW_MASK, UNIT, SHORT)

      ! WRITE COMPLEX MASKED MATRIX MM(n1, n2)

      use common_def, only: c2s, dump_message, i2c
      use matrix, only: write_array

      implicit none

      TYPE(masked_real_matrix_type), INTENT(IN) :: MM
      LOGICAL, OPTIONAL, INTENT(IN)             :: SHOW_MASK, SHORT
      INTEGER, OPTIONAL, INTENT(IN)             :: UNIT
      INTEGER              :: i1, i2, unit_
      LOGICAL              :: show_mask_, short_
      CHARACTER(LEN=400) :: fmt_MM

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN write_masked_real_matrix: &
           &INPUT ISNT ALLOCATED!"

      show_mask_ = .true. ! DEFAULT: SHOW INTEGER MASK
      IF (PRESENT(SHOW_MASK)) show_mask_ = SHOW_MASK
      ! PRINT INTEGER MASK
      IF (show_mask_) THEN
         IF (MM%MASK%n1 <= MM%MASK%n2) THEN
            CALL write_array(MM%MASK%imat, "# "//TRIM(ADJUSTL(MM%title)), &
                             UNIT=UNIT)
         ELSE
            CALL write_array(TRANSPOSE(MM%MASK%imat), "# "// &
                             TRIM(ADJUSTL(MM%title))//" (transposed)", UNIT=UNIT)
         ENDIF
      ELSE
         CALL dump_message(TEXT="# "//TRIM(ADJUSTL(MM%title)), UNIT=UNIT)
      ENDIF

      ! PRINT MATRIX
      IF (ANY(MM%MASK%mat)) THEN
         IF (show_mask_) CALL dump_message(TEXT="# with", UNIT=UNIT)

         short_ = .false. ! DEFAULT: LONG FORMAT
         IF (PRESENT(SHORT)) short_ = SHORT

         IF (short_) THEN
            WRITE (fmt_MM, *) '(2x, a, f0.6)'
         ELSE
            WRITE (fmt_MM, *) '(2x, a, f0.16)'
         ENDIF

         unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
         IF (PRESENT(UNIT)) unit_ = UNIT

         DO i2 = 1, MM%n2
            DO i1 = 1, MM%n1
               IF (MM%MASK%mat(i1, i2)) WRITE (unit_, fmt_MM) '('// &
                  c2s(i2c(i1))//', '//c2s(i2c(i2))//') = '// &
                  c2s(i2c(MM%MASK%imat(i1, i2)))//' = ', MM%mat(i1, i2)
            ENDDO
         ENDDO

      ENDIF
      CALL flush(unit_)
   end subroutine

   subroutine write_raw_masked_real_matrix(MM, UNIT)

      ! WRITE COMPLEX MASKED MATRIX MM(n1, n2)

      use mask_class, only: write_raw_mask

      implicit none

      TYPE(masked_real_matrix_type), INTENT(IN) :: MM
      INTEGER, INTENT(IN)                       :: UNIT

      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN write_raw_masked_real_matrix: &
           &INPUT ISNT ALLOCATED!"
      WRITE (UNIT, *) TRIM(ADJUSTL(MM%title))
      WRITE (UNIT, *) MM%n1, MM%n2
      WRITE (UNIT, *) MM%IS_HERM
      CALL write_raw_mask(MM%MASK, UNIT)
      WRITE (UNIT, *) MM%mat
   end subroutine

   subroutine read_raw_masked_real_matrix(MM, UNIT)

      use mask_class, only: read_raw_mask

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MM
      INTEGER, INTENT(IN)                          :: UNIT
      CHARACTER(LEN=100) :: title
      INTEGER              :: n1, n2
      LOGICAL              :: IS_HERM

      CALL delete_masked_real_matrix(MM)
      READ (UNIT, *) title
      READ (UNIT, *) n1, n2
      READ (UNIT, *) IS_HERM
      CALL new_masked_real_matrix_from_scratch(MM, TRIM(ADJUSTL(title)), n1, &
                                               n2, IS_HERM=IS_HERM)
      CALL read_raw_mask(MM%MASK, UNIT)
      READ (UNIT, *) MM%mat
   end subroutine

   subroutine pad_masked_real_matrix(MMOUT, MMIN)

      use common_def, only: find_rank

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT)        :: MMOUT
      TYPE(masked_real_matrix_type), INTENT(IN), OPTIONAL :: MMIN
      INTEGER :: i1, i2, rankin, rankout

      IF (PRESENT(MMIN)) THEN
         ! FIRST COPY VECTOR OF INDEPENDANT ELEMENTS
         DO rankin = 1, MMIN%MASK%nind
            rankout = find_rank(MMIN%MASK%ivec(rankin), MMOUT%MASK%ivec)
            IF (rankout /= 0) MMOUT%vec(rankout) = MMIN%vec(rankin)
         ENDDO
         ! THEN EXPAND TO FULL MATRIX
         CALL vec2masked_real_matrix(MMOUT)
      ENDIF
      DO i1 = 1, MMOUT%n1
         DO i2 = 1, MMOUT%n2
            IF (MMOUT%MASK%mat(i1, i2)) THEN
               CALL fill_masked_real_matrix(MMOUT, MMOUT%MASK%imat(i1, i2), &
                                            MMOUT%mat(i1, i2))
            ENDIF
         ENDDO
      ENDDO
   end subroutine

   subroutine filter_masked_real_matrix(diag, MM, FILTER)

      use mask_class, only: filter_mask

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: diag
      TYPE(masked_real_matrix_type), INTENT(IN)    :: MM
      LOGICAL, INTENT(IN)                          :: FILTER(:, :)

      CALL delete_masked_real_matrix(diag)
      IF (.NOT. ASSOCIATED(MM%mat)) STOP "ERROR IN filter_masked_real_matrix: &
           &INPUT ISNT ALLOCATED!"
      IF (SIZE(FILTER, 1) /= MM%n1 .OR. SIZE(FILTER, 2) /= MM%n2) STOP "ERROR IN &
           &filter_masked_real_matrix: INCONSISTENT DIMENSIONS!"
      CALL new_masked_real_matrix(diag, MM)
      CALL filter_mask(diag%MASK, MM%MASK, FILTER)
      diag%mat = MERGE(MM%mat, 0.0_DP, FILTER)
   end subroutine

   subroutine slice_masked_real_matrix(MMOUT, MMIN, rbounds, cbounds)

      use mask_class, only: slice_mask
      use common_def, only: c2s, i2c

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MMOUT
      TYPE(masked_real_matrix_type), INTENT(IN)    :: MMIN
      INTEGER, INTENT(IN)                          :: rbounds(2), cbounds(2)
      CHARACTER(LEN=100) :: cslice
      INTEGER              :: n1slice, n2slice
      LOGICAL              :: IS_HERMslice

      CALL delete_masked_real_matrix(MMOUT)
      IF (.NOT. ASSOCIATED(MMIN%mat)) STOP "ERROR IN slice_masked_real_matrix: &
           &INPUT ISNT ALLOCATED!"
      IF (rbounds(1) < LBOUND(MMIN%mat, 1) .OR. rbounds(2) > UBOUND(MMIN%mat, 1)) &
         STOP "ERROR IN slice_masked_real_matrix: INCONSISTENT ROW BOUNDS!"
      IF (cbounds(1) < LBOUND(MMIN%mat, 2) .OR. cbounds(2) > UBOUND(MMIN%mat, 2)) &
           STOP "ERROR IN slice_masked_real_matrix: INCONSISTENT COLUMN &
           &BOUNDS!"
      n1slice = rbounds(2) - rbounds(1) + 1
      n2slice = cbounds(2) - cbounds(1) + 1
      IF (n1slice <= 0) STOP "ERROR IN slice_masked_real_matrix: NON-NULL ROW &
           &DIMENSIONS REQUIRED!"
      IF (n2slice <= 0) STOP "ERROR IN slice_masked_real_matrix: NON-NULL &
           &COLUMN DIMENSIONS REQUIRED!"
      IS_HERMslice = .false.
      IF (ALL(rbounds == cbounds)) IS_HERMslice = MMIN%IS_HERM ! SLICE =
      ! DIAGONAL BLOCK
      cslice = c2s(i2c(rbounds(1)))//"_"//c2s(i2c(rbounds(2)))//"_"// &
               c2s(i2c(cbounds(1)))//"_"//c2s(i2c(cbounds(2)))
      CALL new_masked_real_matrix(MMOUT, TRIM(ADJUSTL(MMIN%title))//"_"// &
                                  TRIM(ADJUSTL(cslice)), n1slice, n2slice, IS_HERM=IS_HERMslice)
      CALL slice_mask(MMOUT%MASK, MMIN%MASK, rbounds, cbounds)

      IF (ASSOCIATED(MMOUT%mat)) DEALLOCATE (MMOUT%mat, STAT=istati)
      IF (ASSOCIATED(MMOUT%vec)) DEALLOCATE (MMOUT%vec, STAT=istati)

      MMOUT%mat => MMIN%mat(rbounds(1):rbounds(2), cbounds(1):cbounds(2))
      MMOUT%vec => MMIN%vec

   end subroutine

   subroutine clean_redundant_imask_real(MM)

      use mask_class, only: mask2vec

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MM(:)
      INTEGER :: iMM, jMM, i1, i2, j1, j2, iind
      INTEGER :: nMM, n1, n2 ! for clarity

      nMM = SIZE(MM)
      n1 = MM(1)%n1
      n2 = MM(1)%n2
      DO iMM = nMM, 2, -1
         DO i1 = 1, n1
            DO i2 = 1, n2
               IF (MM(iMM)%MASK%mat(i1, i2)) THEN
                  iind = MM(iMM)%MASK%imat(i1, i2)
                  DO jMM = 1, iMM - 1
                     DO j1 = 1, n1
                        DO j2 = 1, n2
                           IF (MM(jMM)%MASK%mat(j1, j2) .AND. &
                               MM(jMM)%MASK%imat(j1, j2) == iind) THEN
                              MM(iMM)%MASK%mat(i1, i2) = .false.
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CALL mask2vec(MM(iMM)%MASK) ! update vector of independant elements
      ENDDO
   end subroutine

   subroutine transform_masked_real_matrix(MMOUT, MMIN, STAT)

      ! MMOUT  :  (-1)^F * ( TRANSPOSE(MMIN) - Id ) = (-1)^F * ( MMIN - Id )
      ! I.E.   :  ENFORCE (ANTI-)COMMUTATOR
      ! IF     :  MMIN [ij] = < a[i]*A[j] >
      ! THEN   :  MMOUT[ij] = < A[i]*a[j] >

      use matrix, only: new_id

      implicit none

      TYPE(masked_real_matrix_type), INTENT(INOUT) :: MMOUT
      TYPE(masked_real_matrix_type), INTENT(IN)    :: MMIN
      CHARACTER(LEN=*), INTENT(IN)               :: STAT
      REAL(8) :: Id(MMIN%n1, MMIN%n1)

      IF (MMIN%n1 /= MMOUT%n1 .OR. MMIN%n2 /= MMOUT%n2) STOP "ERROR IN &
           &transform_masked_real_matrix: INCONSISTENT DIMENSIONS!"
      CALL new_Id(Id)
      SELECT CASE (STAT)
      CASE (FERMIONIC)
         MMOUT%mat = -(MMIN%mat - Id)
      CASE (BOSONIC)
         MMOUT%mat = (MMIN%mat - Id)
      END SELECT
      CALL masked_real_matrix2vec(MMOUT)
   end subroutine

end module
