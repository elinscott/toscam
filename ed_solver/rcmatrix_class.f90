MODULE rcmatrix_class

   use genvar, only: DBL

   IMPLICIT NONE

   private

   !----------------------------------------!
   ! GENERIC REAL OR COMPLEX MATRIX TYPE    !
   !----------------------------------------!

   TYPE rcmatrix_type
      INTEGER :: n1 = 0, n2 = 0
#ifdef _complex
      COMPLEX(DBL), POINTER :: rc(:, :) => NULL()
#else
      REAL(DBL),    POINTER :: rc(:, :) => NULL()
#endif
   END TYPE


   TYPE rcmatrix_archive_type ! archive of matrix
   INTEGER                      :: nmat   = 0       ! number of archived matrix
   TYPE(rcmatrix_type), POINTER :: mat(:) => NULL() ! pile of archived matrix
   END TYPE

   INTERFACE new_rcmatrix
      MODULE PROCEDURE new_rcmatrix_from_scratch
      MODULE PROCEDURE new_rcmatrix_from_old
   END INTERFACE

   public :: rcmatrix_type
   public :: new_rcmatrix
   public :: write_raw_rcmatrix

contains

   subroutine new_rcmatrix_from_scratch(MAT, N1, N2)

      implicit none

      TYPE(rcmatrix_type), INTENT(INOUT) :: MAT
      INTEGER, INTENT(IN)                :: N1
      INTEGER, OPTIONAL, INTENT(IN)      :: N2

      MAT%n1 = N1
      IF(PRESENT(N2))THEN
         MAT%n2 = N2
      ELSE
         MAT%n2 = N1
      ENDIF
      IF(MAT%n1 /= 0 .AND. MAT%n2 /= 0)THEN
         ALLOCATE(MAT%rc(MAT%n1, MAT%n2))
         MAT%rc = 0.0_DBL
      ENDIF
   end subroutine

   subroutine new_rcmatrix_from_old(MATOUT, MATIN)

      implicit none

      TYPE(rcmatrix_type), INTENT(INOUT) :: MATOUT
      TYPE(rcmatrix_type), INTENT(IN)    :: MATIN

      IF(MATIN%n1*MATIN%n2 == 0) STOP "ERROR IN new_rcmatrix_from_old: INPUT &
           &ISNT ASSOCIATED!"
      CALL new_rcmatrix_from_scratch(MATOUT, MATIN%n1, MATIN%n2)
      CALL copy_rcmatrix(MATOUT, MATIN)

   end subroutine

   subroutine copy_rcmatrix(MATOUT, MATIN)

      implicit none

      TYPE(rcmatrix_type), INTENT(INOUT) :: MATOUT
      TYPE(rcmatrix_type), INTENT(IN)    :: MATIN

      IF(ANY(SHAPE(MATOUT%rc) /= SHAPE(MATIN%rc))) STOP "ERROR IN copy_mat: &
           &INCONSISTENT DIMENSIONS!"
      MATOUT%n1 = MATIN%n1
      MATOUT%n2 = MATIN%n2
      MATOUT%rc = MATIN%rc
   end subroutine

   subroutine delete_rcmatrix(MAT)

      implicit none

      TYPE(rcmatrix_type), INTENT(INOUT) :: MAT

      IF(ASSOCIATED(MAT%rc)) DEALLOCATE(MAT%rc)
   end subroutine

   subroutine conj_rcmatrix(cmat, mat)

      implicit none

      TYPE(rcmatrix_type), INTENT(INOUT) :: cmat
      TYPE(rcmatrix_type), INTENT(IN)    :: mat

      CALL new_rcmatrix_from_scratch(cmat, mat%n1, mat%n2)
#ifdef _complex
      cmat%rc = CONJG(mat%rc)
#else
      cmat%rc =       mat%rc
#endif
   end subroutine

   subroutine herm_conj_rcmatrix(hcmat, mat)

      implicit none

      TYPE(rcmatrix_type), INTENT(INOUT) :: hcmat
      TYPE(rcmatrix_type), INTENT(IN)    :: mat

      CALL new_rcmatrix_from_scratch(hcmat, mat%n1, mat%n2)
#ifdef _complex
      hcmat%rc = TRANSPOSE(CONJG(mat%rc))
#else
      hcmat%rc = TRANSPOSE(      mat%rc)
#endif
   end subroutine

   function new_rcmatrix_archive(n1, n2, nmat) RESULT(archive)

      implicit none

      INTEGER, INTENT(IN) :: nmat, n1, n2
      TYPE(rcmatrix_archive_type) :: archive
      INTEGER                     :: imat

      archive%nmat = nmat
      IF(nmat > 0)THEN
         ALLOCATE(archive%mat(nmat))
         DO imat = 1, nmat
            CALL new_rcmatrix_from_scratch(archive%mat(imat), n1, n2)
         ENDDO
      ENDIF
   end function

   subroutine delete_rcmatrix_archive(archive)

      implicit none

      TYPE(rcmatrix_archive_type), INTENT(INOUT) :: archive
      INTEGER :: imat

      IF(ASSOCIATED(archive%mat))THEN
         DO imat = 1, archive%nmat
            CALL delete_rcmatrix(archive%mat(imat))
         ENDDO
         DEALLOCATE(archive%mat)
      ENDIF
   end subroutine

   subroutine write_raw_rcmatrix(mat, UNIT, FILEOUT)

      use common_def, only: close_safe, open_safe

      implicit none

      TYPE(rcmatrix_type), INTENT(IN)          :: mat
      INTEGER, INTENT(IN), OPTIONAL            :: UNIT
      CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: FILEOUT
      INTEGER :: unit_

      unit_ = 82547
      IF(PRESENT(UNIT)) unit_ = UNIT

      IF(PRESENT(FILEOUT)) CALL open_safe(unit_, FILEOUT, 'UNKNOWN', 'WRITE')

      WRITE(unit_, *) mat%n1
      WRITE(unit_, *) mat%n2
      WRITE(unit_, *) mat%rc
      CALL flush(unit_)

      IF(PRESENT(FILEOUT)) CALL close_safe(unit_)

   end subroutine

   subroutine read_raw_rcmatrix(mat, UNIT, FILEIN)

      use common_def, only: close_safe, open_safe

      implicit none

      TYPE(rcmatrix_type), INTENT(INOUT)       :: mat
      INTEGER, INTENT(IN), OPTIONAL            :: UNIT
      CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: FILEIN
      INTEGER :: n1, n2, unit_

      unit_ = 82548
      IF(PRESENT(UNIT)) unit_ = UNIT

      IF(PRESENT(FILEIN)) CALL open_safe(unit_, FILEIN, 'UNKNOWN', 'READ')

      READ(unit_, *) n1
      READ(unit_, *) n2
      CALL new_rcmatrix(mat, n1, n2)
      READ(unit_, *) mat%rc

      IF(PRESENT(FILEIN)) CALL close_safe(unit_)

   end subroutine

end module
