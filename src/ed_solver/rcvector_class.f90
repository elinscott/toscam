MODULE rcvector_class

   use genvar, only: DP

   implicit none

   private

   ! GENERIC REAL OR COMPLEX VECTOR TYPE !

   integer :: istati

   TYPE rcvector_type
      INTEGER :: n = 0
#ifdef _complex
      COMPLEX(DP), POINTER :: rc(:) => NULL()
#else
      REAL(DP), POINTER :: rc(:) => NULL()
#endif
   END TYPE

   TYPE rcvector_archive_type
      INTEGER                       :: nvec = 0 ! number of archived vectors
      TYPE(rcvector_type), POINTER  :: vec(:) => NULL() ! pile of archived vectors
   END TYPE

   INTERFACE new_rcvector
      MODULE PROCEDURE new_rcvector_from_scratch
      MODULE PROCEDURE new_rcvector_from_old
   END INTERFACE

   public :: create_fix_initial_vector
   public :: delete_rcvector
   public :: new_rcvector
   public :: norm_rcvector
   public :: rcvector_type
   public :: read_raw_rcvector
   public :: write_raw_rcvector

contains

   subroutine create_fix_initial_vector(initvec, n)

      use random, only: dran_tab

      implicit none

      TYPE(rcvector_type), INTENT(INOUT) :: initvec
      INTEGER, INTENT(IN)                :: n
      INTEGER :: i

      CALL new_rcvector(initvec, n)

#ifdef _complex
      do i = 1, size(initvec%rc)
         initvec%rc(i) = CMPLX(dran_tab(i), dran_tab(i + 20), 8)
      enddo
#else
      do i = 1, size(initvec%rc)
         initvec%rc(i) = dran_tab(i)
      enddo
#endif

   end subroutine

   subroutine new_rcvector_from_scratch(VEC, N)

      implicit none

      TYPE(rcvector_type), INTENT(INOUT) :: VEC
      INTEGER, INTENT(IN)                :: N

      CALL delete_rcvector(VEC)
      VEC%n = N
      IF (N > 0) THEN
         ALLOCATE (VEC%rc(N))
      ENDIF
      VEC%rc = 0.0_DP

   end subroutine

   subroutine new_rcvector_from_old(VECOUT, VECIN)

      implicit none

      TYPE(rcvector_type), INTENT(INOUT) :: VECOUT
      TYPE(rcvector_type), INTENT(IN)    :: VECIN

      CALL delete_rcvector(VECOUT)
      IF (VECIN%n == 0) STOP "ERROR IN new_rcvector_from_old: INPUT ISNT &
           &ASSOCIATED!"
      CALL new_rcvector_from_scratch(VECOUT, VECIN%n)
      CALL copy_rcvector(VECOUT, VECIN)
   end subroutine

   subroutine copy_rcvector(VECOUT, VECIN)

      use common_def, only: create_seg_fault

      implicit none

      TYPE(rcvector_type), INTENT(INOUT) :: VECOUT
      TYPE(rcvector_type), INTENT(IN)    :: VECIN

      IF (SIZE(VECOUT%rc) /= SIZE(VECIN%rc)) then
         write (*, *) "ERROR IN copy_vec: INCONSISTENT DIMENSIONS!"
         call create_seg_fault
         STOP "ERROR IN copy_vec: INCONSISTENT DIMENSIONS!"
      endif
      VECOUT%n = VECIN%n
      VECOUT%rc = VECIN%rc
   end subroutine

   subroutine delete_rcvector(VECIN)

      implicit none

      type(rcvector_type), intent(inout) :: VECIN

      IF (ASSOCIATED(VECIN%rc)) DEALLOCATE (VECIN%rc, STAT=istati)
   end subroutine

   function conj_rcvector(vec) RESULT(cvec)

      implicit none

      TYPE(rcvector_type), INTENT(IN) :: vec
      TYPE(rcvector_type) :: cvec

      CALL new_rcvector_from_scratch(cvec, vec%n)
#ifdef _complex
      cvec%rc = CONJG(vec%rc)
#else
      cvec%rc = vec%rc
#endif
   end function

   function norm_rcvector(vec)

      use globalvar_ed_solver, only: USE_TRANSPOSE_TRICK_MPI
      use mpi_mod, only: mpi_dot_product

      implicit none

      TYPE(rcvector_type), INTENT(IN) :: vec
      REAL(DP) :: norm_rcvector

      norm_rcvector = SQRT(ABS(MPI_DOT_PRODUCT(vec%rc, vec%rc, split= &
                                               USE_TRANSPOSE_TRICK_MPI)))
   end function

   function new_rcvector_archive(n, nvec) RESULT(archive)

      implicit none

      INTEGER, INTENT(IN) :: nvec, n
      TYPE(rcvector_archive_type) :: archive
      INTEGER                     :: ivec

      archive%nvec = nvec
      ALLOCATE (archive%vec(nvec))
      DO ivec = 1, nvec
         CALL new_rcvector_from_scratch(archive%vec(ivec), n)
      ENDDO
   end function

   subroutine delete_rcvector_archive(archive)

      implicit none

      type(rcvector_archive_type), intent(inout) :: archive
      integer :: ivec

      IF (ASSOCIATED(archive%vec)) THEN
         DO ivec = 1, archive%nvec
            CALL delete_rcvector(archive%vec(ivec))
         ENDDO
         DEALLOCATE (archive%vec, STAT=istati)
      ENDIF
   end subroutine

   subroutine write_raw_rcvector(vec, UNIT, FILEOUT)

      use common_def, only: close_safe, open_safe

      implicit none

      TYPE(rcvector_type), INTENT(IN)          :: vec
      INTEGER, INTENT(IN), OPTIONAL            :: UNIT
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: FILEOUT
      INTEGER :: unit_

      unit_ = 82548
      IF (PRESENT(UNIT)) unit_ = UNIT
      IF (PRESENT(FILEOUT)) CALL open_safe(unit_, FILEOUT, 'UNKNOWN', 'WRITE')
      WRITE (unit_, *) vec%n
      WRITE (unit_, *) vec%rc
      CALL flush(unit_)
      IF (PRESENT(FILEOUT)) CALL close_safe(unit_)
   end subroutine

   subroutine read_raw_rcvector(vec, UNIT)

      implicit none

      TYPE(rcvector_type), INTENT(INOUT) :: vec
      INTEGER, INTENT(IN)                :: UNIT
      INTEGER :: n

      CALL delete_rcvector(vec)
      READ (UNIT, *) n
      CALL new_rcvector(vec, n)
      READ (UNIT, *) vec%rc
   end subroutine

end module
