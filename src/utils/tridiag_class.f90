MODULE tridiag_class

   USE common_def
   use genvar, only: DBL, strongstop
   use matrix, only: eigenvector_matrix

   IMPLICIT NONE

   private
   public :: delete_tridiag
   public :: diagonalize_tridiag
   public :: invert_zmtridiag
   public :: new_tridiag
   public :: new_tridiag_from_old
   public :: new_tridiag_from_scratch
   public :: submatrix_tridiag
   public :: tridiag_type

   REAL(DBL), PARAMETER, PRIVATE                 :: zero = 0.0_DBL, one = 1.0_DBL, two = 2.0_DBL, three = 3.0_DBL, four = 4.0_DBL
   LOGICAL, PARAMETER, PRIVATE                 :: F = .FALSE., T = .TRUE.

   ! ! GENERATE TRIDIAGONAL LANCZOS MATRIX RECURSIVELY

   TYPE tridiag_type ! TRI-DIAGONAL LANCZOS MATRIX
      INTEGER            ::          N = 0       ! matrix size
      REAL(DBL), POINTER ::    diag(:) => NULL() ! diagonal    (1..N)
      REAL(DBL), POINTER :: subdiag(:) => NULL() ! subdiagonal (1..N)
   END TYPE

   INTERFACE new_tridiag
      MODULE PROCEDURE new_tridiag_from_scratch
      MODULE PROCEDURE new_tridiag_from_old
   END INTERFACE

CONTAINS

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   SUBROUTINE new_tridiag_from_scratch(tri, N)
      TYPE(tridiag_type), INTENT(INOUT) :: tri
      INTEGER, INTENT(IN)    :: N
      INTEGER                           :: istati

      CALL delete_tridiag(tri)
      tri%N = N
      DEALLOCATE (tri%diag, STAT=istati); DEALLOCATE (tri%subdiag, STAT=istati)
      ALLOCATE (tri%diag(1:N)); ALLOCATE (tri%subdiag(1:N))

      tri%diag = zero
      tri%subdiag = zero

      ! REPRESENTATION = DAGOTTO (RMP 1994)
      !
      !    a(1) b(2) 0    0    0
      !    b(2) a(2) b(3) 0    0
      !    0    b(3) a(3) b(4) 0
      !    0    0    b(4) a(5) b(6)
      !    0    0    0    b(6) a(7)
      !
      ! b(1) is dummy (eases Lanczos recursion)

   END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   SUBROUTINE new_tridiag_from_old(TRIOUT, TRIIN)
      TYPE(tridiag_type), INTENT(INOUT) :: TRIOUT
      TYPE(tridiag_type), INTENT(IN)    :: TRIIN

      IF (TRIIN%N == 0) STOP "ERROR IN new_tridiag_from_old: INPUT ISN'T ALLOCATED!!"
      CALL delete_tridiag(TRIOUT)
      CALL new_tridiag_from_scratch(TRIOUT, TRIIN%N)
      CALL copy_tridiag(TRIOUT, TRIIN)

   END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   SUBROUTINE copy_tridiag(TRIOUT, TRIIN)
      TYPE(tridiag_type), INTENT(INOUT) :: TRIOUT
      TYPE(tridiag_type), INTENT(IN)    :: TRIIN
      TRIOUT%N = TRIIN%N
      TRIOUT%diag = TRIIN%diag
      TRIOUT%subdiag = TRIIN%subdiag
   END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   SUBROUTINE delete_tridiag(TRI)
      TYPE(tridiag_type), INTENT(INOUT) :: TRI
      IF (ASSOCIATED(TRI%diag)) DEALLOCATE (TRI%diag)
      IF (ASSOCIATED(TRI%subdiag)) DEALLOCATE (TRI%subdiag)
   END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   SUBROUTINE submatrix_tridiag(TRIOUT, TRIIN, bnds)

      ! THIS TAKES THE SQUARE SUBMATRIX WITH BOUNDS bnds(1)..bnds(2)
      ! NO ARRAY IS ALLOCATED HERE, THIS MERELY POINTS TO THE SUBMATRIX!

      TYPE(tridiag_type), INTENT(INOUT) :: TRIOUT
      TYPE(tridiag_type), INTENT(IN)    :: TRIIN
      INTEGER, INTENT(IN)    :: bnds(2)
      INTEGER                           :: N

   IF (bnds(1) < LBOUND(TRIIN%diag, 1) .OR. bnds(2) > UBOUND(TRIIN%diag, 1)) STOP "ERROR IN submatrix_tridiag: INCONSISTENT BOUNDS!"

      N = bnds(2) - bnds(1) + 1

      IF (strongstop .and. N <= 1) STOP "ERROR IN submatrix_tridiag: DIMENSION > 1 REQUIRED!"

      CALL new_tridiag(TRIOUT, N)

      TRIOUT%diag = TRIIN%diag(bnds(1):bnds(2))
      TRIOUT%subdiag = TRIIN%subdiag(bnds(1):bnds(2))

   END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   FUNCTION invert_zmtridiag(z, tri) RESULT(zmtrim1)
      COMPLEX(DBL)                   :: zmtrim1
      COMPLEX(DBL), INTENT(IN) :: z
      TYPE(tridiag_type), INTENT(IN) :: tri
      INTEGER                        :: iter
      COMPLEX(DBL)                   :: det_ratio_i, det_ratio_ip1, aa

      !-------------------------------------------------------------------------------------------------!
      ! COMPUTE <0| 1 / ( z - T ) |0> WHERE T IS A TRIDIAGONAL MATRIX AND |0> IS THE FIRST BASIS VECTOR !
      !-------------------------------------------------------------------------------------------------!

      aa = (z - tri%diag(tri%N))

      if (abs(aa) < epsilon(1.d0)) then
         write (*, *) 'error NaN division in zmtridiag, 0 denominator'
         write (*, *) aa
         stop 'critical'
      endif

      det_ratio_ip1 = 1.d0/aa

      DO iter = tri%N - 1, 1, -1
         aa = (z - tri%diag(iter) - det_ratio_ip1*tri%subdiag(iter + 1)**2)
         if (abs(aa) < epsilon(1.d0)) then
            write (*, *) 'error NaN division in zmtridiag, 0 denominator_ '
            write (*, *) aa
            stop 'critical'
         endif
         det_ratio_i = 1.d0/aa
         det_ratio_ip1 = det_ratio_i
      END DO

      if (tri%N == 1) det_ratio_i = det_ratio_ip1
      zmtrim1 = det_ratio_i

   END FUNCTION

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   SUBROUTINE diagonalize_tridiag(Lmatrix, VALP, VECP, EIGENVAL_ONLY)

      ! CONVENIENT WRAPPER OF DSTEDC (DIVIDE&CONQUER/LAPACK) TO DIAGONALIZE THE UPPER-LEFT BLOCK OF SIZE N <= Lmatrix%N
      ! OF THE LANCZOS MATRIX

      TYPE(tridiag_type), INTENT(INOUT) :: Lmatrix
      REAL(DBL), INTENT(INOUT) :: VALP(:), VECP(:, :)
      LOGICAL, OPTIONAL, INTENT(IN)    :: EIGENVAL_ONLY
      INTEGER, ALLOCATABLE            :: IWORK(:)
      REAL(DBL), ALLOCATABLE            :: WORK(:)
      CHARACTER(LEN=1)                  :: FLAG
      INTEGER                           :: N, info
      REAL(DBL)                         :: mat(2, 2)

      N = SIZE(VALP)

      IF (N /= Lmatrix%N .OR. N /= SIZE(VECP, 1) .OR. N /= SIZE(VECP, 2)) then
         write (*, *) "ERROR IN diagonalize: INCONSISTENT DIMENSIONS!"
         write (*, *) 'N,Lmatrix%N : ', N, Lmatrix%N
         write (*, *) 'shape(vecp) : ', shape(VECP)
         STOP
      endif
      FLAG = 'I'  ! COMPUTE EIGENVAL AND EIGENVEC BY DEFAULT
      IF (PRESENT(EIGENVAL_ONLY)) THEN
         IF (EIGENVAL_ONLY) FLAG = 'N'
      ENDIF

      IF (ALLOCATED(WORK)) DEALLOCATE (WORK); IF (ALLOCATED(IWORK)) DEALLOCATE (IWORK)
      ALLOCATE (WORK(1 + 4*N + N**2), IWORK(3 + 5*N)); WORK = 0; IWORK = 0; info = 0

      VALP = 0.d0; VECP = 0.d0

      if (N > 2) then
         CALL DSTEDC(FLAG, N, Lmatrix%diag, Lmatrix%subdiag(2:), VECP, N, WORK, SIZE(WORK), IWORK, SIZE(IWORK), info)
      else
         call small_mat
      endif

      ! REMEMBER STORAGE PATTERN OF TRIDIAG MATRIX: subdiab(1) IS NOT REFERENCED!

      VALP = Lmatrix%diag(1:N)

      IF (info > 0) THEN
         CALL dump_message(TEXT="FAILED TO COMPUTE AN EIGENVALUE [DSTEDC]")
         CALL dump_message(TEXT="WHILE WORKING ON SUBMATRIX "//c2s(i2c(info/(N + 1)))//" x "//c2s(i2c(MOD(info, N + 1))))
      ENDIF

      IF (info < 0) THEN
         CALL dump_message(TEXT=c2s(i2c(info))//"-TH ARGUMENT HAD AN ILLEGAL VALUE [DSTEDC]")
      ENDIF

      DEALLOCATE (WORK); DEALLOCATE (IWORK)

   contains

      !------------------------!
      subroutine small_mat
      if (N == 2) then
         mat(1, 1) = Lmatrix%diag(1)
         mat(2, 2) = Lmatrix%diag(2)
         mat(1, 2) = Lmatrix%subdiag(2)
         mat(2, 1) = mat(1, 2)
         CALL eigenvector_matrix(N, mat=mat, vaps=VALP, eigenvec=VECP)
      elseif (N == 1) then
         VALP = 1.; VECP = 1.
      endif
      end subroutine

      !------------------------!

   END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

END MODULE
