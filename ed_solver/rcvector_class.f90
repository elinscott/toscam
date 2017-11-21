MODULE rcvector_class

  USE masked_matrix_class 
  use random, only : dran_tab

  IMPLICIT NONE

  REAL(DBL),    PARAMETER, PRIVATE  :: zero=0.0_DBL
  LOGICAL,      PARAMETER, PRIVATE  :: F=.FALSE.,T=.TRUE.

! GENERIC REAL OR COMPLEX VECTOR TYPE !

  TYPE rcvector_type
    INTEGER :: n = 0
#ifdef _complex
    COMPLEX(DBL), POINTER :: rc(:) => NULL()
#else
    REAL(DBL),    POINTER :: rc(:) => NULL()
#endif
  END TYPE

  TYPE rcvector_archive_type 
     INTEGER                       ::  nvec   = 0        ! number of archived vectors
     TYPE(rcvector_type), POINTER  ::  vec(:) => NULL()  ! pile of archived vectors
  END TYPE

  INTERFACE new_rcvector
    MODULE PROCEDURE new_rcvector_from_scratch
    MODULE PROCEDURE new_rcvector_from_old
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

 SUBROUTINE create_fix_initial_vector(initvec,n)
 TYPE(rcvector_type), INTENT(INOUT) :: initvec
 INTEGER,             INTENT(IN)    :: n
 INTEGER                            :: i

 CALL new_rcvector(initvec,n)

#ifdef _complex
 do i=1,size(initvec%rc)
  initvec%rc(i)=CMPLX(dran_tab(i),dran_tab(i+20),8)
 enddo
#else
 do i=1,size(initvec%rc)
  initvec%rc(i)=dran_tab(i)
 enddo
#endif

 END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_rcvector_from_scratch(VEC,N)
    TYPE(rcvector_type), INTENT(INOUT) :: VEC
    INTEGER,             INTENT(IN)    :: N
        CALL delete_rcvector(VEC)
        VEC%n = N
        IF(N>0)THEN ; ALLOCATE(VEC%rc(N)) ; ENDIF
        VEC%rc = zero
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_rcvector_from_old(VECOUT,VECIN)
    TYPE(rcvector_type), INTENT(INOUT) :: VECOUT
    TYPE(rcvector_type), INTENT(IN)    :: VECIN
    CALL delete_rcvector(VECOUT)
    IF(VECIN%n==0) STOP "ERROR IN new_rcvector_from_old: INPUT  ISNT ASSOCIATED!"
    CALL new_rcvector_from_scratch(VECOUT,VECIN%n) 
    CALL copy_rcvector(VECOUT,VECIN)
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE copy_rcvector(VECOUT,VECIN)
    TYPE(rcvector_type), INTENT(INOUT) :: VECOUT
    TYPE(rcvector_type), INTENT(IN)    :: VECIN
    IF(SIZE(VECOUT%rc)/=SIZE(VECIN%rc)) then
      write(*,*) "ERROR IN copy_vec: INCONSISTENT DIMENSIONS!"
      call create_seg_fault
      STOP "ERROR IN copy_vec: INCONSISTENT DIMENSIONS!"
    endif
    VECOUT%n  = VECIN%n
    VECOUT%rc = VECIN%rc
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_rcvector(VECIN)
    TYPE(rcvector_type), INTENT(INOUT) :: VECIN
    IF(ASSOCIATED(VECIN%rc)) DEALLOCATE(VECIN%rc,STAT=istati)
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION conj_rcvector(vec) RESULT(cvec)
    TYPE(rcvector_type)              :: cvec
    TYPE(rcvector_type), INTENT(IN)  ::  vec
    CALL new_rcvector_from_scratch(cvec,vec%n)
#ifdef _complex
    cvec%rc = CONJG(vec%rc)
#else
    cvec%rc =       vec%rc 
#endif
  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION norm_rcvector(vec)
    REAL(DBL)                       :: norm_rcvector
    TYPE(rcvector_type), INTENT(IN) :: vec
    norm_rcvector = SQRT(ABS(MPI_DOT_PRODUCT(vec%rc,vec%rc,split=USE_TRANSPOSE_TRICK_MPI)))
  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION new_rcvector_archive(n,nvec) RESULT(archive)
    TYPE(rcvector_archive_type) :: archive
    INTEGER,         INTENT(IN) :: nvec,n
    INTEGER                     :: ivec
    archive%nvec = nvec
    ALLOCATE(archive%vec(nvec))
    DO ivec=1,nvec
     CALL new_rcvector_from_scratch(archive%vec(ivec),n)
    ENDDO
  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_rcvector_archive(archive)
    TYPE(rcvector_archive_type), INTENT(INOUT) :: archive
    INTEGER                                    :: ivec
    IF(ASSOCIATED(archive%vec))THEN
      DO ivec=1,archive%nvec
       CALL delete_rcvector(archive%vec(ivec))
      ENDDO
      DEALLOCATE(archive%vec,STAT=istati)
    ENDIF
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE write_raw_rcvector(vec,UNIT,FILEOUT) 
    TYPE(rcvector_type), INTENT(IN)           :: vec
    INTEGER,             INTENT(IN), OPTIONAL :: UNIT
    CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: FILEOUT
    INTEGER                                   :: unit_
                      unit_ = 82548
    IF(PRESENT(UNIT)) unit_ = UNIT
    IF(PRESENT(FILEOUT)) CALL open_safe(unit_,FILEOUT,'UNKNOWN','WRITE')
    WRITE(unit_,*) vec%n
    WRITE(unit_,*) vec%rc
    CALL flush(unit_)
    IF(PRESENT(FILEOUT))  CALL close_safe(unit_)
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE read_raw_rcvector(vec,UNIT) 
    TYPE(rcvector_type), INTENT(INOUT) :: vec
    INTEGER,             INTENT(IN)    :: UNIT
    INTEGER                            :: n 
    CALL delete_rcvector(vec)
    READ(UNIT,*) n
    CALL new_rcvector(vec,n)
    READ(UNIT,*) vec%rc
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
