MODULE mask_class

  USE common_def, only: i2c, c2s, dump_message
  use genvar,     only: DBL
  use matrix,     only: write_array
  use linalg,     only:

  implicit none 

  private

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !$$ INTEGER/LOGICAL MASK CLASS $$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  CHARACTER(LEN=9), PARAMETER :: SYM  ='SYMMETRIC'
  CHARACTER(LEN=9), PARAMETER :: HERM ='HERMITIC '
  CHARACTER(LEN=9), PARAMETER :: EMPTY='NOSYM    '
  INTEGER,PRIVATE             :: istati 
  
  TYPE mask_type
    !----------------!
    ! FULL 2D MATRIX !
    !----------------!
    CHARACTER(LEN=9) :: SYMMETRY = '\0'
    ! SIZE OF MATRIX
    INTEGER          :: n1=0, n2=0
    ! MASKS
    INTEGER, dimension(:,:), pointer :: imat => NULL()
    LOGICAL, dimension(:,:), pointer :: mat  => NULL()
    !---------------------------------------!
    ! VECTOR OF INDEPENDANT MATRIX ELEMENTS !
    !---------------------------------------!
    ! NUMBER OF INDEPDT ELEMTS
    INTEGER          :: nind=0    
    ! INTEGER MASK
    INTEGER, dimension(:), pointer :: ivec => NULL()
  END TYPE 

  INTERFACE new_mask
    MODULE PROCEDURE new_mask_from_scratch
    MODULE PROCEDURE new_mask_from_old
  END INTERFACE

  public :: build_logical_mask
  public :: HERM
  public :: mask_type
  public :: new_mask
  public :: SYM

CONTAINS


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE order_matrix_elements(tabelem,MASK)
   INTEGER,         INTENT(INOUT) :: tabelem(:,:)
   TYPE(mask_type), INTENT(IN)    :: MASK
   INTEGER                        :: Nc,iind,site1,site2

    IF(SIZE(tabelem,1)/=COUNT(MASK%mat)) then
      write(*,*) "ERROR IN order_matrix_elements: INCONSISTENT DIMENSIONS!"
      STOP
    ENDIF
  !-------------------------------------------------------------------------------!
  ! REORDER MATRIX ELEMENTS FOR PRETTY PRINT: DIAGONAL FIRST + GATHER (i,j) (j,i) !
  !-------------------------------------------------------------------------------!

    Nc      = MASK%n1
    tabelem = 0
    iind    = 0

    DO site1=1,Nc; IF(MASK%mat(site1,site1))THEN
      iind = iind + 1
      tabelem(iind,:) = site1
    ENDIF; ENDDO
    DO site1=1,Nc; DO site2=site1+1,Nc
      IF(MASK%mat(site1,site2))THEN
       iind = iind + 1
       tabelem(iind,:) = (/site1,site2/)
      ENDIF
      IF(MASK%mat(site2,site1))THEN
       iind = iind + 1
       tabelem(iind,:) = (/site2,site1/)
      ENDIF
    ENDDO; ENDDO

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_mask_from_scratch(MASK,n1,n2,IMASK,SYMMETRY) 

    !--------------------------!
    ! CREATE MASK FROM SCRATCH !
    !--------------------------!

    TYPE(mask_type),            INTENT(INOUT) :: MASK
    INTEGER,                    INTENT(IN)    :: n1,n2
    INTEGER,          OPTIONAL, INTENT(IN)    :: IMASK(n1,n2)
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: SYMMETRY
    INTEGER                                   :: iind,i1,i2
    LOGICAL                                   :: is_diag(n1,n1)
    
    CALL delete_mask(MASK)
    
    IF(n1/=n2.AND.(SYMMETRY==HERM.OR.SYMMETRY==SYM)) then
     write(*,*) "ERROR IN new_mask: SYMMETRY = "//TRIM(SYMMETRY)//" SPECIFIED BUT "//c2s(i2c(n1))//"x"//c2s(i2c(n2))//" MATRIX!"
     stop 'termine'
    endif

    MASK%SYMMETRY = EMPTY
    IF(PRESENT(SYMMETRY))THEN
      IF(LEN_TRIM(SYMMETRY)/=0) MASK%SYMMETRY = TRIM(ADJUSTL(SYMMETRY))
    ENDIF
   
    ! FULL MATRIX
    IF(n1*n2==0) then
     write(*,*) 'N1,N2 : ', n1, n2
     write(*,*) "ERROR IN new_mask_from_scratch: WRONG DIMENSIONS!"
     STOP
    endif

    MASK%n1 = n1 ; MASK%n2 = n2

    ALLOCATE(MASK%imat(n1,n2)); ALLOCATE( MASK%mat(n1,n2))

    IF(PRESENT(IMASK))THEN
      MASK%imat = IMASK
      IF(MASK%SYMMETRY/=EMPTY)THEN
       DO i1=1,MASK%n1; DO i2=i1+1,MASK%n2
        IF(MASK%imat(i1,i2)*MASK%imat(i2,i1)/=0) MASK%imat(i2,i1) = 0
       ENDDO; ENDDO
      ENDIF
    ELSE
      CALL ramp(MASK%imat)
      IF(MASK%SYMMETRY/=EMPTY)THEN
      DO i2=1,MASK%n2 
       MASK%imat(i2,i2) = MASK%imat(1,1)
       DO i1=i2+1,MASK%n1
         MASK%imat(i1,i2) = 0
       ENDDO
      ENDDO
      ENDIF
    ENDIF

    CALL build_logical_mask(MASK) 

    CALL new_diag(is_diag,n1)

     IF(MASK%SYMMETRY/=EMPTY)then
      if(size(MASK%imat(:,1))==size(MASK%imat(1,:)))then
       if(ANY((.NOT.is_diag).AND.TRANSPOSE(MASK%imat)*MASK%imat/=0))THEN
        CALL dump_message(TEXT="ERROR IN new_mask: SYMMETRY = "//TRIM(MASK%SYMMETRY)//" INCONSISTENT WITH")
        CALL write_array(MASK%imat," IMASK")
       endif
      endif
     ENDIF

    ! VECTOR
    CALL mask2vec(MASK)

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_mask_from_old(MASKOUT,MASKIN) 
    TYPE(mask_type), INTENT(INOUT) :: MASKOUT
    TYPE(mask_type), INTENT(IN)    :: MASKIN
    IF(.NOT.ASSOCIATED(MASKIN%mat)) then
     write(*,*) "ERROR IN new_mask_from_old: INPUT ISN T ALLOCATED!"
     STOP
    ENDIF
    CALL delete_mask(MASKOUT)
    CALL new_mask_from_scratch(MASKOUT,MASKIN%n1,MASKIN%n2,SYMMETRY=MASKIN%SYMMETRY)
    CALL copy_mask(MASKOUT,MASKIN)
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE copy_mask(MASKOUT,MASKIN) 

    !-------------------------------------------!
    ! CREATE MASKOUT BY COPYING EXISTING MASKIN !
    !-------------------------------------------!

    TYPE(mask_type), INTENT(INOUT) :: MASKOUT
    TYPE(mask_type), INTENT(IN)    :: MASKIN

    IF(.NOT.ASSOCIATED(MASKIN%mat))  then
       write(*,*) "ERROR IN copy_mask: INPUT  ISNT ALLOCATED!"
       STOP
    ENDIF
    IF(.NOT.ASSOCIATED(MASKOUT%mat)) then
       write(*,*) "ERROR IN copy_mask: OUTPUT ISNT ALLOCATED!"
       STOP
    ENDIF
    MASKOUT%SYMMETRY = MASKIN%SYMMETRY

    ! FULL MATRIX
    MASKOUT%n1   = MASKIN%n1
    MASKOUT%n2   = MASKIN%n2
    MASKOUT%imat = MASKIN%imat
    MASKOUT%mat  = MASKIN%mat

    ! VECTOR
    MASKOUT%nind = MASKIN%nind 

!BUG CEDRIC
    if(MASKOUT%nind<1)then
      !write(*,*) 'error, nind is = ', MASKOUT%nind
      !write(*,*) 'n1,n2 = ', MASKOUT%n1,MASKOUT%n2
    endif
!BUG CEDRIC
    IF(ASSOCIATED(MASKIN%ivec).and.MASKOUT%nind>=1)THEN 
      IF(     ASSOCIATED(MASKOUT%ivec)) DEALLOCATE(MASKOUT%ivec,STAT=istati)
      IF(.NOT.ASSOCIATED(MASKOUT%ivec)) ALLOCATE(MASKOUT%ivec(MASKOUT%nind))
      MASKOUT%ivec = MASKIN%ivec 
    ENDIF

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_mask(MASK)
    TYPE(mask_type), INTENT(INOUT) :: MASK
    IF(ASSOCIATED(MASK%imat)) DEALLOCATE(MASK%imat,STAT=istati)
    IF(ASSOCIATED(MASK%mat))  DEALLOCATE(MASK%mat,STAT=istati)
    IF(ASSOCIATED(MASK%ivec)) DEALLOCATE(MASK%ivec,STAT=istati)
!CEDRIC
    !NULLIFY(MASK%imat) 
    !NULLIFY(MASK%mat)
    !NULLIFY(MASK%ivec)
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE write_mask(MASK,UNIT)
    TYPE(mask_type),   INTENT(IN) :: MASK
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    IF(.NOT.ASSOCIATED(MASK%mat)) then
      write(*,*) "ERROR IN write_mask: INPUT ISNT ALLOCATED!"
      STOP
    ENDIF
    CALL write_array(MASK%imat,"Imask",UNIT=UNIT)
    CALL write_array(MASK%mat, " mask",UNIT=UNIT)
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE write_raw_mask(MASK,UNIT)
    TYPE(mask_type), INTENT(IN) :: MASK
    INTEGER,         INTENT(IN) :: UNIT
    IF(.NOT.ASSOCIATED(MASK%mat)) then
     write(*,*) "ERROR IN write_raw_mask: INPUT ISNT ALLOCATED!"
     STOP
    ENDIF
    WRITE(UNIT,*) MASK%n1,MASK%n2
    WRITE(UNIT,*) MASK%SYMMETRY
    WRITE(UNIT,*) MASK%imat
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE read_raw_mask(MASK,UNIT) 
    TYPE(mask_type), INTENT(INOUT) :: MASK
    INTEGER,         INTENT(IN)    :: UNIT
    INTEGER                        :: n1,n2
    CHARACTER(LEN=9)               :: SYMMETRY
    INTEGER, ALLOCATABLE           :: IMASK(:,:)
    CALL delete_mask(MASK)
    READ(UNIT,*) n1,n2
    if(n1==0.or.n2==0) stop 'read_raw_mask n1,n2=0'
    READ(UNIT,*) SYMMETRY
    ALLOCATE(IMASK(n1,n2))
    READ(UNIT,*) IMASK
    CALL new_mask_from_scratch(MASK,n1,n2,IMASK=IMASK,SYMMETRY=SYMMETRY)
    DEALLOCATE(IMASK,STAT=istati)
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE pad_imask(MASK)
    !-------------------------------------------------!
    ! THIS SYMMETRIZES AN INTEGER MASK USING SYMMETRY !
    !-------------------------------------------------!
    TYPE(mask_type), INTENT(INOUT) :: MASK
    INTEGER                        :: i1,i2,sign_cplx
    IF(.NOT.ASSOCIATED(MASK%mat)) then
     write(*,*) "ERROR IN pad_mask: INPUT ISNT ALLOCATED!"
     STOP
    ENDIF
    IF(MASK%n1/=MASK%n2)THEN
      CALL dump_message(TEXT="ERROR in symmetrize_imask: imask must be a SQUARE matrix!")
      CALL write_mask(MASK)
      STOP
    ENDIF

    IF(MASK%SYMMETRY==SYM)  sign_cplx =  1
    IF(MASK%SYMMETRY==HERM) sign_cplx = -1

    DO i1=1,MASK%n1; DO i2=i1+1,MASK%n2
      IF(MASK%imat(i1,i2)*MASK%imat(i2,i1)/=0)THEN
     IF(MASK%imat(i1,i2)/=sign_cplx*MASK%imat(i2,i1))THEN
       CALL dump_message(TEXT="ERROR in symmetrize_imask: INCONSISTENT INPUT IMASK!")
       CALL write_mask(MASK)
       STOP
     ENDIF
      ELSE IF(MASK%imat(i1,i2)/=0.OR. MASK%imat(i2,i1)/=0)THEN
     IF(MASK%imat(i1,i2)==0) MASK%imat(i1,i2) = sign_cplx * MASK%imat(i2,i1)
     IF(MASK%imat(i2,i1)==0) MASK%imat(i2,i1) = sign_cplx * MASK%imat(i1,i2)
      ENDIF
    ENDDO; ENDDO
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE build_logical_mask(MASK,IMASK)
    !---------------------------------------!
    ! BUILD LOGICAL MASK FROM INTEGER MASK  !
    !---------------------------------------!
    TYPE(mask_type),   INTENT(INOUT) :: MASK
    INTEGER, OPTIONAL, INTENT(IN)    :: IMASK(:,:)
    IF(.NOT.ASSOCIATED(MASK%mat)) then
     write(*,*) "ERROR IN build_logical_mask: INPUT ISNT ALLOCATED!"
     STOP
    ENDIF
    IF(PRESENT(IMASK))THEN
      IF(SIZE(IMASK,1)/=MASK%n1.OR.SIZE(IMASK,2)/=MASK%n2) then
        write(*,*) "ERROR IN build_logical_mask: INCONSISTENT DIMENSIONS!"
        STOP
      ENDIF
      MASK%imat = IMASK
    ENDIF
    ! PAD ZERO ELEMTS OF IMASK USING SYMMETRY
    !IF(MASK%SYMMETRY==SYM.OR.MASK%SYMMETRY==HERM) CALL pad_imask(MASK)
    !IF(MASK%SYMMETRY/=EMPTY) CALL pad_imask(MASK)
    ! BUILD LOGICAL MASK
    CALL imask2mask(MASK)
    MASK%nind = COUNT(MASK%mat)
    ! SCRAMBLE IMASK IN A USER-FRIENDLY WAY
    !CALL format_imask(MASK%imat,MASK%mat)
  END SUBROUTINE


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE imask2mask(MASK)
    !--------------------------------------!
    ! BUILD LOGICAL MASK FROM INTEGER MASK !
    !--------------------------------------!
    TYPE(mask_type), INTENT(INOUT) :: MASK
    INTEGER, ALLOCATABLE           :: list_new_elemt(:)
    LOGICAL                        :: new_elemt
    INTEGER                        :: i1,i2,iind,nind
    IF(.NOT.ASSOCIATED(MASK%mat)) then
      write(*,*) "ERROR IN imask2mask: INPUT ISNT ALLOCATED!"
      STOP
    ENDIF
    MASK%mat  = .false.
    MASK%nind = 0
    ALLOCATE(list_new_elemt(MASK%n1*MASK%n2))
    DO i2=1,MASK%n2
      DO i1=1,MASK%n1
     IF(MASK%imat(i1,i2)/=0)THEN
       ! IS THIS A NEW ARRAY ELEMENT?
       new_elemt=.true.
       DO iind=1,MASK%nind
         IF(list_new_elemt(iind)==ABS(MASK%imat(i1,i2)))THEN
           new_elemt=.false.
           EXIT
         ENDIF
       ENDDO
       IF(new_elemt)THEN
         ! IF IT IS NEW THEN MASK=TRUE
         MASK%mat(i1,i2)=.true.
         MASK%nind=MASK%nind+1
         list_new_elemt(MASK%nind)=ABS(MASK%imat(i1,i2))
       ENDIF
     ENDIF
      ENDDO
    ENDDO
    nind = COUNT(MASK%mat)
    IF(MASK%nind/=nind)THEN
      CALL dump_message(TEXT="ERROR IN imask2mask:")
      CALL dump_message(TEXT="nind = "//c2s(i2c(MASK%nind)))
      CALL write_mask(MASK)
      write(*,*) 'critical stop'
      STOP
    ENDIF

    IF(ANY(MASK%mat).AND.MASK%n1==MASK%n2)THEN
      ! NOW WE TRY TO COMPACTIFY THE MASK
      ! IN A SYMMETRIC WAY
      MASK%mat = .false.
      iind_loop: DO iind=1,MASK%nind
     ! FIRST THE UPPER RIGHT PART
     DO i2=1,MASK%n2;DO i1=1,i2
       IF(MASK%imat(i1,i2)==list_new_elemt(iind))THEN
         MASK%mat(i1,i2) = .true.
         CYCLE iind_loop
       ENDIF
     ENDDO;ENDDO
     ! THEN THE LOWER LEFT PART IF NECESSARY
      DO i1=1,MASK%n1;DO i2=1,i1
       IF(MASK%imat(i1,i2)==list_new_elemt(iind))THEN
         MASK%mat(i1,i2) = .true.
         CYCLE iind_loop
       ENDIF
      ENDDO;ENDDO
      ENDDO iind_loop
    ENDIF

  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE mask2vec(MASK)
    !--------------------------------------!
    ! CREATE THE VECTOR OF INDPDT ELEMENTS !
    !--------------------------------------!
    TYPE(mask_type), INTENT(INOUT) :: MASK
    IF(.NOT.ASSOCIATED(MASK%mat)) then
      write(*,*) "ERROR IN mask2vec: INPUT ISNT ALLOCATED!"
      STOP
    ENDIF
    MASK%nind = COUNT(MASK%mat)
    IF(ASSOCIATED(MASK%ivec).AND.SIZE(MASK%ivec)/=MASK%nind) DEALLOCATE(MASK%ivec,STAT=istati)
    IF(.NOT.ASSOCIATED(MASK%ivec).AND.MASK%nind/=0) ALLOCATE(MASK%ivec(MASK%nind))
    IF(MASK%nind/=0) MASK%ivec = PACK(MASK%imat,MASK%mat)
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE filter_mask(diag,MASK,FILTER)
    TYPE(mask_type), INTENT(INOUT) :: diag
    TYPE(mask_type), INTENT(IN)    :: MASK
    LOGICAL,         INTENT(IN)    :: FILTER(:,:)
    IF(.NOT.ASSOCIATED(MASK%mat)) then
      write(*,*) "ERROR IN filter_mask: INPUT ISNT ALLOCATED!"
      STOP
    ENDIF
    IF(SIZE(FILTER,1)/=MASK%n1.OR.SIZE(FILTER,2)/=MASK%n2) then
      write(*,*) "ERROR IN filter_mask: INCONSISTENT DIMENSIONS!"
      STOP
    ENDIF
    CALL delete_mask(diag)
    CALL new_mask(diag,MASK)
    diag%imat = MERGE(MASK%imat,0,FILTER)
    CALL build_logical_mask(diag)
    CALL mask2vec(diag)
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE slice_mask(MASKOUT,MASKIN,rbounds,cbounds) 
    TYPE(mask_type), INTENT(INOUT) :: MASKOUT
    TYPE(mask_type), INTENT(IN)    :: MASKIN
    INTEGER,         INTENT(IN)    :: rbounds(2),cbounds(2)
    INTEGER                        :: n1slice,n2slice
    CHARACTER(LEN=9)               :: SYMslice

    CALL delete_mask(MASKOUT)
    IF(.NOT.ASSOCIATED(MASKIN%mat)) then
     write(*,*) "ERROR IN slice_mask: INPUT ISNT ALLOCATED!"
     STOP
    ENDIF
    IF(rbounds(1)<LBOUND(MASKIN%mat,1).OR.rbounds(2)>UBOUND(MASKIN%mat,1)) then
      write(*,*) "ERROR IN slice_mask: INCONSISTENT ROW    BOUNDS!"
      STOP
    ENDIF
    IF(cbounds(1)<LBOUND(MASKIN%mat,2).OR.cbounds(2)>UBOUND(MASKIN%mat,2)) then
      write(*,*) "ERROR IN slice_mask: INCONSISTENT COLUMN BOUNDS!"
      STOP
    ENDIF
    n1slice = rbounds(2)-rbounds(1)+1
    n2slice = cbounds(2)-cbounds(1)+1
    IF(n1slice<=0) then
       write(*,*) "ERROR IN slice_mask: NON-NULL ROW    DIMENSIONS REQUIRED!"
       STOP
    ENDIF
    IF(n2slice<=0) then
       write(*,*) "ERROR IN slice_mask: NON-NULL COLUMN DIMENSIONS REQUIRED!"
       STOP
    ENDIF
    SYMslice = ''
    IF(ALL(rbounds==cbounds)) SYMslice = MASKIN%SYMMETRY ! SLICE=DIAGONAL BLOCK
    CALL new_mask_from_scratch(MASKOUT,n1slice,n2slice,SYMMETRY=SYMslice)

    IF(ASSOCIATED(MASKOUT%mat))  DEALLOCATE( MASKOUT%mat,STAT=istati)
    MASKOUT%mat  => MASKIN%mat( rbounds(1):rbounds(2),cbounds(1):cbounds(2))
    IF(ASSOCIATED(MASKOUT%imat)) DEALLOCATE(MASKOUT%imat,STAT=istati)
    MASKOUT%imat => MASKIN%imat(rbounds(1):rbounds(2),cbounds(1):cbounds(2))

    CALL mask2vec(MASKOUT)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

END MODULE 



