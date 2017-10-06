MODULE bath_class_vec


  USE bath_class

  IMPLICIT NONE

  REAL(DBL),    PARAMETER, PRIVATE                 :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,      PARAMETER, PRIVATE                 :: F=.FALSE.,T=.TRUE.

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !$$ CONVERT BATH TO VEC AND BACK $$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

CONTAINS


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE bath2vec(BATH)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ CAST BATH INTO REAL VECTOR $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    TYPE(bath_type), INTENT(INOUT) :: BATH
    INTEGER                        :: spin,offset
    LOGICAL, SAVE                  :: first_call = T

    IF(.NOT.ASSOCIATED(BATH%Eb)) STOP "ERROR IN bath2vec: INPUT ISNT ASSOCIATED!"
    
    IF(.NOT.ASSOCIATED(BATH%vec)) ALLOCATE(BATH%vec(BATH%nparam))
    BATH%vec = zero
 
    ! GATHER ALL VECTORS OF INDPDT ELEMTS IN A SINGLE VECTOR OF SIZE nparam

    offset = 0
    DO spin=1,2
      ! HYBRIDIZATION
      CALL add_mm2vec(offset,BATH%Vbc(spin),BATH%vec)  
      CALL add_mm2vec(offset,BATH%PVbc(spin),BATH%vec)
      ! BATH ENERGY
      CALL add_mm2vec(offset, BATH%Eb(spin),BATH%vec)  
    ENDDO

    ! PAIRING ENERGY
    IF(BATH%SUPER) then
        CALL add_mm2vec(offset,BATH%Pb,BATH%vec)  
    endif

    IF(first_call)THEN
      if(.not.ALL_FIRST_CALL) first_call = F
      IF(offset/=BATH%nparam)THEN
        CALL dump_message(TEXT="ERROR IN bath2vec: nparam_verif ="//c2s(i2c(offset))//" /= nparam ="//c2s(i2c(BATH%nparam)))
        write(*,*) 'found... parameters : ', offset
        write(*,*) 'should be           : ', BATH%nparam
        STOP 'total amount of parameters do not match the imposed number of parameter in the bath, stop'
      ENDIF
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

  SUBROUTINE vec2bath(BATH)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ EXTRACT BATH OUT OF REAL VECTOR $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    TYPE(bath_type)                :: BATH
    INTEGER                        :: spin,offset
    LOGICAL, SAVE                  :: first_call = T

    IF(.NOT.ASSOCIATED(BATH%Eb)) STOP "ERROR IN vec2bath: INPUT ISNT ASSOCIATED!"

    ! FIRST EXTRACT THE DIFFERENT VECTORS OF INDPDT ELEMTS FROM THE SINGLE VECTOR

    offset = 0
    DO spin=1,2
      ! HYBRIDIZATION
      CALL extract_mmvec_from_vec(offset,BATH%Vbc(spin),BATH%vec)
      CALL extract_mmvec_from_vec(offset,BATH%PVbc(spin),BATH%vec)
      ! BATH ENERGY
      CALL extract_mmvec_from_vec(offset, BATH%Eb(spin),BATH%vec)
    ENDDO

    IF(BATH%SUPER) CALL extract_mmvec_from_vec(offset,BATH%Pb,BATH%vec)

    ! FILL MISSING SPIN IF NECESSARY
    IF(BATH%Vbc(2)%rc%MASK%nind==0  .AND. BATH%Vbc(1)%rc%MASK%nind/=0)  CALL vec2masked_matrix(BATH%Vbc(2), BATH%Vbc(1))
    IF(BATH%Vbc(1)%rc%MASK%nind==0  .AND. BATH%Vbc(2)%rc%MASK%nind/=0)  CALL vec2masked_matrix(BATH%Vbc(1), BATH%Vbc(2))
    IF(BATH%PVbc(2)%rc%MASK%nind==0 .AND. BATH%PVbc(1)%rc%MASK%nind/=0) CALL vec2masked_matrix(BATH%PVbc(2), BATH%PVbc(1))
    IF(BATH%PVbc(1)%rc%MASK%nind==0 .AND. BATH%PVbc(2)%rc%MASK%nind/=0) CALL vec2masked_matrix(BATH%PVbc(1), BATH%PVbc(2))
    IF(BATH%Eb (2)%rc%MASK%nind==0  .AND. BATH%Eb (1)%rc%MASK%nind/=0)  CALL vec2masked_matrix(BATH%Eb (2), BATH%Eb (1))
    IF(BATH%Eb (1)%rc%MASK%nind==0  .AND. BATH%Eb (2)%rc%MASK%nind/=0)  CALL vec2masked_matrix(BATH%Eb (1), BATH%Eb (2))
 
    IF(first_call)THEN
      if(.not.ALL_FIRST_CALL) first_call = F
      IF(offset/=BATH%nparam)THEN
        CALL dump_message(TEXT="ERROR IN vec2bath: nparam_verif ="//c2s(i2c(offset))//" /= nparam ="//c2s(i2c(BATH%nparam)))
        STOP
      ENDIF
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

  SUBROUTINE add_mm2vec(offset,MM,vec) 

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ APPEND MM%vec TO END OF vec $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    INTEGER,                  INTENT(INOUT) :: offset
    TYPE(masked_matrix_type), INTENT(INOUT) :: MM
    REAL(DBL),                INTENT(INOUT) :: vec(:)
    INTEGER                                 :: nind,ninddiag,nindoffdiag ! for clarity

    IF(.NOT.ASSOCIATED(MM%rc%mat)) STOP "ERROR IN add_mm2vec: INPUT ISNT ALLOCATED!"

    nind = MM%rc%MASK%nind

    !--------------------------------------------------------------------!
    IF(nind/=0)THEN
      CALL masked_matrix2vec(MM)

#ifdef _complex

      ! COMPLEX MATRIX

      IF(.NOT.MM%rc%IS_HERM)THEN 
         ! EVERYBODY IS COMPLEX 
          vec(offset+1:offset+nind)  =  DBLE(MM%rc%vec)    
          offset                     =  offset + nind      
          vec(offset+1:offset+nind)  =  AIMAG(MM%rc%vec)   
          offset                     =  offset + nind      
      ELSE 

         ! DIAGONAL IS REAL/OFF-DIAGONAL IS COMPLEX 
         ninddiag    = MM%rc%MASKdiag%nind
         nindoffdiag = MM%rc%MASKoffdiag%nind

         IF(ninddiag/=0) THEN
           vec(offset+1:offset+ninddiag)     = MM%rc%vecdiag 
           offset                            = offset + ninddiag
         ENDIF

         IF(nindoffdiag/=0) THEN
            vec(offset+1:offset+nindoffdiag) = DBLE(MM%rc%vecoffdiag) 
            offset                           = offset + nindoffdiag
            vec(offset+1:offset+nindoffdiag) = AIMAG(MM%rc%vecoffdiag) 
            offset                           = offset + nindoffdiag
         ENDIF

      ENDIF
#else 
      ! REAL    MATRIX
      vec(offset+1:offset+nind) = MM%rc%vec 
      offset                    = offset + nind
#endif
    ENDIF
    !--------------------------------------------------------------------!


  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE extract_mmvec_from_vec(offset,MM,vec)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ EXTRACT MM%vec FROM END OF vec $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    INTEGER,                  INTENT(INOUT) :: offset
    TYPE(masked_matrix_type), INTENT(INOUT) :: MM
    REAL(DBL),                INTENT(INOUT) :: vec(:) 
    INTEGER                                 :: nind,ninddiag,nindoffdiag,i ! for clarity only

    IF(.NOT.ASSOCIATED(MM%rc%mat)) STOP "ERROR IN extract_mmvec_from_vec: INPUT ISNT ALLOCATED!"

    nind        =    MM%rc%MASK%nind
    ninddiag    =    0
    nindoffdiag =    0

    IF(nind/=0)THEN
     CALL masked_matrix2vec(MM)

#ifdef _complex

      IF(.NOT.MM%rc%IS_HERM)THEN 
        ! EVERYBODY IS COMPLEX 
        MM%rc%vec = CMPLX(vec(offset+1:offset+nind),vec(offset+nind+1:offset+nind*2),8)
        offset    = offset + nind*2
      ELSE 
        ! DIAGONAL IS REAL/OFF-DIAGONAL IS COMPLEX 
        ninddiag    =     MM%rc%MASKdiag%nind
        nindoffdiag =  MM%rc%MASKoffdiag%nind

        if(offset+1>size(vec)) stop 'error extract_mmvec_from_vec size does not match'

        if(size(MM%rc%vec)<ninddiag+nindoffdiag) then
            write(*,*) 'case diag is real, off diag are complex'
            write(*,*) '-------------- ERRRO INFO ---------------'
            write(*,*) 'This error typically happens if in the ed_hybrid input file'
            write(*,*) 'the same variable appears in the diagonal and in the off-diagonal'
            write(*,*) 'part of the mask matrices....check it!'
            write(*,*) 'Reason why, is that for the complex case, the bath Epsilon matrix is hermitian'
            write(*,*) 'and therefore a diag variable is a real number, and an ofddiag one a complex number'
            write(*,*) '-----------------------------------------'
            write(*,*) 'MM%rc%vec size : ', size(MM%rc%vec)
            write(*,*) 'should be      : ', ninddiag+2*nindoffdiag
            write(*,*) 'in diag        : ', ninddiag
            write(*,*) 'off diag       :  ',nindoffdiag
            write(*,*) 'IS HERM        :  ',MM%rc%IS_HERM
            write(*,*) 'EXTRACT MM%vec from vec error size'
            stop 
        endif
        !----------------------------------------------------------------------!
          IF(ninddiag/=0)THEN
             MM%rc%vecdiag = vec(offset+1:offset+ninddiag)
             offset        = offset + ninddiag
          ENDIF
          IF(nindoffdiag>0)THEN
             MM%rc%vecoffdiag = CMPLX(vec(offset+1:offset+nindoffdiag),vec(offset+nindoffdiag+1:offset+nindoffdiag*2),8)
             offset = offset + nindoffdiag*2
             if(.not.associated(MM%rc%vec)) stop 'MM%rc%vec not associated'
          ENDIF
        !----------------------------------------------------------------------!
        call gather_diag_offdiag_vec(MM%rc)

      ENDIF
#else
      ! REAL    MATRIX
      MM%rc%vec = vec(offset+1:offset+nind)
      offset    = offset + nind
#endif
    ENDIF

    CALL vec2masked_matrix(MM)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE write_bathvec(BATH,UNIT)

    TYPE(bath_type),   INTENT(IN) :: BATH
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    CHARACTER(LEN=100), SAVE      :: fmt_vec
    INTEGER                       :: unit_,iparam
    LOGICAL, SAVE                 :: first_call = T

    IF(.NOT.ASSOCIATED(BATH%Eb)) STOP "ERROR IN write_bathvec: INPUT ISNT ASSOCIATED!"

    IF(first_call)THEN
      if(.not.ALL_FIRST_CALL) first_call = F
      WRITE(fmt_vec,'(a)') "(a,"//c2s(i2c(BATH%nparam))//"(f0.6,2X))" 
    ENDIF

    unit_ = log_unit
    IF(PRESENT(UNIT)) unit_ = UNIT
    WRITE(unit_,fmt_vec) "# BATHvec = ",BATH%vec
    CALL flush(unit_)

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
