MODULE green_class

  use Block_lanczos 

  IMPLICIT NONE


  REAL(DBL),  PARAMETER, PRIVATE   :: zero=0.0_DBL,one=1.0_DBL
  LOGICAL,    PARAMETER, PRIVATE   :: F=.FALSE.,T=.TRUE.,ALLOCATE_ALL=.true.


  TYPE green_type 
  ! DYNAMICAL CORRELATIONS OF OPERATORS A,B
    INTEGER                  :: N  = 0 
    INTEGER                  :: Nw = 0 
    CHARACTER(LEN=100)       :: title = '\0'
    LOGICAL                  :: compute(2,2) = RESHAPE((/F,F,F,F/),(/2,2/))  ! independant green's fctns
    TYPE(correl_type)        :: correl (2,2)                                 ! green's fctn in all channels
    TYPE(masked_matrix_type) :: correlstat(2,2)                              ! static green's fctn in all channels particle,hole
#ifdef _complex
    COMPLEX(DBL), POINTER    :: Amean(:,:)=>NULL(),Bmean(:,:)=>NULL()        ! mean of part & hole operators 
#else
    REAL(DBL),    POINTER    :: Amean(:,:)=>NULL(),Bmean(:,:)=>NULL()        ! mean of part & hole operators 
#endif
  END TYPE


  INTERFACE new_green
    MODULE PROCEDURE new_green_rfreq
    MODULE PROCEDURE new_green_ifreq
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

  SUBROUTINE new_green_rfreq(GREEN,compute,title,N,Nw,wmin,wmax,width,STAT,IMASK,AB,force_compute)

    TYPE(green_type),  INTENT(INOUT) :: GREEN
    LOGICAL,           INTENT(IN)    :: compute(2,2)
    CHARACTER(LEN=*),  INTENT(IN)    :: title
    INTEGER,           INTENT(IN)    :: N
    CHARACTER(LEN=9),  INTENT(IN)    :: STAT
    INTEGER,           INTENT(IN)    :: Nw
    REAL(DBL),         INTENT(IN)    :: width,wmin,wmax
    INTEGER, OPTIONAL, INTENT(IN)    :: IMASK(N,N)
    LOGICAL, OPTIONAL, INTENT(IN)    :: AB
    INTEGER                          :: ipm,jpm 
    LOGICAL                          :: IS_HERM,AB_
    LOGICAL,OPTIONAL                 :: force_compute(2,2)

                    AB_ = F
    IF(PRESENT(AB)) AB_ = AB ! TRUE IF <A(z)B>, FALSE IF <A(z)A>
   
    if(messages3) write(log_unit,*) 'delete green function'
    CALL delete_green(GREEN)

    if(messages3) write(log_unit,*) 'done, now define main variables...'

    GREEN%title = TRIM(title(1:MIN(LEN_TRIM(title),100)))
    GREEN%N     = N
    GREEN%Nw    = Nw

    ! ELIMINATE REDUNDANCIES => DEDUCE THEM BY SYMMETRY: <a(z)b> = <A(-z*)B>*

    IF(compute(1,1).OR.compute(2,2))THEN ! COMPUTE ONLY <a(z)*b> BY DEFAULT
      GREEN%compute(1,1) = F
      GREEN%compute(2,2) = T
    ENDIF
    IF(compute(1,2).OR.compute(2,1))THEN ! COMPUTE ONLY <a(z)*B> BY DEFAULT
      GREEN%compute(1,2) = F
      GREEN%compute(2,1) = T
    ENDIF

    if(present(force_compute)) GREEN%compute=force_compute
   
    if(messages3) write(log_unit,*) 'build new arrays' 

   !----------------------------------------------------------------------!
   DO ipm=1,2
    DO jpm=1,2
#ifdef _complex
      IS_HERM = F
#else
      IS_HERM = F
      IF(ipm/=jpm.AND.(.NOT.AB_)) IS_HERM = T
#endif
      ! EQUAL-TIME
      CALL new_masked_matrix(GREEN%correlstat(ipm,jpm),TRIM(GREEN%title)//pm(ipm)//pm(jpm)//"_stat",N,N,IMASK=IMASK,IS_HERM=IS_HERM)
      CALL masked_matrix2vec(GREEN%correlstat(ipm,jpm)) ! dummy to create vec of ind elemts

      ! DYNAMIC 
       IF(GREEN%compute(ipm,jpm).OR.GREEN%compute(3-ipm,3-jpm).OR.ALLOCATE_ALL)THEN 
          CALL new_correl(GREEN%correl(ipm,jpm),TRIM(GREEN%title)//pm(ipm)//pm(jpm),N,Nw,wmin,wmax,width,STAT,IMASK=IMASK,AB=AB)
          CALL correl2vec(GREEN%correl(ipm,jpm)) ! dummy to create vec of ind elemts
       ENDIF

    ENDDO
   ENDDO
   !----------------------------------------------------------------------!

    if(messages3) write(log_unit,*) '...new frequ rfrequ....'
    IF(N/=0)THEN
      if(associated(GREEN%Amean)) deallocate(GREEN%Amean)
      if(associated(GREEN%Bmean)) deallocate(GREEN%Bmean)
      ALLOCATE(GREEN%Amean(N,2),GREEN%Bmean(N,2))
      GREEN%Amean = zero; GREEN%Bmean = zero
    ENDIF
    if(messages3) write(log_unit,*) ' .. allocation  done .. '

  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE new_green_ifreq(GREEN,compute,title,N,Nw,beta,STAT,IMASK,AB,force_compute) 
    TYPE(green_type),  INTENT(INOUT) :: GREEN
    LOGICAL,           INTENT(IN)    :: compute(2,2) 
    CHARACTER(LEN=*),  INTENT(IN)    :: title
    INTEGER,           INTENT(IN)    :: N
    CHARACTER(LEN=9),  INTENT(IN)    :: STAT
    INTEGER,           INTENT(IN)    :: Nw
    REAL(DBL),         INTENT(IN)    :: beta
    INTEGER, OPTIONAL, INTENT(IN)    :: IMASK(N,N)
    LOGICAL, OPTIONAL, INTENT(IN)    :: AB
    INTEGER                          :: ipm,jpm 
    LOGICAL                          :: IS_HERM,AB_
    LOGICAL,OPTIONAL                 :: force_compute(2,2)
    
    AB_ = F
    IF(PRESENT(AB)) AB_ = AB ! TRUE IF <A(z)B>, FALSE IF <A(z)A>

    if(messages3) write(log_unit,*) 'delete green function'

    CALL delete_green(GREEN)
    
    GREEN%title   = TRIM(title(1:MIN(LEN_TRIM(title),100)))
    GREEN%N       = N
    GREEN%Nw      = Nw

    ! ELIMINATE REDUNDANCIES => DEDUCE THEM BY SYMMETRY: <a(z)b> = <A(-z*)B>*

    IF(compute(1,1).OR.compute(2,2))THEN ! COMPUTE ONLY <a(z)*b> BY DEFAULT
      GREEN%compute(1,1) = F
      GREEN%compute(2,2) = T
    ENDIF
    IF(compute(1,2).OR.compute(2,1))THEN ! COMPUTE ONLY <a(z)*B> BY DEFAULT
      GREEN%compute(1,2) = F
      GREEN%compute(2,1) = T
    ENDIF
 
    if(present(force_compute)) GREEN%compute=force_compute
 
    if(messages3) write(log_unit,*) 'build new arrays'
 
   !----------------------------------------------------------------------!
   DO ipm=1,2
    DO jpm=1,2 
#ifdef _complex
      IS_HERM = F
#else
      IS_HERM = F
      IF(ipm/=jpm.AND.(.NOT.AB_)) IS_HERM = T
#endif
      ! EQUAL-TIME
      CALL new_masked_matrix(GREEN%correlstat(ipm,jpm),TRIM(GREEN%title)//pm(ipm)//pm(jpm)//"_stat",N,N,IMASK=IMASK,IS_HERM=IS_HERM)
      CALL masked_matrix2vec(GREEN%correlstat(ipm,jpm)) ! dummy to create vec of ind elemts

      ! DYNAMIC 
      IF(GREEN%compute(ipm,jpm).OR.GREEN%compute(3-ipm,3-jpm).OR.ALLOCATE_ALL)THEN 
         CALL new_correl(GREEN%correl(ipm,jpm),TRIM(GREEN%title)//pm(ipm)//pm(jpm),N,Nw,beta,STAT,IMASK=IMASK,AB=AB)
         CALL correl2vec(GREEN%correl(ipm,jpm)) ! dummy to create vec of ind elemts
      ENDIF

    ENDDO
   ENDDO
   !----------------------------------------------------------------------!

    if(messages3) write(log_unit,*) '...new frequ rfrequ....'
    IF(N/=0)THEN
      if(associated(GREEN%Amean)) deallocate(GREEN%Amean)
      if(associated(GREEN%Bmean)) deallocate(GREEN%Bmean)
      ALLOCATE(GREEN%Amean(N,2),GREEN%Bmean(N,2))
      GREEN%Amean = zero; GREEN%Bmean = zero
    ENDIF
    if(messages3) write(log_unit,*) '...allocation done...'

  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE delete_green(GREEN)
    TYPE(green_type), INTENT(INOUT) :: GREEN
    INTEGER                         :: ipm,jpm

      if(messages3) write(log_unit,*) '...delete correlations contained in green...'

      DO ipm=1,2
       DO jpm=1,2 
       if(messages3) write(log_unit,*) 'delete dynamical part'
        IF(GREEN%compute(ipm,jpm).OR.GREEN%compute(3-ipm,3-jpm)) CALL delete_correl(GREEN%correl(ipm,jpm))
        if(messages3) write(log_unit,*) 'delete equal time correlations'
        CALL delete_masked_matrix(GREEN%correlstat(ipm,jpm))
        if(messages3) write(log_unit,*) 'done......'
       ENDDO
      ENDDO

     if(messages3) write(log_unit,*) '...delete green function s averages...'
     IF(ASSOCIATED(GREEN%Amean)) DEALLOCATE(GREEN%Amean)
     IF(ASSOCIATED(GREEN%Bmean)) DEALLOCATE(GREEN%Bmean)

  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE copy_green(GREENOUT,GREENIN)
    TYPE(green_type), INTENT(IN)    :: GREENIN
    TYPE(green_type), INTENT(INOUT) :: GREENOUT
    INTEGER                         :: ipm,jpm

    GREENOUT%title   = GREENIN%title
    GREENOUT%N       = GREENIN%N
    GREENOUT%Nw      = GREENIN%Nw
    GREENOUT%compute = GREENIN%compute

    DO ipm=1,2  
     DO jpm=1,2
       ! DYNAMIC
        if(messages3) write(log_unit,*) ' -----> copy_green : correl '
        IF(GREENOUT%compute(ipm,jpm).OR.GREENOUT%compute(3-ipm,3-jpm)) then
          CALL copy_correl(GREENOUT%correl(ipm,jpm),GREENIN%correl(ipm,jpm))
        ENDIF
        if(messages3) write(log_unit,*) ' -----> copy_green : mask ' 
       ! EQUAL-TIME
       CALL copy_masked_matrix(GREENOUT%correlstat(ipm,jpm),GREENIN%correlstat(ipm,jpm))
     ENDDO
    ENDDO

    GREENOUT%Amean = GREENIN%Amean ; GREENOUT%Bmean = GREENIN%Bmean

  END SUBROUTINE


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE pad_green(GREENOUT,GREENIN)

    !---------------------------------------------!
    ! PAD REDUNDANT MATRIX ELEMENTS IF NECESSARY  !
    ! USING <a(z)*b> = (-1)^F <A(-z*)*B>*         !
    !---------------------------------------------!

    TYPE(green_type), INTENT(INOUT)        :: GREENOUT
    TYPE(green_type), INTENT(IN), OPTIONAL :: GREENIN
    INTEGER                                :: ipm,jpm,mipm,mjpm,iw,miw

    ! FIRST PAD IN (ipm,jpm) SECTOR USING EXTERNAL GREEN'S FUNCTION IF NECESSARY
 
    DO ipm=1,2; DO jpm=1,2
      IF(PRESENT(GREENIN))THEN
        GREENOUT%Amean = GREENIN%Amean ! WARNING: OVERRIDE!... SO PAD GREENOUT WITH GREENIN BEFORE COMPUTING GREENIN! 
        GREENOUT%Bmean = GREENIN%Bmean ! WARNING: OVERRIDE!... SO PAD GREENOUT WITH GREENIN BEFORE COMPUTING GREENIN! 
        ! EQUAL-TIME
        CALL pad_masked_matrix(GREENOUT%correlstat(ipm,jpm),GREENIN%correlstat(ipm,jpm))
        ! DYNAMIC
        IF(  (GREENOUT%compute(ipm,jpm).OR.GREENOUT%compute(3-ipm,3-jpm))  &
        .AND.( GREENIN%compute(ipm,jpm).OR. GREENIN%compute(3-ipm,3-jpm))) &
        CALL pad_correl(GREENOUT%correl(ipm,jpm),GREENIN%correl(ipm,jpm))
      ELSE
        ! EQUAL-TIME
        CALL pad_masked_matrix(GREENOUT%correlstat(ipm,jpm))
        ! DYNAMIC
        IF(  (GREENOUT%compute(ipm,jpm).OR.GREENOUT%compute(3-ipm,3-jpm))  &
        .AND. (GREENIN%compute(ipm,jpm).OR. GREENIN%compute(3-ipm,3-jpm))) &
        CALL pad_correl(GREENOUT%correl(ipm,jpm))
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
!**************************************************************************


  SUBROUTINE write_green(GREEN)
    TYPE(green_type), INTENT(IN) :: GREEN
    INTEGER                      :: ipm,jpm
    DO ipm=1,2; DO jpm=1,2; IF(GREEN%compute(ipm,jpm).OR.GREEN%compute(3-ipm,3-jpm))THEN
      CALL write_correl(GREEN%correl(ipm,jpm))
    ENDIF; ENDDO; ENDDO
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE glimpse_green(GREEN,UNIT)
    TYPE(green_type), INTENT(IN) :: GREEN
    INTEGER,          OPTIONAL   :: UNIT
    INTEGER                      :: ipm,jpm
    DO ipm=1,2; DO jpm=1,2; IF(GREEN%compute(ipm,jpm).OR.GREEN%compute(3-ipm,3-jpm))THEN
      CALL glimpse_correl(GREEN%correl(ipm,jpm),UNIT)
    ENDIF; ENDDO; ENDDO
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE slice_green(GREENOUT,GREENIN,rbounds,cbounds)
    TYPE(green_type), INTENT(INOUT) :: GREENOUT
    TYPE(green_type), INTENT(IN)    :: GREENIN
    INTEGER,          INTENT(IN)    :: rbounds(2),cbounds(2)
    INTEGER                         :: ipm,jpm,n1slice,n2slice

    CALL delete_green(GREENOUT)
    n1slice = rbounds(2)-rbounds(1)+1 
    n2slice = cbounds(2)-cbounds(1)+1

    IF(n1slice/=n2slice) STOP "ERROR IN slice_green: SQUARE   SLICE      REQUIRED!"
    IF(n1slice<=0)       STOP "ERROR IN slice_green: NON-NULL DIMENSIONS REQUIRED!"

    GREENOUT%title    = TRIM(GREENIN%title)//"_"//c2s(i2c(rbounds(1)))//"_"//c2s(i2c(rbounds(2)))//"_"//c2s(i2c(cbounds(1)))//"_"//c2s(i2c(cbounds(2)))
    GREENOUT%compute  = GREENIN%compute
    GREENOUT%Nw       = GREENIN%Nw
    GREENOUT%N        = n1slice

    DO ipm=1,2 
     DO jpm=1,2
      IF(GREENIN%compute(ipm,jpm).OR.GREENIN%compute(3-ipm,3-jpm))THEN
       ! DYNAMIC
         CALL slice_correl(GREENOUT%correl(ipm,jpm),GREENIN%correl(ipm,jpm),rbounds,cbounds)
       ! EQUAL-TIME
        CALL slice_masked_matrix(GREENOUT%correlstat(ipm,jpm),GREENIN%correlstat(ipm,jpm),rbounds,cbounds)
      ENDIF
     ENDDO 
    ENDDO

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
