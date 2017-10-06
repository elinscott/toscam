
MODULE correl_class

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !$$ CORRELATION CLASS FOR SPIN INDEPENDANT CORRELATIONS $$
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  USE frequency_class
  USE impurity_class  
  USE mesh

  IMPLICIT NONE


  REAL(DBL), PARAMETER, PRIVATE :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,   PARAMETER, PRIVATE :: F=.FALSE.,T=.TRUE.


  TYPE correl_type
    CHARACTER(LEN=100)             :: title = '\0'
    INTEGER                        :: N     = 0 
    INTEGER                        :: Nw    = 0 
    CHARACTER(LEN=9)               :: stat  = '\0'
    TYPE(masked_cplx_matrix_type)  :: MM
    COMPLEX(DBL),    POINTER       :: fctn(:,:,:) => NULL() ! matrix correlation function
    COMPLEX(DBL),    POINTER       :: vec(:,:)    => NULL() ! vector of independant matrix elements
    TYPE(freq_type)                :: freq ! frequency 
  END TYPE


  INTERFACE new_correl
    MODULE PROCEDURE new_correl_from_scratch_rfreq
    MODULE PROCEDURE new_correl_from_scratch_ifreq
    MODULE PROCEDURE new_correl_from_scratch
    MODULE PROCEDURE new_correl_from_old
  END INTERFACE




CONTAINS


!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE correl2vec(CORREL,extmat,vecext)

    !-----------------------------------------------------------!
    ! CAST CORRELATION MATRIX IN VECTOR OF INDEPENDANT ELEMENTS !
    !-----------------------------------------------------------!

    TYPE(correl_type)             :: CORREL
    TYPE(masked_cplx_matrix_type) :: corr
    INTEGER                       :: iw
    complex(8),optional           :: extmat(:,:,:),vecext(:,:)

    IF(.NOT.ASSOCIATED(CORREL%vec)) ALLOCATE(CORREL%vec(CORREL%MM%MASK%nind,CORREL%Nw))
    CALL new_masked_cplx_matrix(corr,CORREL%MM)
    CALL masked_cplx_matrix2vec(corr)
 
    if(present(extmat))then
     if(any(shape(extmat)/=shape(CORREL%fctn))) then
       write(*,*) 'error shape of array in correl2vec, shape extmat :',shape(extmat)
       write(*,*) 'shape correl%func : ', shape(correl%fctn)
       stop
     endif
    endif

   if(present(vecext))then
     if(any(shape(vecext)/=shape(CORREL%vec))) then
       write(*,*) 'error shape of array in correl2vec, shape vecext :',shape(vecext)
       write(*,*) 'shape correl%vec : ', shape(correl%vec)
       stop
     endif
    endif

    DO iw=1,CORREL%Nw
      if(.not.present(extmat))then
       corr%mat = CORREL%fctn(:,:,iw)
      else
       corr%mat = extmat(:,:,iw)
      endif
      CALL masked_cplx_matrix2vec(corr)
      if(.not.present(vecext))then
       CORREL%vec(:,iw) = corr%vec
      else
       vecext(:,iw)     = corr%vec
      endif
    ENDDO

    CALL delete_masked_cplx_matrix(corr)

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE new_correl_from_scratch_ifreq(CORREL,title,N,Nw,beta,STAT,IMASK,AB) 

    !-----------------------------------------------------!
    ! CREATE NEW MATSUBARA FREQ. CORRELATION FROM SCRATCH !
    !-----------------------------------------------------!

    TYPE(correl_type), INTENT(INOUT) :: CORREL
    CHARACTER(LEN=*),  INTENT(IN)    :: title
    INTEGER,           INTENT(IN)    :: N
    CHARACTER(LEN=9),  INTENT(IN)    :: STAT
    INTEGER,           INTENT(IN)    :: Nw
    REAL(DBL),         INTENT(IN)    :: beta
    INTEGER, OPTIONAL, INTENT(IN)    :: IMASK(N,N)
    LOGICAL, OPTIONAL, INTENT(IN)    :: AB

    if(Nw==0) stop 'new correl from scratch ifrequ Nw=0'
    CALL new_correl_from_scratch(CORREL,title,N,Nw,STAT,IMASK,AB=AB)
    CALL new_freq(CORREL%freq,Nw,beta,STAT)

    DEALLOCATE(CORREL%MM%mat,STAT=istati)
    CORREL%MM%mat => CORREL%fctn(:,:,CORREL%freq%iwprint)

  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE new_correl_from_scratch_rfreq(CORREL,title,N,Nw,wmin,wmax,width,STAT,IMASK,AB) 

    !------------------------------------------------!
    ! CREATE NEW REAL FREQ. CORRELATION FROM SCRATCH !
    !------------------------------------------------!

    TYPE(correl_type), INTENT(INOUT) :: CORREL
    CHARACTER(LEN=*),  INTENT(IN)    :: title
    INTEGER,           INTENT(IN)    :: N
    CHARACTER(LEN=9),  INTENT(IN)    :: STAT
    INTEGER,           INTENT(IN)    :: Nw
    REAL(DBL),         INTENT(IN)    :: wmax,width
    INTEGER, OPTIONAL, INTENT(IN)    :: IMASK(N,N)
    LOGICAL, OPTIONAL, INTENT(IN)    :: AB
    REAL(DBL)         ,INTENT(IN)    :: wmin

    if(Nw==0) stop 'new correl from scratch rfrequNw=0'
    CALL new_correl_from_scratch(CORREL,title,N,Nw,STAT,IMASK,AB=AB)
    IF(width>zero)THEN
     CALL new_freq(CORREL%freq,Nw,wmin,wmax,width,RETARDED)
    ELSE
     CALL new_freq(CORREL%freq,Nw,wmin,wmax,width,ADVANCED)
    ENDIF
 
    DEALLOCATE(CORREL%MM%mat,STAT=istati)
    CORREL%MM%mat => CORREL%fctn(:,:,CORREL%freq%iwprint)

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE new_correl_from_scratch(CORREL,title,N,Nw,STAT,IMASK,AB)

  !---------------------------------------------------------------!
  ! CREATE NEW CORRELATION (WITHOUT FREQUENCY ARRAY) FROM SCRATCH !
  !---------------------------------------------------------------!

    TYPE(correl_type), INTENT(INOUT) :: CORREL
    CHARACTER(LEN=*),  INTENT(IN)    :: title
    INTEGER,           INTENT(IN)    :: N,Nw
    CHARACTER(LEN=9),  INTENT(IN)    :: STAT
    INTEGER, OPTIONAL, INTENT(IN)    :: IMASK(N,N)
    LOGICAL, OPTIONAL, INTENT(IN)    :: AB
    LOGICAL                          :: is_diag(N,N),AB_
    INTEGER                          :: imask_(N,N)
    INTEGER                          :: iind,mu,nu

    if(Nw==0) stop 'new correl from scratch Nw=0'
    CALL delete_correl(CORREL)

    CORREL%title = TRIM(ADJUSTL(title))
    CORREL%stat  = STAT
    CORREL%N     = N
    CORREL%Nw    = Nw

    ALLOCATE(CORREL%fctn(N,N,Nw))
    CORREL%fctn  = zero

                    AB_ = F
    IF(PRESENT(AB)) AB_ = AB

!=======================================================================================================!
#ifdef _complex
    IF(PRESENT(IMASK))THEN
      CALL new_diag(is_diag,N)
      IF((.NOT.AB_).AND.ANY((.NOT.is_diag).AND.(IMASK==TRANSPOSE(IMASK)).AND.(IMASK/=0))) &
      CALL dump_message(TEXT="WARNING IN new_correl: CPLX H => "//TRIM(ADJUSTL(CORREL%title))//"/=TRANSPOSE("//TRIM(ADJUSTL(CORREL%title))//") A PRIORI!")
    ENDIF
#else
    IF(PRESENT(IMASK))THEN
      IF(.NOT.AB_.AND.ANY(IMASK/=TRANSPOSE(IMASK))) then
       write(*,*) "ERROR IN new_correl: REAL H => "//TRIM(ADJUSTL(CORREL%title))//"=TRANSPOSE("//TRIM(ADJUSTL(CORREL%title))//") !"
       if(strongstop) stop
      endif
    ENDIF
#endif
!=======================================================================================================!

    CALL new_masked_cplx_matrix(CORREL%MM,CORREL%title,N,N,IMASK=IMASK,IS_HERM=F)

  RETURN
  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE new_correl_from_old(CORRELOUT,CORRELIN) 

    TYPE(correl_type), INTENT(INOUT) :: CORRELOUT
    TYPE(correl_type), INTENT(IN)    :: CORRELIN

    if(CORRELIN%Nw==0) stop 'new correl from old Nw=0'

  !-------------------------------------------------------!
  ! THE AB=T FLAG AVOIDS THE SYMMETRY TESTS IN new_correl !
  !-------------------------------------------------------!

    SELECT CASE(CORRELIN%freq%title)
     CASE(BOSONIC,FERMIONIC)
       CALL new_correl_from_scratch_ifreq(CORRELOUT,CORRELIN%title,CORRELIN%N,CORRELIN%Nw,&
       CORRELIN%freq%beta,CORRELIN%stat,IMASK=CORRELIN%MM%MASK%imat,AB=T)
     CASE(RETARDED,ADVANCED)
       CALL new_correl_from_scratch_rfreq(CORRELOUT,CORRELIN%title,CORRELIN%N,CORRELIN%Nw,CORRELIN%freq%wmin,&
       CORRELIN%freq%wmax,CORRELIN%freq%width,CORRELIN%stat,IMASK=CORRELIN%MM%MASK%imat,AB=T)
    END SELECT

    CALL copy_correl(CORRELOUT,CORRELIN)

  RETURN
  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE copy_correl(CORRELOUT,CORRELIN) 

  !---------------------------------------------------!
  ! CREATE CORRELOUT BY COPYING EXISTING ONE CORRELIN !
  !---------------------------------------------------!

    TYPE(correl_type), INTENT(IN)    :: CORRELIN
    TYPE(correl_type), INTENT(INOUT) :: CORRELOUT
    integer                          :: i1,i2,i(3)

    CORRELOUT%N    = CORRELIN%N
    CALL copy_frequency(CORRELOUT%freq,CORRELIN%freq)
    CORRELOUT%stat = CORRELIN%stat

    if(ASSOCIATED(CORRELIN%fctn)) then
     if(ASSOCIATED(CORRELOUT%fctn)) DEALLOCATE(CORRELOUT%fctn,STAT=istati)
     i=shape(CORRELIN%fctn); ALLOCATE(CORRELOUT%fctn(i(1),i(2),i(3)))
     CORRELOUT%fctn = CORRELIN%fctn
    endif

    if(ASSOCIATED(CORRELOUT%fctn)) then
!CEDRIC BUG
!    IF(ASSOCIATED(CORRELOUT%MM%mat)) DEALLOCATE(CORRELOUT%MM%mat,STAT=istati)
     CORRELOUT%MM%mat => CORRELOUT%fctn(:,:,CORRELOUT%freq%iwprint)
     CORRELOUT%MM%N1 = size(CORRELOUT%fctn,1)
     CORRELOUT%MM%N2 = size(CORRELOUT%fctn,2)
    endif

    IF(ASSOCIATED(CORRELIN%vec))THEN
      IF(ASSOCIATED(CORRELOUT%vec)) DEALLOCATE(CORRELOUT%vec,STAT=istati)
      i(1:2)=shape(CORRELIN%vec); ALLOCATE(CORRELOUT%vec(i(1),i(2))) 
      CORRELOUT%vec = CORRELIN%vec
    ENDIF

  RETURN
  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE delete_correl(CORREL)

    TYPE(correl_type)  :: CORREL

    if(ASSOCIATED(CORREL%MM%mat,CORREL%fctn(:,:,CORREL%freq%iwprint)))then
      NULLIFY(CORREL%MM%mat)
    endif

    if(messages3) write(log_unit,*) '......delete masked matrix.....'
    CALL delete_masked_matrix_(CORREL%MM)
    if(messages3) write(log_unit,*) '......delete correlations......'
    IF(ASSOCIATED(CORREL%fctn)) DEALLOCATE(CORREL%fctn,STAT=istati)
    if(messages3) write(log_unit,*) '......delete vec...............'
    IF(ASSOCIATED(CORREL%vec))  DEALLOCATE(CORREL%vec,STAT=istati)
    if(messages3) write(log_unit,*) '......delete frequencies.......'
    CALL delete_frequency(CORREL%freq)
    if(messages3) write(log_unit,*) '......now return...............'
 
  RETURN
  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE write_raw_correl(CORREL,UNIT)

    ! WRITE CORRELATION TO FILE IN RAW (UNFORMATTED) FORM

    TYPE(correl_type), INTENT(IN) :: CORREL
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    INTEGER                       :: unit_

    ! OPEN OUTPUT FILE
    IF(PRESENT(UNIT))THEN
      unit_ = UNIT
    ELSE
      CALL open_safe(unit_,"./ED_out/"//TRIM(ADJUSTL(CORREL%title))//"_raw.dat","UNKNOWN","WRITE",get_unit=.true.)
    ENDIF

    ! WRITE RAW CORRELATION TO OUTPUT FILE
    WRITE(unit_,*) CORREL%title        
    WRITE(unit_,*) CORREL%N           
    WRITE(unit_,*) CORREL%stat      
    CALL write_raw_masked_cplx_matrix(CORREL%MM,unit_)
    CALL write_raw_freq(CORREL%freq,unit_)
    WRITE(unit_,*) CORREL%fctn

    ! CLOSE OUTPUT FILE
    IF(.NOT.PRESENT(UNIT)) CALL close_safe(unit_)

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************


  SUBROUTINE read_raw_correl(CORREL,UNIT,FILEIN) 

    !--------------------------------------!
    ! CREATE NEW CORRELATION FROM RAW FILE !
    !--------------------------------------!

    TYPE(correl_type),          INTENT(INOUT) :: CORREL
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: FILEIN
    INTEGER,          OPTIONAL, INTENT(IN)    :: UNIT 
    CHARACTER(LEN=100)                        :: title
    CHARACTER(LEN=9)                          :: stat
    INTEGER                                   :: N,Nw,unit_
    TYPE(masked_cplx_matrix_type)             :: MM
    TYPE(freq_type)                           :: FREQ

    if(Nw==0) stop 'read raw correl Nw=0'

    CALL delete_correl(CORREL)
    IF(PRESENT(UNIT))THEN
       unit_ = UNIT
    ELSE
       CALL open_safe(unit_,FILEIN,"UNKNOWN","READ",get_unit=.true.)
    ENDIF

    ! READ RAW CORRELATION FROM INPUT FILE
    READ(unit_,*) title

    IF(PRESENT(FILEIN))THEN 
      CALL dump_message(TEXT="################################")
      CALL dump_message(TEXT="### READ "//TRIM(ADJUSTL(title)))
      CALL dump_message(TEXT="### FROM FILE"//TRIM(ADJUSTL(FILEIN))//" ...")
      CALL dump_message(TEXT="################################")
    ENDIF

    READ(unit_,*) N
    READ(unit_,*) stat
    CALL read_raw_masked_cplx_matrix(MM,unit_)
    CALL read_raw_freq(FREQ,unit_)

    ! WE OVERWRITE THE GENERIC FREQUENCY ACCORDINGLY
    if(FREQ%Nw==0) stop 'read raw corel FREQU%Nw=0'
    SELECT CASE(FREQ%title)
      CASE(FERMIONIC,BOSONIC)
       CALL new_correl_from_scratch_ifreq(CORREL,title,N,FREQ%Nw,FREQ%beta,STAT=stat,IMASK=MM%MASK%imat,AB=T) 
      CASE(RETARDED,ADVANCED)
       CALL new_correl_from_scratch_rfreq(CORREL,title,N,FREQ%Nw,FREQ%wmin,FREQ%wmax,FREQ%width,STAT=stat,IMASK=MM%MASK%imat,AB=T) 
    END SELECT

    READ(unit_,*) CORREL%fctn

    ! CLOSE INPUT FILE
    IF(.NOT.PRESENT(UNIT)) CALL close_safe(unit_)
    CALL delete_frequency(FREQ)
    IF(PRESENT(FILEIN))THEN 
      CALL dump_message(TEXT="################################")
      CALL dump_message(TEXT="###      ... DONE ...         ##")
      CALL dump_message(TEXT="################################")
    ENDIF

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE glimpse_correl(CORREL,UNIT)
    TYPE(correl_type), INTENT(IN) :: CORREL
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    CALL write_array(CORREL%MM%mat,'# Check '//TRIM(ADJUSTL(CORREL%title))//'(iw='//c2s(i2c(CORREL%freq%iwprint))//')',UNIT=UNIT)
  return
  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  subroutine dump_chi0(G,beta)
  implicit none
   TYPE(correl_type)  :: G
   real(8)            :: beta
   INTEGER            :: i1,i2,iw,iw1,iw2,iw_,j,Nw
   COMPLEX(DBL)       :: chi0(size(G%freq%vec)),gg(-size(G%freq%vec):size(G%freq%vec))

    Nw=G%Nw
 
    do j=1,size(G%vec,1)
      chi0=0.0 ; gg=0.0 ;
      do iw=1,Nw
        gg(-iw) = conjg(G%vec(j,iw))
        gg( iw) =       G%vec(j,iw)
      enddo
      DO iw=1,Nw
       do iw_=-Nw+iw+1,Nw-iw-1
         chi0(iw)=chi0(iw)+gg(iw_)*gg(iw_+iw)*2.d0/beta
       enddo
      ENDDO
      !!call plotarray( AIMAG(G%freq%vec(1:Nw/2)),  real(chi0(1:Nw/2)),  'real_part_matsubara_Chi0'//toString(j) )
      !!call plotarray( AIMAG(G%freq%vec(1:Nw/2)),  aimag(chi0(1:Nw/2)),   'im_part_matsubara_Chi0'//toString(j) )
    enddo

  return
  end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE write_correl(CORREL,FILEOUT)

    !----------------------------------------!
    ! WRITE CORREL TO FILE IN FORMATTED FORM !
    !----------------------------------------!

    TYPE(correl_type)         :: CORREL
    CHARACTER(LEN=*),OPTIONAL :: FILEOUT
    CHARACTER(LEN=300)        :: fileout_,fmtcorrel,fmtimask
    INTEGER                   :: UNIT 
    INTEGER                   :: i1,i2,iw

    ! CREATE OUTPUT FILE NAME
                         fileout_ = "./ED_out/"//TRIM(ADJUSTL(CORREL%title))//".dat"
    IF(PRESENT(FILEOUT)) fileout_ = FILEOUT

    ! OPEN OUTPUT FILE
    CALL open_safe(UNIT,fileout_,"UNKNOWN","WRITE", get_unit=.true.)

    ! WRITE HEADER

    CALL dump_message(UNIT=UNIT,TEXT="# INDEPENDANT MATRIX ELEMENTS OF "//TRIM(CORREL%stat)//" CORRELATION "//TRIM(ADJUSTL(CORREL%title)))

    SELECT CASE(CORREL%freq%title)
      CASE(FERMIONIC,BOSONIC)
      if(size(CORREL%freq%vec)>1) then
       WRITE(UNIT,'(a,f0.6)') "# ON MATSUBARA AXIS WITH beta = ",pi2/AIMAG(CORREL%freq%vec(2)-CORREL%freq%vec(1))
      endif
      CASE(RETARDED,ADVANCED)
       WRITE(UNIT,'(a,f0.6)') "# ON ("//TRIM(CORREL%freq%title)//") REAL AXIS WITH BROADENING = ",AIMAG(CORREL%freq%vec(1)) 
    END SELECT

    WRITE(UNIT,'(a,I0)') "# N = ",CORREL%N 
    CALL dump_message(UNIT=UNIT,TEXT="# IMASK ")

    if(CORREL%N>1)then
     WRITE(fmtimask,*) '(',CORREL%N-1,'(a,',CORREL%N,'(I2,x)/),a,',CORREL%N,'(I2,x))'
     WRITE(UNIT,fmtimask) ("# ",(CORREL%MM%MASK%imat(i1,i2),i2=1,CORREL%N),i1=1,CORREL%N)
    else
     WRITE(UNIT,'(a,I5)') "# ",CORREL%MM%MASK%imat(1,1)
    endif

    WRITE(UNIT,'(a,I0)') "# Nw = ",CORREL%Nw 
    WRITE(fmtcorrel,*) '(f0.12,',CORREL%MM%MASK%nind,'(2X,f0.12,2X,f0.12,2X))'

    CALL correl2vec(CORREL)

    SELECT CASE(CORREL%freq%title)
      CASE(FERMIONIC,BOSONIC)
       DO iw=1,CORREL%Nw
         WRITE(UNIT,fmtcorrel) AIMAG(CORREL%freq%vec(iw)),CORREL%vec(:,iw)
       ENDDO
       !!call plotarray( AIMAG(CORREL%freq%vec), real(CORREL%vec),  'real_part_matsubara_'//TRIM(ADJUSTL(CORREL%title)))
       !!call plotarray( AIMAG(CORREL%freq%vec), aimag(CORREL%vec), 'im_part_matsubara_'//TRIM(ADJUSTL(CORREL%title)))
      CASE(RETARDED,ADVANCED)
       DO iw=1,CORREL%Nw
         WRITE(UNIT,fmtcorrel)  DBLE(CORREL%freq%vec(iw)),CORREL%vec(:,iw)
       ENDDO
       !!call plotarray( DBLE(CORREL%freq%vec), real(CORREL%vec), 'real_part_retarded_'//TRIM(ADJUSTL(CORREL%title)))
       !!call plotarray( DBLE(CORREL%freq%vec), aimag(CORREL%vec), 'im_part_retarded_'//TRIM(ADJUSTL(CORREL%title)))
    END SELECT 

    CALL close_safe(UNIT)

  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE read_correl(CORREL,FILEIN)

    !--------------------------------------------!
    ! CREATE NEW CORRELATION FROM FORMATTED FILE !
    !--------------------------------------------!

    TYPE(correl_type), INTENT(INOUT) :: CORREL
    CHARACTER(LEN=*),  INTENT(IN)    :: FILEIN
    LOGICAL                          :: is_sym_
    REAL(DBL)                        :: width,beta
    REAL(DBL), ALLOCATABLE           :: tmp_freq(:)
    INTEGER,   ALLOCATABLE           :: IMASK(:,:)
    INTEGER                          :: UNIT 
    INTEGER                          :: mu,nu,iw,N,Nw
    CHARACTER(LEN=300)               :: header,fmtmask,fmtcorrel,title
    CHARACTER(LEN=9)                 :: FTYPE,STAT
    CHARACTER(LEN=4)                 :: junk

    CALL delete_correl(CORREL)

    ! OPEN INPUT FILE
    CALL open_safe(UNIT,FILEIN,"UNKNOWN","READ",get_unit=.true.)

    ! READ HEADER
    READ(UNIT,'(a)') header 
    READ(header(INDEX(header,"OF")+1:),*) title
    IF(INDEX(header,"FERM")/=0) stat = FERMIONIC
    IF(INDEX(header,"BOSO")/=0) stat =   BOSONIC
    READ(UNIT,'(a)') header 
    IF(INDEX(header,"MATS")/=0)THEN
      FTYPE = MATSUBARA
      READ(header(INDEX(header,"=")+1:),*) beta
    ELSE
      IF(INDEX(header,"RETA")/=0)THEN
       FTYPE = RETARDED
      ELSE
       FTYPE = ADVANCED
      ENDIF
    ENDIF
    READ(header(INDEX(header,"=")+1:),*) width
    READ(UNIT,'(a)') header 
    READ(header(INDEX(header,"=")+1:),*) N
    CALL skip_line(UNIT,1)
    WRITE(fmtmask,*) '(',N-1,'(a4,10x,',N,'(I2,x),3x,',N,'(I2,x)/),a4,10x,',N,'(I2,x),3x,',N,'(I2,x))'
    ALLOCATE(IMASK(N,N))
    READ(UNIT,fmtmask) (junk,(IMASK(mu,nu),nu=1,N),mu=1,N)
    is_sym_ = F
    IF(ALL(IMASK==TRANSPOSE(IMASK))) is_sym_ = T
    READ(UNIT,'(a)') header 
    READ(header(INDEX(header,"=")+1:),*) Nw

    if(Nw==0) stop 'read correl Nw=0'

    ! CREATE CORRELATION
    CALL new_correl_from_scratch(CORREL,title,N,Nw,STAT=STAT,IMASK=IMASK,AB=T)

    ! DUMMY TO CREATE INTEGER MAP imaskvec
    CALL masked_cplx_matrix2vec(CORREL%MM)
    ! DUMMY TO ALLOCATE vec
    CALL correl2vec(CORREL)

    ! READ CORRELATION 
    WRITE(fmtcorrel,*) '(f0.12,',CORREL%MM%MASK%nind,'(2X,f0.12,2X,f0.12,2X))'
    ALLOCATE(tmp_freq(Nw))
    DO iw=1,Nw
      READ(UNIT,fmtcorrel) tmp_freq(iw),CORREL%vec(CORREL%MM%MASK%ivec(:),iw)
    ENDDO
    CALL vec2correl(CORREL)

    ! ADJUST FREQUENCIES
    SELECT CASE(FTYPE)
     CASE(MATSUBARA)
       CALL new_freq(CORREL%freq,Nw,beta,STAT)
     CASE(RETARDED)
       CALL new_freq(CORREL%freq,Nw,MINVAL(tmp_freq),MAXVAL(tmp_freq),width,FTYPE)
    END SELECT
    if(messages3) write(log_unit,*) 'read correl new frequ done'

    ! CLOSE INPUT FILE
    CALL close_safe(UNIT)

    DEALLOCATE(tmp_freq,IMASK,STAT=istati)

  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE vec2correl(CORREL,CORRELEXT) 

    !----------------------------------------------------------!
    ! UNCAST VECTOR OF INDPDT ELEMENTS INTO CORRELATION MATRIX !
    !----------------------------------------------------------!

    TYPE(correl_type), INTENT(INOUT)        :: CORREL
    TYPE(correl_type), INTENT(IN), OPTIONAL :: CORRELEXT
    TYPE(masked_cplx_matrix_type)           :: corr,corrext
    INTEGER                                 :: iw

    IF(PRESENT(CORRELEXT))THEN
      CALL new_masked_cplx_matrix(corr,CORREL%MM)
      CALL new_masked_cplx_matrix(corrext,CORRELEXT%MM)
      CALL masked_cplx_matrix2vec(corrext)
      DO iw=1,CORREL%Nw
         corrext%vec = CORRELEXT%vec(:,iw)
         CALL vec2masked_cplx_matrix(corr,corrext)
         CORREL%fctn(:,:,iw) = corr%mat
      ENDDO
      CALL delete_masked_cplx_matrix(corrext)
    ELSE
        CALL new_masked_cplx_matrix(corr,CORREL%MM)
        CALL masked_cplx_matrix2vec(corr)
        DO iw=1,CORREL%Nw
          corr%vec = CORREL%vec(:,iw)
          CALL vec2masked_cplx_matrix(corr)
          CORREL%fctn(:,:,iw) = corr%mat
        ENDDO
    ENDIF

    CALL delete_masked_cplx_matrix(corr)

  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE diff_correl(diff,CORREL1,CORREL2)
    REAL(DBL),         INTENT(INOUT) :: diff
    REAL(8)                          :: ww_,fit_shift_,sumr,asym_,w2,w1
    TYPE(correl_type), INTENT(IN)    :: CORREL1,CORREL2
    integer                          :: iw,window1_,window2_,tot,boundup,ii,jj,nn,b1,b2
    COMPLEX(DBL)                     :: CORREL1_AV(CORREL1%N,CORREL1%N,CORREL1%Nw)
    logical                          :: asym

    diff=zero; tot=0; asym=abs(lambda_sym_fit)>1.d-3
 
    if(abs(fit_shift)<1.d-10)then
     fit_shift_=abs(CORREL1%freq%vec(1))
    else
     fit_shift_=fit_shift
    endif

    if(window_hybrid2>0) then
       window2_ =  min(window_hybrid2,CORREL1%Nw) 
       window1_ =      window_hybrid
       ww_      = dble(window_weight)
    else
       window1_ = window_hybrid
       window2_ = CORREL1%Nw
       ww_      = 1.d0
    endif

    if(fit_nw<0) then
     boundup=CORREL1%Nw
    else
     boundup=min(fit_nw,CORREL1%Nw)
    endif
    nn = size(CORREL1%fctn,1)

   !--------------------------------------------------------------------!
    if(asym) call average_correlations(CORREL1,CORREL1_AV,average_G>=0,MASK_AVERAGE_SIGMA,boundup)
   !--------------------------------------------------------------------!

   DO iw=1,boundup
     tot  = tot  + 1
     if(iw>window1_.and.iw<window2_)then
        call define_sumr
        diff = diff + ww_ * sumr 
     else
        call define_sumr
        diff = diff + sumr
     endif
   ENDDO

   diff = ABS(diff) / dble( tot * CORREL1%N**2)

  RETURN
 
  contains

 !---------------------------------------!
  subroutine define_sumr
        sumr=0.
        do ii=1,nn
         call get_b1
         do jj=b1,b2
          call get_asym
                     w1=1.d0
          if(ii/=jj) w1=w2
          sumr = sumr + (asym_+ABS(CORREL2%fctn(ii,jj,iw)-CORREL1%fctn(ii,jj,iw)))**weight_expo*w1
         enddo
        enddo
        sumr = sumr/((ABS(CORREL1%freq%vec(iw)-CORREL1%freq%vec(1))+fit_shift_)**fit_weight_power)
  end subroutine
 !---------------------------------------!
  subroutine get_asym
          if(asym) then
           asym_=abs(CORREL1_AV(ii,jj,iw)-CORREL1%fctn(ii,jj,iw))*lambda_sym_fit
          else
           asym_=0.d0
          endif
  end subroutine 
 !---------------------------------------!
  subroutine get_b1
  b1=1
  b2=nn
  w2=1.d0
  if(fmos)then
   b1=ii
   b2=ii
   w2=1.d0
   return
  endif
  if(.not.fast_fit) return
#ifndef _complex
    b1=ii
    b2=nn
    w2=2.d0
    if(.not.superconducting_state)then
     if(ii<=nn/2)then; b2=nn/2; endif
    endif
#else
    b1=1
    b2=nn
    w2=1.d0
    if(.not.superconducting_state)then
     if(ii<=nn/2)then; b2=nn/2  ; endif
     if(ii> nn/2)then; b1=nn/2+1; endif
    endif
#endif
  end subroutine
 !---------------------------------------!

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************


  SUBROUTINE pad_correl(CORRELOUT,CORRELIN)
    TYPE(correl_type), INTENT(INOUT)        :: CORRELOUT
    TYPE(correl_type), INTENT(IN), OPTIONAL :: CORRELIN
    INTEGER                                 :: i1,i2,iw,rankin,rankout,ii1,ii2

    IF(PRESENT(CORRELIN))THEN
      ! FIRST COPY VECTOR OF INDEPENDANT ELEMENTS
      DO rankin=1,CORRELIN%MM%MASK%nind
       rankout = find_rank(CORRELIN%MM%MASK%ivec(rankin),CORRELOUT%MM%MASK%ivec)
       IF(rankout/=0) CORRELOUT%vec(rankout,:) = CORRELIN%vec(rankin,:)
      ENDDO
      ! THEN EXPAND IT TO FULL MATRIX
      CALL vec2correl(CORRELOUT)
    ENDIF

    DO i1=1,CORRELOUT%N; DO i2=1,CORRELOUT%N; IF(CORRELOUT%MM%MASK%mat(i1,i2))THEN
      do ii1=1,size(CORRELOUT%MM%MASK%imat,1)
       do ii2=1,size(CORRELOUT%MM%MASK%imat,2)
        do iw=1,CORRELOUT%Nw
          if(CORRELOUT%MM%MASK%imat(ii1,ii2)==CORRELOUT%MM%MASK%imat(i1,i2)) then
            CORRELOUT%fctn(ii1,ii2,iw) = CORRELOUT%fctn(i1,i2,iw)
          endif
        ENDDO
       ENDDO 
      ENDDO 
    ENDIF; ENDDO; ENDDO

  END SUBROUTINE


!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE slice_correl(CORRELOUT,CORRELIN,rbounds,cbounds)
    TYPE(correl_type), INTENT(INOUT) :: CORRELOUT
    TYPE(correl_type), INTENT(IN)    :: CORRELIN
    INTEGER,           INTENT(IN)    :: rbounds(2),cbounds(2)
    INTEGER                          :: n1slice,n2slice
    CHARACTER(LEN=100)               :: cslice,title
    REAL(DBL)                        :: beta,width,wmax,wmin

    if(CORRELIN%Nw==0) stop 'slice correl Nw=0'

    CALL delete_correl(CORRELOUT)
    ! FIND/CHECK DIMENSIONS
    n1slice = rbounds(2)-rbounds(1)+1
    n2slice = cbounds(2)-cbounds(1)+1
    IF(n1slice/=n2slice) STOP "ERROR IN slice_correl: SQUARE   SLICE      REQUIRED!"
    IF(n1slice<=0)       STOP "ERROR IN slice_correl: NON-NULL DIMENSIONS REQUIRED!"
    ! MAKE TITLE
    title = TRIM(ADJUSTL(CORRELIN%title))
    IF( rbounds(1)/=LBOUND(CORRELIN%fctn,1).OR.rbounds(2)/=UBOUND(CORRELIN%fctn,1).OR. &
      & cbounds(1)/=LBOUND(CORRELIN%fctn,2).OR.cbounds(2)/=UBOUND(CORRELIN%fctn,2))THEN
      cslice = c2s(i2c(rbounds(1)))//"_"//c2s(i2c(rbounds(2)))//"_"//c2s(i2c(cbounds(1)))//"_"//c2s(i2c(cbounds(2)))
      title  = TRIM(ADJUSTL(title))//"_"//TRIM(ADJUSTL(cslice))
    ENDIF

    ! SLICE
    SELECT CASE(CORRELIN%freq%title)
     CASE(FERMIONIC,BOSONIC)
       beta  = pi2 / AIMAG(CORRELIN%freq%vec(2)-CORRELIN%freq%vec(1))
       CALL new_correl_from_scratch_ifreq(CORRELOUT,title,CORRELIN%N,CORRELIN%Nw,beta, &
                                        & STAT=CORRELIN%STAT,IMASK=CORRELIN%MM%MASK%imat,AB=T) 
     CASE(RETARDED,ADVANCED)
       wmin  = CORRELIN%freq%wmin 
       wmax  = CORRELIN%freq%wmax
       width = AIMAG(CORRELIN%freq%vec(1))
       CALL new_correl_from_scratch_rfreq(CORRELOUT,title,CORRELIN%N,CORRELIN%Nw,wmin,wmax,width, &
                                        & STAT=CORRELIN%STAT,IMASK=CORRELIN%MM%MASK%imat,AB=T) 
    END SELECT 

    CALL slice_masked_cplx_matrix(CORRELOUT%MM,CORRELIN%MM,rbounds,cbounds)

    IF(ASSOCIATED(CORRELOUT%fctn))   DEALLOCATE(CORRELOUT%fctn,STAT=istati)
    IF(ASSOCIATED(CORRELOUT%MM%mat)) DEALLOCATE(CORRELOUT%MM%mat,STAT=istati)
    IF(ASSOCIATED(CORRELOUT%vec))    DEALLOCATE(CORRELOUT%vec,STAT=istati)

    CORRELOUT%fctn   => CORRELIN%fctn(rbounds(1):rbounds(2),cbounds(1):cbounds(2),:)
    CORRELOUT%MM%mat => CORRELOUT%fctn(:,:,CORRELOUT%freq%iwprint)
    CORRELOUT%vec    => CORRELIN%vec

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE average_correlations(G,Gm1,offdiag_also,mask,boundup)
  implicit none
  TYPE(correl_type) :: G
  INTEGER           :: iw,siz,i,j,miw,si
  COMPLEX(DBL)      :: Gm1(:,:,:)
  REAL(DBL)         :: rvec(size(Gm1,3))
  LOGICAL           :: offdiag_also
  integer,intent(in):: mask(:,:)
  integer           :: Nw
  integer,optional  :: boundup

  if(present(boundup)) then
    Nw=boundup
  else
    Nw=size(Gm1,3)
  endif

   siz=size(Gm1,1);si=size(Gm1,1)/2

   Gm1(:,:,1:Nw)=G%fctn(:,:,1:Nw)

   !==========================================================!
    SELECT CASE(G%freq%title)
      CASE(FERMIONIC)
       do iw=1,Nw
        Gm1(si+1:siz,si+1:siz,iw) = -CONJG(G%fctn(si+1:siz,si+1:siz,iw))
       enddo
      CASE(BOSONIC)
       do iw=1,Nw
        Gm1(si+1:siz,si+1:siz,iw) =  CONJG(G%fctn(si+1:siz,si+1:siz,iw))
       enddo
      CASE(RETARDED,ADVANCED)
        rvec = real(G%freq%vec)
        Gm1(si+1:siz,si+1:siz,:) = G%fctn(si+1:siz,si+1:siz,:)
        do i=si+1,siz ;do j=si+1,siz;
         call mirror_array(rvec,Gm1(i,j,:))
          SELECT CASE(G%stat)
           CASE(FERMIONIC)
            Gm1(i,j,:)=-CONJG(Gm1(i,j,:))
           CASE(BOSONIC)
            Gm1(i,j,:)= CONJG(Gm1(i,j,:))
           CASE DEFAULT
            stop 'average correlations'
          END SELECT
        enddo ; enddo ;
    END SELECT
    !==========================================================!
    do iw=1,Nw
       call average_matrix(Gm1(:,:,iw),mask,offdiag_also)
    enddo
    !==========================================================!
    SELECT CASE(G%freq%title)
      CASE(FERMIONIC)
       do iw=1,Nw
        Gm1(si+1:siz,si+1:siz,iw) = -CONJG(Gm1(si+1:siz,si+1:siz,iw))
       enddo
      CASE(BOSONIC)
       do iw=1,Nw
        Gm1(si+1:siz,si+1:siz,iw) =  CONJG(Gm1(si+1:siz,si+1:siz,iw))
       enddo
      CASE(RETARDED,ADVANCED)
        rvec = real(G%freq%vec)
        do i=si+1,siz ;do j=si+1,siz;
          call mirror_array(rvec,Gm1(i,j,:))
          SELECT CASE(G%stat)
           CASE(FERMIONIC)
            Gm1(i,j,:)=-CONJG(Gm1(i,j,:))
           CASE(BOSONIC)
            Gm1(i,j,:)= CONJG(Gm1(i,j,:))
           CASE DEFAULT
            stop 'average correlations'
          END SELECT
        enddo; enddo;
    END SELECT
    !==========================================================!

  RETURN
  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE transform_correl(CORRELOUT,CORRELIN)
  implicit none

  !---------------------------------------!
  ! CORRELOUT(z) = (-1)^F * CORRELIN(-z*) !
  !---------------------------------------!

    TYPE(correl_type)            :: CORRELOUT
    TYPE(correl_type),INTENT(IN) :: CORRELIN
    INTEGER                      :: iw,i,j,miw
    COMPLEX(DBL)                 :: vec(size(CORRELOUT%vec,1),size(CORRELOUT%vec,2))
    COMPLEX(DBL)                 :: fctn(size(CORRELOUT%fctn,1),size(CORRELOUT%fctn,2),size(CORRELOUT%fctn,3))
    REAL(DBL)                    :: rvec(size(CORRELOUT%freq%vec))
    COMPLEX(DBL)                 :: vecin(size(CORRELIN%vec,1),size(CORRELIN%vec,2))
    COMPLEX(DBL)                 :: fctnin(size(CORRELIN%fctn,1),size(CORRELIN%fctn,2),size(CORRELIN%fctn,3))
    REAL(DBL)                    :: rvecin(size(CORRELIN%freq%vec))

    if(CORRELIN%Nw==0)                            STOP "transform correl Nw=0"
    IF(CORRELIN%freq%title/=CORRELOUT%freq%title) STOP "ERROR IN transform_correl: INCONSISTENT FREQUENCIES!"
    IF(CORRELIN%stat      /=CORRELOUT%stat)       STOP "ERROR IN transform_correl: INCONSISTENT STATISTICS!"
    IF(CORRELIN%Nw        /=CORRELOUT%Nw)         STOP "ERROR IN transform_correl: INCONSISTENT DIMENSIONS!"

    rvecin = real(CORRELIN%freq%vec)
    vecin  = CORRELIN%vec
    fctnin = CORRELIN%fctn

    !==========================================================!
    SELECT CASE(CORRELIN%freq%title)
      CASE(FERMIONIC)
        CORRELOUT%vec   = -CONJG(vecin )
        CORRELOUT%fctn  = -CONJG(fctnin)
      CASE(BOSONIC)
        CORRELOUT%vec   =  CONJG(vecin )
        CORRELOUT%fctn  =  CONJG(fctnin)
      CASE(RETARDED,ADVANCED)
        !-------------------------------------------------------!
                        call mirror_func_
        !-------------------------------------------------------!
      CASE DEFAULT
        STOP 'ERROR transform_correl case does not exist, critical stop'
     END SELECT
    !==========================================================!

    contains

    !------------------------!
    !------------------------!
    !------------------------!
    !------------------------!

    subroutine dump_graphs
    implicit none
       if(rank/=0) return
       call PGSUBP(1,1)
       do j=1,size(CORRELIN%vec,1)
         !!call plotarray( real(CORRELIN%freq%vec),  real(vecin(j,:)),  & 
         !                real(CORRELIN%freq%vec),  real(CORRELOUT%vec(j,:)) ,'re transfo correl'//toString(j))
         !!call plotarray( real(CORRELIN%freq%vec), aimag(vecin(j,:)),  &
          !             & real(CORRELIN%freq%vec), aimag(CORRELOUT%vec(j,:)) ,'im transfo correl'//toString(j))
       enddo
    end subroutine

    !------------------------!
    !------------------------!
    !------------------------!
    !------------------------!

   subroutine mirror_func_
   implicit none
 
          DO iw=1,CORRELOUT%Nw
           SELECT CASE(CORRELOUT%stat)
              CASE(FERMIONIC)
                CORRELOUT%vec   (:,iw) = - CONJG(vecin   (:,iw))
                CORRELOUT%fctn(:,:,iw) = - CONJG(fctnin(:,:,iw))
              CASE(BOSONIC)
                CORRELOUT%vec   (:,iw) =   CONJG(vecin   (:,iw))
                CORRELOUT%fctn(:,:,iw) =   CONJG(fctnin(:,:,iw))
           END SELECT
          ENDDO

          !-----------!
          ! w --> -w  !
          !-----------!

           rvec = real(CORRELOUT%freq%vec)
           vec  =      CORRELOUT%vec
           fctn =      CORRELOUT%fctn

           do j=1,size(CORRELOUT%vec,1)
              call mirror_array(rvec,vec(j,:))
           enddo
           do i=1,size(fctn,1); do j=1,size(fctn,2)
              call mirror_array(rvec,fctn(i,j,:))
           enddo;enddo

           CORRELOUT%vec  = vec
           CORRELOUT%fctn = fctn

   end subroutine

    !------------------------!
    !------------------------!
    !------------------------!
    !------------------------!

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

END MODULE 
