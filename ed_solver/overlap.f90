MODULE overlap_module


  USE density_matrix


  IMPLICIT NONE

  private

  REAL(DBL), PARAMETER, PRIVATE   :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,   PARAMETER, PRIVATE   :: F=.FALSE.,T=.TRUE.

  !---------------------------------------------------------------------!
  ! COMPUTE THE WEIGHT OF THE GROUND STATE ON CERTAIN REFERENCE VECTORS !
  !---------------------------------------------------------------------!


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
!**************************************************************************
!**************************************************************************


  SUBROUTINE compute_overlap(AIM,beta,GS,vec_list)

    !-----------------------------! 
    ! COMPUTE STATIC CORRELATIONS !
    !-----------------------------!

    TYPE(AIM_type),             INTENT(IN)    :: AIM
    REAL(DBL),                  INTENT(IN)    :: beta
    TYPE(eigensectorlist_type), INTENT(IN)    :: GS
    TYPE(readable_vec_type),    INTENT(INOUT) :: vec_list(:)
    TYPE(eigen_type), POINTER                 :: eigen => NULL()
    INTEGER                                   :: isector,ieigen,iGS,nGS,ivec,iket,AIMrank,IMPstate,BATHstate
    REAL(DBL)                                 :: boltz,Zpart,E0
    COMPLEX(DBL)                              :: coeff
    COMPLEX(DBL), ALLOCATABLE                 :: overlap(:,:),roverlap(:,:)
    REAL(DBL),    ALLOCATABLE                 :: phase_(:)
    INTEGER                                   :: nvec,nkets,kets(100) 
    
    CALL dump_message(TEXT="######################################")
    CALL dump_message(TEXT="### OVERLAP WITH REFERENCE VECTORS ###")
    CALL dump_message(TEXT="######################################")

    !--------------------------------------------!
    ! FIRST WE FIND ALL DISTINCT IMPURITY STATES !
    !--------------------------------------------!

    nvec   = SIZE(vec_list)
    nkets  =  0 
    kets   = -1

    DO ivec=1,nvec
      DO iket=1,vec_list(ivec)%nket 
        IF(find_rank(vec_list(ivec)%state(iket),kets(1:nkets))==0)THEN
          nkets       = nkets + 1
          kets(nkets) = vec_list(ivec)%state(iket)
        ENDIF
      ENDDO
    ENDDO

    !----------------------------------------------------!    
    ! FINALLY COMPUTE THE OVERLAP WITH REFERENCE VECTORS !
    !----------------------------------------------------!

    Zpart   = partition(beta,GS)
    E0      = GSenergy(GS)
    nGS     = nlowest(GS)

    ALLOCATE(overlap(nvec,nGS))
    overlap=zero 
 iGS=0

    DO isector=1,GS%nsector
     DO ieigen=1,GS%es(isector)%lowest%neigen
        iGS = iGS+1
        eigen => GS%es(isector)%lowest%eigen(ieigen)
        boltz = DEXPc(-beta*(eigen%val-E0))/Zpart 
       !---------------------------------------------------------!
         DO AIMrank=1,dimen_func(GS%es(isector)%sector)
          CALL AIM2IMPBATHstate(IMPstate,BATHstate,state_func(AIMrank,GS%es(isector)%sector),AIM%Nc,AIM%Nb) 
          IF(find_rank(IMPstate,kets)/=0)THEN
            DO ivec=1,nvec
              DO iket=1,vec_list(ivec)%nket
               IF(vec_list(ivec)%state(iket)==IMPstate)THEN
                coeff             = vec_list(ivec)%coeff(iket)
                overlap(ivec,iGS) = overlap(ivec,iGS) + conj(coeff)*eigen%vec%rc(AIMrank)*boltz
               ENDIF
              ENDDO
            ENDDO
          ENDIF
         ENDDO
       !---------------------------------------------------------!
     ENDDO
    ENDDO
   
    !-----------------------------------------------------------------------! 
    ! RENORMALIZE WITH OVERALL PHASE FACTOR SUCH THAT FIRST OVERLAP IS REAL !
    !-----------------------------------------------------------------------!

    IF(nvec>=1)THEN
      ALLOCATE(phase_(nGS),roverlap(nvec,nGS))
      DO iGS=1,nGS
        phase_(iGS)     = ATAN2(AIMAG(overlap(1,iGS)),DBLE(overlap(1,iGS)))
        roverlap(:,iGS) = overlap(:,iGS) * MPLX(-phase_(iGS))
      ENDDO
    ENDIF
    
    DO ivec=1,nvec
      write(log_unit,*) '----------------------------------------------------'
      write(log_unit,*) 'projection on state         : ', ivec,TRIM(ADJUSTL(vec_list(ivec)%title))
      write(log_unit,*) 'amplitude in each sectors   : '
      write(log_unit,'(1000f8.5)') (abs(overlap(ivec,iGS)),iGS=1,nGS)
      write(log_unit,*) 'phase in each sectors (pi)  : '
      write(log_unit,'(1000f8.5)') (ATAN2(AIMAG(roverlap(ivec,iGS)),DBLE(roverlap(ivec,iGS)))/pi,iGS=1,nGS)
      write(log_unit,*) '----------------------------------------------------'
    ENDDO

    DEALLOCATE(overlap)
 if(allocated(phase_)) DEALLOCATE(phase_,roverlap)

  END SUBROUTINE 


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
!**************************************************************************
!**************************************************************************

END MODULE 
