MODULE AIM_class

  USE bath_class_hybrid

  IMPLICIT NONE

  !-------------!
  ! GENERIC AIM !
  !-------------!

  TYPE AIM_type
    TYPE(impurity_type), POINTER :: impurity => NULL() 
    TYPE(bath_type),     POINTER :: bath     => NULL()
    INTEGER :: Ns      = 0 ! TOTAL # OF SITES
    INTEGER :: norbs   = 0 ! TOTAL # OF 1-PARTICLE ORBITALS
    INTEGER :: nstates = 0 ! TOTAL # OF STATES
    INTEGER :: Nc      = 0 ! # OF IMPURITY SITES 
    INTEGER :: Nb      = 0 ! # OF BATH     SITES
    !
    ! CONVERSION ARRAY FOR ORBITALS
    !
    INTEGER, POINTER ::  IMPiorb(:,:) => NULL()
    INTEGER, POINTER :: BATHiorb(:,:) => NULL()
  END TYPE 


CONTAINS 



!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE update_AIM_pointer(AIM,impurity,bath)
    TYPE(AIM_type),      INTENT(INOUT)       ::  AIM
    TYPE(impurity_type), INTENT(IN), TARGET  ::  impurity
    TYPE(bath_type),     INTENT(IN), TARGET  ::  bath
    AIM%impurity => impurity
    AIM%bath     => bath
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

  SUBROUTINE new_AIM(AIM,impurity,bath)
    TYPE(AIM_type),      INTENT(INOUT)       ::  AIM
    TYPE(impurity_type), INTENT(IN), TARGET  ::  impurity
    TYPE(bath_type),     INTENT(IN), TARGET  ::  bath

    CALL delete_AIM(AIM)
    CALL update_AIM_pointer(AIM,impurity,bath)

    AIM%Nc      = impurity%Nc
    AIM%Nb      =     bath%Nb
    AIM%Ns      = AIM%Nc + AIM%Nb
    AIM%norbs   = impurity%norbs   + bath%norbs
    AIM%nstates = impurity%nstates * bath%nstates
    
    ! ORDERING OF ORBITALS WITH INCREASING RANK 
    ! = |BATH(site,up),IMP(site,up)> |BATH(site,down),IMP(site,down)>
    
    ! SPIN UP
    ALLOCATE(AIM%BATHiorb(AIM%Nb,2)) 
    CALL ramp(AIM%BATHiorb(:,1),AIM%Nb)
    
    ALLOCATE(AIM%IMPiorb(AIM%Nc,2)) 
    CALL ramp(AIM%IMPiorb(:,1),AIM%Nc)
    AIM%IMPiorb(:,1)  = AIM%IMPiorb(:,1) + bath%Nb
    
    ! SPIN DOWN
    AIM%IMPiorb( :,2) =  AIM%IMPiorb(:,1) + AIM%Ns
    AIM%BATHiorb(:,2) = AIM%BATHiorb(:,1) + AIM%Ns
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

  SUBROUTINE delete_AIM(AIM)
    TYPE(AIM_type), INTENT(INOUT) :: AIM
    if(associated(AIM%BATHiorb)) deallocate(AIM%BATHiorb,STAT=istati)
    if(associated(AIM%IMPiorb )) deallocate(AIM%IMPiorb,STAT=istati)
    AIM%impurity => NULL()
    AIM%bath     => NULL()
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

  SUBROUTINE IMPBATH2AIMstate(AIMstate,IMPstate,BATHstate,Nc,Nb) 

    ! RETURNS THE AIM STATE 
    ! GIVEN   ITS IMPURITY / BATH COMPONENTS
    
    ! REMEMBER ORDERING OF ORBITALS 
    ! | AIMstate> = |BATH(site,up)>| IMP(site,up)> |BATH(site,down)>|IMP(site,down)>
    ! | IMPstate> = | IMP(site,up)>| IMP(site,down)>
    ! |BATHstate> = |BATH(site,up)>|BATH(site,down)>
    
    INTEGER, INTENT(OUT) :: AIMstate  ! AIM      STATE
    INTEGER, INTENT(IN)  :: IMPstate  ! IMPURITY STATE
    INTEGER, INTENT(IN)  :: BATHstate ! BATH     STATE
    INTEGER, INTENT(IN)  :: Nc,Nb     ! IMPURITY / BATH SIZES
    INTEGER              :: IMPstateup, IMPstatedo
    INTEGER              :: BATHstateup,BATHstatedo
    INTEGER              :: AIMstateup, AIMstatedo
    INTEGER              :: nBATHstates_sz,nIMPstates_sz,nAIMstates_sz,ii

    nAIMstates_sz  = 2**(Nc+Nb)
    nBATHstates_sz = 2**    Nb
    nIMPstates_sz  = 2** Nc

    ! RESOLVE up/down COMPONENTS OF BATHstate
    BATHstateup = MOD(BATHstate,nBATHstates_sz)
    BATHstatedo = (BATHstate - BATHstateup) / nBATHstates_sz

    ! RESOLVE up/down COMPONENTS OF  IMPstate
    IMPstateup  = MOD(IMPstate,nIMPstates_sz)
    IMPstatedo  = ( IMPstate -  IMPstateup) / nIMPstates_sz

    ! BUILD AIMstate
    AIMstateup  = BATHstateup +  nBATHstates_sz *  IMPstateup 
    AIMstatedo  = BATHstatedo +  nBATHstates_sz *  IMPstatedo 
    AIMstate    =  AIMstateup +   nAIMstates_sz *  AIMstatedo 
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
 
  SUBROUTINE AIM2IMPBATHstate(IMPstate,BATHstate,AIMstate,Nc,Nb) 
    ! RETURNS THE IMPURITY / BATH COMPONENTS 
    ! OF A GIVEN AIM STATE
    
    ! REMEMBER ORDERING OF ORBITALS 
    ! | AIMstate> = |BATH(site,up)>| IMP(site,up)> |BATH(site,down)>|IMP(site,down)>
    ! | IMPstate> = | IMP(site,up)>| IMP(site,down)>
    ! |BATHstate> = |BATH(site,up)>|BATH(site,down)>
    
    INTEGER, INTENT(OUT) :: IMPstate  ! IMPURITY STATE
    INTEGER, INTENT(OUT) :: BATHstate ! BATH     STATE
    INTEGER, INTENT(IN)  :: AIMstate  ! AIM      STATE
    INTEGER, INTENT(IN)  :: Nc,Nb     ! IMPURITY / BATH SIZES
    INTEGER              :: IMPstateup, IMPstatedo
    INTEGER              :: BATHstateup,BATHstatedo
    INTEGER              :: AIMstateup, AIMstatedo
    INTEGER              :: nBATHstates_sz,nIMPstates_sz,nAIMstates_sz

    nAIMstates_sz  = 2**(Nc+Nb)
    nBATHstates_sz = 2**    Nb
    nIMPstates_sz  = 2** Nc

    ! RESOLVE up/down COMPONENTS OF AIMstate
    AIMstateup  = MOD(  AIMstate, nAIMstates_sz)
    AIMstatedo  = (  AIMstate -  AIMstateup) /  nAIMstates_sz

    ! RESOLVE BATH/IMP COMPONENTS OF AIMstateup
    BATHstateup = MOD(AIMstateup,nBATHstates_sz)
    IMPstateup  = (AIMstateup - BATHstateup) / nBATHstates_sz

    ! RESOLVE BATH/IMP COMPONENTS OF AIMstatedo
    BATHstatedo = MOD(AIMstatedo,nBATHstates_sz)
    IMPstatedo  = (AIMstatedo - BATHstatedo) / nBATHstates_sz

    ! BUILD BATHstate
    BATHstate   = BATHstateup + nBATHstates_sz * BATHstatedo 

    ! BUILD  IMPstate
    IMPstate    =  IMPstateup +  nIMPstates_sz *  IMPstatedo 
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

END MODULE

