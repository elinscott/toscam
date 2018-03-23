MODULE apply_NS

  use apply_C 

  IMPLICIT NONE

  !-----------------------------------------------------!
  ! MODULE TO APPLY NUMBER/SPIN OPERATORS ON EIGENSTATE !
  !-----------------------------------------------------!

  REAL(DBL),PARAMETER ,PRIVATE :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,PARAMETER   ,PRIVATE :: F=.FALSE.,T=.TRUE.
  INTEGER,ALLOCATABLE ,PRIVATE :: IMPiorbupdo(:,:),IMPiorbsz(:,:)
  INTEGER             ,PRIVATE :: Ns,Ns2

  logical,parameter,private :: force_reset_list=.false.

CONTAINS



!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE init_apply_NS(AIM)
    TYPE(AIM_type), INTENT(IN) :: AIM
    ALLOCATE(IMPiorbupdo(SIZE(AIM%IMPiorb,1),SIZE(AIM%IMPiorb,2)))
    IMPiorbupdo(:,1) = AIM%IMPiorb(:,1)
    IMPiorbupdo(:,2) = AIM%IMPiorb(:,1) ! WARNING!
    ALLOCATE(IMPiorbsz  (SIZE(AIM%IMPiorb,1),SIZE(AIM%IMPiorb,2)))
    IMPiorbsz  (:,1) = AIM%IMPiorb(:,1)
    IMPiorbsz  (:,2) = AIM%IMPiorb(:,2)
    Ns  = AIM%Ns
    Ns2 = Ns*2
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
 
  SUBROUTINE Sz_sector(Csec,pm,sector_in) 
    TYPE(sector_type), INTENT(INOUT) :: Csec
    TYPE(sector_type), INTENT(IN)    :: sector_in
    CHARACTER(LEN=1),  INTENT(IN)    :: pm ! not used!
    CALL delete_sector(Csec)
    CALL new_sector(Csec,sector_in)
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE N_sector(Csec,pm,sector_in)
    TYPE(sector_type), INTENT(INOUT) :: Csec
    TYPE(sector_type), INTENT(IN)    :: sector_in
    CHARACTER(LEN=1),  INTENT(IN)    :: pm ! not used !
    CALL delete_sector(Csec)
    CALL new_sector(Csec,sector_in)
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE S_sector(Csec,pm,sector_in) 
    TYPE(sector_type), INTENT(INOUT) :: Csec
    TYPE(sector_type), INTENT(IN)    :: sector_in
    CHARACTER(LEN=1),  INTENT(IN)    :: pm
    INTEGER :: nup,ndo
    CALL delete_sector(Csec)
    IF(ASSOCIATED(sector_in%sz))THEN ! Sz-SECTOR
      ALLOCATE(Csec%sz)
      IF(pm=='+') CALL new_fermion_sector(Csec%sz,npart_func(sector_in)+2,Ns2,SZ=T)
      IF(pm=='-') CALL new_fermion_sector(Csec%sz,npart_func(sector_in)-2,Ns2,SZ=T)
    ELSE IF(ASSOCIATED(sector_in%updo))THEN ! (nup,ndo) SECTOR
      nup = npart_func(sector_in,1)
      ndo = npart_func(sector_in,2)
      ALLOCATE(Csec%updo)
      IF(pm=='+') CALL new_fermion_sector2(Csec%updo,nup+1,ndo-1,Ns)
      IF(pm=='-') CALL new_fermion_sector2(Csec%updo,nup-1,ndo+1,Ns)
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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE apply_Sz(Ces,pm,MASK,es,esrank) 

    !$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ COMPUTE Sz[site]|0> $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$

    TYPE(eigensector_type), INTENT(INOUT) :: Ces
    CHARACTER(LEN=1),       INTENT(IN)    :: pm   ! DUMMY
    LOGICAL,                INTENT(IN)    :: MASK(:)
    TYPE(eigensector_type), INTENT(IN)    :: es
    INTEGER,                INTENT(IN)    :: esrank
    TYPE(eigen_type)                      :: eigen_out
    TYPE(eigen_type), POINTER             :: eigen_in => NULL()
    TYPE(fermion_ket_type)                :: ket_in,ket_in_up,ket_in_do
    INTEGER                               :: istate,iup,ido,site,nup,ndo
    REAL(DBL)                             :: sz
    
    eigen_in => es%lowest%eigen(esrank)
    if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
    CALL new_eigen(eigen_out,dimen_func(Ces%sector))

    eigen_out%val       = eigen_in%val
    eigen_out%converged = eigen_in%converged

  !-----------------------------------------------!
  ! WE ONLY NEED TO CONSIDER A SUBSET OF ORBITALS !
  !-----------------------------------------------!

    DO site=1,SIZE(MASK)
 IF(MASK(site))THEN

      !----------------------------------------------------------!
      ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR !
      !----------------------------------------------------------!

      eigen_out%rank   = site 
      eigen_out%vec%rc = zero

      !-------------------------------------------------------------!
      ! THEN PARSE THE INPUT SECTOR TO APPLY RELEVANT SPIN OPERATOR !
      !-------------------------------------------------------------!

      IF(ASSOCIATED(es%sector%sz))THEN ! Sz-SECTOR
      DO istate=1,es%sector%sz%dimen
 IF(eigen_in%vec%rc(istate)/=zero)THEN
        CALL new_ket(ket_in,es%sector%sz%state(istate),Ns2)
        sz = half * ( ni__(IMPiorbsz(site,1),ket_in) - ( 1 - ni__(IMPiorbsz(site,2),ket_in) ) ) ! NAMBU 
        IF(sz/=zero) eigen_out%vec%rc(istate) = eigen_out%vec%rc(istate) + sz * eigen_in%vec%rc(istate) 
      ENDIF
 ENDDO
      ELSE IF(ASSOCIATED(es%sector%updo))THEN ! (nup,ndo) SECTOR
      DO ido=1,es%sector%updo%down%dimen
        CALL new_ket(ket_in_do,es%sector%updo%down%state(ido),Ns)
        ndo = ni__(IMPiorbupdo(site,2),ket_in_do)
        DO iup=1,es%sector%updo%up%dimen
 istate = rankupdo(iup,ido,es%sector%updo)
 IF(eigen_in%vec%rc(istate)/=zero)THEN 
          CALL new_ket(ket_in_up,es%sector%updo%up%state(iup),Ns)
          nup = ni__(IMPiorbupdo(site,1),ket_in_up)
          IF(nup/=ndo) eigen_out%vec%rc(istate) = eigen_out%vec%rc(istate) + half * ( nup - ndo ) * eigen_in%vec%rc(istate)
        ENDIF
 ENDDO
      ENDDO
      ENDIF

      !-------------------------------------------!
      ! THEN WE UPDATE THE LIST OF OUTPUT VECTORS !
      !-------------------------------------------!

      CALL add_eigen(eigen_out,Ces%lowest)

    ENDIF
 ENDDO

    CALL delete_eigen(eigen_out)

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

  SUBROUTINE apply_S(Ces,pm,MASK,es,esrank) 

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ COMPUTE S^-[site]|0>, S^+[site]|0> $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    TYPE(eigensector_type), INTENT(INOUT) :: Ces
    CHARACTER(LEN=1),       INTENT(IN)    :: pm
    LOGICAL,                INTENT(IN)    :: MASK(:)
    TYPE(eigensector_type), INTENT(IN)    :: es
    INTEGER,                INTENT(IN)    :: esrank
    TYPE(eigen_type)                      :: eigen_out
    TYPE(eigen_type), POINTER             :: eigen_in => NULL()
    TYPE(fermion_ket_type)                :: ket_in,ket_out,ket_in_up,ket_in_do,ket_out_up,ket_out_do
    INTEGER                               :: istate,jstate,iup,ido,jup,jdo,site,IMPstate,BATHstate
    
    eigen_in => es%lowest%eigen(esrank)
    if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
    CALL new_eigen(eigen_out,dimen_func(Ces%sector))

    eigen_out%val       = eigen_in%val
    eigen_out%converged = eigen_in%converged

    ! WE ONLY NEED TO CONSIDER A SUBSET OF ORBITALS

    DO site=1,SIZE(MASK)
 IF(MASK(site))THEN

      ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR

      eigen_out%rank   = site
      eigen_out%vec%rc = zero

      ! THEN PARSE THE INPUT SECTOR TO APPLY RELEVANT SPIN OPERATOR

      IF(ASSOCIATED(es%sector%sz))THEN ! Sz-SECTOR
      DO istate=1,es%sector%sz%dimen
 IF(eigen_in%vec%rc(istate)/=zero)THEN
        CALL new_ket(ket_in,es%sector%sz%state(istate),es%sector%sz%norbs)
        IF(pm=='+') CALL  create_pair(ket_out,IMPiorbsz(site,1),IMPiorbsz(site,2),ket_in)
        IF(pm=='-') CALL destroy_pair(ket_out,IMPiorbsz(site,2),IMPiorbsz(site,1),ket_in)
        IF(.NOT.ket_out%is_nil)THEN
          jstate = Ces%sector%sz%rank(ket_out%state)
          eigen_out%vec%rc(jstate) = eigen_out%vec%rc(jstate) + ket_out%fermion_sign * eigen_in%vec%rc(istate) 
        ENDIF
      ENDIF
 ENDDO
      ELSE IF(ASSOCIATED(es%sector%updo))THEN ! (nup,ndo) SECTOR
      DO ido=1,es%sector%updo%down%dimen
        CALL new_ket(ket_in_do,es%sector%updo%down%state(ido),Ns)
        IF(pm=='+') CALL destroy(ket_out_do,IMPiorbupdo(site,2),ket_in_do)  
        IF(pm=='-') CALL  create(ket_out_do,IMPiorbupdo(site,2),ket_in_do) 
        IF(.NOT.ket_out_do%is_nil)THEN
          DO iup=1,es%sector%updo%up%dimen
 IF(eigen_in%vec%rc(rankupdo(iup,ido,es%sector%updo))/=zero)THEN
            CALL new_ket(ket_in_up,es%sector%updo%up%state(iup),Ns)
            IF(pm=='+') CALL  create(ket_out_up,IMPiorbupdo(site,1),ket_in_up) ! |Ceigen> =  S^+ |eigen>
            IF(pm=='-') CALL destroy(ket_out_up,IMPiorbupdo(site,1),ket_in_up) ! |Ceigen> = -S^- |eigen>
            IF(.NOT.ket_out_up%is_nil)THEN

            !------------------------------------------------! 
            ! WE HAVE TO WRITE S[i]^- = - c[i,up] * C[i,do]  !
            ! TO TAKE CARE BY HAND OF {c[i,up],C[i,do]} = 0  !
            ! SINCE fermion_sign ONLY TAKES CARE OF          !
            ! {c[i,spin],C[j,spin]} = 0                      !
            !------------------------------------------------!

            IF(pm=='-') ket_out_up%fermion_sign = - ket_out_up%fermion_sign

             jup    = Ces%sector%updo%up%  rank(ket_out_up%state) 
             jdo    = Ces%sector%updo%down%rank(ket_out_do%state) 
             istate = rankupdo(iup,ido, es%sector%updo)
             jstate = rankupdo(jup,jdo,Ces%sector%updo)
             eigen_out%vec%rc(jstate) = eigen_out%vec%rc(jstate) + ket_out_up%fermion_sign * ket_out_do%fermion_sign *  eigen_in%vec%rc(istate)
            ENDIF
          ENDIF
 ENDDO
        ENDIF
      ENDDO
      ENDIF

      ! THEN WE UPDATE THE LIST OF OUTPUT VECTORS

      CALL add_eigen(eigen_out,Ces%lowest)
    ENDIF
 ENDDO
    CALL delete_eigen(eigen_out)

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

  SUBROUTINE apply_N(Ces,pm,MASK,es,esrank) 

    !$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ COMPUTE N[site]|0> $$
    !$$$$$$$$$$$$$$$$$$$$$$$$

    TYPE(eigensector_type), INTENT(INOUT) :: Ces
    CHARACTER(LEN=1),       INTENT(IN)    :: pm   ! DUMMY
    LOGICAL,                INTENT(IN)    :: MASK(:)
    TYPE(eigensector_type), INTENT(IN)    :: es
    INTEGER,                INTENT(IN)    :: esrank
    TYPE(eigen_type)                      :: eigen_out
    TYPE(eigen_type), POINTER             :: eigen_in => NULL()
    TYPE(fermion_ket_type)                :: ket_in,ket_in_up,ket_in_do
    INTEGER                               :: istate,iup,ido,nup,ndo,n,site

    eigen_in => es%lowest%eigen(esrank)
    if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
    CALL new_eigen(eigen_out,dimen_func(Ces%sector))

    eigen_out%val       = eigen_in%val
    eigen_out%converged = eigen_in%converged

    ! WE ONLY NEED TO CONSIDER A SUBSET OF ORBITALS
    DO site=1,SIZE(MASK)
 IF(MASK(site))THEN

      ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR

      eigen_out%rank   = site
      eigen_out%vec%rc = zero

      ! THEN PARSE THE INPUT SECTOR TO APPLY RELEVANT SPIN OPERATOR

      IF(ASSOCIATED(es%sector%sz))THEN ! Sz-SECTOR
      DO istate=1,es%sector%sz%dimen
 IF(eigen_in%vec%rc(istate)/=zero)THEN
        CALL new_ket(ket_in,es%sector%sz%state(istate),es%sector%sz%norbs)
        n = ni__(IMPiorbsz(site,1),ket_in) + ( 1 - ni__(IMPiorbsz(site,2),ket_in) ) ! NAMBU 
        IF(n/=zero) eigen_out%vec%rc(istate) = eigen_out%vec%rc(istate) + n * eigen_in%vec%rc(istate) 
      ENDIF
 ENDDO
      ELSE IF(ASSOCIATED(es%sector%updo))THEN ! (nup,ndo) SECTOR
      DO ido=1,es%sector%updo%down%dimen
        CALL new_ket(ket_in_do,es%sector%updo%down%state(ido),Ns)
        ndo = ni__(IMPiorbupdo(site,2),ket_in_do)
        DO iup=1,es%sector%updo%up%dimen
 istate = rankupdo(iup,ido,es%sector%updo)
 IF(eigen_in%vec%rc(istate)/=zero)THEN
          CALL new_ket(ket_in_up,es%sector%updo%up%state(iup),Ns)
          nup = ni__(IMPiorbupdo(site,1),ket_in_up)
          IF(nup+ndo/=0) eigen_out%vec%rc(istate) = eigen_out%vec%rc(istate) + ( nup + ndo ) * eigen_in%vec%rc(istate)
        ENDIF
 ENDDO
      ENDDO
      ENDIF

      ! THEN WE UPDATE THE LIST OF OUTPUT VECTORS

      CALL add_eigen(eigen_out,Ces%lowest)

    ENDIF
 ENDDO

    CALL delete_eigen(eigen_out)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE apply_Id(Ces,pm,MASK,es,esrank)
    TYPE(eigensector_type), INTENT(INOUT) :: Ces
    CHARACTER(LEN=1),       INTENT(IN)    :: pm   ! DUMMY
    LOGICAL,                INTENT(IN)    :: MASK(:)
    TYPE(eigensector_type), INTENT(IN)    :: es
    INTEGER,                INTENT(IN)    :: esrank
    TYPE(eigen_type)                      :: eigen_out
    TYPE(eigen_type), POINTER             :: eigen_in => NULL()
    TYPE(fermion_ket_type)                :: ket_in,ket_in_up,ket_in_do
    INTEGER                               :: istate,iup,ido,nup,ndo,n,site

    eigen_in => es%lowest%eigen(esrank)
    if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
    CALL new_eigen(eigen_out,dimen_func(Ces%sector))

    eigen_out%val       = eigen_in%val
    eigen_out%converged = eigen_in%converged

    ! WE ONLY NEED TO CONSIDER A SUBSET OF ORBITALS
    DO site=1,SIZE(MASK)
 IF(MASK(site))THEN

      ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR

      eigen_out%rank   = site
      eigen_out%vec%rc = zero

      ! THEN PARSE THE INPUT SECTOR TO APPLY RELEVANT SPIN OPERATOR

      IF(ASSOCIATED(es%sector%sz))THEN ! Sz-SECTOR
      DO istate=1,es%sector%sz%dimen
 IF(eigen_in%vec%rc(istate)/=zero)THEN
        eigen_out%vec%rc(istate) = eigen_out%vec%rc(istate) + eigen_in%vec%rc(istate)
      ENDIF
 ENDDO
      ELSE IF(ASSOCIATED(es%sector%updo))THEN ! (nup,ndo) SECTOR
      DO ido=1,es%sector%updo%down%dimen
        DO iup=1,es%sector%updo%up%dimen
 istate = rankupdo(iup,ido,es%sector%updo)
 IF(eigen_in%vec%rc(istate)/=zero)THEN
          eigen_out%vec%rc(istate) = eigen_out%vec%rc(istate) + eigen_in%vec%rc(istate)
        ENDIF
 ENDDO
      ENDDO
      ENDIF

      ! THEN WE UPDATE THE LIST OF OUTPUT VECTORS

      CALL add_eigen(eigen_out,Ces%lowest)

    ENDIF
 ENDDO

    CALL delete_eigen(eigen_out)

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
