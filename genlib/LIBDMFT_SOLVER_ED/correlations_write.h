!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE write_info_correlations(UNIT)
    INTEGER, OPTIONAL, INTENT(IN)  :: UNIT
    INTEGER                        :: unit_

    CALL dump_message(UNIT,TEXT="#############################")
    CALL dump_message(UNIT,TEXT="### IMPURITY CORRELATIONS ###")
    CALL dump_message(UNIT,TEXT="#############################")
                      unit_ = log_unit
    IF(PRESENT(UNIT)) unit_ = UNIT

    WRITE(unit_,'(a,I0)')    "# NUMBER OF LANCZOS ITERATIONS = ",Nitergreenmax
    WRITE(unit_,'(a,E14.7)') "# INVERSE TEMPERATURE          = ",beta 
    WRITE(unit_,'(a,E14.7)') "# BROADENING                   = ",width

    CALL write_array(GNAMBUret%correlstat(2,1)%rc%MASK%imat,"### NAMBU GREEN'S FUNCTION",UNIT=UNIT)
    WRITE(unit_,'(a,f0.6)')  "# REAL FREQ. RANGE = ",MAXVAL(DBLE(GNAMBUret%correl(2,1)%freq%vec))
    CALL write_array( Szret%correlstat(2,1)%rc%MASK%imat,"### < Siz (w) * Sjz >",UNIT=UNIT)
    IF(ANY(Szret%compute))  WRITE(unit_,'(a,f0.6)')  "# REAL FREQ. RANGE = ",MAXVAL(DBLE(    Szret%correl(2,1)%freq%vec))
    CALL write_array(Spmret%correlstat(2,1)%rc%MASK%imat,"### < Si- (w) * Sj+ >",UNIT=UNIT)
    IF(ANY(Spmret%compute)) WRITE(unit_,'(a,f0.6)')  "# REAL FREQ. RANGE = ",MAXVAL(DBLE(   Spmret%correl(2,1)%freq%vec))
    CALL write_array(  Nret%correlstat(2,1)%rc%MASK%imat,"### <  Ni (w) * Nj  >",UNIT=UNIT)
    IF(ANY(Nret%compute))   WRITE(unit_,'(a,f0.6)')  "# REAL FREQ. RANGE = ",MAXVAL(DBLE(     Nret%correl(2,1)%freq%vec))
    CALL write_array( P3ret%correlstat(2,1)%rc%MASK%imat,"### < P3i (w) * P3j >",UNIT=UNIT)
    IF(ANY(P3ret%compute))  WRITE(unit_,'(a,f0.6)')  "# REAL FREQ. RANGE = ",MAXVAL(DBLE(    P3ret%correl(2,1)%freq%vec))
    CALL write_array( P4ret%correlstat(2,1)%rc%MASK%imat,"### < P4i (w) * P4j >",UNIT=UNIT) 
    IF(ANY(P4ret%compute))  WRITE(unit_,'(a,f0.6)')  "# REAL FREQ. RANGE = ",MAXVAL(DBLE(    P4ret%correl(2,1)%freq%vec)) 
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

  SUBROUTINE write_static(AIM,GS,UNIT)
    TYPE(AIM_type),             INTENT(IN) :: AIM
    TYPE(eigensectorlist_type), INTENT(IN) :: GS
    INTEGER, OPTIONAL, INTENT(IN)          :: UNIT
#ifdef _complex
    COMPLEX(DBL)                           :: meanStot2,meanStot,meanDStot 
#else
    REAL(DBL)                              :: meanStot2,meanStot,meanDStot 
#endif
    REAL(DBL), ALLOCATABLE                 :: CHImean(:) 
    REAL(DBL), ALLOCATABLE                 :: dblocc(:),magnet(:) 
    REAL(DBL), ALLOCATABLE                 :: densdens(:,:)
    INTEGER,   ALLOCATABLE                 :: tabpairG(:,:),tabpairF(:,:),tabpairNS(:,:),tabpairP3(:,:),tabpairP4(:,:)
    TYPE(masked_matrix_type)               :: SS,CHIstat,NN,Gdia
    CHARACTER(LEN=200)                     :: fmt_loc,fmt_G,fmt_SS,fmt_NN,fmt_P3,fmt_P4,fmt_F
    INTEGER                                :: unit_,iGS,site1,site2,spin,iP,site,Nc,iind,ipm

                      unit_ = log_unit
    IF(PRESENT(UNIT)) unit_ = UNIT

    CALL dump_message(UNIT=UNIT,TEXT="#")
    CALL dump_message(UNIT=UNIT,TEXT="###############################") 
    CALL dump_message(UNIT=UNIT,TEXT="### EQUAL-TIME CORRELATIONS ###")
    CALL dump_message(UNIT=UNIT,TEXT="###############################")   
    CALL dump_message(UNIT=UNIT,TEXT="#")

    Nc = AIM%Nc
    WRITE(fmt_loc,*) '(a,',Nc,'(f9.6,x))'

    CALL dump_message(UNIT=UNIT,TEXT="###################") 
    CALL dump_message(UNIT=UNIT,TEXT="### LOCAL MEANS ###")
    CALL dump_message(UNIT=UNIT,TEXT="###################")   

    if(allocated(magnet)) deallocate(magnet,dblocc); ALLOCATE(magnet(Nc),dblocc(Nc))

    DO site1=1,Nc
      dens(site1)    =        DBLE(G(1)%correlstat(1,2)%rc%mat(site1,site1) + G(2)%correlstat(1,2)%rc%mat(site1,site1))
      magnet(site1)  = half * DBLE(G(1)%correlstat(1,2)%rc%mat(site1,site1) - G(2)%correlstat(1,2)%rc%mat(site1,site1))
      dblocc(site1)  = half * DBLE(N%correlstat(2,1)%rc%mat(site1,site1) - dens(site1))
    ENDDO

    WRITE(unit_,fmt_loc) "# density       = ",dens 
    WRITE(unit_,fmt_loc) "# double   occ. = ",dblocc
    WRITE(unit_,fmt_loc) "# magnetization = ",magnet*2     
    tot_repulsion=sum(dblocc)
    WRITE(unit_,fmt_loc) "# <U>           = ",tot_repulsion

    call write_array(G(1)%correlstat(1,1)%rc%mat,'G1 11',unit=unit_)
    call write_array(G(1)%correlstat(1,2)%rc%mat,'G1 12',unit=unit_)
    call write_array(G(1)%correlstat(2,1)%rc%mat,'G1 21',unit=unit_)
    call write_array(G(1)%correlstat(2,2)%rc%mat,'G1 22',unit=unit_)
    call write_array(G(2)%correlstat(1,1)%rc%mat,'G2 11',unit=unit_)
    call write_array(G(2)%correlstat(1,2)%rc%mat,'G2 12',unit=unit_)
    call write_array(G(2)%correlstat(2,1)%rc%mat,'G2 21',unit=unit_)
    call write_array(G(2)%correlstat(2,2)%rc%mat,'G2 22',unit=unit_)

    ! BOND OPERATORS !

    if(allocated(tabpairG)) deallocate(tabpairG) ; ALLOCATE(tabpairG(G(1)%correlstat(2,1)%rc%MASK%nind,2))
    CALL order_matrix_elements(tabpairG,G(1)%correlstat(2,1)%rc%MASK)
    CALL dump_message(UNIT=UNIT,TEXT="######################") 
    CALL dump_message(UNIT=UNIT,TEXT="### BOND OPERATORS ###")
    CALL dump_message(UNIT=UNIT,TEXT="######################")   
    CALL new_masked_matrix(Gdia,G(1)%correlstat(1,2))
    Gdia%rc%mat = G(1)%correlstat(1,2)%rc%mat + G(2)%correlstat(1,2)%rc%mat

#ifdef _complex
    CALL dump_message(UNIT=UNIT,TEXT="# link"//&
    "     <C[i,"//ccspin(1)//"]*c[j,"//ccspin(1)//"]>"//&
    "       <C[j,"//ccspin(2)//"]*c[i,"//ccspin(2)//"]>"//&
    "            <C[j]*c[i]>")
#else
    CALL dump_message(UNIT=UNIT,TEXT="# link"//&
    "     <C[i,"//ccspin(1)//"]*c[j,"//ccspin(1)//"]>"//&
    "       <C[j,"//ccspin(2)//"]*c[i,"//ccspin(2)//"]>"//&
    "          <C[j]*c[i]>")
#endif


    DO iind=1,SIZE(tabpairG,1)

      site1 = tabpairG(iind,1) ; site2 = tabpairG(iind,2)

#ifdef _complex
      WRITE(fmt_G,*) '(2x,I0,x,I0,4x,2(2(a,f9.6),a,3x),2x,2(a,f9.6),a)'
      WRITE(unit_,fmt_G)  site1,site2,&
      "(",DBLE(G(1)%correlstat(1,2)%rc%mat(site1,site2)),",",AIMAG(G(1)%correlstat(1,2)%rc%mat(site1,site2)),")",&
      "(",DBLE(G(2)%correlstat(1,2)%rc%mat(site1,site2)),",",AIMAG(G(2)%correlstat(1,2)%rc%mat(site1,site2)),")",&
      "(",DBLE(Gdia%rc%mat(site1,site2)),",",AIMAG(Gdia%rc%mat(site1,site2)),")"
#else
      WRITE(fmt_G,*) '(2x,I0,x,I0,10x,3(f9.6,15x))'
      WRITE(unit_,fmt_G)  site1,site2,&
      G(1)%correlstat(1,2)%rc%mat(site1,site2),G(2)%correlstat(1,2)%rc%mat(site1,site2),Gdia%rc%mat(site1,site2)
#endif
    ENDDO

    ! PAIRING OPERATORS

    IF(SUPER)THEN

      if(allocated(tabpairF)) deallocate(tabpairF) ; ALLOCATE(tabpairF(GF(1)%correlstat(1,1)%rc%MASK%nind,2))
      CALL order_matrix_elements(tabpairF,GF(1)%correlstat(1,1)%rc%MASK)
      CALL dump_message(UNIT=UNIT,TEXT="# link"//&
      "     <C[i,"//ccspin(1)//"]*C[j,"//ccspin(2)//"]>"//&
      "       <c[i,"//ccspin(1)//"]*c[j,"//ccspin(2)//"]>"//&
      "       <C[i,"//ccspin(2)//"]*C[j,"//ccspin(1)//"]>"//&
      "       <c[i,"//ccspin(2)//"]*c[j,"//ccspin(1)//"]>")
      DO iind=1,SIZE(tabpairF,1)
       site1 = tabpairF(iind,1)
       site2 = tabpairF(iind,2)

#ifdef _complex
       WRITE(fmt_F,*) '(2x,I0,x,I0,4x,4(2(a,f9.6),a,4x))'
       WRITE(unit_,fmt_F)  site1,site2,&
       "(",DBLE(GF(1)%correlstat(1,1)%rc%mat(site1,site2)),",",AIMAG(GF(1)%correlstat(1,1)%rc%mat(site1,site2)),")",&
       "(",DBLE(GF(1)%correlstat(2,2)%rc%mat(site1,site2)),",",AIMAG(GF(1)%correlstat(2,2)%rc%mat(site1,site2)),")",&
       "(",DBLE(GF(2)%correlstat(1,1)%rc%mat(site1,site2)),",",AIMAG(GF(2)%correlstat(1,1)%rc%mat(site1,site2)),")",&
       "(",DBLE(GF(2)%correlstat(2,2)%rc%mat(site1,site2)),",",AIMAG(GF(2)%correlstat(2,2)%rc%mat(site1,site2)),")"
#else
       WRITE(fmt_F,*) '(2x,I0,x,I0,10x,4(f9.6,15x))'
       WRITE(unit_,fmt_F)  site1,site2,&
       GF(1)%correlstat(1,1)%rc%mat(site1,site2),GF(1)%correlstat(2,2)%rc%mat(site1,site2),&
       GF(1)%correlstat(2,2)%rc%mat(site2,site1),GF(1)%correlstat(1,1)%rc%mat(site2,site1) ! USE SYMMETRY OF THE REAL CASE
#endif
      ENDDO
    ENDIF

    ! SPIN CORRELATIONS !

    if( donot_compute_holepart_spm )then
      Spm%correlstat(1,2)%rc%mat = conj(Spm%correlstat(2,1)%rc%mat)
    endif

    CALL dump_message(UNIT=UNIT,TEXT="#########################") 
    CALL dump_message(UNIT=UNIT,TEXT="### SPIN CORRELATIONS ###")
    CALL dump_message(UNIT=UNIT,TEXT="#########################")   
#ifdef _complex
    CALL dump_message(UNIT=UNIT,TEXT="# link     <S-[i]*S+[j]>/2        <S+[i]*S-[j]>/2        <Sz[i]*Sz[j]>             <S[i].S[j]>")
#else
    CALL dump_message(UNIT=UNIT,TEXT="# link    <S-[i]*S+[j]>/2    <S+[i]*S-[j]>/2    <Sz[i]*Sz[j]>        <S[i].S[j]>")
#endif
    CALL new_masked_matrix(SS,Sz%correlstat(2,1)) 

    SS%rc%mat = half * (Spm%correlstat(2,1)%rc%mat + Spm%correlstat(1,2)%rc%mat) + Sz%correlstat(2,1)%rc%mat
    if(allocated(tabpairNS)) deallocate(tabpairNS) ; ALLOCATE(tabpairNS(Sz%correlstat(2,1)%rc%MASK%nind,2))
    CALL order_matrix_elements(tabpairNS,Sz%correlstat(2,1)%rc%MASK)

    DO iind=1,Sz%correlstat(2,1)%rc%MASK%nind
      site1 = tabpairNS(iind,1)
      site2 = tabpairNS(iind,2)
#ifdef _complex
      WRITE(fmt_SS,*) '(2x,I0,x,I0,2x,3(2(a,f9.6),a,2x),2x,2(a,f9.6),a)'
      WRITE(unit_,fmt_SS)  site1,site2,&
      "(",half*DBLE(Spm%correlstat(2,1)%rc%mat(site1,site2)),",",half*AIMAG(Spm%correlstat(1,2)%rc%mat(site1,site2)),")",&
      "(",half*DBLE(Spm%correlstat(2,1)%rc%mat(site1,site2)),",",half*AIMAG(Spm%correlstat(1,2)%rc%mat(site1,site2)),")",&
      "(",     DBLE( Sz%correlstat(2,1)%rc%mat(site1,site2)),",",     AIMAG( Sz%correlstat(1,2)%rc%mat(site1,site2)),")",&
      "(",     DBLE(               SS%rc%mat(site1,site2)),",",       AIMAG(                 SS%rc%mat(site1,site2)),")"
#else
      WRITE(fmt_SS,*) '(2x,I0,x,I0,7x,3(f9.6,10x),f9.6)'
      WRITE(unit_,fmt_SS)  site1,site2,&
      half*Spm%correlstat(2,1)%rc%mat(site1,site2),half*Spm%correlstat(1,2)%rc%mat(site1,site2),Sz%correlstat(2,1)%rc%mat(site1,site2),SS%rc%mat(site1,site2)
#endif
    ENDDO
    meanStot2 = SUM(Sz%correlstat(2,1)%rc%mat + half * ( Spm%correlstat(2,1)%rc%mat + Spm%correlstat(1,2)%rc%mat ) )
    meanStot  = SUM(magnet) ! this is true provided <S+> = <S-> = 0 
    meanDStot = SQRT(abs(meanStot2 - meanStot**2))

#ifdef _complex
    WRITE(unit_,'(2(a,f9.6),a)') "# <Stot*Stot> = (",DBLE(meanStot2),",",AIMAG(meanStot2),")"
    WRITE(unit_,'(2(a,f9.6),a)') "# <Stot>      = (",DBLE(meanStot),",",AIMAG(meanStot),")"
    WRITE(unit_,'(2(a,f9.6),a)') "# SQRT(<Stot*Stot> - <Stot>^2) = (",DBLE(meanDStot),",",AIMAG(meanDStot),")"
#else
    WRITE(unit_,'(a,f9.6)') "# <Stot*Stot> = ",meanStot2
    WRITE(unit_,'(a,f9.6)') "# <Stot>      = ",meanStot
    WRITE(unit_,'(a,f9.6)') "# SQRT(<Stot*Stot> - <Stot>^2) = ",meanDStot
#endif

    ! DENSITY CORRELATIONS

    CALL dump_message(UNIT=UNIT,TEXT="############################") 
    CALL dump_message(UNIT=UNIT,TEXT="### DENSITY CORRELATIONS ###")
    CALL dump_message(UNIT=UNIT,TEXT="############################")   
    if(allocated(densdens)) deallocate(densdens) ; ALLOCATE(densdens(Nc,Nc))
    DO site1=1,Nc; DO site2=1,Nc
      densdens(site1,site2) = dens(site1) * dens(site2)
    ENDDO; ENDDO
    CALL new_masked_matrix(NN,N%correlstat(2,1))
    NN%rc%mat = N%correlstat(2,1)%rc%mat - densdens
#ifdef _complex
    CALL dump_message(UNIT=UNIT,TEXT="# link      <n[i]*n[j]>          <n[i]>*<n[j]>      <n[i]*n[j]> - <n[i]>*<n[j]>")
#else
    CALL dump_message(UNIT=UNIT,TEXT="# link     <n[i]*n[j]>       <n[i]>*<n[j]>      <n[i]*n[j]> - <n[i]>*<n[j]>")
#endif
    DO iind=1,NN%rc%MASK%nind
      site1 = tabpairNS(iind,1)
      site2 = tabpairNS(iind,2)
#ifdef _complex
      WRITE(fmt_NN,*) '(2x,I0,x,I0,2x,2(a,f9.6),a,7x,f9.6,10x,2(a,f9.6),a)'
      WRITE(unit_,fmt_NN)  site1,site2,&
      "(",DBLE( N%correlstat(2,1)%rc%mat(site1,site2)),",",AIMAG( N%correlstat(2,1)%rc%mat(site1,site2)),")",&
      densdens(site1,site2),"(",DBLE(NN%rc%mat(site1,site2)),",",AIMAG(NN%rc%mat(site1,site2)),")"
#else
      WRITE(fmt_NN,*) '(2x,I0,x,I0,6x,2(f9.6,11x),4x,f9.6)'
      WRITE(unit_,fmt_NN)  site1,site2,&
      N%correlstat(2,1)%rc%mat(site1,site2),densdens(site1,site2),NN%rc%mat(site1,site2)
#endif
    ENDDO


    !------------------------------------------------------------------------------------------------!
    ! Pijk CORRELATIONS
 
    if(compute_P3)then

    CALL dump_message(UNIT=UNIT,TEXT="#########################") 
    CALL dump_message(UNIT=UNIT,TEXT="### Pijk CORRELATIONS ###")
    CALL dump_message(UNIT=UNIT,TEXT="#########################")   
    if(allocated(CHImean)) deallocate(CHImean); ALLOCATE(CHImean(SIZE(triplets,1)))
    CHImean = zero

#ifdef _complex
    ! CHI
    ! WARNING: sign is different from Wen, Wilczek & Zee cause our Pijk(1) is i->j->k->i
    CHImean = half * AIMAG(P3%Amean(:,1)) 
    CALL dump_message(UNIT=UNIT,TEXT="# (i,j,k)          <Pijk>                   <Pikj>            <Si.(SjxSk)>")    
#else
    CALL dump_message(UNIT=UNIT,TEXT="# (i,j,k)      <Pijk>           <Pikj>        <Si.(SjxSk)>")    
#endif


    DO iP=1,SIZE(triplets,1)

#ifdef _complex
      WRITE(unit_,'(3(a,I0),a,2(2(a,f9.6),a),a,f9.6)') &
      "  (",(triplets(iP,site),",",site=1,2),triplets(iP,3),")  ",&
      ("(",DBLE(P3%Amean(iP,ipm)),",",AIMAG(P3%Amean(iP,ipm)),")    ",ipm=1,2),"  ",CHImean(iP)
#else
      WRITE(unit_,'(3(a,I0),a,2(f9.6,a),f9.6)') &
      "  (",(triplets(iP,site),",",site=1,2),triplets(iP,3),")    ",&
      (P3%Amean(iP,ipm),"        ",ipm=1,2),CHImean(iP)
#endif

    ENDDO

    CALL new_masked_matrix(CHIstat,P3%correlstat(2,1))
    CHIstat%rc%mat=-( P3%correlstat(1,1)%rc%mat + P3%correlstat(2,2)%rc%mat - P3%correlstat(1,2)%rc%mat - P3%correlstat(2,1)%rc%mat ) * quarter**2
    if(allocated(tabpairP3)) deallocate(tabpairP3) ; ALLOCATE(tabpairP3(P3%correlstat(2,1)%rc%MASK%nind,2))
    CALL order_matrix_elements(tabpairP3,P3%correlstat(2,1)%rc%MASK)
#ifdef _complex
    CALL dump_message(UNIT=UNIT,TEXT="# (i,j,k)  (l,m,n)        <Pijk*Plmn>             <Pikj*Plnm>             <Pikj*Plmn>             <Pijk*Plnm>          <Si.(SkxSk)*Sl.(SmxSn)>")    
#else
    CALL dump_message(UNIT=UNIT,TEXT="# (i,j,k)  (l,m,n)     <Pijk*Plmn>          <Pikj*Plnm>          <Pikj*Plmn>          <Pijk*Plnm>          <Si.(SkxSk)*Sl.(SmxSn)>")    
#endif
    DO iind=1,CHIstat%rc%MASK%nind
      site1 = tabpairP3(iind,1)
      site2 = tabpairP3(iind,2)
#ifdef _complex
      WRITE(fmt_P3,*) '(3(a,I0),3(a,I0),a,2x,4(2(a,f9.6),a,3x),3x,2(a,f9.6),a)'
      WRITE(unit_,fmt_P3) &
      "  (",(triplets(site1,site),",",site=1,2),triplets(site1,3),")  (",(triplets(site2,site),",",site=1,2),triplets(site2,3),") ",&
      "(",DBLE(P3%correlstat(1,1)%rc%mat(site1,site2)),",",AIMAG(P3%correlstat(1,1)%rc%mat(site1,site2)),")",&
      "(",DBLE(P3%correlstat(2,2)%rc%mat(site1,site2)),",",AIMAG(P3%correlstat(2,2)%rc%mat(site1,site2)),")",&
      "(",DBLE(P3%correlstat(2,1)%rc%mat(site1,site2)),",",AIMAG(P3%correlstat(2,1)%rc%mat(site1,site2)),")",&
      "(",DBLE(P3%correlstat(1,2)%rc%mat(site1,site2)),",",AIMAG(P3%correlstat(1,2)%rc%mat(site1,site2)),")",&
      "(",DBLE(CHIstat%rc%mat(site1,site2)),",",AIMAG(CHIstat%rc%mat(site1,site2)),")"
#else
      WRITE(fmt_P3,*) '(3(a,I0),3(a,I0),a,2x,4(f9.6,12x),6x,f9.6)'
      WRITE(unit_,fmt_P3) &
      "  (",(triplets(site1,site),",",site=1,2),triplets(site1,3),")  (",(triplets(site2,site),",",site=1,2),triplets(site2,3),")    ",&
      P3%correlstat(1,1)%rc%mat(site1,site2),P3%correlstat(2,2)%rc%mat(site1,site2),&
      P3%correlstat(2,1)%rc%mat(site1,site2),P3%correlstat(1,2)%rc%mat(site1,site2),CHIstat%rc%mat(site1,site2)
#endif
    ENDDO

    endif
    !------------------------------------------------------------------------------------------------!



    !------------------------------------------------------------------------------------------------!
    if(compute_P4)then

    ! Pijkl CORRELATIONS

    CALL dump_message(UNIT=UNIT,TEXT="##########################") 
    CALL dump_message(UNIT=UNIT,TEXT="### Pijkl CORRELATIONS ###")
    CALL dump_message(UNIT=UNIT,TEXT="##########################")   
#ifdef _complex
    CALL dump_message(UNIT=UNIT,TEXT="# (i,j,k,l)           <Pijkl>                  <Pilkj>")    
#else
    CALL dump_message(UNIT=UNIT,TEXT="# (i,j,k,l)      <Pijkl>          <Pilkj>")    
#endif
    DO iP=1,SIZE(quadruplets,1)
#ifdef _complex
      WRITE(unit_,'(4(a,I0),a,2(2(a,f9.6),a))') &
      "  (",(quadruplets(iP,site),",",site=1,3),quadruplets(iP,4),")    ",&
      ("(",DBLE(P4%Amean(iP,ipm)), ",",AIMAG(P4%Amean(iP,ipm)),")    ",ipm=1,2)
#else
      WRITE(unit_,'(4(a,I0),a,f9.6,a,f9.6,a)') &
      "  (",(quadruplets(iP,site),",",site=1,3),quadruplets(iP,4),")    ",&
      (P4%Amean(iP,ipm),"        ",ipm=1,2)
#endif
    ENDDO
    ! REORDER PAIRS OF QUADRUPLETS: DIAGONAL FIRST + GATHER (i,j) (j,i) 
    if(allocated(tabpairP4)) deallocate(tabpairP4) ; ALLOCATE(tabpairP4(P4%correlstat(2,1)%rc%MASK%nind,2))
    CALL order_matrix_elements(tabpairP4,P4%correlstat(2,1)%rc%MASK)
#ifdef _complex
    CALL dump_message(UNIT=UNIT,TEXT="# (i,j,k,l)  (m,n,o,p)        <Pijkl*Pmnop>           <Pilkj*Pmpon>           <Pilkj*Pmnop>           <Pijkl*Pmpon>")    
#else
    CALL dump_message(UNIT=UNIT,TEXT="# (i,j,k,l)  (m,n,o,p)     <Pijkl*Pmnop>          <Pilkj*Pmpon>          <Pilkj*Pmnop>          <Pijkl*Pmpon>")    
#endif
    DO iind=1,P4%correlstat(2,1)%rc%MASK%nind
      site1 = tabpairP4(iind,1)
      site2 = tabpairP4(iind,2)
#ifdef _complex
      WRITE(fmt_P4,*) '(4(a,I0),4(a,I0),a,4(2(a,f9.6),a,3x))'
      WRITE(unit_,fmt_P4) "  (",(quadruplets(site1,site),",",site=1,3),quadruplets(site1,4),&
      ")  (",(quadruplets(site2,site),",",site=1,3),quadruplets(site2,4),")    ",&
      "(",DBLE(P4%correlstat(1,1)%rc%mat(site1,site2)),",",AIMAG(P4%correlstat(1,1)%rc%mat(site1,site2)),")",&
      "(",DBLE(P4%correlstat(2,2)%rc%mat(site1,site2)),",",AIMAG(P4%correlstat(2,2)%rc%mat(site1,site2)),")",&
      "(",DBLE(P4%correlstat(2,1)%rc%mat(site1,site2)),",",AIMAG(P4%correlstat(2,1)%rc%mat(site1,site2)),")",&
      "(",DBLE(P4%correlstat(1,2)%rc%mat(site1,site2)),",",AIMAG(P4%correlstat(1,2)%rc%mat(site1,site2)),")"
#else
      WRITE(fmt_P4,*) '(4(a,I0),4(a,I0),a,3x,4(f9.6,14x))'
      WRITE(unit_,fmt_P4) "  (",(quadruplets(site1,site),",",site=1,3),quadruplets(site1,4),&
      ")  (",(quadruplets(site2,site),",",site=1,3),quadruplets(site2,4),")    ",&
      P4%correlstat(1,1)%rc%mat(site1,site2),P4%correlstat(2,2)%rc%mat(site1,site2),&
      P4%correlstat(2,1)%rc%mat(site1,site2),P4%correlstat(1,2)%rc%mat(site1,site2)
#endif
    ENDDO

    endif
    !------------------------------------------------------------------------------------------------!


    if(allocated(dblocc))    deallocate(dblocc,magnet,densdens)

    if(allocated(CHImean))   deallocate(CHImean)  ; if(allocated(tabpairG))  deallocate(tabpairG)
    if(allocated(tabpairNS)) deallocate(tabpairNS); if(allocated(tabpairP3)) deallocate(tabpairP3)
    if(allocated(tabpairP4)) deallocate(tabpairP4)

    IF(SUPER) then
     if(allocated(tabpairF)) DEALLOCATE(tabpairF)
    ENDIF

    CALL delete_masked_matrix(Gdia); CALL delete_masked_matrix(SS); CALL delete_masked_matrix(NN); CALL delete_masked_matrix(CHIstat)

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
