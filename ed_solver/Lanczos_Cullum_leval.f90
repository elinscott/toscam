MODULE Lanczos_Cullum_leval

  USE Lanczos_Cullum_lesub

  IMPLICIT NONE

  REAL(DBL),PARAMETER, PRIVATE  :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,  PARAMETER, PRIVATE  :: F=.FALSE.,T=.TRUE.

  ! TO SAVE EIGENVALUES
  INTEGER,              PRIVATE :: NG,UNIT1_,UNIT2_,UNIT3_,UNIT4_,UNIT5_
  REAL(8), ALLOCATABLE, PRIVATE :: V2(:)


CONTAINS


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE MATVEC(V1,V2,COEF)
    ! WRAPPER FOR SYMMETRIC LANCZOS ROUTINES BY CULLUM
    REAL(8) :: V1(:)
    REAL(8) :: V2(:)
    REAL(8) :: V2_IN(SIZE(V2))
    REAL(8) :: COEF
    ! COMPUTE V2 = H * V1 - COEF * V2
    if(size(V1)/=size(V2))then
     write(*,*) 'ERROR in MATVEC (Cullum modules), size do not match'
     write(*,*) 'sizes of V1 and V2 : ',size(V1),size(V2)
     stop
    endif
    V2_IN = V2
    CALL Hmult__(size(V1),V2,V1)
    V2 = V2 - COEF * V2_IN
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

  SUBROUTINE find_eigenval_r(FILE1,FILE3,FILE4,FILE5,lowest)
    INTEGER                          ::  ZESEED(ran_siz())
    CHARACTER(LEN=*),     INTENT(IN) ::  FILE1,FILE3,FILE4,FILE5
    TYPE(eigenlist_type)             ::  lowest
    TYPE(eigen_type)                 ::  eigen
    INTEGER                          ::  N_,Neigen_in_window,ieigen,seedsize

    call build_seed_and_shot_in_genran(ZESEED,iseed)

    ! THE SIZE OF THE HAMILTONIAN IS THE ONLY VARYING PARAMETER!

    N_ = dimen_H()
    ALLOCATE(V2(MAX(N_,KMAX_)))

    ! THIS GENERATES ON THE FLY 
    ! THE CORRECT INPUT FILE 'FILE5'

   CALL open_safe(UNIT5_,FILE5,'UNKNOWN',get_unit=.true.)
    if(rank==0) WRITE(UNIT5_,'(a)') "LEVAL INPUT EIGENVALUE COMPUTATION, NO REORTHOGONALIZATION"
    if(rank==0) WRITE(UNIT5_,'(a)') "SYMMETRIC HAMILTONIAN" 
    if(rank==0) WRITE(UNIT5_,'(a)') "           N    KMAX    NMEVS     MATNO"
    if(rank==0) WRITE(UNIT5_,'(4I12)') N_,KMAX_,NMEVS_,MATNO__
    if(rank==0) WRITE(UNIT5_,'(a)') "        SEEDSIZE"
    if(rank==0) WRITE(UNIT5_,'(I12)') SIZE(ZESEED)
    if(rank==0) WRITE(UNIT5_,'(a)') "         SVSEED"
    if(rank==0) WRITE(UNIT5_,*) ZESEED
    if(rank==0) WRITE(UNIT5_,'(a)') "         RHSEED"
    if(rank==0) WRITE(UNIT5_,*) ZESEED+1
    if(rank==0) WRITE(UNIT5_,'(a)') "         MXINIT    MXSTUR" 
    if(rank==0) WRITE(UNIT5_,'(2I12)') MXINIT_,MXSTUR_
    if(rank==0) WRITE(UNIT5_,'(a)') "         ISTART     ISTOP"
    if(rank==0) WRITE(UNIT5_,'(2I12)') ISTART_,ISTOP_
    if(rank==0) WRITE(UNIT5_,'(a)') "           IHIS     IDIST  IWRITE"
    if(rank==0) WRITE(UNIT5_,'(3I12)') IHIS_,IDIST_,IWRITE__
    if(rank==0) WRITE(UNIT5_,'(a)') "         RELTOL (RELATIVE TOLERANCE IN 'COMBINING' GOODEV)" 
    if(rank==0) WRITE(UNIT5_,'(E20.12)') RELTOL__
    if(rank==0) WRITE(UNIT5_,'(a)')"          MB(1)   MB(2)   MB(3)   MB(4)    ...  (ORDERS OF T(1,MEV) )"
    if(rank==0) WRITE(UNIT5_,*) NMEV_
    if(rank==0) WRITE(UNIT5_,'(a)') "          NINT     (NUMBER OF SUB-INTERVALS FOR BISEC)"
    if(rank==0) WRITE(UNIT5_,*) NINT_
    if(rank==0) WRITE(UNIT5_,'(a)') "          LB(1)    LB(2)   LB(3)   LB(4)   ... (INTERVAL LOWER BOUNDS)"
    if(rank==0) WRITE(UNIT5_,*) LB_
    if(rank==0) WRITE(UNIT5_,'(a)') "          UB(1)    UB(2)   UB(3)   UB(4)   ... (INTERVAL UPPER BOUNDS)"
    if(rank==0) WRITE(UNIT5_,*) UB_
   CALL close_safe(UNIT5_)

    ! AND WE RUN LEVAL

    call mpibarrier 
    CALL leval(FILE1=FILE1,FILE3=FILE3,FILE4=FILE4,FILE5=FILE5,SILENT=F)
    call mpibarrier

    ! AND WE CORRECT OUPUT FILES + SAVE OUTPUT EIGENVALUES

    IF(NG==0) STOP "ERROR IN LEVAL: NO GOOD EIGENVALUE FOUND !"
    DO ieigen=1,MIN(Neigen,NG)
      CALL new_eigen(eigen,N_)
      eigen%val       = V2(ieigen)
      eigen%converged = F
      IF(V2(ieigen)<MINVAL(V2(1:NG))+dEmax) eigen%converged = T ! FLAG OUT THE EIGENVALUES OUTSIDE THE WINDOW
      eigen%rank      = ieigen
      CALL add_eigen(eigen,lowest) ! APPEND NEW EIGENPAIR TO LIST
      CALL delete_eigen(eigen)
    ENDDO

    ! CORRECT OUTPUT FILES WITH DESIRED EIGENVALUES IN WINDOW
 
    Neigen_in_window = COUNT(lowest%eigen(:)%converged)
 
    IF(iproc==1)THEN
      ! NG IS THE 1ST NUMBER ON THE 1ST LINE OF FILE3
      CALL system_call("/usr/bin/perl -p -i -e '$n=0;s{(\\d+)}{++$n==1?"//c2s(i2c(Neigen_in_window))//":$1}ige if $.==1' "//TRIM(FILE3))
      ! NG IS THE 3RD NUMBER ON THE 4TH LINE OF FILE4
      CALL system_call("/usr/bin/perl -p -i -e '$n=0;s{(\\d+)}{++$n==3?"//c2s(i2c(Neigen_in_window))//":$1}ige if $.==4' "//TRIM(FILE4))
    ENDIF
 
    DEALLOCATE(V2)
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE leval(FILE1,FILE3,FILE4,FILE5,SILENT)

    !-----LEVAL  (EIGENVALUES OF REAL SYMMETRIC MATRICES)-------------------
    !     CONTAINS MAIN PROGRAM FOR COMPUTING DISTINCT EIGENVALUES OF
    !     A REAL SYMMETRIC MATRIX USING LANCZOS TRIDIAGONALIZATION
    !     WITHOUT REORTHOGONALIZATION.
    !
    !     PFORT VERIFIER IDENTIFIED THE FOLLOWING NONPORTABLE
    !     CONSTRUCTIONS
    !
    !     1.  DATA/MACHEP/ STATEMENT
    !     2.  ALL READ(UNIT5_,*) STATEMENTS (FREE FORMAT)
    !     3.  FORMAT(20A4) USED WITH EXPLANATORY HEADER EXPLAN.
    !     4.  HEXADECIMAL FORMAT (4Z20) USED IN ALPHA/BETA FILES 1 AND 2.
    !-----------------------------------------------------------------------

    LOGICAL,          OPTIONAL, INTENT(IN) :: SILENT
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FILE1,FILE3,FILE4,FILE5
    LOGICAL                                :: silent_
    REAL(8), ALLOCATABLE                   :: V1(:),VS(:)
    REAL(8), ALLOCATABLE                   :: ALPHA(:),BETA(:),LB(:),UB(:)
    REAL(8),          ALLOCATABLE          :: G(:)
    INTEGER,          ALLOCATABLE          :: MP(:),NMEV(:)
    INTEGER,          ALLOCATABLE          :: SVSEED(:),RHSEED(:),SVSOLD(:),IIX(:)
    REAL(8)                                :: BTOL,GAPTOL,TTOL,MACHEP,EPSM,RELTOL
    REAL(8)                                :: SCALE1,SCALE2,SCALE3,SCALE4,BISTOL,CONTOL,MULTOL
    REAL(8)                                :: ONE,ZERO,TEMP,TKMAX,BETAM,BKMIN,T0,T1
    REAL                                   :: EXPLAN(20)
    INTEGER                                :: seedsize,IBMEV,NGM1,IGOOD,IBB,NGOOD,NDIS,IBIST 
    INTEGER                                :: ITEN,NISOM,IWRITO,NM1,MP1,IC,LOOP,IT,IFF,II,NISO
    integer                                :: ICONV,MOLD,MOLD1,ICT,MMB,IPROJ,N,KMAX,NMEVS,MATNO,MEV
    integer                                :: NOLD, MATOLD,ISTART,ISTOP,NINT,IB
    integer                                :: KMAX1,ITEMP
    integer                                :: MXINIT,MXSTUR,J,K,I,IHIS,IDIST,IWRITE

    silent_ = T
    IF(PRESENT(SILENT)) silent_ = SILENT

    IF(PRESENT(FILE1)) then; call open_safe(UNIT1_,FILE1,'UNKNOWN',get_unit=.true.); endif
    IF(PRESENT(FILE3)) then; call open_safe(UNIT3_,FILE3,'UNKNOWN',get_unit=.true.); endif
    IF(PRESENT(FILE4)) then; call open_safe(UNIT4_,FILE4,'UNKNOWN',get_unit=.true.); endif
    IF(PRESENT(FILE5)) then; call open_safe(UNIT5_,FILE5,'UNKNOWN',get_unit=.true.) ;endif

    EPSM = 2.0D0*epsilon(BTOL)

    !-----------------------------------------------------------------------
    !
    !     ARRAYS MUST BE DIMENSIONED AS FOLLOWS:
    !     DIMENSION OF V2 ASSUMES THAT NO MORE THAN KMAX/2 EIGENVALUES
    !     OF THE T-MATRICES ARE BEING COMPUTED IN ANY ONE OF THE
    !     SUB-INTERVALS BEING CONSIDERED.  V2 CONTAINS THE UPPER AND LOWER
    !     BOUNDS FOR EACH T-EIGENVALUE BEING COMPUTED BY BISEC IN ANY ONE
    !     GIVEN INTERVAL.
    !
    !     1.  ALPHA: >= KMAX,   BETA: >= (KMAX+1)
    !     2.  V1:    >= MAX(N,KMAX+1)
    !     3.  V2:    >= MAX(N,KMAX)
    !     4.  VS:    >= KMAX
    !     5.  G:     >= MAX(N,2*KMAX)
    !     6.  MP:    >= KMAX
    !     7.  LB,UB: >= NUMBER OF SUBINTERVALS SUPPLIED TO BISEC.
    !     8.  NMEV:  >= NUMBER OF T-MATRICES ALLOWED.
    !     9.  EXPLAN:  DIMENSION IS 20.
    !
    !     IMPORTANT TOLERANCES OR SCALES THAT ARE USED REPEATEDLY
    !     THROUGHOUT THIS PROGRAM ARE THE FOLLOWING:
    !     SCALED MACHINE EPSILON:  TTOL = TKMAX*EPSM WHERE
    !     EPSM = 2*MACHINE EPSILON AND
    !     TKMAX = MAX(|ALPHA(J)|,BETA(J), J = 1,MEV)
    !     BISEC CONVERGENCE TOLERANCE:  BISTOL = DSQRT(1000+MEV)*TTOL
    !     BISEC T-MULTIPLICITY TOLERANCE:  MULTOL = (1000+MEV)*TTOL
    !     LANCZOS CONVERGENCE TOLERANCE:   CONTOL = BETA(MEV+1)*1.D-10
    !-----------------------------------------------------------------------
    !     OUTPUT HEADER
    IF(.NOT.silent_) WRITE(6,10)
    10 FORMAT(/' LANCZOS PROCEDURE FOR REAL SYMMETRIC MATRICES'/)

    !     SET PROGRAM PARAMETERS
    !     SCALEK ARE USED IN TOLERANCES NEEDED IN SUBROUTINES LUMP,
    !     ISOEV AND PRTEST.  USER MUST NOT MODIFY THESE SCALES.

    IBIST  = 0 ! DON'T STOP AFTER BISEC
    SCALE1 = 5.0D2
    SCALE2 = 5.0D0
    SCALE3 = 5.0D0
    SCALE4 = 1.0D4
    ONE    = 1.0D0
    ZERO   = 0.0D0
    BTOL   = 1.0D-8
    GAPTOL = 1.0D-8
    ICONV  = 0
    MOLD   = 0
    MOLD1  = 1
    ICT    = 0
    MMB    = 0
    IPROJ  = 0

    !-----------------------------------------------------------------------
    !     READ USER-SPECIFIED PARAMETERS FROM INPUT FILE 5 (FREE FORMAT)
 
    !     READ USER-PROVIDED HEADER FOR RUN
    READ(UNIT5_,20) EXPLAN
    IF(.NOT.silent_) WRITE(6,20) EXPLAN
    READ(UNIT5_,20) EXPLAN
    IF(.NOT.silent_) WRITE(6,20) EXPLAN
    20 FORMAT(20A4)

    !     READ ORDER OF MATRICES (N) , MAXIMUM ORDER OF T-MATRIX (KMAX),
    !     NUMBER OF T-MATRICES ALLOWED (NMEVS), AND MATRIX IDENTIFICATION
    !     NUMBERS (MATNO)

    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) N,KMAX,NMEVS,MATNO
  
    ! HERE IS THE PLACE TO ALLOCATE SOME ARRAYS
 
    ALLOCATE(ALPHA(KMAX),BETA(KMAX+1))
    ALLOCATE(V1(MAX(N,KMAX+1)),VS(KMAX))
    ALLOCATE(G(MAX(N,KMAX*2)))
    ALLOCATE(MP(KMAX))
    ALLOCATE(NMEV(NMEVS))

    !     READ SEEDS FOR LANCZS AND INVERR SUBROUTINES (SVSEED AND RHSEED)
    !     READ MAXIMUM NUMBER OF ITERATIONS ALLOWED FOR EACH INVERSE
    !     ITERATION (MXINIT) AND MAXIMUM NUMBER OF STURM SEQUENCES
    !     ALLOWED (MXSTUR)

    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) seedsize

    ! HERE IS THE PLACE TO ALLOCATE SEED ARRAYS

    ALLOCATE(SVSEED(seedsize),RHSEED(seedsize),IIX(seedsize))
    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) SVSEED
    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) RHSEED
    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) MXINIT,MXSTUR
    !
    !     ISTART = (0,1):  ISTART = 0 MEANS ALPHA/BETA FILE IS NOT
    !     AVAILABLE.  ISTART = 1 MEANS ALPHA/BETA FILE IS AVAILABLE ON
    !     FILE 2.
    !     ISTOP = (0,1):  ISTOP = 0 MEANS PROCEDURE GENERATES ALPHA/BETA
    !     FILE AND THEN TERMINATES.  ISTOP = 1 MEANS PROCEDURE GENERATES
    !     ALPHAS/BETAS IF NEEDED AND THEN COMPUTES EIGENVALUES AND ERROR
    !     ESTIMATES AND THEN TERMINATES.

    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) ISTART,ISTOP

    !     IHIS = (0,1):  IHIS = 0 MEANS ALPHA/BETA FILE IS NOT WRITTEN
    !     TO FILE 1.  IHIS = 1 MEANS ALPHA/BETA FILE IS WRITTEN TO FILE 1.
    !     IDIST = (0,1):  IDIST = 0 MEANS DISTINCT T-EIGENVALUES
    !     ARE NOT WRITTEN TO FILE 11.  IDIST = 1 MEANS DISTINCT
    !     T-EIGENVALUES ARE WRITTEN TO FILE 11.
    !     IWRITE = (0,1):  IWRITE = 0 MEANS NO INTERMEDIATE OUTPUT
    !     FROM THE COMPUTATIONS IS WRITTEN TO FILE 6.  IWRITE = 1 MEANS
    !     T-EIGENVALUES AND ERROR ESTIMATES ARE WRITTEN TO FILE 6
    !     AS THEY ARE COMPUTED.

    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) IHIS,IDIST,IWRITE
 
    !     READ IN THE RELATIVE TOLERANCE (RELTOL) FOR USE IN THE
    !     SPURIOUS, T-MULTIPLICITY, AND PRTESTS.
    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) RELTOL

    !     READ IN THE SIZES OF THE T-MATRICES TO BE CONSIDERED.
    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) (NMEV(J), J=1,NMEVS)
    !
    !     READ IN THE NUMBER OF SUBINTERVALS TO BE CONSIDERED.
    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) NINT
    !
    ! HERE IS THE PLACE TO ALLOCATE LB,UB ARRAYS
    !
    ALLOCATE(LB(NINT),UB(NINT))
    !
    !     READ IN THE LEFT-END POINTS OF THE SUBINTERVALS TO BE CONSIDERED.
    !     THESE MUST BE IN ALGEBRAICALLY-INCREASING ORDER
    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) (LB(J), J=1,NINT)
    !
    !     READ IN THE RIGHT-END POINTS OF THE SUBINTERVALS TO BE CONSIDERED.
    !     THESE MUST BE IN ALGEBRAICALLY-INCREASING ORDER
    READ(UNIT5_,20) EXPLAN
    READ(UNIT5_,*) (UB(J), J=1,NINT)
 
    !-----------------------------------------------------------------------
    !     WRITE TO FILE 6, A SUMMARY OF THE PARAMETERS FOR THIS RUN

    IF(.NOT.silent_) WRITE(6,30) MATNO,N,KMAX
    30 FORMAT(/3X,'MATRIX ID',4X,'ORDER OF A',4X,'MAX ORDER OF T'/I12,I14,I18/)
    !
    IF(.NOT.silent_) WRITE(6,40) ISTART,ISTOP
    40 FORMAT(/2X,'ISTART',3X,'ISTOP'/2I8/)
    !
    IF(.NOT.silent_) WRITE(6,50) IHIS,IDIST,IWRITE
    50 FORMAT(/4X,'IHIS',3X,'IDIST',2X,'IWRITE'/3I8/)
    !
 
    IF(.NOT.silent_) WRITE(6,60)
    60 FORMAT(/' SEEDS FOR RANDOM NUMBER GENERATOR'//4X,'LANCZS SEED'/)
    IF(.NOT.silent_) WRITE(6,*) SVSEED
    IF(.NOT.silent_) WRITE(6,61)
    61 FORMAT(/4X,'INVERR SEED'/)
    IF(.NOT.silent_) WRITE(6,*) RHSEED
    !
    IF(.NOT.silent_) WRITE(6,70) (NMEV(J), J=1,NMEVS)
    70 FORMAT(/' SIZES OF T-MATRICES TO BE CONSIDERED'/(6I12))
    !
    IF(.NOT.silent_) WRITE(6,80) RELTOL,GAPTOL,BTOL
    80 FORMAT(/' RELATIVE TOLERANCE USED TO COMBINE COMPUTED T-EIGENVALUES'/E15.3/&
    ' RELATIVE GAP TOLERANCES USED IN INVERSE ITERATION'/E15.3/&
    ' RELATIVE TOLERANCE FOR CHECK ON SIZE OF BETAS'/E15.3/)
    !
    IF(.NOT.silent_) WRITE(6,90) (J,LB(J),UB(J), J=1,NINT)
    90 FORMAT(/' BISEC WILL BE USED ON THE FOLLOWING INTERVALS'/(I6,2E20.6))
    !
    IF (ISTART.EQ.0) GO TO 140
    !
    !     READ IN ALPHA BETA HISTORY
    !
    READ(UNIT2_,100)MOLD,NOLD,SVSOLD,MATOLD
    100 FORMAT(2I6,I12,I8)
    !
    IF (KMAX.LT.MOLD) KMAX = MOLD
    KMAX1 = KMAX + 1
    !
    !     CHECK THAT ORDER N, MATRIX ID MATNO, AND RANDOM SEED SVSEED
    !     AGREE WITH THOSE IN THE HISTORY FILE.  IF NOT PROCEDURE STOPS.
    !
    !ITEMP = (NOLD-N)**2+(MATNO-MATOLD)**2+(SVSEED-SVSOLD)**2
    ITEMP = (NOLD-N)**2+(MATNO-MATOLD)**2+SUM((SVSEED-SVSOLD)**2)
    !
    IF (ITEMP.EQ.0) GO TO 120
    !
    WRITE(6,110)
    110 FORMAT(' PROGRAM TERMINATES'/    ' READ FROM FILE 2 CORRESPONDS TO DIFFERENT MATRIX THAN MATRIX SPECIFIED'/)
    GO TO 640
    !
    120 CONTINUE
    MOLD1 = MOLD+1
    !
    READ(UNIT2_,130)(ALPHA(J), J=1,MOLD)
    READ(UNIT2_,130)(BETA(J), J=1,MOLD1)
    130 FORMAT(4Z20)
    !
    IF (KMAX.EQ.MOLD) GO TO 160
    !
    READ(UNIT2_,130)(V1(J), J=1,N)
    READ(UNIT2_,130)(V2(J), J=1,N)
 
    140 CONTINUE
    IIX = SVSEED

    !-----------------------------------------------------------------------
    CALL LANCZS_R(MATVEC,ALPHA,BETA,V1,V2,G,KMAX,MOLD1,N,IIX)
    !-----------------------------------------------------------------------
 
    KMAX1 = KMAX + 1
 
    IF (IHIS.EQ.0.AND.ISTOP.GT.0) GO TO 160
    !
    WRITE(UNIT1_,150) KMAX,N,MATNO
    150 FORMAT(3I12,' = KMAX,N,MATNO')
    WRITE(UNIT1_,151) seedsize
    151 FORMAT(I12,' = seedsize')
    WRITE(UNIT1_,152)
    152 FORMAT(' SVSEED = ')
    WRITE(UNIT1_,*) SVSEED
    !
    WRITE(UNIT1_,130)(ALPHA(I), I=1,KMAX)
    WRITE(UNIT1_,130)(BETA(I), I=1,KMAX1)
    !
    WRITE(UNIT1_,130)(V1(I), I=1,N)
    WRITE(UNIT1_,130)(V2(I), I=1,N)
    !
    IF (ISTOP.EQ.0) GO TO 540
    !
    160 CONTINUE
    BKMIN = BTOL
    !
    !----------------------------------------------------------------------
    IF(.NOT.silent_) WRITE(6,170)
    170 FORMAT(/' T-MATRICES (ALPHA AND BETA) ARE NOW AVAILABLE'/)
    !-----------------------------------------------------------------------
    !     SUBROUTINE TNORM CHECKS MIN(BETA)/(ESTIMATED NORM(A)) > BTOL .
    !     IF THIS IS VIOLATED IB IS SET EQUAL TO THE NEGATIVE OF THE INDEX
    !     OF THE MINIMAL BETA.  IF(IB < 0) THEN SUBROUTINE TNORM IS
    !     CALLED FOR EACH VALUE OF MEV TO DETERMINE WHETHER OR NOT THERE
    !     IS A BETA IN THE T-MATRIX SPECIFIED THAT VIOLATES THIS TEST.
    !     IF THERE IS SUCH A BETA THE PROGRAM TERMINATES FOR THE USER
    !     TO DECIDE WHAT TO DO.  THIS TEST CAN BE OVER-RIDDEN BY
    !     SIMPLY MAKING BTOL SMALLER, BUT THEN THERE IS THE POSSIBILITY
    !     THAT LOSSES IN THE LOCAL ORTHOGONALITY MAY HURT THE COMPUTATIONS.
    !     BTOL = 1.D-8 IS HOWEVER A CONSERVATIVE CHOICE FOR BTOL.
    !
    !     TNORM ALSO COMPUTES TKMAX = MAX(|ALPHA(K)|,BETA(K), K=1,KMAX).
    !     TKMAX IS USED TO SCALE THE TOLERANCES USED IN THE
    !     T-MULTIPLICITY AND SPURIOUS TESTS IN BISEC. TKMAX IS ALSO USED IN
    !     THE PROJECTION TEST FOR HIDDEN EIGENVALUES THAT HAD 'TOO SMALL'
    !     A PROJECTION ON THE STARTING VECTOR.
    !
    CALL TNORM(ALPHA,BETA,BKMIN,TKMAX,KMAX,IB,SILENT=silent_)
    !-----------------------------------------------------------------------
 
    TTOL = EPSM*TKMAX
 
    !     LOOP ON THE SIZE OF THE T-MATRIX

    180 CONTINUE
    MMB = MMB + 1
    MEV = NMEV(MMB)

    if(KMAX<MEV)then
     write(*,*) 'cullum hleval error, kmax/2 should be bigger than mev'
     write(*,*) 'mev    : ', MEV
     write(*,*) 'Kmax/2 : ', KMAX/2
     stop 'critical'
    endif

    !     IS MEV TOO LARGE ?
    IF(MEV.LE.KMAX) GO TO 200
    WRITE(6,190) MMB, MEV, KMAX
    190 FORMAT(/' TERMINATE PRIOR TO CONSIDERING THE',I6,'TH T-MATRIX'/&
    ' BECAUSE THE SIZE REQUESTED',I6,' IS GREATER THAN THE MAXIMUM SIZE ALLOWED',I6/)
    GO TO 540
    !
    200 MP1 = MEV + 1
    BETAM = BETA(MP1)
    !
    IF (IB.GE.0) GO TO 210
    !
    T0 = BTOL
    !
    !-----------------------------------------------------------------------
    CALL TNORM(ALPHA,BETA,T0,T1,MEV,IBMEV,SILENT=silent_)
    !-----------------------------------------------------------------------
 
    TEMP = T0/TKMAX
    IBMEV = IABS(IBMEV)
    IF (TEMP.GE.BTOL) GO TO 210
    IBMEV = -IBMEV
    GO TO 600

    210 CONTINUE
    IC = MXSTUR-ICT
 
    !-----------------------------------------------------------------------
    !     BISEC LOOP. THE SUBROUTINE BISEC INCORPORATES DIRECTLY THE
    !     T-MULTIPLICITY AND SPURIOUS TESTS. T-EIGENVALUES WILL BE
    !     CALCULATED BY BISEC SEQUENTIALLY ON INTERVALS
    !     (LB(J),UB(J)), J = 1,NINT).
    !
    !     ON RETURN FROM BISEC
    !     NDIS = NUMBER OF DISTINCT EIGENVALUES OF T(1,MEV) ON UNION
    !            OF THE (LB,UB) INTERVALS
    !     VS = DISTINCT T-EIGENVALUES IN ALGEBRAICALLY INCREASING ORDER
    !     MP = MULTIPLICITIES OF THE T-EIGENVALUES IN VS
    !     MP(I) = (0,1,MI), MI>1, I=1,NDIS  MEANS:
    !        (0)  VS(I) IS SPURIOUS
    !        (1)  VS(I) IS T-SIMPLE AND GOOD
    !        (MI) VS(I) IS MULTIPLE AND IS THEREFORE NOT ONLY GOOD BUT
    !             ALSO A CONVERGED GOOD T-EIGENVALUE.
    !
    !
    CALL BISEC(ALPHA,BETA,V1,V2,VS,LB,UB,EPSM,TTOL,MP,NINT,MEV,NDIS,IC,IWRITE,SILENT=silent_)
    IF (IBIST.EQ.1) GO TO 640
    !
    !-----------------------------------------------------------------------
    !
    IF (NDIS.EQ.0) GO TO 620
    !
    !     COMPUTE THE TOTAL NUMBER OF STURM SEQUENCES USED TO DATE
    !     COMPUTE THE BISEC CONVERGENCE AND T-MULTIPLICITY TOLERANCES USED.
    !     COMPUTE THE CONVERGENCE TOLERANCE FOR EIGENVALUES OF A.
    ICT    = ICT + IC
    TEMP   = DBLE(MEV+1000)
    MULTOL = TEMP*TTOL
    TEMP   = DSQRT(TEMP)
    BISTOL = TTOL*TEMP
    CONTOL = BETAM*1.D-10
    !
    !-----------------------------------------------------------------------
    !    SUBROUTINE LUMP 'COMBINES' T-EIGENVALUES THAT ARE 'TOO CLOSE'.
    !    NOTE HOWEVER THAT CLOSE SPURIOUS T-EIGENVALUES ARE NOT AVERAGED
    !    WITH GOOD ONES. HOWEVER, THEY MAY BE USED TO INCREASE THE
    !    MULTIPLICITY OF A GOOD T-EIGENVALUE.
    !
    LOOP = NDIS
    CALL LUMP(VS,RELTOL,MULTOL,SCALE2,MP,LOOP,SILENT=silent_)
    !
    !-----------------------------------------------------------------------
    IF(NDIS.EQ.LOOP) GO TO 230
    IF(.NOT.silent_) WRITE(6,220) NDIS, MEV, LOOP
    220 FORMAT(/I6,' DISTINCT T-EIGENVALUES WERE COMPUTED IN BISEC AT MEV',I6/&
    2X,' LUMP SUBROUTINE REDUCES NUMBER OF DISTINCT EIGENVALUES TO',I6)
 
    230 CONTINUE
    NDIS = LOOP
    BETA(MP1) = BETAM
    !-----------------------------------------------------------------------
    !     THE SUBROUTINE ISOEV LABELS THOSE SIMPLE EIGENVALUES OF T(1,MEV)
    !     WITH VERY SMALL GAPS BETWEEN NEIGHBORING EIGENVALUES OF T(1,MEV)
    !     TO AVOID COMPUTING ERROR ESTIMATES FOR ANY SIMPLE GOOD
    !     T-EIGENVALUE THAT IS TOO CLOSE TO A SPURIOUS EIGENVALUE.
    !     ON RETURN FROM ISOEV, G CONTAINS CODED MINIMAL GAPS
    !     BETWEEN THE DISTINCT EIGENVALUES OF T(1,MEV). (G IS REAL).
    !     G(I) < 0 MEANS MINGAP IS DUE TO LEFT GAP G(I) > 0 MEANS DUE TO
    !     RIGHT GAP. MP(I) = -1 MEANS THAT THE GOOD T-EIGENVALUE IS SIMPLE
    !     AND HAS A VERY SMALL MINGAP IN T(1,MEV) DUE TO A SPURIOUS
    !     EIGENVALUE.  NG = NUMBER OF GOOD T-EIGENVALUES.
    !     NISO = NUMBER OF ISOLATED GOOD T-EIGENVALUES.
    !
    CALL ISOEV(VS,GAPTOL,MULTOL,SCALE1,G,MP,NDIS,NG,NISO,SILENT=silent_)
    !-----------------------------------------------------------------------
    !
    IF(.NOT.silent_) WRITE(6,240)NG,NISO,NDIS
    240 FORMAT(/I6,' GOOD T-EIGENVALUES HAVE BEEN COMPUTED'/&
    I6,' OF THESE ARE T-ISOLATED'/&
    I6,' = NUMBER OF DISTINCT T-EIGENVALUES COMPUTED'/)
    !
    !     DO WE WRITE DISTINCT EIGENVALUES OF T-MATRIX TO FILE 11?
    IF (IDIST.EQ.0) GO TO 280
    !
    WRITE(11,250) NDIS,NISO,MEV,N,SVSEED,MATNO
    250 FORMAT(/4I6,I12,I8,' = NDIS,NISO,MEV,N,SVSEED,MATNO'/)
    !
    WRITE(11,260) (MP(I),VS(I),G(I), I=1,NDIS)
    260 FORMAT(2(I3,E25.16,E12.3))
    !
    WRITE(11,270) NDIS, (MP(I), I=1,NDIS)
    270 FORMAT(/I6,' = NDIS, T-MULTIPLICITIES (0 MEANS  SPURIOUS)'/(20I4))
    !
    280 CONTINUE
    !
    IF (NISO.NE.0) GO TO 310
    !
    if(rank==0) WRITE(UNIT4_,290) MEV
    290 FORMAT(/' AT MEV = ',I6,' THERE ARE NO ISOLATED T-EIGENVALUES'/&
    ' SO NO ERROR ESTIMATES WERE COMPUTED/')
    !
    IF(.NOT.silent_) WRITE(6,300)
    300 FORMAT(/' ALL COMPUTED GOOD T-EIGENVALUES ARE MULTIPLE'/&
    ' THEREFORE ALL SUCH EIGENVALUES ARE ASSUMED TO HAVE CONVERGED')
    !
    ICONV = 1
    GO TO 350
    !
    310 CONTINUE
    !
    !-----------------------------------------------------------------------
    !     SUBROUTINE INVERR COMPUTES ERROR ESTIMATES FOR ISOLATED GOOD
    !     T-EIGENVALUES USING INVERSE ITERATION ON T(1,MEV). ON RETURN
    !     G(J) = MINIMUM GAP IN T(1,MEV) FOR EACH VS(J), J=1,NDIS
    !     G(MEV+I) = BETAM*|U(MEV)| = ERROR ESTIMATE FOR ISOLATED GOOD
    !              T-EIGENVALUES, WHERE I = 1, NISO AND  BETAM = BETA(MEV+1)
    !              U(MEV) IS MEVTH COMPONENT OF THE UNIT EIGENVECTOR OF T
    !              CORRESPONDING TO THE ITH ISOLATED GOOD T-EIGENVALUE.
    !     A NEGATIVE ERROR ESTIMATE MEANS THAT FOR THAT PARTICULAR
    !     EIGENVALUE THE INVERSE ITERATION DID NOT CONVERGE IN <= MXINIT
    !     STEPS AND THAT THE CORRESPONDING ERROR ESTIMATE IS QUESTIONABLE.
    !
    !     V2 CONTAINS THE ISOLATED GOOD T-EIGENVALUES
    !     V1 CONTAINS THE MINGAPS TO THE NEAREST DISTINCT  EIGENVALUE
    !        OF T(1,MEV) FOR EACH ISOLATED GOOD T-EIGENVALUE IN V2.
    !     VS CONTAINS THE NDIS DISTINCT EIGENVALUES OF T(1,MEV)
    !     MP CONTAINS THE CORRESPONDING CODED T-MULTIPLICITIES
    !
    IWRITE = 0
    IT = MXINIT
    CALL INVERR(ALPHA,BETA,V1,V2,VS,EPSM,G,MP,MEV,MMB,NDIS,NISO,N,RHSEED,IT,IWRITE,SILENT=silent_,UNIT4=UNIT4_)
    !
    !-----------------------------------------------------------------------
    !
    !     SIMPLE CHECK FOR CONVERGENCE. CHECKS TO SEE IF ALL OF THE ERROR
    !     ESTIMATES ARE SMALLER THAN CONTOL = BETAM*1.D-10.
    !     IF THIS TEST IS SATISFIED, THEN CONVERGENCE FLAG, ICONV IS SET
    !     TO 1.  TYPICALLY ERROR ESTIMATES ARE VERY CONSERVATIVE.
    !
    IF(.NOT.silent_) WRITE(6,320) CONTOL
    320 FORMAT(/' CONVERGENCE IS TESTED USING THE CONVERGENCE TOLERANCE',E13.4/)
    !
    II = MEV +1
    IFF = MEV+NISO
    DO 330 I = II,IFF
      IF (ABS(G(I)).GT.CONTOL) GO TO 350
      330 CONTINUE
      ICONV = 1
      MMB = NMEVS
      !
      IF(.NOT.silent_) WRITE(6,340) CONTOL
      340 FORMAT(' ALL COMPUTED ERROR ESTIMATES WERE LESS THAN',E15.4/' THEREFORE PROCEDURE TERMINATES'/)
      !
      350 CONTINUE
      !
      !     IF CONVERGENCE IS INDICATED, THAT IS ICONV = 1 ,THEN
      !     THE SUBROUTINE PRTEST IS CALLED TO CHECK FOR ANY CONVERGED
      !     T-EIGENVALUES THAT HAVE BEEN MISLABELLED AS SPURIOUS BECAUSE
      !     THE PROJECTION OF THEIR EIGENVECTOR(S) ON THE STARTING
      !     VECTOR WERE TOO SMALL.
      !     NUMERICAL TESTS INDICATE THAT SUCH EIGENVALUES ARE RARE.
      !     IF FOR SOME REASON MANY OF THESE HIDDEN EIGENVALUES APPEAR
      !     ON SOME RUN, YOU CAN BE CERTAIN THAT SOMETHING IS FOULED UP.
      !
      IF (ICONV.EQ.0) GO TO 480
      !
      !-----------------------------------------------------------------------
      !
      CALL PRTEST (ALPHA,BETA,VS,TKMAX,EPSM,RELTOL,SCALE3,SCALE4,MP,NDIS,MEV,IPROJ,SILENT=silent_)
      !
      !-----------------------------------------------------------------------
      !
      IF(IPROJ.EQ.0) GO TO 470
      !
      IF(IDIST.EQ.1)  WRITE(11,360) IPROJ
      360 FORMAT(' SUBROUTINE PRTEST WANTS TO RELABEL',I6,' SPURIOUS T-EIGENVALUES'/&
      ' WE ACCEPT RELABELLING ONLY IF LAST COMPONENT OF T-EIGENVECTOR IS L.T. 1.D-10'/)
      !
      IIX = RHSEED
      !
      !-----------------------------------------------------------------------
      !
      CALL GENRAN(IIX,G,MEV)
      !
      !-----------------------------------------------------------------------
      !
      ITEN = -10
      NISOM = NISO + MEV
      IWRITO = IWRITE
      IWRITE = 0
      !
      DO 390 J = 1,NDIS
 IF(MP(J).NE.ITEN) GO TO 390
 T0 = VS(J)
 !
 !-----------------------------------------------------------------------
 !
 IT = MXINIT
 CALL INVERM(ALPHA,BETA,V1,V2,T0,TEMP,T1,EPSM,G,MEV,IT,IWRITE,SILENT=silent_)
 !
 !-----------------------------------------------------------------------
 !
 IF(TEMP.LE.1.D-10) GO TO 380
 !     ERROR ESTIMATE WAS NOT SMALL REJECT RELABELLING OF THIS EIGENVALUE
 IF(IDIST.EQ.1)  WRITE(11,370) J,T0,TEMP
 370 FORMAT(/' LAST COMPONENT FOR',I6,'TH EIGENVALUE',E20.12/&
 ' IS TOO LARGE = ',E15.6,' SO DO NOT ACCEPT PRTEST RELABELLING'/)
 MP(J) = 0
 IPROJ = IPROJ - 1
 GO TO 390
 !     RELABELLING ACCEPTED
 380 NISOM = NISOM + 1
 G(NISOM) = BETAM*TEMP
 390 CONTINUE
 IWRITE = IWRITO
 !
 IF(IPROJ.EQ.0) GO TO 430
 IF(.NOT.silent_) WRITE(6,400) IPROJ
 400 FORMAT(/I6,' T-EIGENVALUES WERE RECLASSIFIED AS GOOD.'/&
 ' THESE ARE IDENTIFIED IN FILE 3 BY A T-MULTIPLICITY OF -10'/&
 ' USER SHOULD INSPECT EACH TO MAKE SURE NEIGHBORS HAVE CONVERGED'/)
 !
 IF(IDIST.EQ.1)  WRITE(11,410) IPROJ
 410 FORMAT(/I6,' T-EIGENVALUES WERE RELABELLED AS GOOD'/&
 ' BELOW IS CORRECTED T-MULTIPLICITY PATTERN'/)
 !
 IF(.NOT.silent_) WRITE(6,420) NDIS, (MP(I), I=1,NDIS)
 IF(IDIST.EQ.1)  WRITE(11,420) NDIS, (MP(I), I=1,NDIS)
 420 FORMAT(/I6,' = NDIS, T-MULTIPLICITIES (0 MEANS  SPURIOUS)'/&
 6X, ' (-10) MEANS SPURIOUS T-EIGENVALUE RELABELLED AS GOOD'/(20I4))
 !
 !     RECALCULATE MINGAPS FOR DISTINCT T(1,MEV) EIGENVALUES.
 430 NM1 = NDIS - 1
 G(NDIS) = VS(NM1)-VS(NDIS)
 G(1) = VS(2)-VS(1)
 !
 DO 440 J = 2,NM1
   T0 = VS(J)-VS(J-1)
   T1 = VS(J+1)-VS(J)
   G(J) = T1
   IF (T0.LT.T1) G(J) = -T0
   440 CONTINUE
   IF(IPROJ.EQ.0) GO TO 470
   !     WRITE TO FILE 4 ERROR ESTIMATES FOR THOSE T-EIGENVALUES RELABELLED
   NGOOD = 0
   DO 450 J = 1,NDIS
     IF(MP(J).EQ.0) GO TO 450
     NGOOD = NGOOD + 1
     IF(MP(J).NE.ITEN) GO TO 450
     T0 = VS(J)
     NISO = NISO + 1
     NISOM = MEV + NISO
     if(rank==0) WRITE(UNIT4_,460) NGOOD,T0,G(NISOM),G(J)
     450 CONTINUE
     460 FORMAT(I10,E25.16,2E14.3)
     !
     470 CONTINUE
     !
     !     WRITE THE GOOD T-EIGENVALUES TO FILE 3.  FIRST TRANSFER THEM
     !     TO V2 AND THEIR T-MULTIPLICITIES TO THE CORRESPONDING POSITIONS
     !     IN MP AND COMPUTE THE A-MINGAPS, THE MINIMAL GAPS BETWEEN THE
     !     GOOD T-EIGENVALUES.  THESE GAPS WILL BE PUT IN THE ARRAY G.
     !     SINCE G CURRENTLY CONTAINS THE MINIMAL GAPS BETWEEN THE DISTINCT
     !     EIGENVALUES OF THE T-MATRIX, THESE GAPS WILL FIRST BE
     !     TRANSFERRED TO V1.  NOTE THAT V1<0 MEANS THAT THAT MINIMAL GAP
     !     IN THE T-MATRIX IS DUE TO A SPURIOUS T-EIGENVALUE.
     !     ALL THIS INFORMATION IS PRINTED TO FILE 3
     !
     480 CONTINUE
     !
     NG = 0
     DO 490 I = 1,NDIS
       IF (MP(I).EQ.0) GO TO 490
       NG = NG+1
       MP(NG) = MP(I)
       V2(NG) = VS(I)
       TEMP = G(I)
       TEMP = DABS(TEMP)
       J = I+1
       IF (G(I).LT.ZERO) J = I-1
       IF (MP(J).EQ.0) TEMP = -TEMP
       V1(NG) = TEMP
       490 CONTINUE
       !
       IF(.NOT.silent_) WRITE(6,500)MEV
       500 FORMAT(//' T-EIGENVALUE CALCULATION AT MEV = ',I6,'    IS COMPLETE')
       !
       !     NG = NUMBER OF COMPUTED DISTINCT GOOD T-EIGENVALUES.  NEXT
       !     GENERATE GAPS BETWEEN GOOD T-EIGENVALUES (AMINGAPS) AND PUT THEM
       !     IN G.  G(J) < 0 MEANS THE AMINGAP IS DUE TO THE LEFT-HAND GAP.
       !
       NGM1 = NG - 1
       G(NG) = V2(NGM1)-V2(NG)
       G(1) = V2(2)-V2(1)
       !
       DO 510 J = 2,NGM1
  T0 = V2(J)-V2(J-1)
  T1 = V2(J+1)-V2(J)
  G(J) = T1
  IF (T0.LT.T1) G(J) = -T0
  510 CONTINUE
  !
  !     WRITE GOOD T-EIGENVALUES OUT TO FILE 3.
  !
  !if(rank==0) WRITE(UNIT3_,520)NG,NDIS,MEV,N,SVSEED,MATNO,MULTOL,IB,BTOL
  !520 FORMAT(4I6,I12,I8,' = NG,NDIS,MEV,N,SVEED,MATNO'/&
  !E20.12,I6,E13.4,' = MUTOL,INDEX MINIMAL BETA,BTOL'/&
  !' EV NO',1X,'TMULT',10X,'GOOD EIGENVALUE',7X,'TMINGAP',7X,'AMINGAP')
  if(rank==0) WRITE(UNIT3_,520)NG,NDIS,MEV,N,MATNO
  520 FORMAT(3I6,I12,I8,' = NG,NDIS,MEV,N,MATNO')
  if(rank==0) WRITE(UNIT3_,521) SIZE(SVSEED)
  521 FORMAT(I12,' = seedsize')
  if(rank==0) WRITE(UNIT3_,522)
  522 FORMAT(' SVSEED')
  if(rank==0) WRITE(UNIT3_,*) SVSEED
  if(rank==0) WRITE(UNIT3_,523) MULTOL,IB,BTOL
  523 FORMAT(E20.12,I6,E13.4,' = MUTOL,INDEX MINIMAL BETA,BTOL'/&
  ' EV NO',2X,'TMULT',7X,'GOOD T-EIGENVALUE',7X,'TMINGAP',7X,'AMINGAP')
  !
  if(rank==0) WRITE(UNIT3_,530)(I,MP(I),V2(I),V1(I),G(I), I=1,NG)
  530 FORMAT(2I6,E25.16,2E14.3)
  !
  !     IF CONVERGENCE FLAG ICONV.NE.1 AND NUMBER OF T-MATRICES
  !     CONSIDERED TO DATE IS LESS THAN NUMBER ALLOWED, INCREMENT MEV.
  !     AND LOOP BACK TO 210 TO REPEAT COMPUTATIONS.  RESTORE BETA(MEV+1).
  !
  BETA(MP1) = BETAM
  !
  IF (MMB.LT.NMEVS.AND.ICONV.NE.1) GO TO 180
  !
  !     END OF LOOP ON DIFFERENT SIZE T-MATRICES ALLOWED.
  !
  540 CONTINUE
  !
  IF(ISTOP.EQ.0)THEN;  IF(.NOT.silent_) WRITE(6,550); ENDIF
  550 FORMAT(/' T-MATRICES (ALPHA AND BETA) ARE NOW AVAILABLE, TERMINATE')
  IF (IHIS.EQ.1.AND.KMAX.NE.MOLD) WRITE(UNIT1_,560)
  560 FORMAT(/' ABOVE ARE THE FOLLOWING VECTORS '/&
  '  ALPHA(I), I = 1,KMAX'/&
  '  BETA(I), I = 1,KMAX+1'/&
  ' FINAL TWO LANCZOS VECTORS OF ORDER N FOR I = KMAX,KMAX+1'/&
  ' ALL VECTORS IN THIS FILE HAVE HEX FORMAT 4Z20 '/&
  ' ----- END OF FILE 1 NEW ALPHA, BETA HISTORY---------------'///)
  !
  IF (ISTOP.EQ.0) GO TO 640
  !
  if(rank==0) WRITE(UNIT3_,570)
  570 FORMAT(/' ABOVE ARE COMPUTED GOOD T-EIGENVALUES'/&
  ' NG = NUMBER OF GOOD T-EIGENVALUES COMPUTED'/&
  ' NDIS = NUMBER OF COMPUTED DISTINCT EIGENVALUES OF T(1,MEV)'/&
  ' N = ORDER OF A,  MATNO = MATRIX IDENT'/&
  ' MULTOL = T-MULTIPLICITY TOLERANCE FOR T-EIGENVALUES IN BISEC'/&
  ' TMULT IS THE T-MULTIPLICITY OF GOOD T-EIGENVALUE'/&
  ' TMULT = -1 MEANS SPURIOUS T-EIGENVALUE TOO CLOSE'/&
  ' DO NOT COMPUTE ERROR ESTIMATES FOR SUCH EIGENVALUES'/&
  ' AMINGAP = MINIMAL GAP BETWEEN THE COMPUTED A-EIGENVALUES'/&
  ' AMINGAP .LT. 0. MEANS MINIMAL GAP IS DUE TO LEFT-HAND GAP'/&
  ' TMINGAP= MINIMAL GAP W.R.T.  DISTINCT EIGENVALUES IN T(1,MEV)'/&
  ' TMINGAP .LT. 0. MEANS MINGAP IS DUE TO SPURIOUS T-EIGENVALUE'/&
  ' ----- END OF FILE 3 GOODEIGENVALUES-----------------------'///)
  !
  IF (IDIST.EQ.1) WRITE(11,580)
  580 FORMAT(/' ABOVE ARE THE DISTINCT EIGENVALUES OF T(1,MEV).'/&
  ' THE FORMAT IS      T-MULTIPLICITY    T-EIGENVALUE   TMINGAP'/&
  '        THIS FORMAT IS REPEATED TWICE ON EACH LINE.'/&
  ' T-MULTIPLICITY = -1 MEANS THAT THE SUBROUTINE ISOEV HAS TAGGED'/&
  '   THIS SIMPLE T-EIGENVALUE AS HAVING A VERY CLOSE SPURIOUS'/&
  '   T-EIGENVALUE SO THAT NO ERROR ESTIMATE WILL BE COMPUTED'/&
  '   FOR THAT EIGENVALUE IN SUBROUTINE INVERR.'/&
  ' TMINGAP .LT. 0, TMINGAP IS DUE TO LEFT GAP .GT. 0, RIGHT GAP.'/&
  ' EACH OF THE DISTINCT T-EIGENVALUE TABLES IS FOLLOWED'/&
  ' BY THE T-MULTIPLICITY PATTERN.'/&
  ' NDIS = NUMBER OF COMPUTED DISTINCT EIGENVALUES OF T(1,MEV).'/&
  ' NG = NUMBER OF GOOD T-EIGENVALUES. '/&
  ' NISO = NUMBER OF ISOLATED GOOD T-EIGENVALUES. '/&
  ' NISO ALSO IS THE COUNT OF +1 ENTRIES IN T-MULTIPLICITY PATTERN.'/&
  ' ----- END OF FILE 11 DISTINCT T-EIGENVALUES--------------'///)
  !
  IF(NISO.NE.0.and.rank==0) WRITE(UNIT4_,590)
  590 FORMAT(/' ABOVE ARE THE ERROR ESTIMATES OBTAINED FOR THE ISOLATED GOOD T-EIGENVALUES'/&
  ' OBTAINED VIA INVERSE ITERATION IN THE SUBROUTINE INVERR.'/&
  ' ALL OTHER GOOD T-EIGENVALUES HAVE CONVERGED.'/&
  ' ERROR ESTIMATE = BETAM*ABS(UM)'/&
  ' WHERE BETAM = BETA(MEV+1) AND UM = U(MEV).'/&
  ' U = UNIT EIGENVECTOR OF T WHERE T*U = EV*U AND EV = ISOLATED GOOD T-EIGENVALUE.'/&
  ' TMINGAP = GAP TO NEAREST DISTINCT EIGENVALUE OF T(1,MEV).'/&
  ' TMINGAP .LT. 0. MEANS MINGAP IS DUE TO A LEFT NEIGHBOR.'/&
  ' ERROR ESTIMATE L.T. 0 MEANS INVERSE ITERATION DID NOT CONVERGE'/&
  ' ------ END OF FILE 4 ERRINV -------------------------------'//)
  GO TO 640
  !
  600 CONTINUE
  !
  IBB = IABS(IBMEV)
  IF (IBMEV.LT.0) WRITE(6,610) MEV,IBB,BETA(IBB)
  610 FORMAT(/' PROGRAM TERMINATES BECAUSE MEV REQUESTED = ',I6,' IS .GT',I6/&
  ' AT WHICH AN ABNORMALLY SMALL BETA = ' , E13.4,' OCCURRED'/)
  GO TO 640
  !
  620 IF (NDIS.EQ.0.AND.ISTOP.GT.0) WRITE(6,630)
  630 FORMAT(/' INTERVALS SPECIFIED FOR BISECT DID NOT CONTAIN ANY T-EIGENVALUES'/&
  ' PROGRAM TERMINATES')
  !
  640 CONTINUE
  !
  IF(ALLOCATED(V1))     DEALLOCATE(V1)
  !IF(ALLOCATED(V2))     DEALLOCATE(V2)
  IF(ALLOCATED(VS))     DEALLOCATE(VS)
  IF(ALLOCATED(ALPHA))  DEALLOCATE(ALPHA)
  IF(ALLOCATED(BETA))   DEALLOCATE(BETA)
  IF(ALLOCATED(G))      DEALLOCATE(G)
  IF(ALLOCATED(LB))     DEALLOCATE(LB)
  IF(ALLOCATED(UB))     DEALLOCATE(UB)
  IF(ALLOCATED(MP))     DEALLOCATE(MP)
  IF(ALLOCATED(NMEV))   DEALLOCATE(NMEV)
  IF(ALLOCATED(SVSEED)) DEALLOCATE(SVSEED)
  IF(ALLOCATED(SVSOLD)) DEALLOCATE(SVSOLD)
  IF(ALLOCATED(RHSEED)) DEALLOCATE(RHSEED)
  IF(ALLOCATED(IIX))    DEALLOCATE(IIX)
 
 
  IF(PRESENT(FILE1)) CLOSE(UNIT1_)
  IF(PRESENT(FILE3)) CLOSE(UNIT3_)
  IF(PRESENT(FILE4)) CLOSE(UNIT4_)
  IF(PRESENT(FILE5)) CLOSE(UNIT5_)

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
