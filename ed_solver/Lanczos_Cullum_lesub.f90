MODULE Lanczos_Cullum_lesub

 use Lanczos_Cullum_mult 

IMPLICIT NONE

! HLEVAL/LEVAL PARAMETERS

  INTEGER                      :: ISTART_,ISTOP_,IHIS_,IDIST_,IWRITE__,MATNO__
  INTEGER                      :: KMAX_,MXINIT_,MXSTUR_,NMEVS_,NINT_,MOLD1_
  REAL(8)                      :: RELTOL__
  INTEGER, ALLOCATABLE         :: NMEV_(:)
  REAL(8), ALLOCATABLE         :: LB_(:),UB_(:)

  REAL(DBL),PARAMETER, PRIVATE :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,  PARAMETER, PRIVATE :: F=.FALSE.,T=.TRUE.


CONTAINS

 
  !-----LESUB  (1) REAL SYMMETRIC-----------------------------------------
  !            (2) HERMITIAN MATRICES
  !            (3) FACTORED INVERSES OF REAL SYMMETRIC MATRICES AND
  !            (4) REAL SYMMETRIC GENERALIZED, A*X = EVAL*B*X WHERE
  !                B IS POSITIVE DEFINITE, CHOLESKY FACTOR AVAILABLE
  !
  !     SUBROUTINES    BISEC, INVERR, TNORM, LUMP, ISOEV, PRTEST, AND
  !                    INVERM ARE USED WITH THE LANCZOS EIGENVALUE
  !                    PROGRAMS LEVAL, HLEVAL, LIVAL AND LGVAL. STURMI,
  !                    INVERM, LBISEC, AND TNORM ARE USED WITH THE
  !                    EIGENVECTOR PROGRAMS LEVEC, HLEVEC, LIVEC AND
  !                    LGVEC. 
  !
  !-----COMPUTE T-EIGENVALUES BY BISECTION--------------------------------


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

  SUBROUTINE init_Lanczos_Cullum_leval
    INTEGER :: ievs,sub,meth,m,meth2

    ! SET HLEVAL PARAMETERS TO START FROM SCRATCH
    ISTART_ = 0
    ISTOP_  = 1
    IHIS_   = 1
    IDIST_  = 0
    IWRITE__= 1
    !
    MXINIT_ = 5
    MXSTUR_ = 100000
    MATNO__ = 1
    RELTOL__= tolerance
    !
    KMAX_    = 2*Nitermax

    meth  = 1
    meth2 = 1
    if(allocated(NMEV_)) deallocate(NMEV_)
    if(allocated(LB_)) deallocate(LB_,UB_)


  !############################################################!
   SELECT CASE (meth)
    CASE(1)
     NMEVS_   =  1
     ALLOCATE(NMEV_(1))
     NMEV_(1) =  KMAX_/4
    CASE(2)
     m        =  4
     NMEVS_   =  10
     ALLOCATE(NMEV_(NMEVS_))
     sub      =  INT(DBLE(KMAX_)/DBLE(m*NMEVS_))-1
     NMEV_    =  (/(sub*ievs,ievs=1,NMEVS_-1),KMAX_/m-1/)
   END SELECT

  !############################################################!

   SELECT CASE (meth2)
    CASE(1)
     NINT_    =  7
     ALLOCATE(LB_(NINT_),UB_(NINT_))
     LB_      =  (/ -500.d0, -50.d0,-10.d0,  0.d0, 10.d0,  50.d0,  500.d0 /)
     UB_      =  (/  -50.d0, -10.d0,  0.d0, 10.d0, 50.d0, 500.d0, 1000.d0 /)
    CASE(2)
     NINT_    =  2
     ALLOCATE(LB_(NINT_),UB_(NINT_))
     LB_      =  (/ -500.d0,    0.d0    /)
     UB_      =  (/    0.d0,  500.0d0   /)
   END SELECT


  !############################################################!

   if(maxval(NMEV_)>KMAX_/4)then
    write(*,*) 'NMEVS IS : ', NMEV_
    write(*,*) 'max nmev : ', maxval(NMEV_)
    write(*,*) 'KMAX     : ', KMAX_
    write(*,*) ' the relation KMAX/4>nmev should be satisfied  ' 
    stop
   endif 

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


  SUBROUTINE BISEC(ALPHA,BETA,BETA2,VB,VS,LBD,UBD,EPS,TTOL,MP,NINT,MEV,NDIS,IC,IWRITE,SILENT)
    LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
    LOGICAL                       :: silent_
    REAL(8)                       :: ALPHA(:),BETA(:),BETA2(:),VB(:),VS(:)
    REAL(8)                       :: LBD(:),UBD(:),EPS,EPT,EP0,EP1,TEMP,TTOL
    REAL(8)                       :: ZERO,ONE,HALF,YU,YV,LB,UB,XL,XU,X1,X0,XS,BETAM
    INTEGER                       :: MP(:),IDEF(10),KMAX,MOLD1,N,I,J,IN,MT,MTK,K,NU ,NEV1,IEST,MTI,II,JSTURM,KL,JL,IL,JU,KU 
    INTEGER                       :: NINT,MEV,NDIS,IC,MXSTUR,ISKIP,NA,MD,JIND,MP1,NG,ICT,ISTURM,NEV,IC0  
    INTEGER                       :: KEV,JEV,MPEV,KNE,ITER,KNEW,IWRITE 


                        silent_ = F
    IF(PRESENT(SILENT)) silent_ = SILENT


    !     COMPUTES EIGENVALUES OF T(1,MEV) BY LOOPING INTERNALLY ON THE
    !     USER-SPECIFIED INTERVALS,  (LB(J),UB(J)), J = 1,NINT.  INTERVALS
    !     ARE TREATED AS OPEN ON THE LEFT AND CLOSED ON THE RIGHT.
    !     THE BISEC SUBROUTINE SIMULTANEOUSLY LABELS SPURIOUS T-EIGENVALUES
    !     AND DETERMINES THE T-MULTIPLICITIES OF EACH GOOD T-EIGENVALUE.
    !     SPURIOUS T-EIGENVALUES ARE LABELLED BY A T-MULTIPLICITY = 0.
    !     ANY T-EIGENVALUE WITH A T-MULTIPLICITY >= 1 IS 'GOOD'.
    !
    !     IF IWRITE = 0 THEN MOST OF THE WRITES TO FILE 6 ARE NOT
    !     ACTIVATED.
    !
    !     NOTE THAT PROGRAM ASSUMES THAT NO MORE THAN MMAX/2 EIGENVALUES
    !     OF T(1,MEV) ARE TO BE COMPUTED IN ANY ONE OF THE SUBINTERVALS
    !     CONSIDERED, WHERE MMAX = DIMENSION OF VB SPECIFIED BY THE USER
    !     IN THE MAIN PROGRAM LEVAL.
    !
    !     ON ENTRY
    !     BETA2(J) IS SET = BETA(J)*BETA(J).  THE STORAGE FOR BETA2 COULD
    !     BE ELIMINATED BY RECOMPUTING THE BETA(J)**2 FOR EACH STURM
    !     SEQUENCE.
    !
    !     EPS = 2*MACHEP =  4.4 * 10**-16 ON IBM 3081.
    !     TTOL = EPS*TKMAX WHERE
    !     TKMAX = MAX(|ALPHA(K)|,BETA(K), K=1,KMAX)
    !
    !     ON EXIT
    !     NDIS = TOTAL NUMBER OF COMPUTED DISTINCT EIGENVALUES OF
    !            T(1,MEV) ON THE UNION OF THE (LB,UB) INTERVALS.
    !     VS = COMPUTED DISTINCT EIGENVALUES OF T(1,MEV) IN ALGEBRAICALLY-
    !          INCREASING ORDER
    !     MP = CORRESPONDING T-MULTIPLICITIES OF THESE EIGENVALUES
    !     MP(I) = (0,1,MI), MI>1, I=1,NDIS  MEANS:
    !        (0)  V(I) IS SPURIOUS
    !        (1)  V(I) IS ISOLATED AND GOOD
    !        (MI) V(I) IS MULTIPLE AND HENCE A CONVERGED GOOD T-EIGENVALUE.
    !     IC = TOTAL NUMBER OF STURMS USED
    !
    !     DEFAULTS
    !     ISKIP = 0 INITIALLY. IF DEFAULT OCCURS ON J-TH SUB-INTERVAL, SET
    !             ISKIP=ISKIP+1 AND IDEF(ISKIP) = J
    !             DEFAULTS OCCUR IF THERE ARE NO T-EIGENVALUES  IN THE
    !             SUBINTERVAL SPECIFIED OR IF THE NUMBER
    !             OF STURMS SEQUENCES REQUIRED EXCEEDS MXSTUR.
    !             WHEN A DEFAULT OCCURS THE PROGRAM
    !             SKIPS THE INTERVAL INVOLVED AND GOES ON TO THE NEXT
    !             INTERVAL.
    !
    !-----------------------------------------------------------------------

    !     SPECIFY PARAMETERS

    ZERO   =  0.0D0
    ONE    =  1.0D0
    HALF   =  0.5D0
    MXSTUR =  IC
    NDIS   =  0
    IC     =  0
    ISKIP  =  0
    MP1    =  MEV+1

    if(MP1>SIZE(BETA)) MP1=SIZE(BETA)

    !     SAVE THEN SET BETA(MEV+1) = 0. GENERATE BETA**2

    BETAM     = BETA(MP1)
    BETA(MP1) = ZERO
 
    DO 10 I = 1,MIN(SIZE(BETA2),MP1)
      10 BETA2(I) = BETA(I)*BETA(I)

      !     EP0 IS USED IN T-MULTIPLICITY AND SPURIOUS TESTS
      !     EP1 AND EPS ARE USED IN THE BISEC CONVERGENCE TEST

      TEMP = DBLE(MEV+1000)
      EP0  = TEMP*TTOL
      EP1  = DSQRT(TEMP)*TTOL

      IF(.NOT.silent_) WRITE(6,20)MEV,NINT
      20 FORMAT(/' BISEC CALCULATION'/' ORDER OF T IS',I6/&
      ' NUMBER OF INTERVALS IS',I6/)

      IF(.NOT.silent_) WRITE(6,30) EP0,EP1
      30 FORMAT(/' MULTOL, TOLERANCE USED IN T-MULTIPLICITY AND SPURIOUS TETS = ',E10.3/&
      ' BISTOL, TOLERANCE USED IN BISEC CONVERGENCE TEST =',E10.3/)

      !     LOOP ON THE NINT INTERVALS  (LB(J),UB(J)), J=1,NINT

      DO 430 JIND = 1,NINT
        LB = LBD(JIND)
        UB = UBD(JIND)

 IF(.NOT.silent_) WRITE(6,40)JIND,LB,UB
 40 FORMAT(//1X,'BISEC INTERVAL NO',2X,'LOWER BOUND',2X,'UPPER BOUND'/I18,2E13.5/)

 !     INITIALIZATION AND PARAMETER SPECIFICATION
 !     ICT IS TOTAL STURM COUNT ON (LB,UB)

 NA  = 0
 MD  = 0
 NG  = 0
 ICT = 0
 
 !     START OF T-EIGENVALUE CALCULATIONS
 X1 = UB
 ISTURM = 1
 GO TO 330
 !     FORWARD STURM CALCULATION TO DETERMINE NA = NO. T-EIGENVALUES > UB
 50 NA = NEV
 !
 X1     = LB
 ISTURM = 2
 GO TO 330
 !     FORWARD STURM CALC TO DETERMINE MT = NO. T-EIGENVALUES ON (LB,UB)
 60 CONTINUE
 MT=NEV

 if(MT>size(VB))then
  write(*,*) 'mt bigger then VS, mt vs :', MT,size(VB)
  stop
 endif

 ICT = ICT +2
 !
 IF(.NOT.silent_) WRITE(6,70) MT,NA
 70 FORMAT(/2I6,' = NO. TMEV ON (LB,UB) AND NO. .GT. UB'/)

 !     DEFAULT TEST: IS ESTIMATED NUMBER OF STURMS > MXSTUR?
 IEST = 30*MT
 IF (IEST.LT.MXSTUR) GO TO 90

 IF(.NOT.silent_) WRITE(6,80)
 80 FORMAT(//' ESTIMATED NUMBER OF STURMS REQUIRED EXCEEDS USER LIMIT'/' SKIP THIS SUBINTERVAL')
 GO TO 110

 90 CONTINUE
 
 IF (MT.GE.1) GO TO 120

 IF(.NOT.silent_) WRITE(6,100)
 100 FORMAT(//' THERE ARE NO T-EIGENVALUES ON THIS INTERVAL)'/)

 110 ISKIP = ISKIP+1
 IDEF(ISKIP) = JIND
 GO TO 430

 !     REGULAR CASE.
 120 CONTINUE

 IF (IWRITE.NE.0)THEN; IF(.NOT.silent_) WRITE(6,130); ENDIF
 130 FORMAT(/' DISTINCT T-EIGENVALUES COMPUTED USING BISEC'/13X,'T-EIGENVALUE',2X,'TMULT',3X,'MD',4X,'NG')

 !     SET UP INITIAL UPPER AND LOWER BOUNDS FOR T-EIGENVALUES

 DO 140 I=1,MT

   if(I>SIZE(VB)) then
    write(*,*) 'outside VB boundaries cullum, i vb :',I,size(VB)
    stop
   endif

   VB(I)       = LB
   MTI         = MT + I

   if(MTI>SIZE(VB)) then
     write(*,*) ' outside VB boundaries cullum, mti vb :',MTI,size(VB)
     write(*,*) ' mt i MTI        : ' , MT,I,MTI
     write(*,*) ' NEV             : ',  NEV
     write(*,*) ' 2*MT < size(VB) : ',  2*MT
     stop
   endif

   140 VB(MTI) = UB
  
   !     CALCULATE T-EIGENVALUES FROM LB UP TO UB  K = MT,...,1
   !     MAIN LOOP FOR FINDING KTH T-EIGENVALUE
  
   K      = MT
   150 CONTINUE
   IC0    = 0
   XL     = VB(K)
   MTK    = MT+K
   XU     = VB(MTK)
   ISTURM = 3
   X1     = XU
   IC0    = IC0 + 1
   GO TO 330

   !     FORWARD STURM CALCULATION AT XU
   160 NU=NEV

   !     BISECTION LOOP FOR KTH T-EIGENVALUE. TEST  X1=MIDPOINT OF (XL,XU)

   ISTURM = 4
   170 CONTINUE
   X1     = (XL+XU)*HALF
   XS     = ABS(XL)+ABS(XU)
   X0     = XU-XL
   EPT    = EPS*XS+EP1

   !     EPT IS CONVERGENCE TOLERANCE FOR KTH T-EIGENVALUE

   IF (X0.LE.EPT) GO TO 230
   !
   !     T-EIGENVALUE HAS NOT YET CONVERGED
   !
   IC0 = IC0 + 1
   GO TO 330

   !     FORWARD STURM CALCULATION AT CURRENT T-EIGENVALUE APPROXIMATION.
   180 CONTINUE
 
   !     UPDATE T-EIGENVALUE INTERVAL (XL,XU)
 
   IF (NEV.LT.K) GO TO 190
   !
   !     NUMBER OF T-EIGENVALUES NEV = K
   XL = X1
   GO TO 170
   190 CONTINUE
   !     NUMBER OF T-EIGENVALUES NEV<K
   XU = X1
   NU = NEV

   !     UPDATE OF T-EIGENVALUE BOUNDS
 
   IF (NEV.EQ.0) GO TO 210
 
   DO 200 I = 1,NEV
     200 VB(I) = DMAX1(X1,VB(I))
     210 NEV1 = NEV+1

     DO 220 II = NEV1,K
       I = MT+II
       220 VB(I) = DMIN1(X1,VB(I))
       !
       GO TO 170
       !
       !     END (XL,XU) BISECTION LOOP FOR KTH T-EIGENVALUE ON (LB,UB)
       !     TEST FOR T-MULTIPLICITY AND IF SIMPLE THEN TEST FOR SPURIOUSNESS
       !
       230 CONTINUE
       NDIS     = NDIS+1
       MD       = MD+1
       VS(NDIS) = X1
       !
       JSTURM = 1
       X1     = XL-EP0
       GO TO 370

       !     BACKWARD STURM CALCULATION
       240 KL = KEV
       JL     = JEV
       !
       JSTURM  = 2
       IC0     = IC0 + 2
       X1      = XU+EP0

       GO TO 370
       !     BACKWARD STURM CALCULATION
       250 JU = JEV
       KU = KEV
       !
       !     FOR T(1,MEV)
       !     NU - KU = NO. T-EIGENVALUES ON (XU, XU + EP0)
       !     KL - KU = NO. T-EIGENVALUES ON (XL - EP0, XU + EP0)
       !
       !     FOR T(2,MEV)
       !     JL -JU = NO. T-EIGENVALUES ON (XL - EP0, XU + EP0)
       !
       !     IS THIS A SIMPLE T-EIGENVALUE?
       !
       IF (KL-KU-1.EQ.0) GO TO 290
       !
       !     VS(NDIS) = KTH-T-EIGENVALUE OF (LB,UB) IS MULTIPLE AND HENCE GOOD
       IF (KU.EQ.NU) GO TO 280
       !     CONTINUE TO CHECK FOR T-MULTIPLICITY
       260 CONTINUE
       ISTURM = 5
       X1 = X1+EP0
       IC0 = IC0 + 1
       GO TO 330
       !     FORWARD STURM CALCULATION
       270 KNE = KU-NEV
       KU = NEV
       IF (KNE.NE.0) GO TO 260
       !     SPECIFY T-MULTIPLICITY = MP(NDIS)
       280 MPEV = KL-KU
       KNEW = KU
       GO TO 300
       !     END MULTIPLE CASE
       !
       !     T-EIGENVALUE IS SIMPLE   CHECK IF IT IS SPURIOUS
       290 CONTINUE
       MPEV = 1
       IF (JU.LT.JL) MPEV=0
       KNEW = K-1
       !
       !     X1 >= XU+EP0
       !     SPURIOUS TEST AND T-SIMPLE CASE COMPLETED
       !     START OF NEXT T-EIGENVALUE COMPUTATION
       !
       300 K = KNEW
       MP(NDIS) = MPEV
       IF (MPEV.GE.1) NG = NG + 1
       !
       IF (IWRITE.NE.0)THEN; IF(.NOT.silent_) WRITE(6,310) VS(NDIS),MPEV,MD,NG; ENDIF
       310 FORMAT(E25.16,3I6)
       !
       !     UPDATE STURM COUNT. IC0 = STURM COUNT FOR KTH T-EIGENVALUE
       ICT = ICT + IC0
       !
       !     EXIT TEST FOR K DO LOOP
       !
       IF (K.LE.0) GO TO 410

       !     UPDATE LOWER BOUNDS
       DO 320 I=1,KNEW
  320 VB(I) = DMAX1(X1,VB(I))

  GO TO 150
  !     END OF BISECTION LOOP FOR KTH T-EIGENVALUE


  !     FORWARD STURM CALCULATION
  330 NEV = -NA
  YU      = ONE

  DO 360 I = 1,MEV

    IF (YU.NE.ZERO) GO TO 340
    YV     = BETA(I)/EPS
    GO TO 350
    340 YV = BETA2(I)/YU
    350 YU = X1 - ALPHA(I) - YV
    IF (YU.GE.ZERO) GO TO 360
    NEV    = NEV + 1
    360 CONTINUE

    !     NEV = NUMBER OF T-EIGENVALUES ON (X1,UB)

    GO TO (50,60,160,180,270), ISTURM

    !     BACKWARD STURM CALCULATION FOR T(1,MEV) AND T(2,MEV)
    370 KEV = -NA
    YU  = ONE
 
    DO 400 II = 1,MEV
      I = MP1-II
      IF (YU.NE.ZERO) GO TO 380
      YV = BETA(I+1)/EPS
      GO TO 390
      380 YV = BETA2(I+1)/YU
      390 YU = X1-ALPHA(I)-YV
      JEV = 0
      IF (YU.GE.ZERO) GO TO 400
      KEV = KEV+1
      JEV = 1
      400 CONTINUE
      JEV = KEV-JEV
      !
      GO TO (240,250), JSTURM
      !
      !     KEV = -NA + (NUMBER OF T(1,MEV) EIGENVALUES) > X1
      !     JEV = -NA + (NUMBER OF T(2,MEV) EIGENVALUES) > X1
      !     SET PARAMETERS FOR NEXT INTERVAL
      410 CONTINUE
      IC = ICT+IC
      MXSTUR = MXSTUR-ICT
      !
      IF(.NOT.silent_) WRITE(6,420) JIND,NG,MD
      420 FORMAT(/' T-EIGENVALUE CALCULATION ON INTERVAL',I6,'  IS COMPLETE'/&
      3X,'NO. GOOD',3X,'NO. DISTINCT'/I10,I13)
  
      430 CONTINUE

      !     END LOOP ON THE SUBINTERVALS (LB(J),UB(J)), J=1,NINT
      !     ISKIP OUTPUT

      IF (ISKIP.GT.0)THEN; IF(.NOT.silent_) WRITE(6,440)ISKIP; ENDIF
      440 FORMAT(' BISEC DEFAULTED ON',I3,3X,'INTERVALS'/&
      ' DEFAULTS OCCUR IF AN INTERVAL HAS NO T-EIGENVALUES'/&
      ' OR THE STURM ESTIMATE EXCEEDS THE USER-SPECIFIED LIMIT'/)
      !
      IF (ISKIP.GT.0)THEN; IF(.NOT.silent_) WRITE(6,450)(IDEF(I), I=1,ISKIP); ENDIF
      450 FORMAT(' BISEC DEFAULTED ON INTERVALS'/(10I8))
      !
      !     RESET BETA AT I = MP1
      BETA(MP1) = BETAM

      RETURN
      END subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

      !-----INVERSE ITERATION ON T(1,MEV)-------------------------------------

      SUBROUTINE INVERR(ALPHA,BETA,V1,V2,VS,EPS,G,MP,MEV,MMB,NDIS,NISO,N,IKL,IT,IWRITE,SILENT,UNIT4)

        LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
        LOGICAL                       :: silent_
        REAL(8)                       :: ALPHA(:),BETA(:),V1(:),V2(:),VS(:)
        REAL(8)                       :: X1,U,Z,EST,TEMP,T0,T1,RATIO,SUM,XU,NORM_,TSUM
        REAL(8)                       :: BETAM,EPS,EPS3,EPS4,ZERO_,ONE_
        REAL(8)                       :: G(:),GSUM 
        INTEGER                       :: MP(:),IKL(:),ILL(SIZE(IKL))
        INTEGER                       :: NINT,MEV,NDIS,IC,MXSTUR,ISKIP,NA,MD,JIND,MP1,NG,ICT,ISTURM,NEV,IC0
        INTEGER                       :: KEV,JEV,MPEV,KNE,ITER,KNEW,I,MM1,II,NM1,J,MEVI,ISO,ISPUR,IGOOD,MEVPNI 
        INTEGER                       :: IWRITE,MMB,NISO,N,IT,IB,UNIT4
        REAL                          :: GAP

        silent_ = F
        IF(PRESENT(SILENT)) silent_ = SILENT

        !-----------------------------------------------------------------------
        !     COMPUTES ERROR ESTIMATES FOR COMPUTED ISOLATED GOOD T-EIGENVALUES
        !     IN VS AND WRITES THESE T-EIGENVALUES AND ESTIMATES TO FILE 4.
        !     BY DEFINITION A GOOD T-EIGENVALUE IS ISOLATED IF ITS
        !     CLOSEST T-NEIGHBOR IS ALSO GOOD, OR ITS CLOSEST NEIGHBOR IS
        !     SPURIOUS, BUT THAT NEIGHBOR IS FAR ENOUGH AWAY.  SO
        !     IN PARTICULAR, WE COMPUTE ESTIMATES FOR GOOD T-EIGENVALUES
        !     THAT ARE IN CLUSTERS OF GOOD T-EIGENVALUES.
        !
        !     USES INVERSE ITERATION ON T(1,MEV) SOLVING THE EQUATION
        !     (T - X1*I)V2 = RIGHT-HAND SIDE (RANDOMLY-GENERATED)
        !     FOR EACH SUCH GOOD T-EIGENVALUE X1.
        !
        !     PROGRAM REFACTORS T-X1*I ON EACH ITERATION OF INVERSE ITERATION.
        !     TYPICALLY ONLY ONE ITERATION IS NEEDED PER EIGENVALUE X1.
        !
        !     POSSIBLE STORAGE COMPRESSION
        !     G STORAGE COULD BE ELIMINATED BY REGENERATING THE RANDOM
        !     RIGHT-HAND SIDE ON EACH ITERATION AND PRINTING OUT THE
        !     ERROR ESTIMATES AS THEY ARE GENERATED.
        !
        !     ON ENTRY AND EXIT
        !     MEV = ORDER OF T
        !     ALPHA, BETA CONTAIN THE NONZERO ENTRIES OF THE T-MATRIX
        !     VS = COMPUTED DISTINCT EIGENVALUES OF T(1,MEV)
        !     MP = T-MULTIPLICITY OF EACH T-EIGENVALUE IN VS. MP(I) = -1 MEANS
        !          VS(I) IS A GOOD T-EIGENVALUE BUT THAT IT IS SITTING CLOSE TO
        !          A SPURIOUS T-EIGENVALUE.  MP(I) = 0 MEANS VS(I) IS SPURIOUS.
        !          ESTIMATES ARE COMPUTED ONLY FOR THOSE T-EIGENVALUES
        !          WITH MP(I) = 1. FLAGGING WAS DONE IN SUBROUTINE ISOEV
        !          PRIOR TO ENTERING INVERR.
        !     NISO = NUMBER OF ISOLATED GOOD T-EIGENVALUES CONTAINED IN VS
        !     NDIS =  NUMBER OF DISTINCT T-EIGENVALUES IN VS
        !     IKL = SEED FOR RANDOM NUMBER GENERATOR
        !     EPS = 2. * MACHINE EPSILON
        !
        !     IN PROGRAM:
        !     ITER = MAXIMUM NUMBER OF INVERSE ITERATION STEPS ALLOWED FOR EACH
        !            X1.  ITER = IT ON ENTRY.
        !     G = ARRAY OF DIMENSION AT LEAST MEV + NISO.  USED TO STORE
        !         RANDOMLY-GENERATED RIGHT-HAND SIDE.  THIS IS NOT
        !         REGENERATED FOR EACH X1. G IS ALSO USED TO STORE ERROR
        !         ESTIMATES AS THEY ARE COMPUTED FOR LATER PRINTOUT.
        !     V1,V2 = WORK SPACES USED IN THE FACTORIZATION OF T(1,MEV).
        !     AT THE END OF THE INVERSE ITERATION COMPUTATION FOR X1, V2
        !     CONTAINS THE UNIT EIGENVECTOR OF T(1,MEV) CORRESPONDING TO X1.
        !     V1 AND V2 MUST BE OF DIMENSION AT LEAST MEV.
        !
        !     ON EXIT
        !     G(J) = MINIMUM GAP IN T(1,MEV) FOR EACH VS(J), J=1,NDIS
        !     G(MEV+I) = BETAM*|V2(MEV)| = ERROR ESTIMATE FOR ISOLATED GOOD
        !              T-EIGENVALUES, WHERE I = 1,NISO  AND  BETAM = BETA(MEV+1)
        !              V2(MEV) IS LAST COMPONENT OF THE UNIT EIGENVECTOR OF
        !              T(1,MEV) CORRESPONDING TO ITH ISOLATED GOOD T-EIGENVALUE.
        !
        !     IF FOR SOME X1 IT.GT.ITER THEN THE ERROR ESTIMATE IN G IS MARKED
        !     WITH A - SIGN.
        !
        !     V2 = ISOLATED GOOD T-EIGENVALUES
        !     V1 = MINIMAL T-GAPS FOR THE T-EIGENVALUES IN V2.
        !     THESE ARE CONSTRUCTED FOR WRITE-OUT PURPOSES ONLY AND NOT
        !     NEEDED ELSEWHERE IN THE PROGRAM.
        !-----------------------------------------------------------------------

        !     LABEL OUTPUT FILE  4
        IF (MMB.EQ.1.and.rank==0) WRITE(UNIT4,10)
        10 FORMAT(' INVERSE ITERATION ERROR ESTIMATES'/)
        !
        !     FILE 6 (TERMINAL) OUTPUT OF ERROR ESTIMATES
        IF (IWRITE.NE.0.AND.NISO.NE.0)THEN; IF(.NOT.silent_) WRITE(6,20); ENDIF
        20 FORMAT(/' INVERSE ITERATION ERROR ESTIMATES'/'  JISO',' JDIST',8X,&
        'GOOD T-EIGENVALUE',4X,'BETAM*UM',5X,'TMINGAP')
        !
        !     INITIALIZATION AND PARAMETER SPECIFICATION
        ZERO_ = 0.0D0
        ONE_  = 1.0D0
        NG    = 0
        NISO  = 0
        ITER  = IT
        MP1   = MEV+1
        MM1   = MEV-1
        BETAM = BETA(MP1)
        BETA(MP1) = ZERO_
        !
        !     CALCULATE SCALE AND TOLERANCES
        TSUM = DABS(ALPHA(1))
        DO 30 I = 2,MEV
   30 TSUM = TSUM + DABS(ALPHA(I)) + BETA(I)
   !
   EPS3 = EPS*TSUM
   EPS4 = DBLE(MEV)*EPS3

   !     GENERATE SCALED RANDOM RIGHT-HAND SIDE
   ILL = IKL

   !-----------------------------------------------------------------------
   CALL GENRAN(ILL,G,MEV)
   !-----------------------------------------------------------------------

   GSUM = ZERO_
   DO 40 I = 1,MEV
     40 GSUM = GSUM+ABS(G(I))
     GSUM = EPS4/GSUM
     !
     DO 50 I = 1,MEV
       50 G(I) = GSUM*G(I)
       !
       !     LOOP ON ISOLATED GOOD T-EIGENVALUES IN VS (MP(I) = 1) TO
       !     CALCULATE CORRESPONDING UNIT EIGENVECTOR OF T(1,MEV)
       !
       DO 180 JEV = 1,NDIS
         !
         IF (MP(JEV).EQ.0) GO TO 180
         NG = NG + 1
         IF (MP(JEV).NE.1) GO TO 180
         !
         IT = 1
         NISO = NISO + 1
         X1 = VS(JEV)
         !
         !     INITIALIZE RIGHT HAND SIDE FOR INVERSE ITERATION
         DO 60 I = 1,MEV
    60 V2(I) = G(I)
    !
    !     TRIANGULAR FACTORIZATION WITH NEAREST NEIGHBOR PIVOT
    !     STRATEGY. INTERCHANGES ARE LABELLED BY SETTING BETA < 0.
    !
    70 CONTINUE
    U = ALPHA(1)-X1
    Z = BETA(2)
    !
    DO 90 I = 2,MEV
      IF (BETA(I).GT.DABS(U)) GO TO 80
      !     NO INTERCHANGE
      V1(I-1) = Z/U
      V2(I-1) = V2(I-1)/U
      V2(I) = V2(I)-BETA(I)*V2(I-1)
      RATIO = BETA(I)/U
      U = ALPHA(I)-X1-Z*RATIO
      Z = BETA(I+1)
      GO TO 90
      80 CONTINUE
      !     INTERCHANGE CASE
      RATIO = U/BETA(I)
      BETA(I) = -BETA(I)
      V1(I-1) = ALPHA(I)-X1
      U = Z-RATIO*V1(I-1)
      Z = -RATIO*BETA(I+1)
      TEMP = V2(I-1)
      V2(I-1) = V2(I)
      V2(I) = TEMP-RATIO*V2(I)
      90 CONTINUE
      IF (U.EQ.ZERO_) U = EPS3
      !
      !     SMALLNESS TEST AND DEFAULT VALUE FOR LAST COMPONENT
      !     PIVOT(I-1) = |BETA(I)| FOR INTERCHANGE CASE
      !     (I-1,I+1) ELEMENT IN RIGHT FACTOR = BETA(I+1)
      !     END OF FACTORIZATION AND FORWARD SUBSTITUTION
      !
      !     BACK SUBSTITUTION
      V2(MEV) = V2(MEV)/U
      DO 110 II = 1,MM1
        I = MEV-II
        IF (BETA(I+1).LT.ZERO_) GO TO 100
        !     NO INTERCHANGE
        V2(I) = V2(I)-V1(I)*V2(I+1)
        GO TO 110
        !     INTERCHANGE CASE
        100 BETA(I+1) = -BETA(I+1)
        V2(I) = (V2(I)-V1(I)*V2(I+1)-BETA(I+2)*V2(I+2))/BETA(I+1)
        110 CONTINUE
        !
        !     TESTS FOR CONVERGENCE OF INVERSE ITERATION
        !     IF SUM |V2| COMPS. LE. 1 AND IT. LE. ITER DO ANOTHER INVIT STEP
        !
        NORM_ = DABS(V2(MEV))
        DO 120 II = 1,MM1
          I = MEV-II
          120 NORM_ = NORM_+DABS(V2(I))
          !
          IF (NORM_.GE.ONE_) GO TO 140
          IT = IT+1
          IF (IT.GT.ITER) GO TO 140
          XU = EPS4/NORM_
          !
          DO 130 I = 1,MEV
     130 V2(I) = V2(I)*XU
     !
     GO TO 70
     !     ANOTHER INVERSE ITERATION STEP
     !
     !     INVERSE ITERATION FINISHED
     !     NORMALIZE COMPUTED T-EIGENVECTOR : V2 = V2/||V2||
     140 CONTINUE

     SUM = DBLE(MPI_DOT_PRODUCT(V2(1:MEV),V2(1:MEV)))
     SUM = ONE_/DSQRT(SUM)
 
     DO 150 II = 1,MEV
       150 V2(II) = SUM*V2(II)
       !
       !     SAVE ERROR ESTIMATE FOR LATER OUTPUT
       EST = BETAM*DABS(V2(MEV))
       IF (IT.GT.ITER) EST = -EST
       MEVPNI = MEV + NISO
       G(MEVPNI) = EST
       IF (IWRITE.EQ.0) GO TO 180
       !
       !     FILE 6 (TERMINAL) OUTPUT OF ERROR ESTIMATES.
       IF (JEV.EQ.1) GAP = VS(2) - VS(1)
       IF (JEV.EQ.MEV) GAP = VS(MEV) - VS(MEV-1)
       IF (JEV.EQ.MEV.OR.JEV.EQ.1) GO TO 160
       TEMP = DMIN1(VS(JEV+1)-VS(JEV),VS(JEV)-VS(JEV-1))
       GAP = TEMP
       160 CONTINUE
       !
       IF(.NOT.silent_) WRITE(6,170) NISO,JEV,X1,EST,GAP
       170 FORMAT(2I6,E25.16,2E12.3)
       !
       180 CONTINUE
       !
       !     END ERROR ESTIMATE LOOP ON ISOLATED GOOD T-EIGENVALUES.
       !     GENERATE DISTINCT MINGAPS FOR T(1,MEV).  THIS IS USEFUL AS AN
       !     INDICATOR OF THE GOODNESS OF THE INVERSE ITERATION ESTIMATES.
       !     TRANSFER ISOLATED GOOD T-EIGENVALUES AND CORRESPONDING TMINGAPS
       !     TO V2 AND V1 FOR OUTPUT PURPOSES ONLY.
       !
       NM1 = NDIS - 1
       G(NDIS) = VS(NM1)-VS(NDIS)
       G(1) = VS(2)-VS(1)
       !
       DO 190 J = 2,NM1
         T0 = VS(J)-VS(J-1)
         T1 = VS(J+1)-VS(J)
         G(J) = T1
         IF (T0.LT.T1) G(J)=-T0
         190 CONTINUE
         ISO = 0
         DO 200 J = 1,NDIS
           IF (MP(J).NE.1) GO TO 200
           ISO = ISO+1
           V1(ISO) = G(J)
           V2(ISO) = VS(J)
           200 CONTINUE
           !
           IF(NISO.EQ.0) GO TO 250
           !
           !     ERROR ESTIMATES ARE WRITTEN TO FILE 4
           if(rank==0) WRITE(UNIT4,210)MEV,NDIS,NG,NISO,N,ITER,BETAM
           210 FORMAT(1X,'TSIZE',2X,'NDIS',1X,'NGOOD',2X,'NISO',1X,'ASIZE'/5I12/&
           4X,'MXINIT',5X,'BETAM'/I12,E10.3)
           if(rank==0) WRITE(UNIT4,211) 
           211 FORMAT(4X,'seedsize')
           if(rank==0) WRITE(UNIT4,212) SIZE(IKL)
           212 FORMAT(I12)
           if(rank==0) WRITE(UNIT4,213) 
           213 FORMAT(4X,'RHSEED')
           if(rank==0) WRITE(UNIT4,*) IKL
           if(rank==0) WRITE(UNIT4,214) 
           214 FORMAT(2X,'GOODEVNO',8X,'GOOD T-EIGENVALUE',6X,'BETAM*UM',7X,'TMINGAP')
           !
           ISPUR = 0
           I = 0
           DO 240 J = 1,NDIS
      IF(MP(J).NE.0) GO TO 220
      ISPUR = ISPUR + 1
      GO TO 240
      220 IF(MP(J).NE.1) GO TO 240
      I = I + 1
      MEVI = MEV + I
      IGOOD = J - ISPUR
      if(rank==0) WRITE(UNIT4,230) IGOOD,V2(I),G(MEVI),V1(I)
      230 FORMAT(I10,E25.16,2E14.3)
      240 CONTINUE
      GO TO 270
      !
      250 if(rank==0) WRITE(UNIT4,260)
      260 FORMAT(/' THERE ARE NO ISOLATED T-EIGENVALUES SO NO ERROR ESTIMATE WERE COMPUTED')
      !
      !     RESTORE BETA(MEV+1) = BETAM
      270 BETA(MP1) = BETAM
      RETURN
      END subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

      SUBROUTINE TNORM(ALPHA,BETA,BMIN,TMAX,MEV,IB,SILENT)

        LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
        LOGICAL                       :: silent_
        REAL(8)                       :: ALPHA(:),BETA(:)
        REAL(8)                       :: TMAX,BMIN,BMAX,BSIZE,BTOL
        INTEGER                       :: I,J,MEV,ILOOP,IB 

        silent_ = F
        IF(PRESENT(SILENT)) silent_ = SILENT

        !-----------------------------------------------------------------------
        !     COMPUTE SCALING FACTOR USED IN THE T-MULTIPLICITY, SPURIOUS AND
        !     PRTESTS.   CHECK RELATIVE SIZE OF THE BETA(K), K=1,MEV
        !     AS A TEST ON THE LOCAL ORTHOGONALITY OF THE LANCZOS VECTORS.
        !
        !             TMAX = MAX (|ALPHA(I)|, BETA(I),  I=1,MEV)
        !             BMIN = MIN (BETA(I) I=2,MEV)
        !             BSIZE = BMIN/TMAX
        !             |IB| = INDEX OF MINIMAL(BETA)
        !             IB < 0 IF BMIN/TMAX < BTOL
        !-----------------------------------------------------------------------
        !     SPECIFY PARAMETERS

        IB   = 2
        BTOL = BMIN
        BMIN = BETA(2)
        BMAX = BETA(2)
        TMAX = DABS(ALPHA(1))
  
         write(*,*) BETA
  
         DO 20 I = 2,MEV
           IF (BETA(I).GE.BMIN) GO TO 10
           IB = I
           write(*,*) I
           write(*,*) SIZE(BETA)
           write(*,*) BETA(I)
           write(*,*) BMIN
           if(I>size(BETA)) stop 'error Cumul lanczos, outside bounds'
           BMIN = BETA(I)
           write(*,*) TMAX,I
           writE(*,*) ALPHA(I)
           10 TMAX = DMAX1(TMAX,DABS(ALPHA(I)))
              BMAX = DMAX1(BETA(I),BMAX)
         20 CONTINUE

          TMAX = DMAX1(BMAX,TMAX)

          !     TEST OF LOCAL ORTHOGONALITY USING SCALED BETAS
          BSIZE = BMIN/TMAX
          IF (BSIZE.GE.BTOL) GO TO 40

          !
          !     DEFAULT.  BSIZE IS SMALLER THAN TOLERANCE BTOL SPECIFIED IN MAIN
          !     PROGRAM.  PROGRAM TERMINATES FOR USER TO DECIDE WHAT TO DO
          !     BECAUSE LOCAL ORTHOGONALITY OF THE LANCZOS VECTORS COULD BE
          !     LOST.

          IB = -IB
          IF(.NOT.silent_) WRITE(6,30) MEV
          30 FORMAT(/' BETA TEST INDICATES POSSIBLE LOSS OF LOCAL ORTHOGONALITY OVER 1ST',I6,' LANCZOS VECTORS'/)
          !
          40 CONTINUE
          !
          IF(.NOT.silent_) WRITE(6,50) IB
          50 FORMAT(/' MINIMUM BETA RATIO OCCURS AT',I6,' TH BETA'/)
          !
          IF(.NOT.silent_) WRITE(6,60) MEV,BMIN,TMAX,BSIZE
          60 FORMAT(/1X,'TSIZE',6X,'MIN BETA',5X,'TKMAX',6X,'MIN RATIO'/&
          I6,E14.3,E10.3,E15.3/)
          !
          RETURN
          END subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

          SUBROUTINE LUMP(V1,RELTOL,MULTOL,SCALE2,LINDEX,LOOP,SILENT)
            LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
            LOGICAL                       :: silent_
            INTEGER                       :: NLOOP,J,ICOUNT,I,JI,JF,INDSUM,ISPUR,KK,k,NM1,IDIF
            INTEGER                       :: LOOP,NDIS,NISO,NG
            REAL(8)                       :: V1(:),SUM,RELTOL,MULTOL,THOLD,ZERO,SCALE2
            INTEGER                       :: LINDEX(:)

            silent_ = F
            IF(PRESENT(SILENT)) silent_ = SILENT
            !-----------------------------------------------------------------------
            !     LINDEX(J) = T-MULTIPLICITY OF JTH DISTINCT T-EIGENVALUE
            !     LOOP = NUMBER OF DISTINCT T-EIGENVALUES
            !     LUMP 'COMBINES' COMPUTED 'GOOD' T-EIGENVALUES THAT ARE
            !     'TOO CLOSE'.
            !     VALUE OF RELTOL IS 1.D-10.
            !
            !     IF IN A SET OF T-EIGENVALUES TO BE COMBINED THERE IS AN EIGENVALUE
            !     WITH LINDEX=1, THEN THE VALUE OF THE COMBINED EIGENVALUES IS SET
            !     EQUAL TO THE VALUE OF THAT EIGENVALUE.  NOTE THAT IF A SPURIOUS
            !     T-EIGENVALUE IS TO BE 'COMBINED' WITH A GOOD T-EIGENVALUE, THEN
            !     THIS IS DONE ONLY BY INCREASING THE INDEX, LINDEX, FOR THAT
            !     T-EIGENVALUE.  NUMERICAL VALUES OF SPURIOUS EIGENVALUES ARE NEVER
            !     COMBINE WITH THOSE OF GOOD T-EIGENVALUES.
            !-----------------------------------------------------------------------

            ZERO   = 0.0D0
            NLOOP  = 0
            J      = 0
            ICOUNT = 1
            JI     = 1
            THOLD  = DMAX1(RELTOL*DABS(V1(1)),SCALE2*MULTOL)
            !
            10 J = J+1
            IF (J.EQ.LOOP) GO TO 20
            SUM = DABS(V1(J)-V1(J+1))
            IF (SUM.LT.THOLD) GO TO 60
            20 JF = JI + ICOUNT - 1
            INDSUM = 0
            ISPUR = 0
            !
            DO 30 KK = JI,JF
       IF (LINDEX(KK).NE.0) GO TO 30
       ISPUR = ISPUR + 1
       INDSUM = INDSUM + 1
       30 INDSUM = INDSUM + LINDEX(KK)
       !
       !     IF (JF-JI.GE.1) WRITE(6,40) (V1(KKK), KKK=JI,JF)
       40 FORMAT(/' LUMP LUMPS THE T-EIGENVALUES'/(4E20.13))
       !
       !     COMPUTE THE 'COMBINED' T-EIGENVALUE AND THE RESULTING
       !     T-MULTIPLICITY
       K = JI - 1
       50 K = K+1
       IF (K.GT.JF) GO TO 70
       IF (LINDEX(K) .NE.1) GO TO 50
       NLOOP = NLOOP + 1
       V1(NLOOP) = V1(K)
       GO TO 100
       60 ICOUNT = ICOUNT + 1
       GO TO 10
       !
       !     ALL INDICES WERE 0 OR >1
       70 NLOOP = NLOOP + 1
       IDIF = INDSUM - ISPUR
       IF (IDIF.EQ.0) GO TO 90
       !
       SUM = ZERO
       DO 80 KK = JI,JF
         80 SUM = SUM + V1(KK) * DBLE(LINDEX(KK))
         !
         V1(NLOOP) = SUM/DBLE(IDIF)
         GO TO 100
         90 V1(NLOOP) = V1(JI)
         100 LINDEX(NLOOP) = INDSUM
         IDIF = INDSUM - ISPUR
         IF (IDIF.EQ.0.AND.ISPUR.EQ.1) LINDEX(NLOOP) = 0
         IF (J.EQ.LOOP) GO TO 110
         ICOUNT = 1
         JI= J+1
         THOLD = DMAX1(RELTOL*DABS(V1(JI)),SCALE2*MULTOL)
         IF (JI.LT.LOOP) GO TO 10
         NLOOP = NLOOP + 1
         V1(NLOOP)= V1(JI)
         LINDEX(NLOOP) = LINDEX(JI)
         110 CONTINUE
         !
         !     ON RETURN V1 CONTAINS THE DISTINCT T-EIGENVALUES
         !     LINDEX CONTAINS THE CORRESPONDING T-MULTIPLICITIES
         !
         LOOP = NLOOP
         RETURN
         END subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

         SUBROUTINE ISOEV(VS,GAPTOL,MULTOL,SCALE1,G,MP,NDIS,NG,NISO,SILENT)
           LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
           LOGICAL                       :: silent_
           REAL(8)                       :: VS(:),T0,T1,MULTOL,GAPTOL,SCALE1,TEMP
           REAL(8)                       :: G(:)
           REAL                          :: GAP
           INTEGER                       :: MP(:),NM1,I,J,NSIGMA,LOOP,NDIS,NISO,NG

           silent_ = F
           IF(PRESENT(SILENT)) silent_ = SILENT

           !-----------------------------------------------------------------------
           !     GENERATE DISTINCT TMINGAPS AND USE THEM TO LABEL THE ISOLATED
           !     GOOD T-EIGENVALUES THAT ARE VERY CLOSE TO SPURIOUS ONES.
           !     ERROR ESTIMATES WILL NOT BE COMPUTED FOR THESE T-EIGENVALUES.
           !
           !     ON ENTRY AND EXIT
           !     VS CONTAINS THE COMPUTED DISTINCT EIGENVALUES OF T(1,MEV)
           !     MP CONTAINS THE CORRESPONDING T-MULTIPLICITIES
           !     NDIS = NUMBER OF DISTINCT EIGENVALUES
           !     GAPTOL = RELATIVE GAP TOLERANCE SET IN MAIN
           !
           !     ON EXIT
           !     G CONTAINS THE TMINGAPS.
           !     G(I) < 0 MEANS MINGAP IS DUE TO LEFT GAP
           !     MP(I) IS NOT CHANGED EXCEPT THAT  MP(I)=-1, IF MP(I)=1,
           !     TMINGAP WAS TOO SMALL AND DUE TO A SPURIOUS T-EIGENVALUE.
           !
           !     IF MP(I)=-1 THAT SIMPLE GOOD T-EIGENVALUE WILL BE SKIPPED
           !     IN THE SUBSEQUENT ERROR ESTIMATE COMPUTATIONS IN INVERR
           !     THAT IS, WE COMPUTE ERROR ESTIMATES ONLY FOR THOSE GOOD
           !     T-EIGENVALUES WITH MP(I)=1.
           !-----------------------------------------------------------------------

           !     CALCULATE MINGAPS FOR DISTINCT T(1,MEV) EIGENVALUES.

           NM1 = NDIS - 1
           G(NDIS) = VS(NM1)-VS(NDIS)
           G(1) = VS(2)-VS(1)
           !
           DO 10 J = 2,NM1
             T0 = VS(J)-VS(J-1)
             T1 = VS(J+1)-VS(J)
             G(J) = T1
             IF (T0.LT.T1) G(J) = -T0
             10 CONTINUE
             !
             !     SET MP(I)=-1 FOR SIMPLE GOOD T-EIGENVALUES WHOSE MINGAPS  ARE
             !     'TOO SMALL' AND DUE TO SPURIOUS T-EIGENVALUES.
             !
             NISO = 0
             NG = 0
             DO 20 J = 1,NDIS
        IF (MP(J).EQ.0) GO TO 20
        NG = NG+1
        IF (MP(J).NE.1) GO TO 20
        !     VS(J) IS NEXT SIMPLE GOOD T-EIGENVALUE
        NISO = NISO + 1
        I = J+1
        IF(G(J).LT.0.0) I = J-1
        IF(MP(I).NE.0) GO TO 20
        GAP  = ABS(G(J))
        T0   = DMAX1(SCALE1*MULTOL,GAPTOL*DABS(VS(J)))
        TEMP = T0
        IF (GAP.GT.TEMP) GO TO 20
        MP(J) = -MP(J)
        NISO  = NISO-1
        20 CONTINUE
        !
        RETURN
        END subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

        SUBROUTINE PRTEST(ALPHA,BETA,TEIG,TKMAX,EPSM,RELTOL,SCALE3,SCALE4,TMULT,NDIST,MEV,IPROJ,SILENT)
          LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
          LOGICAL :: silent_
          REAL(8) :: ALPHA(:), BETA(:),TEIG(:),SIGMA(10)
          REAL(8) :: EPSM,RELTOL,PRTOL,TKMAX,LRATIO,URATIO
          REAL(8) :: EPS,EPS1,BETAM,LBD,UBD,SIG,YU,YV,LRATS,URATS
          REAL(8) :: ZERO,ONE,TEN,BISTOL,SCALE3,SCALE4,AEV,TEMP
          INTEGER :: TMULT(:),ISIGMA(10),I,J,NSIGMA,IFIN,MF,ML,MEVL,IL,MEV1L,MEVU,MEV1U,NEV1,K
          INTEGER ::  MM1,NDIST1,ICOUNT,MEVUS,MEVLS,IPROJ, NDIST,MEV

          silent_ = F

          IF(PRESENT(SILENT)) silent_ = SILENT

          !-----------------------------------------------------------------------
          !     AFTER CONVERGENCE HAS BEEN ESTABLISHED, SUBROUTINE PRTEST
          !     TESTS COMPUTED EIGENVALUES OF T(1,MEV) THAT HAVE BEEN LABELLED
          !     SPURIOUS TO DETERMINE IF ANY EIGENVALUES OF A HAVE BEEN
          !     MISSED BY LANCZOS PROCEDURE.  AN EIGENVALUE WITH A VERY SMALL
          !     PROJECTION ON THE STARTING VECTOR (< SINGLE PRECISION)
          !     CAN BE MISSED BECAUSE IT IS ALSO AN EIGENVALUE OF T(2,MEV) TO
          !     WITHIN THE SQUARE OF THIS ORIGINAL PROJECTION.
          !     OUR EXPERIENCE IS THAT SUCH SMALL PROJECTIONS OCCUR ONLY
          !     VERY INFREQUENTLY.
          !
          !     THIS SUBROUTINE IS CALLED ONLY AFTER CONVERGENCE HAS BEEN
          !     ESTABLISHED. ONCE CONVERGENCE HAS BEEN OBSERVED ON THE
          !     OTHER EIGENVALUES THEN ONE CAN EXPECT TO ALSO HAVE CONVERGENCE
          !     ON ANY SUCH HIDDEN EIGENVALUES.(IF THERE ARE ANY).  THIS
          !     PROCEDURE CONSIDERS ONLY SPURIOUS T-EIGENVALUES AND ONLY THOSE
          !     SPURIOUS T-EIGENVALUES THAT ARE ISOLATED FROM GOOD T-EIGENVALUES.
          !     FOR EACH SUCH T-EIGENVALUE IT DOES 2 STURM SEQUENCES
          !     AND A FEW SCALAR MULTIPLICATIONS.  UPON RETURN TO MAIN
          !     PROGRAM ERROR ESTIMATES WILL BE COMPUTED FOR ANY EIGENVALUES
          !     THAT HAVE BEEN LABELLED AS 'HIDDEN'.  SUCH T-EIGENVALUES
          !     WILL BE RELABELLED AS 'GOOD' ONLY IF THESE ERROR ESTIMATES
          !     ARE SUFFICIENTLY SMALL.
          !-----------------------------------------------------------------------

          ZERO = 0.0D0
          ONE  = 1.0D0
          TEN  = 10.0D0
          PRTOL = 1.D-6
          TEMP = DBLE(MEV+1000)
          TEMP = DSQRT(TEMP)
          BISTOL = TKMAX*EPSM*TEMP
          NSIGMA = 4
          SIGMA(1) = TEN*TKMAX
          !
          DO 10 J = 2,NSIGMA
            10 SIGMA(J) = TEN*SIGMA(J-1)
            !
            IFIN = 0
            MF = 1
            ML = MEV
            BETAM = BETA(MF)
            BETA(MF) = ZERO
            IPROJ = 0
            J = 1
            !
            IF (TMULT(1).NE.0) GO TO 110
            !
            AEV = DABS(TEIG(1))
            TEMP = PRTOL*AEV
            EPS1 = DMAX1(TEMP,SCALE4*BISTOL)
            TEMP = RELTOL*AEV
            EPS  = DMAX1(TEMP,SCALE3*BISTOL)
            IF (TEIG(2)-TEIG(1).LT.EPS1.AND.TMULT(2).NE.0) GO TO 110
            !
            20 LBD = TEIG(J) - EPS
            UBD = TEIG(J) + EPS
            MEVL = 0
            IL = 0
            YU = ONE
            !
            DO 50 I=MF,ML
              IF (YU.NE.ZERO) GO TO 30
              YV = BETA(I)/EPSM
              GO TO 40
              30 YV = BETA(I)*BETA(I)/YU
              40 YU = ALPHA(I)-LBD-YV
              IF (YU.GE.ZERO) GO TO 50
              !     MEVL INCREMENTED
              MEVL = MEVL + 1
              IL = I
              50 CONTINUE
              !
              LRATIO = YU
              MEV1L = MEVL
              IF (IL.EQ.ML) MEV1L=MEVL-1
              !
              !     MEVL = NUMBER OF EVS OF T(1,MEV) WHICH ARE < LBD
              !     MEV1L = NUMBER OF EVS OF T(1,MEV-1) WHICH ARE < LBD
              !     LRATIO = DET(T(1,MEV)-LBD)/DET(T(1,MEV-1)-LBD):
              !
              MEVU = 0
              IL = 0
              YU = ONE
              !
              DO 80 I=MF,ML
         IF (YU.NE.ZERO) GO TO 60
         YV = BETA(I)/EPSM
         GO TO 70
         60 YV = BETA(I)*BETA(I)/YU
         70 YU = ALPHA(I)-UBD-YV
         IF (YU.GE.ZERO) GO TO 80
         !     MEVU INCREMENTED
         MEVU = MEVU + 1
         IL = I
         80 CONTINUE
         !
         URATIO = YU
         MEV1U = MEVU
         IF (IL.EQ.ML) MEV1U=MEVU-1
         !
         !     MEVU = NUMBER OF EVS OF T(MEV) WHICH ARE < UBD
         !     MEV1U = NUMBER OF EVS OF T(MEV-1) WHICH ARE < UBD
         !     URATIO = DET(TM-UBD)/DET(T(M-1)-UBD): TM=T(MF,ML)
         !
         NEV1 = MEV1U-MEV1L
         !
         DO 90 K=1,NSIGMA
           SIG = SIGMA(K)
           LRATS = LRATIO-SIG
           URATS = URATIO-SIG
           !     NOTE THE INCREMENT IS ON NUMBER OF EVALUES OF T(M-1)
           MEVLS = MEV1L
           IF (LRATS.LT.0.) MEVLS=MEV1L+1
           MEVUS = MEV1U
           IF (URATS.LT.0.) MEVUS=MEV1U+1
           ISIGMA(K) = MEVUS - MEVLS
           90 CONTINUE
           !
           ICOUNT = 0
           DO 100 K=1,NSIGMA
             100 IF (ISIGMA(K).EQ.1) ICOUNT=ICOUNT + 1
             !
             IF (ICOUNT.LT.2.OR.NEV1.EQ.0) GO TO 110
             TMULT(J) = -10
             IPROJ=IPROJ+1
             !
             110 J=J+1
             !
             IF (J.GE.NDIST) GO TO 120
             IF (TMULT(J).NE.0) GO TO 110
             !
             AEV = DABS(TEIG(J))
             TEMP = PRTOL*AEV
             EPS1 = DMAX1(TEMP,SCALE4*BISTOL)
             TEMP = RELTOL*AEV
             EPS  = DMAX1(TEMP,SCALE3*BISTOL)
             !
             IF (TEIG(J)-TEIG(J-1).LT.EPS1.AND.TMULT(J-1).NE.0) GO TO 110
             IF (TEIG(J+1)-TEIG(J).LT.EPS1.AND.TMULT(J+1).NE.0) GO TO 110
             !
             GO TO 20
             !
             120 IF (IFIN.EQ.1) GO TO 130
             IF (TMULT(NDIST).NE.0) GO TO 130
             !
             AEV  = DABS(TEIG(NDIST))
             TEMP = PRTOL*AEV
             EPS1 = DMAX1(TEMP,SCALE4*BISTOL)
             TEMP = RELTOL*AEV
             EPS  = DMAX1(TEMP,SCALE3*BISTOL)
             !
             NDIST1=NDIST -1
             TEMP = TEIG(NDIST)-TEIG(NDIST1)
             IF (TEMP.LT.EPS1.AND.TMULT(NDIST1).NE.0) GO TO 130
             IFIN = 1
             !
             GO TO 20
             !
             130 BETA(MF) = BETAM
             !
             RETURN
             END subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
 
             SUBROUTINE STURMI(ALPHA,BETA,X1,TOLN,EPSM,MMAX,MK1,MK2,IC,IWRITE,SILENT)
               LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
               LOGICAL :: silent_
               REAL(8) :: ALPHA(:),BETA(:)
               REAL(8) :: EPSM,X1,TOLN,EVL,EVU,BETA2
               REAL(8) :: U1,U2,V1,V2,ZERO,ONE
               INTEGER :: I,IC,ICD,IC0,IC1,IC2,MK1,MK2,MMAX,IWRITE

                                   silent_ = F
               IF(PRESENT(SILENT)) silent_ = SILENT

               !-----------------------------------------------------------------------
               !
               !     FOR ANY EIGENVALUE OF A THAT HAS CONVERGED AS AN EIGENVALUE
               !     OF THE T-MATRICES THIS SUBROUTINE CALCULATES
               !     THE SMALLEST SIZE OF THE T-MATRIX, T(1,MK1) DEFINED
               !     BY THE ALPHA AND BETA ARRAYS SUCH THAT MK1.LE.MMAX
               !     AND THE INTERVAL (X1-TOLN,X1+TOLN) CONTAINS AT LEAST ONE
               !     EIGENVALUE OF T(1,MK1). IT ALSO CALCULATES MK2 <= MMAX
               !     AS THE SMALLEST SIZE T-MATRIX (IF ANY) SUCH THAT THIS INTERVAL
               !     CONTAINS AT LEAST TWO EIGENVALUES OF T(1,MK2).
               !     IF NO T-MATRIX OF ORDER < MMAX SATISFIES THIS REQUIREMENT
               !     THEN MK2 IS SET EQUAL TO MMAX.  THE EIGENVECTOR PROGRAM
               !     USES THESE VALUES TO DETERMINE AN APPROPRIATE 1ST GUESS AT
               !     AN APPROPRIATE SIZE T-MATRIX FOR THE EIGENVALUE X1.
               !
               !     ON EXIT IC = NUMBER OF EIGENVALUES OF T(1,MK2) IN THIS INTERVAL
               !
               !     STURMI REGENERATES THE QUANTITIES BETA(I)**2 EACH TIME IT IS
               !     CALLED, OBVIOUSLY FOR THE PRICE OF ANOTHER VECTOR OF LENGTH
               !     MMAX THIS GENERATION COULD BE DONE ONCE IN THE MAIN
               !     PROGRAM BEFORE THE LOOP ON THE CALLS TO SUBROUTINE STURMI.
               !
               !     IF ANY OF THE EIGENVALUES BEING CONSIDERED WERE MULTIPLE
               !     AS EIGENVALUES OF THE USER-SPECIFIED MATRIX, THEN
               !     THIS SUBROUTINE COULD BE MODIFIED TO COMPUTE ADDITIONAL
               !     SIZES MKJ, J = 3, ...  WHICH COULD THEN BE USED IN THE
               !     MAIN LANCZOS EIGENVECTOR PROGRAM TO COMPUTE ADDITIONAL
               !     EIGENVECTORS CORRESPONDING TO THESE MULTIPLE EIGENVALUES.
               !     THE MAIN PROGRAM PROVIDED DOES NOT INCLUDE THIS OPTION.
               !
               !-----------------------------------------------------------------------

               !     INITIALIZATION OF PARAMETERS

               MK1     = 0
               MK2     = 0
               ZERO    = 0.0D0
               ONE     = 1.0D0
               BETA(1) = ZERO
               EVL     = X1-TOLN
               EVU     = X1+TOLN
               U1      = ONE
               U2      = ONE

               IC0 = 0
               IC1 = 0
               IC2 = 0
  
    !     MAIN LOOP FOR CALCULATING THE SIZES MK1,MK2
          DO 60 I = 1,MMAX
          BETA2 = BETA(I)*BETA(I)
          IF (U1.NE.ZERO) GO TO 10
          V1 = BETA(I)/EPSM
          GO TO 20
          10 V1 = BETA2/U1
          20 U1 = EVL - ALPHA(I) - V1
          IF (U1.LT.ZERO) IC1 = IC1+1
          IF (U2.NE.ZERO) GO TO 30
          V2 = BETA(I)/EPSM
          GO TO 40
          30 V2 = BETA2/U2
          40 U2 = EVU - ALPHA(I) - V2
          IF (U2.LT.ZERO) IC2 = IC2+1

    !     TEST FOR CHANGE IN NUMBER OF T-EIGENVALUES ON (EVL,EVU)
          ICD = IC1-IC2
          IC = ICD-IC0
          IF (IC.GE.1) GO TO 50
          GO TO 60
          50 CONTINUE
          IF (IC0.EQ.0) MK1 = I
          IC0 = IC0+1
          IF (IC0.GT.1) GO TO 70
          60 CONTINUE
 
          I = I-1
          IF (IC0.EQ.0) MK1 = MMAX
          70 MK2 = I
          IC = ICD
 
          IF (IWRITE.EQ.1)THEN; IF(.NOT.silent_) WRITE(6,80) X1,MK1,MK2,IC; ENDIF
          80 FORMAT(' EVAL =',E20.12,' MK1 =',I6,' MK2 =',I6,' IC =',I3/)
    
          RETURN
          END subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************




          !-----START OF INVERM---------------------------------------------------

          SUBROUTINE INVERM(ALPHA,BETA,V1,V2,X1,ERROR,ERRORV,EPS,G,MEV,IT,IWRITE,SILENT)
            LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
            LOGICAL                       :: silent_
            REAL(8)                       :: ALPHA(:),BETA(:),V1(:),V2(:)
            REAL(8)                       :: X1,U,Z,TEMP,RATIO,SUM,XU,NORM,TSUM,BETAM
            REAL(8)                       :: EPS,EPS3,EPS4,ERROR,ERRORV,ZERO,ONE
            REAL(8)                       :: G(:),GSUM 
            INTEGER                       :: MP1,NDIST1,MEVLS,MEVUS,ITER,I,J,K,II,NEV,NA,JSTURM,JM,MEV,MM1,IWRITE,IT 

            silent_ = F

            IF(PRESENT(SILENT)) silent_ = SILENT

            !-----------------------------------------------------------------------
            !
            !     COMPUTES T-EIGENVECTORS FOR ISOLATED GOOD T-EIGENVALUES X1
            !     USING INVERSE ITERATION ON T(1,MEV(X1)) SOLVING EQUATION
            !     (T - X1*I)V2 = RIGHT-HAND SIDE (RANDOMLY-GENERATED) .
            !     PROGRAM REFACTORS T- X1*I ON EACH ITERATION OF INVERSE ITERATION.
            !     TYPICALLY ONLY ONE ITERATION IS NEEDED PER T-EIGENVALUE X1.
            !
            !     IF IWRITE = 1 THEN THERE ARE EXTENDED WRITES TO FILE 6 (TERMINAL)
            !
            !     ON ENTRY G CONTAINS A REAL*4 RANDOM VECTOR WHICH WAS GENERATED
            !     IN MAIN PROGRAM.
            !
            !     ON ENTRY AND EXIT
            !     MEV = ORDER OF T
            !     ALPHA, BETA CONTAIN THE DIAGONAL AND OFFDIAGONAL ENTRIES OF T.
            !     EPS = 2. * MACHINE EPSILON
            !
            !     IN PROGRAM:
            !     ITER = MAXIMUM NUMBER STEPS ALLOWED FOR INVERSE ITERATION
            !     ITER = IT ON ENTRY.
            !     V1,V2 = WORK SPACES USED IN THE FACTORIZATION OF T(1,MEV).
            !     V1 AND V2 MUST BE OF DIMENSION AT LEAST MEV.
            !
            !     ON EXIT
            !     V2 = THE UNIT EIGENVECTOR OF T(1,MEV) CORRESPONDING TO X1.
            !     ERROR = |V2(MEV)| = ERROR ESTIMATE FOR CORRESPONDING
            !             RITZ VECTOR FOR X1.
            !
            !     ERRORV = || T*V2 - X1*V2 || = ERROR ESTIMATE ON T-EIGENVECTOR.
            !     IF IT.GT.ITER THEN ERRORV = -ERRORV
            !     IT = NUMBER OF ITERATIONS ACTUALLY REQUIRED
            !-----------------------------------------------------------------------


            !     INITIALIZATION AND PARAMETER SPECIFICATION

            ONE  = 1.0D0
            ZERO = 0.0D0
            ITER = IT
            MP1 = MEV+1
            MM1 = MEV-1
            BETAM = BETA(MP1)
            BETA(MP1) = ZERO
            !
            !     CALCULATE SCALE AND TOLERANCES
            TSUM = DABS(ALPHA(1))
            DO 10 I = 2,MEV
              10 TSUM = TSUM + DABS(ALPHA(I)) + BETA(I)
              !
              EPS3 = EPS*TSUM
              EPS4 = DBLE(MEV)*EPS3
              !
              !     GENERATE SCALED RANDOM RIGHT-HAND SIDE
              GSUM = ZERO
              DO 20 I = 1,MEV
                20 GSUM = GSUM+ABS(G(I))
                GSUM = EPS4/GSUM
                !
                !     INITIALIZE RIGHT HAND SIDE FOR INVERSE ITERATION
                DO 30 I = 1,MEV
           30 V2(I) = GSUM*G(I)
           IT = 1
           !
           !     CALCULATE UNIT EIGENVECTOR OF T(1,MEV) FOR ISOLATED GOOD
           !     T-EIGENVALUE X1.
           !
           !     TRIANGULAR FACTORIZATION WITH NEAREST NEIGHBOR PIVOT
           !     STRATEGY. INTERCHANGES ARE LABELLED BY SETTING BETA < 0.
           !
           40 CONTINUE
           U = ALPHA(1)-X1
           Z = BETA(2)
           !
           DO 60 I=2,MEV
             IF (BETA(I).GT.DABS(U)) GO TO 50
             !     NO PIVOT INTERCHANGE
             V1(I-1) = Z/U
             V2(I-1) = V2(I-1)/U
             V2(I) = V2(I)-BETA(I)*V2(I-1)
             RATIO = BETA(I)/U
             U = ALPHA(I)-X1-Z*RATIO
             Z = BETA(I+1)
             GO TO 60
             !     PIVOT INTERCHANGE
             50 CONTINUE
             RATIO = U/BETA(I)
             BETA(I) = -BETA(I)
             V1(I-1) = ALPHA(I)-X1
             U = Z-RATIO*V1(I-1)
             Z = -RATIO*BETA(I+1)
             TEMP = V2(I-1)
             V2(I-1) = V2(I)
             V2(I) = TEMP-RATIO*V2(I)
             60 CONTINUE
             !
             IF (U.EQ.ZERO) U=EPS3
             !
             !     SMALLNESS TEST AND DEFAULT VALUE FOR LAST COMPONENT
             !     PIVOT(I-1) = |BETA(I)| FOR INTERCHANGE CASE
             !     (I-1,I+1) ELEMENT IN RIGHT FACTOR = BETA(I+1)
             !     END OF FACTORIZATION AND FORWARD SUBSTITUTION
             !
             !     BACK SUBSTITUTION
             V2(MEV) = V2(MEV)/U
             DO 80 II = 1,MM1
               I = MEV-II
               IF (BETA(I+1).LT.ZERO) GO TO 70
               !     NO PIVOT INTERCHANGE
               V2(I) = V2(I)-V1(I)*V2(I+1)
               GO TO 80
               !     PIVOT INTERCHANGE
               70 BETA(I+1) = -BETA(I+1)
               V2(I) = (V2(I)-V1(I)*V2(I+1)-BETA(I+2)*V2(I+2))/BETA(I+1)
               80 CONTINUE
               !
               !
               !     TESTS FOR CONVERGENCE OF INVERSE ITERATION
               !     IF SUM |V2| COMPS. LE. 1 AND IT. LE. ITER DO ANOTHER INVIT STEP
               !
               NORM = DABS(V2(MEV))
               DO 90 II = 1,MM1
                 I = MEV-II
                 90 NORM = NORM+DABS(V2(I))
                 !
                 !     IS DESIRED GROWTH IN VECTOR ACHIEVED ?
                 !     IF NOT, DO ANOTHER INVERSE ITERATION STEP UNLESS NUMBER ALLOWED IS
                 !     EXCEEDED.
                 IF (NORM.GE.ONE) GO TO 110
                 !
                 IT=IT+1
                 IF (IT.GT.ITER) GO TO 110
                 !
                 XU = EPS4/NORM
                 DO 100 I=1,MEV
            100 V2(I) = V2(I)*XU
            !
            GO TO 40
            !
            !     NORMALIZE COMPUTED T-EIGENVECTOR : V2 = V2/||V2||
            !
            110 CONTINUE
 
            SUM = DBLE(MPI_DOT_PRODUCT(V2(1:MEV),V2(1:MEV)))
            SUM = ONE/DSQRT(SUM)

            DO 120 II = 1,MEV
              120 V2(II) = SUM*V2(II)
              !
              !     SAVE ERROR ESTIMATE FOR LATER OUTPUT
              ERROR = DABS(V2(MEV))
              !
              !     GENERATE ERRORV = ||T*V2 - X1*V2||.
              V1(MEV) = ALPHA(MEV)*V2(MEV)+BETA(MEV)*V2(MEV-1)-X1*V2(MEV)
              DO 130 J = 2,MM1
                JM = MP1 - J
                V1(JM) = ALPHA(JM)*V2(JM) + BETA(JM)*V2(JM-1) + BETA(JM+1)*V2(JM+1) - X1*V2(JM)
                130 CONTINUE

                V1(1) = ALPHA(1)*V2(1) + BETA(2)*V2(2) - X1*V2(1)
                ERRORV = DBLE(MPI_DOT_PRODUCT(V1(1:MEV),V1(1:MEV)))
                ERRORV = DSQRT(ERRORV)
                IF (IT.GT.ITER) ERRORV = -ERRORV
                IF (IWRITE.EQ.0) GO TO 150

                !     FILE 6 (TERMINAL) OUTPUT OF ERROR ESTIMATES.
                IF(.NOT.silent_) WRITE(6,140) MEV,X1,ERROR,ERRORV
                140 FORMAT(2X,'TSIZE',15X,'EIGENVALUE',11X,'U(M)',9X,'ERRORV'/I6,E25.16,2E15.5)
                !
                !     RESTORE BETA(MEV+1) = BETAM
                150 CONTINUE
                BETA(MP1) = BETAM
            RETURN
            END subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************




                !-----START OF LBISEC---------------------------------------------------

                SUBROUTINE LBISEC(ALPHA,BETA,EPSM,EVAL,EVALN,LB,UB,TTOL,M,NEVT,SILENT)
                  LOGICAL, OPTIONAL, INTENT(IN) :: SILENT
                  LOGICAL :: silent_
                  REAL(8) :: ALPHA(:),BETA(:),X0,X1,XL,XU,YU,YV,LB,UB
                  REAL(8) :: EPSM,EP1,EVAL,EVALN,EVD,EPT
                  REAL(8) :: ZERO,ONE,HALF,TTOL,TEMP
                  INTEGER :: NA,JSTURM,NEV,I,M,NEVT 

                  silent_ = F

                  IF(PRESENT(SILENT)) silent_ = SILENT
                  !-----------------------------------------------------------------------
                  !     SPECIFY PARAMETERS
                  ZERO = 0.0D0
                  HALF = 0.5D0
                  ONE  = 1.0D0
                  XL = LB
                  XU = UB
                  !
                  !     EP1 = DSQRT(1000+M)*TTOL     TTOL = EPSM*TKMAX
                  !     TKMAX = MAX(|ALPHA(K)|,BETA(K), K= 1,KMAX)
                  !
                  TEMP = DBLE(1000+M)
                  EP1 = DSQRT(TEMP)*TTOL
                  !
                  NA = 0
                  X1 = XU
                  JSTURM = 1
                  GO TO 60
                  !     FORWARD STURM CALCULATION
                  10 NA = NEV
                  X1 = XL
                  JSTURM = 2
                  GO TO 60
                  !     FORWARD STURM CALCULATION
                  20 NEVT = NEV
                  !
                  !     WRITE(6,30) M,EVAL,NEVT,EP1
                  30 FORMAT(/3X,'TSIZE',23X,'EV',9X/I8,E25.16/I6,' = NUMBER OF T(1,M) EIGENVALUES ON TEST INTERVAL'/&
                  E12.3,' =  CONVERGENCE TOLERANCE'/)
                  !
                  IF (NEVT.NE.1) GO TO 120
                  !
                  !     BISECTION LOOP
                  JSTURM = 3
                  40 X1 = HALF*(XL+XU)
                  X0 = XU-XL
                  EPT = EPSM*(DABS(XL) + DABS(XU)) + EP1
                  !     CONVERGENCE TEST
                  IF (X0.LE.EP1) GO TO 100
                  GO TO 60
                  !     FORWARD STURM CALCULATION
                  50 CONTINUE
                  IF(NEV.EQ.0) XU = X1
                  IF(NEV.EQ.1) XL = X1
                  GO TO 40
                  !     NEV = NUMBER OF T-EIGENVALUES OF T(1,M) ON (X1,XU)
                  !     THERE IS EXACTLY ONE T-EIGENVALUE OF T(1,M) ON (XL,XU)
                  !
                  !     FORWARD STURM CALCULATION
                  60 NEV = -NA
                  YU = ONE
                  DO 90 I = 1,M
             IF (YU.NE.ZERO) GO TO 70
             YV = BETA(I)/EPSM
             GO TO 80
             70 YV = BETA(I)*BETA(I)/YU
             80 YU = X1 - ALPHA(I) - YV
             IF (YU.GE.ZERO) GO TO 90
             NEV = NEV+1
             90 CONTINUE
             GO TO (10,20,50), JSTURM
             !
             100 CONTINUE
             !
             EVALN = X1
             EVD = DABS(EVALN-EVAL)
             !     WRITE(6,110) EVALN,EVAL,EVD
             110 FORMAT(/20X,'EVALN',21X,'EVAL',6X,'CHANGE'/2E25.16,E12.3/)
             !
             120 CONTINUE

        RETURN
        END SUbROUTINE


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************



   END MODULE 
