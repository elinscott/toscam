module splines

  use smooth_data
  use genvar
  use random
  use geometry
  use tools_fit
  use mpirout
  use derivative_noise
  use splines2, only: spline_overhauser_val
  use sorting, only: qsort_array, qsort_adj_array, group_data_rrr

  private
  public :: resampleit
  public :: resampleit___b

 !===============================================================!
   type spline
     real,dimension(:),allocatable   :: om,Frc,w,cc,tt,smooth
     real,dimension(:,:),allocatable :: der
     integer                         :: nn,kk,N,initialized=0
     integer                         :: order,nest,ierr
   end type
 !===============================================================!
! 
!  !---------------------------------------------------!
!   INTERFACE fourier_filter
!       MODULE PROCEDURE  fourier_filter_,fourier_filter__ 
!   END INTERFACE
 !---------------------------------------------------!
  INTERFACE spline_inter_extrapolate
      MODULE PROCEDURE spline_inter_extrapolate__, spline_inter_extrapolate_, &
          spline_inter_extrapolate_array__,spline_inter_extrapolate_array_
  END INTERFACE
 !---------------------------------------------------!
!   INTERFACE derivate_spline
!       MODULE PROCEDURE derivate_splinea,derivate_splineb,derivate_splinebr,derivate_splined
!   END INTERFACE
 !---------------------------------------------------!
  INTERFACE evaluate_spline
      MODULE PROCEDURE evaluate_splinea,evaluate_splineb,evaluate_splinebr
  END INTERFACE
 !---------------------------------------------------!
!   INTERFACE smoothit
!       MODULE PROCEDURE smoothit_,smoothit__,smoothit___,smoothitb_____
!   END INTERFACE
!  !---------------------------------------------------!
!   INTERFACE derivateit
!       MODULE PROCEDURE derivateit_,derivateit__,derivateit___,derivateit____
!   END INTERFACE
 !---------------------------------------------------!
  INTERFACE resampleit
      MODULE PROCEDURE resampleit_b, resampleit__b, resampleit___b, &
            resampleit_matrix, resampleit_matrix_, resampleit_xonly_r, resampleit_xonly_c
  END INTERFACE
 !---------------------------------------------------!
! 
!   INTERFACE
!     subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
!      implicit none
!       integer ::  n,k,nu,m,ier
!       real(4) ::  t(n),c(n),x(m),y(m),wrk(n)
!       integer ::  i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn
!     end subroutine
!   END INTERFACE

contains
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!   subroutine check_x_y(x,y)
!   implicit none
!   integer :: i 
!   real(8) :: x(:),y(:)
!   
!     do i=1,size(x)-1
!      if(x(i)>=x(i+1))then
!       write(*,*) 'error in smoothit, x not monotonic : ', x(i),x(i+1)
!       stop
!      endif
!     enddo
!   
!   return
!   end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! !--------------------------------------------------------------------!
! !  Extracted from "A practical guide to splines," 1st edition,       !
! !  Applied Mathematical Sciences 27, Carl de Boor, Springer, 1978.   !
! !--------------------------------------------------------------------!
! 
!       SUBROUTINE CUBSPL ( TAU, C, N, IBCBEG, IBCEND)
! !     N = NUMBER OF DATA POINTS. ASSUMED TO BE .GE. 2.
! !     (TAU(I), C(1,I), I=1,...,N) = ABSCISSAE AND ORDINATES OF THE
! !     DATA POINTS. TAU IS ASSUMED TO BE STRICTLY INCREASING.
! !     IBCBEG, IBCEND = BOUNDARY CONDITION INDICATORS, AND
! !     C(2,1), C(2,N) = BOUNDARY CONDITION INFORMATION. SPECIFICALLY,
! !        IBCBEG = 0  MEANS NO BOUNDARY CONDITION AT TAU(1) IS GIVEN.
! !           IN THIS CASE, THE NOT-A-KNOT CONDITION IS USED, I.E. THE
! !           JUMP IN THE THIRD DERIVATIVE ACROSS TAU(2) IS FORCED TO
! !           ZERO, THUS THE FIRST AND THE SECOND CUBI! POLYNOMIAL PIECES
! !           ARE MADE TO COINCIDE.)
! !        IBCBEG = 1  MEANS THAT THE SLOPE AT TAU(1) IS MADE TO EQUAL
! !           C(2,1), SUPPLIED BY INPUT.
! !        IBCBEG = 2  MEANS THAT YHE SECOND DERIVATIVE TAU(1) IS 
! !           MADE TO EQUAL C(2,1), SUPPLIED BY INPUT.
! !        IBCEND = 0, 1, OR 2 HAS ANALOGOUS MEANING CONCERNING THE
! !           BOUNDARY CONDITION AT TAU(N), WITH THE ADDITIONAL INFOR-
! !           MATION TAKEN FROM C(2,N).
! !     C(J,I), J=1,...,4; I=1,...,L (= N-1) = THE POLYNOMIAL COEFFICIENTS 
! !        OF THE CUBI! INTERPOLATING SPLINE WITH INTERIOR KNOTS (OR
! !        JOINTS) TAU(2),...,TAU(N-1). PRECISELY, IN THE 
! !        INTERVAL (TAU(I), TAU(I+1)), THE SPLINE F IS GIVEN BY
! !           F(X) = C(1,I)+H*(C(2,I)+H*(C(3,I)+H*C(4,I)/3.)/2.)
! !        WHERE H = X - TAU(I). THE FUNCTION PROGRAM *PPVALU* MAY BE 
! !        USED TO EVALUATE F OR ITS DERIVATIVES FROM TAU, C, L = N-1,
! !        AND K=4.
!       IMPLICIT NONE
!       INTEGER IBCBEG, IBCEND, N, I, J, L, M
!       DOUBLE PRECISION C(4,N), TAU(N), DIVDF1, DIVDF3, DTAU, G
! !     A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SPLOPES S(I) OF
! !     F AT TAU(I), I=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS 
! !     ELIMINATION, WITH S(I) ENDING UP IN C(2,I), ALL I.
! !     C(3,.) AND C(4,.) ARE USED INITIALLY FOR TEMPORARY STORAGE.
!       L = N - 1
! !     COMPUTE FIRST DIFFERENCES OF TAU SEQUENCE AND STORE IN C(3,.). ALSO,
! !     COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN C(4,.).
!       DO 10 M = 2, N
!          C(3,M) = TAU(M) - TAU(M-1)
!    10    C(4,M) = (C(1,M) - C(1,M-1))/C(3,M)
! !     CONSTRUCT FIRST EQUATION FROM THE BOUNDARY CONDITION, OF THE FORM
! !             C(4,1)*S(1) + C(3,1)*S(2) = C(2,1)
!       IF (IBCBEG-1)                     11,15,16
!    11 IF (N .GT. 2)                     GO TO 12
! !     NO CONDITION AT LEFT END AND N = 2.
!       C(4,1) = 1.D0
!       C(3,1) = 1.D0
!       C(2,1) = 2.D0*C(4,2)
!                                         GO TO 25
! !     NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
!    12 C(4,1) = C(3,3)
!       C(3,1) = C(3,2) + C(3,3)
!       C(2,1) = ((C(3,2)+2.D0*C(3,1))*C(4,2)*C(3,3)+C(3,2)**2*C(4,3))/C(3,1)
!                                         GO TO 19
! !     SLOPE PRESCRIBED AT LEFT END.
!    15 C(4,1) = 1.D0
!       C(3,1) = 0.D0
!                                         GO TO 18
! !     SECOND DERIVATIVE PRESCRIBED AT LEFT END.
!    16 C(4,1) = 2.D0
!       C(3,1) = 1.D0
!       C(2,1) = 3.D0*C(4,2) - C(3,2)/2.D0*C(2,1)
!    18 IF (N .EQ. 2)                     GO TO 25
! !  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESP. EQUATIONS AND CAR-
! !  RY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE M-TH
! !  EQUATION READS    C(4,M)*S(M) + C(3,M)*S(M+1) = C(2,M).
!    19 DO 20 M=2,L
!          G = -C(3,M+1)/C(4,M-1)
!          C(2,M) = G*C(2,M-1) + 3.D0*(C(3,M)*C(4,M+1)+C(3,M+1)*C(4,M))
!    20    C(4,M) = G*C(3,M-1) + 2.D0*(C(3,M) + C(3,M+1))
! !     CONSTRUCT LAST EQUATION FROM THE SECOND BOUNDARY CONDITION, OF THE FORM
! !           (-G*C(4,N-1))*S(N-1) + C(4,N)*S(N) = C(2,N)
! !     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
! !     SUBSTITUTION, SINCE C ARRAY HAPPENS TO BE SET UP JUST RIGHT FOR IT
! !     AT THIS POINT.
!       IF (IBCEND-1)                     21,30,24
!    21 IF (N .EQ. 3 .AND. IBCEND .EQ. 0) GO TO 22
! !     NOT-A-KNOT AND N .GE. 3, AND EITHER N .GT. 3 OR ALSO NOT-A-KNOT AT
! !     LEFT END POINT.
!       G = C(3,N-1) + C(3,N)
!       C(2,N) = ((C(3,N)+2.D0*G)*C(4,N)*C(3,N-1) + C(3,N)**2*(C(1,N-1)-C(1,N-2))/C(3,N-1))/G
!       G = -G/C(4,N-1)
!       C(4,N) = C(3,N-1)
!                                         GO TO 29
! !     EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND NOT-A-
! !     KNOT AT LEFT END POINT).
!    22 C(2,N) = 2.D0*C(4,N)
!       C(4,N) = 1.D0
!                                         GO TO 28
! !     SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
!    24 C(2,N) = 3.D0*C(4,N) + C(3,N)/2.D0*C(2,N)
!       C(4,N) = 2.D0
!                                         GO TO 28
!    25 IF (IBCEND-1)                     26,30,24
!    26 IF (IBCBEG .GT. 0)                GO TO 22
! !     NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
!       C(2,N) = C(4,N)
!                                         GO TO 30
!    28 G = -1.D0/C(4,N-1)
! !  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
!    29 C(4,N) = G*C(3,N-1) + C(4,N)
!       C(2,N) = (G*C(2,N-1) + C(2,N))/C(4,N)
! !  CARRY OUT BACK SUBSTITUTION
!    30 DO 40 J=L,1,-1
!    40    C(2,J) = (C(2,J) - C(3,J)*C(2,J+1))/C(4,J)
! !  GENERATE CUBIC COEFFICIENTS IN EACH INTERVAL, I.E., THE DERIV.S
! !  AT ITS LEFT ENDPOINT, FROM VALUE AND SLOPE AT ITS ENDPOINTS.
!       DO 50 I=2,N
!          DTAU = C(3,I)
!          DIVDF1 = (C(1,I) - C(1,I-1))/DTAU
!          DIVDF3 = C(2,I-1) + C(2,I) - 2.D0*DIVDF1
!          C(3,I-1) = 2.D0*(DIVDF1 - C(2,I-1) - DIVDF3)/DTAU
!    50    C(4,I-1) = (DIVDF3/DTAU)*(6.D0/DTAU)
!          RETURN
!    END subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! !--------------------------------------------------------------------!
! !  Extracted from "A practical guide to splines," 1st edition,       !
! !  Applied Mathematical Sciences 27, Carl de Boor, Springer, 1978.   !
! !--------------------------------------------------------------------!
! 
!    SUBROUTINE INTERV ( XT, LXT, X, LEFT, MFLAG )
! !  COMPUTES LEFT = MAX( I , 1 .LE. I .LE. LXT .AND. XT(I) .LE. X ).
! !
! !  XT.....A REAL SEQUENCE, OF LENGHT LXT, ASSUMED TO BE NONDECREASING
! !  LXT.....NUMBER OF TERMS IN THE SEQUENCE XT.
! !  X.....THE POINT WHOSE LOCATION WITH RESPECT TO THE SEQUENCE XT IS 
! !        TO BE DETERMINED.
! !  LEFT, MFLAG.....BOTH INTEGERS, WHOSE VALUE IS
! !
! !   1     -1      IF               X .LT. XT(1)
! !   I      0      IF   XT(I)  .LE. X .LT. XT(I+1)
! !  LXT     1      IF  XT(LXT) .LE. X
! !         IN PARTICULAR, MFLAG = 0 IS THE 'USUAL' CASE. MFLAG .NE. 0
! !         INDICATES THAT X LIES OUTSIDE THE HALFOPEN INTERVAL
! !         XT(1) .LE. Y .LT. XT(LXT). THE ASYMETRIC TREATMENT OF THE 
! !         INTERVAL IS DUE TO THE DECISION TO MAKE ALL PP FUNCTIONS CONT-
! !         INUOUS FROM THE RIGHT.
! !  THE PROGRAM IS DESIGNED TO BE EFFICIENT IN THE COMMON SITUATION THAT
! !  IT IS CALLED REPEATEDLY, WITH X TAKEN FROM AN INCREASING OR DECREA-
! !  SING SEQUENCE. THIS WILL HAPPEN, E.G., WHEN A PP FUNCTION IS TO BE
! !  GRAPHED. THE FIRST GUESS FOR LEFT IS THEREFORE TAKEN TO BE THE VAL-
! !  UE RETURNED AT THE PREVIOUS CALL AND STORED IN THE L O C A L VARIA-
! !  BLE  ILO . A FIRST CHECK ASCERTAINS THAT  ILO .LT. LXT (THIS IS NEC-
! !  ESSARY SINCE THE PRESENT CALL MAY HAVE NOTHING TO DO WITH THE PREVI- 
! !  OUS CALL). THEN, IF  XT(ILO) .LE. X .LT. XT(ILO+1), WE SET  LEFT = 
! !  ILO  AND ARE DONE AFTER JUST THREE COMPARISONS.
! !     OTHERWISE, WE REPEATEDLY DOUBLE THE DIFFERENCE  ISTEP = IHI - ILO
! !  WHILE ALSO MOVING  ILO  AND  IHI  IN THE DIRECTION OF  X, UNTIL
! !                      XT(ILO) .LE. X .LT. XT(IHI) ,
! !  AFTER WHICH WE USE BISECTION TO GET, IN ADDITION, ILO+1 = IHI .
! !  LEFT = ILO  IS THE RETURNED.
! 
!       IMPLICIT NONE
!       INTEGER  LEFT, LXT,MFLAG, IHI, ILO, ISTEP, MIDDLE
!       REAL(8) X, XT(LXT)
!       DATA ILO /1/
! !     SAVE ILO  (A VALID FORTRAN STATEMENT IN THE NEW 1977 STANDARD)
!       IHI = ILO + 1
!       IF (IHI .LT. LXT)                 GO TO 20
!          IF (X .GE. XT(LXT))            GO TO 110
!          IF (LXT .LE. 1)                GO TO 90
!          ILO = LXT - 1
!          IHI = LXT
!   20  IF (X .GE. XT(IHI))               GO TO 40
!       IF (X .GE. XT(ILO))               GO TO 100
! !     NOW X .LT. XT(ILO). DECREASE  ILO  TO CAPTURE X .
!       ISTEP = 1
!   31     IHI = ILO
!          ILO = IHI - ISTEP
!          IF (ILO .LE. 1)                GO TO 35
!          IF (X .GE. XT(ILO))            GO TO 50
!          ISTEP = ISTEP*2
!                                         GO TO 31
!   35  ILO = 1
!       IF (X .LT. XT(1))                 GO TO 90
!                                         GO TO 50
! !              **** NOW X .GE. XT(IHI). INCREASE  IHI  TO CAPTURE  X . 
!   40  ISTEP = 1
!   41     ILO = IHI
!          IHI = ILO + ISTEP
!          IF (IHI .GE. LXT)              GO TO 45
!          IF (X .LT. XT(IHI))            GO TO 50
!          ISTEP = ISTEP*2
!                                         GO TO 41
!   45  IF (X .GE. XT(LXT))               GO TO 110
!       IHI = LXT
! !           **** NOW  XT(ILO) .LE. X .LT. XT(IHI). NARROW THE INTERVAL.
!   50  MIDDLE = (ILO + IHI)/2
!       IF (MIDDLE .EQ. ILO)              GO TO 100
! !     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1 .
!       IF (X .LT. XT(MIDDLE))            GO TO 53
!          ILO = MIDDLE
!                                         GO TO 50
!   53     IHI = MIDDLE
!                                         GO TO 50
!   90  MFLAG = -1
!       LEFT = 1
!                                         RETURN
!  100  MFLAG = 0
!       LEFT = ILO
!                                         RETURN
!  110  MFLAG = 1
!       LEFT = LXT
!                                         RETURN
!     END SUBROUTINE
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! !---------------------------------------------------------------------!
! !  Extracted from "A practical guide to splines," 1st edition,        !
! !  Applied Mathematical Sciences 27, Carl de Boor, Springer, 1978.    !
! !---------------------------------------------------------------------!
! 
!    REAL(8) FUNCTION PPVALU (BREAK, COEF, L, K, X, JDERIV)
!    IMPLICIT NONE
! !  CALCULATES VALUE AT X OF JDERIV-TH DERIVATIVE OF PP FCT FROM PP-REPR
! !  TO BE EVALUATED. SPECIFICALLY, THE J-TH DERIVATIVE OF F IS      
! !  GIVEN BY
! !  (D**J)F(X) = COEF(J+1,I) + H*(COEF(J+2,I) + H*( ... (COEF(K-1,I) +
! !                             + H*COEF(K,I)/K-J-I))/(K-J-2) ... )/2)/1
! !  WITH  H = X - BREAK(I),  AND
! !  I = MAX( 1 , MAX( J , BREAK(J) .LE. X , 1 .LE. J .LE. L ) ).
! !  X.....THE POINT AT WHICH TO EVALUATE.
! !  JDERIV.....INTEGER GIVING THE ORDER OF THE DERIVATIVE TO BE EVALUAT-
! !  ED. ASSUMED TO BE ZERO OR POSITIVE.
! !  PPVALU.....THE VALUE OF THE (JDERIV)-TH DERIVATIVE OF F AT X.
! !  THE INTERVAL INDEX I, APPROPRIATE FOR X, IS FOUND THROUGHT A  
! !  CALL TO INTERV. THE FORMULA ABOVE FOR THE JDERIV-YH DERIVATIVE
! !  OF F IS THEN EVALUATED (BY NESTED MULTIPLICATION).
!       INTEGER JDERIV, K, L, I, M, NDUMMY
!       REAL(8) BREAK(L), COEF(K,L), X, FMMJDR, H
!       PPVALU = 0.D0
!       FMMJDR = K - JDERIV
! !     DERIVATIVES OF ORDER K OR HIGHER ARE IDENTICALLY ZERO.
!       IF (FMMJDR .LE. 0.D0)           GO TO 99
! !
! !              FIND INDEX I OF LARGEST BREAKPOINT TO THE LEFT OF X.
! !     CALL INTERV ( BREAK, L, X, I, NDUMMY )
! !     EVALUATE JDERIV-TH DERIVATIVE OF I-TH POLYNOMIAL PIECE AT X.
!       H = X - BREAK(I)
!       DO 10 M=K,JDERIV+1,-1
!          PPVALU = (PPVALU/FMMJDR)*H + COEF(M,I)
!   10     FMMJDR = FMMJDR - 1.D0
!   99  RETURN
!       END function
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  !-------------------------------------------------------------------------------------!
!  ! interpolate a G^qmc(-Lfak:Lfak) to G^smooth(-L:L) using the cubic-spline algorithm. !
!  !-------------------------------------------------------------------------------------!
! 
!       subroutine interp(gtmp,g,Lfak,Lfak1,L)
!       implicit real*8(a-h,o-z)
!       integer L,Lfak,Lfak1
!       double precision gtmp(-Lfak:Lfak),g(-L:L)
!       double precision xa(Lfak1),ya(4,Lfak1)
!       do 10 i=1,Lfak1
!          xa(i)=float(i-1)/float(Lfak)
!          ya(1,i) = gtmp(i-1)
!  10   continue
!       call CUBSPL(xa,ya,Lfak1,0,0)
!       do 20 i=1,L
!          x=float(i)/float(L)
!          g(i) = PPVALU(xa,ya,Lfak,4,x,0)
!  20   continue
!       g(0)=gtmp(0)
!       do 40 i=1,L
!          g(-i)=-g(L-i)
!  40   continue
!       return
!       end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
! subroutine get_derivative_noise_test
! implicit none
!     integer  :: i,k,nphas,win,npoint,vois,j
!     parameter(nphas=3,npoint=50)
!     REAL(8)  :: dxdt(nphas+1,npoint),x(npoint),tab(npoint)
! 
!       vois=3; win=3
! 
!       do i=1,npoint
!          x(i)   =  dble(i)/dble(npoint)*2.d0*pi
!          tab(i) = cos(x(i))
!          write(15,*) x(i),tab(i)
!       enddo
! 
!       write(*,*) ' get derivative noise test ' 
!       call get_derivative_noise_rout(x,tab,win,vois,nphas,dxdt)
!       write(*,*) 'plot result derivative noise, npoint = ', npoint
!       do i=1,npoint
!        write(16,'(10f8.3)') x(i),(dxdt(j,i),j=1,nphas+1)
!       enddo
! 
! stop
! end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!      !-------------------------------------!
! 
! subroutine fourier_filter_(xx,ttab,vois)
! implicit none
!    integer     :: vois,k
!    REAL(8)     :: xx(:),i(size(xx)),j(size(xx))
!    complex(8)  :: ttab(:)
!       i=real(ttab)
!       j=aimag(ttab) 
!       call fourier_filter__(xx,i,vois)
!       call fourier_filter__(xx,j,vois)
!       ttab=i+imi*j
! end subroutine
! 
!       !-------------------------------------!
! 
! subroutine fourier_filter__(xx,ttab,vois)
! implicit none
!     integer  :: i,k,npoint,vois,j 
!     parameter(npoint=2**11)
!     REAL(8)  :: x(npoint),tab(npoint)
!     REAL(8)  :: xx(:),ttab(:)
!       x=[ ( xx(1) + (xx(size(xx)) - xx(1)) * dble(i-1)/dble(npoint-1) , i=1,npoint  )  ]
!       call resampleit(xx,ttab,x,tab,0.d0)
!       call fourier_filter_noise_(npoint,tab,vois)
!       call resampleit(x,tab,xx,ttab,0.d0)
! end subroutine 
! 
!       !-------------------------------------!
! 
! subroutine get_derivative_noise_rout(xx,ttab,win,vois,nphas,ddxdt)
! implicit none
!     integer  :: i,k,nphas,win,npoint,vois,j
!     parameter(npoint=2**8)
!     REAL(8)  :: dxdt(nphas+1,npoint),x(npoint),tab(npoint)
!     REAL(8)  :: xx(:),ttab(:),ddxdt(nphas+1,size(xx))
!       x=[ ( xx(1) + (xx(size(xx)) - xx(1)) * dble(i-1)/dble(npoint-1) , i=1,npoint  )  ]
!       call resampleit(xx,ttab,x,tab,0.d0)
!       call get_derivative_noise(npoint,x,tab,nphas,win,vois,dxdt)
!       call resampleit(x,tab,xx,ttab,0.d0)
!       do j=1,nphas+1
!        call resampleit(x,dxdt(j,:),xx,ddxdt(j,:),0.d0)
!       enddo
! end subroutine
! 
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  real(8) function inflexion_point(y,x)
!  implicit none
!  integer,parameter  ::  mm=1000
!  integer            ::  i
!  real(8)            ::  x(:),y(:),xx(mm),yy(mm)
! 
!   xx(1:mm)=(/(minval(x)+dble(i-1)/dble(mm-1)*(maxval(x)-minval(x)),i=1,mm)/)
!   call derivateit(x,y,xx,yy,2,0.d0)
!   inflexion_point = xx( minloci(  abs(yy)    ) )
! 
!  end function
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  subroutine interpolate_array_(kx,ky,arrayin,xin,yin,arrayout,xout,yout,smooth)
!  implicit none
!   integer,parameter                :: nnx=70,nny=70
!   real(8)                          :: smooth,arrayin(:),xin(:),yin(:),xout(:),yout(:),arrayout(:,:)
!   real(4)                          :: tx(nnx),ty(nny),w(size(xin)),xb,xe,yb,ye,s,eps
!   integer                          :: j,ier,iopt,siz1,siz2,rr,ss,k
!   integer                          :: mx,my,ne,km,nx,ny,m,kx,ky,nxest,nyest
!   integer                          :: kwrk,nmax,i,u,v,bx,by,b1,b2,lwrk1,lwrk2
!   real(4)                          :: xr(size(xout)),yr(size(yout)),zr(size(xout)*size(yout))
!   real(4)                          :: fp,c((nnx-kx-1)*(nny-ky-1))
!   real,dimension(:),allocatable    :: work1,work2
!   integer,dimension(:),allocatable :: iwrk
! 
!   nxest=nnx
!   nyest=nny
!   m=size(arrayin(:))
!   iopt=-1
!   w=1.0
!   s=real(smooth)
!   eps=0.000001
!   nmax=max(nxest,nyest)
!   xb=minval(xin)-0.001
!   xe=maxval(xin)+0.001
!   yb=minval(yin)-0.001
!   ye=maxval(yin)+0.001
!   u = nxest-kx-1
!   v = nyest-ky-1 
!   km = max(kx,ky)+1
!   ne = max(nxest,nyest)
!   bx = kx*v+ky+1
!   by = ky*u+kx+1
!   if(bx.le.by) then
!      b1 = bx; b2 = b1+v-ky 
!    endif
!   if(bx.gt.by) then
!      b1 = by; b2 = b1+u-kx 
!   endif
! 
!   lwrk1 = 10+u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
!   lwrk2 = u*v*(b2+1)+b2
!   kwrk = 10+m+(nxest-2*kx-1)*(nyest-2*ky-1)
! 
!   tx(1:nnx)=(/( dble(i-1)/dble(nnx)*(xe-xb)+xb,i=1,nnx)/)
!   ty(1:nny)=(/( dble(i-1)/dble(nny)*(ye-yb)+yb,i=1,nny)/)
!   nx=nnx
!   ny=nny
! 
!   write(*,*) 'iopt,smooth : ', iopt,s
!   write(*,*) 'nxest,nyest : ',nxest,nyest
!   write(*,*) 'm,ne        : ', m,ne
!   write(*,*) 'kx,ky,km    : ', kx,ky,km
!   write(*,*) 'u,v         : ', u,v
!   write(*,*) 'b1,b2       : ', b1,b2
!   write(*,*) 'eps         : ', eps
!   write(*,*) 'lwrk1,lwrk2,kwrk : ', lwrk1,lwrk2,kwrk
!  
!   write(*,*) 'test on entries...'
!   write(*,*)  m>=(kx+1)*(ky+1)
!   write(*,*)  nmax>=nxest
!   write(*,*)  nmax>=nyest
!   write(*,*)  lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
!   write(*,*)  kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
! 
! 
!   allocate(work1(lwrk1),work2(lwrk2),iwrk(kwrk))
! 
!   write(*,*) 'surfit start...'
!   call surfit(iopt,m,real(xin),real(yin),real(arrayin),w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
!      &  nmax,eps,nx,tx,ny,ty,c,fp,work1,lwrk1,work2,lwrk2,iwrk,kwrk,ier)
!   write(*,*) 'surfit done...'
!   deallocate(work1,work2,iwrk)
! 
!   write(*,*) 'ERROR : ', ier
! 
!   mx=size(xout)
!   my=size(yout)
!   lwrk1 = 10+mx*(kx+1)+my*(ky+1)
!   kwrk  = mx+my+10
!   write(*,*) 'lwrk1,kwrk : ', lwrk1,kwrk
!   allocate(work1(lwrk1),iwrk(kwrk))        
!   write(*,*)'define inputs...'
! 
!   xr=real(xout)
!   yr=real(yout)
!   write(*,*) '...bispev start...'
!   write(*,*) 'x y knots number : ', nx,ny
!   write(*,*) '      max number : ', nnx,nny
!   call bispev(tx,nx,ty,ny,c,kx,ky,xr,mx,yr,my,zr,work1,lwrk1,iwrk,kwrk,ier)
!   write(*,*) '...done....'
! 
!   do i=1,mx
!    do j=1,my
!     arrayout(i,j)=zr(my*(i-1)+j)
!    enddo
!   enddo
! 
!   write(*,*) 'output done...'
!   deallocate(work1,iwrk)
!   write(*,*) 'return'
!  return
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  subroutine repair_curves_(x,y,yerr,kk,iii1,mismatch)
!  implicit none
!  real(8)      :: x(:),y(:),yerr(:),vv(2),xx,yy,dxdt(2,size(x))
!  integer      :: kk,i,j,kkk,z,ok,ii,ii1,iii1,cc,l,lll,ll
!  real(8)      :: xxx1,xxx2,der(2)
!  real(8)      :: a,deriv,deriv_tail
!  logical      :: crossing
!  real(8)      :: mismatch
! 
!  ! iii1 : 0.9 * matsu
!  ! kk   :       matsu
!  ! mismatch : in derivative to join the curve
! 
!  call get_derivative_noise_rout(x(iii1:kk),y(iii1:kk),4,4,1,dxdt(1:2,iii1:kk))
!  deriv=dxdt(2,iii1+3)
!  write(*,*) 'joining curves.....'
!  do i=iii1+8,size(x)-2
!   deriv_tail= (y(i+1)-(y(i-1)))/(x(i+1)-x(i-1))
!   if(abs(deriv_tail-deriv)<mismatch) exit
!  enddo
!  write(*,*) 'matching of derivative gave : ', deriv_tail-deriv,mismatch,i,size(x)
!  !linear interpolaiton
!    do l=iii1+3,i
!        y(l) = (y(i)-y(iii1)) * dble(l-iii1)/dble(i-iii1) + y(iii1)
!    enddo
! 
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  subroutine repair_curves(x,y,yerr,kk,iii1,mismatch)
!  implicit none
!  real(8)          :: x(:),y(:),yerr(:),vv(2),xx,yy
!  integer          :: kk,i,j,kkk,z,ok,ii,ii1,iii1,cc
!  real(8)          :: xxx1,xxx2,der(2)
!  real(8)          :: a
!  logical          :: crossing
!  real(8),optional :: mismatch
! 
!    if(kk>2.and.present(mismatch))then
!       a=(maxval(y)-minval(y))*mismatch
!       if(abs(y(kk+1)-y(kk-1))<a.and.maxval(yerr(kk-10:kk+10))<2.*a) return
!    endif
! 
!    ii1=max(2,iii1); xxx1=x(kk); xxx2=x(kk+1)
!    der=0.; cc=6
!    if(kk-cc-1<1) cc=kk-2
!    do i=0,cc
!       der(1)=der(1)+(x(max(2,kk-i))-x(max(1,kk-i-1)))/dble(cc)
!       der(2)=der(2)+(y(max(2,kk-i))-y(max(1,kk-i-1)))/dble(cc)
!    enddo
!    der=der/1.75
!    if(kk-3>0)then
!       y(kk)=(y(kk)+y(kk-1)+y(kk-2)+y(kk-3))/4.d0
!    endif
! 
!    call doit
!    if(ok==1) return
! 
!  return
! 
!  contains
! 
!   !------------------!
!   !------------------!
!   !------------------!
!   !------------------!
! 
!  subroutine doit
!  implicit none
!    call cross_or_not(der)
!    kkk=minloci(abs(xx-x))
!    ok=0
!    if(crossing.and.kkk>kk.and.abs(kkk-kk)<100)then
!     write(*,*) 'match_curves_success'
!     call join_curves
!     ok=1
!    else
!     kkk=min(kk+180,size(x))
!     call join_curves
!     ok=1
!    endif
! 
!  end subroutine
! 
!   !------------------!
!   !------------------!
!   !------------------!
!   !------------------!
! 
!  subroutine cross_or_not(der)
!  implicit none
!  real(8) :: der(2),v1(2,2),v2(2,2)
!  integer :: i
!   v1(1,1)=x(kk)
!   v1(1,2)=y(kk)
!   v1(2,1)=x(kk)+1000.*der(1)
!   v1(2,2)=y(kk)+1000.*der(2)
!   crossing=.false.
!   do i=kk+6,size(x)-1
!    v2(1,1)=x(i)
!    v2(1,2)=y(i)
!    v2(2,1)=x(i+1)
!    v2(2,2)=y(i+1)
!    call intersect(v1,v2,crossing,xx,yy,1.d-5)
!    if(crossing) exit
!   enddo
!  end subroutine
! 
!   !------------------!
!   !------------------!
!   !------------------!
!   !------------------!
! 
!   subroutine join_curves
!   implicit none
!   integer :: z,i,j
!   logical :: crossing
! 
!    vv(1) = (x(kkk)-x(kk))/dble(kkk-kk)
!    vv(2) = (y(kkk)-y(kk))/dble(kkk-kk)
!    write(*,*)  'join curves'
!    write(*,*)  'x(kkk),x(kk),y(kkk),y(kk) : ', x(kkk),x(kk),y(kkk),y(kk)
!    write(*,*)  'kk  : ', kk,kkk
!    write(*,*)  'vv  : ', vv
!    write(*,*)  'size(y),size(x):',size(y),size(x)
!    z=0
! 
!    do j=kk+1,kkk
!      z=z+1
!      y(j)=y(kk)+vv(2)*dble(z)
!    enddo
! 
!   call bendit__(ii1,kkk,x,y,yerr,0.3d0)
! 
!   return
!   end subroutine
! 
!   !------------------!
!   !------------------!
!   !------------------!
!   !------------------!
! 
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
!  
!  subroutine test_spline
!  implicit none
!  real(8)      :: x(100),y(100),xx(1000),dh(1000),yy(1000),yyy(3,1000)
!  integer      :: i,j,k
!  type(spline) :: rspline
! 
!    call rand_init(iseed=1234)
!    do i=1,100
!      x(i)=dble(i)/100.d0*2.d0*pi
!      y(i)=cos(x(i))+(-1.d0+2.d0*drand1())/4.
!      write(44,*) x(i),y(i)
!    enddo
! 
!   call init_spline(rspline,x,y,3)
!   xx=(/( 2.d0*pi*dble(i-1)/999.d0 , i =1,1000 )/)
! 
!   call evaluate_spline(rspline,xx,yy)
!   call derivate_spline(rspline,xx,yyy,3)
! 
!   call smoothit(x,y,1.6d0)
!   do i=1,100
!    write(48,*) x(i),y(i)
!   enddo
! 
!   do i=1,100
!     y(i)=cos(x(i))+(-1.d0+2.d0*drand1())/4.d0
!   enddo
! 
!   call derivateit(x,y,y,1,1.6d0)
!   do i=1,100
!    write(49,*) x(i),y(i)
!   enddo
! 
!   do i=1,1000
!    write(45,*) xx(i),yy(i)
!    write(46,'(4f12.4)') xx(i),(yyy(k,i),k=1,3)
!   enddo
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  !--------------------!
! 
!  subroutine derivateit__(x,yy,yyy,k,ss)
!  implicit none
!  real(8)      :: x(:),yy(:,:),yyy(:,:),ss
!  integer      :: i,j,k
!  type(spline) :: rspline
!   do i=1,size(yy(:,1))
!    call derivateit_(x,yy(i,:),yyy(i,:),k,ss)
!   enddo
!  end subroutine
! 
!  !--------------------!
! 
!  subroutine derivateit_(x,y,yy,k,ss)
!  implicit none
!  real(8)      :: x(:),y(:),yy(:),ss
!  integer      :: i,j,k
!  type(spline) :: rspline
!    call init_spline(rspline,x,y,5,ss)
!    call derivate_spline(rspline)
!    yy=rspline%der(min(rspline%order,k+1),:)
!    call kill_spline(rspline)
!  end subroutine
! 
!  !--------------------!
! 
!  subroutine derivateit___(x,y,xx,yy,k,ss)
!  implicit none
!  real(8)      :: x(:),y(:),xx(:),yy(:),ss
!  integer      :: i,j,k
!  type(spline) :: rspline
!     call init_spline(rspline,x,y,5,ss)
!     call derivate_splined(rspline,xx,yy,k)
!     call kill_spline(rspline)
!  end subroutine
! 
!  !--------------------!
! 
!  subroutine derivateit____(x,y,x_,y_,k,ss)
!  implicit none
!  real(8)      :: x(:),y(:),x_,y_,xx(1),yy(1),ss
!  integer      :: i,j,k
!  type(spline) :: rspline
!     xx=x_
!     call init_spline(rspline,x,y,5,ss)
!     call derivate_splined(rspline,xx,yy,k)
!     call kill_spline(rspline)
!     y_=yy(1)
!  end subroutine
! 
!  !--------------------!
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  subroutine bendit_(i,j,x,y,yerr,ss)
!  implicit none
!  real(8)      :: x(:),y(:),ss,yerr(:)
!  real(4)      :: ys(size(x)),temp(size(x)),yp(size(x)),errors(size(x))
!  real(4)      :: td(size(x)),tsd1(size(x)),hd(size(x)),hsd1(size(x)),hsd2(size(x))
!  real(4)      :: rd(size(x)),rsd1(size(x)),rsd2(size(x)),v(size(x))
!  integer      :: i,j,k,ierr
!  real(4)      :: sss,aa,bb,dd
! 
!   dd=100000.
!   aa=0.
!   bb=100000.
! 
!   errors(1:i)=dd
!   errors(i+1:j)=aa
!   errors(j+1:size(x))=bb
! 
!   sss=2.d0/float(size(x))
!   call curvss (size(x),real(x),real(y),errors,1,real(ss),sss, &
!              & ys,yp,1.,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ierr)
!   if(ierr==1.or.ierr==0) y=ys
! 
!  return
!  end subroutine
! 
!  !--------------------!
! 
!  subroutine bendit__(i,j,x,y,yerr,ss)
!  implicit none
!  real(8)      :: x(:),y(:),yerr(:),ss,errors(size(x)),yout(size(y))
!  integer      :: i,j,k,NDEG,STATUS,IERR
!  real(8)      :: eps,A(7*size(x)+3)
!  real(4)      :: sss,aa,bb,cc,dd
! 
!   eps =  -1.
!   dd  =  1000.
!   aa  =  3. 
!   bb  =  1000.
!   errors(1:i)=dd
!   errors(i+1:j)=aa
!   errors(j+1:size(x))=bb 
!   if(messages3) write(*,*) 'CALL TO SLATEC IN BENDING FUNCTION'
!   call DPOLFT(size(x),x,y,errors,3,NDEG,eps,yout,IERR,A,STATUS)
!   if(IERR<=1) then
!    y=yout 
!   else
!    write(*,*) 'BEND IT IERR : ', IERR
!   endif
! 
!  end subroutine
! 
!  !--------------------!
! 
!  subroutine bendit(i,j,x,y,yerr,ss)
!  implicit none
!  real(8)      :: x(:),y(:),yerr(:),ss,errors(size(x))
!  integer      :: i,j,k
!  type(spline) :: rspline
!  real(4)      :: aa,bb,cc,dd
! 
!   dd=  0.0001
!   aa=  maxval(abs(y))/5.
!   bb=  0.0001
! 
!   errors(1:i)        = dd
!   errors(i:j/2)      = (/( dd + (aa-dd)*dble(k-i)/dble(j/2-1) , k=i,j/2 )/)
!   errors(j/2+1:j)    = aa
!   errors(j+1:size(x))= bb
!   
!   call init_spline(rspline,x,y,3,ss,errors=errors)
!   call evaluate_spline(rspline)
!   y=rspline%smooth
!   call kill_spline(rspline)
! 
!  return
!  end subroutine
! 
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

   !---------------------------------------!

 subroutine spline_inter_extrapolate_array__(xin,yin,x,y)

 implicit none
 real(8) :: xin(:),yin(:),x(:),y(:)
 real(8) :: y_arr(1)
 integer :: i,j,k,l,siz,ndim,ninter
  siz=size(xin(:))
  ! ebl: converting to be compatible with spline overhauser
  do i=1,size(x)
   call spline_overhauser_val ( 1, siz, xin, [yin], x(i), y_arr)
   y(i) = y_arr(1)
  enddo
 end subroutine

   !---------------------------------------!

 subroutine spline_inter_extrapolate_array_(xin,yin,x,y)
 implicit none
 real(4) :: xin(:),yin(:),x(:),y(:)
 real(8) :: yo(1)
 integer :: i,j,k,l,siz,ndim,ninter
  siz=size(xin(:))
  do i=1,size(x)
   call spline_overhauser_val ( 1, siz, dble(xin), dble(yin), dble(x(i)), yo)
   y(i)=yo(1)
  enddo
 end subroutine

   !---------------------------------------!

 subroutine spline_inter_extrapolate__(xin,yin,x,y)
 implicit none
 real(8) :: xin(:),yin(:),x,y, y_arr(1)
 integer :: i,j,k,l,siz,ndim,ninter
  siz=size(xin(:))
  call spline_overhauser_val ( 1, siz, xin, [yin], x, y_arr)  
  y = y_arr(1)
 end subroutine

   !---------------------------------------!

 subroutine spline_inter_extrapolate_(xin,yin,x,y)
 implicit none
 real(4) :: xin(:),yin(:),x,y
 real(8) :: yo(1)
 integer :: i,j,k,l,siz,ndim,ninter
  siz=size(xin(:))
  call spline_overhauser_val ( 1, siz, dble(xin), [dble(yin)], dble(x), yo)
  y=yo(1)
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
! 
!  subroutine bezier_it(xin,t,x)
!  implicit none
!  real(8) :: xin(:),x(1,1),t,tt(1),xx(1,size(xin))
!  integer :: i,j,k,l,siz,ndim,ninter
! 
!   tt=t
!   xx(1,:)=xin
!   siz=size(xin(:))
!   ndim=1
!   ninter=(siz-1)/3
!   call spline_bezier_val (ndim, ninter, xx, 1 , tt, x )
! 
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  !--------------------!
!  !--------------------!
! 
!  subroutine smoothit____(x,y,ss,errors)
!  implicit none
!  real(8)      :: x(:),y(:),ss,errors(:)
!  real(4)      :: ys(size(x)),temp(size(x)),yp(size(x)),errors2(size(x))
!  real(4)      :: td(size(x)),tsd1(size(x)),hd(size(x)),hsd1(size(x)),hsd2(size(x))
!  real(4)      :: rd(size(x)),rsd1(size(x)),rsd2(size(x)),v(size(x))
!  integer      :: i,j,k,ierr
!  real(4)      :: sss
! 
!   if(maxval(abs(x))<1.d-13.or.maxval(abs(y))<1.d-13)then
!    return
!   endif
!   call check_x_y(x,y)
! 
!   sss=2.d0/float(size(x))
!   errors2=abs(errors)
!   call curvss (size(x),real(x),real(y),errors2,1,real(ss),sss, &
!              & ys,yp,1.,td,tsd1,hd,hsd1,hsd2,rd,rsd1,rsd2,v,ierr)
!   if(ierr==1.or.ierr==0) y=ys
! 
!  return
!  end subroutine
! 
!  !--------------------!
!  !--------------------!
! 
!  subroutine smoothitb_____(x,y,errors,degree)
!  implicit none
!  real(8)          :: x(:),y(:),ss,errors(:),errors2(size(errors)),yout(size(y))
!  integer          :: i,j,k,NDEG,STATUS,IERR,degree_
!  real(8)          :: eps,A(7*size(x)+3),aa,rrr
!  integer,optional :: degree
! 
!   if(maxval(abs(x))<1.d-13.or.maxval(abs(y))<1.d-13)then
!    return
!   endif
!   call check_x_y(x,y)
! 
!   eps=-1.d0
!   rrr=max(1.d-2,maxval(abs(errors)))
!   errors2=abs(errors/rrr)
!   where(errors2<1.d-2) errors2=1.d-2
!   errors2=1.d0/(errors2)
!   NDEG=5
!   if(constantarray(y)) return
!   if(present(degree))then
!    degree_=degree
!   else
!    degree_=min(size(x),max(5,size(x)/3-3))
!   endif
!   if(messages3) then
!   write(*,*) 'CALLING SLATER DPOLFT with MAX Y,ERR : ', maxval(abs(y)),maxval(abs(errors2))
!   write(*,*) 'WITH DEGREE = ' , degree_
!   write(*,*) 'NUMBER OF POINTS = ', size(x)
!   endif
!   call DPOLFT(size(x),x,y,errors2,degree_,NDEG,eps,yout,IERR,A,STATUS)
!   if(IERR>1)then
!     write(*,*) 'IERR smoothing with polynomial : ', IERR
!     if(strongstop) stop 'smooth routine error, not a problem the smoothing was simply not applied'
!   else
!     write(*,*) 'SUCESSFULL SMOOTHING'
!     y=yout
!   endif
! 
!  return
!  end subroutine
! 
!  !--------------------!
!  !--------------------!
! 
!  subroutine smoothitbb_____(x,y,errors,fac,kkk)
!  implicit none
!  integer          :: fac
!  real(8)          :: x(:),y(:),errors(:),x2(kkk),y2(kkk),errors2(kkk),dmax
!  real(8)          :: x3(size(x)),y3(size(y)),errors3(size(errors)),back(size(x))
!  integer          :: i,ii,order(size(errors)),kkk,ikil,itot,kkk2
! 
!   if(maxval(abs(x))<1.d-13.or.maxval(abs(y))<1.d-13)then
!    return
!   endif
!   call check_x_y(x,y)
! 
!  if(maxval(abs(errors))<1.d-7)then
!   errors=0.d0
!   errors(kkk/2:kkk)=0.001
!  endif
!  if(kkk+1>size(x))then
!   stop 'error smoothitbb_____ input kkk bigger than array size'  
!  endif
! 
!  x2(1:kkk)      = x(1:kkk)
!  y2(1:kkk)      = y(1:kkk)
!  errors2(1:kkk) = errors(1:kkk)
! 
!  call qsort_array(errors2,order)
!  call qsort_adj_array(x2,order)
!  call qsort_adj_array(y2,order)
!  
!  if(fac==0) then
!   ii   = kkk
!  else
!   ii   = kkk/fac
!  endif
!  ikil = kkk - ii
!  itot = size(x)-ikil
! 
!  x3(1:ii)      = x2(1:ii)
!  y3(1:ii)      = y2(1:ii)
!  errors3(1:ii) = errors2(1:ii)
!  
!  !==================!
!  if(fac>1)then
!     kkk2=ii + max(2,(kkk-ii)/2) -1
!     do i=ii+1,kkk2
!        x3(i) = dble(i-ii) / dble(kkk2-ii-1) * (x(kkk+1)-x2(ii))  +  x2(ii)
!        y3(i) = dble(i-ii) / dble(kkk2-ii-1) * (y(kkk+1)-y2(ii))  +  y2(ii)
!     enddo 
!     itot=size(x)
!     ii=kkk2
!     errors3(ii+1:ii+kkk2)=0.03
!  else
!     itot=size(x)
!     ii=kkk
!  endif
!  !==================!
! 
!  x3(ii+1:itot)      = x(ii+1:size(x))
!  y3(ii+1:itot)      = y(ii+1:size(x))
!  errors3(ii+1:itot) = errors(ii+1:size(x))
! 
!  call qsort_array(x3(1:itot),order(1:itot))
!  call qsort_adj_array(y3(1:itot),order(1:itot))
!  call qsort_adj_array(errors3(1:itot),order(1:itot))
!  where(x3>x(kkk)) errors3=0.d0
!  call smoothitb_____(x3(1:itot),y3(1:itot),errors3(1:itot))
!  back=y
!  call resampleit(x3(1:itot),y3(1:itot),x(:),y(:),0.d0)
!  y(kkk+1:size(x))=back(kkk+1:size(x))
! 
!  return
!  end subroutine
! 
!  !--------------------!
!  !--------------------!
!  !--------------------!
!  !--------------------!
! 
!  subroutine smoothitfft(x,y,win)
!  implicit none
!  real(8)          :: x(:),y(:),win
!  integer          :: i,j
!  real(4)          :: xx(2048),yy(2048)
! 
!   if(maxval(abs(x))<1.d-13.or.maxval(abs(y))<1.d-13)then
!    return
!   endif
!   call check_x_y(x,y)
!  
!   write(*,*) 'smooth fft start'
!   j=2048-2*NINT(win)
!   xx(1:j)=(/(real(x(1)+dble(i-1)/dble(j-1)*(x(size(x))-x(1))),i=1,j)/)
!   call resampleit(real(x),real(y),xx(1:j),yy(1:j),0.d0)
!   call SMOOFT(yy,j,real(win),2048)
!   write(*,*) 'smooth done'
!   call resampleit(dble(xx(1:j)),dble(yy(1:j)),x,y,0.d0) 
!   write(*,*) 'resampling done, now return'
! 
!  return
!  end subroutine
! 
!  !--------------------!
!  !--------------------!
! 
!  subroutine smoothitred(x,y,ss,errors,fac,kkk)
!  implicit none
!  integer          :: fac
!  real(8)          :: x(:),y(:),errors(:),x2(kkk),y2(kkk),errors2(kkk),dmax
!  real(8)          :: x3(size(x)),y3(size(y)),errors3(size(errors)),ss
!  integer          :: degree,i,ii,order(size(errors)),kkk,ikil,itot
!  integer          :: ddeg
!  type(spline)     :: rspline
! 
!   if(maxval(abs(x))<1.d-13.or.maxval(abs(y))<1.d-13)then
!    return
!   endif
!   call check_x_y(x,y)
! 
!  x2(1:kkk)=x(1:kkk)
!  y2(1:kkk)=y(1:kkk)
!  errors2(1:kkk)=errors(1:kkk)
!  call qsort_array(errors2,order)
!  call qsort_adj_array(x2,order)
!  call qsort_adj_array(y2,order)
! 
!  ii   = kkk/fac
!  ikil = kkk - ii
!  itot = size(x)-ikil
! 
!  x3(1:ii)=x2(1:ii)
!  y3(1:ii)=y2(1:ii)
!  errors3(1:ii)=errors2(1:ii)
!  
!  x3(ii+1:itot)=x(kkk+1:size(x))
!  y3(ii+1:itot)=y(kkk+1:size(x))
!  errors3(ii+1:itot)=errors(kkk+1:size(x))
! 
!  call qsort_array(x3(1:itot),order(1:itot))
!  call qsort_adj_array(y3(1:itot),order(1:itot))
!  call qsort_adj_array(errors3(1:itot),order(1:itot))
! 
!   call init_spline(rspline,x3(1:ii),y3(1:ii),3,ss,errors3(1:ii))
!   call evaluate_spline(rspline)
!   y3(1:ii)=rspline%smooth(1:ii)
!   call kill_spline(rspline)
!   call resampleit(x3(1:itot),y3(1:itot),x(:),y(:),0.d0)
!  return
!  end subroutine
! 
!  !--------------------!
!  !--------------------!
! 
!  subroutine smoothit___(x,y,ss,errors)
!  implicit none
!  real(8)      :: x(:),y(:),ss,errors(:)
!  integer      :: i,j,k
!  type(spline) :: rspline
! 
!   if(maxval(abs(x))<1.d-13.or.maxval(abs(y))<1.d-13)then
!    return
!   endif
!   call check_x_y(x,y)
!   call init_spline(rspline,x,y,3,ss,errors)
!   call evaluate_spline(rspline)
!   y=rspline%smooth
!   call kill_spline(rspline)
!  end subroutine
! 
!  !--------------------!
!  !--------------------!
! 
!  subroutine smoothit__(x,yy,ss)   
!  implicit none
!  real(8)      :: x(:),yy(:,:),ss
!  integer      :: i,j,k
!  type(spline) :: rspline
!   do i=1,size(yy(:,1))
!    call smoothit_(x,yy(i,:),ss)
!   enddo
!  end subroutine
! 
!  !--------------------!
!  !--------------------!
! 
!  subroutine smoothit_(x,y,ss)
!  implicit none
!  real(8)       :: x(:),y(:),ss
!  integer      :: i,j,k
!  type(spline) :: rspline
!   if(maxval(abs(x))<1.d-13.or.maxval(abs(y))<1.d-13)then
!    return
!   endif
!   call check_x_y(x,y)
!   call init_spline(rspline,x,y,3,ss)
!   call evaluate_spline(rspline)
!   y=rspline%smooth
!   call kill_spline(rspline)
!  end subroutine
! 
!  !--------------------!
!  !--------------------!
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  subroutine prolong_fix_deriv(x,y,xin,yin)
!  implicit none
!  integer      :: ns
!  real(8)      :: step,deriv,xnew,ynew
!  real(8)      :: xin(:),yin(:),x(size(xin)+2),y(size(yin)+2)
!   
!   ns    =  size(xin)
!   step  =  xin(2)-xin(1)
!   deriv = (yin(2)-yin(1))/step
!   xnew  =  xin(1)-step
!   ynew  =  yin(1)-step*deriv
!   x(1)=xnew
!   y(1)=ynew
! 
!   x(2:ns+1)=xin
!   y(2:ns+1)=yin
! 
!   step  =  xin(ns)-xin(ns-1)
!   deriv = (yin(ns)-yin(ns-1))/step
!   xnew  =  xin(ns)+step
!   ynew  =  yin(ns)+step*deriv
!   x(ns+2)=xnew
!   y(ns+2)=ynew
! 
!  end subroutine
! 
  !--------------------------------!

 subroutine prolong_fix_deriv_r(x,y,xin,yin)
 implicit none
 integer   :: ns
 real(4)   :: step,deriv,xnew,ynew
 real(4)   :: xin(:),yin(:),x(size(xin)+2),y(size(yin)+2)

  ns    =  size(xin)
  step  =  xin(2)-xin(1)
  deriv = (yin(2)-yin(1))/step
  xnew  =  xin(1)-step
  ynew  =  yin(1)-step*deriv
  x(1)=xnew
  y(1)=ynew
  x(2:ns+1)=xin
  y(2:ns+1)=yin
  step  =  xin(ns)-xin(ns-1)
  deriv = (yin(ns)-yin(ns-1))/step
  xnew  =  xin(ns)+step
  ynew  =  yin(ns)+step*deriv
  x(ns+2)=xnew
  y(ns+2)=ynew

 end subroutine

  !--------------------------------!

 subroutine prolong_fix_deriv_c(x,y,xin,yin)
 implicit none
 integer    :: ns
 real(8)    :: step,xnew
 real(8)    :: xin(:),x(size(xin)+2)
 complex(8) :: yin(:),y(size(yin)+2),ynew,deriv

  ns    =  size(xin)
  step  =  xin(2)-xin(1)
  deriv = (yin(2)-yin(1))/step
  xnew  =  xin(1)-step
  ynew  =  yin(1)-step*deriv
  x(1)=xnew
  y(1)=ynew
  x(2:ns+1)=xin
  y(2:ns+1)=yin
  step  =  xin(ns)-xin(ns-1)
  deriv = (yin(ns)-yin(ns-1))/step
  xnew  =  xin(ns)+step
  ynew  =  yin(ns)+step*deriv
  x(ns+2)=xnew
  y(ns+2)=ynew

 end subroutine

  !--------------------------------!

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

!  subroutine resampleit_matrix___(xin,yin,xx,yy,ss)
!  implicit none
!  real(8)    :: xin(:),xx(:),trace,diag_(size(xin)),vlam,w(size(xin)),sy(size(xin))
!  complex(8) :: yin(:,:,:,:)
!  complex(8) :: yy(:,:,:,:)
!  real(8)    :: ss,tt(size(yin(1,1,1,:))),ttr(size(yin(1,1,1,:)))
!  integer    :: i,j,k,u(4),ierr,v(3)
! 
!  if(size(xin)==size(xx))then
!   if(maxval(abs(xin-xx))<1.d-12) then
!    yy=yin
!    return
!   endif
!  endif
!  if(maxval(abs(yin))<1.d-15) then
!    yy=0.d0
!    return
!  endif
! 
!  v(1)=2; v(2)=3; v(3)=1; w=1.d0; yy=0.; u=shape(yin)
!  do k=1,u(3)
!   do i=1,u(1)
!    do j=1,u(2)
!     call css(-2000.d0,xin,real(yin(i,j,k,:)),w,sy,trace,diag_,vlam,size(xx),xx,ttr,v,0,ierr)
!     call css(-2000.d0,xin,aimag(yin(i,j,k,:)),w,sy,trace,diag_,vlam,size(xx),xx,tt,v,0,ierr)
!     yy(i,j,k,:)=CMPLX(ttr,tt,kind=8)
!    enddo
!   enddo
!  enddo
! 
!  return
!  end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit_matrix_(xin,yin,xx,yy,ss)
 implicit none
 real(8)      :: xin(:),xx(:)
 complex(8)   :: yin(:,:,:,:)
 complex(8)   :: yy(:,:,:,:)
 real(8)      :: ss,tt(size(yin(1,1,1,:))),ttr(size(yin(1,1,1,:)))
 integer      :: i,j,k,u(4)
 type(spline) :: rspline,ispline

 if(size(xin)==size(xx))then
 if(maxval(abs(xin-xx))<1.d-12) then
  yy=yin
  return
 endif
 endif
 if(maxval(abs(yin))<1.d-15) then
   yy=0.d0
   return
 endif

 yy=0.; u=shape(yin)
 do k=1,u(3)
  do i=1,u(1)
   do j=1,u(2)
     if(maxval(abs(yin(i,j,k,:)))>1.d-11)then
      call init_spline(rspline,xin,real(yin(i,j,k,:)),3,ss)
      call init_spline(ispline,xin,aimag(yin(i,j,k,:)),3,ss)
      call evaluate_spline(rspline,xx,ttr)
      call evaluate_spline(ispline,xx,tt)
      yy(i,j,k,:)=CMPLX(ttr,tt,kind=8)
     endif
   enddo
  enddo
 enddo

 call kill_spline(rspline)
 call kill_spline(ispline)

 return
 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit_matrix__(xin,yin,xx,yy,ss)
 implicit none
 real(8)      :: xin(:),xx(:)
 complex(8)   :: yin(:,:,:,:)
 complex(8)   :: yy(:,:,:,:)
 real(8)      :: ss,tt,ttr
 real(8)      :: tt_arr(1),ttr_arr(1)
 integer      :: i,j,k,u(4),l
 type(spline) :: rspline,ispline

 if(size(xx)==size(xin))then
 if(maxval(abs(xin-xx))<1.d-12) then
  yy=yin
  return
 endif
 endif
 if(maxval(abs(yin))<1.d-15) then
   yy=0.d0
   return
 endif

 yy=0.; u=shape(yin)
 do k=1,u(3)
 do i=1,u(1)
 do j=1,u(2)
  do l=1,size(xx)
   call spline_overhauser_val ( 1, size(xin), xin, [real(yin(i,j,k,:))], xx(l), ttr_arr)
   call spline_overhauser_val ( 1, size(xin), xin, [aimag(yin(i,j,k,:))], xx(l), tt_arr)
   ttr = ttr_arr(1)
   tt = tt_arr(1)
   yy(i,j,k,l)=CMPLX(ttr,tt,kind=8)
  enddo
 enddo
 enddo
 enddo

 return
 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit_matrix(xin,yin,xx,yy,ss)
 implicit none
 real(8)      :: xin(:),xx(:)
 complex(8)   :: yin(:,:,:)
 complex(8)   :: yy(:,:,:)
 real(8)      :: ss,tt(size(xx)),ttr(size(xx))
 integer      :: i,j,k,u(3)
 type(spline) :: rspline,ispline

  if(maxval(abs(xin))<1.d-10)then
    write(*,*) 'resample matrix error, xin is 0'
    stop
  endif

  if(size(xx)==size(xin))then
   if(maxval(abs(xin-xx))<1.d-12) then
    if(size(xin)/=size(xx))then
     write(*,*) 'strange... scales xin and x are the same, but different sizes'
     write(*,*) 'size  xin : ', size(xin)
     write(*,*) 'size xout : ', size(xx)
     stop
    endif
    yy=yin
    return
   endif
  endif

  u=shape(yin)
  do i=1,u(1)
   do j=1,u(2)
    if(maxval(abs(yin(i,j,:)))>1.d-11)then
      call init_spline(rspline,xin,real(yin(i,j,:)),3,ss)
      call init_spline(ispline,xin,aimag(yin(i,j,:)),3,ss)
      call evaluate_spline(rspline,xx,ttr)
      call evaluate_spline(ispline,xx,tt)
      yy(i,j,:)=CMPLX(ttr,tt,kind=8)
    endif
   enddo
  enddo
  call kill_spline(rspline)
  call kill_spline(ispline)

 end subroutine


  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
! 
!  subroutine linear_interpolation(xin,yin,xx,yy)
!  real(8)      :: xin(:),yin(:)
!  real(8)      :: xx(:),yy(:),fac
!  integer      :: i,j,k
!    do i=1,size(xx)
!        j=minloci(abs(xx(i)-xin))
!        if(j<size(yin))then
!         fac=(xx(i)-xin(j))/(xin(j+1)-xin(j))
!         yy(i) = (yin(j+1)-yin(j)) * fac + yin(j)
!        else
!         yy(i) = yin(j)
!        endif
!    enddo
!  end subroutine
! 
!   !--------------------------------!
!   !--------------------------------!
!   !--------------------------------!
!   !--------------------------------!
!   !--------------------------------!
! 
!  subroutine resampleit__(xin,yin,xx,yy,ss)
!  implicit none
!  real(8)      :: xin(:),yin(:),x(size(xin)+2),y(size(yin)+2)
!  real(8)      :: xx(:),yy(:),ss
!  integer      :: i,j,k
!  type(spline) :: rspline
! 
!   if(size(xin)==size(xx))then
!   if(maxval(abs(xin-xx))<1.d-12) then
!    yy=yin
!    return
!   endif
!   endif
!   if(maxval(abs(yin))<1.d-15) then
!    yy=0.d0
!    return
!   endif
! 
!   call prolong_fix_deriv(x,y,xin,yin)
!   call init_spline(rspline,x,y,4,ss)
!   call evaluate_splineb(rspline,xx,yy)
!   call kill_spline(rspline)
! 
!  return
!  end subroutine
! 
!   !--------------------------------!
!   !--------------------------------!
!   !--------------------------------!
!   !--------------------------------!
!   !--------------------------------!
! 
!  subroutine resampleit_lin_(xin,yin,xx,yy)
!  implicit none
!  real(8)      :: xin(:),yin(:)
!  real(8)      :: xx(:),yy(:),ss,a,b,vv
!  integer      :: i,j,k,iout
! 
!   if(size(xx)==size(xin))then
!   if(maxval(abs(xin-xx))<1.d-12) then
!    yy=yin
!    return
!   endif
!   endif
!  if(maxval(abs(yin))<1.d-15) then
!    yy=0.d0
!    return
!  endif
! 
!   a  = minval(xin)
!   b  = maxval(xin)
!   j  = size(xx)
!   xx = (/( a + (b-a) * dble(i-1)/dble(j-1),i=1,j )/)
! 
!   do i=1,size(xx)
!    call findin1Dmesh_(xin,size(xin),xx(i),iout)
!    if(iout+1>size(xin) ) yy(i)=yin(size(yin))
!    if(iout<1           ) yy(i)=yin(1)
!    if(iout>=1.and.iout+1<=size(xin))then
!     vv    = (yin(iout+1)-yin(iout))/(xin(iout+1)-xin(iout))
!     yy(i) =  yin(iout) + (xx(i)-xin(iout))*vv
!    endif
!   enddo
! 
!   return
!   contains
!      !------------------------------------!
!      subroutine findin1Dmesh_(xin,n,xx,iout)
!      implicit none
!       integer :: n,iout,jjj
!       real(8) :: xin(n),xx
!       
!       iout=0
!       do jjj=1,n-1
!        if(xx>=xin(jjj).and.xx<xin(jjj+1))then
!          iout=jjj
!          exit
!        endif
!       enddo
!       if(xx>xin(n-1)) iout=n
!       if(xx<xin(1))   iout=0
!    
!      end subroutine
!      !------------------------------------!
!  end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit_(xin,yin,xx,yy,ss)
 implicit none
 real(4)      :: xin(:),yin(:),x(size(xin)+2),y(size(yin)+2),xx(:),yy(:)
 real(8)      :: ss
 integer      :: i,j,k
 type(spline) :: rspline

  if(size(xin)==size(xx))then
  if(maxval(abs(xin-xx))<1.d-12) then
   yy=yin
   return
  endif
  endif
 if(maxval(abs(yin))<1.d-15) then
   yy=0.d0
   return
 endif

  call prolong_fix_deriv_r(x,y,xin,yin)
  call init_spline(rspline,dble(x),dble(y),4,ss)
  call evaluate_splinebr(rspline,xx,yy)
  call kill_spline(rspline)

 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit___(xin,yin,xx,yy,ss)
 implicit none
 complex(8)   :: yin(:),y(size(yin)+2),yy(:)
 real(8)      :: xin(:),x(size(xin)+2),xx(:)
 real(8)      :: ss,yyy(size(yy))
 integer      :: i,j,k
 type(spline) :: rspline

  if(size(xin)==size(xx))then
  if(maxval(abs(xin-xx))<1.d-12) then
   yy=yin
   return
  endif
  endif
 if(maxval(abs(yin))<1.d-15) then
   yy=0.d0
   return
 endif

  call prolong_fix_deriv_c(x,y,xin,yin)
  call init_spline(rspline,x,real(y),4,ss)
  call evaluate_splineb(rspline,xx,yyy)
  yy=yyy
  call kill_spline(rspline)
  call init_spline(rspline,x,aimag(y),4,ss)
  call evaluate_splineb(rspline,xx,yyy)
  yy=yy+imi*yyy
  call kill_spline(rspline)

 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit__b(x,y,xx,yy,ss)
 implicit none
 real(8)      :: x(:),y(:)
 real(8)      :: xx(:),yy(:),ss
 integer      :: i,j,k
 type(spline) :: rspline

  if(size(x)==size(xx))then
  if(maxval(abs(x-xx))<1.d-12) then
   yy=y
   return
  endif
  endif

 if(maxval(abs(y))<1.d-15) then
   yy=0.d0
   return
 endif

  call init_spline(rspline,x,y,4,ss)
  call evaluate_splineb(rspline,xx,yy)
  call kill_spline(rspline)

 return
 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit_b(x,y,xx,yy,ss)
 implicit none
 real(4)      :: x(:),y(:),xx(:),yy(:)
 real(8)      :: ss
 integer      :: i,j,k
 type(spline) :: rspline

  if(size(x)==size(xx))then
  if(maxval(abs(x-xx))<1.d-12) then
   yy=y
   return
  endif
  endif
 if(maxval(abs(y))<1.d-15) then
   yy=0.d0
   return
 endif

  call init_spline(rspline,dble(x),dble(y),4,ss)
  call evaluate_splinebr(rspline,xx,yy)
  call kill_spline(rspline)

 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit___b(x,y,xx,yy,ss)
 implicit none
 complex(8)   :: y(:),yy(:)
 real(8)      :: x(:),xx(:)
 real(8)      :: ss,yyy(size(yy))
 integer      :: i,j,k
 type(spline) :: rspline

  if(size(x)==size(xx))then
  if(maxval(abs(x-xx))<1.d-12) then
   yy=y
   return
  endif
  endif
 if(maxval(abs(y))<1.d-15) then
   yy=0.d0
   return
 endif

  call init_spline(rspline,x,real(y),4,ss)
  call evaluate_splineb(rspline,xx,yyy)
  yy=yyy
  call kill_spline(rspline)
  call init_spline(rspline,x,aimag(y),4,ss)
  call evaluate_splineb(rspline,xx,yyy)
  yy=yy+imi*yyy
  call kill_spline(rspline)

 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit_xonly_r(x,y,xx,ss)
 implicit none
 real(8)    :: y(:)
 real(8)    :: yy(size(y))
 real(8)    :: x(:),xx(:),ss
 integer    :: i,j,k
  if(size(x)==size(xx))then
  if(maxval(abs(x-xx))<1.d-12) then
   return
  endif
  endif
 if(maxval(abs(y))<1.d-15) then
   return
 endif

  if(size(xx)/=size(x)) stop 'resampleit_xonly : shape of arrays differ'    
  call spline_inter_extrapolate(x,y,xx,yy) 
  y=yy

 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

 subroutine resampleit_xonly_c(x,y,xx,ss)
 implicit none
 complex(8) :: y(:)
 real(8)    :: yy(size(y)),yyi(size(y))
 real(8)    :: x(:),xx(:),ss
 integer    :: i,j,k
  if(size(x)==size(xx))then
  if(maxval(abs(x-xx))<1.d-12) then
   return
  endif
  endif
 if(maxval(abs(y))<1.d-15) return

  if(size(xx)/=size(x)) stop 'resampleit_xonly : shape of arrays differ'
  call spline_inter_extrapolate( x, real(y),xx,yy  )
  call spline_inter_extrapolate( x,aimag(y),xx,yyi )
  y=yy+imi*yyi

 end subroutine

  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!
  !--------------------------------!

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine init_spline(rspline,omi,Frc,kkk2,s2,errors)
  implicit none
    !--------------------------------------------------------------!
    ! omi : mesh x, Frc : array to fit , kk : order of the splines
    ! nn  : number of knots
    ! tt  : knots
    ! cc  : coefficient of spline
    ! w   : weights
    !--------------------------------------------------------------!
  real(8)               :: omi(:),Frc(size(omi))
  type(spline)          :: rspline 
  integer               :: kk,kkk,i,j,kkk2
  real(8),optional      :: s2
  real(8),optional      :: errors(size(omi))
  real(4)               :: s
  integer               :: ier,nn,lwrk,k,tryit
  real(4)               :: fp,xl,xr

  rspline%ierr=0
  kk=size(omi)

  if(rspline%initialized>0) then
     if(testing) write(*,*) 'danger :  spline already initialized'
     goto 21
  endif
  rspline%initialized=1

  if(allocated(rspline%om))     deallocate(rspline%om)
  if(allocated(rspline%Frc))    deallocate(rspline%Frc)
  if(allocated(rspline%w ))     deallocate(rspline%w)
  if(allocated(rspline%cc ))    deallocate(rspline%cc)
  if(allocated(rspline%tt ))    deallocate(rspline%tt)
  if(allocated(rspline%der))    deallocate(rspline%der)
  if(allocated(rspline%smooth)) deallocate(rspline%smooth)

  kkk           = kkk2
  rspline%N     = kk
  rspline%order = kkk
  rspline%nest  = max(kk+kkk+1,2*kkk+2)
  if(kkk>5) then
    if(testing) write(*,*) 'error spline of too high order'
    rspline%ierr=1
    return
  endif
  if(kk<=kkk)    kkk=kk-1

  allocate(rspline%om(kk))
  allocate(rspline%Frc(kk))
  allocate(rspline%cc(rspline%nest))
  allocate(rspline%tt(rspline%nest))
  allocate(rspline%der(kkk,kk))
  allocate(rspline%smooth(kk))
  allocate(rspline%w(kk))


  21 continue

  rspline%om=0
  rspline%Frc=0
  rspline%cc=0
  rspline%tt=0
  rspline%der=0
  rspline%smooth=0
  rspline%w=0


  rspline%om=omi
  rspline%Frc=Frc
  rspline%w=(/( 1.d0,i=1,kk )/)

  if(present(errors))then
   rspline%w=abs(errors)
   rspline%w=(10.d0*rspline%w/maxval(rspline%w))**2.d0
   where(rspline%w<1.d-2) rspline%w=1.d-2
   rspline%w=1.d0/rspline%w
  endif

  tryit=0
  20 continue

  if(present(s2)) then
    call define_spline(s2=s2)
  else
    call define_spline
  endif

  if(ier==10) then
   if(tryit==0)then
    if(testing) write(*,*) 'error init spline, try to sort it out....'
    call sortitout
    tryit=1
    goto 20
   else
    rspline%ierr=1
    write(*,*) 'error init spline, nope, nothing possible....'
    if(strongstop) stop
   endif
  endif


 return
  
 contains
 
   !-----------------------!
   !-----------------------!
   !-----------------------!

 subroutine define_spline(s2) 
  implicit none
  real(8),optional  :: s2
  real(4)           :: wrk(rspline%N*(rspline%order+3)+rspline%nest*(7+4*rspline%order)+10),tt(rspline%nest),cc(rspline%nest)
  integer           :: iwrk(rspline%nest)

   rspline%ierr=0
   lwrk=size(wrk)
   if(present(s2))then
    s=real(s2)
   else
    s=0.
   endif
   if(rspline%order>=rspline%N) then
    write(*,*) 'spline not enough points'
    write(*,*) 'number of points : ' , rspline%N
    write(*,*) ' order of splines: ', rspline%order
    rspline%ierr=1
    return
   endif
   if(rspline%nest<=2*rspline%order+2) then
     write(*,*) 'spline nest too small'
     rspline%ierr=1
     return
   endif
   xl=rspline%om(1)-1.e-7
   xr=rspline%om(rspline%N)+1.e-7
   call curfit(0,rspline%N,rspline%om,rspline%Frc,rspline%w,xl,xr, &
            & rspline%order,s,rspline%nest,rspline%nn,rspline%tt,rspline%cc,fp,wrk,lwrk,iwrk,ier)

 end subroutine
 
   !-----------------------!
   !-----------------------!
   !-----------------------!

 subroutine sortitout
 implicit none
 integer :: order(rspline%N)
 real(8) :: iiFrc(rspline%N),iiom(rspline%N)
  iiFrc=rspline%Frc
  iiom=rspline%om
  call qsort_array(iiom,order)
  call qsort_adj_array(iiFrc,order)
  k=rspline%N
  call group_data_rrr(k,iiom,iiFrc,iiom,k)
  call kill_spline(rspline)
  rspline%ierr=0
  kk=k
  if(allocated(rspline%om))     deallocate(rspline%om)
  if(allocated(rspline%Frc))    deallocate(rspline%Frc)
  if(allocated(rspline%w ))     deallocate(rspline%w)
  if(allocated(rspline%cc ))    deallocate(rspline%cc)
  if(allocated(rspline%tt ))    deallocate(rspline%tt)
  if(allocated(rspline%der))    deallocate(rspline%der)
  if(allocated(rspline%smooth)) deallocate(rspline%smooth)
  kkk           = kkk2
  rspline%N     = kk
  rspline%order = kkk
  rspline%nest  = max(kk+kkk+1,2*kkk+2)
  if(kk<=kkk)    kkk=kk-1
  allocate(rspline%om(kk))
  allocate(rspline%Frc(kk))
  allocate(rspline%cc(rspline%nest))
  allocate(rspline%tt(rspline%nest))
  allocate(rspline%der(kkk,kk))
  allocate(rspline%smooth(kk))
  allocate(rspline%w(kk))
  rspline%om=0
  rspline%Frc=0
  rspline%cc=0
  rspline%tt=0
  rspline%der=0
  rspline%smooth=0
  rspline%w=0
  rspline%om=iiom
  rspline%Frc=iiFrc
  rspline%w=(/( 1.d0,i=1,kk )/)
 end subroutine

   !-----------------------!
   !-----------------------!
   !-----------------------!

 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine kill_spline(rspline)
  implicit none
  type(spline)    :: rspline
  if(allocated(rspline%om))     deallocate(rspline%om)
  if(allocated(rspline%Frc))    deallocate(rspline%Frc)
  if(allocated(rspline%w ))     deallocate(rspline%w)
  if(allocated(rspline%cc ))    deallocate(rspline%cc)
  if(allocated(rspline%tt ))    deallocate(rspline%tt)
  if(allocated(rspline%der))    deallocate(rspline%der)
  if(allocated(rspline%smooth)) deallocate(rspline%smooth)
  rspline%initialized=0
end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine evaluate_splineb(rspline,xx,yy)
 implicit none
 type(spline)      :: rspline
 integer           :: kk,nn,iii,mm,siz
 real(8)           :: xx(:),yy(:)
 real(4)           :: yy2(size(yy)),xxx(size(xx))
 integer           :: i,ier

   if(rspline%ierr>0) return
   mm=size(xx(:))
   kk=rspline%order
   xxx(1:size(xx))=xx
   call splev(rspline%tt,rspline%nn,rspline%cc,rspline%order,xxx,yy2,mm,ier)
   yy=yy2(1:size(yy))

 return
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine evaluate_splinebr(rspline,xx,yy)
 implicit none
 type(spline)      :: rspline
 integer           :: kk,nn,iii,mm,siz
 real(4)           :: xx(:),yy(:)
 integer           :: i,ier

   if(rspline%ierr>0) return
   mm=size(xx(:))
   kk=rspline%order
   call splev(rspline%tt,rspline%nn,rspline%cc,rspline%order,xx,yy,mm,ier)

 return
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine evaluate_splinea(rspline)
 implicit none
 type(spline)      :: rspline
 integer           :: kk,nn,iii
 integer           :: i,ier

   !---------------------------------------------------------------------!
   !  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
   !  a spline s(x) of degree k, given in its b-spline representation.
   !  calling sequence:
   !     call splev(t,n,c,k,x,y,m,ier)
   !  input parameters:
   !    t    : array,length n, which contains the position of the knots.
   !    n    : integer, giving the total number of knots of s(x).
   !    c    : array,length n, which contains the b-spline coefficients.
   !    k    : integer, giving the degree of s(x).
   !    x    : array,length m, which contains the points where s(x) must
   !           be evaluated.
   !    m    : integer, giving the number of points where s(x) must be
   !           evaluated.
   !  output parameter:
   !    y    : array,length m, giving the value of s(x) at the different
   !           points.
   !    ier  : error flag
   !      ier = 0 : normal return
   !      ier =10 : invalid input data (see restrictions)
   !---------------------------------------------------------------------!

  if(rspline%ierr>0) return
  call splev(rspline%tt,rspline%nn,rspline%cc,rspline%order, &
              & rspline%om,rspline%smooth,rspline%N,ier)

 return
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
! 
!  subroutine derivate_splinea(rspline)
!  implicit none
!  type(spline)      :: rspline
!  integer           :: kk,nn,iii,mm
!  integer           :: i,ier
!  real(4)           :: wrk(size(rspline%tt))
! 
!    !------------------------------------------------------------------------!
!    !  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
!    !  the derivative of order nu of a spline s(x) of degree k,given in
!    !  its b-spline representation.
!    !  calling sequence:
!    !     call splder(t,n,c,k,nu,x,y,m,wrk,ier)
!    !  input parameters:
!    !    t    : array,length n, which contains the position of the knots.
!    !    n    : integer, giving the total number of knots of s(x).
!    !    c    : array,length n, which contains the b-spline coefficients.
!    !    k    : integer, giving the degree of s(x).
!    !    nu   : integer, specifying the order of the derivative. 0<=nu<=k
!    !    x    : array,length m, which contains the points where the deriv-
!    !           ative of s(x) must be evaluated.
!    !    m    : integer, giving the number of points where the derivative
!    !           of s(x) must be evaluated
!    !    wrk  : real array of dimension n. used as working space.
!    !  output parameters:
!    !    y    : array,length m, giving the value of the derivative of s(x)
!    !           at the different points.
!    !    ier  : error flag
!    !      ier = 0 : normal return
!    !      ier =10 : invalid input data (see restrictions)
!    !------------------------------------------------------------------------!
! 
!    ier=0;wrk=0
!    if(rspline%ierr>0) return
! 
!     do iii=0,rspline%order-1
!      call splder(rspline%tt,                           &
!               & rspline%nn,                            &
!               & rspline%cc,                            &
!               & rspline%order,                         &
!               & iii,                                   &
!               & rspline%om,                            &
!               & rspline%der(iii+1,1:size(rspline%om)), &
!               & size(rspline%om),wrk,ier)
!     enddo
! 
!  return
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  subroutine derivate_splineb(rspline,xx,yy,mu)
!  implicit none
!  type(spline)      :: rspline
!  integer           :: kk,nn,iii,mm
!  integer           :: mu
!  real(8)           :: xx(:),yy(:,:)
!  real(4)           :: yy2(size(yy(:,1)),size(yy(1,:)))
!  integer           :: i,ier
!  real(4)           :: wrk(size(rspline%tt))
! 
!    if(rspline%ierr>0) return
!    mm=size(xx(:))
!     do iii=0,mu-1
!        call splder(rspline%tt,rspline%nn,rspline%cc,rspline%order,iii,real(xx),yy2(iii+1,:),mm,wrk,ier)
!     enddo
!     yy=yy2
! 
!  return
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  subroutine derivate_splined(rspline,xx,yy,mu)
!  implicit none
!  type(spline)      :: rspline
!  integer           :: kk,nn,mm
!  integer           :: mu
!  real(8)           :: xx(:),yy(:)
!  real(4)           :: yy2(size(yy(:)))
!  integer           :: i,ier
!  real(4)           :: wrk(size(rspline%tt))
!    if(rspline%ierr>0) return
!    mm=size(xx)
!    call splder(rspline%tt,rspline%nn,rspline%cc,rspline%order,mu,real(xx),yy2,mm,wrk,ier)
!    yy=yy2
!  return
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! 
!  subroutine derivate_splinebr(rspline,xx,yy,mu)
!  implicit none
!  type(spline)      :: rspline
!  integer           :: kk,nn,iii,mm
!  integer           :: mu
!  real(4)           :: xx(:),yy(:,:)
!  integer           :: i,ier
!  real(4)           :: wrk(size(rspline%tt))
!    if(rspline%ierr>0) return
!    mm=size(xx(:))
!     do iii=0,mu-1
!        call splder(rspline%tt,rspline%nn,rspline%cc,rspline%order,iii,xx,yy(iii+1,:),mm,wrk,ier)
!     enddo
!  end subroutine
! 
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!
! !*****************************************************!

end module
