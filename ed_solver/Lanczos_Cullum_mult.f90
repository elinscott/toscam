MODULE Lanczos_Cullum_mult

  use correlations
  use namelistmod
  
  IMPLICIT NONE


  logical,private :: verbose=.false.

CONTAINS



  !-----HLEMULT----HERMITIAN MATRICES-------------------------------------
  !     CONTAINS SUBROUTINE LANCZS

  !-----LANCZS-COMPUTE THE LANCZOS TRIDIAGONAL MATRICES-------------------
  !     GRAM-SCHMIDT ORTHOGONALIZATION 
  !     REQUIRES EXTRA VECTOR VS IN LANCZS.  

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE LANCZS_C(MATVEC,V1,V2,VS,ALPHA,BETA,GR,GC,G,KMAX,MOLD1,N,IIX)
    COMPLEX(8) :: V1(:),V2(:),ZEROC,TEMP,VS(:)
    REAL(8)    :: ALPHA(:), BETA(:), BATA, SUM, ONE, ZERO
    REAL(8)    :: GR(:),GC(:)
    REAL(8)    :: G(:)
    INTEGER    :: IIX(:), IIL(SIZE(IIX)),KMAX,MOLD1,N,I,J,IN

    INTERFACE 
      SUBROUTINE MATVEC(V1,V2,COEF)
        COMPLEX(8)       :: V1(:)
        COMPLEX(8)       :: V2(:)
        REAL(8)          :: COEF
      END SUBROUTINE 
    END INTERFACE

    if(SIZE(VS)/=SIZE(V1).or.SIZE(V1)/=SIZE(V2)) stop 'error LANCZS_C shape of arrays do not match'

    if(N<2) then
      write(*,*) 'danger, diagonalize hilbert space with N<2'
      stop
    endif

    ZERO  = 0.D0
    ONE   = 1.D0
    ZEROC = CMPLX(ZERO,ZERO,8)

          IF(MOLD1.GT.1) GO TO 50

    !     ALPHA/BETA GENERATION STARTS AT I = 1
    !     MOLD1 = 1 SET V1 = 0. AND V2 = RANDOM UNIT VECTOR

    IIL=IIX

    !-----------------------------------------------------------------------
    CALL GENRAN(IIL,G,N)
    !-----------------------------------------------------------------------
   
    if(verbose) write(*,*) '========> start LEVAL with N = : ', N
 
    DO 10 I = 1,N
      10 GR(I) = G(I)

      !-----------------------------------------------------------------------
      CALL GENRAN(IIL,G,N)
      !-----------------------------------------------------------------------

      DO 20 I = 1,N
      20 GC(I) = G(I)
        !
      DO 30 I = 1,N
      30 V2(I) = CMPLX(GR(I),GC(I),8)
 
          !-----------------------------------------------------------------------
          SUM = DBLE(MPI_DOT_PRODUCT(V2(1:N),V2(1:N)))
          !-----------------------------------------------------------------------
 
          if(sum<epsilonr) then
            write(*,*) 'LANCZOS COMPLEX NEGATIVE SQR OR ZERO DIVISION (a): ', sum
            write(*,*) ' N : ', N
            call write_first_el(10,V1)
            call write_first_el(10,V2)
            stop 'error termine'
          endif

          SUM = ONE/DSQRT(SUM)

            DO 40 I = 1,N
               V1(I) = ZEROC
            40 V2(I) = V2(I)*SUM
  
           BETA(1) = ZERO

          ! ALPHA BETA GENERATION LOOP
            50 CONTINUE

            !########################################################! 
            DO 80 I=MOLD1,KMAX
              SUM = ZERO

              !-----------------------------------------------------------------------
              !     MATVEC(V2,VS,SUM) CALCULATES  VS = A*V2 - SUM*VS
              CALL MATVEC(V2(1:N),VS(1:N),SUM)
              if(verbose) write(*,*) 'leval, build VS : ', VS(1:10)
              SUM = DBLE(MPI_DOT_PRODUCT(V2(1:N),VS(1:N)))
              !-----------------------------------------------------------------------
 
              ALPHA(I) = SUM
              BATA     = BETA(I)

              DO 60 J=1,N
                60 V1(J) = (VS(J)-BATA*V1(J)) - SUM*V2(J)

                !-----------------------------------------------------------------------
                SUM = DBLE(MPI_DOT_PRODUCT(V1(1:N),V1(1:N)))
                !-----------------------------------------------------------------------

               IN = I+1

               if(sum<epsilonr) then
                 write(*,*) 'LANCZOS COMPLEX NEGATIVE SQR OR ZERO DIVISION (b): ', sum,epsilonr
                 write(*,*) ' N : ',N
                 call write_first_el(10,V1)
                 call write_first_el(10,V2)
                 call write_first_el(10,VS) 
                 write(*,*) norme(V1),norme(V2),norme(VS)
                 write(*,*) 'bata :',bata
                 write(*,*) 'dot : ',MPI_DOT_PRODUCT(V1(1:N),V1(1:N))
                 write(*,*) 'I MOLD1 : ', I, MOLD1
                 stop 'error termine'
               endif

               BETA(IN) = DSQRT(SUM)
               SUM = ONE/BETA(IN)

                DO 70 J=1,N
                   TEMP = SUM*V1(J)
                   V1(J) = V2(J)
                70 V2(J) = TEMP
               80 CONTINUE

            !########################################################!



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

 
  SUBROUTINE LANCZS_R(MATVEC,ALPHA,BETA,V1,V2,G,KMAX,MOLD1,N,IIX)
    REAL(8) ::  ALPHA(:),BETA(:),V1(:),V2(:),SUM,TEMP,ONE,ZERO
    REAL(8) ::  G(:)
    INTEGER ::  IIX(:),IIL(SIZE(IIX)),KMAX,MOLD1,N,I,J,IN

    INTERFACE
      SUBROUTINE MATVEC(V1,V2,COEF)
        REAL(8) :: V1(:)
        REAL(8) :: V2(:)
        REAL(8) :: COEF
      END SUBROUTINE MATVEC
    END INTERFACE

    !-----------------------------------------------------------------------
    ZERO = 0.D0
    ONE  = 1.D0
    !
    IF(MOLD1.GT.1)GO TO 30
    !
    !    ALPHA/BETA GENERATION STARTS AT I = 1
    !    MOLD1 = 1 SET V1 = 0. AND V2 = RANDOM UNIT VECTOR
    BETA(1) = ZERO
    IIL=IIX
    !
    !-----------------------------------------------------------------------
    CALL GENRAN(IIL,G,N)
    !-----------------------------------------------------------------------
    !
    DO 10 I = 1,N
      10 V2(I) = G(I)

      !-----------------------------------------------------------------------
      SUM = MPI_DOT_PRODUCT(V2(1:N),V2(1:N))
      !-----------------------------------------------------------------------

      SUM = ONE/DSQRT(SUM)
      DO 20 I = 1,N
 V1(I) = ZERO
 20 V2(I) = V2(I)*SUM
 !
 !    ALPHA BETA GENERATION LOOP
 30 CONTINUE
 !
 DO 60 I=MOLD1,KMAX
   SUM = BETA(I)
   !    MATVEC(V2,V1,SUM) CALCULATES  V1 = A*V2 - SUM*V1

   !-----------------------------------------------------------------------
   CALL MATVEC(V2(1:N),V1(1:N),SUM)
   !-----------------------------------------------------------------------

   !-----------------------------------------------------------------------
   SUM = MPI_DOT_PRODUCT(V1(1:N),V2(1:N))
   !-----------------------------------------------------------------------
 
   ALPHA(I) = SUM
   DO 40 J=1,N
   40 V1(J) = V1(J)-SUM*V2(J)

     !-----------------------------------------------------------------------
     SUM = MPI_DOT_PRODUCT(V1(1:N),V1(1:N))
     !-----------------------------------------------------------------------

     IN = I+1
     BETA(IN) = DSQRT(SUM)
     SUM = ONE/BETA(IN)
     DO 50 J=1,N
       TEMP = SUM*V1(J)
       V1(J) = V2(J)
       50 V2(J) = TEMP
       60 CONTINUE
       !
       !     END ALPHA, BETA GENERATION LOOP
       !
   RETURN
   END subroutine



   END MODULE 
