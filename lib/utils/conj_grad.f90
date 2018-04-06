MODULE conj_grad

  USE common_def, only: c2s, dump_message, i2c
  use genvar,     only: DBL, log_unit

  implicit none

  LOGICAL,              private :: verbose=.false.

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

  SUBROUTINE conjugate_gradient(funct,x,n,fmin,step,fmax,mode,ncall,iprint)
    USE common_def
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ CONJUGUATE GRADIENT MINIMIZATION OF FUNTION 'funct' $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    INTERFACE
      SUBROUTINE funct(val,n,x)
        USE common_def
        use genvar, only: DBL
        REAL(DBL), INTENT(OUT) :: val
        INTEGER,   INTENT(IN)  :: n
        REAL(DBL), INTENT(IN)  :: x(n)
      END SUBROUTINE
    END INTERFACE
    REAL(DBL), INTENT(OUT)   :: fmin
    REAL(DBL), INTENT(INOUT) :: x(n)      ! initial vector to evaluate funct.
    INTEGER,   INTENT(IN)    :: ncall     ! max. # of calls to funct.
    INTEGER,   INTENT(IN)    :: n         ! # of parameters
    REAL(DBL), INTENT(IN)    :: step,fmax ! step/tolerance
    INTEGER,   INTENT(IN)    :: iprint,mode
    REAL(DBL)                :: g(n),xm(n)
    REAL(DBL)                :: w(n*4)
    REAL(DBL)                :: h(n*(n+1)/2)
    INTEGER                  :: iexit
    INTEGER                  :: n4,np1,nm1,nnp1s2
    INTEGER                  :: is,iu,iv,ib,idiff,i_,j_,ij,i1
    INTEGER                  :: jk,ik,k,itn,icall,link,int_
    REAL(DBL)                :: f_,df,dfn,z,zz,dmin,gs0,afmax
    REAL(DBL)                :: alpha,ff,tot,f1,f2,gys,dgs,sig
    f_=0.
    n4  = n * 4
    np1 = n + 1
    nm1 = n - 1
    nnp1s2 = (n*np1)/2
    dfn = - 0.5_DBL
    xm = ABS(x) + 1.D-10
    is = n
    iu = n
    iv = n + n
    ib = iv + n
    idiff = 1
    iexit = 0
    SELECT CASE(mode)
      CASE(1)
        !
        h = 0.0_DBL
        DO i_=1,n
          h(nnp1s2-i_*(i_+1)/2+1) = 1.0_DBL
        END DO
      CASE(2)
        !
        ij = 1
        DO i_ = 2,n
          z = h(ij)
          IF(z<=0.0_DBL)THEN
            fmin = f_
            RETURN
          END IF
          ij = ij + 1
          i1 = ij
          DO j_ = i_,n
            zz = h(ij)
            h(ij) = h(ij) / z
            jk = ij
            ik = i1
            DO k = i_,j_
              jk = jk + np1 - k
              h(jk) = h(jk) - h(ik) * zz
              ik = ik + 1
            END DO
            ij = ij + 1
          END DO
        END DO
        IF(h(ij)<=0.0_DBL)THEN
          fmin = f_
          RETURN
        END IF
    END SELECT
    ij = np1
    dmin = h(1)
    DO i_=2,n
      IF(h(ij)<dmin) dmin = h(ij)
      ij = ij + np1 - i_
    END DO
    IF(dmin<=0.0_DBL)THEN
      fmin = f_
      RETURN
    END IF
    z = f_
    itn = 0
    CALL funct(f_,n,x)
    icall = 1
    df = dfn
    IF(dfn==0.0_DBL) df = f_ - z
    IF(dfn<0.0_DBL)  df = ABS(df*f_)
    IF(df<=0.0_DBL)  df = 1.0_DBL
    17 CONTINUE
    w(1:n) = x(1:n)
    link = 1
    IF(idiff-1) 100, 100, 110 ! COMPUTE GRADIENT/SECOND ORDER DERIVATIVE
    18 CONTINUE
    IF(icall>=ncall) GOTO 90
    20 CONTINUE
    IF(iprint/=0.AND.MOD(itn,iprint)==0)THEN
      CALL dump_message(TEXT="#################")
      CALL dump_message(TEXT="# iteration     = "//c2s(i2c(itn+1))//" (number of calls to fctn = "//c2s(i2c(icall))//")")
      WRITE(log_unit,'(a,f0.12)')                        "#          f(x) = ",f_
      IF(iprint>0)THEN
        WRITE(log_unit,'(a,'//c2s(i2c(n))//'(f9.6,2x))') "# at point   x  = ",(x(i_),i_=1,n)
        WRITE(log_unit,'(a,'//c2s(i2c(n))//'(f9.6,2x))') "# gradient g(x) = ",(g(i_),i_=1,n)
      END IF
    END IF
    itn = itn + 1
    w(1) = -g(1)
    DO i_ = 2,n
      ij = i_
      i1 = i_ - 1
      z = -g(i_)
      DO j_ = 1, i1
        z = z - h(ij) * w(j_)
        ij = ij + n - j_
      END DO
      w(i_) = z
    END DO
    w(is+n) = w(n) / h(nnp1s2)
    ij = nnp1s2
    DO i_ = 1,nm1
      ij = ij - 1
      z = 0.0_DBL
      DO j_ = 1,i_
        z = z + h(ij) * w(is+np1-j_)
        ij = ij - 1
      END DO
      w(is+n-i_) = w(n-i_) / h(ij) - z
    END DO
    z   = 0.0_DBL
    gs0 = SUM(g(1:n)*w(is+1:is+n))
    DO i_ = 1,n
      IF(z*xm(i_)<ABS(w(is+i_))) z = ABS(w(is+i_)) / xm(i_)
    END DO
    afmax = fmax / z
    iexit = 2
    IF(gs0>=0.0_DBL) GOTO 92
    alpha = - 2.0_DBL * df / gs0
    IF(alpha>1.0_DBL) alpha = 1.0_DBL
    ff = f_
    tot = 0.0_DBL
    int_ = 0
    iexit = 1
    30  CONTINUE
    IF(icall>=ncall) GOTO 90
    w(1:n) = x(1:n) + alpha * w(is+1:is+n)
    CALL funct(f1,n,w)
    icall = icall + 1
    IF(f1>=f_) GOTO 40
    f2 = f_
    tot = tot + alpha
    32  CONTINUE
    x(1:n) = w(1:n)
    f_ = f1
    IF(int_-1) 35, 49, 50
    35  CONTINUE
    IF(icall>=ncall) GOTO 90
    w(1:n) = x(1:n) + alpha * w(is+1:is+n)
    CALL funct(f1,n,w)
    icall = icall + 1
    IF(f1>=f_) GOTO 50
    IF((f1+f2>=f_+f_).AND.(7.0d0*f1+5.0d0*f2>12.0d0*f_)) int_ = 2
    tot = tot + alpha
    alpha = 2.0_DBL * alpha
    GOTO 32
    40  CONTINUE
    IF(alpha<afmax) GOTO 92
    IF(icall>=ncall) GOTO 90
    alpha = 0.5_DBL * alpha
    w(1:n) = x(1:n) + alpha * w(is+1:is+n)
    CALL funct(f2,n,w)
    icall = icall + 1
    IF(f2>=f_) GOTO 45
    tot = tot + alpha
    f_ = f2
    x(1:n) = w(1:n)
    GOTO 49
    45  CONTINUE
    z = 0.1d0
    IF(f1+f_>f2+f2) z = 1.0_DBL + 0.5_DBL * (f_ - f1) / (f_ + f1 - f2 - f2)
    IF(z<0.1d0) z = 0.1d0
    alpha = z * alpha
    int_ = 1
    GOTO 30
    49  CONTINUE
    IF(tot<afmax) GOTO 92
    50  CONTINUE
    alpha = tot
    w(1:n)       = x(1:n)
    w(ib+1:ib+n) = g(1:n)
    link = 2
    IF(idiff-1) 100, 100, 110 ! COMPUTE GRADIENT/SECOND ORDER DERIVATIVE
    54  CONTINUE
    IF(icall>=ncall) GOTO 90
    w(1:n) = w(ib+1:ib+n)
    gys = SUM(g(1:n)*w(is+1:is+n))
    df = ff - f_
    dgs = gys - gs0
    IF(dgs<=0.0_DBL) GOTO 20
    link = 1
    IF(dgs+alpha*gs0>0.0_DBL)THEN
      zz = alpha / (dgs - alpha * gs0)
      z = dgs * zz - 1.0_DBL
      w(iu+1:iu+n) = z * w(1:n) + g(1:n)
      sig = 1.0_DBL / (zz * dgs * dgs)
    ELSE
      w(iu+1:iu+n) = g(1:n) - w(1:n)
      sig = 1.0_DBL / (alpha * dgs)
    END IF
    GOTO 70
    60  CONTINUE
    link = 2
    w(iu+1:iu+n) = w(1:n)
    IF(dgs+alpha*gs0>0.0_DBL)THEN
      sig = - zz
    ELSE
      sig = 1.0_DBL / gs0
    END IF
    70  CONTINUE
    w(iv+1) = w(iu+1)
    DO i_ = 2,n
      ij = i_
      i1 = i_ - 1
      z = w(iu+i_)
      DO j_ = 1,i1
        z = z - h(ij) * w(iv+j_)
        ij = ij + n - j_
      END DO
      w(iv+i_) = z
    END DO
    ij = 1
    DO i_ = 1,n
      z = h(ij) + sig * w(iv+i_) * w(iv+i_)
      IF(z<=0.0_DBL)  z = dmin
      IF(z<dmin) dmin = z
      h(ij) = z
      w(ib+i_) = w(iv+i_) * sig / z
      sig = sig - w(ib+i_) * w(ib+i_) * z
      ij = ij + np1 - i_
    END DO
    ij = 1
    DO i_ = 1,nm1
      ij = ij + 1
      i1 = i_ + 1
      DO j_ = i1,n
        w(iu+j_) = w(iu+j_) - h(ij) * w(iv+i_)
        h(ij) = h(ij) + w(ib+i_) * w(iu+j_)
        ij = ij + 1
      END DO
    END DO
    GOTO (60, 20), link
    90  CONTINUE
    iexit = 3
    GOTO 94
    92  CONTINUE
    IF(idiff==2) GOTO 94
    idiff = 2
    GOTO 17
    94  CONTINUE
    fmin = f_
    IF(iprint/=0)THEN
      CALL dump_message(TEXT="#################")
      CALL dump_message(TEXT="#################")
      CALL dump_message(TEXT="# STOPPED iexit ="//c2s(i2c(iexit)))
      CALL dump_message(TEXT="# iteration     = "//c2s(i2c(itn+1))//" (number of calls to fctn = "//c2s(i2c(icall))//")")
      WRITE(log_unit,'(a,f0.12)')                      "#          f(x) = ",f_
      WRITE(log_unit,'(a,'//c2s(i2c(n))//'(f9.6,2x))') "# at point   x  = ",(x(i_),i_=1,n)
      WRITE(log_unit,'(a,'//c2s(i2c(n))//'(f9.6,2x))') "# gradient g(x) = ",(g(i_),i_=1,n)
      CALL dump_message(TEXT="#################")
      CALL dump_message(TEXT="#################")
    END IF
    RETURN
    100 CONTINUE
    ! COMPUTE GRADIENT g OF FUNCTION funct
    CALL gradient (g(1:n),funct,w(1:n),step*xm(1:n),f_,icall)
    GOTO (18, 54), link
    110 CONTINUE
    ! COMPUTE SECOND ORDER DERIVATIVE g OF FUNCTION funct
    CALL gradient2(g(1:n),funct,w(1:n),step*xm(1:n),icall)
    GOTO (18, 54), link
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

  SUBROUTINE gradient(grad,funct,x,step,fx,ncall)
    USE common_def
    INTERFACE
      SUBROUTINE funct(val,n,xx)
        USE common_def
        use genvar, only: DBL
        REAL(DBL), INTENT(OUT) :: val
        INTEGER,   INTENT(IN)  :: n
        REAL(DBL), INTENT(IN)  :: xx(n)
      END SUBROUTINE
    END INTERFACE
    REAL(DBL), INTENT(INOUT) :: grad(:)
    REAL(DBL), INTENT(IN)    :: x(:),step(:)
    REAL(DBL), INTENT(IN),    OPTIONAL :: fx    ! value of function at x
    INTEGER,   INTENT(INOUT), OPTIONAL :: ncall ! number of calls to funct
    REAL(DBL) :: y(SIZE(x)),fy,myfx 
    INTEGER   :: icomp,ncomp
    !
    ! COMPUTE grad(funct)[i](x) = ( f(x+step(i)) - f(x) ) / step(i) 
    ! GRADIENT OF funct at x ALONG 'i'th DIRECTION
    !
    IF(SIZE(grad)/=SIZE(x).OR.SIZE(step)/=SIZE(x)) STOP "ERROR IN gradient: INCONSISTENT SIZES!"
    ncomp = SIZE(x)
    DO icomp=1,ncomp
      y = x
      y(icomp) = y(icomp) + step(icomp)
      CALL funct(fy,ncomp,y)
      IF(PRESENT(fx))THEN ! f(x) was provided
        grad(icomp) = (fy-fx) / step(icomp)
      ELSE ! f(x) was NOT provided => compute it
        CALL funct(myfx,ncomp,x)
        grad(icomp) = (fy-myfx) / step(icomp)
      END IF
    END DO
    IF(PRESENT(ncall))THEN
      IF(PRESENT(fx))THEN
        ncall = ncall + ncomp
      ELSE
        ncall = ncall + ncomp*2
      END IF
    END IF
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

 
  SUBROUTINE gradient2(grad2,funct,x,step,ncall)
    USE common_def
    INTERFACE
      SUBROUTINE funct(val,n,xx)
        USE common_def
        use genvar, only: DBL
        REAL(DBL), INTENT(OUT) :: val
        INTEGER,   INTENT(IN)  :: n
        REAL(DBL), INTENT(IN)  :: xx(n)
      END SUBROUTINE
    END INTERFACE
    REAL(DBL), INTENT(INOUT) :: grad2(:)
    REAL(DBL), INTENT(IN)    :: x(:),step(:)
    INTEGER,   INTENT(INOUT), OPTIONAL :: ncall ! number of calls to funct
    REAL(DBL) :: y(SIZE(x)),fyp,fym 
    INTEGER   :: icomp,ncomp
    ! COMPUTE SECOND ORDER DERIVATIVE g(i) ALONG i'th DIRECTION
    IF(SIZE(grad2)/=SIZE(x).OR.SIZE(step)/=SIZE(x)) STOP "ERROR IN gradient: INCONSISTENT SIZES!"
    ncomp = SIZE(x)
    DO icomp=1,ncomp
      y = x
      y(icomp) = y(icomp) + step(icomp)
      CALL funct(fyp,ncomp,y)
      y = x
      y(icomp) = y(icomp) - step(icomp)
      CALL funct(fym,ncomp,y)
      grad2(icomp) = (fyp-fym) / (2.0_DBL*step(icomp))
    END DO
    IF(PRESENT(ncall)) ncall = ncall + ncomp*2
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

  
  ! PROGRAM: minimize
  ! TYPE   : subroutine
  ! PURPOSE: conjugent gradient search
  ! I/O    : 
  ! VERSION: 30-Sep-95
  ! COMMENT: This is a most reliable conjugent gradient routine! It has
  !          served us well for many years, and is capable to cope with
  !          a very large number of variables. Unfortunately, we don't
  !          know who wrote this routine (original name: 'va10a'), and 
  !          we find it very obscure.
  !          Don't worry, it works just fine.
 


  subroutine minimize(funct,n,x,f,g,h,w,dfn,xm,hh,eps,mode,maxfn,iprint,iexit)

    INTERFACE
      SUBROUTINE funct(val,n,x)
        USE common_def
        use genvar, only: DBL
        REAL(DBL), INTENT(OUT) :: val
        INTEGER,   INTENT(IN)  :: n
        REAL(DBL), INTENT(IN)  :: x(n)
      END SUBROUTINE
    END INTERFACE

    INTEGER                        :: n,mode,maxfn,iprint,iexit
    DOUBLE PRECISION, dimension(*) :: x, g, h, w, xm
    DOUBLE PRECISION               :: f,dfn,eps,hh
    INTEGER                        :: i,j,np,n1,nn,is,iu,iv,ib,idiff,ij,i1,jk,ik,k,itn,ifn,link,int
    DOUBLE PRECISION               :: f1,f2,z,zz,dmin,df,gs0,aeps,alpha,ff,tot,gys,dgs,sig

    if (iprint .ne. 0 .and. verbose ) write (6,1000)
    1000 format (' entry into minimize')
    np = n + 1
    n1 = n - 1
    nn=(n*np)/2
    is = n
    iu = n
    iv = n + n
      ib = iv + n
      idiff = 1
      iexit = 0
      if (mode .eq. 3) go to 15
      if (mode .eq. 2) go to 10
      ij = nn + 1
      do 5 i = 1, n
      do 6 j = 1, i
      ij = ij - 1
   6  h(ij) = 0.0_DBL
   5  h(ij) = 1.0_DBL
      go to 15
  10  continue
      ij = 1
      do 11 i = 2, n
      z = h(ij)
      if (z .le. 0.0_DBL) return
      ij = ij + 1
      i1 = ij
      do 11 j = i, n
      zz = h(ij)
      h(ij) = h(ij) / z
      jk = ij
      ik = i1
      do 12 k = i, j
      jk = jk + np - k
      h(jk) = h(jk) - h(ik) * zz
      ik = ik + 1
  12  continue
      ij = ij + 1
  11  continue
      if (h(ij) .le. 0.0_DBL) return
  15  continue
      ij = np
      dmin = h(1)
      do 16 i = 2, n
      if (h(ij) .ge. dmin) go to 16
      dmin = h(ij)
  16  ij = ij + np - i
      if (dmin .le. 0.0_DBL) return
      z = f
      itn = 0
      call funct(f,n,x)
      ifn = 1
      df = dfn
      if (dfn .eq. 0.0_DBL) df = f - z
      if (dfn .lt. 0.0_DBL) df = abs (df * f)
      if (df .le. 0.0_DBL) df = 1.0_DBL
  17  continue
      do 19 i = 1, n
      w(i) = x(i)
  19  continue
      link = 1
      if (idiff - 1) 100, 100, 110
  18  continue
      if (ifn .ge. maxfn) go to 90
  20  continue
      if (iprint .eq. 0) go to 21
      if (mod (itn, iprint) .ne. 0) go to 21
      if(verbose)  write (6,1001) itn, ifn
1001  format (1x,'itn = ',i5,' ifn = ',i5)
      if(verbose) write (6,1002) f
1002  format (1x,'f = ',e15.7)
      if (iprint .lt. 0) go to 21
      if(verbose) write (6,1003) (x(i), i = 1, n)
!***
!***
1003  format (1x,'x = ',4e15.7 / (5x, 4e15.7))
      if(verbose) write (6,1004) (g(i), i = 1, n)
1004  format (1x,'g = ',4e15.7 / (5x, 4e15.7))
!**
!***
  21  continue
      itn = itn + 1
      w(1) = -g(1)
      do 22 i = 2, n
      ij = i
      i1 = i - 1
      z = -g(i)
      do 23 j = 1, i1
      z = z - h(ij) * w(j)
      ij = ij + n - j
  23  continue
  22  w(i) = z
      w(is+n) = w(n) / h(nn)
      ij = nn
      do 25 i = 1, n1
      ij = ij - 1
      z = 0.0_DBL
      do 26 j = 1, i
      z = z + h(ij) * w(is+np-j)
      ij = ij - 1
  26  continue
  25  w(is+n-i) = w(n-i) / h(ij) - z
      z   = 0.0_DBL
      gs0 = 0.0_DBL
      do 29 i = 1, n
      if (z * xm(i) .ge. abs (w(is+i))) go to 28
      z = abs (w(is+i)) / xm(i)
  28  gs0 = gs0 + g(i) * w(is+i)
  29  continue
      aeps = eps / z
      iexit = 2
      if (gs0 .ge. 0.0_DBL) go to 92
      alpha = -2.0_DBL * df / gs0
      if (alpha .gt. 1.0_DBL) alpha = 1.0_DBL
      ff = f
      tot = 0.0_DBL
      int = 0
      iexit = 1
  30  continue
      if (ifn .ge. maxfn) go to 90
      do 31 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  31  continue
      call funct(f1,n,w)
      ifn = ifn + 1
      if (f1 .ge. f) go to 40
      f2 = f
      tot = tot + alpha
  32  continue
      do 33 i = 1, n
      x(i) = w(i)
  33  continue
      f = f1
      if (int - 1) 35, 49, 50
  35  continue
      if (ifn .ge. maxfn) go to 90
      do 34 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  34  continue
      call funct(f1,n,w)
      ifn = ifn + 1
      if (f1 .ge. f) go to 50
      if ((f1 + f2 .ge. f + f) .and. (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
      tot = tot + alpha
      alpha = 2.0_DBL * alpha
      go to 32
  40  continue
      if (alpha .lt. aeps) go to 92
      if (ifn .ge. maxfn) go to 90
      alpha = 0.5_DBL * alpha
      do 41 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  41  continue
      call funct(f2,n,w)
      ifn = ifn + 1
      if (f2 .ge. f) go to 45
      tot = tot + alpha
      f = f2
      do 42 i = 1, n
      x(i) = w(i)
  42  continue
      go to 49
  45  continue
      z = 0.1d0
      if (f1 + f .gt. f2 + f2) &
      z = 1.0_DBL + 0.5_DBL * (f - f1) / (f + f1 - f2 - f2)
      if (z .lt. 0.1d0) z = 0.1d0
      alpha = z * alpha
      int = 1
      go to 30
  49  continue
      if (tot .lt. aeps) go to 92
  50  continue
      alpha = tot
      do 56 i = 1, n
      w(i) = x(i)
      w(ib+i) = g(i)
  56  continue
      link = 2
      if (idiff - 1) 100, 100, 110
  54  continue
      if (ifn .ge. maxfn) go to 90
      gys = 0.0_DBL
      do 55 i = 1, n
      w(i) = w(ib+i)
      gys = gys + g(i) * w(is+i)
  55  continue
      df = ff - f
      dgs = gys - gs0
      if (dgs .le. 0.0_DBL) go to 20
      link = 1
      if (dgs + alpha * gs0 .gt. 0.0_DBL) go to 52
      do 51 i = 1, n
      w(iu + i) = g(i) - w(i)
  51  continue
      sig = 1.0_DBL / (alpha * dgs)
      go to 70
  52  continue
      zz = alpha / (dgs - alpha * gs0)
      z = dgs * zz - 1.0_DBL
      do 53 i = 1, n
      w(iu+i) = z * w(i) + g(i)
  53  continue
      sig = 1.0_DBL / (zz * dgs * dgs)
      go to 70
  60  continue
      link = 2
      do 61 i = 1, n
      w(iu+i) = w(i)
  61  continue
      if (dgs + alpha * gs0 .gt. 0.0_DBL) go to 62
      sig = 1.0_DBL / gs0
      go to 70
  62  continue
      sig = -zz
  70  continue
      w(iv+1) = w(iu+1)
      do 71 i = 2, n
      ij = i
      i1 = i - 1
      z = w(iu+i)
      do 72 j = 1, i1
      z = z - h(ij) * w(iv+j)
      ij = ij + n - j
  72  continue
      w(iv+i) = z
  71  continue
      ij = 1
      do 75 i = 1, n
      z = h(ij) + sig * w(iv+i) * w(iv+i)
      if (z .le. 0.0_DBL) z = dmin
      if (z .lt. dmin) dmin = z
      h(ij) = z
      w(ib+i) = w(iv+i) * sig / z
      sig = sig - w(ib+i) * w(ib+i) * z
      ij = ij + np - i
  75  continue
      ij = 1
      do 80 i = 1, n1
      ij = ij + 1
      i1 = i + 1
      do 80 j = i1, n
      w(iu+j) = w(iu+j) - h(ij) * w(iv+i)
      h(ij) = h(ij) + w(ib+i) * w(iu+j)
      ij = ij + 1
  80  continue
      go to (60, 20), link
  90  continue
      iexit = 3
      go to 94
  92  continue
      if (idiff .eq. 2) go to 94
      idiff = 2
      go to 17
  94  continue
      if (iprint .eq. 0) return
      if(verbose) write (6,1005) itn, ifn, iexit
1005  format (1x,'itn = ',i5, ' ifn = ',i5,' iexit = ',i5)
      if(verbose) then
        write (6,1002) f
        write (6,1003) (x(i), i = 1, n)
        write (6,1004) (g(i), i = 1, n)
      endif
      return
 100  continue
      if(verbose) write(6,*) '#######################'
      if(verbose) write(6,*) '### compute gradient at'
      if(verbose) write(6,*) '  x =',w(1:n)        
      if(verbose) write(6,*) 'f(x)=',f
      do 101 i = 1, n
      if(verbose) write(6,*)
      if(verbose) write(6,*) 'direction',i
      z = hh * xm(i)
      w(i) = w(i) + z
      if(verbose) write(6,*) '  x+epsilon =',w(1:n)        
      call funct(f1,n,w)
      if(verbose) write(6,*) 'f(x+epsilon)=',f1
      g(i) = (f1 - f) / z
      if(verbose) write(6,*) 'grad(f,i)(x)=',g(i)
      w(i) = w(i) - z
      if(verbose) write(6,*)
      101  continue
      if(verbose) write(6,*) '### ... done computing gradient at'
      if(verbose) write(6,*) '  x =',w(1:n)        
      if(verbose) write(6,*) 'f(x)=',f
      if(verbose) write(6,*) 'g(x)=',g(1:n)
      ifn = ifn + n
      go to (18, 54), link
 110  continue
      if(verbose) write(6,*) '###############################'
      if(verbose) write(6,*) '### compute gradientgradient at'
      if(verbose) write(6,*) '  x =',w(1:n) 
      do 111 i = 1, n
      if(verbose) write(6,*)
      if(verbose) write(6,*) 'direction',i
      z = hh * xm(i)
      w(i) = w(i) + z
      if(verbose) write(6,*) '  x+epsilon =',w(1:n)        
      call funct(f1,n,w)
      if(verbose) write(6,*) 'f(x+epsilon)=',f1
      w(i) = w(i) - z - z
      if(verbose) write(6,*) '  x-epsilon =',w(1:n)        
      call funct(f2,n,w)
      if(verbose) write(6,*) 'f(x-epsilon)=',f2
      g(i) = (f1 - f2) / (2.0_DBL * z)
      if(verbose) write(6,*) 'gradgrad(f,i)(x)=',g(i)
      w(i) = w(i) + z
 111  continue
      if(verbose) write(6,*) '### ... done computing gradientgradient at'
      if(verbose) write(6,*) '  x =',w(1:n)        
      if(verbose) write(6,*) 'g(x)=',g(1:n)
      ifn = ifn + n + n
      go to (18, 54), link
      end subroutine

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


