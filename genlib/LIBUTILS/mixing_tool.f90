
module mixing

contains

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

subroutine Simple_Mix(F, x, alpha, nsize)
  IMPLICIT NONE
  !Passed variables
  real(8), intent(in)    :: F(nsize)
  real(8), intent(inout) :: x(nsize)
  real(8), intent(in)    :: alpha
  INTEGER, intent(in)    :: nsize
  INTEGER                :: i
  do i=1,nsize
     x(i) = x(i) + F(i)*alpha
  enddo
end subroutine Simple_Mix

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

subroutine Broyden_Mix(filename, dim, F, x, restart0, alpha0, alpha)
  IMPLICIT NONE
  !Passed variables
  CHARACTER*100, intent(in)      :: filename
  INTEGER, intent(in)            :: dim
  real(8),  intent(in)           :: F(dim)
  real(8),  intent(inout)        :: x(dim)
  LOGICAL, intent(in)            :: restart0
  real(8),  intent(in)           :: alpha0
  real(8),  intent(in), optional :: alpha
  real(8)                        :: xold(dim), Fold(dim), Jt(dim,dim)
  real(8)                        :: dx(dim), dF(dim), Jdx(dim)
  INTEGER                        :: ipiv(dim)
  LOGICAL                        :: file_exists, restart
  INTEGER,PARAMETER              :: fh=1001
  INTEGER                        :: ios, i
  real(8)                        :: nrm, fct, alpha_

  restart = .true.
  if (.NOT.restart0) then
     INQUIRE(file=filename,exist=file_exists)
     if (file_exists) then
        open(fh,file=filename,status='old',form='unformatted')
        read(fh,iostat=ios) xold, Fold, Jt
        close(fh)
        if (ios.EQ.0) then
           restart = .false.
        else
           print *, 'READING WAS NOT SUCCESSFUL'
        endif
     endif
  endif
  if (restart) print *, 'RESTARTING!'
  !Update Jacoobian Jt
  IF (restart) THEN
     Jt=0
     DO i=1,size(Jt,1)
        Jt(i,i)=1/alpha0
     ENDDO
  ELSE
     dx = x-xold
     dF = F-Fold
     nrm = 1/sqrt(dot_product(dx,dx))
     dx = dx*nrm
     dF = dF*nrm
     CALL ProductMV(Jdx,Jt,dx,dim,dim)
     Jdx = Jdx+dF
     fct = -1.
     CALL TensorProduct(Jt,Jdx,dx,fct,dim,dim)
  ENDIF
  !Store data for next itteration
  open(fh,file=filename,status='unknown',form='unformatted')
  write(fh) x, F, Jt
  close(fh)
  !Find size of a step
  alpha_ = 1
  if (present(alpha)) alpha_ = alpha
  if (restart) alpha_ = alpha0
  !Find new x
  dx = F
  if (.not.restart) CALL SolveSOLA(Jt,dx,ipiv,dim)
  x = x + alpha_*dx
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
!**************************************************************************

subroutine Vanderbilt_Mix(filename, dim, F, x, restart0, alpha0, alpha)
  IMPLICIT NONE

  CHARACTER*100, intent(in)     :: filename
  INTEGER, intent(in)           :: dim
  real(8), intent(in)           :: F(dim)
  real(8), intent(inout)        :: x(dim)
  LOGICAL, intent(in)           :: restart0
  real(8), intent(in)           :: alpha0
  real(8), intent(in), optional :: alpha
  real(8)                       :: xold(dim), Fold(dim), Jt(dim,dim)
  real(8)                       :: dx(dim), dF(dim), Gamma(dim,dim), Phi(dim,dim)
  INTEGER                       :: ipiv(dim)
  LOGICAL                       :: file_exists, restart
  INTEGER,PARAMETER             :: fh=1001
  INTEGER                       :: ios, i
  real(8)                       :: nrm, fct, wm, wghs0, c, alpha_

  wm = 1
  restart = .true.
  INQUIRE(file=filename,exist=file_exists)
  if (.NOT.restart0) then
     if (file_exists) then
        open(fh,file=filename,status='old',form='unformatted')
        read(fh,iostat=ios) xold, Fold, Gamma, Phi
        close(fh)
        if (ios.EQ.0) then
           dx = x-xold
           nrm = sqrt(dot_product(dx,dx))
           if (nrm.GT.1e-12) restart = .false.
        else
           print *, 'READING WAS NOT SUCCESSFUL'
        endif
     endif
  endif
  if (restart) print *, 'RESTARTING!'
  !Update Jacoobian Jt
  IF (restart) THEN
     Jt=0
     DO i=1,size(Jt,1)
        Jt(i,i)=1/alpha0
     ENDDO
     wghs0=0.01
     Phi = 0
     DO i=1,dim
        Phi(i,i) = 1/wghs0**2
     ENDDO
     Gamma = Jt*(wghs0**2)
  ELSE
     dx = x-xold
     dF = F-Fold
     nrm = 1/sqrt(dot_product(dx,dx))
     dx = dx*(nrm*wm)
     dF = dF*(nrm*wm)
     fct = -1.
     CALL TensorProduct(Gamma,dF,dx,fct,dim,dim);
     CALL ProductMV(xold,Phi,dx,dim,dim)
     c = 1/(1+dot_product(dx,xold));
     CALL TensorProduct(Phi,xold,xold,-c,dim,dim);
     CALL ProductMM(Jt,Gamma,Phi,dim,dim,dim)
  ENDIF
  !Store data for next itteration
  open(fh,file=filename,status='unknown',form='unformatted')
  write(fh) x, F, Gamma, Phi
  close(fh)

  !Find size of a step

  alpha_ = 1
  if (present(alpha)) alpha_ = alpha
  if (restart) alpha_ = alpha0

  !Find new x
  dx = F
  if (.not.restart) CALL SolveSOLA(Jt,dx,ipiv,dim)
  x = x + alpha_*dx

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
!**************************************************************************

subroutine Johnson_Mix(xnew, filename, F, x, restart0, alpha0, alpha, rdim)
  IMPLICIT NONE
  !Passed variables
  INTEGER, intent(in)           :: rdim
  real(8), intent(out)          :: xnew(rdim)
  CHARACTER*500, intent(in)     :: filename
  real(8), intent(in)           :: F(rdim)
  real(8), intent(in)           :: x(rdim)
  LOGICAL, intent(in)           :: restart0
  real(8), intent(in)           :: alpha0
  real(8), intent(in)           :: alpha
  INTEGER,PARAMETER             :: fh=1001
  real(8)                       :: xold(rdim), Fold(rdim), Gt(rdim,rdim)
  real(8)                       :: dx(rdim), dF(rdim), Gamma(rdim,rdim), Phi(rdim,rdim)
  INTEGER                       :: ipiv(rdim)
  LOGICAL                       :: file_exists, restart
  INTEGER                       :: ios, i
  real(8)                       :: nrm, fct, wm, wghs0, c, alpha_

  wm = 1
  restart = .true.
  if (.NOT.restart0) then
     INQUIRE(file=TRIM(filename),exist=file_exists)
     if (file_exists) then
        open(fh,file=TRIM(filename),status='old',form='unformatted')
        read(fh,iostat=ios) xold, Fold, Gamma, Phi
        close(fh)
        if (ios.EQ.0) then
           dF = F-Fold
           nrm = sqrt(dot_product(dF,dF))
           if (nrm.GT.1e-18) then 
              restart = .false.
           else
              print *, 'dF.dF < 1e-18'
           endif
        else
           print *, 'READING WAS NOT SUCCESSFUL'
        endif
     else
        print *, 'file does not exists!'
     endif
  endif
  if (restart) print *, 'RESTARTING!'


  !Update inverse Jacoobian Gt
  IF (restart) THEN
     Gt=0
     DO i=1,size(Gt,1)
        Gt(i,i)=alpha0
     ENDDO
     wghs0=0.01
     Phi = 0
     DO i=1,rdim
        Phi(i,i) = 1/wghs0**2
     ENDDO
     Gamma = Gt*(wghs0**2)
  ELSE
     dx = x-xold
     dF = F-Fold

     nrm = 1/sqrt(dot_product(dF,dF))

     dx = dx*(nrm*wm)
     dF = dF*(nrm*wm)
     fct = -1.d0
     CALL TensorProduct(Gamma,dx,dF,fct,rdim,rdim);
     CALL ProductMV(Fold,Phi,dF,rdim,rdim)
     c = 1.d0/(1+dot_product(dF,Fold));
     CALL TensorProduct(Phi,Fold,Fold,-c,rdim,rdim);
     CALL ProductMM(Gt,Gamma,Phi,rdim,rdim,rdim)
  ENDIF

  !Store data for next itteration
  open(fh,file=TRIM(filename),status='unknown',form='unformatted')
  write(fh) x, F, Gamma, Phi
  close(fh)

  !Find size of a step
  alpha_ = alpha
  if (restart) alpha_ = alpha0
  !Find new x
  dx = F
  if (.not.restart) CALL ProductMV(dx,Gt,F,rdim,rdim)
  !    print *, 'alpha_=', alpha_
  xnew = x + alpha_*dx

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
!**************************************************************************

subroutine SolveSOLA(A, b, ipiv, ndim)
  IMPLICIT NONE
  real(8), intent(inout) :: A(ndim,ndim)
  real(8), intent(inout) :: b(ndim)
  INTEGER, intent(inout) :: ipiv(ndim)
  INTEGER, intent(in)    :: ndim
  INTEGER                :: info
  CALL dgetrf(ndim,ndim,A,ndim,ipiv,info)
  IF (info.NE.0) print *, 'ERROR in dgetrf. Code=', info
  CALL dgetrs('N',ndim,1,A,ndim,ipiv,b,ndim,info)
  IF (info.NE.0) print *, 'ERROR in dgetrs. Code=', info
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
!**************************************************************************

subroutine TensorProduct(A,x,y,alpha,m,n)
  IMPLICIT NONE ! Computes A += x * y
  real(8), intent(out) :: A(m,n)
  real(8), intent(in)  :: x(m)
  real(8), intent(in)  :: y(n)
  INTEGER, intent(in)  :: m, n
  real(8), intent(in)  :: alpha
  CALL dger(m,n,alpha,x,1,y,1,A,m)
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
!**************************************************************************

subroutine ProductMV(y,A,x,m,n)
  IMPLICIT NONE
  real(8), intent(out) :: y(m)
  real(8), intent(in)  :: A(m,n)
  real(8), intent(in)  :: x(n)
  INTEGER, intent(in)  :: m, n
  real(8)              :: alpha, beta
  alpha = 1.
  beta = 0.
  CALL dgemv('N',m,n,alpha,A,m,x,1,beta,y,1)
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
!**************************************************************************

subroutine ProductMM(C,A,B,m,n,k)
  IMPLICIT NONE
  real(8), intent(out) :: C(m,n)
  real(8), intent(in)  :: A(m,k)
  real(8), intent(in)  :: B(k,n)
  INTEGER :: m,n,k
  real(8) :: alpha, beta
  alpha = 1.
   beta = 0
  CALL dgemm('N','N',m,n,k,alpha,A,m,B,k,beta,C,m)
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
!**************************************************************************

end module
