
module kramers_kronig_mod

 !--------------!
   use genvar
   use splines
   use linalg
 !--------------!

 INTERFACE principal_part
  MODULE PROCEDURE principal_part_,principal_part__,principal_part___
 END INTERFACE

 INTERFACE kramers_kronig
  MODULE PROCEDURE kramers_kronig_,kramers_kronig__
 END INTERFACE


 contains

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

 subroutine matsum_(gr, gi, iom, x0, wb)
 IMPLICIT NONE
  real(8), intent(out) :: gr(:), gi(:)
  real(8), intent(in)  :: iom(:),x0(:),wb(:)
  INTEGER              :: nom, nx
  real(8)              :: omn, w1, sumr, sumi,dh
  INTEGER              :: n,k

  nom=size(gr)
  nx=size(x0)
  gr=0.d0 ; gi=0.d0

  do n=1,nom
     omn  = iom(n)
     sumr = 0.d0 ; sumi = 0.d0
     do k=1,nx
        if(k<nx)then
         dh=x0(k+1)-x0(k)
        else
         dh=x0(nx)-x0(nx-1)
        endif
        w1   = wb(k) / ( omn**2 + x0(k)**2 )
        sumr = sumr + w1 * x0(k) * dh
        sumi = sumi + w1 * omn   * dh
     enddo
     gr(n) = -sumr ; gi(n) = -sumi
  enddo

 return
 end subroutine

  !-------------------------------------------!

 subroutine matsum(gr, gi, iom, x0, dh, wb, nom, nx)
 IMPLICIT NONE
  real(8), intent(out) :: gr(nom), gi(nom)
  real(8), intent(in)  :: iom(nom),x0(nx),dh(nx),wb(nx)
  INTEGER, intent(in)  :: nom, nx
  real(8)              :: omn, w1, sumr, sumi
  INTEGER              :: n,k

  gr=0.d0 ; gi=0.d0
  do n=1,nom
     omn  = iom(n)
     sumr = 0.d0 ; sumi = 0.d0
     do k=1,nx
        w1   = wb(k) / ( omn**2 + x0(k)**2 )
        sumr = sumr + w1 * x0(k) * dh(k)
        sumi = sumi + w1 * omn   * dh(k)
     enddo
     gr(n) = -sumr ; gi(n) = -sumi
  enddo

 return
 end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

SUBROUTINE kramarskronig(gr, fi, om, i0, nom)
  IMPLICIT NONE
  real(8), intent(out):: gr
  real(8), intent(in) :: fi(nom)
  real(8), intent(in) :: om(nom)
  INTEGER, intent(in) :: i0, nom
  real(8) :: sm
  real(8) :: dh(nom)
  real(8) :: S0, x
  INTEGER :: j, i0f

  if(nullarray(om).or.constantarray(fi)) return

  i0f = i0+1
  S0  = fi(i0f)
  x   = om(i0f)

  dh(1) = 0.5*(om(2)-om(1))
  do j=2,nom-1
     dh(j) = 0.5*(om(j+1)-om(j-1))
  enddo
  dh(nom) = 0.5*(om(nom)-om(nom-1))

  sm=0
  do j=1,i0f-2
     sm = sm + (fi(j)-S0)*dh(j)/(om(j)-x)
  enddo

  if (i0f>1)   sm = sm + (fi(i0f-1)-S0)*(dh(i0f-1)+0.5*dh(i0f))/(om(i0f-1)-x)

  if (i0f<nom) sm = sm + (fi(i0f+1)-S0)*(dh(i0f+1)+0.5*dh(i0f))/(om(i0f+1)-x)

  do j=i0f+1,nom
     sm = sm + (fi(j)-S0)*dh(j)/(om(j)-x)
  enddo

  if (.NOT.(i0f .EQ. 1 .OR. i0f.EQ.nom)) then
     sm = sm + S0*log(abs((om(nom)-x)/(x-om(1))))
  endif

  gr = sm/pi

  RETURN

END SUBROUTINE 


!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

subroutine test_kramers_kronig
implicit none
integer :: i,j,k,l,m
real(8) :: x(1000),f(1000),cauchy(1000),cauchy2(1000),cauchy3(1000),cauchy4(1000)

 do i=1,1000
  x(i)=-10. + 20./999.*dble(i-1) 
  f(i)= EXP(-x(i)**2.d0)
 enddo
 
 do i=1,1000
  call principal_part(x,f,x(i),-10.d0,10.d0,cauchy(i))
  call principal_part(x(i),-10.d0,10.d0,cauchy2(i),func) 
 enddo
 call principal_part(x,f,x,cauchy3)

 call kramers_kronig_(x,f,cauchy4,realtoim=.true.)

 stop


contains

 real(8) function func(x)
   real(8),intent(in) :: x 
   func = EXP(-x**2 )
 end function

end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

! Kramers - Kronig relation :

subroutine KramsKron(Frc, om, F0, wb, x0, dhx)
  IMPLICIT NONE
  complex(8), intent(out)  :: Frc(:)
  real(8), intent(in)      :: om(:)
  real(8), intent(in)      :: F0(:)
  real(8), intent(in)      :: wb(:)
  real(8), intent(in)      :: x0(:)
  real(8), intent(in)      :: dhx(:)
  INTEGER                  :: nom, nx
  INTEGER                  :: j, k
  real(8)                  :: dd, omj, wj, sums, Fre, Fri

  nx=size(x0) ; nom=size(F0)
 
  do j=1,nom
     omj = om(j); wj  = F0(j)
     sums=0.d0
     do k=1,nx
       dd   = ( omj-x0(k) ) 
       if(abs(dd)<1.d-5) dd=1.d-5
       sums = sums + (wb(k)-wj)*dhx(k) / dd
     enddo
     dd     =  omj-x0(1)
     if(abs(dd)<1.d-5) dd=1.d-5
     Fre    =  sums - wj * dlog(abs((x0(nx)-omj)/dd))
     Fri    = -pi*F0(j)
     Frc(j) =  cmplx(Fre,Fri,kind=8)
  enddo

 return
 end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

subroutine kramers_kronig__(xout,yout,funcint,realtoim)
implicit none
integer           ::  npoints,i,j,k,l
real(8)           ::  C,a,b,intg,rerror
real(8)           ::  xout(:),dh0(size(xout)),yout(:)
logical           ::  plotdata,realtoim

 interface
   real(8) function funcint(x)
   implicit none
    real(8),intent(in) :: x
    end function
 end interface

 plotdata=.false.
 npoints=size(xout)

 call center_function(4,100.d0,0.008d0,xout,dh0,yout,funcint)
 a=minval(xout); b=maxval(xout)
 do i=1,npoints 
  call principal_part__(xout(i),a,b,yout(i),funcint)
 enddo
 if(realtoim)then
  yout=yout/( -pi )
 else
  yout=yout/(  pi)
 endif

return
end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

subroutine kramers_kronig_(xin,yin,yout,realtoim)
implicit none
integer           ::  npoints,nin,ns
real(8)           ::  C,a,b,intg,rerror
real(8)           ::  xin(:),yin(:),yout(:),xcopy(size(xin))
real(8)           ::  i,j,k,l
logical           ::  plotdata,realtoim

 if(nullarray(xin).or.nullarray(yin)) return

 plotdata=.false.
 xcopy=xin+0.0001
 ns=size(xin)

 xcopy(1)  = (xcopy(1)    + xcopy(2)  ) /2.d0
 xcopy(ns) = (xcopy(ns-1) + xcopy(ns) ) /2.d0

 call principal_part___(xin,yin,xcopy,yout) 

 if(realtoim)then
  yout=yout/( -pi )
 else
  yout=yout/(  pi)
 endif

return
end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

subroutine principal_part__(C2,a,b,intg,funcint)
implicit none
real(8)           ::  C,C2,a,b,intg,rerror
real(8)           ::  i,j,k,l
integer,parameter ::  LIMIT=1000000
integer           ::  NEVAL,IER,LAST,IWORK(LIMIT)
real(8)           ::  WORK(4*LIMIT+2)
logical           ::  plotdata

interface
  real(8) function funcint(x)
  implicit none
   real(8),intent(in) :: x
  end function
end interface

   if(abs(C2-a)<100.*error .or.abs(C2-b)<100.*error )then
    C=C2+100.*error
   else
    C=C2
   endif
   call DQAWC (funcint, a, b, C, 0.00001d0, 0.0001d0, intg, rerror, NEVAL, IER, size(IWORK), SIZE(WORK), LAST, IWORK, WORK)

   if(IER>0.and.strongstop) stop 'error in principal_part'

return
end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

subroutine principal_part_(om,im,C2,a,b,intg)
implicit none
real(8)           ::  om(:),im(:),C,C2,a,b,intg,rerror
real(8)           ::  i,j,k,l
type(spline)      ::  splineim
integer,parameter ::  LIMIT=1000000
integer           ::  NEVAL,IER,LAST,IWORK(LIMIT)
real(8)           ::  WORK(4*LIMIT+2)
logical           ::  plotdata

if(abs(C2-a)<100.*error .or.abs(C2-b)<100.*error )then
 C=C2+100.*error
else
 C=C2
endif

!-------------------------------------------------------!
! QAWC adaptive integration for Cauchy principal values !
!-------------------------------------------------------!

!            The routine calculates an approximation result to a
!            Cauchy principal value I = INTEGRAL of F*W over (A,B)
!            (W(X) = 1/((X-C), C.NE.A, C.NE.B), hopefully satisfying
!            following claim for accuracy
!            ABS(I-RESULT).LE.MAX(EPSABE,EPSREL*ABS(I)).
!
!        Computation of a Cauchy principal value
!        Standard fortran subroutine
!
!        PARAMETERS
!         ON ENTRY
!            F      - Real
!                     Function subprogram defining the integrand
!                     Function F(X). The actual name for F needs to be
!                     declared E X T E R N A L in the driver program.
!            A      - Real
!            B      - Real
!                     limit of integration
!            C      - Parameter in the weight function, C.NE.A, C.NE.B.
!                     If C = A or C = B, the routine will end with
!                     IER = 6 .
!            EPSABS - Real
!                     Absolute accuracy requested
!            EPSREL - Real
!                     Relative accuracy requested
!         ON RETURN
!            RESULT - Real
!                     Approximation to the integral
!            ABSERR - Real
!                     Estimate or the modulus of the absolute error,
!                     Which should equal or exceed ABS(I-RESULT)
!            NEVAL  - Integer
!                     Number of integrand evaluations
!
!            IER    - Integer
!                     IER = 0 Normal and reliable termination of the
!                             routine. It is assumed that the requested
!                             accuracy has been achieved.
!                     IER.GT.0 Abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. It is assumed that the
!                             requested accuracy has not been achieved.
!            ERROR MESSAGES
!                     IER = 1 Maximum number of subdivisions allowed
!                             has been achieved. One can allow more sub-
!                             divisions by increasing the value of LIMIT
!                             (and taking the according dimension
!                             adjustments into account). However, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties.
!                             If the position of a local difficulty
!                             can be determined (e.g. SINGULARITY,
!                             DISCONTINUITY within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling
!                             appropriate integrators on the subranges.
!                         = 2 The occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 Extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 The input is invalid, because
!                             C = A or C = B or
!                             (EPSABS.LE.0 and
!                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
!                             or LIMIT.LT.1 or LENW.LT.LIMIT*4.
!                             RESULT, ABSERR, NEVAL, LAST are set to
!                             zero.  Except when LENW or LIMIT is
!                             invalid, IWORK(1), WORK(LIMIT*2+1) and
!                             WORK(LIMIT*3+1) are set to zero, WORK(1)
!                             is set to A and WORK(LIMIT+1) to B.
!
!         DIMENSIONING PARAMETERS
!            LIMIT - Integer
!                    Dimensioning parameter for IWORK
!                    LIMIT determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (A,B), LIMIT.GE.1.
!                    If LIMIT.LT.1, the routine will end with IER = 6.
!
!           LENW   - Integer
!                    Dimensioning parameter for WORK
!                    LENW must be at least LIMIT*4.
!                    If LENW.LT.LIMIT*4, the routine will end with
!                    IER = 6.
!
!            LAST  - Integer
!                    On return, LAST equals the number of subintervals
!                    produced in the subdivision process, which
!                    determines the number of significant elements
!                    actually in the WORK ARRAYS.
!
!         WORK ARRAYS
!            IWORK - Integer
!                    Vector of dimension at least LIMIT, the first K
!                    elements of which contain pointers
!                    to the error estimates over the subintervals,
!                    such that WORK(LIMIT*3+IWORK(1)), ... ,
!                    WORK(LIMIT*3+IWORK(K)) form a decreasing
!                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2),
!                    and K = LIMIT+1-LAST otherwise
!
!            WORK  - Real
!                    Vector of dimension at least LENW
!                    On return
!                    WORK(1), ..., WORK(LAST) contain the left
!                     end points of the subintervals in the
!                     partition of (A,B),
!                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
!                     the right end points,
!                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
!                     the integral approximations over the subintervals,
!                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
!                     contain the error estimates.
!-------------------------------------------------------!
 
     call init_spline(splineim,om,im,5,0.0001d0)
     plotdata=.false.
     write(*,*) a,b,C,LIMIT
     call DQAWC (funcint, a, b, C, 0.00001d0, 0.0001d0, intg, rerror, NEVAL, IER, size(IWORK), SIZE(WORK), LAST, IWORK, WORK)
     if(IER>0.and.strongstop) stop 'error in principal_part'
     call kill_spline(splineim)

return

contains

 !---------------------!
 !---------------------!
 !---------------------!
 !---------------------!

  real(8) function funcint(x) 
  implicit none
   real(8),intent(in)  :: x
   real(8)             :: xx(1),yy(1)  
   xx=x
   call evaluate_spline(splineim,xx,yy)
   funcint=yy(1) 
  end function

 !---------------------!
 !---------------------!
 !---------------------!
 !---------------------!

end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

subroutine principal_part___(om,im,omout,imout)
implicit none
real(8)            ::  om(:),im(:),omout(:),imout(:)
real(8)            ::  C,C2,a,b,intg,rerror
integer(4)         ::  i,j,k,l
type(spline)       ::  splineim
integer,parameter  ::  LIMIT=1000000
integer            ::  NEVAL,IER,LAST,IWORK(LIMIT)
real(8)            ::  WORK(4*LIMIT+2)
logical            ::  plotdata

    a=minval(om); b=maxval(om)
    if(nullarray(om).or.nullarray(im))return

    call init_spline(splineim,om,im,5,0.0001d0)
    plotdata=.false.

    do i=1,size(omout)
     C=omout(i) 
     if(abs(C-a)<100.*error) C=C+100.*error
     if(abs(C-b)<100.*error) C=C-100.*error

     call DQAWC (funcint, a, b, C, 0.00001d0, 0.0001d0, intg, rerror, NEVAL, IER, size(IWORK), SIZE(WORK), LAST, IWORK, WORK)
     if(IER>0.and.strongstop) stop 'error in principal_part' 
     imout(i)=intg
    enddo
    call kill_spline(splineim)

return

contains

 !---------------------!
 !---------------------!
 !---------------------!
 !---------------------!
 
  real(8) function funcint(x)
  implicit none
   real(8),intent(in) :: x
   real(8)            :: xx(1),yy(1)
   xx=x
   call evaluate_spline(splineim,xx,yy)
   funcint=yy(1)
  end function

 !---------------------!
 !---------------------!
 !---------------------!
 !---------------------!

end subroutine


!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

end module
