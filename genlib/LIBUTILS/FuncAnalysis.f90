module FuncAnalysis
use genvar

implicit none
contains  
 
!--------------------------------------------------------------!
 !subroutine bissect (a,b,fct,pasmax,eps,racine)              !
 !function longueur (fctx,fcty,ta,tb,pasmax,eps)              !
 !function minimum (a,b,fct,pasmax,eps)                       !
 !subroutine simpson (a,b,fct,pasmax,eps,somme)               !
 !subroutine troispoints (a,b,fct,pasmax,eps,somme)           !
 !subroutine rk4 (fct, xi, yi, xf, yf, x, y, n, itmax, eps)   ! 
 !subroutine rk4_multi(fct,xi,yi,xf,yf,x,y,n,Nvar,itmax,eps)  !
 !subroutine simptab (y, h, n, sum)                           !
!--------------------------------------------------------------!


!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!

subroutine bissect3(a,b,fct,nstep,racine,decide,critere)
implicit none
real(8),intent (in)  :: a,b
real(8),intent (out) :: racine
real(8)              :: critere
integer              :: nzero,zeros(nstep)

interface
 function fct (x)
  real(8),intent (in) :: x
  real(8)             :: fct
 end function 
 function decide (x)
  real(8),intent (in) :: x
  real(8)             :: decide
 end function 
end interface

integer              :: nstep
real(8)              :: meshx(nstep),funcmesh(nstep),funcfin(nstep)
integer              :: i,j,kk,u(1)

 meshx    = (/( a + (b-a)*dble(i-1)/dble(nstep-1), i=1,nstep )/)
 funcmesh = (/( fct(meshx(i)) , i=1,nstep )/)

 nzero=0
 do i=1,nstep

   u = minloc(abs(funcmesh))

   if(nzero==0.or.funcmesh(u(1))<critere)then
     nzero=nzero+1
     zeros(nzero)=u(1)
     funcmesh(u(1))=1.d20
     funcfin(nzero)=decide(meshx(u(1)))
   else 
     exit
   endif 

 enddo

 u = minloc(funcfin(1:nzero))
 racine=meshx(zeros(u(1)))

end subroutine



subroutine bissect2(a,b,fct,nstep,racine)
implicit none
real(8),intent (in)  :: a,b
real(8),intent (out) :: racine
interface
 function fct (x)
  real(8),intent (in) :: x
  real(8)             :: fct
 end function fct
end interface
integer              :: nstep
real(8)              :: meshx(nstep),funcmesh(nstep)
integer              :: i,j,kk,u(1)

 meshx    = (/( a + (b-a)*dble(i-1)/dble(nstep-1), i=1,nstep )/)
 funcmesh = (/( fct(meshx(i)) , i=1,nstep )/)
 u = minloc(abs(funcmesh))
 racine=meshx(u(1))

end subroutine

!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
 
subroutine bissect(a,b,fct,pasmax,eps,racine)
implicit none
real(8),intent (in)  :: a,b,eps
real(8)              :: aa,bb
integer              :: pasmax
real(8),intent (out) :: racine
interface
 function fct (x)
  real(8),intent (in) :: x
  real(8)             :: fct
 end function fct
end interface
real(8)              :: xp,xn,xm
integer              :: pas

 aa=fct(a); bb=fct(b)

 if(aa*bb>0) then
   write(*,*) 'Bissect : conditions aux extremites !'
   write(*,*) 'f(a),f(b):', aa,bb
   write(*,*) 'a,b : ',a,b
   stop
 endif

 if(aa<0.)then
  xn=a;xp=b
 else
  xn=b;xp=a
 endif

 pas = 1
 xm  = (a+b)/2.

 do
  if (pas >= pasmax) then
   write(*,*) 'Bissect : pas assez d iterations !'
   write(*,*) 'racine is ', xm 
   write(*,*) 'function is : ', fct(xm)
   racine=xm
   return
  endif
  
   aa=fct(xm)
   if (aa < -eps) then
    xn = xm
   else if ( aa > eps ) then
    xp = xm
   else if ( abs(aa) <eps ) then
    racine = xm
   endif
  
   xm = (xn+xp)/2.
  
   if(abs((xm-xp)/xm)<eps) then
    racine = xm
    exit
   endif
  
   pas = pas+1
 enddo

 return
 end subroutine bissect


!***********************************************************! 
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
 
 
function longueur (fctx,fcty,ta,tb,pasmax,eps)
implicit none
real(8),intent (in)    :: ta,tb,eps
integer,intent (in)   :: pasmax
real(8) :: longueur
interface
  function fctx (x)
   real(8),dimension (:),intent (in)  :: x
   real(8),dimension(size(x)):: fctx
  end function
  function fcty (x)
   real(8),dimension(:),intent (in):: x
   real(8),dimension (size(x))  :: fcty
  end function
end interface
real(8),dimension (:),allocatable :: t,xt,yt,delx,dely,dels
integer                          :: pas,subdiv,i
real(8)                           :: newlong,oldlong,err

oldlong = 0.0 ; pas = 1
!
do
  pas = pas + 1
  subdiv = 2**(pas-1)
  allocate (t (subdiv+1),xt (subdiv+1),yt (subdiv+1))
  allocate (delx (subdiv),dely (subdiv),dels(subdiv))
  t = (/ (ta + (i-1)*((tb-ta)/subdiv),i = 1,subdiv+1) /)
  xt = fctx (t)
  yt = fcty (t)
  delx = xt(2:subdiv+1) - xt(1:subdiv)
  dely = yt(2:subdiv+1) - yt(1:subdiv)
  dels = sqrt (delx*delx + dely*dely)
  
  newlong = sum (dels)
  deallocate (t,xt,yt,delx,dely,dels)
  err = newlong - oldlong
  
  if ( err <= eps ) exit
  if ( pasmax <= pas ) then
 write (*,*) 'Longueur : nombre d iterations insuffisant'
 stop 'termine'
  endif
  if ( oldlong > newlong ) then
 write (*,*) 'Longueur : incoherence !'
 stop 'termine'
  endif
  
  oldlong = newlong
enddo

longueur = newlong

 write (*,*)
 write (*,*) 'La longueur est :  ',longueur
 write (*,*)
 write (*,*) 'Erreur absolue :  ',err
 write (*,*)
 write (*,*) 'Nombre d iterations :  ',pas
 write (*,*)

return
end function longueur
 

!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
 
function minimum(a,b,fct,pasmax,eps,epsy)
implicit none
real(8),intent (in)  :: a,b,eps
real(8),optional     :: epsy
integer,intent (in) :: pasmax
real(8)              :: minimum,ra,rb,rd,rg
interface
  function fct(x)
   real(8) :: x
   real(8) :: fct
  end function
end interface
real(8)              :: xa,xb,xg,xd
integer             :: pas

xa=a;xb=b;
xg=xa+(xb-xa)/3.d0; xd=xb-(xb-xa)/3.d0
pas=1;ra=fct(xa);rb=fct(xb);

 do
  if (pas==pasmax )then
   write(*,*) 'Minimum : danger,pas assez d iterations !'
   minimum=(xa+xb)/2.d0
   write(*,*) 'minimum : ', minimum
   write(*,*) 'fmin    : ', fct(minimum)
   return
  endif

  !ra=fct(xa); rb=fct(xb); 
   rd=fct(xd); rg=fct(xg)

  if(ra<=rg)then
   xb=xg;rb=rg
  elseif(rd>=rg)then
   xb=xd;rb=rd
  elseif(rb>=rd)then
   xa=xg;ra=rg
  elseif(rb<rd)then
   xa=xd;ra=rd
  endif
  
  xg=xa+(xb-xa)/3.d0; xd=xb-(xb-xa)/3.d0
 
  if (abs((xb-xa)/xa)<eps) then
   minimum = (xa + xb)/2.d0
   exit 
  endif
  if(present(epsy))then
   if (abs(ra)<epsy.and.abs(rb)<epsy) then
    minimum = (xa + xb)/2.d0
    exit
   endif
  endif
  
  pas=pas+1
 enddo

 if(messages4) then
  write (*,*) '###############################'
  write (*,*) 'Le minimum est      :  ',minimum
  write (*,*) 'Erreur relative     :  ',abs((xb-xa)/xa)
  write (*,*) 'Erreur absolue      :  ',abs(xb-xa)
  write (*,*) 'Nombre d iterations :  ',pas-1
  write (*,*) '###############################'
 endif

return
end function minimum

!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!

subroutine simpson (a,b,fct,pasmax,eps,somme)
implicit none
real(8),intent (in)  :: a,b,eps
real(8),intent (out) :: somme
integer,intent (in) :: pasmax
interface
  function fct (x)
   real(8),intent (in)  :: x
   real(8)              :: fct
  end function
end interface
integer                           :: n,pas,i
real(8)                           :: h,sigma1,sigma2
real(8)                           :: sigma4,newsum,oldsum
real(8),dimension (:),allocatable :: x4,f4

if ( b <= a ) then
  write (*,*) 'Simpson : bornes d integrations !'
  stop 'termine'
endif

pas=1; n=1; h=b-a
sigma1 = fct(a) + fct(b)
sigma2 = 0.0
sigma4 = fct( (a+b)/2.0 )
newsum = h*( sigma1 + 2.0*sigma2 + 4.0*sigma4 ) / 6.0

 do
  pas = pas + 1
  if ( pas >= pasmax ) then
   write (*,*) 'Simpson : pas assez d iterations !'
   stop 'termine'
  endif
  
  oldsum = newsum
  n = 2*n
  h = h / 2.0
  allocate ( x4 (n),f4 (n) )
  x4 = (/(a - (h/2.0) + i*h, i=1,n)/)
  f4 = (/(fct (x4(i)),i=1,n)/)
  sigma2 = sigma2 + sigma4
  sigma4 = sum (f4)
  deallocate ( x4,f4 )
  newsum = h*( sigma1 + 2.0*sigma2 + 4.0*sigma4 ) / 6.0
  
  if ( eps >= abs(((newsum-oldsum)/15.0)/newsum) ) exit
 enddo

somme = newsum + (newsum-oldsum)/15.0

end subroutine
 

!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
 
subroutine troispoints (a,b,fct,pasmax,eps,somme)
implicit none
real(8),intent(in)  :: a,b,eps
integer,intent(in) :: pasmax
real(8),intent(out) :: somme
interface
  function fct (x)
   real(8),intent (in)  :: x
   real(8)              :: fct
  end function
end interface
real(8),dimension (:),allocatable :: x,fx
real(8)                           :: h,sigma1,sigma2
real(8)                           :: newsum,oldsum
integer                          :: pas,i,n

h = b-a
sigma1 = fct ( (a+b)/2.0 )
sigma2 = fct ( (3.0*a + b)/4.0 ) + fct ( (a+3.0*b)/4.0 )
oldsum = h*( 2.0*sigma2 - sigma1)/ 3.0
pas = 1; n = 1

do
  if ( pas == pasmax ) then
   write (*,*) 'Trois points : pas assez d iterations !'
   stop 'termine'
  endif
  
  n = 2*n; h = h/2.0; sigma1 = sigma2
  allocate (x (2*n),fx (2*n) )
  x = (/ (a+(2*i-1)*(h/4.0),i=1,2*n)/)
  fx = (/ (fct (x(i)),i=1,2*n)/)
  sigma2 = sum (fx)
  deallocate (x,fx)
  newsum = h*(2.0*sigma2 - sigma1) / 3.0
  
  if ( abs(((newsum-oldsum)/15.0)/newsum) < eps ) then
   somme = newsum
   exit
  endif
  
  oldsum = newsum
  pas = pas + 1
enddo
end subroutine troispoints
 
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
 
subroutine rk4 (fct, xi, yi, xf, yf, x, y, n, itmax, eps)
implicit none
integer, intent(in)             :: itmax
real(8), intent(in)              :: xi, yi, xf, eps
real(8), intent(out)             :: yf
integer, intent(out)            :: n
real(8), pointer, dimension (:)  :: x , y ,xt, yt 
interface
 function fct (x,y)
  real(8) :: x,y
  real(8) :: fct
 end function fct
end interface
integer                         ::  i, npas, it
real(8)                          ::  t1,t2,t3,t4,xc,yc,yf_old,h,err

npas=5; it=-1

 do 
 
  it=it+1
 
  if(it>itmax) then
   print *, 'rk4: pas de convergence itmax trop petit'
   stop
  endif

  npas=npas*2; h= (xf - xi)/ npas  
  
  allocate (xt(npas+1),yt(npas+1))
  xt(1)=xi; yt(1)=yi
  xc=xi; yc=yi
  
 do i=1,npas
  t1 = h*fct(xc,yc)
  t2 = h*fct(xc+(h/2.),yc + (t1)/2)
  t3 = h*fct(xc+(h/2.),yc + (t2)/2)
  t4 = h*fct(xc+ h,yc + t3)
  xc= xc + h
  yc= yc + (1./6.)*(t1 +2*t2 + 2*t3 + t4)
  xt(i+1)=xc
  yt(i+1)=yc
 enddo

 yf=yc
 print *, 'yc vaut:', yf

  if(it>=1) then
  err=abs(  (yf-yf_old)/(15.*yf))
 
  if(err<= eps) then
   yf= ( 16.*yf - yf_old) / 15.
   allocate (x(npas+1),y(npas+1))
   x=xt
   y=yt
   n=npas
   exit
  endif
 
  endif
  
 deallocate(xt,yt)
 yf_old=yf
  
 enddo

end subroutine rk4  

!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!

subroutine rk4_multi(fct,xi,yi,xf,yf,x,y,n,Nvar,itmax,eps)
implicit none
integer, intent (in)                :: itmax, Nvar
real(8), intent(in)                  :: xi, xf
integer, intent (out)               :: n
real(8), dimension (:), intent (in)  :: yi, eps 
real(8), dimension (:), intent (out) :: yf  
real(8), pointer, dimension (:)      :: x  
real(8), pointer, dimension (:,:)    :: y  

interface
 function fct (x, y)
  real(8), intent (in)                :: x
  real(8), dimension (:), intent (in) :: y  
  real(8), dimension (size(y))        :: fct  
 end function fct 
end interface

integer                               :: it,j,npas
real(8)                                :: h,xc
real(8), dimension (:), allocatable    :: err,yc,yf_old
real(8), dimension (:), allocatable    :: xt
real(8), dimension (:,:), allocatable  :: yt 
real(8), dimension (:,:), allocatable  :: t 
logical, dimension (:), allocatable   :: masque


allocate(err(1:Nvar),yc(1:Nvar),yf_old(1:Nvar),masque(1:Nvar))

it = -1; npas= 8

 do
  it = it + 1

  if (it > itmax) then
   write (*,*) 'rk4: pas de convergence'
   stop
  endif

 npas = npas*2 
 
 h = (xf - xi) / npas
 
 allocate ( xt(1:npas+1 ) , yt(1:Nvar,1:npas+1) )
 
 allocate (t(1:Nvar,1:4))

 xt(1) = xi
 yt(:,1) = yi(:)

 xc=xi
 yc=yi(:)
 
 do j =  1, npas
  t(:,1) = h * fct( xc, yc )
  t(:,2) = h * fct( xc + 0.5*h , yc + t(:,1)/2. )
  t(:,3) = h * fct( xc + 0.5*h , yc + t(:,2)/2. )
  t(:,4) = h * fct( xc, yc + t(:,3) )
  
  yc = yc + (t(:,1) + 2*t(:,2) + 2*t(:,3) + t(:,4)) / 6.
  xc=xc+h
  
  xt(j+1)  = xc
  yt(:,j+1)= yc
 enddo
 
 yc(:) = yt(:,npas+1)
 
 if (it >= 2) then
 where (eps<0)  err = abs((yc - yf_old) / (15.*yc))  
 where (eps>0)  err = abs ((yc - yf_old) / (15.))
  
  masque = .false.
  WHERE (err < abs(eps)) masque = .true.
  if (all(masque)) then
   n = npas
   yf(:) = (16.*yc - yf_old)/15. 
   allocate ( x(1:n+1) , y(1:Nvar, 1:n+1) )
   x = xt
   y = yt
   deallocate (xt,yt,t)
   exit
  endif
 endif

 yf_old = yc
 deallocate (xt,yt,t)
enddo

deallocate(err,yc,yf_old,masque)
return
end subroutine rk4_multi


!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!

subroutine simptab (y, h, n, sum) 
 !Integrale avec entree de tableaux!
  implicit none
  integer, intent (in)                :: n
  real(8), intent (in)                :: h
  real(8), intent (in), dimension (:) :: y
  real(8), intent (out)               :: sum
  integer                             :: i
  real(8)                             :: wgt
  sum = y(1)
  wgt = 2.0
! h= (xf - xi)/ n
  do i = 2, n
  if (wgt == 2.0) then
   wgt = 4.0
  else
   wgt = 2.0
  endif
  sum = sum + wgt*y(i)
  enddo
  sum = sum + y(n+1)
  sum = h*sum/3.0
end subroutine 
 
 
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
!***********************************************************!
 
end module 
