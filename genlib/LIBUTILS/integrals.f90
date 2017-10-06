module integrals

 use linalg
 use splines
 use mesh
 use FuncAnalysis

contains

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine spline_and_integrate_0_to_infty(xx,ff,fac,integ,logmesh,funct)
 implicit none
 real(8) :: xx(:),ff(:),integ
 integer :: i,j,k,l,fac,siz,siz_,sizb
 real(8) :: xxb(size(xx)+1),ffb(size(ff)+1)
 real(8) :: xx_(fac*size(xx)),ff_(fac*size(ff))
 logical :: logmesh

 interface
  real(8) function funct(x)
  real(8) x
  end function
 end interface
 optional :: funct

  siz  = size(xx)
  siz_ = size(xx_)

  if(size(ff)/=siz) stop 'error spline and integrate : sizes' 

  if(fac>1.and.abs(xx(1))>1.d-4) then
   sizb=siz+1
   xxb(1)=0.0001d0
   call spline_inter_extrapolate(xx(1:5),ff(1:5),xxb(1),ffb(1))
   ffb(2:siz+1) = ff(1:siz)
   xxb(2:siz+1) = xx(1:siz)
  else
   sizb=siz
   ffb(1:siz) = ff(1:siz)
   xxb(1:siz) = xx(1:siz)
  endif 

 if(.not.logmesh)then
  call build1Dmesh(xx_(1:siz_),siz_,xxb(1),xxb(sizb))
 else
  call build1DmeshLOG(xx_(1:siz_),siz_,xxb(1),xxb(sizb))
 endif

  if(fac>1.or.logmesh)then
    call resampleit(xxb(1:sizb),ffb(1:sizb),xx_,ff_,0.0d0) 
  else
    xx_=xx
    ff_=ff
  endif

  if(present(funct)) then
   do i=1,siz_
     ff_(i)=ff_(i)*funct(xx_(i))
   enddo
  endif

 if(.not.logmesh.and..not.fac==1)then
  call simptab (ff_, xx_(2)-xx_(1), size(xx_), integ)
 else
  integ=0.d0
  do i=1,siz_-1
   integ=integ+ ff_(i) * (xx_(i+1)-xx_(i)) 
  enddo 
 endif

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
!*****************************************************!

end module
