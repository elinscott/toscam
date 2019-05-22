
 module mesh

    !-------------!
    use geometry
    use linalg
    use genvar
    use splines, only: resampleit
    use sorting
    use StringManip
    !-------------!

    implicit none
    private
    public :: bin_function
    public :: build1Dmesh
    public :: findin1Dmesh
    public :: findin1Dmesh_d
    public :: mirror_array
    public :: mirror_arrayb__

    INTERFACE build1Dmesh
       MODULE PROCEDURE build1Dmesh_r, build1Dmesh_d, build1Dmesh_dr, build1Dmesh_rd, build1Dmesh_i
    END INTERFACE

    ! INTERFACE build2Dmesh
    !  MODULE PROCEDURE build2Dmesh_r,build2Dmesh_d
    ! END INTERFACE

    INTERFACE findin1Dmesh
       MODULE PROCEDURE findin1Dmesh_r, findin1Dmesh_d
    END INTERFACE

    ! INTERFACE find2Dmesh
    !  MODULE PROCEDURE find2Dmesh_r,find2Dmesh_d
    ! END INTERFACE

    ! INTERFACE find_threshold
    !  MODULE PROCEDURE find_threshold_,find_threshold__,find_threshold___
    ! END INTERFACE

    ! INTERFACE density_distr
    !  MODULE PROCEDURE density_distr_,density_distr__
    ! END INTERFACE

    INTERFACE mirror_array
       MODULE PROCEDURE mirror_array_, mirror_array__, mirror_arrayb_, mirror_arrayb__
    END INTERFACE

 contains
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine mirror_array_sym(xin,yin)
!  implicit none
!  complex(8) :: yin(:),y(size(yin))
!  real(8)    :: xin(:),x(size(xin))
!  integer    :: Nw,iw,miw
!
!   Nw=size(x)
!   DO iw=1,Nw
!    miw=Nw-(iw-1); x(iw)=xin(miw); y(iw)=yin(miw)
!   ENDDO
!   yin=y; xin=x
!
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine renormalize_scale_color_plot(outk,cutoff,cutoff_plot_ak_slope,LOG_SCALE,DISCONTINUED)
!  implicit none
!  real(8)          :: outk(:,:),ar,cutoff,cutoff_plot_ak_slope,cutoff_
!  integer          :: i,j,np,nl
!  logical          :: LOG_SCALE
!  logical,optional :: DISCONTINUED
!
!  np=size(outk,1); nl=size(outk,2)
!
!  if(present(DISCONTINUED))then
!    cutoff_=0.25; ar=cutoff_*4.
!    do i=1,np
!     do j=1,nl
!      if(outk(i,j)>cutoff_)then
!       outk(i,j)=cutoff_+min(1.d0,((outk(i,j))-(cutoff_))/((ar)-(cutoff_)))*0.3*cutoff_
!      endif
!     enddo
!    enddo
!   else
!    ar=maxval(outk)
!     do i=1,np
!      do j=1,nl
!       if(outk(i,j)>cutoff)then
!        if(LOG_SCALE)then
!         outk(i,j)=cutoff + LOG(1.d0+outk(i,j)-cutoff)*cutoff_plot_ak_slope
!        else
!         outk(i,j)=cutoff + (outk(i,j)-cutoff)/(ar-cutoff)*cutoff_plot_ak_slope
!        endif
!       endif
!      enddo
!     enddo
!  endif
!
!  return
!  end subroutine
!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

    subroutine mirror_arrayb_(x, y, xo, yo)
       implicit none
       real(8)           :: x(:), y(:), xx(size(x)), yy(size(y)), xo(:), yo(:)
       xx = x
       yy = y
       call mirror_array_(xx, yy)
       call resampleit(xx, yy, xo, yo, 0.d0)
    end subroutine

    subroutine mirror_arrayb__(x, y, xo, yo)
       implicit none
       complex(8)           :: y(:), yo(:), yy(size(y))
       real(8)              :: x(:), xo(:), xx(size(x))
       xx = x
       yy = y
       call mirror_array__(xx, yy)
       call resampleit(xx, yy, xo, yo, 0.d0)
    end subroutine

    !---------------------!

    subroutine mirror_array_(x, y)
       implicit none
       real(8)           :: x(:), y(:), xx(size(x)), yy(size(y))
       real(8), parameter :: epsilonr = 1.d-5
       integer           :: i, j, i1, i2, ii1, ii2

       xx = x; yy = y
       call reverse_array(xx)
       call reverse_array(yy)
       xx = -xx

       if (xx(1) < x(1) - epsilonr) then
          i1 = fastsearchreal(x(1), xx)
          if (i1 == 0) i1 = 1
       else
          i1 = 1
       endif
       if (xx(size(xx)) > x(size(x)) + epsilonr) then
          i2 = fastsearchreal(x(size(x)), xx)
          if (i2 == 0) i2 = size(x)
       else
          i2 = size(x)
       endif

       if (x(1) < xx(1) - epsilonr) then
          ii1 = fastsearchreal(xx(1), x)
          if (ii1 == 0) ii1 = 1
       else
          ii1 = 1
       endif
       if (x(size(x)) > xx(size(xx)) + epsilonr) then
          ii2 = fastsearchreal(xx(size(xx)), x)
          if (ii2 == 0) ii2 = size(xx)
       else
          ii2 = size(xx)
       endif

       y = 0.d0
       call resampleit(xx(i1:i2), yy(i1:i2), x(ii1:ii2), y(ii1:ii2), 0.d0)

    end subroutine

    !--------------------------!
    !--------------------------!
    !--------------------------!
    !--------------------------!

    subroutine mirror_array__(x, y)
       implicit none
       real(8)           :: x(:), xx(size(x))
       complex(8)        :: y(:), yy(size(y))
       real(8), parameter :: epsilonr = 1.d-5
       integer           :: i, j, i1, i2, ii1, ii2

       xx = x; yy = y
       call reverse_array(xx)
       call reverse_array(yy)
       xx = -xx

       if (xx(1) < x(1) - epsilonr) then
          i1 = fastsearchreal(x(1), xx)
          if (i1 == 0) i1 = 1
       else
          i1 = 1
       endif
       if (xx(size(xx)) > x(size(x)) + epsilonr) then
          i2 = fastsearchreal(x(size(x)), xx)
          if (i2 == 0) i2 = size(x)
       else
          i2 = size(x)
       endif

       if (x(1) < xx(1) - epsilonr) then
          ii1 = fastsearchreal(xx(1), x)
          if (ii1 == 0) ii1 = 1
       else
          ii1 = 1
       endif
       if (x(size(x)) > xx(size(xx)) + epsilonr) then
          ii2 = fastsearchreal(xx(size(xx)), x)
          if (ii2 == 0) ii2 = size(xx)
       else
          ii2 = size(xx)
       endif

       y = 0.d0
       call resampleit(xx(i1:i2), yy(i1:i2), x(ii1:ii2), y(ii1:ii2), 0.d0)

    end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!
!  subroutine derivate_mesh_x(m,n,x,arr,derx)
!  implicit none
!  integer:: n,m
!  real(8) :: arr(m,n),derx(m,n),x(n),xx(n)
!  integer :: i,j,k,l
!  do i=1,m
!   xx=x
!   call derivateit(x(1:n),arr(i,1:n),xx(1:n),derx(i,1:n),1,0.0d0)
!  enddo
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine derivate_mesh_y(m,n,x,arr,derx)
!  implicit none
!  integer:: n,m
!  real(8) :: arr(m,n),derx(m,n),x(m),xx(m)
!  integer :: i,j,k,l
!  do i=1,n
!   xx=x
!   call derivateit(x(1:m),arr(1:m,i),xx(1:m),derx(1:m,i),1,0.0d0)
!  enddo
!  end subroutine
!
!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

    subroutine bin_function(gg, ii, m, maxpos, dd1, dd2)
       implicit none
       integer          :: m, ii(m)
       real(8)          :: gg(:), param(m), d2, d1
       real(8), optional :: maxpos, dd1, dd2
       integer          :: i, j, k, l, iout

       d1 = minval(gg)
       d2 = maxval(gg)
       if (present(dd1)) d1 = dd1
       if (present(dd2)) d2 = dd2

       call build1Dmesh(param(1:m), m, d1 - 1.d-10, d2 + 1.d-10)
       ii = 0

       do i = 1, size(gg)
          call findin1Dmesh(param, m, gg(i), iout)
          if (iout == 0) then
             write (*, *) 'trouble in bin function, not found :', gg(i), d1, d2
             write (*, *) ' m ', m
             write (*, *) 'param : ', param
             stop
          endif
          ii(iout) = ii(iout) + 1
       enddo

       if (present(maxpos)) then
          j = maxloci(ii)
          maxpos = param(j)
       endif

    end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!
!  subroutine combinemesh(x1,f1,x2,f2)
!  implicit none
!  real(8) :: x1(:),f1(:),x2(:),f2(:)
!  integer :: i,j,k,l,m,siz1,siz2,iout
!   siz1=size(x1)
!   siz2=size(x2)
!   do i=1,siz2
!     call findin1Dmesh(x1,siz1,x2(i),iout)
!     if(iout>0) f1(iout)=f1(iout)+f2(i)
!   enddo
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine convert2Dmesh_to_1d(x,y,array,xx,yy,arrayy)
!  implicit none
!   real(8)       :: x(:),y(:),array(:,:)
!   real(8)       :: xx(:),yy(:),arrayy(:)
!   integer       :: i,j,nn
!   nn=0
!   do i=1,size(x)
!    do j=1,size(y)
!     nn=nn+1
!     xx(nn)=x(i)
!     yy(nn)=y(j)
!     arrayy(nn)=array(i,j)
!    enddo
!   enddo
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine convert1Dmesh_to_2d(x,y,array,xx,yy,arrayy,nn)
!  implicit none
!   real(8)       :: x(:),y(:),array(:),dist(nn,nn)
!   real(8)       :: xx(nn),yy(nn),arrayy(nn,nn),xxx,yyy,min,max,miny,maxy,step,vratio,xproj,yproj,scalv
!   integer       :: i,j,iout,jout,nn,nnn,np,count(nn,nn),u(2),v(2)
!
!       np=size(x)
!       min=minval(x)
!       max=maxval(x)
!       step=(max-min)/dble(nn)
!       min=min - step/2.
!       max=max - step/2.
!       miny=minval(y)
!       maxy=maxval(y)
!       step=(maxy-miny)/dble(nn)
!       miny=miny - step/2.
!       maxy=maxy - step/2.
!       call build1Dmesh(xx,nn,min,max)
!       call build1Dmesh(yy,nn,miny,maxy)
!       arrayy=0.d0
!       count=0
!       do i=1,np
!        xxx=x(i)
!        yyy=y(i)
!        call findin1Dmesh(xx,nn,xxx,iout)
!        call findin1Dmesh(yy,nn,yyy,jout)
!        count(iout,jout)=count(iout,jout)+1
!        arrayy(iout,jout)=arrayy(iout,jout)+array(i)
!       enddo
!
!       where(count>0) arrayy=arrayy/dble(count)
!
!       do iout=1,nn
!       do jout=1,nn
!        if(count(iout,jout)==0)then
!
!          do i=1,nn
!           do j=1,nn
!            dist(i,j)=abs(i-iout)+abs(j-jout)
!            if(count(i,j)==0) dist(i,j)=1.d20
!           enddo
!          enddo
!          u=minloc(dist)
!          dist(u(1),u(2))=1.d20
!          v=minloc(dist)
!          vratio =  sqrt ( ( xx(v(1))-xx(u(1)) )**2 + ( yy(v(2))-yy(u(2)) )**2 )
!
!          if(abs(vratio)>1.d-8)then
!            scalv  = (xx(v(1))-xx(u(1)))/vratio * xx(iout) + (yy(v(2))-yy(u(2)))/vratio * yy(iout)
!            xproj  = scalv*(xx(v(1))-xx(u(1)))/vratio
!            yproj  = scalv*(yy(v(2))-yy(u(2)))/vratio
!            vratio = sqrt(  ((xproj-x(u(1)))**2 + (yproj-y(u(2)))**2)  ) / vratio
!            if(vratio>1.d0)vratio=1.d0
!            if(vratio<0.d0)vratio=0.d0
!            arrayy(iout,jout) = arrayy(u(1),u(2)) + (arrayy(v(1),v(2))-arrayy(u(1),u(2)))*vratio
!          else
!            write(*,*) 'resample mesh, 0 vratio occured'
!            arrayy(iout,jout) = arrayy(u(1),u(2))
!          endif
!
!        endif
!       enddo
!       enddo
!
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine area_contained_inside_max_2d_array(x,y,array,area,nn,fourfoldsym,thres2,polyout,nbordout)
!  implicit none
!  real(8)          :: x(:),y(:),array(:),area,xx(nn),yy(nn),arrayy(nn,nn)
!  integer          :: nn,nbord,points(nn*nn,2),order(4*nn*nn),i,j,k,l
!  real(8)          :: slice,threshold,polygon(4*nn*nn,2),angles(4*nn*nn),mat(4,2,2)
!  integer,optional :: fourfoldsym,nbordout
!  real(4)          :: arear
!  real(8),optional :: thres2,polyout(4*nn*nn,2)
!
!  if(nn==0) then
!   write(*,*) 'area contain inside max 2d array , error, nn=0'
!   return
!  endif
!  if(size(array)<1)then
!   write(*,*) 'area contain inside max 2d array , error, nn=0'
!   return
!  endif
!
!  call convert1Dmesh_to_2d(x,y,array,xx,yy,arrayy,nn)
!
!  arrayy=arrayy-minval(arrayy)
!  slice=maxval(arrayy)
!
!                      threshold=abs(slice/10.d0)
!  if(present(thres2)) threshold=abs(slice)/abs(thres2)
!
!  call point2d_inside_threshold_mesh(nn,arrayy,slice,threshold,points,nbord)
!
!  write(*,*) '=============================================='
!  write(*,*) 'threshold                          : ', threshold
!  write(*,*) 'slice                              : ', slice
!  write(*,*) 'area contained inside max boundary : ', nbord
!  write(*,*) 'total points                       : ', nn*nn
!  write(*,*) '=============================================='
!
!  do i=1,nbord
!   polygon(i,1)=xx(points(i,1))
!   polygon(i,2)=yy(points(i,2))
!  enddo
!
!  mat(2,1,1)= 0.
!  mat(2,2,1)=-1.
!  mat(2,1,2)= 1.
!  mat(2,2,2)= 0.
!  do i=3,4
!   mat(i,:,:)=MATMUL(mat(i-1,:,:),mat(2,:,:))
!  enddo
!
!  if(present(fourfoldsym))then
!   do i=2,4
!    do j=nbord*(i-1)+1,nbord*i
!     polygon(j,1:2)=MATMUL(mat(i,:,:),polygon(j-(nbord*(i-1)),1:2))
!    enddo
!   enddo
!   nbord=4*nbord
!  endif
!
!  angles(1:nbord)=(/( ANGLE(polygon(j,1:2)), j=1,nbord )/)
!  call qsort_array(angles(1:nbord),order(1:nbord))
!  call qsort_adj_array(polygon(1:nbord,1),order(1:nbord))
!  call qsort_adj_array(polygon(1:nbord,2),order(1:nbord))
!
!  call polygon_area_2d(nbord,real(polygon(1:nbord,1)),real(polygon(1:nbord,2)),arear)
!  area=arear
!
!  if(present(polyout))then
!   nbordout=nbord
!   polyout=polygon
!  endif
!
!  return
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine point2d_inside_threshold_mesh(nn,array,slice,threshold,points,nbord)
!  implicit none
!  integer :: nn,nbord,i,j
!  real(8) :: threshold,array(nn,nn),slice
!  integer :: points(nn*nn,2)
!
!  points=0; nbord=0
!  do i=1,nn
!  do j=1,nn
!   if( abs( array(i,j)-slice ) < threshold ) then
!    nbord=nbord+1
!    points(nbord,1)=i
!    points(nbord,2)=j
!   endif
!  enddo
!  enddo
!
!  return
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  real(8) function weighted_sum(x,weight)
!  implicit none
!  real(8),intent(in) :: x(:),weight(:)
!  integer :: i,j,k,l
!
!  weighted_sum=0.d0
!  do i=1,size(x)
!   weighted_sum=weighted_sum + weight(i)*x(i)
!  enddo
!  weighted_sum=weighted_sum/sum(weight)
!
!  end function
!
!   !-----------------------------------------!
!
!  real(8) function sum_with_error_bars(x,err)
!  implicit none
!  real(8),intent(in) :: x(:),err(:)
!  real(8) :: weight(size(x))
!  integer :: i,j,k,l
!  weight=1.d0/1.d-5
!  where(abs(err)>1.d-5) weight=1.d0/err
!
!  sum_with_error_bars=0.d0
!  do i=1,size(x)
!   sum_with_error_bars=sum_with_error_bars + weight(i)*x(i)
!  enddo
!  sum_with_error_bars=sum_with_error_bars/sum(weight)
!
!  end function
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!      SUBROUTINE SMOOTH(a,ns)
!      implicit none
!      INTEGER  :: i,j,i1,i2,ii,ns,nw
!      real(8)  :: a(:),a1(size(a))
!
!       nw=size(a)
!
!       DO i=1,nw
!          i1=i-ns
!          IF (i1<1) i1=1
!          i2=i+ns
!          IF (i2>nw) i2=nw
!          a1(i)=0.
!          ii=0
!          DO j=i1,i2
!             ii=ii+1
!             a1(i)=a1(i)+a(j)
!          enddo
!          a1(i)=a1(i)/dble(ii)
!       enddo
!       a=a1
!
!      RETURN
!      END subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
! subroutine build_random_mesh_in_polygon(polygonx,polygony,nmesh,meshx,meshy,borders)
! implicit none
! real(8)          :: meshx(nmesh),meshy(nmesh),xl,xr,yl,yr,rx,ry,polygonx(:),polygony(:)
! real(8)          :: xx,yy,largex,largey
! integer          :: nmesh,NN
! integer          :: i,j,k,l,m
! logical          :: inside
! real(8)          :: point(2),polygon(2,size(polygonx))
! logical,optional :: borders
!
!  NN=size(polygonx)
!  xl=minval(polygonx(:))
!  xr=maxval(polygonx(:))
!  yl=minval(polygony(:))
!  yr=maxval(polygony(:))
!  largex=xr-xl
!  largey=yr-yl
!
!  polygon(1,:)=polygonx
!  polygon(2,:)=polygony
!
!  k=0
!  do
!   rx=drand1(); ry=drand1()
!   xx=xl+rx*largex; yy=yl+ry*largey
!   point(1)=xx; point(2)=yy
!   call polygon_contains_point_2d_(NN,polygon,point,inside)
!   if(inside)then
!     k=k+1
!     meshx(k)=xx; meshy(k)=yy
!   endif
!   if(k==nmesh)exit
!  enddo
!
! if(present(borders))then
!  if(nmesh>size(polygonx))then
!  do i=1,size(polygonx)
!   meshx(i)=polygonx(i)
!   meshy(i)=polygony(i)
!  enddo
!  endif
! endif
!
! return
! end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!   SUBROUTINE density_distr_(discr,Elevel,dens)
!    implicit none
!    real(8)                 :: Elevel(:)
!    INTEGER                :: i,j,k,l,m,jj,jjj,discr
!    real(8)                 :: min,max,resolution
!    real(8)                 :: dens(discr,2),value
!
!    min=minval(Elevel)+2.;max=maxval(Elevel)-2.
!    call build1Dmesh(dens(:,1),discr,min,max)
!
!    resolution=abs((max-min))/100.
!
!    do j=1,discr
!     dens(j,2) = 0.
!      do i=1,size(Elevel)
!       value = sqrt(resolution) / pi / ( (Elevel(i)-dens(j,1))**2 + resolution**2)
!       dens(j,2)=dens(j,2)+value
!      enddo
!    enddo
!    dens(:,2)=dens(:,2) / integrate_non_unif_array(dens(:,1),dens(:,2))
!
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!   SUBROUTINE density_distr__(discr,Elevel,dens)
!    implicit none
!    real(8)                 :: Elevel(:,:)
!    INTEGER                :: i,j,k,l,m,jj,jjj,discr
!    real(8)                 :: min,max,resolution
!    real(8)                 :: dens(discr,2),value
!
!    min=minval(Elevel)-2.;max=maxval(Elevel)+2.
!    call build1Dmesh(dens(:,1),discr,min,max)
!    resolution=abs((max-min))/100.
!
!    do j=1,discr
!     dens(j,2) = 0.
!      do i=1,size(Elevel(:,1))
!      do k=1,size(Elevel(1,:))
!       value = sqrt(resolution) / pi / ( (Elevel(i,k)-dens(j,1))**2 + resolution**2)
!       dens(j,2)=dens(j,2)+value
!      enddo
!      enddo
!    enddo
!    dens(:,2)=dens(:,2) / integrate_non_unif_array(dens(:,1),dens(:,2))
!
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  integer function find_threshold__(x,f,seuil)
!  implicit none
!  real  :: x(:),f(:),seuil,dist
!  integer :: i,j,k
!
!  dist=-1000.
!  k=maxloci(abs(x))
!  do i=1,size(x)
!   if(abs(f(i))>seuil.and.abs(x(i))>dist) then
!     k=i
!     dist=abs(x(i))
!   endif
!  enddo
!
!  find_threshold__=k
!
!  end function
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  integer function find_threshold___(x,f,seuil)
!  implicit none
!  complex(8)  :: f(:)
!  real(8)  :: seuil,x(:),dist
!  integer :: i,j,k
!
!  dist=-1000.
!  k=maxloci(abs(x))
!  do i=1,size(x)
!   if(abs(f(i))>seuil.and.abs(x(i))>dist) then
!     k=i
!     dist=abs(x(i))
!   endif
!  enddo
!
!  find_threshold___=k
!
!  end function
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  integer function find_threshold_(x,f,seuil)
!  implicit none
!  real(8)  :: x(:),seuil,f(:),dist
!  integer :: i,j,k
!
!  dist=-1000.
!  k=maxloci(abs(x))
!  do i=1,size(x)
!   if(abs(f(i))>seuil.and.abs(x(i))>dist) then
!     k=i
!     dist=abs(x(i))
!   endif
!  enddo
!
!  find_threshold_=k
!
!  end function
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!   subroutine normalize_non_unif_array(x,f)
!   implicit none
!   real(8)  :: x(:),f(:),sum
!   integer :: i,n
!    sum=integrate_non_unif_array(x,f)
!    f=f/sum
!   end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!   real(8) function integrate_non_unif_array(x,f,cutoff)
!   implicit none
!   real(8)          :: x(:),f(:),sum,cutoff1
!   integer          :: i,n
!   real(8),optional :: cutoff
!
!      sum=0.d0; n=size(x)
!      if(size(x)/=size(f)) stop 'error integrate_non_unif_array'
!
!      if(present(cutoff))then
!       cutoff1=cutoff
!      else
!       cutoff1=x(n)+1.d-15
!      endif
!
!      do i=2,n
!        if(x(i)<=cutoff1) then
!         sum = sum + abs(x(i)-x(i-1)) * (f(i)+f(i-1)) /2.d0
!        else
!         exit
!        endif
!      enddo
!      integrate_non_unif_array=sum
!   return
!   end function
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
! subroutine create_log_mesh(M1,M2,l,T,omi,oms)
! implicit none
! integer :: M1,M2,inp,last
! real(8)  :: omi(M1+M2+1),oms(M1+M2+1)
! real(8)  :: T,alpha
! integer :: i,j,k,l,m,in
!
! oms=0
! do i=1,M1
!  oms(i)=omi(i)
! enddo
!
! alpha=log( (dble(size(omi)-1)/dble(M1)) ) / dble(M2-1)
! inp = int(0.5*( oms(M1-1)/(pi*T)-1.d0 ))
!
! l=0
! do i=1,M2
!  in = INT(dble(M1)*exp(alpha*dble(i-1))+0.5  )
!  if(in/=inp) then
!    l=l+1
!    oms(M1+l)=dble(2*in+1)*(pi*T)
!  endif
!  inp=in
! enddo
!
! last=INT(0.5*(omi(size(omi)) / (pi*T) -1.d0 ))
! if(inp/=last.and.M1+l<=M1+M2)  then
!  l=l+1
!  oms(M1+l)=dble(2*last+1)*(pi*T)
! endif
!
! return
! end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
! subroutine FILL_MESH(zin)
! implicit none
! real(8),intent(inout)  :: zin(:,:)
! real(8)                :: back(size(zin(:,1)),size(zin(1,:)))
! integer               :: i,j,k,l,m,n,siz1,siz2,u(1)
! real(8)                :: rr,step
! integer               :: rrr(size(zin(:,1))*size(zin(1,:)),2)
! real(8)                :: v1(2),v2(2)
! real(8)                :: dist(size(zin(:,1))*size(zin(1,:)))
!
! !------------------------------------!
! back=0.;k=0
! siz1=size(zin(:,1))
! siz2=size(zin(1,:))
!  do i=1,siz1
!   do j=1,siz2
!    if(abs(zin(i,j))>1.d-10) then
!     k=k+1
!     rrr(k,1)=i
!     rrr(k,2)=j
!    endif
!   enddo
!  enddo
!  if(messages)write(*,*) 'il y a k points dans zin : ',k
! !------------------------------------!
!
! do i=1,siz1
!  do j=1,siz2
!  !************************!
!   if(abs(zin(i,j))<1.d-10)then
!    dist=1.d20
!    v1(1)=dble(i)
!    v1(2)=dble(j)
!    do l=1,k
!     v2(1)=dble(rrr(l,1))
!     v2(2)=dble(rrr(l,2))
!     dist(l)=norme(v1-v2)
!    enddo
!    u=minloc(dist)
!    back(i,j)=zin(rrr(u(1),1),rrr(u(1),2))
!   endif
!  !************************!
!  enddo
! enddo
!
! !------------------------------------!
!
! zin=zin+back
!
! return
! end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
! subroutine duplicate_mesh(xin,zin,T1,T2,xout,yout,zout)
! implicit none
! real(8)                            :: T1(3),T2(3),g1(3),g2(3)
! real(8)                            :: xin(:,:),zin(:)
! real(8),dimension(:),allocatable   :: xout,yout,zinp
! real(8),dimension(:,:),allocatable :: zout,xinp
! integer                            :: Ninp,NX,NY
!
!   !--------------------------------------------------------------!
!   ! periodise la grille originale                                !
!   ! trouve ou sont les points grille originale dans la nouvelle  !
!   ! point i = point le plus proche de la grille originale        !
!   !--------------------------------------------------------------!
!
!     !==================================================================!
!     if(allocated(xinp)) deallocate(xinp,zinp,xout,yout,zout)
!     call periodic_mesh(xin,zin,T1,T2,Ninp)
!     allocate(xinp(Ninp,3),zinp(Ninp))
!     call periodic_mesh(xin,zin,T1,T2,Ninp,xinp,zinp)
!     write(*,*) ' Ninp : ', Ninp
!     !==================================================================!
!     g1=0.;g2=0.;g1(1)=0.1;g2(2)=0.1;
!     call convert_meshvectors_to_matrix(g1,g2,xinp,zinp,NX,NY)
!     write(*,*) 'NX,NY : ', NX,NY
!     allocate(xout(NX),yout(NY),zout(NX,NY))
!     xout=0.;yout=0.;zout=0.
!     call convert_meshvectors_to_matrix(g1,g2,xinp,zinp,NX,NY,xout,yout,zout,.true.)
!     !==================================================================!
!     call FILL_MESH(zout)
!     !==================================================================!
! return
! end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!
!  subroutine build2DmeshGEN(param,mesh,g1,g2,base1,base2)
!  implicit none
!  real(8)  :: param(:,:),param2(size(param(:,1)),size(param(1,:)))
!  real(8)  :: base1(:),base2(:)
!  real(8)  :: min1,max1,min2,max2,g1(2),g2(2)
!  integer :: mesh,i,j,k,count
!   call build2Dmesh(param2,mesh,0.d0,1.d0,0.d0,1.d0)
!   count=0
!   do i=1,mesh
!    do j=1,mesh
!     count=count+1
!     param(count,1)=param2(count,1)*g1(1)+param2(count,2)*g2(1)
!     param(count,2)=param2(count,1)*g1(2)+param2(count,2)*g2(2)
!     if(i==2.and.j==1) base2=param(count,:)
!     if(i==1.and.j==2) base1=param(count,:)
!    enddo
!   enddo
!  return
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!      !-----------------------------!
!
!   subroutine find2Dmesh_(param,xx,iout)
!   implicit none
!   integer  :: iout(2),i,j,k
!   real(8)  :: dist
!   real(8)  :: param(:,:),xx(:),temp1(size(param,1)),temp2(size(param,2))
!     temp1=param(:,1)
!     temp2=param(1,:)
!     iout(1)=fastsearchreal(xx(1),temp1(:))
!     iout(2)=fastsearchreal(xx(2),temp2(:))
!     k=min(iout(1)+1,size(param,1))
!     dist=abs(xx(1)-param(iout(1),1))
!     if(abs(xx(1)-param(k,1))<dist) iout(1)=k
!     k=max(iout(1)-1,1)
!     if(abs(xx(1)-param(k,1))<dist) iout(1)=k
!     k=min(iout(2)+1,size(param,2))
!     dist=abs(xx(2)-param(1,iout(2)))
!     if(abs(xx(2)-param(1,k))<dist) iout(2)=k
!     k=max(iout(2)-1,1)
!     if(abs(xx(2)-param(1,k))<dist) iout(2)=k
!   return
!   end subroutine
!
!      !-----------------------------!
!
!  subroutine find2Dmesh_d(param,xx,iout)
!  implicit none
!  integer  :: iout(2),i,j
!  real(8)  :: param(:,:),xx(:),temp1(size(param,1)),temp2(size(param,2))
!    temp1=param(:,1)
!    temp2=param(1,:)
!    iout(1)=fastsearchreal(xx(1),temp1(:))
!    iout(2)=fastsearchreal(xx(2),temp2(:))
!  return
!  end subroutine
!
!      !-----------------------------!
!
!  subroutine find2Dmesh_r(param,xx,iout)
!  implicit none
!  integer  ::  iout(2),i,j
!  real(4)  ::  param(:,:),xx(:),temp1(size(param,1)),temp2(size(param,2))
!    temp1=param(:,1)
!    temp2=param(1,:)
!    iout(1)=fastsearchreal(xx(1),temp1(:))
!    iout(2)=fastsearchreal(xx(2),temp2(:))
!  return
!  end subroutine
!
!      !-----------------------------!
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine build2Dmesh_d(param,mesh,min1,max1,min2,max2)
!  implicit none
!  real(8)  :: param(:,:),min1,max1,min2,max2
!  integer :: mesh,i,j,k,count
!   count=0
!   do j=1,mesh
!    do i=1,mesh
!     count=count+1
!     param(count,1)= min1+dble(i-1)/dble(mesh-1)*(max1-min1)
!     param(count,2)= min2+dble(j-1)/dble(mesh-1)*(max2-min2)
!    enddo
!   enddo
!  return
!  end subroutine
!
!  subroutine build2Dmesh_r(param,mesh,min1,max1,min2,max2)
!  implicit none
!  real(4) :: param(:,:),min1,max1,min2,max2
!  integer :: mesh,i,j,k,count
!   count=0
!   do j=1,mesh
!    do i=1,mesh
!     count=count+1
!     param(count,1)= min1+dble(i-1)/dble(mesh-1)*(max1-min1)
!     param(count,2)= min2+dble(j-1)/dble(mesh-1)*(max2-min2)
!    enddo
!   enddo
!  return
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine build1DmeshLOG(param,mesh,min,max)
!  real(8) :: param(:),min,max
!  integer :: mesh
!  integer :: i
!   if(mesh<=1.or.abs(min-max)<epsilonr) then
!    param=min
!    return
!   endif
!   param(1:mesh)=(/( min+ (max-min)*Func(i),i=1,mesh)/)
!  contains
!   real(8) function Func(i)
!   integer :: i
!     Func = (1.-LOG( dble(mesh+1-i) ) / LOG ( dble(mesh) ))**2.
!   end function
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine build1Dadapt(xin,yin,param,min1,max1,kkk,cutoff,strong)
!  implicit none
!  real(8)          :: param(:),xin(:),yin(:),min1,max1,derivative(size(param)),mmax,step
!  integer          :: nbin,tot
!  real(8)          :: dist,xout(size(param)),dd,bounda
!  integer          :: mesh
!  integer          :: i,j,kkk2,isi
!  integer,optional :: kkk
!  real(8),optional :: cutoff
!  logical,optional :: strong
!
!   nbin=max(2,size(param)/10)
!
!   if(size(xin)<=1.or.abs(min1-max1)<epsilonr) then
!    param=min1
!    return
!   endif
!
!   if(present(kkk))then
!     kkk2=kkk
!   else
!     kkk2=1
!   endif
!
!   if(present(cutoff)) then
!     bounda=cutoff
!   else
!     bounda=1.d-3
!   endif
!
!   derivative=0.
!   call build1Dmesh(xout(1:nbin+1),nbin+1,min1,max1)
!   call derivateit___(xin,yin,xout(1:nbin+1),derivative(1:nbin+1),kkk2,0.d0)
!
!   do i=1,nbin
!    derivative(i)=(derivative(i)+derivative(i+1))/2.
!   enddo
!
!   mmax=maxval(abs(derivative(1:nbin)))
!   derivative(1:nbin)=abs(derivative(1:nbin))/mmax
!   where(derivative<abs(bounda)) derivative=bounda
!   if(present(strong)) derivative(1:nbin)=derivative(1:nbin)**2
!
!   dd=sum(derivative(1:nbin))
!   do i=1,nbin
!    derivative(i)=MAX(2.d0,derivative(i)/dd*dble(size(param)))
!   enddo
!   derivative(1:nbin)=NINT(derivative(1:nbin))
!   tot=NINT(sum(derivative(1:nbin)))-size(param)
!
!   isi=sign(1,tot)
!   if(tot/=0)then
!   do
!    do j=1,nbin
!     if(derivative(j)>2.01.or.isi<0) then
!      derivative(j)=derivative(j)-dble(1*isi)
!      tot=tot-1*isi
!      if(tot==0) exit
!     endif
!    enddo
!    if(tot==0) exit
!   enddo
!   endif
!
!   if(abs(sum(derivative(1:nbin))-dble(size(param)))>1.d-3) then
!     write(*,*) 'nbin           : ', nbin
!     write(*,*) 'xin            : ', xin(1:10)
!     write(*,*) 'yin            : ', yin(1:10)
!     write(*,*) '    derivative : ', derivative(1:nbin)
!     write(*,*) 'sum derivative : ', sum(derivative(1:nbin))
!     write(*,*) 'size param     : ', size(param)
!     write(*,*) 'diff           : ', tot
!     stop
!   endif
!
!   tot=0
!   do i=1,nbin
!      j=NINT(derivative(i))
!      if(j<2) then
!       write(*,*) '<2 elements in bin number : ', i
!       write(*,*) ' derivative(i) = ', derivative(i)
!       stop 'critical'
!      endif
!      step=(xout(i+1)-xout(i))/dble(2*j)
!      if(i==1)then
!       call build1Dmesh(param(tot+1:tot+j),j,xout(i),xout(i+1)-step)
!      elseif(i>1.and.i<nbin)then
!       call build1Dmesh(param(tot+1:tot+j),j,xout(i)+step,xout(i+1)-step)
!      elseif(i==nbin)then
!       call build1Dmesh(param(tot+1:tot+j),j,xout(i)+step,xout(i+1))
!      endif
!
!      tot=tot+j
!   enddo
!
!   if(tot/=size(param))then
!    write(*,*) 'build adapt, tot different from size(param) : ', tot, size(param)
!    stop
!   endif
!
!   if(abs(param(1)-min1)>1.d-3)then
!     write(*,*) 'error build mesh 1d adapt, first point does not match : '
!     write(*,*) 'param(1:3)  : ',  param(1:3)
!     write(*,*) 'min1        : ',  min1
!     write(*,*) 'derivative1 : ',  derivative(1)
!     write(*,*) 'nbin        : ',  nbin
!     stop 'critical'
!   endif
!
!   if(abs(param(size(param))-max1)>1.d-3)then
!     write(*,*) 'error build mesh 1d adapt, last point does not match : '
!     write(*,*) 'max1                  : ', max1
!     write(*,*) 'xout(nbin+1)          : ', xout(nbin+1)
!     write(*,*) 'param(size(param))    : ', param(size(param))
!     write(*,*) 'param last            : ', param(size(param)-10:size(param))
!     stop 'critical'
!   endif
!
!  return
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine build1Dadapt_resample(xin,yin,param,yout,min,max)
!  implicit none
!  real(8) :: param(:),xin(:),yin(:),yre(size(param)),min,max,derivative(size(param))
!  real(8) :: dist,xout(size(param)),param2(size(param)),yout(:),yout2(size(yout))
!  integer :: mesh
!  integer :: i
!    call build1Dadapt(xin,yin,param2,min,max)
!    call resampleit(xin,yin,param2,yout2,0.d0)
!    param=param2
!    yout=yout2
!  return
!  end subroutine
!
!    !---------------------------------!
!
!  subroutine build1Dadapt_resampleb(xin,yin,kkk,inp)
!  implicit none
!  real(8) :: xin(:),yin(:)
!  real(8) :: min,max,derivative(size(xin)),param(size(xin))
!  real(8) :: yout(size(yin)),xtemp(inp),ytemp(inp)
!  integer :: mesh,inp
!  integer :: i,kkk
!
!    min=minval(xin); max=maxval(xin)
!    call build1Dmesh(xtemp,inp,min,max)
!    call resampleit(xin,yin,xtemp,ytemp,0.d0)
!    call build1Dadapt(xtemp,ytemp,param,min,max,kkk)
!    call resampleit(xtemp,ytemp,param,yout,0.d0)
!
!    xin=param
!    yin=yout
!
!  return
!  end subroutine
!
!    !---------------------------------!
!
!  subroutine build1Dadapt_resamplec(xin,yin,peakr)
!  implicit none
!  real(8) :: xin(:),yin(:),yback(size(yin)),peakr
!  real(8) :: min1,max1
!  real(8) :: xtemp(size(xin)),ytemp(size(xin)),xtemp2(size(xin)),ytemp2(size(xin))
!  integer :: mesh,inp,peak1,peak2
!  integer :: i,kkk,l1,l2
!  real(8) :: large,center1,center2,a1,a2,b1,b2
!
!    inp=size(xin)
!    min1=minval(xin); max1=maxval(xin)
!    call build1Dmesh(xtemp,inp,min1,max1)
!    call resampleit(xin,yin,xtemp,ytemp,0.d0)
!
!    yback=ytemp
!    peak1=maxloci(abs(yback))
!    yback(peak1-inp/20:peak1+inp/20)=0.
!    peak2=maxloci(abs(yback))
!
!    if(ytemp(peak2)>ytemp(peak1))then
!     peakr=xtemp(peak2)
!    else
!     peakr=xtemp(peak1)
!    endif
!
!    large=abs(max1-min1)/80.
!    center1=xtemp(peak1)
!    center2=xtemp(peak2)
!    if(center1>center2) call swap(center1,center2)
!    if(center1>center2) stop 'error'
!
!    a1=center1-large
!    if(a1<min1) a1=min1
!    a2=center1+large
!    b1=center2-large
!    if(a2>b1) then
!     call swap(a2,b1)
!     a2=a2-large/10.
!     b1=b1+large/10.
!    endif
!    b2=center2+large
!    if(b2>max1) b2=max1
!
!    l1=1;l2=inp/10
!    call build1Dmesh(xtemp2(l1:l2),l2-l1+1,min1,a1+epsilonr)
!    l1=inp/10+1;l2=2*inp/5
!    call build1Dmesh(xtemp2(l1:l2),l2-l1+1,a1+epsilonr,a2)
!    l1=2*inp/5+1;l2=3*inp/5
!    call build1Dmesh(xtemp2(l1:l2),l2-l1+1,a2+epsilonr,b1)
!    l1=3*inp/5+1;l2=9*inp/10
!    call build1Dmesh(xtemp2(l1:l2),l2-l1+1,b1+epsilonr,b2)
!    l1=9*inp/10+1;l2=inp
!    call build1Dmesh(xtemp2(l1:l2),l2-l1+1,b2+epsilonr,max1+2.*epsilonr)
!
!    call resampleit(xtemp,ytemp,xtemp2,ytemp2,0.d0)
!    xin=xtemp2; yin=ytemp2
!
!  return
!  end subroutine
!
!    !---------------------------------!
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine build1DmeshLOGsym(param,meshtot,acenter,min,max)
!  implicit none
!  real(8) :: param(:),min,max,acenter
!  integer :: mesh,meshtot
!  integer :: i
!
!   if(min*max>0.)then
!    call build1DmeshLOG(param,meshtot,min,max)
!    return
!   endif
!
!   mesh=meshtot/2;
!   if(mesh<=1.or.abs(min-max)<epsilonr) then
!     param=min
!     return
!   endif
!
!   call build1DmeshLOG(param(1:mesh),mesh,abs(acenter),-min)
!   call reverse_array(param(1:mesh))
!   param=-param
!   call build1DmeshLOG(param(mesh+1:meshtot),meshtot-mesh,abs(acenter),max)
!
!  end subroutine
!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

    subroutine build1Dmesh_i(param, mesh, min, max)
       integer :: param(:), min, max
       integer :: mesh
       integer :: i
       real(8) :: param2(mesh), min2, max2
       if (mesh <= 1 .or. min == max) then
          param = min
          return
       endif
       min2 = dble(min)
       max2 = dble(max)
       call build1Dmesh_d(param2, mesh, min2, max2)
       param = NINT(param2)
    end subroutine

    subroutine build1Dmesh_d(param, mesh, min, max)
       real(8) :: param(:), min, max
       integer :: mesh
       integer :: i
       if (mesh <= 1 .or. abs(min - max) < epsilonr) then
          param = min
          return
       endif
       param(1:mesh) = (/(min + dble(i - 1)/dble(mesh - 1)*(max - min), i=1, mesh)/)
    end subroutine

    subroutine build1Dmesh_r(param, mesh, min, max)
       real :: param(:), min, max
       integer :: mesh
       integer :: i
       if (mesh <= 1 .or. abs(min - max) < 1.e-6) then
          param = min
          return
       endif
       param(1:mesh) = (/(min + dble(i - 1)/dble(mesh - 1)*(max - min), i=1, mesh)/)
    end subroutine

    subroutine build1Dmesh_dr(param, mesh, min, max)
       real(8)  :: param(:)
       real(4) :: min, max
       integer :: mesh
       integer :: i
       if (mesh <= 1 .or. abs(min - max) < 1.e-6) then
          param = min
          return
       endif
       param(1:mesh) = (/(min + dble(i - 1)/dble(mesh - 1)*(max - min), i=1, mesh)/)
    end subroutine

    subroutine build1Dmesh_rd(param, mesh, min, max)
       real(4) :: param(:)
       real(8)  :: min, max
       integer :: mesh
       integer :: i
       if (mesh <= 1 .or. abs(min - max) < epsilonr) then
          param = min
          return
       endif
       param(1:mesh) = (/(min + dble(i - 1)/dble(mesh - 1)*(max - min), i=1, mesh)/)
    end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

    subroutine findin1Dmesh_d(param, mesh, xx, iout)
       integer :: iout, mesh, i
       real(8) :: param(:), xx
       iout = fastsearchreal(xx, param)
    end subroutine

    subroutine findin1Dmesh_r(param, mesh, xx, iout)
       integer :: iout, mesh, i
       real :: param(:), xx
       iout = fastsearchreal(xx, param)
    end subroutine

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!
!  subroutine periodic_mesh(xin,zin,T1,T2,Nout,xout,zout)
!  implicit none
!  !-------------------------------------------!
!  ! T1,T2 : periodic vectors for translations !
!  !-------------------------------------------!
!  integer                         :: Nout
!  real(8),intent(in)              :: xin(:,:)
!  real(8)                         :: xinb(size(xin(:,1)),size(xin(1,:)))
!  real(8),intent(in)              :: zin(:)
!  real(8),dimension(:,:),optional :: xout
!  real(8),dimension(:),optional   :: zout
!  integer                         :: siz1,siz2,i,j,k,l,m,jj,count
!  real(8)                         :: a1,a2,b1,b2,T1(3),T2(3)
!  real(8)                         :: T1b(2),T2b(2),v(3)
!  logical                         :: out
!    xinb=xin
!    call shift_to_corner_mesh(xinb,T1b(1:2),T2b(1:2))
!    call findNmesh(xinb,T1b,T2b,Nout,T1(1:2),T2(1:2))
!    if(present(xout))then
!     call findNmesh(xinb,T1b(1:2),T2b(1:2),Nout,T1(1:2),T2(1:2),zout,zin,xout)
!    endif
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine convert_meshvectors_to_matrix(g1,g2,xin,zin,Nx,Ny,xout,yout,zout,fill)
!  implicit none
!  real(8)           :: xin(:,:),zin(:),g1(:),g2(:),base(2,2)
!  integer           :: Nx,Ny,N1,N2
!  real(8),optional  :: xout(Nx),yout(Ny),zout(Nx,Ny)
!  integer           :: i,j,k,l,m,n,siz1
!  integer           :: decomp(size(xin(:,1)),2)
!  real(8)           :: decompr(size(xin(:,1)),2)
!  integer           :: x1,y1
!  logical,optional  :: fill
!
!  base(1,:)= g1(1:2);base(2,:)= g2(1:2)
!  siz1=size(xin(:,1))
!  write(*,*) 'size in convert_meshvectors : ' , siz1
!  if(siz1/=size(zin(:))) then
!    write(*,*) 'error routine convert_meshvectors : ', siz1,size(zin)
!    stop
!  endif
!
!  write(*,*) ' CONVERT TO 2D PERIODIC MESH , SIZE : ', size(zin(:)),size(xin(:,1))
!  do i=1,siz1
!   call decomposevec(decompr(i,:),xin(i,1:2),base)
!  enddo
!  decomp=NINT(decompr)
!
!  if(testing)then
!   if(maxval(abs(decompr-dble(decomp)))>1.d-1) then
!    write(*,*) 'base1 : ',g1
!    write(*,*) 'base2 : ',g2
!    write(*,*) 'coordonnee : '
!    do i=1,siz1
!     write(*,*) '---------------------------'
!     write(*,*) 'initial : ',decompr(i,:)
!     write(*,*) 'decomp  : ', decomp(i,:)
!    enddo
!    stop 'error in convert_meshvectors_to_matrix (module mesh) decomposition not integer'
!   endif
!  endif
!
!  decomp=decomp+1
!  Nx=maxval(decomp(:,1))
!  Ny=maxval(decomp(:,2))
!
!  if(present(xout))then
!   xout=0.;yout=0.;zout=0.
!   if(present(fill))then
!    do i=1,Nx
!     xout(i)=dble(i-1)*g1(1)
!    enddo
!    do i=1,Ny
!     yout(i)=dble(i-1)*g2(2)
!    enddo
!    if(abs(g1(2))>1.d-5.or.abs(g2(1))>1.d-5) &
! & write(*,*)'MIGHT BE A PROBLEM IN ROUTINE convert_meshvector... in mesh.f90'
!   endif
!
!   do i=1,siz1
!    x1=decomp(i,1);y1=decomp(i,2)
!    if(messages3.and.abs( xin(i,1) - xout(x1) )>1.d-2) then
!     write(*,*) 'error in convert mesh'
!     write(*,*) 'x1,xin(i,1),xout(x1) : ', x1,xin(i,1),xout(x1)
!    endif
!    if(messages3.and.abs( xin(i,2) - yout(x1) )>1.d-2) then
!     write(*,*) 'error in convert mesh'
!     write(*,*) 'y1,xin(i,2),yout(y1) : ', y1,xin(i,2),yout(y1)
!    endif
!    xout(x1)=xin(i,1)
!    yout(y1)=xin(i,2)
!    if(messages3.and.abs(zout(x1,y1))>1.d-3) then
!     write(*,*) 'x,y : ', xin(i,1),xin(i,2)
!     write(*,*) i,x1,y1
!     write(*,*) zout(x1,y1)
!     write(*,*) 'DANGER, error for zout'
!    endif
!    zout(x1,y1)=zin(i)
!
!    if(messages4) write(*,*) 'step  : ' ,i
!    if(messages4) write(*,*) 'xin   : ' ,xin(i,1:2)
!    if(messages4) write(*,*) 'x1,y1 : ' ,x1,y1
!    if(messages4) write(*,*) 'zin   : ' ,zin(i)
!    if(messages4) write(*,*) 'zout  : ' ,zout(x1,y1)
!   enddo
!  endif
!
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine convert_meshmatrix_to_vectors(xin,yin,zin,Nx,Ny,xout,zout)
!  implicit none
!  real(8) :: xin(:),yin(:),zin(:,:)
!  integer :: i,j,k,l,m
!  integer :: x1,y1,Nx,Ny
!  integer :: xout(Nx*Ny,3),zout(Nx*Ny)
!  i=0
!  do x1=1,Nx
!  do y1=1,Ny
!  i=i+1
!   xout(i,1)=xin(x1)
!   xout(i,2)=yin(y1)
!   xout(i,3)=0.
!   zout(i)=zin(x1,y1)
!  enddo
!  enddo
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
! subroutine shift_to_corner_mesh(xin,T1b,T2b)
! implicit none
! real(8),intent(inout) :: xin(:,:)
! real(8)               :: a1,a2,b1,b2,T1b(2),T2b(2)
!  a1=minval(xin(:,1))
!  a2=maxval(xin(:,1))
!  b1=minval(xin(:,2))
!  b2=maxval(xin(:,2))
!  xin(:,1)=xin(:,1)-a1
!  a2=a2-a1
!  xin(:,2)=xin(:,2)-b1
!  b2=b2-b1
!  T1b=0.;T2b=0.
!  T1b(1)=a2
!  T2b(2)=b2
! end subroutine
!
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
! subroutine findNmesh(xin,T1b,T2b,Nout,T1,T2,zout,zin,xout)
! implicit none
! real(8),optional :: T1(2),T2(2)
! real(8),optional :: zout(:),zin(:)
! real(8),optional :: xout(:,:)
! real(8)          :: xin(:,:)
! logical          :: out
! integer          :: Nout,i,j,k,l,jj,siz1
! real(8)          :: T1b(2),T2b(2),TT1(2),TT2(2),v(3)
! integer          :: a1,a2
!  siz1=size(xin(:,1))
!  if(present(T1)) then
!   a1=1;TT1=T1;TT2=T2
!  else
!   a1=0;TT1=0;TT2=0
!  endif
!  Nout=0
!  if(present(zout)) zout=0.
!  if(present(xout)) xout=0.
!  do i=1,siz1
!   do j=-a1,a1
!    do jj=-a1,a1
!      v(1)=xin(i,1)+j*TT1(1)+jj*TT2(1)
!      v(2)=xin(i,2)+j*TT1(2)+jj*TT2(2)
!      v(3)=xin(i,3)
!      call inout(v(1),v(2),T1b,T2b,out)
!      if(.not.out) then
!       Nout=Nout+1
!       if(present(zout)) zout(Nout)=zin(i)
!       if(present(xout)) xout(Nout,:)=v(:)
!      endif
!    enddo
!   enddo
!  enddo
! return
! end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
!  subroutine shift_down_min_of_array(array,percent)
!  implicit none
!  real(8) :: percent, array(:,:)
!  integer :: u(2)
!  real(8) :: amin,amax
!    amin=minval(array)
!    amax=maxval(array)
!    u=minloc(array)
!    array(u(1),u(2)) = array(u(1),u(2))-percent*(amax-amin)
!  end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
!
! subroutine rescale_array_moins1_1(array)
! implicit none
! real(8) :: array(:,:)
! integer :: i,j,k,l,m
! real(8) :: r1,r2,t,u,v,ar,am
!    ar=minval(array)
!    array=array-ar
!    am=maxval(array)
!    array=array*2.d0/am-1.d0
! end subroutine
!
!  !-----------------!
!
! subroutine rescale_array_atanh(array)
! implicit none
! real(8) :: array(:,:)
!    call rescale_array_moins1_1(array)
!    array=atanh(array)
! end subroutine
!
!  !-----------------!
!
! subroutine rescale_array_tanh(array)
! implicit none
! real(8) :: array(:,:)
!    call rescale_array_moins1_1(array)
!    array=tanh(array)
! end subroutine
!
!  !-----------------!
!
! subroutine rescale_array_tan(array)
! implicit none
! real(8) :: array(:,:)
!    call rescale_array_moins1_1(array)
!    array=tanh(array*pi)
! end subroutine
!
!  !-----------------!
!
! subroutine rescale_array_cos(array)
! implicit none
! real(8) :: array(:,:)
!    call rescale_array_moins1_1(array)
!    array=tanh((array+1)*pi/2.-pi)
! end subroutine
!
!  !-----------------!
!
! subroutine rescale_array_atan(array)
! implicit none
! real(8) :: array(:,:)
! integer :: i,j,k,l,m
!    call rescale_array_moins1_1(array)
!    array=atan(array)
! end subroutine
!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!
! !********************************************************************!

 end module
