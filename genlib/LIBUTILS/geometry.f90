
module geometry

  use sorting
  use random
  use linalg
  use matrix
  use genvar
  use geomlib

  INTERFACE inout
     MODULE PROCEDURE inout_,inout__
  END INTERFACE

contains


!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  real(8) FUNCTION Distance(v0, v1, indim)
    IMPLICIT NONE
    complex(8), intent(in) :: v0(indim), v1(indim)
    INTEGER, intent(in)    :: indim
    INTEGER                :: i, j
    real(8)                :: ds
    ds = 0.0
    DO i=1,indim
       ds = ds + abs(v0(i)-v1(i))
    END DO
    Distance = ds/indim
  END FUNCTION 

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

function periodist_vec_open_bc(x,y,T1,T2,T3,open)
implicit none
integer :: i,j,ii,jj,i_
integer :: m,k,l
real(8)  :: x(3),y(3),yy(3),T1(3),T2(3),T3(3)
real(8)  :: min,minb,periodist_vec_open_bc
logical :: open

 min=100000.
 do i=-1,1
  do j=-1,1
   do i_=-1,1

    if(.not.open.or.(i==0.and.j==0.and.i_==0))then
      yy=y+dble(i)*T1+dble(j)*T2+dble(i_)*T3
      minb=norme(yy-x)
      if(minb<min)then
       l=i
       m=j
       min=minb
      endif
    endif

   enddo
  enddo
 enddo
 periodist_vec_open_bc=min

return
end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

function periodist_vec(x,y,T1,T2,T3)
implicit none
integer :: i,j,ii,jj,i_
integer :: m,k,l
real(8) :: x(3),y(3),yy(3),T1(3),T2(3),T3(3)
real(8) :: min,minb,periodist_vec

 min=100000.
 do i=-1,1
  do j=-1,1
   do i_=-1,1
      yy=y+dble(i)*T1+dble(j)*T2+dble(i_)*T3
      minb=norme(yy-x)
      if(minb<min)then
       l=i
       m=j
       min=minb
      endif
   enddo
  enddo
 enddo
 periodist_vec=min

return
end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

 subroutine rapproche_point(A,B,scale)
 implicit none
 real(8)   :: A(:),B(:),scale,V(size(A)),W(size(A)),C(size(A)),D(size(A))
  V=B-A
  C=A+V*(1.-scale)/2. 
  D=B-V*(1.-scale)/2.
  A=C
  B=D
 end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

 subroutine dist_point_to_polygon(dx,T1,T2,disto)
 implicit none
 real(8) :: T1(2),T2(2),disto(4),dx(2)
 real :: dist(4)
  call line_seg_point_dist_2d(0.,0.,real(T1(1)),real(T1(2)),real(dx(1)),real(dx(2)),dist(1))
  call line_seg_point_dist_2d(0.,0.,real(T2(1)),real(T2(2)),real(dx(1)),real(dx(2)),dist(2))
  call line_seg_point_dist_2d(real(T1(1)),real(T1(2)), real(T1(1)+T2(1)),real(T1(2)+T2(2)),&
                            & real(dx(1)),real(dx(2)),dist(3))
  call line_seg_point_dist_2d(real(T2(1)),real(T2(2)), real(T2(1)+T1(1)),real(T2(2)+T1(2)),&
                            & real(dx(1)),real(dx(2)),dist(4))
 disto=dist
 return
 end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  !---------------------------------------------------!
  ! length of a polygon from point 1 to NN<=nb points !
  !---------------------------------------------------!

 real(8) function lengthpoly(curvein,NN)
 real(8) :: curvein(NN,3)
 integer :: i,j,k,l,NN
  lengthpoly=0
  do i=1,NN
   if(i+1/=NN)then
    lengthpoly=lengthpoly + norme( curvein(i+1,:)-curvein(i,:) )
   else
    lengthpoly=lengthpoly + norme( curvein(1,:)-curvein(NN,:) )
   endif
  enddo

 return
 end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  !------------------------------------------------------------------------------------!
  ! given V1,V2,V3 orientated triangle, and a distance B, func give the distance from  !
  ! V2 (in the way to V3) such that the point lies at a distance B from V1 and is      !
  ! on the line V3-V2                                                                  !
  !------------------------------------------------------------------------------------!

 real(8) function generaltriangle(V1,V2,V3,B)
 implicit none
  real(8) :: V1(3),V2(3),V3(3),B
  real    :: vV1(3),vV2(3),vV3(3),vA
  real(8) :: chi,A,t1(3),t2(3),C,D
  real(8) :: theta
  C=norme(V3-V1)
  t1=V3-V1
  t2=V3-V2
  chi=ACOS( scalprod(t1,t2)/norme(t1)/norme(t2) )
  vV1=V1
  vV2=V2
  vV3=V3
  call line_exp_point_dist_3d ( vV2(1),vV2(2),vV2(3),vV3(1),vV3(2),vV3(3),vV1(1),vV1(2),vV1(3),vA)
  A=vA
  theta=ASIN(A/B)
  D=(C-cos(theta-chi)*B)/cos(chi)
  D=norme(V3-V2)-D
  generaltriangle=D
 
 return
 end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

      subroutine dir1(costh,sinth,cospsi,sinpsi,rot)
      implicit none
      real(8) costh,sinth,cospsi,sinpsi,cosom,sinom
      real(8) rot(3,3)
      rot(1,1)=costh*cospsi
      rot(1,2)=costh*sinpsi
      rot(1,3)=-sinth
      rot(2,1)=-sinpsi
      rot(2,2)=cospsi
      rot(2,3)=0.0d0
      rot(3,1)=cospsi*sinth
      rot(3,2)=sinpsi*sinth
      rot(3,3)=costh
      return
      end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

      subroutine dir2(cosom,sinom,rot)
      implicit none
       real(8) :: cosom,sinom
       real(8) :: rot(3,3)
       rot(1,1)=cosom
       rot(1,2)=-sinom
       rot(1,3)=0.0d0
       rot(2,1)=sinom
       rot(2,2)=cosom
       rot(2,3)=0.0d0
       rot(3,1)=0.0d0
       rot(3,2)=0.0d0
       rot(3,3)=1.0d0
      return
      end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
 
    !---------------------------------------------!
    ! nsiti    = nbre de points du polygone       !
    ! nf       = nombre de points pour rotations  !
    ! cross1   = number of crossing               !
    ! crossref = nbre crossing demande            !
    !---------------------------------------------!

   subroutine krank(nf2,nn,z1,zn1,cross,nocrossing,error)
   implicit none

      integer              :: nn,idum,nocrossing,i,j,k,k1,nf2,error
      real(8),intent(in)   :: z1(nn,3)
      real(8),intent(out)  :: zn1(nn,3) 
      real(8)              :: po(0:nn-1,3),po1(0:nn-1,3),ro(3,3)
      real(8)              :: cos_the,sin_the,cos_psi,sin_psi,cos_ome,sin_ome
      real(8)              :: ome, rd, som(3),check
      real(8)              :: xx,yy,zz,dx,dy,dz,uu,vv,ww, du,dv,dw,xi,yi
      real(8)              :: aa,bb,cc,a,b,s,t, so,co,coe,delta,coi,soi
      integer              :: nm,nf, fa1, fa2, cross,cross1

      error=0
      cross1=0
34    continue
      nm=int(drand1()*float(nn)+1.)
      nf=int(drand1()*float(nf2)+2.)
      ome=(0.5-drand1())*(pi/2.)       
      if(nf>nn-2) goto 34

      do i=0,nn-1
        do j=1,3
          po(i,j)=z1(mod(nm+i-1,nn)+1,j)-z1(nm,j)
        enddo
      enddo            

      rd=norme(po(nf,:))
      cos_the=po(nf,3)/rd
      sin_the=sqrt(1.0-cos_the**2.)
      if(sin_the<0.0001)goto 33
      cos_psi=po(nf,1)/(rd*sin_the)
      if(po(nf,1)==0.0)goto 33
      sin_psi=po(nf,2)/(rd*sin_the)
      cos_ome=dcos(ome)
      sin_ome=dsin(ome)

      call dir1(cos_the,sin_the,cos_psi,sin_psi,ro) 
      do i=0,nn-1
         do k=1,3
            som(k)=0.
            do k1=1,3
               som(k)=som(k)+ro(k,k1)*po(i,k1)
            enddo
         enddo
         po(i,:)=som(:)
      enddo

      cross=0
      do i=0,nf-1

         xx=po(i,1)
         yy=po(i,2)
         zz=po(i,3)
         dx=po(i+1,1)-po(i,1)
         dy=po(i+1,2)-po(i,2)
         dz=po(i+1,3)-po(i,3)

         do j=nf,nn-1
            if((j.ne.i+1).and.(j-nn+1.ne.i)) then

               uu=po(j,1)
               vv=po(j,2)
               ww=po(j,3)
               du=po(mod(j+1,nn),1)-po(j,1)
               dv=po(mod(j+1,nn),2)-po(j,2)
               dw=po(mod(j+1,nn),3)-po(j,3)

               a=(zz-ww)/dw
               b=dz/dw

               cc=-xx**2.0-yy**2.0+uu**2.0+vv**2.0 &
     &              +(du**2.0+dv**2.0)*a**2.0 +2.0*a*(uu*du+vv*dv)
               bb=-2.0*(xx*dx+yy*dy)+2.0*a*b*(du**2.0+dv**2.0) & 
     &              +2.0*b*(uu*du+vv*dv)
               aa=-(dx**2.0+dy**2.0)+(du**2.0+dv**2.0)*b**2.0

               coe=1.0d0
               do k=1,2

                  coe=-1.0*coe
                  delta=bb**2.0-4.0*aa*cc
                  if (delta.ge.0.0d0) then
                     s=(-bb+coe*sqrt(delta))/2.0/aa
                     t=a+s*b
                  else   
                     goto 98
                  endif
                  
                  delta=(xx+s*dx)**2.0+(yy+s*dy)**2.0

                  co=((uu+t*du)*(xx+s*dx)+(vv+t*dv)*(yy+s*dy))/delta
                  so=((xx+s*dx)*(vv+t*dv)-(uu+t*du)*(yy+s*dy))/delta
                  
                  if ( (s.gt. 0.0d0).and.(s.le.1.0d0) ) then
                  if ( (t.ge.-0.0d0).and.(t.le.1.0d0) ) then
                  if ( (co.le.1.0d0).and.(co.ge.cos_ome)) then
                  if (  so*sin_ome.gt.0.0d0)   then
                     cross1=cross1+1
                     fa1=i
                     fa2=j
                     coi=co
                     soi=so
                     xi=xx+s*dx
                     yi=yy+s*dy
                  endif
                  endif
                  endif
                  endif
                  
 98               continue
               enddo
            endif
         enddo
      enddo

      if(cross1/=0.and.nocrossing==1) then
        !zn1=0
        !zn1=z1
        cross=cross1
        error=1
        return
      endif

      do i=0,nn-1
         do k=1,3
            po1(i,k)=po(i,k)
         enddo
      enddo

      call dir2(cos_ome,sin_ome,ro)
      do i=0,nf
         do k=1,3     
         som(k)=0.
            do k1=1,3
               som(k)=som(k)+ro(k,k1)*po(i,k1)
            enddo
         enddo
         po(i,:)=som(:)
      enddo

 33   continue

      do i=0,nn-1
        zn1(i+1,:)=po(i,:)+z1(nm,:)
      enddo

      cross=cross1
      
    return
    end subroutine


!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

 subroutine dist3D(Aa,Bb,Cc,D_,dP)
 implicit none

 real :: u(3),v(3),w(3),a,b,c,d,e,Dd,sc,sN,sD,tc,tN,tD, &
          & S1(2,3),S2(2,3),dP,Aa(3),Bb(3),Cc(3),D_(3)

    S1(1,:)=Aa
    S1(2,:)=Bb
    S2(1,:)=Cc 
    S2(2,:)=D_
    u=S1(2,:)-S1(1,:)
    v=S2(2,:)-S2(1,:)
    w=S1(1,:)-S2(1,:)
    a = scalprod(u,u)      
    b = scalprod(u,v)
    c = scalprod(v,v)     
    d = scalprod(u,w)
    e = scalprod(v,w)
    Dd = a*c - b*b      
    sD = Dd    
    tD = Dd     
    if (D < 1.d-6) then
        sN = 0.0
        sD = 1.0
        tN = e
        tD = c
    else 
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if (sN < 0.0) then
            sN = 0.0
            tN = e
            tD = c
        elseif (sN > sD) then
            sN = sD
            tN = e + b
            tD = c
        endif
    endif
    if (tN < 0.0) then
        tN = 0.0
        if (-d < 0.0) then
            sN = 0.0
        elseif (-d > a) then
            sN = sD
        else 
            sN = -d
            sD = a
        endif
    elseif (tN > tD) then
        tN = tD
        if ((-d + b) < 0.0) then
            sN = 0
        elseif ((-d + b) > a) then
            sN = sD
        else 
            sN = (-d + b)
            sD = a
        endif
    endif
    sc = sN / sD
    tc = tN / tD
    dP = norme((w + (sc*u) - (tc*v)))

   return 
   end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  function PointInTriangle2(p,a,b,c,eps2)
  implicit none
   real(8),dimension(3),intent(in) :: p,a,b,c
   real(8),intent(in)              :: eps2
   logical                         :: PointInTriangle2

    if(SameSide(p,a,b,c,eps2).and.  &
    &  SameSide(p,b,a,c,eps2).and.  &
    &  SameSide(p,c,a,b,eps2)) then
     PointInTriangle2= .true.
    else
     PointInTriangle2= .false.
    endif
    if(norme(p-a)<1.d-9) PointInTriangle2= .false.
    if(norme(p-b)<1.d-9) PointInTriangle2= .false.
    if(norme(p-c)<1.d-9) PointInTriangle2= .false.

  return
  end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  function PointInTriangle(kk,ll,mm,p,a,b,c,eps2)
  implicit none
  real(8),dimension(3),intent(in) :: p,a,b,c
  real(8),intent(in)              :: eps2
  integer                         :: kk,ll,mm
  logical                         :: PointInTriangle

    if(SameSide(p,a,b,c,dble(kk)*eps2).and. &
    &  SameSide(p,b,a,c,dble(ll)*eps2).and. &
    &  SameSide(p,c,a,b,dble(mm)*eps2)) then
     PointInTriangle= .true.
    else
     PointInTriangle= .false.
    endif

  return
  end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  function perp(v)
  implicit none
  real(8),intent(in) :: v(3)
  real(8)            :: perp(3)
    perp(1)= v(2)
    perp(2)=-v(1)
    perp(3)= 0.d0
  end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

    function SameSide(p1,p2,a,b,eps2)
    implicit none
    real(8),dimension(3),intent(in) :: p1,p2,a,b
    real(8),intent(in)              :: eps2
    real(8)                         :: cp1(3),cp2(3)
    logical                        :: SameSide

    cp1 = vecprod(b-a, p1-a)
    cp2 = vecprod(b-a, p2-a)
    if(scalprod(cp1,cp2)>=eps2) then
     SameSide= .true.
    else
     SameSide= .false.
    endif

   return
   end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

   subroutine knotjoin(tab_knot,nInter,points,chaos,n_pas,npasreel,ok)
   implicit none
   integer                                   :: temp,temp2,nx,ny,nz,nmax,n
   integer                                   :: i,j,k,compteur,tab,u,l,boucle,total,delta1,delta2,lastcount
   integer,intent(in)                        :: nInter,n_pas
   integer,intent(out)                       :: npasreel
   logical                                   :: correct,test,test2,test3
   logical, intent(out)                      :: ok
   integer, dimension(6)                     :: direc
   integer, dimension (n_pas,3),intent(out)  :: tab_knot
   real(8),dimension(:,:), allocatable       :: deplacement,deplacement_temp,deplacement_temp_temp
   real(8), intent(in)                       :: chaos
   real(8), dimension (3)                    :: vect_temp, pos_temp,tab_mov,A,B
   real(8), dimension (nInter,3), intent(in) ::  points
   real(8)                                   :: z,sum,norm,prod,hasard,long, normebefore,normeafter
   
   total=0
   l=0
   lastcount=0
   ok=.true.
   
   do while (l<=nInter-1)
   l=l+1
   test2=.true.
   u=0
   if(lastcount>5) then
   ok=.false.
   exit
   endif
   if(l<nInter) then
   A(:)=points(l,:)
   B(:)=points(l+1,:)
   else
   A(:)=points(nInter,:)
   B(:)=points(1,:)
   endif
   
   nmax=0
   do i=1,3
   nmax=nmax+ABS(A(i)-B(i))
   enddo
   nmax=nmax*6
   boucle=0
   do while(test2)
   u=0
   do i = total+1, n_pas
   tab_knot(i,:) = vnull(:)
   enddo
   
   do i=1,total-1
   
   vect_temp=dble(tab_knot(i,:)-tab_knot(i+1,:))
   hasard= norme(vect_temp)
   
   
   if(hasard/=1.) then
   write(*,*) 'PROG PRINCIPAL: ce n est pas un noeud'
   write(*,*) 'Au point a relier numero : ', l
   write(*,*) 'Appuyer sur une touche'
   stop
   endif
   enddo
   
   
   tab_knot(total+1,:)=Int(A(:))
   tab_mov(:)=vnull(:)
   compteur = 1
   k=0
   temp2=0
   test=.true.
   boucle=boucle+1
   if(boucle>10) then
   test2=.false.
   l=0
   total=0
   lastcount=lastcount+1
   endif
   
   if(test2) then
   do while (test) 
   norm=drand1()
   tab_mov(:) = vnull
   if(compteur>temp2) then
   do i=1,6
    direc(i)=0
   enddo
   sum=0  
   else
   sum=0
   do i=1,6
    sum=sum+direc(i)
   enddo
   
   if(sum==6) then
   k=k+1
   do i = total+2, n_pas
   tab_knot(i,:) = vnull(:)
   enddo
   compteur=1
   temp2=0
   sum=0
   do i=1,6
   direc(i)=0
   enddo
   endif
   endif
   
   if ( norm < 1./6. ) then
   tab_mov(1) =1.
   if (compteur == temp2) direc(1)=1
   endif
   if ( norm >= 1./6. .AND. norm < 2./6.) then
   tab_mov(2) =1.
   if (compteur == temp2) direc(2)=1
   endif
   if ( norm >= 2./6. .AND. norm < 3./6.) then
   tab_mov(3) =1.
   if (compteur == temp2) direc(3)=1
   endif
   if ( norm >= 3./6. .AND. norm < 4./6.) then
   tab_mov(3) = - 1.
   if (compteur == temp2) direc(4)=1
   endif
   if ( norm >= 4./6. .AND. norm < 5./6.) then
   tab_mov(2) = - 1.
   if (compteur == temp2) direc(5)=1
   endif
   if ( norm >= 5./6. .AND. norm < 1.) then
   tab_mov(1) = - 1.
   if (compteur == temp2) direc(6)=1
   endif
   
   if( total+compteur >= n_pas-2) then
   if(mod(boucle,10)==0 .AND. mod(lastcount,4)==0) write(*,*) &
   & 'DERNIER POINT REUSSI: ',l,' /EN TOUT: ', nInter+1
   test=.false.
   endif
   
   vect_temp(:) = tab_mov(:) + tab_knot(total+compteur,:)
   correct = .true.
   
   do i=1,nInter
   if(vect_temp(1)==points(i,1) .AND. vect_temp(2)==points(i,2) &
    & .AND. vect_temp(3)==points(i,3)) correct=.false.
   enddo
   
   if(vect_temp(1)==B(1).AND.vect_temp(2)==B(2) &
    & .AND.vect_temp(3)==B(3))then
   test=.false.
   test2=.false.
   total=total+compteur
   endif
   
   if(test) then
   if(total+compteur-1 >1 )then
    do j = 1, total+compteur-1 
   temp = 0
   do i=1,3
   if (vect_temp(i) == tab_knot(j,i)) temp=temp+1
    if (temp==3) correct = .false.
   enddo
   enddo
   endif
   
   normebefore=norme(B(:)-tab_knot(total+compteur,:))
   normeafter= norme(B(:)-vect_temp)
   hasard=drand1()
   if(normeafter>normebefore .AND. hasard>chaos) correct=.false.
   
   temp2=compteur  
   
   if(correct) then
   compteur = compteur + 1
   tab_knot(total+compteur,:) = Int(vect_temp)
   endif
   endif
   enddo
   endif
   enddo
   enddo
   
   npasreel=total
   
   end subroutine
   
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
     
   subroutine cross2d(n_pas,tabx,taby,tabtot_x,tabtot_y,nInter,err,Alexander,signes,Inter)
   implicit none
   integer,intent(in)                            ::  n_pas,nInter
   real(8), dimension(3)                         ::  v1,v2,v3
   real(8), dimension(n_pas,2,2)                 ::  segmentproj
   real(8), dimension(n_pas)  ,intent(in)        ::  tabx,taby
   real(8), dimension(:,:),allocatable           ::  Idisorder,Iorder
   real(8)                                       ::  x,y
   real(8), dimension(2)                         ::  A
   real(8), dimension(nInter,2),intent(out)      ::  Inter
   real(8), intent(in)                           ::  err
   real(8), dimension (n_pas+nInter),intent(out) ::  tabtot_x,tabtot_y
   integer  , dimension(nInter),intent(in)       ::  signes
   integer  , dimension(nInter)                  ::  internum
   integer  , dimension(nInter,2)                ::  Inter2
   integer  , dimension(nInter/2)                ::  Interarc
   integer  , dimension(nInter/2,2)              ::  arc
   integer                                       ::  i,j,l,k,t,u,v,sign,b,Epsilon
   integer  , dimension(:), allocatable          ::  Jorder, Jdisorder
   integer  , dimension(n_pas)                   ::  numbcross
   logical                                       ::  crossing
   character(4), dimension(nInter/2,nInter/2)    ::  Alexandertot
   character(4), dimension(nInter/2-1,nInter/2-1),intent(out) ::  Alexander
   
   do i=1,n_pas
   segmentproj(i,1,1)=tabx(i)
   segmentproj(i,1,2)=taby(i)
   if(i<n_pas)then
   segmentproj(i,2,1)=tabx(i+1)
   segmentproj(i,2,2)=taby(i+1)
   else
   segmentproj(i,2,1)=tabx(1)
   segmentproj(i,2,2)=taby(1)
   endif
   enddo
    
   do i=1,n_pas
    numbcross(i)=0
   enddo
   
   l=0
   do i=1,n_pas
   k=0
   do j=1,n_pas
   crossing=.FALSE.
   call intersect(segmentproj(i,:,:),segmentproj(j,:,:),crossing,x,y,err)
   if (crossing) then 
   k=k+1
   endif
   enddo
   numbcross(i)=k
   l=l+k
   enddo
     
   if(l/=nInter) write(23,*) 'Routine cross2d: TEST nInter ERROR'
   
   l=0
   u=0
   v=0
   
   do i=1,n_pas
   
   allocate(Idisorder(numbcross(i),2),Iorder(numbcross(i),2),&
   & Jdisorder(numbcross(i)),Jorder(numbcross(i)))
    
   v=v+1
   
   tabtot_x(v)=tabx(i)
   tabtot_y(v)=taby(i)
   
   if (numbcross(i)>1) then 
   k=0
   do j=1,n_pas
   crossing=.FALSE.
   call intersect(segmentproj(i,:,:),segmentproj(j,:,:),crossing,x,y, err)
   if (crossing) then
   k=k+1
   Idisorder(k,1)=x
   Idisorder(k,2)=y
   Jdisorder(k)=j
   endif
   enddo
   A(1)=segmentproj(i,1,1)
   A(2)=segmentproj(i,1,2)
   call ordertab(2,A,Idisorder,Iorder, Jdisorder, Jorder, numbcross(i))
   do t=1,k
   Inter(u+t,:)=Iorder(t,:)
   Inter2(u+t,1)=i
   Inter2(u+t,2)=Jorder(t)
   tabtot_x(v+t)=Iorder(t,1)
   tabtot_y(v+t)=Iorder(t,2)
   enddo
   u=u+k
   v=v+k
   endif
   if (numbcross(i) == 1 ) then
   do j=1,n_pas
   crossing=.FALSE.
   call intersect(segmentproj(i,:,:),segmentproj(j,:,:),crossing,x,y,err)
   if (crossing) then
   Iorder(1,1)=x
   Iorder(1,2)=y
   Jorder(1)=j
   Inter2(u+1,1)=i
   Inter2(u+1,2)=j
   Inter(u+1,:)=Iorder(1,:)
   tabtot_x(v+1)=Iorder(1,1)
   tabtot_y(v+1)=Iorder(1,2)
   endif
   enddo
   u=u+1
   v=v+1
   endif
   
   deallocate(Idisorder,Iorder,Jdisorder,Jorder)
    
   enddo
   
   do i=1,nInter
   internum(i)=i
   enddo
   
   k=1
   do i=1,nInter
   do j=i+1,nInter
   if(ABS(Inter(j,1)-Inter(i,1))<err .AND. ABS(Inter(j,2)-Inter(i,2))<err ) then
   internum(j)=k
   internum(i)=k
   if(signes(i)/=-signes(j)) write(23,*) &
     &'ERROR routines cross,signes pas coherents'
   k=k+1
   endif
   enddo
   enddo
   
   
   k=1
   
   do i=1,nInter 
   if(signes(i)==-1 .AND. i<nInter) then
   j=i
   do while(j<nInter)
   j=j+1
   if( signes(j)==-1 ) then
   k=k+1
   arc(k,1)=i
   arc(k,2)=j
   j=nInter+1
   endif
   enddo 
   endif
   if(signes(i)==1 ) Interarc(internum(i))=k
   enddo
   
   
   if(signes(nInter)==1) then
   do i=1,nInter
   if(signes(i)==-1) then
   k=0
   else
   k=k+1
   endif
   enddo
   endif
   
   if(signes(nInter)==-1) k=0
   
   if(k>0) then
   do i=1,k
    Interarc(internum(nInter-i+1))=1
   enddo
   endif
   
   arc(1,1)=arc(nInter/2,2)
   arc(1,2)=arc(2,1)  
   
   do i=1,nInter/2
   do j=1,2
   arc(i,j) = internum(arc(i,j))
   enddo
   enddo
   
   
   do i=1,nInter/2
   do j=1,nInter/2
   Alexandertot(i,j)='0'
   enddo
   enddo
   Epsilon=0
   l=0
   do i=1,nInter
   if(signes(i)==1) then
   l=l+1
   k=0
   do j=1,nInter/2
   if (arc(j,1)==internum(i)) then
   t=j
   k=k+1
   endif
   if (arc(j,2)==internum(i)) then
   u=j
   k=k+1
   endif
   
   
   if(k>2) write(23,*)'ERREUR: routine cross'
   enddo
   
   v1(1)=segmentproj(Inter2(i,1),2,1) - segmentproj(Inter2(i,1),1,1) 
   v1(2)=segmentproj(Inter2(i,1),2,2) - segmentproj(Inter2(i,1),1,2)
   v1(3)=0
   v2(1)=segmentproj(Inter2(i,2),2,1) - segmentproj(Inter2(i,2),1,1) 
   v2(2)=segmentproj(Inter2(i,2),2,2) - segmentproj(Inter2(i,2),1,2)
   v2(3)=0
   v3=vecprod(v1,v2)
   if(v3(3)==0.) write(23,*) 'Error cross2d: right left'
   if(v3(3)<0.) then
   b=u
   u=t
   t=b
   Epsilon=Epsilon-1
   else
   Epsilon=Epsilon+1
   endif
   b=Interarc(internum(i))
   
   if(b/=u .AND. b/=t)  Alexandertot(internum(i),b)= '1-t'
   if(b==u .AND. b/=t)  Alexandertot(internum(i),b)= '-t'
   if(b==t .AND. b/=u)  Alexandertot(internum(i),b)= '1'
   if(b==t .AND. b==u)  Alexandertot(internum(i),b)= '0'
   if(u==t .AND. u/=b)  Alexandertot(internum(i),u)= '-1+t'
   if(u/=t .AND. u/=b)  Alexandertot(internum(i),u)= '-1'  
   if(u/=t .AND. t/=b)  Alexandertot(internum(i),t)= 't' 
    
   endif
   enddo
   
   if(nInter/2>1) then
   do i=1,nInter/2-1
   do j=1,nInter/2-1
    Alexander(i,j)=Alexandertot(i,j)
   enddo
   enddo
   else
   write(23,*) 'Polynome d Alexander trivial: P=1' 
   endif
   
   end subroutine cross2d
   
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
     
   subroutine cross3d(n_pas,vect_1,vect_2,vect_3,tab_knot_2d_x,&
       & tab_knot_2d_y,tab_knot,tabtot_x,tabtot_y,nInter,err,Alexander)
   
   implicit none

   integer,intent(in)                                           ::  n_pas,nInter
   real(8), dimension (3,3)                                     ::  P
   real(8), dimension (3),intent(in)                            ::  vect_2,vect_3,vect_1
   real(8), dimension (3)                                       ::  v1,v2,v3
   real(8), dimension (n_pas,2,3)                               ::  segment
   real(8), dimension (n_pas,2,2)                               ::  segmentproj
   integer, dimension (n_pas,3),intent(in)                      ::  tab_knot
   real(8), dimension (n_pas)  ,intent(in)                      ::  tab_knot_2d_x,tab_knot_2d_y
   real(8), dimension (n_pas+nInter),intent(out)                ::  tabtot_x,tabtot_y
   real(8), dimension (:,:),allocatable                         ::  Idisorder,Iorder
   real(8)                                                      ::  x,y
   real(8), dimension(3)                                        ::  A,Iplan
   real(8), dimension(nInter,3)                                 ::  Inter
   real(8), intent(in)                                          ::  err
   integer  , dimension(nInter)                                 ::  signes
   integer  , dimension(nInter)                                 ::  internum
   integer  , dimension(nInter,2)                               ::  Inter2
   integer  , dimension(nInter/2)                               ::  Interarc
   integer  , dimension(nInter/2,2)                             ::  arc
   integer                                                      ::  i,j,l,k,t,u,v,sign,b,Epsilon
   integer  , dimension(:), allocatable                         ::  Jorder, Jdisorder
   integer  , dimension(n_pas)                                  ::  numbcross
   logical                                                      ::  crossing
   character(4), dimension(nInter/2,nInter/2)                   ::  Alexandertot
   character(4), dimension(nInter/2-1,nInter/2-1),intent(out)   ::  Alexander
   
   P(:,1)=vect_1(:)
   P(:,2)=vect_2(:)
   P(:,3)=vect_3(:)
   
   do i=1,n_pas
   segmentproj(i,1,1)=tab_knot_2d_x(i)
   segmentproj(i,1,2)=tab_knot_2d_y(i)
   if(i<n_pas)then
   segmentproj(i,2,1)=tab_knot_2d_x(i+1)
   segmentproj(i,2,2)=tab_knot_2d_y(i+1)
   else
   segmentproj(i,2,1)=tab_knot_2d_x(1)
   segmentproj(i,2,2)=tab_knot_2d_y(1)
   endif
   enddo
    
   do i=1,n_pas
   segment(i,1,:)=dble(tab_knot(i,:))
   if(i<n_pas)then
   segment(i,2,:)=dble(tab_knot(i+1,:))
   else
   segment(i,2,:)=dble(tab_knot(1,:))
   endif
   enddo
   
   do i=1,n_pas
   k=0
   do j=1,n_pas
   crossing=.FALSE.
   call intersect(segmentproj(i,:,:),segmentproj(j,:,:),crossing,x,y,err)
   if (crossing) then 
   k=k+1
   endif
   enddo
   numbcross(i)=k
   enddo
     
   l=0
   u=0
   v=0
   
   do i=1,n_pas
   
   allocate(Idisorder(numbcross(i),3),Iorder(numbcross(i),3),&
   & Jdisorder(numbcross(i)),Jorder(numbcross(i)))
    
   v=v+1
   
   tabtot_x(v)=tab_knot_2d_x(i)
   tabtot_y(v)=tab_knot_2d_y(i)
   
   if (numbcross(i)>1) then 
   k=0
   do j=1,n_pas
   crossing=.FALSE.
   call intersect(segmentproj(i,:,:),segmentproj(j,:,:),crossing,x,y,err)
   if (crossing) then
   k=k+1
   Idisorder(k,1)=x
   Idisorder(k,2)=y
   Idisorder(k,3)=0.
   Jdisorder(k)=j
   endif
   enddo
   A(1)=segmentproj(i,1,1)
   A(2)=segmentproj(i,1,2)
   A(3)=0.
   call ordertab(3,A,Idisorder,Iorder, Jdisorder, Jorder, numbcross(i))
   do t=1,k
   Inter(u+t,:)=Iorder(t,:)
   Inter2(u+t,1)=i
   Inter2(u+t,2)=Jorder(t)
   tabtot_x(v+t)=Iorder(t,1)
   tabtot_y(v+t)=Iorder(t,2)
   enddo
   u=u+k
   v=v+k
   endif
   if (numbcross(i) == 1 ) then
   do j=1,n_pas
   crossing=.FALSE.
   call intersect(segmentproj(i,:,:),segmentproj(j,:,:),crossing,x,y,err)
   if (crossing) then
   Iorder(1,1)=x
   Iorder(1,2)=y
   Iorder(1,3)=0.
   Jorder(1)=j
   Inter2(u+1,1)=i
   Inter2(u+1,2)=j
   Inter(u+1,:)=Iorder(1,:)
   tabtot_x(v+1)=Iorder(1,1)
   tabtot_y(v+1)=Iorder(1,2)
   endif
   enddo
   u=u+1
   v=v+1
   endif
   
   do j=1,numbcross(i)
    Iplan(:)=MATMUL(P,Iorder(j,:))
    call updown(segment(i,:,:),segment(Jorder(j),:,:),Iplan,vect_3,err,sign)
    l=l+1
    signes(l)=sign
   enddo
   deallocate(Idisorder,Iorder,Jdisorder,Jorder)
   enddo
    
   do i=1,nInter
   internum(i)=i
   enddo
   
   k=1
   do i=1,nInter
   do j=i+1,nInter
   if(ABS(Inter(j,1)-Inter(i,1))<err .AND. ABS(Inter(j,2)-Inter(i,2))<err &
   & .AND. ABS(Inter(j,3)-Inter(i,3))<err) then
   internum(j)=k
   internum(i)=k
   if(signes(i)/=-signes(j)) write(23,*) &
   & 'ERROR routines cross, signes pas coherents'
   k=k+1
   endif
   enddo
   enddo
   
   k=1
   do i=1,nInter 
   if(signes(i)==-1 .AND. i<nInter) then
   j=i
   do while(j<nInter)
   j=j+1
   if( signes(j)==-1 ) then
   k=k+1
   arc(k,1)=i
   arc(k,2)=j
   j=nInter+1
   endif
   enddo 
   endif
   if(signes(i)==1 ) Interarc(internum(i))=k
   enddo
   if(signes(nInter)==1) then
   do i=1,nInter
   if(signes(i)==-1) then
   k=0
   else
   k=k+1
   endif
   enddo
   endif
   
   if(signes(nInter)==-1) k=0
   
   if(k>0) then
   do i=1,k
    Interarc(internum(nInter-i+1))=1
   enddo
   endif
   
   arc(1,1)=arc(nInter/2,2)
   arc(1,2)=arc(2,1) 
   
   do i=1,nInter/2
   do j=1,2
   arc(i,j) = internum(arc(i,j))
   enddo
   enddo
   
   do i=1,nInter/2
   do j=1,nInter/2
   Alexandertot(i,j)='0'
   enddo
   enddo
   
   Epsilon=0
   l=0
   do i=1,nInter
   if(signes(i)==1) then
   l=l+1
   k=0
   do j=1,nInter/2
   if (arc(j,1)==internum(i)) then
   t=j
   k=k+1
   endif
   if (arc(j,2)==internum(i)) then
   u=j
   k=k+1
   endif
   if(k>2) write(23,*) 'ERREUR: routine cross'
   enddo
   v1(1)=segmentproj(Inter2(i,1),2,1) - segmentproj(Inter2(i,1),1,1) 
   v1(2)=segmentproj(Inter2(i,1),2,2) - segmentproj(Inter2(i,1),1,2)
   v1(3)=0
   v2(1)=segmentproj(Inter2(i,2),2,1) - segmentproj(Inter2(i,2),1,1) 
   v2(2)=segmentproj(Inter2(i,2),2,2) - segmentproj(Inter2(i,2),1,2)
   v2(3)=0
   v3=vecprod(v1,v2)
   if(v3(3)<0.) then
   b=u
   u=t
   t=b
   Epsilon=Epsilon+1
   else
   Epsilon=Epsilon-1
   endif
   b=Interarc(internum(i))
   if(b/=u .AND. b/=t)  Alexandertot(internum(i),b)='1-t'
   if(b==u .AND. b/=t)  Alexandertot(internum(i),b)= '-t'
   if(b==t .AND. b/=u)  Alexandertot(internum(i),b)= '1'
   if(b==t .AND. b==u)  Alexandertot(internum(i),b)= '0'
   if(u==t .AND. u/=b)  Alexandertot(internum(i),u)= '-1+t'
   if(u/=t .AND. u/=b)  Alexandertot(internum(i),u)= '-1'  
   if(u/=t .AND. t/=b)  Alexandertot(internum(i),t)= 't' 
   endif
   enddo
    
   if(nInter/2>1) then
   do i=1,nInter/2-1
   do j=1,nInter/2-1
     Alexander(i,j)=Alexandertot(i,j)
   enddo
   enddo
   else
   write(23,*) 'Polynome d Alexander trivial: P=1' 
   endif  
     
   end subroutine 
  
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
     
   subroutine updown(sgm1, sgm2, Iplan, e3, err, sign)
   implicit none
   integer, intent(out)                 ::  sign
   real(8), dimension(2,3), intent(in)  ::  sgm1,sgm2
   real(8), dimension(3), intent(in)    ::  Iplan,e3
   real(8), intent(in)                  ::  err
   real(8)                              ::  k1,k2,temp
   integer                              ::  k
   integer,dimension(2)                 ::  notnullv
   integer,dimension(2,2)               ::  nullv
   real(8), dimension(3)                ::  vec,test
   real(8), dimension(2,3)              ::  Itab,vect,Itabverify
   real(8), dimension(2,2,3)            ::  sgm
   real(8), dimension(2)                ::  Icompare
   integer                              ::  i,j
   
   sgm(1,:,:)=sgm1
   sgm(2,:,:)=sgm2
   
   do i=1,2
   vect(i,:)=sgm(i,2,:) - sgm(i,1,:)
   enddo
    
   temp=scalprod(Iplan,e3)
   if(temp>err) write(23,*) &
   & 'ERREUR1: routine updown, Iplan n est pas dans le plan de proj'
   
   do j=1,2
   do i=1,3
     if (ABS(vect(j,i)) >0) notnullv(j)=i
   enddo
   enddo
   
    do j=1,2
   k=0
   do i=1,3
    if (i /= notnullv(j)) then
   k=k+1
   nullv(j,k)=i
   endif
   enddo
    enddo
   
   do i=1,2
   
   vec(:)=Iplan(:)-sgm(i,1,:)
    
   if(ABS(e3(nullv(i,1)))  > 0) then
   k1= -vec(nullv(i,1)) / e3(nullv(i,1))
   else
   if( e3(nullv(i,2)) == 0 ) write(23,*) &
   &'ERREUR2: routine updown, e3 pas adapte, 2 composantes nulles'
   k1= -vec(nullv(i,2)) / e3(nullv(i,2))
   endif
   
   Itab(i,:)=Iplan(:) + k1*e3(:)
   
   temp=scalprod(Itab(i,:),e3)
   Icompare(i)=temp
   
    k2=(Itab(i,notnullv(i))-sgm(i,1,notnullv(i)))/vect(i,notnullv(i))
    Itabverify(i,:)=sgm(i,1,:) + k2*vect(i,:)
   
   if(ABS(Itabverify(i,1)-Itab(i,1)) > err .or. &
   & ABS(Itabverify(i,2)-Itab(i,2)) > err .or. &
   & ABS(Itabverify(i,3)-Itab(i,3)) > err) then
   write(23,*) 'EREUR3: routine updown, err certainement trop petit'
    endif
   
    enddo
   
    do i=1,2
     test(:)= Itab(i,:) - sgm(i,1,:)
     call vecnorm(test,test)
     vec=vecprod(test,vect(i,:))
     temp=norme(vec)
   
     if (temp>err) write(23,*) &
     & 'ERROR4: routine updown, err certainement trop petit'
   
     test(:)= Itab(i,:) - Iplan(:)
     call vecnorm(test,test)
     vec=vecprod(test,e3)
     temp=norme(vec)
     if (temp>err) write(23,*) &
     & 'ERROR5: routine updown, err certainement trop petit'
     enddo
   
    If (Icompare(1)>Icompare(2)) then 
     sign=1
    else
     sign=-1
    endif
   
    If (Icompare(1)==Icompare(2)) then
    sign=0
    write(23,*) 'ERROR6: routine updown, ne passe ni dessus ni dessous'
    endif
 end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

subroutine intersect(sgm1,sgm2,crossing,x,y,err)
implicit none
  real(8), intent(out)               ::  x, y
  real(8), intent(in)                ::  err
  real(8), dimension(2)              ::  vect_dir
  integer                           ::  c,d,z,e,f,vec,up1,up2
  real(8), dimension(2)              ::  a,b
  real(8)                            ::  k
  real(8), dimension(2,2),intent(in) ::  sgm1,sgm2
  real(8), dimension(2,2,2)          ::  sgm
  logical,intent(out)               ::  crossing
    
    c=0.; x=0.; y=0.
    
     if(sgm1(1,1)>sgm1(2,1)) then
      sgm(1,1,:)=sgm1(2,:)
      sgm(1,2,:)=sgm1(1,:)
     else
      sgm(1,1,:)=sgm1(1,:)
      sgm(1,2,:)=sgm1(2,:)
     endif

     if(sgm2(1,1)>sgm2(2,1)) then
      sgm(2,1,:)=sgm2(2,:)
      sgm(2,2,:)=sgm2(1,:)
     else
      sgm(2,1,:)=sgm2(1,:)
      sgm(2,2,:)=sgm2(2,:)
     endif
    
    up1=-1; up2=-1
    if(sgm1(1,2)<sgm1(2,2))           up1=1
    if(DABS(sgm1(1,2)-sgm1(2,2))<err) up1=0
    if(sgm2(1,2)<sgm2(2,2))           up2=1
    if(DABS(sgm2(1,2)-sgm2(2,2))<err) up2=0
     
    do e = 1, 2
    do f = 1, 2
    if (sgm(1,e,1) == sgm(2,f,1) .AND. sgm(1,e,2) == sgm(2,f,2)) then
    c = c + 1
    crossing=.False.
    endif
    enddo
    enddo
    
    if (c/=0) crossing=.FALSE.
    
    if (c==0) then  
    d = 0
    z = 0
    vec=0
    do e=1,2
    vect_dir(:) = sgm(e,2,:) - sgm(e,1,:)
    if ( ABS(vect_dir(2)) < err ) then
    b(e) = sgm(e,1,2)
    d=d+1
    endif
    if ( ABS(vect_dir(1)) < err ) then
    z = z + 1
    x = sgm(e,1,1)
    vec=e-1
    vec=(1+(-1)**(vec))/2 +1
    a(vec)=(sgm(vec,2,2)-sgm(vec,1,2))/(sgm(vec,2,1)-sgm(vec,1,1))
    endif
    if ( ABS(vect_dir(2)) >= err .AND. ABS(vect_dir(1)) >= err ) then
    k = - ( sgm(e,1,1) / vect_dir(1) )
    b(e) = sgm(e,1,2) + k * vect_dir(2)
    endif
    enddo
    
    if(d==2) crossing =.FALSE.
    if(z==2) crossing =.FALSE.
    
    if( d<2 .AND. z==0 ) then
    do e=1,2
    a(e) = ( sgm(e,2,2) - sgm(e,1,2) ) / ( sgm(e,2,1) - sgm(e,1,1) )
    enddo
    if(DABS(a(1)-a(2))>0) then
    x = (b(2)-b(1))/(a(1)-a(2))
    y = a(1)*x + b(1)
    crossing=.true.
    else
    crossing=.false.
    endif
     endif
    
    if( z==0 .AND. d<2) then
    if( x>=sgm(1,1,1) .AND. x<=sgm(1,2,1) .AND. &
      & x>=sgm(2,1,1) .AND. x<=sgm(2,2,1) ) then 
     crossing = .TRUE.
    else
    crossing = .FALSE.
    endif
    endif
     
    if( z == 1 .AND. d < 2 ) then
    if( x >= sgm(vec,1,1) .AND. x <= sgm(vec,2,1) ) then 
    crossing=.TRUE.
    y = a(vec)*x + b(vec)
    else
    crossing= .FALSE. 
    endif
    endif
    
    do f=1,2
    do e=1,2
    if(ABS(x-sgm(f,e,1) )<=err .and. ABS( y-sgm(f,e,2) )<=err) crossing=.FALSE.
    enddo
    enddo
    endif  
    
    
    if(crossing .AND. d==0) then
    if(up1== 1 )then
    if(y<sgm1(1,2).or.y>sgm1(2,2)) crossing=.false.
     endif
    if(up1==-1 )then
    if(y>sgm1(1,2).or.y<sgm1(2,2)) crossing=.false.
    endif
     if(up2== 1 )then
    if(y<sgm2(1,2).or.y>sgm2(2,2)) crossing=.false.
     endif
    if(up2==-1 )then
    if(y>sgm2(1,2).or.y<sgm2(2,2)) crossing=.false.
    endif
    endif 
    
    if(crossing .AND. d==1) then
    if(up1==0 .AND. up2== 1)then
     if(y<sgm2(1,2).or.y>sgm2(2,2)) crossing=.false.   
    endif
    if(up1==0 .AND. up2==-1)then
     if(y>sgm2(1,2).or.y<sgm2(2,2)) crossing=.false.   
    endif
    if(up2==0 .AND. up1== 1)then
     if(y<sgm1(1,2).or.y>sgm1(2,2)) crossing=.false.   
    endif
    if(up2==0 .AND. up1==-1)then
     if(y>sgm1(1,2).or.y<sgm1(2,2)) crossing=.false.   
    endif
    endif
    
return
end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

subroutine crosscount(n_pas,nInter,tab_knot_2d_x,tab_knot_2d_y,err)
implicit none
integer, intent(in)                 ::  n_pas
real(8),intent(in)                  ::  err
real(8)                             ::  x,y
real(8),dimension(n_pas,2,3)        ::  segment
real(8),dimension(n_pas,2,2)        ::  segmentproj
real(8),dimension(n_pas),intent(in) ::  tab_knot_2d_x,tab_knot_2d_y
integer                             ::  i,j
logical                             ::  crossing
integer,intent(out)                 ::  nInter

  do i=1,n_pas
  segmentproj(i,1,1)=tab_knot_2d_x(i)
  segmentproj(i,1,2)=tab_knot_2d_y(i)
  if(i<n_pas)then
    segmentproj(i,2,1)=tab_knot_2d_x(i+1)
    segmentproj(i,2,2)=tab_knot_2d_y(i+1)
  else
    segmentproj(i,2,1)=tab_knot_2d_x(1)
    segmentproj(i,2,2)=tab_knot_2d_y(1)
  endif
 enddo
 
 nInter=0
 do i=1,n_pas
  do j=1,n_pas
   crossing=.FALSE.
   call intersect(segmentproj(i,:,:),segmentproj(j,:,:),crossing,x,y, err)
   if (crossing) nInter=nInter+1
  enddo
 enddo

end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

subroutine basegen(mat_plan, vect_1,vect_2,vect_3 )
implicit none
 real(8),dimension(3,2),intent(in) :: mat_plan
 real(8)                           :: sign
 real(8),dimension(3),intent(out)  :: vect_1,vect_2,vect_3
 real(8),dimension(3)              :: temp

 temp(:)=vnull(:)
 call vecnorm(mat_plan(:,1),vect_1)
 call vecnorm(mat_plan(:,2),vect_2)
 vect_3=vecprod(vect_1,vect_2)
 call vecnorm(vect_3,vect_3)
 vect_2=vecprod(vect_3,vect_1)
 temp=vecprod(vect_1,vect_2)
 sign=scalprod(vect_3,temp)

 if (sign >=1.+1d-10 .or.sign<= 1.-1d-10) write (23,*) &
                  & 'routine basegen: base non directe !'
end subroutine basegen

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
 
subroutine project(tab_knot, n_pas, vect_1, vect_2, vect_3, x, y)
implicit none
integer, intent(in)                     ::  n_pas
integer                                 ::  i
real(8), dimension(n_pas),intent(out)   ::  x,y 
real(8), dimension(n_pas,3),intent(in)  ::  tab_knot
real(8), dimension (2,2)                ::  mat_temp1,mat_temp2
real(8), dimension (2,3)                ::  mat_plan_trans,mat_temp3
real(8), dimension (n_pas,3)            ::  knot_proj
real(8)                                 ::  det,temp,sign
real(8), dimension(3)                   ::  temp2
real(8), dimension(3),intent(in)        ::  vect_1,vect_2,vect_3
real(8), dimension(3,3)                 ::  mat_proj
real(8), dimension(3,2)                 ::  mat_plan

 mat_plan(:,1)  = vect_1(:)
 mat_plan(:,2)  = vect_2(:)
 mat_plan_trans = TRANSPOSE (mat_plan)
 mat_temp1      = MATMUL (mat_plan_trans, mat_plan)

 det = mat_temp1(1, 1)*mat_temp1(2, 2) - mat_temp1(1, 2)*mat_temp1(2, 1)

 if(det == 0) write(23,*) 'routine project: vecteurs du plan colineaires'

 mat_temp2(1, 1) = ( 1.d0 / det) * mat_temp1(2, 2)
 mat_temp2(2, 2) = ( 1.d0 / det) * mat_temp1(1, 1)
 mat_temp2(1, 2) = (-1.d0 / det) * mat_temp1(1, 2)
 mat_temp2(2, 1) = (-1.d0 / det) * mat_temp1(2, 1)

 mat_temp3 = MATMUL (mat_temp2, mat_plan_trans)
 mat_proj  = MATMUL (mat_plan, mat_temp3)

 do i = 1, n_pas
  knot_proj(i,:) = MATMUL (mat_proj, tab_knot(i,:))
 enddo

 do i=1,n_pas
   x(i)=scalprod(knot_proj(i,:),vect_1(:))
   y(i)=scalprod(knot_proj(i,:),vect_2(:))
   temp=scalprod(knot_proj(i,:),vect_3(:))
   if(abs(temp)>1d-6) &
   write(23,*) 'routine projec : projection sur e3 non nulle !'
 enddo

end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  subroutine rotateTRI(face,A,B,vect_1,vect_2,vect_3,faceout,Aout,Bout)
    implicit none
    real(8) :: A(3),B(3),face(3,3),veca(3),vecb(3),vecc(3),vect_1(3),vect_2(3),vect_3(3)
    real(8) :: faceout(3,3),Aout(3),Bout(3),mat(3,3),det
    complex(8) :: det2
    integer :: colin,pdet
    colin=0
    mat(:,1)=vect_1(:)
    mat(:,2)=vect_2(:)
    mat(:,3)=vect_3(:)
    call invmat(3,mat,det2,det,pdet)
    if(abs(det)<1.d-12) then
     colin=1
     stop
     return
    endif
    Aout = matmul(mat,A)
    Bout = matmul(mat,B)
    faceout(1,:)=matmul(mat,face(1,:))
    faceout(2,:)=matmul(mat,face(2,:))
    faceout(3,:)=matmul(mat,face(3,:))
   end subroutine


!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  subroutine test_inout_3d
  implicit none
  integer,parameter :: nn=10000
  real(8) :: x(nn,3),b1(3),b2(3),b3(3)
  logical :: out
  integer :: ii
  
  b1=0;b2=0;b3=0

  b1(1)=1.
  b1(2)=0.1
  b1(3)=0.05
  b1=b1/norme(b1)

  b2(1)=0.05
  b2(2)=1.
  b3(3)=0.3
  b2=b2/norme(b2)

  b3(1)=0.3
  b3(2)=0.2
  b3(3)=1.
  b3=b3/norme(b3)
  
  do ii=1,nn
   x(ii,:)=           drand1()*b1
   x(ii,:)= x(ii,:) + drand1()*b2
   x(ii,:)= x(ii,:) + drand1()*b3
   call inout__(x(ii,1),x(ii,2),x(ii,3),b1,b2,b3,out,err=0.d0)
   if(out) write(140,*)  x(ii,:)
  enddo
  stop 'done'

  end subroutine


!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

   !-------------------------------!
   ! inout d'un parallelepipede 3D !
   !-------------------------------!

   subroutine inout__(dx,dy,dz,bbound1,bbound2,bbound3,out,err)
   implicit none
   real(8)          :: aa(3),bb(3),cc(3),pp(3),v1(3),v2(3),v3(3),v4(3),v_(3),dx,dy,dz
   real(8)          :: bbound1(3),bbound2(3),bbound3(3),corner(8,3)
   integer          :: kk,ll,mm,i1,i2,i3,i4
   logical          :: out,inside
   real(8),optional :: err

        pp(1)=dx; pp(2)=dy; pp(3)=dz
        v1=0.d0; v2=bbound1; v3=bbound2; v4=bbound3; v_=v2+v3+v4; v_=v_/norme(v_)*1.d-3
        v1=v1-v_ 
        v2=v2-v_
        v3=v3-v_
        v4=v4-v_

        if(present(err)) then
           pp(1)=pp(1)+err*1.001
           pp(2)=pp(2)+err*1.005
           pp(3)=pp(3)+err*1.003
        else
           pp(1)=pp(1)+0.000000000001313212
           pp(2)=pp(2)+0.000000000001641235
           pp(3)=pp(3)+0.000000000001422121
        endif

        inside= .false.
        corner(1,:)=v1
        corner(2,:)=v2
        corner(3,:)=v3
        corner(4,:)=v4
        corner(5,:)=v4+(v2-v1)
        corner(6,:)=v4+(v3-v1)
        corner(7,:)=v4+(v2-v1)+(v3-v1)
        corner(8,:)=v3+(v2-v1)
        inside =inside .or. & 
               tetra(v1                ,v2           ,v3          ,v4                   ,pp)  .or. &
             & tetra(v3+(v2-v1)        ,v2           ,v3          ,v3+(v2-v1)+(v4-v1)   ,pp)  .or. &
             & tetra(v1+(v3-v1)+(v4-v1),v4           ,v2+(v4-v1)  ,v3                   ,pp)  .or. &
             & tetra(v2+(v4-v1)        ,v2           ,v4          ,v3+(v4-v1)           ,pp)  .or. &
             & tetra(v3+(v4-v1)        ,v3           ,v2+(v4-v1) ,v2+(v3-v1)+(v4-v1)    ,pp)  .or. &
             & tetra(v4                ,v3           ,v2         ,v2+(v3-v2)            ,pp)  .or. &
             & tetra(v2                ,v3           ,v2+(v4-v1) ,v2+(v3-v1)+(v4-v1)    ,pp)
               inside=tetra(corner(1,:),corner(2,:),corner(3,:),corner(6,:),pp).or.inside 
               inside=tetra(corner(1,:),corner(2,:),corner(4,:),corner(6,:),pp).or.inside

       out=.not.inside

   return
   contains

    logical function tetra(v1,v2,v3,v4,pp)
    implicit none
     real(8) :: v1(3),v2(3),v3(3),v4(3),pp(3)
     real    :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x, y, z
      x1=v1(1); y1=v1(2); z1=v1(3);
      x2=v2(1); y2=v2(2); z2=v2(3);
      x3=v3(1); y3=v3(2); z3=v3(3);
      x4=v4(1); y4=v4(2); z4=v4(3);
      x=pp(1) ;  y=pp(2);  z=pp(3);
      call tetra_contains_point_3d( x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x, y, z, tetra)
    end function


   end subroutine

     !=====================================!
     !=====================================!
     !=====================================!
     !=====================================!
     !=====================================!
     !=====================================!
     !=====================================!
     !=====================================!

   !-------------------------------!
   ! inout d'un parallelepipede 2D !
   !-------------------------------!

    subroutine inout_(dx,dy,bbound1,bbound2,out,err,shift)
    implicit none
    real(8)          :: aa(3),bb(3),cc(3),gg(3),pp(3),bound1(2),bound2(2),dx,dy
    real(8)          :: bbound1(2),bbound2(2)
    integer          :: kk,ll,mm
    logical          :: triang1,triang2
    logical          :: out
    real(8),optional :: err
    real(8),optional :: shift(3)

        gg=0.d0; pp=0.d0; pp(1)=dx; pp(2)=dy

        bound1=bbound1; bound2=bbound2

        if(present(shift))then
         pp=pp-shift
         bound1(1)=bound1(1)-shift(1)
         bound2(2)=bound2(2)-shift(2)
        endif

        if(present(err)) then
         gg(1)=-err*1.123123123312
         gg(2)=-err*1.05324234234
        else
         gg(1)=-0.0000000013132121323
         gg(2)=-0.0000000016412351234
        endif

        !triangle 1
        aa=0.; bb=0.; cc=0.
        bb(1)=bound1(1)
        bb(2)=bound1(2)
        cc(1)=bound2(1)
        cc(2)=bound2(2)
        kk=-1 !avec le cote AB
        ll=-1 !avec le cote BC
        mm=-1 !avec le cote CA
        aa=aa+gg; bb=bb+gg; cc=cc+gg
        triang1=PointInTriangle(kk,ll,mm,pp,aa,bb,cc,1.d-11)

        !triangle 2
        aa=0.; bb=0.; cc=0.
        aa(1)=bound1(1)
        aa(2)=bound1(2)
        bb(1)=bound1(1)+bound2(1)
        bb(2)=bound1(2)+bound2(2)
        cc(1)=bound2(1)
        cc(2)=bound2(2)
        kk=1
        ll=1
        mm=1
        aa=aa+gg; bb=bb+gg; cc=cc+gg
        triang2=PointInTriangle(kk,ll,mm,pp,aa,bb,cc,1.d-11)

        out=.not.(triang1.or.triang2)

    end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

      SUBROUTINE POLMAT(RPOL,THAXI,PHAXI,ALFA)
      REAL THAXI,PHAXI,ALFA
      REAL RPOL(3,3)
      REAL SINT,SINTQ,SINP,SINPQ,SINA
      REAL COST,COSTQ,COSP,COSPQ,COSA,EMCOSA
      SINT = SIN(THAXI)
      COST = COS(THAXI)
      SINP = SIN(PHAXI)
      COSP = COS(PHAXI)
      SINA = SIN(ALFA)
      COSA = COS(ALFA)
      EMCOSA = 1.0-COSA
      SINTQ = SINT*SINT
      COSTQ = COST*COST
      SINPQ = SINP*SINP
      COSPQ = COSP*COSP
      RPOL(1,1) =  COSA+COSPQ*SINTQ*EMCOSA
      RPOL(2,1) =  COST*SINA+SINP*COSP*SINTQ*EMCOSA
      RPOL(3,1) = -SINP*SINT*SINA+SINT*COST*COSP*EMCOSA
      RPOL(1,2) = -COST*SINA+SINP*COSP*SINTQ*EMCOSA
      RPOL(2,2) =  COSA+SINPQ*SINTQ*EMCOSA
      RPOL(2,3) = -COSP*SINT*SINA+SINP*SINT*COST*EMCOSA
      RPOL(1,3) =  SINP*SINT*SINA+SINT*COST*COSP*EMCOSA
      RPOL(3,2) =  COSP*SINT*SINA+COST*SINT*SINP*EMCOSA
      RPOL(3,3) =  COSA+COSTQ*EMCOSA
      RETURN
      END subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

      SUBROUTINE CHNAX(THAXI3,PHAXI3,THAXI2,PHAXI2,THAXI1,PHAXI1)
      REAL THAXI1,PHAXI1,PHAXI2,THAXI2,PHAXI3,THAXI3,PI
      PARAMETER (PI=3.14159265359)
      THAXI3 = THAXI3 - PI*0.32
      PHAXI3 = PHAXI3 + PI*0.28
      THAXI2 = THAXI2 + PI*0.18
      PHAXI2 = PHAXI2 - PI*0.14
      THAXI1 = THAXI1 - PI*0.12
      PHAXI1 = PHAXI1 + PI*0.08
      RETURN
      END subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

      SUBROUTINE OFFSET (XOFF,YOFF,ZOFF)
      REAL XOFF,YOFF,ZOFF
      WRITE(*,*)' Rotation with shifting'
      XOFF=-0.0002
      YOFF=+0.0004
      ZOFF=-0.0002
      RETURN
      END subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

      SUBROUTINE ROTX (X, Y, Z, U, V, W, S, C)
      U = X
      V = Y * C + Z * S
      W = - Y * S + Z * C
      END SUBROUTINE ROTX

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

      SUBROUTINE ROTY(X,Y,Z,U,V,W,S,C)
      U=X*C-Z*S
      V=Y
      W=X*S+Z*C
      END subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

      SUBROUTINE ROTZ(X,Y,Z,U,V,W,S,C)
      U=X*C+Y*S
      V=-X*S+Y*C
      W=Z
      END subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

                                                            
      SUBROUTINE VCOPY (X, Y, N) 
      REAL X ( * ), Y ( * ) 
      DO I = 1, N 
       Y (I) = X (I) 
      enddo
      END SUBROUTINE VCOPY                          
      
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

                                                                 
      SUBROUTINE VRFILL (X, A, N) 
      REAL X ( * ) 
      DO I = 1, N 
        X (I) = A 
      enddo
      END SUBROUTINE VRFILL   

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
 
  subroutine giration (a,longueur,radius)
     implicit none
     integer,intent(in)::longueur
     integer::d,e
     real(8),dimension(:,:),intent(in)::a
     real(8)::compt,rgs
     real(8),intent(out)::radius
     compt=0.d0
     do d=1,longueur-1
       do e=d+1,longueur
         rgs=scalprod(a(e,:)-a(d,:),a(e,:)-a(d,:))
         compt=compt+rgs
       enddo
     enddo
     radius=compt/(longueur**2)
  end subroutine
 
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

function intersection(pt_a1,pt_a2,pt_b1,pt_b2)
implicit none

real(8),dimension(2),intent(in)::pt_a1,pt_a2,pt_b1,pt_b2
real(8),dimension(2)::intersection,vect
real(8),dimension(2,2)::matrice
real(8)::det_2


matrice=reshape((/ pt_a2(2)-pt_a1(2),pt_b2(2)-pt_b1(2), &
   & pt_a1(1)-pt_a2(1),pt_b1(1)-pt_b2(1) /),(/2,2/))
                   
vect=(/pt_a1(1)*pt_a2(2)-pt_a2(1)*pt_a1(2),pt_b1(1)* & 
              & pt_b2(2)-pt_b2(1)*pt_b1(2)/)

det_2=determinant(matrice)

if (abs(det_2)<1.e-5) then

        intersection=(/9999.,9999./)
else
        intersection(1)=determinant(reshape( & 
    & (/vect(1),vect(2),matrice(1,2),matrice(2,2)/),(/2,2/)))/det_2
        intersection(2)=determinant(reshape( & 
    & (/matrice(1,1),matrice(2,1),vect(1),vect(2)/),(/2,2/)))/det_2
endif
        
end function intersection


!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

function writheLB(r,long)
implicit none
integer,intent(in)::long
real(8),dimension(long,3),intent(in)::r
real(8),dimension(3,long)::p_vect,u
real(8),dimension(3,long,long)::rseg
real(8),dimension(long,long)::norm
real(8)::writheLB,dir_writhe,writhe,norm_temp,normp,normi,normj,tes1
integer::i,j
real(8),dimension(long)::phi,khi,eta
real(8),dimension(3)::vect_temp
real(8),dimension(2)::pt_inter

dir_writhe=0.d0
writhe=0.d0
rseg=0.d0

do i=1,long-1
        do j=i+1,long
                rseg(:,i,j)=r(j,:)-r(i,:)
                norm(i,j)=norme(rseg(:,i,j)) 
        enddo
enddo
     
do i=1,long-1
        u(:,i)=rseg(:,i,i+1)/norm(i,i+1)
enddo
u(:,long)=-rseg(:,1,long)/norm(1,long)
       
do i=2,long
        p_vect(:,i)=vecprod(u(:,i-1),u(:,i))
        tes1=norme(p_vect(:,i))
        if (tes1>1.d-5) then
                p_vect(:,i)=p_vect(:,i)/tes1
        else
                print *,'bug p_vect'
        endif
enddo
p_vect(:,1)=vecprod(u(:,long),u(:,1))
tes1=norme(p_vect(:,1))
if (tes1>1.d-5) then
        p_vect(:,1)=p_vect(:,1)/tes1
else
        print *,'bug p_vect'
endif

do i=1,long-3
        do j=i+2,long-1
                pt_inter=intersection(r(i,1:2),r(i+1,1:2),r(j,1:2),r(j+1,1:2))
                norm_temp=norme(pt_inter-r(i,1:2))
                normi=norme(r(i,1:2)-r(i+1,1:2))
                normj=norme(r(j,1:2)-r(j+1,1:2))
                if (norm_temp.le.normi) then
                        norm_temp=norme(pt_inter-r(i+1,1:2))
                        if (norm_temp<normi) then
                                norm_temp=norme(pt_inter-r(j,1:2))
                                if (norm_temp.le.normj) then
                                        norm_temp=norme(pt_inter-r(j+1,1:2))
                                        if (norm_temp<normj) then
              vect_temp(:)=vecprod(rseg(:,i,j),u(:,j))
              dir_writhe=dir_writhe+sign(1.0d0,scalprod(u(:,i),vect_temp))
                                        endif
                                endif
                        endif
                endif
        enddo
enddo
do j=2,long-2
        pt_inter=intersection(r(long,1:2),r(1,1:2),r(j,1:2),r(j+1,1:2))
        norm_temp=norme(pt_inter-r(long,1:2))
        normi=norme(r(long,1:2)-r(1,1:2))
        normj=norme(r(j,1:2)-r(j+1,1:2))
        if (norm_temp.le.normi) then
                norm_temp=norme(pt_inter-r(1,1:2))
                if (norm_temp.le.normi) then
                        norm_temp=norme(pt_inter-r(j,1:2))
                        if (norm_temp.le.normj) then
                                norm_temp=norme(pt_inter-r(j+1,1:2))
                                if (norm_temp.le.normj) then
             vect_temp(:)=vecprod(rseg(:,long,j),u(:,j))
             dir_writhe=dir_writhe+sign(1.0d0,scalprod(u(:,long),vect_temp))
                                endif
                        endif
                endif
        endif
enddo                                          

do i=1,long-1
        phi(i)=acospi(u(1:2,i),(/1.d0,0.d0/))
        eta(i)=acos(p_vect(3,i+1))
        if (eta(i)>pi/2.) then
                eta(i)=pi-eta(i)
                khi(i)=acospi(-p_vect(1:2,i+1),(/1.d0,0.d0/))
        else
                khi(i)=acospi(p_vect(1:2,i+1),(/1.d0,0.d0/))
        endif
enddo
phi(long)=acospi(u(1:2,long),(/1.d0,0.d0/))
eta(long)=acos(p_vect(3,1))
if (eta(long)>pi/2.) then
        eta(long)=pi-eta(long)
        khi(long)=acospi(-p_vect(1:2,1),(/1.d0,0.d0/))
else
        khi(long)=acospi(p_vect(1:2,1),(/1.d0,0.d0/))
endif

do i=1,long-1
        writhe=writhe+asin(sin(eta(i))*sin(phi(i+1)-khi(i)))/2.d0/pi
        writhe=writhe-asin(sin(eta(i))*sin(phi(i)-khi(i)))/2.d0/pi
enddo
writhe=writhe+asin(sin(eta(long))*sin(phi(1)-khi(long)))/2.d0/pi
writhe=writhe-asin(sin(eta(long))*sin(phi(long)-khi(long)))/2.d0/pi

writheLB=writhe+dir_writhe

end function writheLB     

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

subroutine acn_compteur(a,long,sample_rate,acn_value)
implicit none
integer::i,j,k,l,m,long,number_cross,ii,ij,ik,in
integer,intent(in)::sample_rate
real(8),intent(inout)::acn_value
real(8),dimension(long,3),intent(in)::a
real(8),dimension(long,3)::noeud,b,backb
real(8),dimension(3,3)::matrice_pi
real(8)::theta,phi,func,vv(3)
real(8)::xi1,xi2,yi1,yi2,xj1,xj2,yj1,yj2,x1,x2,y1,y2,x11,x22,y11,y22,droitex,droitey


number_cross=0
do k=1,sample_rate

     if(k>1)then
       do i=1,long
         b(i,:)=backb(i,:) 
       enddo
     endif

     vv(:)=RANVECSPHERE()
     call azimutangle(vv,theta,phi)

     matrice_pi(1,1)= dsin(theta)*dcos(phi)
     matrice_pi(2,1)=-dsin(phi)
     matrice_pi(3,1)=-dcos(theta)*dcos(phi)
     matrice_pi(1,2)= dsin(theta)*dsin(phi)
     matrice_pi(2,2)= dcos(phi)
     matrice_pi(3,2)=-dcos(theta)*dsin(phi)
     matrice_pi(1,3)= dcos(theta)
     matrice_pi(2,3)= 0.d0
     matrice_pi(3,3)= dsin(theta)

     if(k.eq.1)then
          do i=1,long
               do j=1,3
                    noeud(i,j)=0.d0
                    do m=1,3
                         noeud(i,j)=noeud(i,j)+(matrice_pi(j,m)*a(i,m))
                    enddo
               enddo
               backb(i,:)=noeud(i,:)
          enddo
     else
          do i=1,long
               do j=1,3
                    noeud(i,j)=0.d0
                    do m=1,3
                         noeud(i,j)=noeud(i,j)+(matrice_pi(j,m)*b(i,m))
                    enddo
               enddo
               backb(i,:)=noeud(i,:)
          enddo
     endif

     do i=1,long-3
          xi1=noeud(i,1)
          yi1=noeud(i,2)
          xi2=noeud(i+1,1)
          yi2=noeud(i+1,2)
          if(xi1.gt.xi2) then
               x1=xi2
               x2=xi1
          else
               x1=xi1
               x2=xi2
          endif
               
          if(yi1.gt.yi2) then
               y1=yi2
               y2=yi1
          else
               y1=yi1
               y2=yi2
          endif
               
          do j=i+2,long-1
               xj1=noeud(j,1)
               yj1=noeud(j,2)
               xj2=noeud(j+1,1)
               yj2=noeud(j+1,2)
               droitex=(((yi2-yi1)/(xi2-xi1))*xi1-((yj2-yj1)/(xj2-xj1))*xj1+yj1-yi1)&
                      &/(((yi2-yi1)/(xi2-xi1))-((yj2-yj1)/(xj2-xj1)))
               droitey=(((xi2-xi1)/(yi2-yi1))*yi1-((xj2-xj1)/(yj2-yj1))*yj1+xj1-xi1) &
                      & /(((xi2-xi1)/(yi2-yi1))-((xj2-xj1)/(yj2-yj1)))
               if(xj1.gt.xj2) then
                    x11=xj2
                    x22=xj1
               else
                    x11=xj1
                    x22=xj2
               endif
                    
               if(yj1.gt.yj2) then
                    y11=yj2
                    y22=yj1
               else
                    y11=yj1
                    y22=yj2
               endif
                    
               if(x1.lt.droitex.and.droitex.lt.x2.and.x11.lt.droitex.and.droitex.lt.x22) then
                    if(y1.lt.droitey.and.droitey.lt.y2.and.y11.lt.droitey.and.droitey.lt.y22) then
                         number_cross=number_cross+1
                    endif
               endif
          enddo
          
          if(i>=2) then
               xj1=noeud(long,1)
               yj1=noeud(long,2)
               xj2=noeud(1,1)
               yj2=noeud(1,2)
               droitex=(((yi2-yi1)/(xi2-xi1))*xi1-((yj2-yj1)/(xj2-xj1))*xj1+yj1-yi1) &
                    & /(((yi2-yi1)/(xi2-xi1))-((yj2-yj1)/(xj2-xj1)))
               droitey=(((xi2-xi1)/(yi2-yi1))*yi1-((xj2-xj1)/(yj2-yj1))*yj1+xj1-xi1) &
                    & /(((xi2-xi1)/(yi2-yi1))-((xj2-xj1)/(yj2-yj1)))
               if(xj1.gt.xj2) then
                    x11=xj2
                    x22=xj1
               else
                    x11=xj1
                    x22=xj2
               endif
                    
               if(yj1.gt.yj2) then
                    y11=yj2
                    y22=yj1
               else
                    y11=yj1
                    y22=yj2
               endif
                    
               if(x1.lt.droitex.and.droitex.lt.x2.and.x11.lt.droitex.and.droitex.lt.x22) then
                    if(y1.lt.droitey.and.droitey.lt.y2.and.y11.lt.droitey.and.droitey.lt.y22) then
                         number_cross=number_cross+1
                    endif
               endif
          endif
     enddo
enddo

acn_value=dble(number_cross)/dble(sample_rate)

end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

subroutine acn_compteur2(a,long,sample_rate,acn_value)
implicit none
integer::i,j,k,l,m,long,number_cross,ii,ij,ik,in,NN,count
integer,intent(in)::sample_rate
real(8) :: vv(3),vperp(3),plan(3,2),vv2(3),vect_1(3),vect_2(3),vect_3(3)
real(8),dimension(:),allocatable ::acx,acy
real(8),intent(inout)::acn_value
real(8),dimension(long,3),intent(in)::a

!CALCUL DE L'ACN

    allocate(acx(long),acy(long))
    acn_value=0.

    do count=1,sample_rate
     vv(:)=RANVECSPHERE()
     vv2(:)=RANVECSPHERE()
     plan(1,1) = vv(1)
     plan(2,1) = vv(2)
     plan(3,1) = vv(3)
     plan(1,2) = vv2(1)
     plan(2,2) = vv2(2)
     plan(3,2) = vv2(3)

    call basegen(plan,vect_1,vect_2,vect_3)
    call project(a, long, vect_1, vect_2, vect_3, acx, acy )
    call crosscount(long,NN,acx,acy,1.d-10)
    acn_value=acn_value+float(NN)/2.
    enddo
    deallocate(acx,acy)
    acn_value=acn_value/float(sample_rate)

!FIN CALCUL

end subroutine


!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************


function acospi(vect1,vect2)
implicit none
real(8)::acospi,temp,temp1,norm1,norm2,res
real(8),intent(in),dimension(2)::vect1,vect2
 temp=vect1(1)*vect2(2)-vect1(2)*vect2(1)
 norm1=sqrt(scalprod(vect1,vect1))
 norm2=sqrt(scalprod(vect2,vect2))
 if(norm1<1.d-6.or.norm2<1.d-6) then
  temp1=1.d0
 else
  temp1=scalprod(vect1/norm1,vect2/norm2)
 endif
 if(abs(temp)>1.d-7) then
  if(temp1>1.) then
     res=0.
   else    
     if(temp1<-1.d0) then
       res=-4.d0*atan(1.0d0)*sign(1.0d0,temp)
     else
       res=-acos(temp1)*sign(1.0d0,temp)
     endif
  endif
 else
   res=0.d0
 endif
 acospi=res
end function acospi


!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

    subroutine azimutangle(vv,theta,phi)
    implicit none
    real(8),intent(in)::vv(3)
    real(8),intent(out)::theta,phi
    phi=0.
    theta=0.
    if(vv(1)/=0.) phi=ATAN(vv(2)/vv(1))
    if(vv(1)<=0.) phi=-phi
    if(vv(1)**2+vv(2)**2/=0.) &
    & theta=ATAN(vv(3)/sqrt(vv(1)**2+vv(2)**2))
    end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

   integer function pointcircle(i,N)
   integer :: N
    pointcircle=i
    if(i>N) pointcircle=i-N
    if(i<0) pointcircle=i+N
   end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

   integer function distcircle(i1,j1,N)
   integer :: i,j,k,l,N,i1,j1
    i=i1
    j=j1
    distcircle = abs(j-i)
    l=abs(j-(i+N))
    if(l < distcircle) distcircle=l
    l=abs(j-(i-N))
    if(l < distcircle) distcircle=l
   end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  subroutine pointdist(rayon,rr,start,point,nclose)
   implicit none
   real(8)           :: rr(:,:),rayon,dtemp
   integer          :: start,point,a1,i,j,k,l
   integer          :: tab(size(rr(:,1))*size(rr(:,1)),2)
   integer          :: count,kk,jj
   integer,optional :: nclose

   a1=size(rr(:,1));tab=0;point=0;start=0;jj=0;tab=0

   if(present(nclose)) then
    nclose=0
   endif

   do j=1,a1
    i=j
    do
     if(i+1>a1) exit
     i=i+1
     dtemp=norme(rr(i,:)-rr(j,:))

     if(dtemp<rayon-1.d-6.and.distcircle(i,j,a1)>3) then
       jj=jj+1
       tab(jj,1)=i 
       tab(jj,2)=j 
     else
       i=i+max(0,INT(dtemp-rayon))
     endif
     if(i>=a1) exit 
    enddo
   enddo

   i=0
   if(jj>0) then
     i=INT(drand1()*jj)+1
     if(i>jj) i=jj 
     point=tab(i,1)
     start=tab(i,2)
   else
     start=0
     point=0
   endif

   dtemp=norme(rr(point,:)-rr(start,:))
   if(dtemp>rayon+1.d-4) stop 'error routine pointdist'
   if(size(rr(1,:))/=3.and.size(rr(:,1))/=2) write(*,*) 'aie aie strange thing in pointdist'
  
   if(present(nclose)) then
    nclose=jj
   endif

  return
  end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  !vv define the rotation, vector to send to ez
  !vvb vector to be rotated

  function rotatevtoez(vv)
  implicit none
  real(8)  :: vv(3),d1,d2,d3
  real(8)  :: A(3,3),norm,norm2,rotatevtoez(3,3)
  integer  :: i

  d1=vv(1)
  d2=vv(2)
  d3=vv(3)

  norm=norme(vv)

  norm2=sqrt(d1**2+d2**2)

  A=0. 

  if(abs(norm2)>1.d-7.and.abs(norm)>1.d-7)then

   A(1,:) =  (/ d1*d3/norm2, d2*d3/norm2, -norm2 /)
   A(1,:) =  A(1,:)/norm

   A(2,:) =  (/ d2, -d1 , 0.d0  /)
   A(2,:) = -A(2,:) / norm2

   A(3,:) =  (/ d1, d2, d3  /)
   A(3,:) =  A(3,:) / norm

  else
   
   do i=1,3
    A(i,i)=1.d0
   enddo

  endif
  
  rotatevtoez=A

  return
  end function

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

subroutine writhesub(n_pas,knot2,writhe)
implicit none
integer                      ::  n_pas
real(8), dimension (3)       :: v1,v11,v2,v22,r12,r34,n1,n2,n3,n4,r23,r13,r24,r14
real(8), dimension (n_pas,3) ::  knot2
real(8)                      ::  x,y,writhe,omega,d1,d2,d3,d4
integer                      ::  i,j,l,k,t,u,v,sign,b

    writhe=0.d0
    do i=3,n_pas+1
     do j=2,i-1
      if(i/=j)then
      if(i/=n_pas+1) then
       v1 =knot2(i,:)
       v11=knot2(i-1,:)
      else
       v1 =knot2(1,:)
       v11=knot2(n_pas,:)
      endif
      if(j/=n_pas+1) then
       v2 =knot2(j,:)
       v22=knot2(j-1,:)
      else
       v2 =knot2(1,:)
       v22=knot2(n_pas,:)
      endif
      !1=v11 2=v1 3=v22 4=v2
      r12=v1 - v11
      r34=v2 - v22
      r14=v2 - v11
      r24=v2 - v1
      r23=v22- v1
      r13=v22- v11

      n1=vecprod(r13,r14)
      n2=vecprod(r14,r24)
      n3=vecprod(r24,r23)
      n4=vecprod(r23,r13)

      if(norme(n1)>1.d-12) n1=n1/norme(n1)
      if(norme(n2)>1.d-12) n2=n2/norme(n2)
      if(norme(n3)>1.d-12) n3=n3/norme(n3)
      if(norme(n4)>1.d-12) n4=n4/norme(n4)

      d1=scalprod(n1,n2)
      d2=scalprod(n2,n3)
      d3=scalprod(n3,n4)
      d4=scalprod(n4,n1)
      if(d1>1.d0)  d1= 1.d0
      if(d1<-1.d0) d1=-1.d0
      if(d2>1.d0)  d2= 1.d0
      if(d2<-1.d0) d2=-1.d0
      if(d3>1.d0)  d3= 1.d0
      if(d3<-1.d0) d3=-1.d0
      if(d4>1.d0)  d4= 1.d0
      if(d4<-1.d0) d4=-1.d0
      if(dabs(d1)>1.+1.d-5) write(*,*) 'error'
      if(dabs(d2)>1.+1.d-5) write(*,*) 'error'
      if(dabs(d3)>1.+1.d-5) write(*,*) 'error'
      if(dabs(d4)>1.+1.d-5) write(*,*) 'error'

      omega=DASIN(d1)+DASIN(d2)+DASIN(d3)+DASIN(d4)

      writhe=writhe+2.d0*omega/4.d0/dacos(-1.d0)*DSIGN(1.d0,scalprod(vecprod(r34,r12),r13))
      endif

     enddo
    enddo

    return
    end subroutine

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

end module
