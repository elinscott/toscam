module minimization_wrapping_module

  IMPLICIT NONE 
  PRIVATE

  PUBLIC                           :: minimize_func_wrapper,dspace_rot,drand1,rotatevtoez,norme,rotate_our_notation
  PUBLIC                           :: rearrange_columns_to_identity,det
  REAL(8),  PARAMETER, PRIVATE     :: zero=0.0_8,one=1.0_8,two=2.0_8,three=3.0_8,four=4.0_8,half=0.5d0
  LOGICAL,    PARAMETER, PRIVATE   :: F=.FALSE.,T=.TRUE.

contains

subroutine rearrange_columns_to_identity(lsize,mat,diagdens)
implicit none
integer                   :: lsize
real(8)                   :: diagdenstemp(lsize),diagdens(lsize),mat(lsize,lsize),mat_temp(lsize,lsize)
integer                   :: i
integer                   :: uu(1),uuu(lsize),uuu_not_placed(lsize),j,k
real(8)                   :: dist(lsize)

 uuu=0
 uuu_not_placed=0
 k=0

 do i=1,lsize
  uu=maxloc(abs(mat(:,i)))
  j=uu(1)
  if(uuu(j)==0)then
    uuu(j)=i
  else
    k=k+1
    uuu_not_placed(k)=i
  endif
 enddo

 k=0
 do i=1,lsize
  if(uuu(i)==0)then
    k=k+1
    uuu(i)=uuu_not_placed(k)
  endif
 enddo

 mat_temp=mat
 diagdenstemp=diagdens
 do i=1,lsize
   mat(:,i)=mat_temp(:,uuu(i))
   diagdens(i)=diagdenstemp(uuu(i))
 enddo

return
end subroutine




real(8) function norme(vv)
implicit none
real(8) :: vv(:)
  norme=sqrt(sum(abs(vv)**2))
end function

function rotatevtoez(vv)
  implicit none
  real(8)  :: vv(3),d1,d2,d3
  real(8)  :: A(3,3),norm,norm2,rotatevtoez(3,3)
  integer  :: i

  d1=vv(1); d2=vv(2); d3=vv(3); norm=norme(vv)

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



function drand1()
implicit none
real(8) :: drand1
real(4) :: r
 call random_number(r)
 drand1=dble(r)
end function

     !===========================!
     ! ONETEP IS WITH :
     ! 1 IS m=-2 : dxy
     ! 2 IS m=-1 : dzy
     ! 3 IS m= 0 : d3z2-r2
     ! 4 IS m= 1 : dxz
     ! 5 IS m= 2 : dx2-y2
     !===========================!
     ! THIS ROUTINE WITH NOTATIONS :
     !     yz
     !     zx
     !     xy
     !   x2-y2
     !   3z2-1
     !===========================!

function dspace_rot(rot_vec)
implicit none
integer :: i,j,k,aa(5)
real(8) :: dspace_rot(5,5)
real(8) :: Rd(5,5),rot_vec(3,3),cmp_rot(5,5)

        CALL Rotate_d(Rd,rot_vec)

        aa(1)=3
        aa(2)=1
        aa(3)=5
        aa(4)=2
        aa(5)=4

        do i=1,5
         do j=1,5
          cmp_rot(i,j)=Rd(aa(i),aa(j))
         enddo
        enddo

       dspace_rot=cmp_rot

end function

SUBROUTINE rotate_our_notation(Rd,R)
 IMPLICIT NONE
 real(8), intent(in)  :: R(3,3)
 real(8), intent(out) :: Rd(5,5)
 real(8)              :: s3,r1,r2,r3,r4,r5,r6,r7,r8,r9

   s3=sqrt(3.d0)

   r1=R(1,1)
   r2=R(1,2)
   r3=R(1,3)
   r4=R(2,1)
   r5=R(2,2)
   r6=R(2,3)
   r7=R(3,1)
   r8=R(3,2)
   r9=R(3,3)

   Rd(1,1)=r2*r4+r1*r5
   Rd(1,2)=r3*r5+r2*r6
   Rd(1,3)=s3*r3*r6
   Rd(1,4)=r3*r4+r1*r6
   Rd(1,5)=2.d0*r1*r4+r3*r6

   Rd(2,1)=r5*r7+r4*r8
   Rd(2,2)=r6*r8+r5*r9
   Rd(2,3)=s3*r6*r9
   Rd(2,4)=r6*r7+r4*r9
   Rd(2,5)= 2.*r4*r7+r6*r9

   Rd(3,1)= s3*r7*r8
   Rd(3,2)= s3*r8*r9
   Rd(3,3)= (3.* r9**2 -1.d0)/2.d0
   Rd(3,4)= s3*r7*r9
   Rd(3,5)= s3*(2.*r7**2 + r9**2 - 1.d0)/2.d0

   Rd(4,1)= r2*r7+r1*r8
   Rd(4,2)=  r3*r8+r2*r9
   Rd(4,3)=  s3*r3*r9
   Rd(4,4)= r3*r7+r1*r9
   Rd(4,5)= 2.*r1*r7 +r3*r9

   Rd(5,1)=r1*r2-r4*r5
   Rd(5,2)=r2*r3-r5*r6
   Rd(5,3)=s3*(r3**2-r6**2)/2.d0
   Rd(5,4)=r1*r3-r4*r6
   Rd(5,5)=(2.*r1**2 + r3**2 - 2.*r4**2 - r6**2)/2.d0

 end subroutine

function det(mat)
implicit none
real(8),intent(in)::mat(:,:)
real(8),dimension(3)::detc
real(8)::det
   detc(1)=mat(1,1)*mat(2,2)*mat(3,3)-mat(1,1)*mat(3,2)*mat(2,3)
   detc(2)=mat(2,1)*mat(3,2)*mat(1,3)-mat(2,1)*mat(1,2)*mat(3,3)
   detc(3)=mat(3,1)*mat(1,2)*mat(2,3)-mat(3,1)*mat(2,2)*mat(1,3)
   det=sum(detc)
end function


SUBROUTINE Rotate_d(Rd, R)

! Rotation of the d electron cubic harmonic orbitals. Real space
! rotation matrix is R(3,3). The cubic wave functions are defined as follows:

! b1 = 2yz/r^2 * a
! b2 = 2xz/r^2 * a
! b3 = 2xy/r^2 * a
! b4 = (x^2-y^2)/r^2 * a
! b5 = (3z^2-r^2)/r^2 * a

! Rotation changes them as follows

! R.b1 =   (R(2,2)*R(3,3) + R(2,3)*R(3,2))*b1 +
!          (R(2,1)*R(3,3) + R(2,3)*R(3,1))*b2
!        + (R(2,1)*R(3,2) + R(2,2)*R(3,1))*b3
!        + (R(2,1)*R(3,1) - R(2,2)*R(3,2))*b4 +
!          sqrt(3.)*R(2,3)*R(3,3)*b5

  IMPLICIT NONE
  real(8), intent(in)  :: R(3,3)
  real(8), intent(out) :: Rd(5,5)

  Rd(1,1) = R(2,2)*R(3,3) + R(2,3)*R(3,2);
  Rd(1,2) = R(2,1)*R(3,3) + R(2,3)*R(3,1);
  Rd(1,3) = R(2,1)*R(3,2) + R(2,2)*R(3,1);
  Rd(1,4) = R(2,1)*R(3,1) - R(2,2)*R(3,2);
  Rd(1,5) = sqrt(3.)*R(2,3)*R(3,3);
  Rd(2,1) = R(1,2)*R(3,3) + R(1,3)*R(3,2);
  Rd(2,2) = R(1,1)*R(3,3) + R(1,3)*R(3,1);
  Rd(2,3) = R(1,1)*R(3,2) + R(1,2)*R(3,1);
  Rd(2,4) = R(1,1)*R(3,1) - R(1,2)*R(3,2);
  Rd(2,5) = sqrt(3.)*R(1,3)*R(3,3);
  Rd(3,1) = R(1,2)*R(2,3) + R(1,3)*R(2,2);
  Rd(3,2) = R(1,1)*R(2,3) + R(1,3)*R(2,1);
  Rd(3,3) = R(1,1)*R(2,2) + R(1,2)*R(2,1);
  Rd(3,4) = R(1,1)*R(2,1) - R(1,2)*R(2,2);
  Rd(3,5) = sqrt(3.)*R(1,3)*R(2,3);
  Rd(4,1) = R(1,2)*R(1,3) - R(2,2)*R(2,3);
  Rd(4,2) = R(1,1)*R(1,3) - R(2,1)*R(2,3);
  Rd(4,3) = R(1,1)*R(1,2) - R(2,1)*R(2,2);
  Rd(4,4) = 0.5*(R(1,1)**2+R(2,2)**2) -0.5*(R(1,2)**2 + R(2,1)**2);
  Rd(4,5) = 0.5*sqrt(3.)*(R(1,3)**2-R(2,3)**2);
  Rd(5,1) = sqrt(3.)*R(3,2)*R(3,3);
  Rd(5,2) = sqrt(3.)*R(3,1)*R(3,3);
  Rd(5,3) = sqrt(3.)*R(3,1)*R(3,2);
  Rd(5,4) = 0.5*sqrt(3.)*(R(3,1)**2-R(3,2)**2);
  Rd(5,5) = -0.5 + 1.5*R(3,3)**2;

return
END SUBROUTINE




!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine minimize_civelli (funct, n, x, f, g, h, w, dfn, xm, hh, eps, mode, maxfn, iprint, iexit)
implicit double precision (a-h,o-z)
integer :: n,mode,iprint,iexit,i,ifn,ij,ik,int,is,itn,iu,iv,jk,k,link,n1,np,i1,ib,idiff,j,nn
integer :: maxfn
dimension x(*), g(*), h(*), w(*), xm(*)
INTERFACE
SUBROUTINE funct(n,x,val)
REAL(8), INTENT(OUT) :: val
INTEGER,   INTENT(IN)  :: n
REAL(8), INTENT(IN)  :: x(n)
END SUBROUTINE
END INTERFACE

if (iprint .ne. 0) write (6,1000)
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
6  h(ij) = zero
5  h(ij) = one
go to 15
10  continue
ij = 1
do 11 i = 2, n
z = h(ij)
if (z .le. zero) then
 write(*,*) 'less than zero'
 return
endif
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
if (h(ij) .le. zero) then
 write(*,*) 'less than zero'
 return
endif
15  continue
ij = np
dmin = h(1)
do 16 i = 2, n
if (h(ij) .ge. dmin) go to 16
dmin = h(ij)
16  ij = ij + np - i
if (dmin .le. zero) then
 write(*,*) 'less than zero'
 return
endif
z = f
itn = 0
call funct (n, x, f)
ifn = 1
df = dfn
if (dfn .eq. zero) df = f - z
if (dfn .lt. zero) df = abs (df * f)
if (df .le. zero) df = one
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
!write (6,1001) itn, ifn
!1001  format (1x,'itn = ',i5,' ifn = ',i5)
!write (6,1002) f
!1002  format (1x,'f = ',e15.7)
if (iprint .lt. 0) go to 21
!write (6,1003) (x(i), i = 1, n)
!1003  format (1x,'x = ',4e15.7 / (5x, 4e15.7))
!write (6,1004) (g(i), i = 1, n)
!1004  format (1x,'g = ',4e15.7 / (5x, 4e15.7))
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
z = zero
do 26 j = 1, i
z = z + h(ij) * w(is+np-j)
ij = ij - 1
26  continue
25  w(is+n-i) = w(n-i) / h(ij) - z
z = zero
gs0 = zero
do 29 i = 1, n
if (z * xm(i) .ge. abs (w(is+i))) go to 28
z = abs (w(is+i)) / xm(i)
28  gs0 = gs0 + g(i) * w(is+i)
29  continue
aeps = eps / z
iexit = 2
if (gs0 .ge. zero) go to 92
alpha = -two * df / gs0
if (alpha .gt. one) alpha = one
ff = f
tot = zero
int = 0
iexit = 1
30  continue
if (ifn .ge. maxfn) go to 90
do 31 i = 1, n
w(i) = x(i) + alpha * w(is+i)
31  continue
call funct (n, w, f1)
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
call funct (n, w, f1)
ifn = ifn + 1
if (f1 .ge. f) go to 50
if ((f1 + f2 .ge. f + f) .and. (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
tot = tot + alpha
alpha = two * alpha
go to 32
40  continue
if (alpha .lt. aeps) go to 92
if (ifn .ge. maxfn) go to 90
alpha = half * alpha
do 41 i = 1, n
w(i) = x(i) + alpha * w(is+i)
41  continue
call funct (n, w, f2)
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
if (f1 + f .gt. f2 + f2)  z = one + half * (f - f1) / (f + f1 - f2 - f2)
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
gys = zero
do 55 i = 1, n
w(i) = w(ib+i)
gys = gys + g(i) * w(is+i)
55  continue
df = ff - f
dgs = gys - gs0
if (dgs .le. zero) go to 20
link = 1
if (dgs + alpha * gs0 .gt. zero) go to 52
do 51 i = 1, n
w(iu + i) = g(i) - w(i)
51  continue
sig = one / (alpha * dgs)
go to 70
52  continue
zz = alpha / (dgs - alpha * gs0)
z = dgs * zz - one
do 53 i = 1, n
w(iu+i) = z * w(i) + g(i)
53  continue
sig = one / (zz * dgs * dgs)
go to 70
60  continue
link = 2
do 61 i = 1, n
w(iu+i) = w(i)
61  continue
if (dgs + alpha * gs0 .gt. zero) go to 62
sig = one / gs0
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
if (z .le. zero) z = dmin
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
if (iprint .eq. 0) then
 write(*,*) 'iprint'
 return
endif
!write (6,1005) itn, ifn, iexit
!1005  format (1x,'itn = ',i5, ' ifn = ',i5,' iexit = ',i5)
!write (6,1002) f
!write (6,1003) (x(i), i = 1, n)
!write (6,1004) (g(i), i = 1, n)
write(*,*) '100 return'
return
100  continue
do 101 i = 1, n
z = hh * xm(i)
w(i) = w(i) + z
call funct (n, w, f1)
g(i) = (f1 - f) / z
w(i) = w(i) - z
101  continue
ifn = ifn + n
go to (18, 54), link
110  continue
do 111 i = 1, n
z = hh * xm(i)
w(i) = w(i) + z
call funct (n, w, f1)
w(i) = w(i) - z - z
call funct (n, w, f2)
g(i) = (f1 - f2) / (two * z)
w(i) = w(i) + z
111  continue
ifn = ifn + n + n
go to (18, 54), link
end subroutine


!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE minimize_func_wrapper(func_,test,nnn,Niter_search_max_, &
&                                     dist_min,dist_max,search_step)
  implicit none
    INTEGER                        :: iparam,ii,nnn,Niter_search_max,Niter_search_max_
    INTEGER                        :: iexit,ifault,i,j
    real(8)                        :: dist_min,dist_max,search_step
    INTEGER                        :: mode,iprint
    REAL(8)                        :: dfn
    real(8),allocatable            :: g(:),xtemp(:),xprmt(:),hess(:),w(:),WORK(:)
    real(8),allocatable            :: var(:),stepvec(:)
    real(8)                        :: test(:)
    real(8)                        :: dist_init

    INTERFACE
     SUBROUTINE func_(dist,n,vec)
        REAL(8), INTENT(OUT)   :: dist
        INTEGER,   INTENT(IN)  :: n
        REAL(8), INTENT(IN)    :: vec(n)
     END SUBROUTINE
    END INTERFACE


    write(*,*) 'starting minimization'

   allocate(g(nnn),xtemp(nnn),xprmt(nnn),hess(nnn**2),w(4*nnn),WORK(max(5*nnn+2,nnn*(nnn+7)/2)))
   allocate(var(nnn),stepvec(nnn))

   if(size(test)/=nnn)then
    write(*,*) 'ERROR - minimize_func_wrapper, dimension of test do no match the argument'
    write(*,*) 'size of test : ', size(test)
    write(*,*) 'nnn          : ', nnn
    stop
   endif

write(*,*) 'starting minimization'

   Niter_search_max=Niter_search_max_

   write(*,*) '============================================='
   write(*,*) '---- MINIMIZATION WRAPPER -----'
   write(*,*) ' NPARAM         : ', nnn
   write(*,*) ' size test      : ', size(test)
   write(*,*) ' NITER          : ', Niter_search_max
   write(*,*) ' dist max       : ', dist_max
   write(*,*) ' dist min       : ', dist_min
   write(*,*) '============================================='

   write(*,*) 'starting minimization'

   stepvec   =  0.05d0
   iprint    =  1
   mode      =  1
   dist_min  =  0.d0
   hess      =  0.d0
   g         =  0.d0
   xprmt     =  0.d0
   w         =  0.d0
   var       =  0.d0

  !---------------------------------------------------------------------------------------------!
   call init_minimize
   CALL minimize_civelli(distance_func______,nnn,xtemp,dist_min,g,hess,w,dfn,& 
&                          xprmt,search_step,dist_max,mode,Niter_search_max,iprint,iexit)
        test = xtemp(1:nnn)
  !---------------------------------------------------------------------------------------------!

     call func_(dist_min,nnn,test);
     write(*,*) ' ------- '
     write(*,*) ' value at minimum is : ', dist_min
     write(*,*) ' ...end....'
     write(*,*) '============================================='
     write(*,*) '============================================='
     write(*,*) '============================================='

   deallocate(g,xtemp,xprmt,hess,w,work,var,stepvec)

  contains

   !---------------------!

  SUBROUTINE distance_func______(n,x,f)
    REAL(8), INTENT(OUT) :: f
    INTEGER, INTENT(IN)  :: n
    REAL(8), INTENT(IN)  :: x(n)
    call func_(f,n,x)
  END SUBROUTINE

   !---------------------!

  subroutine init_minimize
  implicit none
        iexit        = 0
        dfn          = -half
        xtemp        =  zero
        xtemp(1:nnn) =  test
        xprmt        =  zero
        xprmt(1:nnn) =  ABS(test) + 1.d-12
        write(*,*) 'calling initial distance'
        call func_(dist_init,nnn,test)
        write(*,*) 'initial distance is : ', dist_init
  end subroutine

   !---------------------!

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

end module















