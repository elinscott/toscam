module linalg 

 use genvar
 use geomlib

!---------------------------------------------------!
 INTERFACE norm
   MODULE PROCEDURE norm_cplx_vec
   MODULE PROCEDURE norm_real_vec
 END INTERFACE
!---------------------------------------------------!
 INTERFACE conj
   MODULE PROCEDURE conj_cplx
   MODULE PROCEDURE conj_real
 END INTERFACE
!---------------------------------------------------!
 INTERFACE ramp
   MODULE PROCEDURE ramp_vec
   MODULE PROCEDURE ramp_mat
 END INTERFACE
!---------------------------------------------------!
INTERFACE outerprod
 MODULE PROCEDURE outerprod_r,outerprod_c,outerprod_rs,outerprod_cs
END INTERFACE
 !---------------------------------------------------!
INTERFACE reverse_array
 MODULE PROCEDURE reverse_array_,reverse_array__,reverse_array___,reverse_array____
END INTERFACE
 !---------------------------------------------------!
INTERFACE printmatrice
 MODULE PROCEDURE printmatrice__,printmatrice_
END INTERFACE
 !---------------------------------------------------!
INTERFACE put_and_shift
 MODULE PROCEDURE put_and_shift_,put_and_shift__,put_and_shift___,put_and_shift____
END INTERFACE 
 !---------------------------------------------------!
INTERFACE maxontop_array
 MODULE PROCEDURE maxontop_array_d,maxontop_array_r,maxontop_array_d_,maxontop_array_r_,maxontop_array_d__,maxontop_array_r__
END INTERFACE 
 !---------------------------------------------------!
INTERFACE nullarray
 MODULE PROCEDURE nullarray_d,nullarray_r,nullarray_d_,nullarray_r_,nullarray_d__,nullarray_r__
END INTERFACE 
 !---------------------------------------------------!
INTERFACE constantarray
 MODULE PROCEDURE constantarray_d,constantarray_r,constantarray_d_,constantarray_r_,constantarray_d__,constantarray_r__
END INTERFACE 
 !---------------------------------------------------!
INTERFACE matrixinverse
 MODULE PROCEDURE matrixinverse_r,matrixinverse_c,matrixinverse_qr,matrixinverse_qc,matrixinverse_rs,matrixinverse_cs
END INTERFACE
 !---------------------------------------------------!
INTERFACE same_array
 MODULE PROCEDURE same_array_r,same_array_rr,same_array_c,same_array_cc
END INTERFACE
 !---------------------------------------------------!
INTERFACE dble
 MODULE PROCEDURE dble_r,dble_q
END INTERFACE
 !---------------------------------------------------!
INTERFACE aimag
 MODULE PROCEDURE aimag_r,aimag_q
END INTERFACE
 !---------------------------------------------------!
INTERFACE cabs1
 MODULE PROCEDURE cabs1_r,cabs1_q
END INTERFACE
 !---------------------------------------------------!
INTERFACE csign1
 MODULE PROCEDURE csign1_r,csign1_q
END INTERFACE
 !---------------------------------------------------!
INTERFACE DEXPc
 MODULE PROCEDURE DEXPc_r,DEXPc_q,DEXPc_rr
END INTERFACE
 !---------------------------------------------------!
INTERFACE ZSUM
 MODULE PROCEDURE ZSUM_,ZSUM__ 
END INTERFACE
 !---------------------------------------------------!
INTERFACE ZAXPYb
 MODULE PROCEDURE ZAXPY__,ZAXPY_
END INTERFACE
 !---------------------------------------------------!
INTERFACE ZDOT
 MODULE PROCEDURE ZDOT_,ZDOT__
END INTERFACE
 !---------------------------------------------------!
INTERFACE norme
 MODULE PROCEDURE dnorm,inorm,snorm,cnorm
END INTERFACE
 !---------------------------------------------------!
INTERFACE scalprod
 MODULE PROCEDURE idot_,ddot_,sdot_,cdot_
END INTERFACE
 !---------------------------------------------------!
INTERFACE vecprod
 MODULE PROCEDURE dcross,icross,scross,ccross
END INTERFACE
 !---------------------------------------------------!
INTERFACE swap
 module procedure swapi,swaps,swapr,swapa,swapc,swapqc,swapqr,zswap_,zswap__,swapcs
END INTERFACE
 !---------------------------------------------------!
INTERFACE errorsmallenough
 module procedure error_r,error_c
END INTERFACE
 !---------------------------------------------------!
INTERFACE packarray
 module procedure pack_r,pack_c,pack_i,pack_s,unpack_r,unpack_c,unpack_i,unpack_s
END INTERFACE
 !---------------------------------------------------!
INTERFACE packarray_
 module procedure pack_r_,pack_i_,pack_s_,unpack_r_,unpack_i_,unpack_s_,pack_c_,unpack_c_
END INTERFACE
 !---------------------------------------------------!
INTERFACE minloc_array
 module procedure minloc_with_test_c,minloc_with_test_r
END INTERFACE
 !---------------------------------------------------!
INTERFACE nambu_pack
 module procedure nambu_pack_,nambu_pack__
END INTERFACE
 !---------------------------------------------------!
INTERFACE prod
 module procedure prod_c,prod_r,prod_i
END INTERFACE
 !---------------------------------------------------!
INTERFACE percent
 module procedure per_c,per_c_vec,per_c_mat
END INTERFACE
 !---------------------------------------------------!
INTERFACE array_to_int
 module procedure array_to_int_1,array_to_int_2
END INTERFACE
 !---------------------------------------------------!
INTERFACE int_to_array
 module procedure int_to_array_1,int_to_array_2
END INTERFACE
 !---------------------------------------------------!
INTERFACE erase_divergence
 module procedure rerase_divergence_vec,rerase_divergence_mat,   &
                & cerase_divergence_vec,cerase_divergence_mat,   & 
                & ierase_divergence_vec,ierase_divergence_mat,   &
                & ierase_divergence_scal,cerase_divergence_scal, &
                & rerase_divergence_scal,rrerase_divergence_vec, &
                & rrerase_divergence_mat,rrerase_divergence_scal, &
                & rerase_divergence_mat___,rerase_divergence_mat__, &
                & rerase_divergence_mat_,cerase_divergence_mat___,cerase_divergence_mat__, &
                & cerase_divergence_mat_,rerase_divergence_mat____

END INTERFACE
 !---------------------------------------------------!
INTERFACE erase_neg
 MODULE PROCEDURE erase_negr,erase_negd,erase_negr_,erase_negd_,erase_negr__,erase_negd__,erase_negr___,erase_negd___
END INTERFACE
 !---------------------------------------------------!
INTERFACE PHASE
   MODULE PROCEDURE PHASE_cel,PHASE_rel,PHASE_iel,PHASE_rel2
END INTERFACE
 !---------------------------------------------------!
INTERFACE ANGLE
   MODULE PROCEDURE ANGLE_c,ANGLE_c2,ANGLE_cb,ANGLE_r,ANGLE_rs
END INTERFACE
 !---------------------------------------------------!
interface
 logical function disnan(x)
  real(8) :: x
 end function
 logical function disinf(x)
  real(8) :: x
 end function
 logical function isnan(x)
  real(4) :: x
 end function
 logical function isinf(x)
  real(4) :: x
 end function
end interface
 !---------------------------------------------------!
INTERFACE maxloci
 module procedure maxloci_c,maxloci_r,maxloci_i,maxloci_rr
END INTERFACE
 !---------------------------------------------------!
INTERFACE minloci
 module procedure minloci_c,minloci_r,minloci_i,minloci_rr
END INTERFACE
 !---------------------------------------------------!
INTERFACE ISNAN_TEST 
 module procedure ISNANR,ISNANC,ISNANS 
END INTERFACE
 !---------------------------------------------------!
INTERFACE ISINF_TEST 
 module procedure ISINFR,ISINFC,ISINFS
END INTERFACE
 !---------------------------------------------------!
INTERFACE too_large
 module procedure too_large_r,too_large_i,too_large_c,&
                & too_large_c_vec,too_large_c_mat
END INTERFACE
 !---------------------------------------------------!
INTERFACE pol
 MODULE PROCEDURE polr,polc,pold,poli
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (.DP.)
 MODULE PROCEDURE dirprod
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (.QN.)
 MODULE PROCEDURE qnumber
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (<)
   MODULE PROCEDURE comparearray_s,comparearray_sc
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (>)
   MODULE PROCEDURE comparearray_g,comparearray_gc
END INTERFACE
 !---------------------------------------------------!
INTERFACE none_zero_term_in_there
 MODULE PROCEDURE none_zero_term_in_there_i, &
               &  none_zero_term_in_there_r, &
               &  none_zero_term_in_there_c
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (.eqt.)
   MODULE PROCEDURE vecegali,vecegalr,vecegalis,vecegalrs
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (.bgt.)
   MODULE PROCEDURE vecbigi,vecbigr,vecbigis,vecbigrs
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (.smt.)
   MODULE PROCEDURE vecsmi,vecsmr,vecsmis,vecsmrs
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (+)
   MODULE PROCEDURE add_r,add_c
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (-)
   MODULE PROCEDURE sub_r,sub_c
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (*)
   MODULE PROCEDURE mult_r,mult_c
END INTERFACE
 !---------------------------------------------------!
INTERFACE OPERATOR (/)
   MODULE PROCEDURE div_r,div_c
END INTERFACE
 !---------------------------------------------------!
INTERFACE write_first_el
   MODULE PROCEDURE write_first_el_,write_first_el__
END INTERFACE
 !---------------------------------------------------!
INTERFACE floor_
   MODULE PROCEDURE floor_r,floor_c
END INTERFACE
 !---------------------------------------------------!
TYPE arrayr
 real(8),dimension(:),allocatable :: subarray
END TYPE
 !---------------------------------------------------!
TYPE arrayi
 integer,dimension(:),allocatable :: subarray
END TYPE
 !---------------------------------------------------!
TYPE arrayc
 complex(8),dimension(:),allocatable :: subarray
END TYPE
 !---------------------------------------------------!


 real(8),  private  :: prec=1.d-5
 real(8)            :: vnull(3)=(/0.d0,0.d0,0.d0/)


contains

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION norm_cplx_vec(vec) RESULT(norm)
    REAL(DBL)    :: norm
    COMPLEX(DBL) :: vec(:)
    norm = SQRT(DBLE(DOT_PRODUCT(vec,vec)))
  END FUNCTION


  FUNCTION norm_real_vec(vec) RESULT(norm)
    REAL(DBL) :: norm
    REAL(DBL) :: vec(:)
    norm = SQRT(      DOT_PRODUCT(vec,vec))
  END FUNCTION


  ELEMENTAL FUNCTION conj_real(r) RESULT(cc)
    REAL(DBL)             :: cc
    REAL(DBL), INTENT(IN) :: r
    cc = r
  END FUNCTION


  ELEMENTAL FUNCTION conj_cplx(c) RESULT(cc)
    COMPLEX(DBL)             :: cc
    COMPLEX(DBL), INTENT(IN) :: c
    cc = CONJG(c)
  END FUNCTION



  SUBROUTINE ramp_vec(zeramp,n)
    INTEGER, INTENT(INOUT)        :: zeramp(:)
    INTEGER, INTENT(IN), OPTIONAL :: n
    INTEGER                       :: i,nn
    nn = SIZE(zeramp)
    IF(PRESENT(n)) nn = n
    zeramp(1:nn) = (/(i,i=1,nn)/)
  END SUBROUTINE


  SUBROUTINE ramp_mat(zeramp,n1,n2)
    INTEGER, INTENT(INOUT)        :: zeramp(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: n1,n2
    INTEGER                       :: i,nn1,nn2
    nn1 = SIZE(zeramp,1)
    IF(PRESENT(n1)) nn1 = n1
    nn2 = SIZE(zeramp,2)
    IF(PRESENT(n2)) nn2 = n2
    zeramp(1:nn1,1:nn2) = RESHAPE((/(i,i=1,nn1*nn2)/),(/nn1,nn2/))
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 integer function last_non_zero_array_element(arr)
 implicit none
  integer :: arr(:)
  integer :: i
  last_non_zero_array_element=0
  do i=1,size(arr)
   if(arr(i)/=0) then
    last_non_zero_array_element=i
   endif
  enddo
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 integer(8) function long_sum(x)
 integer(4) :: x(:)
  long_sum=0
  do i=1,size(x)
   long_sum=long_sum+x(i)
  enddo
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine write_first_el__(n,v)
implicit none
 integer    :: n
 real(8) :: v(:)
 write(*,*)
 write(*,*) '####################################'
 write(*,*) v(1:min(n,size(v)))
 write(*,*) '####################################'
end subroutine

subroutine write_first_el_(n,v)
implicit none
 integer    :: n
 complex(8) :: v(:)
 write(*,*) '####################################'
 write(*,*) v(1:min(n,size(v)))
 write(*,*) '####################################'
end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  !------------------------!

   real(16) function q_zsum(n,z,inc)
   implicit none
   integer     :: inc,n,i
   complex(16) :: z(n)
    q_zsum=0.d0
    do i=1,n,inc
     q_zsum = q_zsum + abs(real(z(i)))+abs(aimag(z(i)))
    enddo
   end function

  !------------------------!

   subroutine  q_zscal(n,za,zx,incx)
      complex(16) za,zx(*)
      integer i,incx,ix,n
      if( n<=0 .or. incx<=0 )return
      if(incx==1)go to 20
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end subroutine

  !------------------------!

      integer function q_zamax(n,zx,incx)
      complex(16) zx(*)
      real(16)    smax
      integer i,incx,ix,n
      q_zamax = 0
      if( n<1 .or. incx<=0 )return
      q_zamax = 1
      if(n==1)return
      if(incx==1)go to 20
      ix = 1
      smax = cabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(zx(ix))<=smax) go to 5
         q_zamax = i
         smax = cabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
   20 smax = cabs1(zx(1))
      do 30 i = 2,n
         if(cabs1(zx(i))<=smax) go to 30
         q_zamax = i
         smax = cabs1(zx(i))
   30 continue
      return
      end function
   !---------------------------!

      integer function q_damax(n,zx,incx)
      real(16) zx(*)
      real(16)    smax
      integer i,incx,ix,n
      q_damax = 0
      if( n<1 .or. incx<=0 )return
      q_damax = 1
      if(n==1)return
      if(incx==1)go to 20
      ix = 1
      smax = abs(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(zx(ix))<=smax) go to 5
         q_damax = i
         smax = abs(zx(ix))
    5    ix = ix + incx
   10 continue
      return
   20 smax = abs(zx(1))
      do 30 i = 2,n
         if(abs(zx(i))<=smax) go to 30
         q_damax = i
         smax = abs(zx(i))
   30 continue
      return
      end function

   !---------------------------!

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 pure real(8) function floor_r(rr,i)
 implicit none
 integer,intent(in)  :: i
 integer             :: ii
 real(16)            :: dd,pp
 real(8),intent(in)  :: rr 
  dd=1.d0
  do ii=1,i
   dd=10._16*dd
  enddo
  pp =  ANINT ( dd*rr ) / dd
  floor_r = pp
 end function

   !-----------------------------!

 pure complex(8) function floor_c(rr,i)
 implicit none
 integer,intent(in)    :: i
 integer               :: ii
 complex(8),intent(in) :: rr
  floor_c= floor_r(real(rr),i) + imi * floor_r(aimag(rr),i)
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 logical function upper_triang_mat(vec)
 implicit none
  integer :: vec(2)
  upper_triang_mat=vec(1)<=vec(2)
 end function

 logical function lower_triang_mat(vec)
 implicit none
  integer :: vec(2)
  lower_triang_mat=vec(1)>=vec(2)
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine nambu_unpack(mat)
  implicit none
  complex(8) :: mat(:,:)
  integer    :: i,j,k,l,siz,si
  
   siz = size(mat,1)
   si  = siz/2
   if(siz/=size(mat,2)) stop 'nambu_unpack matrix is not square'
   mat(si+1:siz,si+1:siz) = - transpose(mat(si+1:siz,si+1:siz))
   mat(1:si,si+1:siz)     =   conjg(mat(1:si,si+1:siz))

  end subroutine

    !-----------------------------!

  subroutine nambu_pack__(mat)
  implicit none
  complex(8) :: mat(:,:)
  integer    :: i,j,k,l,siz,si

   siz = size(mat,1)
   si  = siz/2
   if(siz/=size(mat,2)) stop 'nambu_unpack matrix is not square'
   mat(si+1:siz,si+1:siz) = - transpose(mat(si+1:siz,si+1:siz))
   mat(1:si,si+1:siz)     =   conjg(mat(1:si,si+1:siz))

  end subroutine

    !-----------------------------!

  subroutine nambu_pack_(mat)
  implicit none
  complex(8) :: mat(:)
  integer    :: i,j,k,l,siz,si
   siz = size(mat,1); si  = siz/2
   mat(si+1:siz) = -mat(si+1:siz)
  end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  real(8) function qsign(quad)
  implicit none
  real(16) :: quad
   if(quad>0.) then
    qsign=1.d0
   else
    qsign=-1.d0
   endif
  end function 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

         !-------------------------!

      subroutine zaxpy__(n,za,zx,incx,zy,incy)
      implicit none
      complex(16) :: zx(n),zy(n),za
      integer     :: n,incx,ix,iy,i,incy
 
      if(n.le.0)return
      if (cabs1(za)<1.d-30) return
      if (incx.eq.1.and.incy.eq.1) go to 20

      ix = 1 ; iy = 1

      if(incx.lt.0) ix = ( -n + 1 ) * incx + 1
      if(incy.lt.0) iy = ( -n + 1 ) * incy + 1

      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue

      return

   20 do 30 i = 1, n
        zy(i) = zy(i) + za*zx(i)
   30 continue

      return
      end subroutine

         !-------------------------!

      subroutine zaxpy_(n,za,zx,incx,zy,incy)
      implicit none
      complex(8)  ::  zx(n),zy(n),za
      complex(16) ::  tt
      integer     :: n,incx,incy,i,ix,iy

      if(n.le.0)return
      if (cabs1(za) < 1.d-16) return
      if (incx.eq.1.and.incy.eq.1) go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

   20 do 30 i = 1,n
        tt=zy(i)+za*zx(i)
        if(abs(tt)<MAX_REAL)then
         zy(i)=tt
        else
         zy(i)=MAX_REAL
        endif
   30 continue

      return
      end subroutine

         !-------------------------!

      COMPLEX(16) FUNCTION ZDOT__(N,ZX,INCX,ZY,INCY)
      implicit none
      INTEGER     :: INCX,INCY,N
      COMPLEX(16) :: ZX(*),ZY(*)
      COMPLEX(16) :: ZTEMP
      INTEGER     :: I,IX,IY

      ZTEMP = (0.0d0,0.0d0)
      ZDOT__ = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZTEMP = ZTEMP + CONJG(ZX(IX))*ZY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      ZDOT__ = ZTEMP
      RETURN
   20 DO 30 I = 1,N
          ZTEMP = ZTEMP + CONJG(ZX(I))*ZY(I)
   30 CONTINUE
      ZDOT__ = ZTEMP

      RETURN
      END function

         !-------------------------!

      COMPLEX(8) FUNCTION ZDOT_(N,ZX,INCX,ZY,INCY)
      implicit none
      INTEGER        ::  INCX,INCY,N
      complex(8)     ::  ZX(*),ZY(*)
      COMPLEX(16)    ::  tt
      INTEGER        ::  I,IX,IY
      ZDOT_ = 0.0d0
      tt=0. 
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          tt = tt + CONJG(ZX(IX))*ZY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE 

      if(abs(tt)<MAX_REAL)then
       ZDOT_ = tt
      else
       ZDOT_ = MAX_REAL
      endif
      RETURN

   20 DO 30 I = 1,N
          tt = tt + CONJG(ZX(I))*ZY(I)
   30 CONTINUE

      if(abs(tt)<MAX_REAL)then
      ZDOT_ = tt
      else
      ZDOT_ = MAX_REAL
      endif

      RETURN
      END function

         !-------------------------!

      REAL(8) FUNCTION ZSUM__(N,ZX,J) 
      implicit none
      INTEGER          :: N
      COMPLEX(8)       :: ZX(*)
      INTEGER          :: I,IX,J
      REAL(8),SAVE     :: t
      real(16)         :: tt

       if(J>0) then
        ZSUM__=t
        return
       endif
       tt=0.
       DO I = 1,N
           tt = tt + ABS(real(ZX(I),kind=8)) + ABS(aimag(ZX(I)))
       enddo 
       if(abs(tt)<MAX_REAL)then
        ZSUM__=tt
       else
        ZSUM__= MAX_REAL
       endif
       t=ZSUM__
      END function

         !-------------------------!

      FUNCTION ZSUM_(N,ZX,J)
      INTEGER        ::  N,J
      COMPLEX(16)    :: ZX(*)
      REAL(16)       :: ZSUM_
      REAL(16),SAVE  :: t
      INTEGER I,IX
       if(J>0) then
        ZSUM_=t 
        return
       endif
       ZSUM_ = 0.
       DO I = 1,N
           ZSUM_ = ZSUM_ + ABS(real(ZX(I),kind=16)) + ABS(aimag(ZX(I)))
       enddo 
       t=ZSUM_
      END function

         !-------------------------!

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
!**************************************************************************

 integer function invcell(a1,a2,a3)
 implicit none
 integer :: a1,a2,a3,cadran

!-------!
! 1 3 2 !
! 4 5 6 !
! 8 7 9 !
!-------!

    if(a1==-1.and.a2==-1) cadran=8
    if(a1==-1.and.a2== 0) cadran=4
    if(a1==-1.and.a2== 1) cadran=1
    if(a1== 0.and.a2==-1) cadran=7
    if(a1== 0.and.a2== 0) cadran=5
    if(a1== 0.and.a2== 1) cadran=3
    if(a1== 1.and.a2==-1) cadran=9
    if(a1== 1.and.a2== 0) cadran=6
    if(a1== 1.and.a2== 1) cadran=2

    if(a3== 1) cadran =  cadran+9
    if(a3==-1) cadran =  cadran+18

    invcell=cadran

 return
 end function


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine cells(cadran,aa,bb,cc)
 integer   :: cadran
 real(8)   :: aa,bb,cc

 aa=0.d0;bb=0.d0;cc=0.d0

 SELECT CASE (cadran)
 !----------------------------!
 CASE(2)
  aa=1.d0 ;  bb=1.d0
 CASE(3)
  aa=0.d0 ;  bb=1.d0
 CASE(4)
  aa=-1.d0;  bb= 0.d0
 CASE(5)
  aa=0.d0 ;  bb=0.d0
 CASE(6)
  aa=1.d0 ;  bb=0.d0
 CASE(7)
  aa= 0.d0;  bb=-1.d0
 CASE(8)
  aa=-1.d0;  bb=-1.d0
 CASE(9)
  aa= 1.d0;  bb=-1.d0
 CASE(1)
  aa=-1.d0;  bb= 1.d0
 !----------------------------!
 CASE(20)
  aa=1.d0 ;  bb=1.d0  ; cc=-1.d0
 CASE(21)
  aa=0.d0 ;  bb=1.d0  ; cc=-1.d0
 CASE(22)
  aa=-1.d0;  bb=0.d0  ; cc=-1.d0
 CASE(23)
  aa=0.d0 ;  bb=0.d0  ; cc=-1.d0
 CASE(24)
  aa=1.d0 ;  bb=0.d0  ; cc=-1.d0
 CASE(25)
  aa= 0.d0;  bb=-1.d0 ; cc=-1.d0
 CASE(26)
  aa=-1.d0;  bb=-1.d0 ; cc=-1.d0
 CASE(27)
  aa= 1.d0;  bb=-1.d0 ; cc=-1.d0
 CASE(19)
  aa=-1.d0;  bb= 1.d0 ; cc=-1.d0
 !----------------------------!
 CASE(11)
  aa=1.d0 ;  bb=1.d0  ; cc=1.d0
 CASE(12)
  aa=0.d0 ;  bb=1.d0  ; cc=1.d0
 CASE(13)
  aa=-1.d0;  bb= 0.d0 ; cc=1.d0
 CASE(14)
  aa=0.d0 ;  bb=0.d0  ; cc=1.d0
 CASE(15)
  aa=1.d0 ;  bb=0.d0  ; cc=1.d0
 CASE(16)
  aa= 0.d0;  bb=-1.d0 ; cc=1.d0
 CASE(17)
  aa=-1.d0;  bb=-1.d0 ; cc=1.d0
 CASE(18)
  aa= 1.d0;  bb=-1.d0 ; cc=1.d0
 CASE(10)
  aa=-1.d0;  bb= 1.d0 ; cc=1.d0
 !----------------------------!
 END SELECT

 return
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 !--------------------------!

 function realvec_to_complex(vec)
 implicit none
 real(8)    :: vec(:)
 complex(8) :: realvec_to_complex(size(vec)/2) 
 integer    :: i
  if(mod(size(vec),2)/=0) stop 'error realvec_to_complex :: vec size s odd' 
  do i=1,size(vec)/2
    realvec_to_complex(i)=CMPLX(vec(2*i-1),vec(2*i),kind=8)
  enddo
 end function

 !--------------------------!

 function complex_to_realvec(vec)
 implicit none
 complex(8) :: vec(:)
 real(8)    :: complex_to_realvec(2*size(vec))
 integer    :: i
   do i=1,size(vec)
     complex_to_realvec(2*i-1) =  real(vec(i))
     complex_to_realvec(2*i  ) = aimag(vec(i))
   enddo
 end function

 !--------------------------!

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
!**************************************************************************
!**************************************************************************

 !--------------------------!

 subroutine put_and_shift____(g,ii,xx)
 implicit none
 integer :: g(:),xx
 integer :: ii,i,j,k,l,siz
 siz=size(g)
 if(ii+1<=siz)then
  g(ii+1:siz)=g(ii:siz-1)
  g(ii)=xx
 endif
 end subroutine

 !--------------------------!

 subroutine put_and_shift___(g,ii,xx)
 implicit none
 complex(8) :: g(:),xx
 integer :: ii,i,j,k,l,siz
 siz=size(g)
 if(ii+1<=siz)then
  g(ii+1:siz)=g(ii:siz-1)
  g(ii)=xx
 endif
 end subroutine
 !--------------------------!

 subroutine put_and_shift__(g,ii,xx)
 implicit none
 real(4) :: g(:),xx
 integer :: ii,i,j,k,l,siz
 siz=size(g)
 if(ii+1<=siz)then
  g(ii+1:siz)=g(ii:siz-1)
  g(ii)=xx
 endif
 end subroutine

 !--------------------------!

 subroutine put_and_shift_(g,ii,xx)
 implicit none
 real(8) :: g(:),xx
 integer :: ii,i,j,k,l,siz
 siz=size(g)
 if(ii+1<=siz)then
  g(ii+1:siz)=g(ii:siz-1)
  g(ii)=xx 
 endif
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************

 integer function modi(j,kk)
 implicit none
 !1,2,3....kk 1 2 3... kk
 integer :: j,kk
  modi= mod(j,kk)  
  if(modi==0) modi=kk
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 real(8) function minloc_parabola(a)
 implicit none
 real(8) :: a(3)
  minloc_parabola=-a(2)/2.d0/a(3)
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 real(8) function minval_parabola(a)
 implicit none
  real(8) :: a(3)
  minval_parabola=pol(a,minloc_parabola(a)) 
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine shift_lin_pos_parabola(a)
 implicit none
 real(8) :: a(3)
   a(2) = sqrt(abs(a(1)*a(3)*4))
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 function root_equation(a)
 implicit none
 real(8) :: root_equation(2)
 real(8) :: a(3),dd

  if(abs(a(1))>error)then
   dd=abs(a(2)**2-4.d0*a(1)*a(3))
   if(abs(dd)>error)then
    root_equation(1)=(-a(2)+sqrt(dd))/2.d0/a(1)
    root_equation(2)=(-a(2)-sqrt(dd))/2.d0/a(1)
   endif
  endif

 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 logical function maxontop_array_d(x)
 implicit none
 real(8) :: x(:)
  maxontop_array_d=abs(maxval(abs(x))-maxval(x))<epsilonr
 end function

 logical function constantarray_d(x)
 implicit none
 real(8) :: x(:)
  constantarray_d=abs(maxval(x)-minval(x))<1.d-5
 end function

 logical function nullarray_d(x)
 implicit none
 real(8) :: x(:)
  nullarray_d=maxval(abs(x))<epsilonr
 end function

 logical function maxontop_array_r(x)
 implicit none
 real(4) :: x(:)
  maxontop_array_r=abs(maxval(abs(x))-maxval(x))<epsilonr
 end function

 logical function constantarray_r(x)
 implicit none
 real(4) :: x(:)
  constantarray_r=abs(maxval(x)-minval(x))<1.d-5
 end function

 logical function nullarray_r(x)
 implicit none
 real(4) :: x(:)
  nullarray_r=maxval(abs(x))<epsilonr
 end function

 logical function maxontop_array_d_(x)
 implicit none
 real(8) :: x(:,:)
  maxontop_array_d_=abs(maxval(abs(x))-maxval(x))<epsilonr
 end function

 logical function constantarray_d_(x)
 implicit none
 real(8) :: x(:,:)
  constantarray_d_=abs(maxval(x)-minval(x))<1.d-5
 end function

 logical function nullarray_d_(x)
 implicit none
 real(8) :: x(:,:)
  nullarray_d_=maxval(abs(x))<epsilonr
 end function

 logical function maxontop_array_r_(x)
 implicit none
 real(4) :: x(:,:)
  maxontop_array_r_=abs(maxval(abs(x))-maxval(x))<epsilonr
 end function

 logical function constantarray_r_(x)
 implicit none
 real(4) :: x(:,:)
  constantarray_r_=abs(maxval(x)-minval(x))<1.d-5
 end function

 logical function nullarray_r_(x)
 implicit none
 real(4) :: x(:,:)
  nullarray_r_=maxval(abs(x))<epsilonr
 end function

 logical function maxontop_array_d__(x)
 implicit none
 real(8) :: x(:,:,:)
  maxontop_array_d__=abs(maxval(abs(x))-maxval(x))<epsilonr
 end function

 logical function constantarray_d__(x)
 implicit none
 real(8) :: x(:,:,:)
  constantarray_d__=abs(maxval(x)-minval(x))<1.d-5
 end function

 logical function nullarray_d__(x)
 implicit none
 real(8) :: x(:,:,:)
  nullarray_d__=maxval(abs(x))<epsilonr
 end function

 logical function maxontop_array_r__(x)
 implicit none
 real(4) :: x(:,:,:)
  maxontop_array_r__=abs(maxval(abs(x))-maxval(x))<epsilonr
 end function

 logical function constantarray_r__(x)
 implicit none
 real(4) :: x(:,:,:)
  constantarray_r__=abs(maxval(x)-minval(x))<1.d-5
 end function

 logical function nullarray_r__(x)
 implicit none
 real(4) :: x(:,:,:)
  nullarray_r__=maxval(abs(x))<epsilonr
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

logical function same_array_rr(mat1,mat2)
implicit none
 real(8)  :: mat1(:,:),mat2(:,:)
 integer  :: i 

 interface
  integer function same_address(d1,d2)
   real(8) :: d1,d2
  end function
 end interface

 i=same_address(mat1(1,1),mat2(1,1))
 if(i==1) then
  same_array_rr=.true.
 else
  same_array_rr=.false.
 endif
return
end function

   !-----------------------------!

logical function same_array_r(mat1,mat2)
implicit none
 real(4)  :: mat1(:,:),mat2(:,:)
 integer  :: i
 interface
  integer function same_address_r(d1,d2)
   real(4) :: d1,d2
  end function
 end interface

 i=same_address_r(mat1(1,1),mat2(1,1))
 if(i==1) then
  same_array_r=.true.
 else
  same_array_r=.false.
 endif
return
end function

   !-----------------------------!

logical function same_array_cc(mat1,mat2)
implicit none
 complex(4)  :: mat1(:,:),mat2(:,:)
 integer     :: i
 interface
  integer function same_address_r(d1,d2)
   complex(4) :: d1,d2
  end function
 end interface

 i=same_address_r(mat1(1,1),mat2(1,1))
 if(i==1) then
  same_array_cc=.true.
 else
  same_array_cc=.false.
 endif
return
end function

   !-----------------------------!

logical function same_array_c(mat1,mat2)
implicit none
 complex(8)  :: mat1(:,:),mat2(:,:)
 integer  :: i
 interface
  integer function same_address(d1,d2)
   complex(8) :: d1,d2
  end function
 end interface

 i=same_address(mat1(1,1),mat2(1,1))
 if(i==1) then
  same_array_c=.true.
 else
  same_array_c=.false.
 endif
return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine reverse_array____(vec)
 implicit none
 real(8) :: vec(:)
 integer :: i,j,mid,k    
 k=size(vec)
 if(mod(k,2)==0)then
  mid=k/2          
 else
  mid=(k+1)/2
 endif
 do i=1,mid     
  call swap(vec(i),vec(k-i+1))
 enddo
 return 
 end subroutine

 subroutine reverse_array_(vec)
 implicit none
 complex(8) :: vec(:)
 integer    :: i,j,mid,k
 k=size(vec)
 if(mod(k,2)==0)then
  mid=k/2
 else
  mid=(k+1)/2
 endif
 do i=1,mid
  call swap(vec(i),vec(k-i+1))
 enddo
 return
 end subroutine

 subroutine reverse_array__(vec)
 implicit none
 real(4)    :: vec(:)
 integer    :: i,j,mid,k
 k=size(vec)
 if(mod(k,2)==0)then
  mid=k/2
 else
  mid=(k+1)/2
 endif
 do i=1,mid
  call swap(vec(i),vec(k-i+1))
 enddo
 return
 end subroutine

 subroutine reverse_array___(vec)
 implicit none
 integer(4) :: vec(:)
 integer    :: i,j,mid,k
 k=size(vec)
 if(mod(k,2)==0)then
  mid=k/2
 else
  mid=(k+1)/2
 endif
 do i=1,mid
  call swap(vec(i),vec(k-i+1))
 enddo
 return
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine remove_linear_contribution(x,f)
 implicit none
 real(8)   ::  x(:),f(:),slope,v(6),off
 integer   ::  i,j,k,l,m,siz

  siz=size(x)
  v= (/( (f(i+1)-f(i))/(x(i+1)-x(i)), i=siz-6,siz-1 )/) / 6.d0 
  slope=sum(v)
  v= (/( f(i), i=siz-6,siz-1 )/) / 6.d0
  off=sum(v)
  f = (/( f(i) - (off+slope*x(i)) , i=1,siz )/)

 return
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine erase_negr(vec)
 implicit none
  real :: vec(:)
  where(vec<0.)vec=0.
 end subroutine

 subroutine erase_negr_(vec)
 implicit none
  real :: vec(:,:)
  where(vec<0.)vec=0.
 end subroutine

 subroutine erase_negr__(vec)
 implicit none
  real :: vec(:,:,:)
  where(vec<0.)vec=0.
 end subroutine

 subroutine erase_negr___(vec)
 implicit none
  real :: vec(:,:,:,:)
  where(vec<0.)vec=0.
 end subroutine

 subroutine erase_negd(vec)
 implicit none
  real(8) :: vec(:)
  where(vec<0.)vec=0.
 end subroutine

 subroutine erase_negd_(vec)
 implicit none
  real(8) :: vec(:,:)
  where(vec<0.)vec=0.
 end subroutine

 subroutine erase_negd__(vec)
 implicit none
  real(8) :: vec(:,:,:)
  where(vec<0.)vec=0.
 end subroutine

 subroutine erase_negd___(vec)
 implicit none
  real(8) :: vec(:,:,:,:)
  where(vec<0.)vec=0.
 end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine erase_noise_(x,f,ferr,aa)
 implicit none
 real(8)          :: x(:),f(:),ferr(:),aa
 real(8)          :: u,v,w,maxerr,thres,ddm,dd2
 real(4)          :: ax,ay,bx,by,cx,cy,tempx,tempy,a,b,c
 integer          :: i,j,k,ierror
 logical          :: mask(size(f))

 mask=.true.; maxerr=maxval(ferr); thres=maxerr*aa

 if(maxval(ferr)<1.d-4) return

  if(size(x)<3) return
  do i=3,(size(ferr)-1)
   if(ferr(i)>thres.and.ferr(i-1)<ferr(i).and.ferr(i-2)<ferr(i).and.ferr(i+1)<ferr(i)) then
     f(i)=dumpi(i)
     mask(i)=.false.
     ferr(i)=(ferr(i-2)+ferr(i-1)+ferr(i+1))/3.d0
   endif
  enddo 
  do i=3,(size(ferr)-2)
   if(ferr(i)>thres.and.ferr(i-1)<ferr(i).and.ferr(i+2)<ferr(i).and.ferr(i+1)<ferr(i)) then
    if(mask(i+1).and.mask(i-1).and.mask(i+2))then
     f(i)=dumpi_for(i)
     mask(i)=.true.
     ferr(i)=(ferr(i-2)+ferr(i-1)+ferr(i+1))/3.d0
    endif
   endif
  enddo

  do i=2,size(x)-1
   if(mask(i-1).and.mask(i).and.mask(i+1))then
    if(ferr(i)<ferr(i-1).and.ferr(i)<ferr(i+1).and.ferr(i)>thres)then
     f(i)=(f(i-1)+f(i+1))/2.d0
     ferr(i)=(ferr(i-1)+ferr(i+1))/2.d0
     mask(i)=.false.
    endif 
   endif
  enddo


 contains

  !----------!
  !----------!

 real(8) function dumpi(i)
 integer :: i
   ax=x(i-2)
   ay=f(i-2)
   bx=x(i-1)
   by=f(i-1)
   cx=x(i+1)
   cy=f(i+1)
   call parabola_ex2( ax, ay, bx, by, cx , cy, tempx, tempy, a,b,c, ierror )
   dumpi=c + b*x(i) + a*(x(i)**2.d0)
 end function

 real(8) function dumpi_for(i)
 integer :: i
   ax=x(i-1)
   ay=f(i-1)
   bx=x(i+1)
   by=f(i+1)
   cx=x(i+2)
   cy=f(i+2)
   call parabola_ex2( ax, ay, bx, by, cx , cy, tempx, tempy, a,b,c, ierror )
   dumpi_for=c + b*x(i) + a*(x(i)**2.d0)
 end function


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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 function polr(P,x)
 implicit none
   real    :: polr,P(:),x
   integer :: i
  polr=P(1)
  do i=2,size(P)
   polr=polr+P(i)*(x**float(i-1))
  enddo
 end function

 function pold(P,x)
 implicit none
   real(8) :: pold,P(:),x
   integer :: i
  pold=P(1)
  do i=2,size(P)
   pold=pold+P(i)*(x**dble(i-1))
  enddo
 end function

 function poli(P,x)
 implicit none
   integer :: poli,P(:),x
   integer :: i
  poli=P(1)
  do i=2,size(P)
   poli=poli+P(i)*(x**(i-1))
  enddo
 end function

 function polc(P,x)
 implicit none
   complex(8) :: polc,P(:),x
   integer    :: i
  polc=P(1)
  do i=2,size(P)
   polc=polc+P(i)*(x**dble(i-1))
  enddo
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine print_machine_limitations
 implicit none
  write(*,*) ' radix       :  ',  radix ( 1 )
  write(*,*) ' digits      :  ',  digits ( 1.d0 )
  write(*,*) ' maxexponent :  ',  maxexponent (1.d0)
  write(*,*) ' minexponent :  ',  minexponent (1.d0)
  write(*,*) ' tiny        :  ',  tiny (1.d0)
  write(*,*) ' huge        :  ',  huge (1.d0)
  write(*,*) ' epsilon     :  ',  epsilon (1.d0)
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

integer function maxloci_c(mat)
complex(8) :: mat(:)
integer    :: u(1)
 u=maxloc(abs(mat))
 maxloci_c=u(1)
end function

integer function maxloci_r(mat)
real(8)     :: mat(:)
integer    :: u(1)
 u=maxloc(mat)
 maxloci_r=u(1)
end function

integer function maxloci_rr(mat)
real(4)    :: mat(:)
integer    :: u(1)
 u=maxloc(mat)
 maxloci_rr=u(1)
end function

integer function maxloci_i(mat)
integer(4)  :: mat(:)
integer    :: u(1)
 u=maxloc(mat)
 maxloci_i=u(1)
end function


integer function minloci_c(mat)
complex(8) :: mat(:)
integer    :: u(1)
 u=minloc(abs(mat))
 minloci_c=u(1)
end function

integer function minloci_r(mat)
real(8)     :: mat(:)
integer    :: u(1)
 u=minloc(mat)
 minloci_r=u(1)
end function

integer function minloci_rr(mat)
real(4)    :: mat(:)
integer    :: u(1)
 u=minloc(mat)
 minloci_rr=u(1)
end function

integer function minloci_i(mat)
integer(4)  :: mat(:)
integer    :: u(1)
 u=minloc(mat)
 minloci_i=u(1)
end function

     !-------------------------------!

real(8) function maxvalcondrr(YR,XR,a1,a2)
implicit none
real(8) :: YR(:),XR(:),a1,a2,minr
integer :: i,j,k
minr=-1.d20
do i=1,size(XR)
 if(XR(i)<=a2.and.XR(i)>=a1)then
   if(YR(i)>minr)then
     minr=YR(i)
   endif
 endif
enddo
maxvalcondrr=minr
end function

     !-------------------------------!

real(8) function minvalcondrr(YR,XR,a1,a2)
implicit none
real(8) :: YR(:),XR(:),a1,a2,minr
integer :: i,j,k
minr=1.d20
do i=1,size(XR)
 if(XR(i)<=a2.and.XR(i)>=a1)then
   if(YR(i)<minr)then
     minr=YR(i)
   endif
 endif
enddo
minvalcondrr=minr
end function

     !-------------------------------!

real(4) function maxvalcondr(YR,XR,a1,a2)
implicit none
real(4) :: YR(:),XR(:),a1,a2,minr
integer :: i,j,k
minr=-1.e10
do i=1,size(XR)
 if(XR(i)<=a2.and.XR(i)>=a1)then
   if(YR(i)>minr)then
     minr=YR(i)
   endif
 endif
enddo
maxvalcondr=minr
end function

     !-------------------------------!

real(4) function minvalcondr(YR,XR,a1,a2)
implicit none
real(4) :: YR(:),XR(:),a1,a2,minr
integer :: i,j,k
minr=1.e10
do i=1,size(XR)
 if(XR(i)<=a2.and.XR(i)>=a1)then
   if(YR(i)<minr)then
     minr=YR(i)
   endif
 endif
enddo
minvalcondr=minr
end function

     !-------------------------------!

real(4) function maxvalcondr_(YR,XR,a1,a2)
implicit none
real(4) :: YR(:,:),XR(:,:),a1,a2,minr
integer :: i,j,k
logical :: test

minr=-1.e10
do i=1,size(XR(1,:))
 test=.true.
 do j=1,size(XR(:,1))
 if(XR(j,i)>a2.or.XR(j,i)<a1)then
  test=.false.
  exit
 endif
 enddo
 if(test)then
   if(maxval(YR(:,i))>minr)then
     minr=maxval(YR(:,i))
   endif
 endif
enddo

maxvalcondr_=minr

end function

     !-------------------------------!

real(4) function minvalcondr_(YR,XR,a1,a2)
implicit none
real(4) :: YR(:,:),XR(:,:),a1,a2,minr
integer :: i,j,k
logical :: test

minr=1.e10
do i=1,size(XR(1,:))
 test=.true.
 do j=1,size(XR(:,1)) 
 if(XR(j,i)>a2.or.XR(j,i)<a1)then
  test=.false.
  exit
 endif
 enddo
 if(test)then
   if(minval(YR(:,i))<minr)then
     minr=minval(YR(:,i))
   endif
 endif
enddo

minvalcondr_=minr

end function

     !-------------------------------!

logical function none_zero_term_in_there_i(tab)
implicit none
integer :: i,j,k,l,m,siz1,tab(:)
 none_zero_term_in_there_i=.false.
 do i=1,size(tab)
  if(abs(tab(i))>0)then
   none_zero_term_in_there_i=.true.
   return
  endif
 enddo
end function

     !-------------------------------!

logical function none_zero_term_in_there_r(tab)
implicit none
integer :: i,j,k,l,m,siz1
real(8)  :: tab(:)
 none_zero_term_in_there_r=.false.
 do i=1,size(tab)
  if(abs(tab(i))>1.d-4)then
   none_zero_term_in_there_r=.true.
   return
  endif
 enddo
end function

     !-------------------------------!

logical function none_zero_term_in_there_c(tab)
implicit none
integer     :: i,j,k,l,m,siz1
complex(8)  :: tab(:)
 none_zero_term_in_there_c=.false.
 do i=1,size(tab)
  if(abs(tab(i))>1.d-4)then
   none_zero_term_in_there_c=.true.
   return
  endif
 enddo
end function

     !-------------------------------!

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
!**************************************************************************
!**************************************************************************

elemental complex(8) function csqrt(rin)
implicit none
real(8),intent(in)     :: rin
     if(rin>-0.000001d0)Then
       csqrt=sqrt(abs(rin))
     else
       csqrt=cmplx(0.d0,1.d0,kind=8)*sqrt(abs(rin))
     endif
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine rep_re(xx,aa)
implicit none
complex(8) :: xx,yy
real(8)     :: aa
 yy=cmplx(0.d0,1.d0,kind=8)*aimag(xx)
 yy=yy + aa
 xx=yy
end subroutine

 !---------!

subroutine rep_im(xx,aa)
implicit none
complex(8) :: xx,yy
real(8)     :: aa
 yy=real(xx)
 yy=yy + aa*cmplx(0.d0,1.d0,kind=8)
 xx=yy
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
!**************************************************************************
!**************************************************************************

real(8) function theta_func(rr,a)
implicit none
real(8) :: rr
real(8),optional :: a

if(.not.present(a))then
if(rr<0.)then
 theta_func=1.
else
 theta_func=0.
endif
else
if(rr<a)then
 theta_func=1.
else
 theta_func=0.
endif
endif

end function

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
!**************************************************************************
!**************************************************************************

function per_c(rr)
implicit none
complex(8)  :: rr
real(8)      :: per_c,ss,tt
 per_c=0.
 ss=abs(aimag(rr))
 tt=abs(real(rr))
 if(ss>1.d-12) per_c=ss/tt
return
end function

 !---------------!
 !---------------!
 !---------------!

function per_c_vec(rr)
implicit none
complex(8)  :: rr(:)
integer     :: i
real(8)      :: per_c_vec(size(rr)),ss,tt
 do i=1,size(rr)
  per_c_vec(i)=0.
  ss=abs(aimag(rr(i)))
  tt=abs(real(rr(i)))
  if(ss>1.d-12) per_c_vec(i)=ss/tt
 enddo
return
end function

 !---------------!
 !---------------!
 !---------------!

function per_c_mat(rr)
implicit none
complex(8)  :: rr(:,:)
integer     :: i,siz1,siz2,j
real(8)      :: per_c_mat(size(rr(:,1)),size(rr(1,:))),ss,tt

siz1=size(rr(:,1))
siz2=size(rr(1,:))

 do i=1,siz1
 do j=1,siz2
  per_c_mat(i,j)=0.
  ss=abs(aimag(rr(i,j)))
  tt=abs(real(rr(i,j)))
  if(ss>1.d-12) per_c_mat(i,j)=ss/tt
 enddo
 enddo

return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function prod_r(rr)
implicit none
real(8)  :: rr(:,:)
real(8)  :: prod_r
integer :: s1,s2,i,j,k,l,m

prod_r=1.
s1=size(rr(:,1)); s2=size(rr(1,:))

do i=1,s1
do j=1,s2
 prod_r=prod_r * (rr(i,j))
enddo
enddo

return
end function

 !---------------!
 !---------------!
 !---------------!

function prod_i(rr)
implicit none
integer  :: rr(:,:)
integer  :: prod_i
integer    :: s1,s2,i,j,k,l,m

prod_i=1
s1=size(rr(:,1)); s2=size(rr(1,:))

do i=1,s1
 do j=1,s2
  prod_i=prod_i * (rr(i,j))
 enddo
enddo

return
end function

 !---------------!
 !---------------!
 !---------------!

function prod_c(rr)
implicit none
complex(8)  :: rr(:,:)
complex(8)  :: prod_c
integer     :: s1,s2,i,j,k,l,m

 prod_c=1.
 s1=size(rr(:,1));s2=size(rr(1,:))

 do i=1,s1
  do j=1,s2
   prod_c=prod_c * (rr(i,j))
  enddo
 enddo

return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function minloc_with_test_r(array,arraytest,value)
implicit none
real(8)  :: array(:),arraytemp(size(array))
integer :: arraytest(:),value,minloc_with_test_r
integer :: i,j,k,l,m,n,u(1),siz1
 siz1=size(array);arraytemp=array;k=0
 do 
  u=minloc(arraytemp)
  if(arraytest(u(1))==value)then
   exit
  else
   k=k+1   
   arraytemp(u(1))=1.d20
  endif 
  if(k==siz1) then
   minloc_with_test_r=1
   return 
  endif
 enddo
 minloc_with_test_r=u(1)
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function minloc_with_test_c(array,arraytest,value)
implicit none
complex(8)  :: array(:),arraytemp(size(array))
integer     :: arraytest(:),value,minloc_with_test_c
integer     :: i,j,k,l,m,n,u(1),siz1
 siz1=size(array);arraytemp=array;k=0
 do
  u=minloc(real(arraytemp))
  if(arraytest(u(1))==value)then
   exit
  else
   k=k+1
   arraytemp(u(1))=1.d20
  endif
  if(k==siz1) then
   minloc_with_test_c=1
   return
  endif
 enddo
 minloc_with_test_c=u(1)
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

logical function too_large_r(a)
real(8) :: a
too_large_r=.false.
 if(a>MAX_REAL)then
  too_large_r=.true.  
 endif
end function

logical function too_large_c(a)
complex(8) :: a
too_large_c=.false.
 if(real(a)>MAX_REAL.or.aimag(a)>MAX_REAL)then
  too_large_c=.true.
 endif
end function

logical function too_large_c_vec(a)
complex(8) :: a(:)
too_large_c_vec=.false.
 if(maxval(abs(real(a)))>MAX_REAL.or.maxval(abs(aimag(a)))>MAX_REAL)then
  too_large_c_vec=.true.
 endif
end function

logical function too_large_c_mat(a)
complex(8) :: a(:,:)
 too_large_c_mat=.false.
 if(maxval(abs(real(a)))>MAX_REAL.or.maxval(abs(aimag(a)))>MAX_REAL)then
  too_large_c_mat=.true.
 endif
end function

logical function too_large_i(a)
integer :: a
too_large_i=.false.
 if(i>MAX_INT)then
  too_large_i=.true.
 endif
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 logical function error_r(b,a,threshold)
 implicit none
 !b=observable
 !a=error
 !threshold : in [percent]
 real(8),intent(in) :: a,b
 real(8),intent(in) :: threshold
 real(8) :: ratio

  error_r=.true.

  if(abs(b)<1.d-9) then
   error_r=.false.
   return
  endif

  ratio=abs(a)/abs(b)*100.d0
  if(ratio>threshold) error_r=.false.

  return
 end function 

 logical function error_c(b,a,threshold)
 implicit none
 !b=observable
 !a=error
 !threshold : in [percent]
 complex(8),intent(in) :: a,b
 real(8),intent(in) :: threshold
 real(8) :: ratio

  error_c=.true.

  if(abs(b)<1.d-9) then
   error_c=.false.
   return
  endif

  ratio=abs(a)/abs(b)*100.d0
  if(ratio>threshold) error_c=.false.

  return
 end function

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

 function pack_r(mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,vv(2)
 real(8),intent(in) :: mat(:,:)
 real(8)            :: pack_r(size(mat(:,1))*size(mat(1,:)))
  siz1=size(mat,1); siz2=size(mat,2)
  count=0
  do i=1,siz1 
   do j=1,siz2
    count=count+1
    pack_r(count)=mat(i,j)
   enddo
  enddo
 end function

 subroutine pack_r_(pack,mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,vv(2)
 real(8),intent(in) :: mat(:,:)
 real(8),intent(out):: pack(size(mat(:,1))*size(mat(1,:)))
  siz1=size(mat,1); siz2=size(mat,2)
  count=0
  do i=1,siz1
   do j=1,siz2
    count=count+1
    pack(count)=mat(i,j)
   enddo
  enddo
 end subroutine

 function pack_c(mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer               :: siz1,siz2,i,j,k,l,m,count,vv(2)
 complex(8),intent(in) :: mat(:,:)
 complex(8)            :: pack_c(size(mat(:,1))*size(mat(1,:)))
  siz1=size(mat(:,1))
  siz2=size(mat(1,:))
  count=0
  do i=1,siz1 
   do j=1,siz2
    count=count+1
    pack_c(count)=mat(i,j)
   enddo
  enddo
 end function

 subroutine pack_c_(pack,mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,vv(2)
 complex(8),intent(in) :: mat(:,:)
 complex(8),intent(out):: pack(size(mat(:,1))*size(mat(1,:)))
  siz1=size(mat,1); siz2=size(mat,2)
  count=0
  do i=1,siz1
   do j=1,siz2
    count=count+1
    pack(count)=mat(i,j)
   enddo
  enddo
 end subroutine


  function pack_i(mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,vv(2)
 integer,intent(in) :: mat(:,:)
 integer            :: pack_i(size(mat(:,1))*size(mat(1,:)))
  siz1=size(mat(:,1))
  siz2=size(mat(1,:))
  count=0
  do i=1,siz1 
   do j=1,siz2
    count=count+1
    pack_i(count)=mat(i,j)
   enddo
  enddo
 end function

 subroutine pack_i_(pack,mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer                :: siz1,siz2,i,j,k,l,m,count,vv(2)
 integer(4),intent(in)  :: mat(:,:)
 integer(4),intent(out) :: pack(size(mat(:,1))*size(mat(1,:)))
  siz1=size(mat,1); siz2=size(mat,2)
  count=0
  do i=1,siz1
   do j=1,siz2
    count=count+1
    pack(count)=mat(i,j)
   enddo
  enddo
 end subroutine


 function pack_s(mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer         :: siz1,siz2,i,j,k,l,m,count,vv(2)
 real,intent(in) :: mat(:,:)
 real(4)         :: pack_s(size(mat(:,1))*size(mat(1,:)))
  siz1=size(mat(:,1))
  siz2=size(mat(1,:))
  count=0
  do i=1,siz1 
   do j=1,siz2
    count=count+1
    pack_s(count)=mat(i,j)
   enddo
  enddo
 end function

 subroutine pack_s_(pack,mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,vv(2)
 real(4),intent(in) :: mat(:,:)
 real(4),intent(out):: pack(size(mat(:,1))*size(mat(1,:)))
  siz1=size(mat,1); siz2=size(mat,2)
  count=0
  do i=1,siz1
   do j=1,siz2
    count=count+1
    pack(count)=mat(i,j)
   enddo
  enddo
 end subroutine


 function unpack_r(mat,ligne)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,ligne,vv(2)
 real(8),intent(in) :: mat(:)
 real(8)            :: unpack_r(size(mat)/ligne,ligne)
  siz1=size(mat(:))
  count=0
  do i=1,siz1 
    vv=k_to_ij(ligne,i)
    unpack_r(vv(1),vv(2))=mat(i)
  enddo
  return
 end function


 subroutine unpack_r_(unpack,mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,ligne,vv(2)
 real(8),intent(in) :: mat(:)
 real(8),intent(out):: unpack(:,:)
  ligne=size(unpack,2)
  if(mod(size(mat),ligne)/=0) stop 'unpack error, size do not match'
  if(size(unpack,1)/=size(mat)/ligne) then
      write(*,*) size(unpack,1),size(mat),ligne
      stop ' unpack_r error, size problem'
  endif
  siz1=size(mat(:))
  count=0
  do i=1,siz1
    vv=k_to_ij(ligne,i)
    unpack(vv(1),vv(2))=mat(i)
  enddo
  return
 end subroutine


 function unpack_c(mat,ligne)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer               :: siz1,siz2,i,j,k,l,m,count,ligne,vv(2)
 complex(8),intent(in) :: mat(:)
 complex(8)            :: unpack_c(size(mat)/ligne,ligne)
  siz1=size(mat(:))
  count=0
  do i=1,siz1 
    vv=k_to_ij(ligne,i)
    unpack_c(vv(1),vv(2))=mat(i)
  enddo
  return
 end function

 subroutine unpack_c_(unpack,mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,ligne,vv(2)
 complex(8),intent(in) :: mat(:)
 complex(8),intent(out):: unpack(:,:)
  ligne=size(unpack,2)
  if(mod(size(mat),ligne)/=0) stop 'unpack error, size do not match'
  if(size(unpack,1)/=size(mat)/ligne) then
      write(*,*) size(unpack,1),size(mat),ligne
      stop ' unpack_c error, size problem'
  endif
  siz1=size(mat(:))
  count=0
  do i=1,siz1
    vv=k_to_ij(ligne,i)
    unpack(vv(1),vv(2))=mat(i)
  enddo
  return
 end subroutine


 function unpack_s(mat,ligne)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer         :: siz1,siz2,i,j,k,l,m,count,ligne,vv(2)
 real,intent(in) :: mat(:)
 real            :: unpack_s(size(mat)/ligne,ligne)
  siz1=size(mat(:))
  count=0
  do i=1,siz1 
    vv=k_to_ij(ligne,i)
    unpack_s(vv(1),vv(2))=mat(i)
  enddo
  return
 end function

 subroutine unpack_s_(unpack,mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,ligne,vv(2)
 real(4),intent(in) :: mat(:)
 real(4),intent(out):: unpack(:,:)
  ligne=size(unpack,2)
  if(mod(size(mat),ligne)/=0) stop 'unpack error, size do not match'
  if(size(unpack,1)/=size(mat)/ligne) then
      write(*,*) size(unpack,1),size(mat),ligne 
      stop ' unpack_s error, size problem'
  endif
  siz1=size(mat(:))
  count=0
  do i=1,siz1
    vv=k_to_ij(ligne,i)
    unpack(vv(1),vv(2))=mat(i)
  enddo
  return
 end subroutine

 function unpack_i(mat,ligne)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,ligne,vv(2)
 integer,intent(in) :: mat(:)
 integer            :: unpack_i(size(mat)/ligne,ligne)
  siz1=size(mat(:))
  count=0
  do i=1,siz1
    vv=k_to_ij(ligne,i)
    unpack_i(vv(1),vv(2))=mat(i)
  enddo
  return
 end function

 subroutine unpack_i_(unpack,mat)
 !pack matrice ligne par ligne dans vecteur
 implicit none
 integer            :: siz1,siz2,i,j,k,l,m,count,ligne,vv(2)
 integer(4),intent(in) :: mat(:)
 integer(4),intent(out):: unpack(:,:)
  ligne=size(unpack,2)
  if(mod(size(mat),ligne)/=0) stop 'unpack error, size do not match'
  if(size(unpack,1)/=size(mat)/ligne) then
    write(*,*) size(unpack,1),size(mat),ligne  
    stop ' unpack_i error, size problem'
  endif
  siz1=size(mat(:))
  count=0
  do i=1,siz1
    vv=k_to_ij(ligne,i)
    unpack(vv(1),vv(2))=mat(i)
  enddo
  return
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

logical function ISNANR(x)
 real(8) :: x
 ISNANR=disnan(x)
end function

logical function ISNANC(x)
 complex(8) :: x
 ISNANC=disnan(REAL(x))
 if(ISNANC) return
 ISNANC=disnan(aimag(x))
end function

logical function ISINFR(x)
 real(8) :: x
 ISINFR=disinf(x)
end function

logical function ISINFC(x)
 complex(8) :: x
 ISINFC=disinf(REAL(x))
 if(ISINFC) return
 ISINFC=disinf(aimag(x))
end function

logical function ISNANS(x)
 real(4) :: x
 ISNANS=isnan(x)
end function

logical function ISINFS(x)
 real(4) :: x
 ISINFS=isinf(x)
end function


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function isqrt(ii)
integer :: ii,isqrt
 isqrt=NINT ( SQRT( dble(ii) ) )
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function comparearray_s(array,rrr)
implicit none
type(arrayc),intent(in) :: array
logical :: comparearray_s
real(8),intent(in) :: rrr
if(real(array%subarray(1))<rrr) then
 comparearray_s=.true.
else
 comparearray_s=.false.
endif
return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function comparearray_g(array,rrr)
implicit none
type(arrayc),intent(in) :: array
logical :: comparearray_g
real(8),intent(in) :: rrr
if(real(array%subarray(1))>rrr) then
 comparearray_g=.true.
else
 comparearray_g=.false.
endif
return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function comparearray_sc(array,rrr)
implicit none
type(arrayc),intent(in) :: array
logical :: comparearray_sc
complex(8),intent(in) :: rrr
if(real(array%subarray(1))<real(rrr)) then
 comparearray_sc=.true.
else
 comparearray_sc=.false.
endif
return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function comparearray_gc(array,rrr)
implicit none
type(arrayc),intent(in) :: array
logical :: comparearray_gc
complex(8),intent(in) :: rrr
if(real(array%subarray(1))>real(rrr)) then
 comparearray_gc=.true.
else
 comparearray_gc=.false.
endif
return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function mult_c(array,rrr)
implicit none
type(arrayc),intent(in) :: array
type(arrayc) :: mult_c
complex(8),intent(in) :: rrr
 mult_c%subarray=array%subarray(:)*rrr
return
end function

!************************************************
!************************************************

function div_c(array,rrr)
implicit none
type(arrayc),intent(in) :: array
type(arrayc) :: div_c
complex(8),intent(in) :: rrr
 div_c%subarray=array%subarray(:)/rrr
return
end function

!************************************************
!************************************************

function add_c(array,rrr)
implicit none
type(arrayc),intent(in) :: array
type(arrayc) :: add_c
complex(8),intent(in) :: rrr
 add_c%subarray=array%subarray(:)+rrr
return
end function

!************************************************
!************************************************

function sub_c(array,rrr)
implicit none
type(arrayc),intent(in) :: array
type(arrayc) :: sub_c
complex(8),intent(in) :: rrr
 sub_c%subarray=array%subarray(:)-rrr
return
end function

!************************************************
!************************************************

function mult_r(array,rrr)
implicit none
type(arrayc),intent(in) :: array
type(arrayc)            :: mult_r
real(8),intent(in)      :: rrr
 mult_r%subarray=array%subarray(:)*rrr
return
end function

!************************************************
!************************************************

function div_r(array,rrr)
implicit none
type(arrayc),intent(in) :: array
type(arrayc) :: div_r
real(8),intent(in) :: rrr
 div_r%subarray=array%subarray(:)/rrr
return
end function

!************************************************
!************************************************

function add_r(array,rrr)
implicit none
type(arrayc),intent(in) :: array
type(arrayc) :: add_r
real(8),intent(in) :: rrr
 add_r%subarray=array%subarray(:)+rrr
return
end function

!************************************************
!************************************************

function sub_r(array,rrr)
implicit none
type(arrayc),intent(in) :: array
type(arrayc) :: sub_r
real(8),intent(in) :: rrr
 sub_r%subarray=array%subarray(:)-rrr
return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

      SUBROUTINE ZSWAP__(N,ZX,INCX,ZY,INCY)
      INTEGER INCX,INCY,N
      COMPLEX(16) ZX(*),ZY(*)
      COMPLEX(16) ZTEMP
      INTEGER I,IX,IY
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZTEMP = ZX(IX)
          ZX(IX) = ZY(IY)
          ZY(IY) = ZTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
          ZTEMP = ZX(I)
          ZX(I) = ZY(I)
          ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END subroutine

!*****************

      SUBROUTINE ZSWAP_(N,ZX,INCX,ZY,INCY)
      INTEGER     :: INCX,INCY,N
      COMPLEX(8)  :: ZX(*),ZY(*)
      COMPLEX(8)  :: ZTEMP
      INTEGER     :: I,IX,IY

      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZTEMP = ZX(IX)
          ZX(IX) = ZY(IY)
          ZY(IY) = ZTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
   20 DO 30 I = 1,N
          ZTEMP = ZX(I)
          ZX(I) = ZY(I)
          ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END subroutine

!*****************
elemental subroutine swapi(a,b)
integer,intent(inout) :: a,b
integer :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swaps(a,b)
real(4),intent(inout) :: a,b
real(4) :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swapr(a,b)
real(8),intent(inout) :: a,b
real(8) :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swapqr(a,b)
real(16),intent(inout) :: a,b
real(16) :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swapqc(a,b)
complex(16),intent(inout) :: a,b
complex(16) :: c
c=a
a=b
b=c
end subroutine
!*****************
subroutine swapa(a,b)
character*(*),intent(inout) :: a,b
character*20000:: c
integer :: maxab
maxab=LEN_TRIM(a)
if(LEN_TRIM(b)>maxab) maxab=LEN_TRIM(b)
c(1:maxab)=a(1:maxab)
a(1:maxab)=b(1:maxab)
b(1:maxab)=c(1:maxab)
end subroutine
!*****************
elemental subroutine swapc(a,b)
complex(8),intent(inout) :: a,b
complex(8) :: c
c=a
a=b
b=c
end subroutine
!*****************
elemental subroutine swapcs(a,b)
complex(4),intent(inout) :: a,b
complex(4) :: c
c=a
a=b
b=c
end subroutine
!*****************

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine rerase_divergence_scal(vec)
implicit none
real(8)::vec
integer ::i
if(disnan(vec).or.disinf(vec)) vec=0
end subroutine

subroutine rrerase_divergence_scal(vec)
implicit none
real::vec
integer ::i
if(disnan(dble(vec)).or.disinf(dble(vec))) vec=0
end subroutine

subroutine ierase_divergence_scal(vec)
implicit none
integer::vec
integer ::i
if(disnan(dble(vec)).or.disinf(dble(vec))) vec=0
end subroutine

subroutine cerase_divergence_scal(vec)
implicit none
complex(8)::vec
integer ::i
if(disnan(real(vec)).or.disnan(aimag(vec)).or.disinf(real(vec)).or.disinf(aimag(vec))) vec=0
end subroutine

subroutine rerase_divergence_vec(vec)
implicit none
real(8)::vec(:)
integer ::i
do i=1,size(vec)
 if(disnan(vec(i)).or.disinf(vec(i))) vec(i)=0
enddo
end subroutine

subroutine rrerase_divergence_vec(vec)
implicit none
real::vec(:)
integer ::i
do i=1,size(vec)
 if(disnan(dble(vec(i))).or.disinf(dble(vec(i)))) vec(i)=0
enddo
end subroutine

subroutine rerase_divergence_mat(mat)
implicit none
real(8)::mat(:,:)
integer ::i,j
do i=1,size(mat(:,1))
 do j=1,size(mat(1,:))
  if(disnan(mat(i,j)).or.disinf(mat(i,j))) mat(i,j)=0
 enddo
enddo
end subroutine

subroutine rrerase_divergence_mat(mat)
implicit none
real::mat(:,:)
integer ::i,j
do i=1,size(mat(:,1))
 do j=1,size(mat(1,:))
  if(disnan(dble(mat(i,j))).or.disinf(dble(mat(i,j)))) mat(i,j)=0
 enddo
enddo
end subroutine

subroutine ierase_divergence_vec(vec)
implicit none
integer::vec(:)
integer ::i
do i=1,size(vec)
 if(disnan(dble(vec(i))).or.disinf(dble(vec(i)))) vec(i)=0
enddo
end subroutine

subroutine ierase_divergence_mat(mat)
implicit none
integer::mat(:,:)
integer ::i,j
do i=1,size(mat(:,1))
 do j=1,size(mat(1,:))
  if(disnan(dble(mat(i,j))).or.disinf(dble(mat(i,j)))) mat(i,j)=0
 enddo
enddo
end subroutine

subroutine cerase_divergence_vec(vec)
implicit none
complex(8)::vec(:)
integer ::i
do i=1,size(vec)
 if(disnan(real(vec(i))).or.disnan(aimag(vec(i))).or.disinf(real(vec(i))).or.disinf(aimag(vec(i)))) vec(i)=0.
enddo
end subroutine

subroutine cerase_divergence_mat(mat)
implicit none
complex(8)::mat(:,:)
integer ::i,j
do i=1,size(mat(:,1))
 do j=1,size(mat(1,:))
  if(disnan(real(mat(i,j))).or.disnan(aimag(mat(i,j))).or.disinf(real(mat(i,j))).or.disinf(aimag(mat(i,j)))) mat(i,j)=0.
 enddo
enddo
end subroutine

subroutine rerase_divergence_mat_(mat)
implicit none
real(8) :: mat(:,:,:)
integer :: i1,i2,i3,i4,i5,i6,s(3)
s=shape(mat)
do i1=1,s(1)
 do i2=1,s(2)
  do i3=1,s(3)
   if(disnan(mat(i1,i2,i3)).or.disinf(mat(i1,i2,i3))) mat(i1,i2,i3)=0.
  enddo
 enddo
enddo
return
end subroutine

subroutine rerase_divergence_mat__(mat)
implicit none
real(8) :: mat(:,:,:,:)
integer :: i1,i2,i3,i4,i5,i6,s(4)
s=shape(mat)
do i1=1,s(1)
 do i2=1,s(2)
  do i3=1,s(3)
   do i4=1,s(4)
   if(disnan(mat(i1,i2,i3,i4)).or.disinf(mat(i1,i2,i3,i4))) mat(i1,i2,i3,i4)=0.
   enddo
  enddo
 enddo
enddo
return
end subroutine

subroutine rerase_divergence_mat___(mat)
implicit none
real(8) :: mat(:,:,:,:,:)
integer :: i1,i2,i3,i4,i5,i6,s(5)
s=shape(mat)
do i1=1,s(1)
 do i2=1,s(2)
  do i3=1,s(3)
   do i4=1,s(4)
    do i5=1,s(5)
     if(disnan(mat(i1,i2,i3,i4,i5)).or.disinf(mat(i1,i2,i3,i4,i5))) mat(i1,i2,i3,i4,i5)=0.
    enddo
   enddo
  enddo
 enddo
enddo 
return
end subroutine

subroutine rerase_divergence_mat____(mat)
implicit none
real(8) :: mat(:,:,:,:,:,:)
integer :: i1,i2,i3,i4,i5,i6,s(6)
s=shape(mat)
do i1=1,s(1)
 do i2=1,s(2)
  do i3=1,s(3)
   do i4=1,s(4)
    do i5=1,s(5)
     do i6=1,s(6)
      if(disnan(mat(i1,i2,i3,i4,i5,i6)).or.disinf(mat(i1,i2,i3,i4,i5,i6))) mat(i1,i2,i3,i4,i5,i6)=0.
     enddo
    enddo
   enddo
  enddo
 enddo
enddo
return
end subroutine

subroutine cerase_divergence_mat_(mat)
implicit none
complex(8) :: mat(:,:,:)
integer :: i1,i2,i3,i4,i5,i6,s(3)
s=shape(mat)
do i1=1,s(1)
 do i2=1,s(2)
  do i3=1,s(3)
   if(ISNAN_TEST(mat(i1,i2,i3)).or.ISINF_TEST(mat(i1,i2,i3))) mat(i1,i2,i3)=0.
  enddo
 enddo
enddo
return
end subroutine

subroutine cerase_divergence_mat__(mat)
implicit none
complex(8) :: mat(:,:,:,:)
integer :: i1,i2,i3,i4,i5,i6,s(4)
s=shape(mat)
do i1=1,s(1)
 do i2=1,s(2)
  do i3=1,s(3)
   do i4=1,s(4)
   if(ISNAN_TEST(mat(i1,i2,i3,i4)).or.ISINF_TEST(mat(i1,i2,i3,i4))) mat(i1,i2,i3,i4)=0.
   enddo
  enddo
 enddo
enddo
return
end subroutine

subroutine cerase_divergence_mat___(mat)
implicit none
complex(8) :: mat(:,:,:,:,:)
integer :: i1,i2,i3,i4,i5,i6,s(5)
s=shape(mat)
do i1=1,s(1)
 do i2=1,s(2)
  do i3=1,s(3)
   do i4=1,s(4)
    do i5=1,s(5)
     if(ISNAN_TEST(mat(i1,i2,i3,i4,i5)).or.ISINF_TEST(mat(i1,i2,i3,i4,i5))) mat(i1,i2,i3,i4,i5)=0.
    enddo
   enddo
  enddo
 enddo
enddo 
return
end subroutine

subroutine cerase_divergence_mat____(mat)
implicit none
complex(8) :: mat(:,:,:,:,:,:)
integer :: i1,i2,i3,i4,i5,i6,s(6)
s=shape(mat)
do i1=1,s(1)
 do i2=1,s(2)
  do i3=1,s(3)
   do i4=1,s(4)
    do i5=1,s(5)
     do i6=1,s(6)
      if(ISNAN_TEST(mat(i1,i2,i3,i4,i5,i6)).or.ISINF_TEST(mat(i1,i2,i3,i4,i5,i6))) mat(i1,i2,i3,i4,i5,i6)=0.
     enddo
    enddo
   enddo
  enddo
 enddo
enddo
return
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

 pure real(8) function dble_r(zdum)
 implicit none
  complex(8),intent(in) :: zdum
  dble_r=real(zdum)
 end function

 pure real(8) function aimag_r(zdum)
 implicit none
  complex(8),intent(in) :: zdum
  aimag_r=(0.0d0,-1.0d0)*zdum
 end function

 pure real(8) function cabs1_r(zdum)
 implicit none
 complex(8),intent(in) :: zdum
 real(16) :: tt
  tt=dabs(dble(zdum)) + dabs(aimag(zdum))
  if(abs(tt)>MAX_REAL)then
   cabs1_r=MAX_REAL
   return
  endif
  cabs1_r = tt
 end function

 pure real(8) function csign1_r(zdum1,zdum2)
 implicit none
  complex(8),intent(in) :: zdum1,zdum2
  csign1_r = cabs1(zdum1)*(zdum2/cabs1(zdum2))
 end function

 pure real(8) function DEXPc_r(rr)
 implicit none
  real(8),intent(in) :: rr 
  if(rr<MAX_EXP) then
   if(rr<MIN_EXP)then
    DEXPc_r=0.d0
   else
    DEXPc_r=EXP(rr)
   endif
  else
    DEXPc_r=EXP(MAX_EXP)
  endif
 end function

 pure real(8) function DEXPc_rr(rr)
 implicit none
  real(4),intent(in) :: rr
  if(rr<MAX_EXP_r) then
   if(rr<MIN_EXP_r)then
    DEXPc_rr=0.d0
   else
    DEXPc_rr=EXP(rr)
   endif
  else
    DEXPc_rr=EXP(MAX_EXP_r)
  endif
 end function

  !-------------------------!

 pure real(16) function dble_q(zdum)
 implicit none
  complex(16),intent(in) :: zdum
  dble_q=real(zdum)
 end function

 pure real(16) function aimag_q(zdum)
 implicit none
  complex(16),intent(in) :: zdum
  aimag_q=real(-1.0d0*zdum*imi)
 end function

 pure real(16) function cabs1_q(zdum)
 implicit none
 complex(16),intent(in) :: zdum
  cabs1_q = abs(dble_q(zdum)) + abs(aimag_q(zdum))
 end function

 pure real(16) function csign1_q(zdum1,zdum2)
 implicit none
  complex(16),intent(in) :: zdum1,zdum2
  csign1_q = cabs1_q(zdum1)*(zdum2/cabs1_q(zdum2))
 end function

 pure real(16) function DEXPc_q(rr)
 implicit none
  real(16),intent(in) :: rr
  if(rr<MAX_EXP_QUAD) then
   if(rr<MIN_EXP_QUAD)then
    DEXPc_q=0.d0
   else
    DEXPC_q=EXP(rr)
   endif
  else
    DEXPc_q=EXP(MAX_EXP_QUAD)
  endif
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

elemental real function PHASE_rel2(dd)
implicit none
real,intent(in) :: dd 
  if(abs(dd)>1.e-8) then
   PHASE_rel2=dd/abs(dd)
  else
   PHASE_rel2=1.
  endif
end function

elemental function PHASE_iel(dd)
implicit none
integer,intent(in)    :: dd
integer               :: i,s1,s2,j
integer               :: PHASE_iel
  if(abs(dd)>0) then
   PHASE_iel=dd/abs(dd)
  else
   PHASE_iel=1
  endif
end function

elemental function PHASE_rel(dd)
implicit none
real(8),intent(in)     :: dd
integer               :: i,s1,s2,j
real(8)                :: PHASE_rel
  if(abs(dd)>1.d-13) then
   PHASE_rel=dd/abs(dd)
  else
   PHASE_rel=1.d0
  endif
end function

elemental function PHASE_cel(dd)
implicit none
complex(8),intent(in) :: dd
integer               :: i,s1,s2,j
complex(8)            :: PHASE_cel
  if(abs(dd)>1.d-13) then
   PHASE_cel=dd/abs(dd)
  else
   PHASE_cel=1.d0
  endif
end function

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
!**************************************************************************
!**************************************************************************

function ANGLE_r(dd)
implicit none
real(8)     :: dd(2)
integer    :: i,s1,s2,j
real(8)     :: ANGLE_r

  ANGLE_r=0.

  if(abs(dd(2))>1.d-8) then
   if(abs(dd(1))<1.d-8)then
    if(dd(2)>0.) ANGLE_r=pi/2.
    if(dd(2)<0.) ANGLE_r=-pi/2.
    return
   endif
   ANGLE_r=ATAN(dd(2)/dd(1))
   if(dd(1)<0.) ANGLE_r=ANGLE_r+pi
  else
   if(norme(dd)>1.d-8)then
    if(dd(1)>0.) ANGLE_r =0.
    if(dd(1)<=0.)ANGLE_r= pi
   endif
  endif

 do 
  if(ANGLE_r>pi)  ANGLE_r=ANGLE_r-2.*pi
  if(ANGLE_r<-pi) ANGLE_r=ANGLE_r+2.*pi
  if(abs(ANGLE_r)<pi+error) exit
 enddo

end function

 !----------------------------!

function ANGLE_rs(dd)
implicit none
real(4)     :: dd(2)
integer    :: i,s1,s2,j
real(8)     :: ANGLE_rs

  ANGLE_rs=0.
  if(abs(dd(2))>1.d-8) then
   if(abs(dd(1))<1.d-8)then
    if(dd(2)>0.) ANGLE_rs=pi/2.
    if(dd(2)<0.) ANGLE_rs=-pi/2.
    return
   endif 
   ANGLE_rs=ATAN(dd(2)/dd(1))
   if(dd(1)<0.) ANGLE_rs=ANGLE_rs+pi
  else
   if(norme(dd)>1.d-8)then
    if(dd(1)>0.) ANGLE_rs =0.
    if(dd(1)<=0.)ANGLE_rs= pi
   endif
  endif

  do
    if(ANGLE_rs>pi)  ANGLE_rs=ANGLE_rs-2.*pi
    if(ANGLE_rs<-pi) ANGLE_rs=ANGLE_rs+2.*pi
    if(abs(ANGLE_rs)<pi+error) exit
  enddo

end function

 !----------------------------!

function ANGLE_c(dd)
implicit none
complex(8) :: dd(:)
integer    :: i,s1,s2,j
real(8)     :: ANGLE_c(size(dd(:))),vv(2)
s1=size(dd(:))
 ANGLE_c=0.
 do i=1,s1
  vv(1)=real(dd(i))
  vv(2)=aimag(dd(i))
  ANGLE_c(i)=ANGLE_r(vv)
 enddo
end function

 !----------------------------!

function ANGLE_c2(dd)
implicit none
complex(8) :: dd(:,:)
integer    :: i,s1,s2,j
real(8)     :: ANGLE_c2(size(dd(:,1)),size(dd(1,:))),vv(2)
s1=size(dd(:,1))
s2=size(dd(1,:))
 ANGLE_c2=0.d0
 do i=1,s1
  do j=1,s2
   vv(1)=real(dd(i,j))
   vv(2)=aimag(dd(i,j))
   ANGLE_c2(i,j)=ANGLE_r(vv)
  enddo
 enddo
end function

 !----------------------------!

function ANGLE_cb(dd)
implicit none
complex(8) :: dd
integer    :: i,s1,s2,j
real(8)     :: ANGLE_cb,vv(2)
 ANGLE_cb=0.
 vv(1)=real(dd)
 vv(2)=aimag(dd)
 ANGLE_cb=ANGLE_r(vv)
end function

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
!**************************************************************************
!**************************************************************************

elemental function MPLX(dd)
implicit none
complex(8)            :: MPLX
real(8),intent(in)    :: dd
  MPLX=CMPLX(cos(dd),sin(dd),kind=8)
end function

    !------------------!

elemental function MPLXamp(dd)
implicit none
complex(8)            :: MPLXamp
real(8),intent(in)    :: dd
  MPLXamp=CMPLX(1.d0,dd,kind=8) 
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function int_to_array_1(n,inti)
implicit none
integer,intent(in) :: n,inti
integer :: i
integer :: int_to_array_1(n)
int_to_array_1=0
do i=1,n
 if(ibits(inti,i-1,1)/=0) int_to_array_1(i)=1 
enddo
end function

 !---------------------------!

function int_to_array_2(n,m,inti)
implicit none
integer,intent(in) :: n,inti,m
integer :: i,j
integer :: int_to_array_2(n)
int_to_array_2=0
j=0
if(m==0)return
do i=1,n
 if(ibits(inti,i-1,1)/=0) then
  int_to_array_2(i)=1
  j=j+1
  if(j==m)return
 endif
enddo
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function array_to_int_1(n,arr)
implicit none
integer :: array_to_int_1
integer,intent(in) :: n
integer :: i
integer,intent(in) :: arr(n)
 array_to_int_1=0
 do i=1,n
  if(arr(i)/=0) then 
    array_to_int_1=ibset(array_to_int_1,i-1)
  endif
 enddo
end function

  !---------------!

function array_to_int_2(n,m,arr)
implicit none
integer :: array_to_int_2
integer,intent(in) :: n,m
integer :: i,j
integer,intent(in) :: arr(n)
 array_to_int_2=0
 if(m==0)return
 j=0
 do i=1,n
  if(arr(i)/=0) then
    array_to_int_2=ibset(array_to_int_2,i-1)
    j=j+1
    if(j==m)return
  endif
 enddo
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function ij_to_k(N,i,j)
implicit none
 !i=1....N range en lignes
 integer :: ij_to_k,N,i,j
 ij_to_k=(i-1)*N+j
end function

 !-----------------!

function k_to_ij(N,k)
integer :: k_to_ij(2),N,m,i,j
k_to_ij(2)=mod(k,N)
if(k_to_ij(2)==0)k_to_ij(2)=N
k_to_ij(1)=(k-k_to_ij(2))/N+1
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function sta(i)
implicit none
integer :: sta(1),i
  sta(1)=i
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine printmatrice__(AA,format)
implicit none
complex(8) :: AA(:,:)
integer :: i,k,sizea,sizeb
character*(*) :: format
 sizeb=size(AA(1,:))
 sizea=size(AA(:,1))
 writE(*,*) '=========== MATRICE =============='
 do i=1,sizea
  write(*,'('//format//')') (AA(i,k),k=1,sizeb)
 enddo
 write(*,*) '=================================='
end subroutine

       !-----------------------------!

subroutine printmatrice_(AA,format)
implicit none
real(8) :: AA(:,:)
integer :: i,k,sizea,sizeb
character*(*) :: format 
 sizeb=size(AA(1,:))
 sizea=size(AA(:,1))
 writE(*,*) '=========== MATRICE =============='
 do i=1,sizea
  write(*,'('//format//')') (AA(i,k),k=1,sizeb)
 enddo
 write(*,*) '=================================='
end subroutine

       !-----------------------------!

subroutine check_symmetric(errmess,A,nodiag)
implicit none
complex(8)       :: A(:,:),B(size(A,1),size(A,2))
integer          :: qq,i,j,i1,i2
character*(*)    :: errmess
real(8)          :: norm_max
logical,optional :: nodiag

 B=A; qq=size(B(:,1)); norm_max=maxval(abs(B))
 if(present(nodiag))then
  do i=1,size(A,1)
   B(i,i)=0.
  enddo
 endif

if(maxval(abs(B-TRANSPOSE(B)))/max(1.d-8,norm_max)>0.03) then
 if(messages2)then
 write(*,*) '====================================================='
 write(*,*)  errmess
 write(*,*) 'max element      : ', norm_max
 write(*,*) 'B-Trans(conj(B)) : ', maxval(abs(B-CONJG(TRANSPOSE(B))))
 write(*,*) ' element         : ', maxloc(abs(B-CONJG(TRANSPOSE(B))))
 do i=1,3 !qq
  do j=1,3 !qq
   write(*,*) i,j,abs(B(i,j)-B(j,i))
  enddo
 enddo
 if(qq<8)then
    do i1=1,qq
    WRITE(*,'(100f10.4)',err=10) (DBLE(B(i1,i2)),AIMAG(B(i1,i2)),i2=1,qq)
    enddo
 endif
 write(*,*) '====================================================='
 endif
 10 continue 
 if(strongstop) stop
endif

end subroutine

       !-----------------------------!

logical function is_symmetric(A,nodiag)
implicit none
complex(8)       :: A(:,:),B(size(A,1),size(A,2))
integer          :: qq,i,j,i1,i2
real(8)          :: norm_max
logical,optional :: nodiag
B=A; qq=size(B(:,1)); norm_max=maxval(abs(B))
if(present(nodiag))then; do i=1,size(A,1); B(i,i)=0.; enddo; endif
is_symmetric=.not.(maxval(abs(B-TRANSPOSE(B)))/max(1.d-8,norm_max)>0.03).or.maxval(abs(B-TRANSPOSE(B)))<1.d-14
end function

       !-----------------------------!

subroutine check_hermitian(errmess,A,nodiag)
implicit none
complex(8)       :: A(:,:),B(size(A,1),size(A,2))
integer          :: qq,i,j,i1,i2
character*(*)    :: errmess
real(8)          :: norm_max
logical,optional :: nodiag

 B=A; qq=size(B(:,1)); norm_max=maxval(abs(B))

 if(present(nodiag))then
  do i=1,size(A,1)
   B(i,i)=0.
  enddo
 endif

if(maxval(abs(B-CONJG(TRANSPOSE(B))))/max(1.d-8,norm_max)>0.03) then 
 if(messages2)then
 write(*,*) '====================================================='
 write(*,*)  errmess
 write(*,*) 'max element      : ', norm_max
 write(*,*) 'B-Trans(conj(B)) : ', maxval(abs(B-CONJG(TRANSPOSE(B))))
 write(*,*) ' element         : ', maxloc(abs(B-CONJG(TRANSPOSE(B))))
 do i=1,3 !qq
  do j=1,3 !qq
   write(*,*) i,j,abs(B(i,j)-CONJG(B(j,i)))
  enddo
 enddo
 if(qq<8)then
    do i1=1,qq
    WRITE(*,'(100f10.4)',err=10) (DBLE(B(i1,i2)),AIMAG(B(i1,i2)),i2=1,qq)
    enddo
 endif
 write(*,*) '====================================================='
 endif
 10 continue
 if(strongstop) stop
endif

end subroutine

       !-----------------------------!

subroutine check_squew(errmess,B)
implicit none
complex(8)    :: B(:,:)
real(8)       :: maxvalB
integer       :: qq,i,j
character*(*) :: errmess

qq=size(B(:,1)); maxvalB=maxval(abs(B))

if(maxval(abs(B+TRANSPOSE(B)))/maxvalB>0.1) then
if(messages3) then
 write(*,*) 'error matrix not squew'
 write(*,*) errmess
 write(*,*) 'B+Trans(B) : ', maxval(abs(B+TRANSPOSE(B)))
 write(*,*) 'maxval B   : ', maxvalB
 do i=1,10 !qq
  do j=1,10 !qq
   write(*,*) i,j,B(i,j)
  enddo
 enddo
 endif
 if(strongstop) stop
endif

end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine vecnorm(x,y)
implicit none
real(8),dimension(:),intent(in) ::  x
real(8),dimension(:),intent(out)::  y
real(8) :: norm
integer:: j
 norm=norme(x)
 y=x/norm
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function matrixchange(g,h,i)
implicit none
real(8)                          :: o,p
real(8),dimension(3),intent(in)  :: g,h,i
real(8),dimension(3,3)           :: matrix,matrixchange

 o=-g(3)*h(2)*i(1)+g(2)*h(3)*i(1)+g(3)*h(1)*i(2)& 
 & -g(1)*h(3)*i(2)-g(2)*h(1)*i(3)+g(1)*h(2)*i(3)

  p=g(3)*h(2)*i(1)-g(2)*h(3)*i(1)-g(3)*h(1)*i(2) & 
 & +g(1)*h(3)*i(2)+g(2)*h(1)*i(3)-g(1)*h(2)*i(3)

  if(abs(o)<epsilonr) o=epsilonr
  if(abs(p)<epsilonr) p=epsilonr

  matrix(1,1)=(-h(3)*i(2)+h(2)*i(3))/o
  matrix(1,2)=(-h(3)*i(1)+h(1)*i(3))/p
  matrix(1,3)=( h(2)*i(1)-h(1)*i(2))/p
  matrix(2,1)=(-g(3)*i(2)+g(2)*i(3))/p
  matrix(2,2)=(-g(3)*i(1)+g(1)*i(3))/o
  matrix(2,3)=(-g(2)*i(1)+g(1)*i(2))/p
  matrix(3,1)=( g(3)*h(2)-g(2)*h(3))/p
  matrix(3,2)=(-g(3)*h(1)+g(1)*h(3))/p
  matrix(3,3)=(-g(2)*h(1)+g(1)*h(2))/o

  matrixchange=matrix

end function 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function matrixrotation(f)
implicit none
real(8),intent(in)::f
real(8),dimension(3,3)::rotation,matrixrotation

    rotation(1,1)= 1.
    rotation(1,2)= 0.
    rotation(1,3)= 0.
    rotation(2,1)= 0.
    rotation(2,2)= dcos(f)
    rotation(2,3)= dsin(f)
    rotation(3,1)= 0.
    rotation(3,2)=-dsin(f)
    rotation(3,3)= dcos(f)

    matrixrotation=rotation

end function 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function matrixinverse_qr(j,q)
implicit none
real(16)                           :: q
real(16),dimension(3,3),intent(in) :: j
real(16),dimension(3,3)            :: matrixinverse_qr

    q=-j(3,1)*j(2,2)*j(1,3)+j(2,1)*j(3,2)*j(1,3)+j(3,1)*j(1,2)*j(2,3) &
    & -j(1,1)*j(3,2)*j(2,3)-j(2,1)*j(1,2)*j(3,3)+j(1,1)*j(2,2)*j(3,3)

    if(abs(q)<epsilonr) q=epsilonr

    matrixinverse_qr(1,1)=(-j(3,2)*j(2,3)+j(2,2)*j(3,3))/q
    matrixinverse_qr(2,1)=(j(3,1)*j(2,3)-j(2,1)*j(3,3))/q
    matrixinverse_qr(3,1)=(-j(3,1)*j(2,2)+j(2,1)*j(3,2))/q
    matrixinverse_qr(1,2)=(j(3,2)*j(1,3)-j(1,2)*j(3,3))/q
    matrixinverse_qr(2,2)=(-j(3,1)*j(1,3)+j(1,1)*j(3,3))/q
    matrixinverse_qr(3,2)=(j(3,1)*j(1,2)-j(1,1)*j(3,2))/q
    matrixinverse_qr(1,3)=(-j(2,2)*j(1,3)+j(1,2)*j(2,3))/q
    matrixinverse_qr(2,3)=(j(2,1)*j(1,3)-j(1,1)*j(2,3))/q
    matrixinverse_qr(3,3)=(-j(2,1)*j(1,2)+j(1,1)*j(2,2))/q

return
end function

   !-----------------------------------------------!

function matrixinverse_qc(j,q)
implicit none
complex(16)                           :: q
complex(16),dimension(3,3),intent(in) :: j
complex(16),dimension(3,3)            :: matrixinverse_qc

    q=-j(3,1)*j(2,2)*j(1,3)+j(2,1)*j(3,2)*j(1,3)+j(3,1)*j(1,2)*j(2,3) &
    & -j(1,1)*j(3,2)*j(2,3)-j(2,1)*j(1,2)*j(3,3)+j(1,1)*j(2,2)*j(3,3)

    if(abs(q)<epsilonr) q=epsilonr

    matrixinverse_qc(1,1)=(-j(3,2)*j(2,3)+j(2,2)*j(3,3))/q
    matrixinverse_qc(2,1)=(j(3,1)*j(2,3)-j(2,1)*j(3,3))/q
    matrixinverse_qc(3,1)=(-j(3,1)*j(2,2)+j(2,1)*j(3,2))/q
    matrixinverse_qc(1,2)=(j(3,2)*j(1,3)-j(1,2)*j(3,3))/q
    matrixinverse_qc(2,2)=(-j(3,1)*j(1,3)+j(1,1)*j(3,3))/q
    matrixinverse_qc(3,2)=(j(3,1)*j(1,2)-j(1,1)*j(3,2))/q
    matrixinverse_qc(1,3)=(-j(2,2)*j(1,3)+j(1,2)*j(2,3))/q
    matrixinverse_qc(2,3)=(j(2,1)*j(1,3)-j(1,1)*j(2,3))/q
    matrixinverse_qc(3,3)=(-j(2,1)*j(1,2)+j(1,1)*j(2,2))/q

return
end function

   !-----------------------------------------------!

function matrixinverse_r(j,q)
implicit none
real(8)                           :: q
real(8),dimension(3,3),intent(in) :: j
real(8),dimension(3,3)            :: matrixinverse_r

    q=-j(3,1)*j(2,2)*j(1,3)+j(2,1)*j(3,2)*j(1,3)+j(3,1)*j(1,2)*j(2,3) &
    & -j(1,1)*j(3,2)*j(2,3)-j(2,1)*j(1,2)*j(3,3)+j(1,1)*j(2,2)*j(3,3)

    if(abs(q)<epsilonr) q=epsilonr

    matrixinverse_r(1,1)=(-j(3,2)*j(2,3)+j(2,2)*j(3,3))/q
    matrixinverse_r(2,1)=(j(3,1)*j(2,3)-j(2,1)*j(3,3))/q
    matrixinverse_r(3,1)=(-j(3,1)*j(2,2)+j(2,1)*j(3,2))/q
    matrixinverse_r(1,2)=(j(3,2)*j(1,3)-j(1,2)*j(3,3))/q
    matrixinverse_r(2,2)=(-j(3,1)*j(1,3)+j(1,1)*j(3,3))/q
    matrixinverse_r(3,2)=(j(3,1)*j(1,2)-j(1,1)*j(3,2))/q
    matrixinverse_r(1,3)=(-j(2,2)*j(1,3)+j(1,2)*j(2,3))/q
    matrixinverse_r(2,3)=(j(2,1)*j(1,3)-j(1,1)*j(2,3))/q
    matrixinverse_r(3,3)=(-j(2,1)*j(1,2)+j(1,1)*j(2,2))/q

return
end function

   !-----------------------------------------------!

function matrixinverse_rs(j,q)
implicit none
real(4)                           :: q
real(4),dimension(3,3),intent(in) :: j
real(4),dimension(3,3)            :: matrixinverse_rs

    q=-j(3,1)*j(2,2)*j(1,3)+j(2,1)*j(3,2)*j(1,3)+j(3,1)*j(1,2)*j(2,3) &
    & -j(1,1)*j(3,2)*j(2,3)-j(2,1)*j(1,2)*j(3,3)+j(1,1)*j(2,2)*j(3,3)

    if(abs(q)<1.d-8) q=1.d-8

    matrixinverse_rs(1,1)=(-j(3,2)*j(2,3)+j(2,2)*j(3,3))/q
    matrixinverse_rs(2,1)=(j(3,1)*j(2,3)-j(2,1)*j(3,3))/q
    matrixinverse_rs(3,1)=(-j(3,1)*j(2,2)+j(2,1)*j(3,2))/q
    matrixinverse_rs(1,2)=(j(3,2)*j(1,3)-j(1,2)*j(3,3))/q
    matrixinverse_rs(2,2)=(-j(3,1)*j(1,3)+j(1,1)*j(3,3))/q
    matrixinverse_rs(3,2)=(j(3,1)*j(1,2)-j(1,1)*j(3,2))/q
    matrixinverse_rs(1,3)=(-j(2,2)*j(1,3)+j(1,2)*j(2,3))/q
    matrixinverse_rs(2,3)=(j(2,1)*j(1,3)-j(1,1)*j(2,3))/q
    matrixinverse_rs(3,3)=(-j(2,1)*j(1,2)+j(1,1)*j(2,2))/q

return
end function

   !-----------------------------------------------!

function matrixinverse_c(j,q)
implicit none
complex(8)                           :: q
complex(8),dimension(3,3),intent(in) :: j
complex(8),dimension(3,3)            :: matrixinverse_c

    q=-j(3,1)*j(2,2)*j(1,3)+j(2,1)*j(3,2)*j(1,3)+j(3,1)*j(1,2)*j(2,3) &
    & -j(1,1)*j(3,2)*j(2,3)-j(2,1)*j(1,2)*j(3,3)+j(1,1)*j(2,2)*j(3,3)

    if(abs(q)<epsilonr) q=epsilonr

    matrixinverse_c(1,1)=(-j(3,2)*j(2,3)+j(2,2)*j(3,3))/q
    matrixinverse_c(2,1)=(j(3,1)*j(2,3)-j(2,1)*j(3,3))/q
    matrixinverse_c(3,1)=(-j(3,1)*j(2,2)+j(2,1)*j(3,2))/q
    matrixinverse_c(1,2)=(j(3,2)*j(1,3)-j(1,2)*j(3,3))/q
    matrixinverse_c(2,2)=(-j(3,1)*j(1,3)+j(1,1)*j(3,3))/q
    matrixinverse_c(3,2)=(j(3,1)*j(1,2)-j(1,1)*j(3,2))/q
    matrixinverse_c(1,3)=(-j(2,2)*j(1,3)+j(1,2)*j(2,3))/q
    matrixinverse_c(2,3)=(j(2,1)*j(1,3)-j(1,1)*j(2,3))/q
    matrixinverse_c(3,3)=(-j(2,1)*j(1,2)+j(1,1)*j(2,2))/q

return
end function

   !-----------------------------------------------!

function matrixinverse_cs(j,q)
implicit none
complex(4)                           :: q
complex(4),dimension(3,3),intent(in) :: j
complex(4),dimension(3,3)            :: matrixinverse_cs

    q=-j(3,1)*j(2,2)*j(1,3)+j(2,1)*j(3,2)*j(1,3)+j(3,1)*j(1,2)*j(2,3) &
    & -j(1,1)*j(3,2)*j(2,3)-j(2,1)*j(1,2)*j(3,3)+j(1,1)*j(2,2)*j(3,3)

    if(abs(q)<epsilonr) q=epsilonr

    matrixinverse_cs(1,1)=(-j(3,2)*j(2,3)+j(2,2)*j(3,3))/q
    matrixinverse_cs(2,1)=(j(3,1)*j(2,3)-j(2,1)*j(3,3))/q
    matrixinverse_cs(3,1)=(-j(3,1)*j(2,2)+j(2,1)*j(3,2))/q
    matrixinverse_cs(1,2)=(j(3,2)*j(1,3)-j(1,2)*j(3,3))/q
    matrixinverse_cs(2,2)=(-j(3,1)*j(1,3)+j(1,1)*j(3,3))/q
    matrixinverse_cs(3,2)=(j(3,1)*j(1,2)-j(1,1)*j(3,2))/q
    matrixinverse_cs(1,3)=(-j(2,2)*j(1,3)+j(1,2)*j(2,3))/q
    matrixinverse_cs(2,3)=(j(2,1)*j(1,3)-j(1,1)*j(2,3))/q
    matrixinverse_cs(3,3)=(-j(2,1)*j(1,2)+j(1,1)*j(2,2))/q

return
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function determinant(mat)
implicit none
real(8),intent(in)::mat(:,:)
real(8),dimension(3)::detc
real(8)::determinant
 if(size(mat(1,:))==3)then
   detc(1)=mat(1,1)*mat(2,2)*mat(3,3)-mat(1,1)*mat(3,2)*mat(2,3)
   detc(2)=mat(2,1)*mat(3,2)*mat(1,3)-mat(2,1)*mat(1,2)*mat(3,3)
   detc(3)=mat(3,1)*mat(1,2)*mat(2,3)-mat(3,1)*mat(2,2)*mat(1,3)
   determinant=sum(detc)
  endif
  if(size(mat(1,:))==2) determinant=mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2)
end function

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
!**************************************************************************
!**************************************************************************

function vecegali(v3,v4)
 implicit none
 logical :: vecegali
 integer,intent(in),dimension(:)  :: v3,v4
 integer ::i
 vecegali=.true.
 if(size(v3)/=size(v4))then
  vecegali=.false.
  return
 endif
 do i=1,size(v3)
  if(v3(i)/=v4(i))then
   vecegali=.false.
   return
  endif
 enddo
end function

!**************************************************************************

function vecegalr(vv3,vv4)
 implicit none
 logical :: vecegalr
 real(8),intent(in),dimension(:)  ::vv3,vv4
 integer ::i
 vecegalr=.true.
 if(size(vv3)/=size(vv4))then
  vecegalr=.false.
  return
 endif
 if(norme(vv3-vv4)>prec) vecegalr=.false.
end function

!**************************************************************************

function vecbigr(v3,v4)
 implicit none
 logical :: vecbigr
 real(8),intent(in),dimension(:) :: v3,v4
 vecbigr=.true.
 if(size(v3)/=size(v4))then
  vecbigr=.false.
  return
 endif
  if(norme(v3)<norme(v4))then
   vecbigr=.false.
   return
  endif
end function

!**************************************************************************

 function vecbigi(v3,v4)
 implicit none
 logical :: vecbigi
 integer,intent(in),dimension(:) :: v3,v4
 vecbigi=.true.
 if(size(v3)/=size(v4))then
  vecbigi=.false.
  return
 endif
  if(norme(v3)<norme(v4))then
   vecbigi=.false.
   return
  endif
end function

!**************************************************************************

function vecsmr(v3,v4)
 implicit none
 logical :: vecsmr
 real(8),intent(in),dimension(:) :: v3,v4
 vecsmr=.not.(vecbigr(v3,v4))
end function

!**************************************************************************

function vecsmi(vv3,vv4)
implicit none
logical :: vecsmi
integer,intent(in),dimension(:) :: vv3,vv4
 vecsmi=.not.(vecbigi(vv3,vv4))
end function

!**************************************************************************

function vecegalis(v3,v4)
implicit none
logical :: vecegalis
integer,intent(in),dimension(:)  :: v3
integer,intent(in) :: v4
integer::i
 vecegalis=.true.
 do i=1,size(v3)
  if(v3(i)/=v4)then
   vecegalis=.false.
   return
  endif
 enddo
end function

!**************************************************************************

function vecegalrs(vv3,vv4)
implicit none
logical :: vecegalrs
real(8),intent(in),dimension(:)  :: vv3
real(8),intent(in) :: vv4
integer::i
 vecegalrs=.true.
 if(norme(vv3-vv4)>prec) vecegalrs=.false.
end function

!**************************************************************************

function vecbigrs(v3,v4)
implicit none
logical :: vecbigrs
real(8),intent(in),dimension(:) :: v3
real(8),intent(in) :: v4
 vecbigrs=.true.
  if(norme(v3)<v4)then
   vecbigrs=.false.
   return
  endif
end function

!**************************************************************************

function vecbigis(v3,v4)
implicit none
logical :: vecbigis
integer,intent(in),dimension(:) :: v3
integer,intent(in) :: v4
 vecbigis=.true.
  if(norme(v3)<v4)then
   vecbigis=.false.
   return
  endif
end function

!**************************************************************************

function vecsmrs(vv3,vv4)
implicit none
logical vecsmrs
real(8),intent(in),dimension(:) :: vv3
real(8),intent(in) :: vv4
 vecsmrs=.not.(vecbigrs(vv3,vv4))
end function

!**************************************************************************

function vecsmis(vv3,vv4) 
implicit none
logical :: vecsmis
integer,intent(in),dimension(:) :: vv3
integer,intent(in) :: vv4
 vecsmis=.not.(vecbigis(vv3,vv4))
end function

!**************************************************************************

function cdot_(x,y)
implicit none
integer                            :: j
complex(8),dimension(:),intent(in) :: x,y
complex(8)                         :: cdot_
 cdot_=0.
 do j=1,size(x)
  cdot_ = cdot_ + conjg(x(j))*y(j)
 enddo
end function

!**********************************************************************************

function ddot_(xx,yy)
implicit none
integer                         ::  j
real(8),dimension(:),intent(in) ::  xx,yy
real(8)                         ::  ddot_
 ddot_=0.
 do j=1,size(xx)
  ddot_=ddot_+xx(j)*yy(j)
 enddo
end function

!**********************************************************************************

function sdot_(x,y)
implicit none
integer ::j
real(4),dimension(:),intent(in) ::  x,y
real(8) :: sdot_
 sdot_=0.
 do j=1,size(x)
  sdot_=sdot_+x(j)*y(j)
 enddo
end function

!**********************************************************************************

function idot_(u,v)
implicit none
integer,dimension(:),intent(in) ::u,v
real(8) :: idot_
integer :: j
 idot_=0.
 do j=1,size(u)
  idot_=idot_+dble(u(j)*v(j))
 enddo
end function

!*********************************************************************************

function icross(x1,y1)
implicit none
integer, dimension(3),intent(in) ::  x1,y1
real(8),dimension(3)::  icross
  icross(1) = x1(2)*y1(3) -x1(3)*y1(2)
  icross(2) = x1(3)*y1(1) -x1(1)*y1(3)
  icross(3) = x1(1)*y1(2) -x1(2)*y1(1)
end function

!**********************************************************************************

function dnorm(xx)
implicit none
real(8), dimension(:),intent(in) :: xx 
real(8) :: dnorm
integer :: j
 dnorm=0.
if(enable_mpi_dot)then
!$OMP PARALLEL PRIVATE(j), SHARED(xx), REDUCTION(+:dnorm)
!$OMP DO
 do j=1,size(xx)
  dnorm=dnorm+xx(j)**2
 enddo
!$OMP END DO
!$OMP END PARALLEL
else
 do j=1,size(xx)
  dnorm=dnorm+xx(j)**2
 enddo
endif
 dnorm=SQRT(dnorm)
end function

!**********************************************************************************

function norme_squ(xx)
implicit none
real(8), dimension(:),intent(in) :: xx
real(8) :: norme_squ
integer :: j
 norme_squ=0.
if(enable_mpi_dot)then
!$OMP PARALLEL PRIVATE(j), SHARED(xx), REDUCTION(+:norme_squ)
!$OMP DO
 do j=1,size(xx)
  norme_squ=norme_squ+xx(j)**2
 enddo
!$OMP END DO
!$OMP END PARALLEL
else
 do j=1,size(xx)
  norme_squ=norme_squ+xx(j)**2
 enddo
endif
end function

!**********************************************************************************

function inorm(xx)
implicit none
integer,dimension(:),intent(in)::xx
real(8) ::  inorm
integer :: j
 inorm=0.
 do j=1,size(xx)
  inorm=inorm+dble(xx(j))**2
 enddo
 inorm=SQRT(inorm)
end function

!**********************************************************************************

function cnorm(xx)
implicit none
complex(8),dimension(:),intent(in)::xx
real(8) ::  cnorm
integer ::  j
 cnorm=0.
 do j=1,size(xx)
  cnorm=cnorm+abs(xx(j))**2
 enddo
 cnorm=SQRT(cnorm)
end function

!**********************************************************************************

function snorm(x)
implicit none
real(4), dimension(:),intent(in) ::  x
real(8):: snorm
integer :: j
 snorm=0.
 do j=1,size(x)
  snorm=snorm+ x(j)**2
 enddo
 snorm=SQRT(snorm)
end function

!**********************************************************************************

function dcross(x,y)
implicit none
real(8), dimension(3),intent(in) ::  x,y
real(8), dimension(3)            ::  dcross
 dcross(1) = x(2)*y(3) -x(3)*y(2)
 dcross(2) = x(3)*y(1) -x(1)*y(3)
 dcross(3) = x(1)*y(2) -x(2)*y(1)
end function

!**********************************************************************************

function scross(xx,yy)
implicit none
real(4), dimension(3),intent(in) ::  xx,yy
real(8), dimension(3)            ::  scross
 scross(1) = xx(2)*yy(3) -xx(3)*yy(2)
 scross(2) = xx(3)*yy(1) -xx(1)*yy(3)
 scross(3) = xx(1)*yy(2) -xx(2)*yy(1)
end function

!**********************************************************************************

function ccross(x,y)
implicit none
complex(8), dimension(3),intent(in) ::  x,y
complex(8), dimension(3)            :: ccross
 ccross(1) = ( x(2)*conjg(y(3)) -x(3)*conjg(y(2)) )
 ccross(2) = ( x(3)*conjg(y(1)) -x(1)*conjg(y(3)) )
 ccross(3) = ( x(1)*conjg(y(2)) -x(2)*conjg(y(1)) )
return
end function

!**********************************************************************************




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
!**************************************************************************


SUBROUTINE eigsys(ham, olm, zek, evl, evr, ndim)

!!!----------------------------------------------------------------------!!!
!!! This routine solves the generalized right/left eigenvalue problems   !!!
!!!       H.Ev^{R} = E O.Ev^{R}    and    Ev^{L}H = E Ev^{L}.O           !!!
!!! ( . denotes matrix multiplication )  where H is a general matrix and !!!
!!! O is a Hermitian positive definite matrix.                           !!!
!!! The problem is solved in several steps:                              !!!
!!!                                                                      !!!
!!! 1) Cholesky decompose O    :  O = U^{+}.U                            !!!
!!! 2) Invert Cholesky factor  :  B = U^{-1}                             !!!
!!! 3) Transform Hamiltonian   :  H' = B^{+}.H.B                         !!!
!!! 4) Solve eigenvalue problem:  H'.A^{R} = EA^{R} & A^{L}.H' = EA^{L}  !!!
!!! 5) Normalize eigenvectors  :  A^{L}.A^{R} = 1                        !!!
!!! 6) Transform eigenvectors  :  Ev^{R} = B.A^{R} & Ev^{L} = A^{L}.B^{+}!!!
!!! 7) Reorder eigenvalues                                               !!!
!!!                                                                      !!!
!!! Variable contents on entry:                                          !!!
!!! ===========================                                          !!!
!!!                                                                      !!!
!!! ham        : Hamiltonian                                             !!!
!!! olm        : Overlap matrix                                          !!!
!!! ndim       : Dimension of matrices                                   !!!
!!! fh_stderr  : Standard error file handle (Unit number)                !!!
!!!                                                                      !!!
!!! Variable contents on exit:                                           !!!
!!! ==========================                                           !!!
!!!                                                                      !!!
!!! ham        : Right eigenvectors as columns                           !!!
!!! olm        : Left eigenvectors as rows                               !!!
!!! zek        : Eigenvalues                                             !!!
!!!----------------------------------------------------------------------!!!

  IMPLICIT NONE

  complex(8), INTENT(in)   :: ham(ndim,ndim)   ! Hamiltonian
  complex(8), INTENT(in)   :: olm(ndim,ndim)   ! Overlap matrix
  complex(8), INTENT(out)  :: zek(ndim)        ! Vector of eigenvalues 
  complex(8), INTENT(out)  :: evl(ndim,ndim)   ! left eigenvectors
  complex(8), INTENT(out)  :: evr(ndim,ndim)   ! right eigenvectors
  INTEGER, INTENT(in)      :: ndim             ! Dimension of hamiltonian
  real(8), PARAMETER       :: smalleps = 1e-5
  CHARACTER(1), PARAMETER  :: uplo  = "U"      ! Job parameter for zpotrf/ztrtri
  CHARACTER(1), PARAMETER  :: diag  = "N"      ! Job parameter for ztrtri
  CHARACTER(1), PARAMETER  :: jobvl = "V"      ! Jop parameter for zgeev
  CHARACTER(1), PARAMETER  :: jobvr = "V"      ! Job parameter for zgeev
  CHARACTER(1), PARAMETER  :: transn  = "N"    ! Job parameter for ztrmm
  CHARACTER(1), PARAMETER  :: transc  = "C"    ! Job parameter for ztrmm
  CHARACTER(1), PARAMETER  :: left   = "L"     ! Job parameter for ztrmm
  CHARACTER(1), PARAMETER  :: right  = "R"     ! Job parameter for ztrmm
  INTEGER                  :: ierr             ! Error parameter for lapack
  INTEGER                  :: irow             ! Loop index for rows
  INTEGER                  :: icol             ! Loop index for columns
  INTEGER                  :: p,q,r            ! Loop index 
  complex(8)               :: ctmp             ! Temporary variable 
  complex(8)               :: alpha_one        ! For ztrmm
  complex(8)               :: umat(ndim,ndim)  ! Cholesky factor
  complex(8)               :: bmat(ndim,ndim)  ! Inverse Cholesky factor
  complex(8)               :: htmp(ndim,ndim)  ! Transformed Hamiltonian
  complex(8)               :: scaler(ndim)     ! Array of normalization parameters
  complex(8)               :: cworkvec(8*ndim) ! Work array for zgeev
  real(8)                  :: rworkvec(8*ndim) ! Work array for zgeev

  !========== Step 1, Cholesky decompose O, O = U^{+}U ==========
  umat = olm
  CALL zpotrf(uplo,ndim,umat,ndim,ierr)
  IF(ierr.NE.0) THEN 
     WRITE(6,*) 'Error code of zpotrf ', ierr
     WRITE(0,*) 'Error in dia_gho! Stopp the code!'
  ENDIF
  !---------- umat is now upper triangular, lower
  !---------- triangular part is not referenced

  !========== Step 2, Invert Cholesky matrix U, B = U^{-1} ==========
  bmat = umat
  CALL ztrtri(uplo,diag,ndim,bmat,ndim,ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'Error code of ztrtri ', ierr
     WRITE(0,*) 'Error in dia_gho! Stopp the code!'
  ENDIF
  !---------- bmat is now upper triangular, lower
  !---------- triangular part is not referenced

  !========== Step 3, Compute new Hamiltonian, H' = B^{+}HB ==========
  alpha_one = 1.0
  htmp = ham
  CALL ztrmm(right, uplo, transn, diag, ndim, ndim, alpha_one, bmat, ndim, htmp, ndim)
  CALL ztrmm(left,  uplo, transc, diag, ndim, ndim, alpha_one, bmat, ndim, htmp, ndim)

  !========== Step 4, Solve eigenvalue problem, H'v' = Ev' ==========
  CALL zgeev(jobvl,jobvr,ndim,htmp,ndim,zek,evl,ndim,evr,ndim,cworkvec,8*ndim,rworkvec,ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'Error code of zgeev ', ierr
     WRITE(0,*) 'Error in dia_gho! Stopp the code!'
  ENDIF
  ! transpose left eigenvectors
  evl = conjg(TRANSPOSE(evl))

  !========== Step 5, Make degenerate eigenvectors orthogonal
  DO q=1,ndim
     DO p=1,q-1
        IF (abs(zek(p)-zek(q)).LT.smalleps .AND. abs(scalprod(evl(p,:),evr(:,q),ndim)) .GT.smalleps) THEN
           evr(:,q) = evr(:,q) - scalprod(evl(p,:),evr(:,q),ndim)/scalprod(evl(p,:),evr(:,p),ndim) * evr(:,p)
        ENDIF
     ENDDO
     DO p=1,q-1
        IF (abs(zek(p)-zek(q)).LT.smalleps .AND. abs(scalprod(evl(q,:),evr(:,p),ndim)) .GT.smalleps) THEN
           evl(q,:) = evl(q,:) - scalprod(evl(q,:),evr(:,p),ndim)/scalprod(evl(p,:),evr(:,p),ndim) * evl(p,:)
        ENDIF
     ENDDO
  ENDDO

  !========= Step 6, Normalize eigenvectors
  DO p = 1,ndim
     ctmp = 0.d0
     DO q = 1,ndim
        ctmp = ctmp+evl(p,q)*evr(q,p)
     ENDDO
     scaler(p) = SQRT(ctmp)
  ENDDO
  DO p = 1,ndim
     evl(p,:) = evl(p,:)/scaler(p)
     evr(:,p) = evr(:,p)/scaler(p)
  ENDDO

  !========== Step 7, Transform eigenvectors and store ==========
  CALL ztrmm(left,  uplo, transn, diag, ndim, ndim, alpha_one, bmat, ndim, evr, ndim)
  CALL ztrmm(right, uplo, transc, diag, ndim, ndim, alpha_one, bmat, ndim, evl, ndim)

  RETURN

CONTAINS
  complex(8) FUNCTION scalprod(a,b,ndim)
    IMPLICIT NONE
    complex(8) :: a(:), b(:)
    INTEGER    :: ndim
    INTEGER    :: i
    scalprod = 0.0
    DO i=1,ndim
       scalprod = scalprod + a(i)*b(i)
    ENDDO
  END FUNCTION scalprod

END SUBROUTINE eigsys


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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

SUBROUTINE eigvals(ham, olm, zek, ndim)

!!!----------------------------------------------------------------------!!!
!!! This routine solves the generalized right/left eigenvalue problems   !!!
!!!       H.Ev^{R} = E O.Ev^{R}    and    Ev^{L}H = E Ev^{L}.O           !!!
!!! ( . denotes matrix multiplication )  where H is a general matrix and !!!
!!! O is a Hermitian positive definite matrix.                           !!!
!!! The problem is solved in several steps:                              !!!
!!!                                                                      !!!
!!! 1) Cholesky decompose O    :  O = U^{+}.U                            !!!
!!! 2) Invert Cholesky factor  :  B = U^{-1}                             !!!
!!! 3) Transform Hamiltonian   :  H' = B^{+}.H.B                         !!!
!!! 4) Solve eigenvalue problem:  H'.A^{R} = EA^{R} & A^{L}.H' = EA^{L}  !!!
!!! 5) Normalize eigenvectors  :  A^{L}.A^{R} = 1                        !!!
!!! 6) Transform eigenvectors  :  Ev^{R} = B.A^{R} & Ev^{L} = A^{L}.B^{+}!!!
!!! 7) Reorder eigenvalues                                               !!!
!!!                                                                      !!!
!!! Variable contents on entry:                                          !!!
!!! ===========================                                          !!!
!!!                                                                      !!!
!!! ham        : Hamiltonian                                             !!!
!!! olm        : Overlap matrix                                          !!!
!!! ndim       : Dimension of matrices                                   !!!
!!! fh_stderr  : Standard error file handle (Unit number)                !!!
!!!                                                                      !!!
!!! Variable contents on exit:                                           !!!
!!! ==========================                                           !!!
!!!                                                                      !!!
!!! ham        : Right eigenvectors as columns                           !!!
!!! olm        : Left eigenvectors as rows                               !!!
!!! zek        : Eigenvalues                                             !!!
!!!----------------------------------------------------------------------!!!

  IMPLICIT NONE

  complex(8), INTENT(in)  :: ham(ndim,ndim)   ! Hamiltonian
  complex(8), INTENT(in)  :: olm(ndim,ndim)   ! Overlap matrix
  complex(8), INTENT(out) :: zek(ndim)        ! Vector of eigenvalues 
  INTEGER, INTENT(in)     :: ndim             ! Dimension of hamiltonian
  real(8), PARAMETER      :: smalleps = 1e-5
  CHARACTER(1),PARAMETER  :: uplo  = "U"      ! Job parameter for zpotrf/ztrtri
  CHARACTER(1),PARAMETER  :: diag  = "N"      ! Job parameter for ztrtri
  CHARACTER(1),PARAMETER  :: jobvl = "N"      ! Jop parameter for zgeev
  CHARACTER(1),PARAMETER  :: jobvr = "N"      ! Job parameter for zgeev
  CHARACTER(1),PARAMETER  :: transn  = "N"    ! Job parameter for ztrmm
  CHARACTER(1),PARAMETER  :: transc  = "C"    ! Job parameter for ztrmm
  CHARACTER(1),PARAMETER  :: left   = "L"     ! Job parameter for ztrmm
  CHARACTER(1),PARAMETER  :: right  = "R"     ! Job parameter for ztrmm
  INTEGER                 :: ierr             ! Error parameter for lapack
  INTEGER                 :: irow             ! Loop index for rows
  INTEGER                 :: icol             ! Loop index for columns
  INTEGER                 :: p,q,r            ! Loop index 
  complex(8)              :: ctmp             ! Temporary variable 
  complex(8)              :: alpha_one        ! For ztrmm
  complex(8)              :: umat(ndim,ndim)  ! Cholesky factor
  complex(8)              :: bmat(ndim,ndim)  ! Inverse Cholesky factor
  complex(8)              :: htmp(ndim,ndim)  ! Transformed Hamiltonian
  complex(8)              :: scaler(ndim)     ! Array of normalization parameters
  complex(8)              :: cworkvec(8*ndim) ! Work array for zgeev
  real(8)                 :: rworkvec(8*ndim) ! Work array for zgeev
  complex(8)              :: evl(1,1)         ! left eigenvectors
  complex(8)              :: evr(1,1)         ! right eigenvectors

  !========== Step 1, Cholesky decompose O, O = U^{+}U ==========
  umat = olm
  CALL zpotrf(uplo,ndim,umat,ndim,ierr)
  IF(ierr.NE.0) THEN 
     WRITE(6,*) 'Error code of zpotrf ', ierr
     WRITE(0,*) 'Error in dia_gho! Stopp the code!'
  ENDIF
  !---------- umat is now upper triangular, lower
  !---------- triangular part is not referenced

  !========== Step 2, Invert Cholesky matrix U, B = U^{-1} ==========
  bmat = umat
  CALL ztrtri(uplo,diag,ndim,bmat,ndim,ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'Error code of ztrtri ', ierr
     WRITE(0,*) 'Error in dia_gho! Stopp the code!'
  ENDIF
  !---------- bmat is now upper triangular, lower
  !---------- triangular part is not referenced

  !========== Step 3, Compute new Hamiltonian, H' = B^{+}HB ==========
  alpha_one = 1.0
  htmp = ham
  CALL ztrmm(right, uplo, transn, diag, ndim, ndim, alpha_one, bmat, ndim, htmp, ndim)
  CALL ztrmm(left,  uplo, transc, diag, ndim, ndim, alpha_one, bmat, ndim, htmp, ndim)

  !========== Step 4, Solve eigenvalue problem, H'v' = Ev' ==========
  CALL zgeev(jobvl,jobvr,ndim,htmp,ndim,zek,evl,ndim,evr,ndim,cworkvec,8*ndim,rworkvec,ierr)
  IF(ierr.NE.0) THEN
     WRITE(6,*) 'Error code of zgeev ', ierr
     WRITE(0,*) 'Error in dia_gho! Stopp the code!'
  ENDIF
  RETURN

END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

      SUBROUTINE CINV (N, A, IDIM, R, IFAIL) 

!     ******************************************************************
!                 REPLACES A BY ITS INVERSE.                                        
!     ******************************************************************

      IMPLICIT real(8) (A-H,O-Z)
      complex(8)  A(IDIM,N)
      INTEGER     N, IDIM
      INTEGER     IFAIL
      integer     R(N)
      real(8)     T1,T2,T3
      complex(8)  ONE,DET,TEMP,S,C11,C12,C13,C21,C22,C23,C31,C32,C33
      CHARACTER*6 NAME  

      DATA NAME/'CINV'/,KPRNT/0/    
      DATA ONE/(1.D0,0.D0)/
                                                                       
!  TEST FOR PARAMETER ERRORS.                                           
                                                                       
      IF ( (N.LT.1) .OR. (N.GT.IDIM) ) GOTO 7 
                                                                       
!  TEST FOR N.LE.3.                                                     
                                                                       
      IF (N.GT.3) GOTO 6 
      IFAIL = 0 
      IF (N.LT.3) GOTO 4 
                                                                       
!  N=3 CASE.                                                            
                                                                       
!     COMPUTE COFACTORS.                                                
      C11 = A (2, 2) * A (3, 3) - A (2, 3) * A (3, 2) 
      C12 = A (2, 3) * A (3, 1) - A (2, 1) * A (3, 3) 
      C13 = A (2, 1) * A (3, 2) - A (2, 2) * A (3, 1) 
      C21 = A (3, 2) * A (1, 3) - A (3, 3) * A (1, 2) 
      C22 = A (3, 3) * A (1, 1) - A (3, 1) * A (1, 3) 
      C23 = A (3, 1) * A (1, 2) - A (3, 2) * A (1, 1) 
      C31 = A (1, 2) * A (2, 3) - A (1, 3) * A (2, 2) 
      C32 = A (1, 3) * A (2, 1) - A (1, 1) * A (2, 3) 
      C33 = A (1, 1) * A (2, 2) - A (1, 2) * A (2, 1) 
      T1 = ABS (REAL (A (1, 1) ) ) + ABS (aimag (A (1, 1) ) ) 
      T2 = ABS (REAL (A (2, 1) ) ) + ABS (aimag (A (2, 1) ) ) 
      T3 = ABS (REAL (A (3, 1) ) ) + ABS (aimag (A (3, 1) ) ) 
                                                                       
      IF (T1.GE.T2) GOTO 1 
      IF (T3.GE.T2) GOTO 2 
      TEMP = A (2, 1) 
      DET = C13 * C32 - C12 * C33 
      GOTO 3 
    1 IF (T3.GE.T1) GOTO 2 
      TEMP = A (1, 1) 
      DET = C22 * C33 - C23 * C32 
      GOTO 3 
    2 TEMP = A (3, 1) 
      DET = C23 * C12 - C22 * C13 
!     SET ELEMENTS OF INVERSE IN A.                                     
    3 IF (REAL (DET) .EQ.0.D0.AND.aimag (DET) .EQ.0.D0) GOTO 8 
      S = TEMP / DET 
      A (1, 1) = S * C11 
      A (1, 2) = S * C21 
      A (1, 3) = S * C31 
      A (2, 1) = S * C12 
      A (2, 2) = S * C22 
      A (2, 3) = S * C32 
      A (3, 1) = S * C13 
      A (3, 2) = S * C23 
      A (3, 3) = S * C33 
      RETURN 
                                                                       
    4 IF (N.LT.2) GOTO 5 
                                                                       
!  N=2 CASE BY CRAMERS RULE.                                            
                                                                       
      DET = A (1, 1) * A (2, 2) - A (1, 2) * A (2, 1) 
      IF (REAL (DET) .EQ.0.D0.AND.aimag (DET) .EQ.0.D0) GOTO 8 
      S = ONE / DET 
      C11 = S * A (2, 2) 
      A (1, 2) = - S * A (1, 2) 
      A (2, 1) = - S * A (2, 1) 
      A (2, 2) = S * A (1, 1) 
      A (1, 1) = C11 
      RETURN 
                                                                       
!  N=1 CASE.                                                            
                                                                       
    5 IF (REAL (A (1, 1) ) .EQ.0.D0.AND.aimag (A (1, 1) ) .EQ.0.D0)    &
      GOTO 8                                                            
      A (1, 1) = ONE / A (1, 1) 
      RETURN 
                                                                       
!  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.                          
                                                                       
    6 CALL CFACT (N, A, IDIM, R, IFAIL, DET, JFAIL) 
      IF (IFAIL.NE.0) RETURN 
      CALL CFINV (N, A, IDIM, R) 
      RETURN 
                                                                       
!  ERROR EXITS.                                                         
                                                                       
    7 IFAIL = + 1 
      RETURN 
                                                                       
    8 IFAIL = - 1 
      RETURN 
                                                                       
      END SUBROUTINE
                                                                        
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
                                                                        
      SUBROUTINE CFINV (N, A, IDIM, IR) 
      IMPLICIT REAL (8)(A - H, O - Z) 
      INTEGER       :: IR (2 * N) 
      complex(8)    :: A (IDIM, N), X, Y, TI 
      CHARACTER(6)  :: HNAME 
      complex(8)    ::  ZERO, S31, S32, S33, S34, DC, DOTF 

      DOTF (X, Y, S31) = X * Y + S31 
      DATA ZERO / (0.D0, 0.D0) / 
      DATA HNAME / ' CFINV' / 

      IF (IDIM.GE.N.AND.N.GT.0) GOTO 310 
      WRITE (0, * ) 'ERROR IN SUBROUTINE CFINV STOP' 
      return
  310 IF (N.EQ.1) RETURN 
      A (2, 1) = - A (2, 2) * DOTF (A (1, 1), A (2, 1), ZERO) 
      A (1, 2) = - A (1, 2) 
      IF (N.EQ.2) GOTO 330 
      DO 314 I = 3, N 
         IM2 = I - 2 
         DO 312 J = 1, IM2 
            S31 = ZERO 
            S32 = A (J, I) 
            DO 311 K = J, IM2 
               S31 = DOTF (A (K, J), A (I, K), S31) 
               S32 = DOTF (A (J, K + 1), A (K + 1, I), S32) 
  311       END DO 
            A (I, J) = - A (I, I) * DOTF (A (I - 1, J), A (I, I - 1),   &
            S31)                                                        
            A (J, I) = - S32 
  312    END DO 
         A (I, I - 1) = - A (I, I) * DOTF (A (I - 1, I - 1), A (I, I -  &
         1), ZERO)                                                      
         A (I - 1, I) = - A (I - 1, I) 
  314 END DO 
  330 NM1 = N - 1 
      DO 335 I = 1, NM1 
         NMI = N - I 
         DO 332 J = 1, I 
            S33 = A (I, J) 
            DO 331 K = 1, NMI 
               S33 = DOTF (A (I + K, J), A (I, I + K), S33) 
  331       END DO 
            A (I, J) = S33 
  332    END DO 
         DO 334 J = 1, NMI 
            S34 = ZERO 
            DO 333 K = J, NMI 
               S34 = DOTF (A (I + K, I + J), A (I, I + K), S34) 
  333       END DO 
            A (I, I + J) = S34 
  334    END DO 
  335 END DO 
      NXCH = IR (N) 
      IF (NXCH.EQ.0) RETURN 
      DO 342 M = 1, NXCH 
         K = NXCH - M + 1 
         IJ = IR (K) 
         I = IJ / 4096 
         J = MOD (IJ, 4096) 
         DO 341 K = 1, N 
            TI = A (K, I) 
            A (K, I) = A (K, J) 
            A (K, J) = TI 
  341    END DO 
  342 END DO 
      RETURN 
      END SUBROUTINE
                                                                        
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
                                                                        
      SUBROUTINE CFACT (N, A, IDIM, IR, IFAIL, DET, JFAIL) 
      IMPLICIT REAL (8)(A - H, O - Z) 
      INTEGER      :: IR(2*N), IPAIRF 
      complex(8)   :: A (IDIM, N), DET, ZERO, ONE, X, Y, TF 
      real(8)      :: G1, G2 
      real(8)      :: PIVOTF, P, Q, SIZEF, T 
      CHARACTER(6) :: HNAME 
      complex(8)   :: S11, S12, DC, DOTF 

      DATA G1, G2 / 1.D-19, 1.D19 / 
      DATA HNAME / ' CFACT' / 
      DATA ZERO, ONE / (0.D0, 0.D0), (1.D0, 0.D0) / 
      DATA NORMAL, IMPOSS / 0, - 1 / 
      DATA JRANGE, JOVER, JUNDER / 0, + 1, - 1 / 

      DOTF (X, Y, S11) = X * Y + S11
      IPAIRF (J, K) = J * 2**12 + K
      PIVOTF (X) = DMAX1 (ABS (REAL (X) ), ABS (aimag (X) ) )
      SIZEF (X) = DMAX1 (ABS (REAL (X) ), ABS (aimag (X) ) )

      IF (IDIM.GE.N.AND.N.GT.0) GOTO 110 
      WRITE (0, * ) 'ERROR IN SUBROUTINE CFACT STOP' 
      return
  110 IFAIL = NORMAL 
      JFAIL = JRANGE 
      NXCH = 0 
      DET = ONE 
      DO 144 J = 1, N 
  120    K = J 
         P = PIVOTF (A (J, J) ) 
         IF (J.EQ.N) GOTO 122 
         JP1 = J + 1 
         DO 121 I = JP1, N 
            Q = PIVOTF (A (I, J) ) 
            IF (Q.LE.P) GOTO 121 
            K = I 
            P = Q 
  121    END DO 
         IF (K.NE.J) GOTO 123 
  122    IF (P.GT.0.) GOTO 130 
         DET = ZERO 
         IFAIL = IMPOSS 
         JFAIL = JRANGE 
         RETURN 
  123    DO 124 L = 1, N 
            TF = A (J, L) 
            A (J, L) = A (K, L) 
            A (K, L) = TF 
  124    END DO 
         NXCH = NXCH + 1 
         IR (NXCH) = IPAIRF (J, K) 
  130    DET = DET * A (J, J) 
         A (J, J) = ONE / A (J, J) 
         T = SIZEF (DET) 
         IF (T.LT.G1) GOTO 132 
         IF (T.GT.G2) GOTO 133 
  131    IF (J.EQ.N) GOTO 144 
         GOTO 140 
  132    DET = ZERO 
         IF (JFAIL.EQ.JRANGE) JFAIL = JUNDER 
         GOTO 131 
  133    DET = ONE 
         IF (JFAIL.EQ.JRANGE) JFAIL = JOVER 
         GOTO 131 
  140    JM1 = J - 1 
         JP1 = J + 1 
         DO 143 K = JP1, N 
            S11 = - A (J, K) 
            S12 = - A (K, J + 1) 
            IF (J.EQ.1) GOTO 142 
            DO 141 I = 1, JM1 
               S11 = DOTF (A (I, K), A (J, I), S11) 
               S12 = DOTF (A (I, J + 1), A (K, I), S12) 
  141       END DO 
  142       A (J, K) = - S11 * A (J, J) 
            A (K, J + 1) = - DOTF (A (J, J + 1), A (K, J), S12) 
  143    END DO 
  144 END DO 
  150 IF (MOD (NXCH, 2) .NE.0) DET = - DET 
      IF (JFAIL.NE.JRANGE) DET = ZERO 
      IR (N) = NXCH 
      RETURN 
      END SUBROUTINE

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_MM(C,A,B,transa,transb, na1, na2, nb1, nb2, nc1, nc2)
  IMPLICIT NONE
  complex(8), intent(out) :: C(nc1,nc2)
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: B(nb1,nb2)
  CHARACTER,  intent(in)  :: transa
  CHARACTER,  intent(in)  :: transb
  INTEGER,    intent(in)  :: na1, na2, nb1, nb2, nc1, nc2
  ! local variables
  complex(8) :: alpha_, beta_
  INTEGER :: m, n, k, kp
  alpha_ = 1.
  beta_ = 0
  IF (transa.EQ.'N' .AND. transb.EQ.'N') THEN
     m = na1
     k = na2
     n = nb2
     kp= nb1
  ELSEIF (transa.EQ.'C' .AND. transb.EQ.'N') THEN
     m = na2
     k = na1
     n = nb2
     kp= nb1
  ELSEIF (transa.EQ.'N' .AND. transb.EQ.'C') THEN
     m = na1
     k = na2
     n = nb1
     kp= nb2
  ELSEIF (transa.EQ.'C' .AND. transb.EQ.'C') THEN
     m = na2
     k = na1
     n = nb1
     kp= nb2
  ELSE
     print *, 'The call to ZProduct_MM was wrong!'
     return
  ENDIF
  IF (kp.NE.k) print *, 'Error in ProductMM, sizes not correct!'
  IF (nc1.NE.m) print *, 'Error in ProductMM, sizes not correct!'
  IF (nc2.NE.n) print *, 'Error in ProductMM, sizes not correct!'
  CALL zgemm(transa,transb,m,n,k,alpha_,A,na1,B,nb1,beta_,C,nc1)
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_sum_MM(C,A,B,transa,transb, na1, na2, nb1, nb2, nc1, nc2)
  IMPLICIT NONE
  complex(8), intent(out) :: C(nc1,nc2)
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: B(nb1,nb2)
  CHARACTER,  intent(in)  :: transa
  CHARACTER,  intent(in)  :: transb
  INTEGER,    intent(in)  :: na1, na2, nb1, nb2, nc1, nc2
  ! local variables
  complex(8) :: alpha_, beta_
  INTEGER :: m, n, k, kp
  alpha_ = 1.
  beta_ = 1.
  IF (transa.EQ.'N' .AND. transb.EQ.'N') THEN
     m = na1
     k = na2
     n = nb2
     kp= nb1
  ELSEIF (transa.EQ.'C' .AND. transb.EQ.'N') THEN
     m = na2
     k = na1
     n = nb2
     kp= nb1
  ELSEIF (transa.EQ.'N' .AND. transb.EQ.'C') THEN
     m = na1
     k = na2
     n = nb1
     kp= nb2
  ELSEIF (transa.EQ.'C' .AND. transb.EQ.'C') THEN
     m = na2
     k = na1
     n = nb1
     kp= nb2
  ELSE
     print *, 'The call to ZProduct_MM was wrong!'
     return
  ENDIF
  IF (kp.NE.k) print *, 'Error in ProductMM, sizes not correct!'
  IF (nc1.NE.m) print *, 'Error in ProductMM, sizes not correct!'
  IF (nc2.NE.n) print *, 'Error in ProductMM, sizes not correct!'
  CALL zgemm(transa,transb,m,n,k,alpha_,A,na1,B,nb1,beta_,C,nc1)
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZTransform(A,U,trans,isize)
  IMPLICIT NONE
  complex(8), intent(inout) :: A(isize,isize)
  complex(8), intent(in)    :: U(isize,isize)
  CHARACTER,  intent(in)    :: trans
  INTEGER,    intent(in)    :: isize
  complex(8) :: temp(isize,isize)  
  IF (trans.EQ.'N') THEN
     CALL Zproduct_MM(temp,A,U,'N','N',isize,isize,isize,isize,isize,isize)
     CALL Zproduct_MM(A,U,temp,'C','N',isize,isize,isize,isize,isize,isize)
  ELSEIF (trans.EQ.'C') THEN
     CALL Zproduct_MM(temp,A,U,'N','C',isize,isize,isize,isize,isize,isize)
     CALL Zproduct_MM(A,U,temp,'N','N',isize,isize,isize,isize,isize,isize)
  ELSE     
     print *, 'The call to ZTransform was wrong!'
     return
  ENDIF
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_NN(C, A, B, na1, na2, nb2)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: B(na2,nb2)
  complex(8), intent(out) :: C(na1,nb2)
  INTEGER,    intent(in)  :: na1, na2, nb2
  ! local variables
  complex(8) :: alpha_, beta_
  alpha_ = 1.
  beta_ = 0
  CALL zgemm('N','N',na1,nb2,na2,alpha_,A,na1,B,na2,beta_,C,na1)
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_NC(C, A, B, na1, na2, nb1)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: B(nb1,na2)
  complex(8), intent(out) :: C(na1,nb1)
  INTEGER,    intent(in)  :: na1, na2, nb1
  ! local variables
  complex(8) :: alpha_, beta_
  alpha_ = 1.
  beta_ = 0
  CALL zgemm('N','C',na1,nb1,na2,alpha_,A,na1,B,nb1,beta_,C,na1)
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_CN(C, A, B, na1, na2, nb2)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: B(na1,nb2)
  complex(8), intent(out) :: C(na2,nb2)
  INTEGER,    intent(in)  :: na1, na2, nb2
  ! local variables
  complex(8) :: alpha_, beta_
  alpha_ = 1.
  beta_ = 0
  CALL zgemm('C','N',na2,nb2,na1,alpha_,A,na1,B,na1,beta_,C,na2)
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_CC(C, A, B, na1, na2, nb1)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: B(nb1,na1)
  complex(8), intent(out) :: C(na2,nb1)
  INTEGER,    intent(in)  :: na1, na2, nb1
  ! local variables
  complex(8) :: alpha_, beta_
  alpha_ = 1.
  beta_ = 0
  CALL zgemm('C','C',na2,nb1,na1,alpha_,A,na1,B,nb1,beta_,C,na2)
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine seigval(w, vr, A, isize)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(isize,isize)
  real(8), intent(out)    :: w(isize)
  complex(8), intent(out) :: vr(isize,isize)
  INTEGER, intent(in)     :: isize
  complex(8), allocatable :: work(:)
  real(8), allocatable    :: rwork(:)
  INTEGER, allocatable    :: iwork(:)
  INTEGER                 :: lwork, lrwork, liwork, info

  lwork  = 2*isize + isize*isize + 5
  lrwork = 5*isize + 2*isize*isize + 5
  liwork = 3 + 5*isize + 5

  ALLOCATE(work(lwork), rwork(lrwork), iwork(liwork))

  vr = A
  CALL ZHEEVD('V', 'U', isize, vr, isize, w, work, lwork, rwork, lrwork, iwork, liwork, info)
  
  IF (info .NE. 0) THEN
     WRITE(0,*) 'ERROR in ZHEEVD inside seigval'
  ENDIF

  DEALLOCATE(work, rwork, iwork)
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine zeigval(zek, A, isize)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(isize,isize)
  complex(8), intent(out) :: zek(isize)
  INTEGER, intent(in)     :: isize
  complex(8), allocatable :: work(:), At(:,:)
  real(8), allocatable    :: rwork(:)
  INTEGER                 :: lwork, lrwork, info
  real(8), PARAMETER      :: smalleps = 1e-5
  complex(8)              :: evl, evr

  lwork  = 3*isize + 1
  lrwork = 2*isize
  
  ALLOCATE(work(lwork), rwork(lrwork), At(isize,isize))

  At = A
  CALL zgeev('N','N',isize,At,isize,zek,evl,isize,evr,isize,work,lwork,rwork,info)

  IF (info .NE. 0) THEN
     WRITE(0,*) 'ERROR in ZGEEV inside seigval'
  ENDIF
  DEALLOCATE(At, work, rwork)
END subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine zeigsys(zek, evl, evr, A, isize)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(isize,isize)
  complex(8), intent(out) :: zek(isize)
  complex(8), intent(out) :: evl(isize,isize)
  complex(8), intent(out) :: evr(isize,isize)
  INTEGER, intent(in)     :: isize
  ! temporaries
  complex(8), allocatable :: work(:), At(:,:)
  real(8), allocatable     :: rwork(:)
  INTEGER    :: lwork, lrwork, info
  real(8), PARAMETER       :: smalleps = 1e-5
  complex(8), allocatable :: scaler(:)
  INTEGER    :: q, p
  complex(8) :: ctmp
  lwork  = 3*isize + 1
  lrwork = 2*isize
  
  ALLOCATE(work(lwork), rwork(lrwork))
  ALLOCATE(At(isize,isize))
  ALLOCATE(scaler(isize))

  At = A
  CALL zgeev('V','V',isize,At,isize,zek,evl,isize,evr,isize,work,lwork,rwork,info)

  IF (info .NE. 0) THEN
     WRITE(0,*) 'ERROR in ZGEEV inside seigval'
  ENDIF
  ! transpose left eigenvectors
  evl = conjg(TRANSPOSE(evl))

  ! Maybe this is not really necessary casue zgeev gives already orthogonal eigensystem
  !========== Step 5, Make degenerate eigenvectors orthogonal
  DO q=1,isize
     DO p=1,q-1
        IF (abs(zek(p)-zek(q)).LT.smalleps .AND. abs(scalprod(evl(p,:),evr(:,q),isize)) .GT.smalleps) THEN
           evr(:,q) = evr(:,q) - scalprod(evl(p,:),evr(:,q),isize)/scalprod(evl(p,:),evr(:,p),isize) * evr(:,p)
        ENDIF
     ENDDO
     DO p=1,q-1
        IF (abs(zek(p)-zek(q)).LT.smalleps .AND. abs(scalprod(evl(q,:),evr(:,p),isize)) .GT.smalleps) THEN
           evl(q,:) = evl(q,:) - scalprod(evl(q,:),evr(:,p),isize)/scalprod(evl(p,:),evr(:,p),isize) * evl(p,:)
        ENDIF
     ENDDO
  ENDDO
  !========= Step 6, Normalize eigenvectors
  DO p = 1,isize
     ctmp = 0.d0
     DO q = 1,isize
        ctmp = ctmp+evl(p,q)*evr(q,p)
     ENDDO
     scaler(p) = SQRT(ctmp)
  ENDDO
  DO p = 1,isize
     evl(p,:) = evl(p,:)/scaler(p)
     evr(:,p) = evr(:,p)/scaler(p)
  ENDDO

  !========== Deallocate dynamic arrays ==========
  DEALLOCATE(scaler)
  DEALLOCATE(At)
  DEALLOCATE(work, rwork)

CONTAINS

  complex(8) FUNCTION scalprod(a,b,ndim)
    IMPLICIT NONE
    complex(8) :: a(:), b(:)
    INTEGER    :: ndim
    INTEGER    :: i
    scalprod = 0.0
    DO i=1,ndim
       scalprod = scalprod + a(i)*b(i)
    ENDDO
  END FUNCTION scalprod

end subroutine 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_MD(C, A, D, na1, na2)
  IMPLICIT NONE
  complex(8), intent(out) :: C(na1,na2)
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: D(na2)
  INTEGER,    intent(in)  :: na1, na2
  ! local variables
  INTEGER :: i, j
  DO i=1,na1
     DO j=1,na2
        C(j,i) = A(j,i)*D(i)
     ENDDO
  ENDDO
end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine check_causality(E, small_positive, isize)
  IMPLICIT NONE
  complex(8), intent(inout) :: E(isize)
  real(8), intent(in)        :: small_positive
  INTEGER, intent(in)       :: isize  
  INTEGER :: i
  real(8)  :: a
  DO i=1,isize
     IF (aimag(E(i)).GT.-small_positive) E(i) = cmplx(dble(E(i)),-small_positive,kind=8)
  ENDDO
end subroutine 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine check_causality2(sig, small_positive, isize)
  IMPLICIT NONE
  complex(8), intent(inout) :: sig(isize,isize)
  real(8), intent(in)        :: small_positive
  INTEGER, intent(in)       :: isize  
  INTEGER :: i, j
  real(8)  :: a
  DO i=1,isize
     IF (aimag(sig(i,i)).GT.-small_positive) sig(i,i) = cmplx(dble(sig(i,i)),-small_positive,kind=8)
  ENDDO
end subroutine 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_ADA(C, A, D, B, na1, na2, nb2)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: D(na2)
  complex(8), intent(in)  :: B(na2,nb2)
  complex(8), intent(out) :: C(na1,nb2)
  INTEGER,    intent(in)  :: na1, na2, nb2
  ! local variables
  INTEGER    :: i, j
  complex(8) :: tmp(na1,na2)
  
  DO j=1,na2
     DO i=1,na1
        tmp(i,j) = A(i,j)*D(j)
     ENDDO
  ENDDO

  CALL ZProduct_NN(C, tmp, B, na1, na2, nb2)

end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

subroutine ZProduct_ADAt(C, A, D, B, na1, na2, nb2, temp)
  IMPLICIT NONE
  complex(8), intent(in)  :: A(na1,na2)
  complex(8), intent(in)  :: D(na2)
  complex(8), intent(in)  :: B(na2,nb2)
  complex(8), intent(out) :: C(na1,nb2)
  INTEGER,    intent(in)  :: na1, na2, nb2
  complex(8), intent(inout) :: temp(na1,na2)
  ! local variables
  INTEGER    :: i, j
  complex(8) :: tmp(na1,na2)
  
  DO j=1,na2
     DO i=1,na1
        tmp(i,j) = A(i,j)*D(j)
     ENDDO
  ENDDO

  CALL ZProduct_NN(C, tmp, B, na1, na2, nb2)

end subroutine 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION dirprod(p,q)
    REAL(DBL), INTENT(IN) :: p(:,:)
    REAL(DBL), INTENT(IN) :: q(:,:)
    REAL(DBL)             :: dirprod(SIZE(p,1)*SIZE(q,1),SIZE(p,2)*SIZE(q,2))
    INTEGER               :: p1,p2,q1,q2
    INTEGER               :: i1,i2
    p1=SIZE(p,1)
    p2=SIZE(p,2)
    q1=SIZE(q,1)
    q2=SIZE(q,2)
    DO i1 = 1, p1
      DO i2 = 1, p2 
        dirprod((i1-1)*q1+1:i1*q1,(i2-1)*q2+1:i2*q2) = p(i1,i2)*q
      END DO
    END DO
  END FUNCTION

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE dirp(pq,p,q)
    REAL(DBL), INTENT(IN)    :: p(:,:)
    REAL(DBL), INTENT(IN)    :: q(:,:)
    REAL(DBL), INTENT(INOUT) :: pq(:,:)
    INTEGER                  :: p1,p2,q1,q2
    INTEGER                  :: i1,i2
    p1=SIZE(p,1)
    p2=SIZE(p,2)
    q1=SIZE(q,1)
    q2=SIZE(q,2)
    DO i1 = 1, p1
      DO i2 = 1, p2 
        pq((i1-1)*q1+1:i1*q1,(i2-1)*q2+1:i2*q2) = p(i1,i2)*q
      END DO
    END DO
  END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE dirpsum(pq,p,q)
    REAL(DBL), INTENT(IN)    :: p(:,:)
    REAL(DBL), INTENT(IN)    :: q(:,:)
    REAL(DBL), INTENT(INOUT) :: pq(:,:)
    INTEGER                  :: p1,p2,q1,q2
    INTEGER                  :: i1,i2
    p1=SIZE(p,1)
    p2=SIZE(p,2)
    q1=SIZE(q,1)
    q2=SIZE(q,2)
    DO i1 = 1, p1
       DO i2 = 1, p2 
        pq((i1-1)*q1+1:i1*q1,(i2-1)*q2+1:i2*q2) = pq((i1-1)*q1+1:i1*q1,(i2-1)*q2+1:i2*q2)+p(i1,i2)*q
       END DO
    END DO
  END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION qnumber(p,q)
    INTEGER, INTENT(IN) :: p(:)
    INTEGER, INTENT(IN) :: q(:)
    INTEGER :: qnumber(SIZE(p)*SIZE(q))
    INTEGER :: p1,q1
    INTEGER :: i1
    p1=SIZE(p)
    q1=SIZE(q)
    DO i1 = 1, p1
       qnumber((i1-1)*q1+1:i1*q1) = p(i1)+q
    END DO
  END FUNCTION

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION outerprod_rs(a,b)
    REAL(4), DIMENSION(:), INTENT(in) :: a,b
    REAL(4), DIMENSION(SIZE(a),SIZE(b)) :: outerprod_rs
    outerprod_rs=SPREAD(a,dim=2,ncopies=SIZE(b))* SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION

  FUNCTION outerprod_r(a,b)
    REAL(8), DIMENSION(:), INTENT(in) :: a,b
    REAL(8), DIMENSION(SIZE(a),SIZE(b)) :: outerprod_r
    outerprod_r=SPREAD(a,dim=2,ncopies=SIZE(b))* SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION outerprod_cs(a,b)
    COMPLEX(4), DIMENSION(:), INTENT(in) :: a,b
    COMPLEX(4), DIMENSION(SIZE(a),SIZE(b)) :: outerprod_cs
    outerprod_cs=SPREAD(a,dim=2,ncopies=SIZE(b))*SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION

  FUNCTION outerprod_c(a,b)
    COMPLEX(8), DIMENSION(:), INTENT(in) :: a,b
    COMPLEX(8), DIMENSION(SIZE(a),SIZE(b)) :: outerprod_c
    outerprod_c=SPREAD(a,dim=2,ncopies=SIZE(b))*SPREAD(b,dim=1,ncopies=SIZE(a))
  END FUNCTION 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

        SUBROUTINE ISCAL(largo,escal,vector,inc)
        integer(4) largo,inc
        integer(4) escal,vector(*)
        do i=1,largo,inc
                vector(i)=escal*vector(i)
        end do
        return
        end subroutine

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

end module
