
module tetraedron_method_tool

  INTEGER :: debug_count

contains

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

complex(8) FUNCTION lv(x,y)

!!!-------------------------------------------------------!!!
!!! We define the Lambin-Vigneron function as follows     !!!
!!!  lv(x,y) = z( 1 - z*Log[ 1 + 1/z ] )  where z = x/y   !!!
!!! or                                                    !!!
!!!  lv(x,y) = ( 1 - Log[ 1 + w ]/w ) /w  where w = y/x.  !!!
!!! For small z or w we use the series expansions:        !!!
!!! lv(x,y) =  z(1+z(log(z)-z(1-z(1/2-z(1/3-z(1/4-z ...   !!!
!!! lv(x,y) =  1/2-w(1/3-w(1/4-w(1/5- ...                 !!!
!!!-------------------------------------------------------!!!

  complex(8), intent(in) :: x
  complex(8), intent(in) :: y
  real(8),     PARAMETER :: eps = 1.d-1
  complex(8)             :: z
  complex(8)             :: w
  complex(8)             :: i2
  complex(8)             :: i3
  complex(8)             :: i4
  complex(8)             :: i5
  complex(8)             :: i6
  complex(8)             :: i7
  complex(8)             :: i8


  i2 = cmplx(1.d0/2.d0,8)
  i3 = cmplx(1.d0/3.d0,8)
  i4 = cmplx(1.d0/4.d0,8)
  i5 = cmplx(1.d0/5.d0,8)
  i6 = cmplx(1.d0/6.d0,8)
  i7 = cmplx(1.d0/7.d0,8)
  i8 = cmplx(1.d0/8.d0,8)

  IF(ABS(y).GT.ABS(x)) THEN
     z = x/y
     lv = z*( 1.d0 - z* (LOG(x+y) - LOG(x)))

     !IF(ABS(z).GT.eps) THEN
     !   lv = z*( 1.d0 - z*LOG( 1 + 1/z ) )
     !ELSE
     !   lv = z*( 1.d0 + z*( LOG(z) - z*( 1.d0 - z*( i2 - z*( i3 - &
     !        z*( i4 - z*( i5 - z*( i6 - z*( i7 - z*i8 ))))))))) 
     !ENDIF
  ELSE
     w = y/x
     IF(ABS(w).GT.eps) THEN
        !lv = ( 1.d0 - LOG( 1 + w )/w )/w
        lv = (1.d0 - (LOG(x+y)-LOG(x))/w )/w
     ELSE
        lv = i2 - w*( i3 - w*( i4 - w*( i5 - w*( i6 - w*(i7 - w*i8 )))))
     ENDIF
  ENDIF

RETURN
END FUNCTION 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
 
SUBROUTINE arbek_optimized( z, e, res )
  IMPLICIT NONE
  complex(8), intent(in)   :: z
  complex(8), intent(in)   :: e(1:4)
  complex(8), intent(out)  :: res(1:4)
  complex(8)               :: lv
  INTEGER                  :: i,j
  complex(8)               :: zmie(1:4)
  complex(8)               :: edif(1:4,1:4)
  complex(8)               :: u(1:4), v1, v2, v3, v4
  complex(8              ) :: p12,p13,p14,p21,p23,p24,p31,p32,p34,p41,p42,p43

  DO i = 1,4
     zmie(i) = z - e(i)
     DO j = 1,4
        edif(j,i) = e(j)-e(i)
     ENDDO
  ENDDO

  u(1) = zmie(1)/(          edif(2,1)*edif(3,1)*edif(4,1))
  u(2) = zmie(2)/(edif(1,2)          *edif(3,2)*edif(4,2))
  u(3) = zmie(3)/(edif(1,3)*edif(2,3)          *edif(4,3))
  u(4) = zmie(4)/(edif(1,4)*edif(2,4)*edif(3,4)          )

  p21 = edif(2,1)*lv(zmie(2),edif(2,1))
  p31 = edif(3,1)*lv(zmie(3),edif(3,1))
  p41 = edif(4,1)*lv(zmie(4),edif(4,1))
  p32 = edif(3,2)*lv(zmie(3),edif(3,2))
  p42 = edif(4,2)*lv(zmie(4),edif(4,2))
  p43 = edif(4,3)*lv(zmie(4),edif(4,3))

  v1 = 1/zmie(1)
  v2 = 1/zmie(2)
  v3 = 1/zmie(3)
  v4 = 1/zmie(4)

  p12 = zmie(1)*(1-zmie(1)*v2*(1-v2*p21))
  p13 = zmie(1)*(1-zmie(1)*v3*(1-v3*p31))
  p14 = zmie(1)*(1-zmie(1)*v4*(1-v4*p41))
  p23 = zmie(2)*(1-zmie(2)*v3*(1-v3*p32))
  p24 = zmie(2)*(1-zmie(2)*v4*(1-v4*p42))
  p34 = zmie(3)*(1-zmie(3)*v4*(1-v4*p43))

  res(1) =           - u(2)*p21 - u(3)*p31 - u(4)*p41
  res(2) = -u(1)*p12            - u(3)*p32 - u(4)*p42
  res(3) = -u(1)*p13 - u(2)*p23            - u(4)*p43
  res(4) = -u(1)*p14 - u(2)*p24 - u(3)*p34
END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE ek2idn( z, e, res )
  complex(8), intent(in)  :: z
  complex(8), intent(in)  :: e(1:4)
  complex(8), intent(out) :: res(1:4)
  complex(8)              :: lv
  INTEGER                 :: i,j
  complex(8)              :: zmie(2:4)
  complex(8)              :: edif(2:4,2:4)
  complex(8)              :: head(2:4)
  complex(8)              :: haus(3:4)
  complex(8)              :: kopf(3:4)

  DO i = 2,4
     zmie(i) = z - e(i)
     DO j = 2,4
        IF(j==i) CYCLE
        edif(i,j) = e(i)-e(j)
     ENDDO
  ENDDO
  head(2) = zmie(2)/(edif(3,2)*edif(4,2))
  head(3) = zmie(3)/(edif(2,3)*edif(4,3))
  head(4) = zmie(4)/(edif(2,4)*edif(3,4))
  haus(3) = zmie(3)/edif(2,3)**2
  haus(4) = zmie(4)/edif(2,4)**2
  kopf(3) = zmie(3)/(edif(4,2)*edif(2,3))
  kopf(4) = zmie(4)/(edif(3,2)*edif(2,4))

  res(2)  = 0.5d0*head(2)  &
       + head(3)*lv(zmie(3),edif(3,2)) &
       + head(4)*lv(zmie(4),edif(4,2)) 

  res(3)  = head(2) &
       + (2*kopf(3)-haus(4))*lv(zmie(2),edif(2,3)) &
       + haus(4)*lv(zmie(4),edif(4,3)) 

  res(4)  = head(2) &
       + (2*kopf(4)-haus(3))*lv(zmie(2),edif(2,4)) &
       + haus(3)*lv(zmie(3),edif(3,4)) 
  res(1)  = res(2)
RETURN
END SUBROUTINE

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE ek22idn( z, e, res )
  complex(8), intent(in)  :: z
  complex(8), intent(in)  :: e(1:4)
  complex(8), intent(out) :: res(1:4)
  complex(8)              :: lv
  INTEGER                 :: i,j
  complex(8)              :: zmie(2:3)
  complex(8)              :: edif(2:3,2:3)
  complex(8)              :: head(2:3)

  DO i = 2,3
     zmie(i) = z-e(i)
     DO j = 2,3
        IF(j==i) CYCLE
        edif(i,j) = e(i)-e(j)
     ENDDO
  ENDDO

  head(2) = zmie(2)/edif(3,2)**2
  head(3) = zmie(3)/edif(2,3)**2
  res(2)  = 0.5d0*head(2)+head(3)-3*head(2)*lv(zmie(3),edif(3,2))
  res(3)  = 0.5d0*head(3)+head(2)-3*head(3)*lv(zmie(2),edif(2,3))
  res(1)  = res(2)
  res(4)  = res(3)
  RETURN
END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE ek3idn( z, e, res )
  complex(8), intent(in)  :: z
  complex(8), intent(in)  :: e(1:4)
  complex(8), intent(out) :: res(1:4)
  complex(8)              :: lv
  INTEGER                 :: i,j
  complex(8)              :: zmie(3:4)
  complex(8)              :: edif(3:4,3:4)
  complex(8)              :: head(3:4)
  complex(8)              :: c1
  complex(8)              :: c5

  c1 = cmplx(1.d0/3.d0,8)
  c5 = cmplx(5.d0/6.d0,8)
  DO i = 3,4
     zmie(i) = z-e(i)
     DO j = 3,4
        IF(j==i) CYCLE
        edif(i,j) = e(i)-e(j)
     ENDDO
  ENDDO
  head(3) = zmie(3)/edif(4,3)**2
  head(4) = zmie(4)/edif(3,4)**2
  res(3)  = c1*head(3)-c5*head(4) + head(4)*lv(zmie(4),edif(4,3))
  res(4)  = 0.5d0*head(3)+head(4)-3*head(3)*lv(zmie(4),edif(4,3))
  res(1)  = res(3)
  res(2)  = res(3)
  RETURN
END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE ek4idn( z, e, res )
  complex(8), intent(in)  :: z
  complex(8), intent(in)  :: e(1:4)
  complex(8), intent(out) :: res(1:4)
  res(:) = 0.25d0/(z-e(1))
  RETURN
END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION ccu(zj, eij, ekj, elj)
    IMPLICIT NONE
    real(8)              :: ccu
    real(8), intent(in)  :: zj, eij, ekj, elj
    ccu = 0
    IF (zj.GT.0 .AND. zj.LT.eij) ccu = 0.25*zj**4/(eij**2*ekj*elj)
    RETURN
  END FUNCTION

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION ccv(zj, eij, ekj)
    IMPLICIT NONE
    real(8)              :: ccv
    real(8), intent(in)  :: zj, eij, ekj
    ccv = 0
    IF (zj.GT.0 .AND. zj.LT.eij) ccv = 0.25*zj**3/(eij**2*ekj)*(4-zj/ekj-2*zj/eij)
    RETURN
  END FUNCTION 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION ccw(zj, eij)
    IMPLICIT NONE
    real(8)              :: ccw
    real(8), intent(in)  :: zj, eij
    ccw = 0
    IF (zj.GT.0 .AND. zj.LT.eij) ccw = 0.25*(zj/eij)**3*(4-3*zj/eij)
    RETURN
  END FUNCTION

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE INIT(z, en, x, e, zi, eij)
    IMPLICIT NONE
    complex(8)  :: z, en(1:4)
    real(8)     :: x, e(1:4), zi(1:4), eij(1:4,1:4)
    INTEGER     :: i, j
    x = real(z)
    DO i=1,4
       e(i) = real(en(i))
    ENDDO
    DO i = 1,4
       zi(i) = x - e(i)
       DO j = 1,4
          eij(j,i) = e(j)-e(i)
       ENDDO
    ENDDO
  END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE arbek_int(z, en, res)
    IMPLICIT NONE
    complex(8) :: z
    complex(8) :: en(1:4)
    complex(8) :: res(1:4)
    real(8)    :: ccu
    real(8)    :: x, e(1:4), zi(1:4), eij(1:4,1:4)

    CALL INIT(z, en, x, e, zi, eij)
    ! Main algorithm
    IF (x.GE.e(1)) THEN 
       res(1)=0.25
    ELSE
       res(1)=0
       res(1) = res(1) + ccu(zi(2), eij(1,2), eij(3,2), eij(4,2))
       res(1) = res(1) + ccu(zi(3), eij(1,3), eij(2,3), eij(4,3))
       res(1) = res(1) + ccu(zi(4), eij(1,4), eij(2,4), eij(3,4))
    ENDIF
    IF (x.GE.e(2)) THEN
       res(2)=0.25
    ELSE
       res(2)=0
       res(2) = res(2) + ccu(zi(3), eij(2,3), eij(4,3), eij(1,3))
       res(2) = res(2) + ccu(zi(4), eij(2,4), eij(3,4), eij(1,4))
       res(2) = res(2) + ccu(zi(1), eij(2,1), eij(3,1), eij(4,1))
    ENDIF
    IF (x.GE.e(3)) THEN
       res(3)=0.25
    ELSE
       res(3)=0
       res(3) = res(3) + ccu(zi(4), eij(3,4), eij(1,4), eij(2,4))
       res(3) = res(3) + ccu(zi(1), eij(3,1), eij(4,1), eij(2,1))
       res(3) = res(3) + ccu(zi(2), eij(3,2), eij(4,2), eij(1,2))
    ENDIF
    IF (x.GE.e(4)) THEN
       res(4)=0.25
    ELSE
       res(4)=0
       res(4) = res(4) + ccu(zi(1), eij(4,1), eij(2,1), eij(3,1))
       res(4) = res(4) + ccu(zi(2), eij(4,2), eij(1,2), eij(3,2))
       res(4) = res(4) + ccu(zi(3), eij(4,3), eij(1,3), eij(2,3))
    ENDIF
  END SUBROUTINE

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE ek2idn_int(z, en, res)
    IMPLICIT NONE
    complex(8) :: z
    complex(8) :: en(1:4)
    complex(8) :: res(1:4)
    real(8)    :: ccu
    real(8)    :: ccv
    real(8)    :: x, e(1:4), zi(1:4), eij(1:4,1:4)

    CALL INIT(z, en, x, e, zi, eij)
    ! Main algorithm
    IF (x.GE.e(1)) THEN 
       res(1)=0.25
    ELSE
       res(1)=0
       res(1) = res(1) + ccu(zi(3), eij(1,3), eij(2,3), eij(4,3))
       res(1) = res(1) + ccu(zi(4), eij(1,4), eij(2,4), eij(3,4))
    ENDIF
    res(2) = res(1)
    IF (x.GE.e(3)) THEN
       res(3)=0.25
    ELSE
       res(3)=0
       res(3) = res(3) + ccu(zi(4), eij(3,4), eij(1,4), eij(2,4))
       res(3) = res(3) + ccv(zi(1), eij(3,1), eij(4,1))
    ENDIF
    IF (x.GE.e(4)) THEN
       res(4)=0.25
    ELSE
       res(4)=0
       res(4) = res(4) + ccv(zi(1), eij(4,1), eij(3,1))
       res(4) = res(4) + ccu(zi(3), eij(4,3), eij(1,3), eij(2,3))
    ENDIF
  END SUBROUTINE

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE ek22idn_int(z, en, res)
    IMPLICIT NONE
    complex(8) :: z
    complex(8) :: en(1:4)
    complex(8) :: res(1:4)
    real(8)    :: ccw
    real(8)    :: x, e(1:4), zi(1:4), eij(1:4,1:4)

    CALL INIT(z, en, x, e, zi, eij)
    ! Main algorithm
    IF (x.GE.e(1)) THEN 
       res(1)=0.25
    ELSE
       res(1) = ccw(zi(3), eij(1,3))
    ENDIF
    res(2) = res(1)
    IF (x.GE.e(3)) THEN
       res(3)=0.25
    ELSE
       res(3) = ccw(zi(1), eij(3,1))
    ENDIF
    res(4) = res(3)
  END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE ek3idn_int(z, en, res)
    IMPLICIT NONE
    complex(8) :: z
    complex(8) :: en(1:4)
    complex(8) :: res(1:4)
    real(8)    :: x, e(1:4), zi(1:4), eij(1:4,1:4), w

    CALL INIT(z, en, x, e, zi, eij)
    ! Main algorithm
    IF (x.GE.e(1)) THEN 
       res(1)=0.25
    ELSE
       res(1) = 0
       IF (zi(4).GT.0 .AND. zi(4).LT.eij(1,4)) res(1) = 0.25*(zi(4)/eij(1,4))**4
    ENDIF
    res(2) = res(1)
    res(3) = res(1)
    IF (x.GE.e(4)) THEN
       res(4)=0.25
    ELSE
       res(4)=0
       w = zi(1)/eij(4,1)
       IF (zi(1).GT.0 .AND. zi(1).LT.eij(4,1)) res(4) = 0.25*w**2*(6-8*w+3*w*2)
    ENDIF
  END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE ek4idn_int(z, en, res)
    IMPLICIT NONE
    complex(8) :: z
    complex(8) :: en(1:4)
    complex(8) :: res(1:4)
    res(1)=0
    IF (real(z).GE.real(en(1))) res(1)=0.25
    res(2) = res(1)
    res(3) = res(1)
    res(4) = res(1)
  END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

!***********************************************************************************
!*** Tetrahedron for complex eigenvalues - integral over frequency and momentum  ***
!***********************************************************************************

  SUBROUTINE terms(a, b, E42, term1, term2, term3)

  !************************************************************************************************
  ! Taylor expansion of the following function:
  ! g(a,b)/(x*E42) - g(a-x,b+x)/(x*(E42-x)) = term1 + term2*x + term3*x^2
  ! where
  ! g(a,b) == int_lv(a,b) = 1/4*y^2*(u^3 + 1/2 u^2 - u - u^4 Log[1 + 1/u] + Log[1 + u]) where u=x/y
  !*************************************************************************************************

    IMPLICIT NONE
    complex(8), intent(out) :: term1, term2, term3
    complex(8), intent(in)  :: a, b, E42
    complex(8)              :: denom_1, nom1_1, nom2_1
    complex(8)              :: lnab, lnab_a, lna
    complex(8)              :: denom_2, nom1_2, nom2_2, nom3_2
    complex(8)              :: denom_3, nom1_3, nom2_3
    complex(8)              :: ab, a2, b2, a3, b3, b4, b5, E42_2, E42_3, E42_4
    complex(8)              :: v0, v1, v2, v3, term1x, term2x, term3x
    real(8),  PARAMETER     :: smallb = 1e-3

    ! Only few logarithms can appear (only those that have correct brach cuts):
    ! log(a), log(a+b)
    ! For example, log(b) should not appear

    lnab = log(a+b)
    lna = log(a)
    lnab_a = lnab-lna
    ab = a*b
    a2 = a**2
    b2 = b**2
    a3 = a2*a
    b3 = b2*b
    b4 = b3*b
    b5 = b4*b
    E42_2 = E42**2
    E42_3 = E42_2*E42
    E42_4 = E42_3*E42

    IF (abs(b).LT.smallb) THEN
       v1 = a*(-3*a + 10*E42)/(12*E42_2)
       v2 = (8*a - 11*E42 - 12*E42*lna)/(24*E42_2)
       v3 = -(5*a + 28*E42 + 20*a*lna)/(80*a*E42_2)
       term1 = v1 + v2*b + v3*b**2

       v1 = -(12*a2 - 40*a*E42 + 31*E42_2 + 12*E42_2*lna)/(48*E42_3)
       v2 = (40*a - 55*E42 - 12*E42_2/a - 60*E42*lna)/(120*E42_3)
       v3 = -(5 + 28*E42/a - 2*E42_2/a2 + 20*lna)/(80*E42_3)
       term2 = v1 + v2*b + v3*b**2

       v0 = 1/(12*E42)
       v1 = (-60*a2/E42_4 + 200*a/E42_3     - 155/E42_2     + 12/(E42*a)  -  60*lna/E42_2)/240.
       v2 = ( 40* a/E42_4 - 55   /E42_3     - 12/(E42_2*a)  - 2 /(E42*a2) -  60*lna/E42_3)/120.
       v3 = (-35   /E42_4 - 196  /(E42_3*a) + 14/(E42_2*a2) + 4 /(E42*a3) - 140*lna/E42_4)/560.
       term3 = v0/b + v1 + v2*b + v3*b**2

       return
    ENDIF

    denom_1 = 8.
    nom1_1 = a*(-a/E42_2 + 6*a/(b*E42) - 2*a2/(b*E42_2) + 4*a2/(b2*E42) + 2*b/(E42_2) + 4/E42)
    nom2_1 =  2*a3*(a/(b2*E42_2) - 2*(a/(b3*E42_2)+2/(b2*E42_2))*E42)*lnab_a - 2*b*(b/E42_2 + 2/E42)*lnab-2*b/E42
    term1 = (nom1_1 + nom2_1)/denom_1

    denom_2 = 8. 
    nom1_2  = a*( 2*(b/E42-1)*(b/E42+3)/(b*E42) - 2*a2*(1/(b*E42_3) - 2/(b2*E42_2) + &
             & 3/(b3*E42)) - a*(1/E42_3 - 6/(b*E42_2) + 13/(b2*E42)))
    nom2_2 = 2*a2*(a2/(b2*E42_3) - 2*a*(a/b + 2)/(b2*E42_2) + (3*a2/b2 + 8*a/b + 6)/(b2*E42))*lnab_a
    nom3_2 = -2*(b/E42 + 1)**2*lnab/(E42) -2*(b/E42 + 1.5)/E42
    term2 = (nom1_2 + nom2_2 + nom3_2)/denom_2
    
    denom_3 = 24.
    nom1_3 = 6*(a/b2)*(a3/E42_4 - 2*a2*(a/b+2)/E42_3 + a*(3*a2/b2+8*a/b+6)/E42_2 - 4*(a/b+1)**3/E42)*(lnab_a)
    nom2_3 = ( -3*a*(2*a2/b + a - 2*b)/E42_4 +  6*a*(2*a2/b2 + 3*a/b + 2)/E42_3 - 3*a*(3*a/b+2)*(2*a/b+3)/(b*E42_2) &
          & + 4*(a/b+1)*(6*a2/b2+9*a/b+2)/(b*E42) - 6*(b/E42+1)**2*(lnab)/E42_2 - (6*b2/E42_2 + 9*b/E42 + 2)/(b*E42))
    term3 = (nom1_3 + nom2_3)/denom_3
    
  END SUBROUTINE

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE termse(a, b, term1, term2, term3)

  !************************************************************************************************
  ! Taylor expansion of the following function:
  ! g(a,b)/(x3*x4) + g(a-x3,b+x3)/(x3*(x3-x4)) + g(a-x4,b+x4)/(x4*(x4-x3))= term1 + term2*(x3+x4) + term3*x3*x4
  ! where
  ! g(a,b) == int_lv(a,b) = 1/4*y^2*(u^3 + 1/2 u^2 - u - u^4 Log[1 + 1/u] + Log[1 + u]) where u=x/y
  !*************************************************************************************************

    IMPLICIT NONE
    complex(8), intent(out) :: term1, term2, term3
    complex(8), intent(in)  :: a, b
    complex(8)              :: lna, lnab_a
    complex(8)              :: u, u2, u3 !, ab, a2, b2, a3, b3, a4, b4, b5, b6
    real(8),  PARAMETER     :: smallb = 1e-3

    lna = log(a)    
    lnab_a = log(a+b)-lna

    u = a/b
    u2 = u*u
    u3 = u2*u
    term1 = ( (6*u3  + 13*u2 + 6*u  + 3) + 2*lna + 2*(1-3*u)*(1+u)**3*lnab_a)/(8.)
    term2 = (-(12*u3 + 30*u2 + 22*u + 3) - 12*u*(u+1)**3*lnab_a)/(12*b)
    term3 = ((2*u+1)*(30*u2+66*u+37) - 12*(u+1)**3*(5*u+1)*lnab_a)/(48*b**2) 
  END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION int_lv(x,y)
    ! u=y/x
    ! int_lv = 0.25*y**2*(u**3+0.5*u**2-u-u**4*ln(1+1/u)+ln(1+u))
    IMPLICIT NONE
    complex(8)             :: int_lv
    complex(8), intent(in) :: x, y
    complex(8)             :: u, v
    complex(8)             :: lnxy_x, lnxy
    real(8)                :: ay, ax
    real(8), PARAMETER     :: small = 1e-3

    u = x/y
    v = y/x
    lnxy_x = LOG(x+y)-LOG(x)
    lnxy = LOG(x+y)
    ay = abs(y)
    ax = abs(x)
    IF (ay.LT.small .and. ax.GT.small) THEN
       int_lv = 0.25*x**2 - x*y/3. + 0.25*y**2*(1./4.-v*(1./5.-v*(1./6.-v*(1./7.-v*(1/8.))))) + 0.25*y**2*lnxy
       RETURN
    ENDIF
    !int_lv = 0.25*y**2(u**3+0.5*u**2-u-u**4*(lnxy-lnx)+lnxy-lny)
    int_lv = 0.25*y**2*(u*(-1. + u*(0.5 + u*(1 - u*(lnxy_x)))) + lnxy)
    RETURN

  END FUNCTION 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION jsign(r)
    IMPLICIT NONE
    real(8) :: jsign
    real(8), intent(in) :: r
    jsign = -1
    IF (r.GT.0) jsign=1
    return
  END FUNCTION

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION zar_case(case, om, E1, E2, E3, E4, prnt) !a, x2, x3, x4, prnt)

    ! a=om-E2
    ! x2 = E2-E1
    ! x3 = E3-E2
    ! x4 = E4-E2
    ! Computes function
    ! g(a,x2)/(x3*x4) + g(a-x3,x2+x3)/(x3*(x3-x4)) + g(a-x4,x2+x4)/(x4*(x4-x3))
    ! which corresponds to the tetrahedron integral

    IMPLICIT NONE

    complex(8)             :: zar_case
    INTEGER, intent(in)    :: case
    complex(8), intent(in) :: om, E1, E2, E3, E4 !a, x2, x3, x4
    LOGICAL, intent(in)    :: prnt
    complex(8)             :: int_lv
    complex(8)             :: a, x2, x3, x4
    complex(8)             :: w1,w2,w3, term1,term2,term3
    real(8),PARAMETER      :: pi = 3.14159266!358979323844
    real(8), PARAMETER     :: Dw = 4.
    real(8)                :: omDw, maxaimag
    
    maxaimag = MAX(MAX(MAX(abs(aimag(E1)),abs(aimag(E2))),abs(aimag(E3))),abs(aimag(E4)))
    !print *, 'maxaimag=', maxaimag

    omDw = real(om)+Dw+maxaimag*10
    IF (omDw<dble(E1) .AND. omDw<dble(E2) .AND. omDw<dble(E3) .AND. omDw<dble(E4)) THEN ! empty
       zar_case = cmplx(0,0.25*pi,8)
       RETURN 
    ENDIF
    omDw = dble(om)-Dw-maxaimag*10
    IF (omDw>dble(E1) .AND. omDw>dble(E2) .AND. omDw>dble(E3) .AND. omDw>dble(E4)) THEN ! full
       zar_case = 0
       RETURN
    ENDIF

    a = om-E2
    x2 = E2-E1
    x3 = E3-E2
    x4 = E4-E2
    
    IF (case.EQ.1) THEN    ! valid when min(x2,x3,x4)>>1e-4 or x2<x3<x4 or x2<x4<x3
       w1 = int_lv(a,x2)
       w2 = int_lv(a-x3,x2+x3)
       w3 = int_lv(a-x4,x2+x4)
       zar_case = w1/(x3*x4)+w2/(x3*(x3-x4))+w3/(x4*(x4-x3))
       RETURN
    ELSEIF (case.EQ.2) THEN ! valid when x3<<x2<x4
       CALL terms(a, x2, x4, term1, term2, term3)
       w3 = int_lv(a-x4,x2+x4)
       zar_case = term1 + term2*x3 + term3*x3**2 + w3/(x4*(x4-x3))
       RETURN
    ELSEIF (case.EQ.3) THEN ! valid when x4<<x2<x3
       CALL terms(a, x2, x3, term1, term2, term3)
       w2 = int_lv(a-x3,x2+x3)
       zar_case = term1 + term2*x4 + term3*x4**2 + w2/(x3*(x3-x4))
       RETURN
    ELSEIF (case.EQ.4) THEN ! valid when |x4-x3|<<|x2+x3|<|x3|
       CALL terms(a-x3, x2+x3, -x3, term1, term2, term3)
       w1 = int_lv(a,x2)
       zar_case = term1 + term2*(x4-x3) + term3*(x4-x3)**2 + w1/(x3*x4)
       RETURN
    ELSEIF (case.EQ.5) THEN ! valid when (x3,x4) << x2
       CALL termse(a, x2, term1, term2, term3)
       zar_case = term1 + term2*(x3+x4) + term3*x3*x4 
       RETURN
    ELSE
       print *, 'THIS CASE DOES NOT EXIST'
    ENDIF
  END FUNCTION 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
 
  FUNCTION zarbek_r1(om1,om2,E1,E2,E3,E4,prnt)

    ! Computes the function
    ! g(a,x2)/(x3*x4) + g(a-x3,x2+x3)/(x3*(x3-x4)) + g(a-x4,x2+x4)/(x4*(x4-x3))
    ! and takes care of various singular limits. Calls above functions 
    ! which give Taylor expansion in various limits

    IMPLICIT NONE
    complex(8)               :: zarbek_r1
    complex(8), intent(in)   :: om1,om2,E1,E2,E3,E4
    LOGICAL, intent(in)      :: prnt
    complex(8)               :: int_lv, zar_case!1, zar_case2, zar_case3, zar_case4, zar_case5
    real(8), PARAMETER       :: smallE = 1e-4
    real(8), PARAMETER       :: smallEi = 1e-3
    complex(8)               :: w1, w2, w3, xi, term1_a, term2_a, term3_a, term1_b, term2_b, term3_b
    complex(8)               :: zarbek_r1x, E32, E42, v1, v2, v3
    real(8)                  :: Estr(4)
    INTEGER                  :: iesrt(4), i

    debug_count = debug_count + 1

    Estr(1) = abs(E3-E2) ! x3
    Estr(2) = abs(E4-E2) ! x4
    Estr(3) = abs(E3-E4) ! x3-x4
    Estr(4) = abs(E2-E1) ! x2
    DO i=1,4
       iesrt(i)=i
    ENDDO
    CALL ssort(Estr, iesrt, 4)

    IF (prnt) THEN
       print *, '-------------------------'
       print *, debug_count, om2-E2, om1-E2, E2-E1, E3-E2, E4-E2
       print *, iesrt
    ENDIF

    IF (abs(E3-E2)<smallEi .AND. abs(E4-E2)<smallEi) THEN
       IF (iesrt(1)+iesrt(2).LE.4) THEN
          zarbek_r1 = zar_case(5, om2, E1, E2, E3, E4, prnt) - zar_case(5, om1, E1, E2, E3, E4, prnt)
          if (prnt) print *, 'v5', debug_count, zarbek_r1
          RETURN
       ENDIF
       IF (Estr(1)/Estr(4)<smallEi .AND. Estr(1)/Estr(2)<smallEi .AND. iesrt(1).EQ.1) THEN
          zarbek_r1 = zar_case(2, om2, E1, E2, E3, E4, prnt)-zar_case(2, om1, E1, E2, E3, E4, prnt)
          if (prnt) print *, 'v2', debug_count, zarbek_r1
          RETURN
       ENDIF
       IF (Estr(2)/Estr(4)<smallEi .AND. Estr(2)/Estr(1)<smallEi .AND. iesrt(1).EQ.2) THEN
          zarbek_r1 = zar_case(3, om2, E1, E2, E3, E4, prnt)-zar_case(3, om1, E1, E2, E3, E4, prnt)
          if (prnt) print *, 'v3', debug_count, zarbek_r1
          RETURN
       ENDIF
       IF (Estr(3)/Estr(4)<smallEi .AND. Estr(3)/Estr(1)<smallEi .AND. iesrt(1).EQ.3) THEN
          zarbek_r1 = zar_case(4, om2, E1, E2, E3, E4, prnt)-zar_case(4, om1, E1, E2, E3, E4, prnt)
          if (prnt) print *, 'v4', debug_count, zarbek_r1
          RETURN
       ENDIF
       IF (iesrt(1).EQ.4) THEN
          zarbek_r1 = zar_case(1, om2, E1, E2, E3, E4, prnt)-zar_case(1, om1, E1, E2, E3, E4, prnt)
          if (prnt) print *, 'v1', debug_count, zarbek_r1
          RETURN
       ENDIF
    ENDIF
    IF (Estr(1)/Estr(4)<smallEi .AND. Estr(1)/Estr(2)<smallEi .AND. iesrt(1).EQ.1) THEN       
       zarbek_r1 = zar_case(2, om2, E1, E2, E3, E4, prnt)-zar_case(2, om1, E1, E2, E3, E4, prnt)
       if (prnt) print *, 'v2', zarbek_r1
      RETURN
    ENDIF
    IF (Estr(2)/Estr(4)<smallEi .AND. Estr(2)/Estr(1)<smallEi .AND. iesrt(1).EQ.2) THEN
       zarbek_r1 = zar_case(3, om2, E1, E2, E3, E4, prnt)-zar_case(3, om1, E1, E2, E3, E4, prnt)
       if (prnt) print *, 'v3', debug_count, zarbek_r1
       RETURN
    ENDIF
    IF (Estr(3)/Estr(4)<smallEi .AND. Estr(3)/Estr(1)<smallEi .AND. iesrt(1).EQ.3) THEN
       zarbek_r1 = zar_case(4, om2, E1, E2, E3, E4, prnt)-zar_case(4, om1, E1, E2, E3, E4, prnt)
       if (prnt) print *, 'v4', debug_count, zarbek_r1
       RETURN
    ENDIF
    
    zarbek_r1 = zar_case(1, om2, E1, E2, E3, E4, prnt)-zar_case(1, om1, E1, E2, E3, E4, prnt)

    if (prnt) print *, 'v1', debug_count, zarbek_r1, zar_case(1, om2, E1, E2, E3, E4, prnt), zar_case(1, om1, E1, E2, E3, E4, prnt)
    RETURN
  END FUNCTION 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE zarbek_int(om1, om2, en, res)
    IMPLICIT NONE
    complex(8), intent(in)  :: om1, om2
    complex(8), intent(in)  :: en(1:4)
    complex(8), intent(out) :: res(1:4)
    complex(8)              :: zarbek_r1
    res(1) = zarbek_r1(om1,om2, en(1), en(2), en(3), en(4), .FALSE.)
    res(2) = zarbek_r1(om1,om2, en(2), en(3), en(4), en(1), .FALSE.)
    res(3) = zarbek_r1(om1,om2, en(3), en(4), en(1), en(2), .FALSE.)
    res(4) = zarbek_r1(om1,om2, en(4), en(1), en(2), en(3), .FALSE.)
  END SUBROUTINE 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE zarbek_int_debug(om1, om2, en, res)
    IMPLICIT NONE
    complex(8), intent(in)  :: om1, om2
    complex(8), intent(in)  :: en(1:4)
    complex(8), intent(out) :: res(1:4)
    complex(8)              :: zarbek_r1

    res(1) = zarbek_r1(om1,om2, en(1), en(2), en(3), en(4), .TRUE.)
    res(2) = zarbek_r1(om1,om2, en(2), en(3), en(4), en(1), .TRUE.)
    res(3) = zarbek_r1(om1,om2, en(3), en(4), en(1), en(2), .TRUE.)
    res(4) = zarbek_r1(om1,om2, en(4), en(1), en(2), en(3), .TRUE.)
  END SUBROUTINE

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  FUNCTION MINIMUM(x)
    IMPLICIT NONE
    real(8)                :: MINIMUM
    complex(8), intent(in) :: x(4)
    real(8)                :: xmin
    INTEGER                :: i
    xmin = real(x(1))
    DO i=2,4
       IF (real(x(i)).LT.xmin) xmin=real(x(i))
    ENDDO
    MINIMUM = xmin
    RETURN
  END FUNCTION

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
  
  FUNCTION MAXIMUM(x)
    IMPLICIT NONE
    real(8)                :: MAXIMUM
    complex(8), intent(in) :: x(4)
    real(8)                :: xmax
    INTEGER                :: i
    xmax = real(x(1))
    DO i=2,4
       IF (real(x(i)).GT.xmax) xmax=real(x(i))
    ENDDO
    MAXIMUM = xmax
    RETURN
  END FUNCTION 

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************


end module
