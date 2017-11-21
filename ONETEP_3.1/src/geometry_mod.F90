! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   The subroutines in this file were written by Chris-Kriton Skylaris
!
!   February 2000
!   
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE GEOMETRY

  use CONSTANTS, only: DP

  IMPLICIT NONE 

  private

  TYPE POINT
     real(kind=DP) :: X,Y,Z
  END TYPE POINT


  INTERFACE OPERATOR(+)
     MODULE PROCEDURE ADD_POINTS
  END INTERFACE

  INTERFACE OPERATOR(-)
     MODULE PROCEDURE SUBTRACT_POINTS
  END INTERFACE

  INTERFACE OPERATOR(*)
     MODULE PROCEDURE SCALE_POINT_REAL
     MODULE PROCEDURE SCALE_POINT_INT
  END INTERFACE

  INTERFACE OPERATOR(.CROSS.)
     MODULE PROCEDURE CROSS_PRODUCT
  END INTERFACE

  INTERFACE OPERATOR(.DOT.)
     MODULE PROCEDURE DOT_PROD
  END INTERFACE


  public :: POINT

  public :: geometry_distance
  public :: local_displacement
  public :: magnitude
  public :: unit_vector
  public :: geometry_magnitude

  public :: operator(+), operator(-), operator(*)
  public :: operator(.CROSS.) 
  public :: operator(.DOT.) 

CONTAINS  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TYPE (POINT) FUNCTION UNIT_VECTOR(A)

    IMPLICIT NONE

    TYPE(POINT), intent(in) :: A

    UNIT_VECTOR=(1.0_DP/(MAGNITUDE(A)))*A


  END FUNCTION UNIT_VECTOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  FUNCTION geometry_MAGNITUDE(A)

    use constants, only: DP
    IMPLICIT NONE

    real(kind=DP) :: geometry_magnitude

    TYPE(POINT), intent(in) :: A

    geometry_MAGNITUDE=SQRT(A%X**2+A%Y**2+A%Z**2)

  END FUNCTION GEOMETRY_MAGNITUDE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  FUNCTION MAGNITUDE(A)

    use constants, only: DP
    IMPLICIT NONE

    real(kind=DP) :: magnitude

    TYPE(POINT), intent(in) :: A

    MAGNITUDE=SQRT(A%X**2+A%Y**2+A%Z**2)

  END FUNCTION MAGNITUDE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  FUNCTION geometry_DISTANCE(point1,point2)

    use constants, only: DP
    IMPLICIT NONE

    real(kind=DP) :: geometry_distance

    type(POINT), intent(in) :: point1,point2


    geometry_DISTANCE=sqrt((point1%x-point2%x)**2 & 
         + (point1%y-point2%y)**2 + (point1%z-point2%z)**2 )

  END FUNCTION GEOMETRY_DISTANCE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  TYPE(POINT) FUNCTION LOCAL_DISPLACEMENT(&
       A1_UNIT,A2_UNIT,A3_UNIT,LOCAL_1,LOCAL_2,LOCAL_3)

    IMPLICIT NONE 

    TYPE(POINT),   intent(in) :: A1_UNIT,A2_UNIT,A3_UNIT

    real(kind=DP), intent(in) :: LOCAL_1,LOCAL_2,LOCAL_3


    LOCAL_DISPLACEMENT=LOCAL_1*A1_UNIT+LOCAL_2*A2_UNIT+LOCAL_3*A3_UNIT



  END FUNCTION LOCAL_DISPLACEMENT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  TYPE(POINT) FUNCTION NORMAL_UNIT_VECTOR(A,B)

    IMPLICIT NONE

    TYPE(POINT), intent(in) :: A,B
    TYPE(POINT)             :: NORMAL

    ! cks: NORMAL_UNIT_VECTOR = AxB/|AxB|

    NORMAL%X=A%Y*B%Z-A%Z*B%Y
    NORMAL%Y=A%Z*B%X-A%X*B%Z
    NORMAL%Z=A%X*B%Y-A%Y*B%X


    NORMAL_UNIT_VECTOR=UNIT_VECTOR(NORMAL)


  END FUNCTION NORMAL_UNIT_VECTOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! cks: for definition of new operators

  TYPE(POINT) FUNCTION ADD_POINTS(A,B)

    IMPLICIT NONE

    TYPE(POINT), INTENT(IN) :: A,B

    ADD_POINTS%X=A%X+B%X
    ADD_POINTS%Y=A%Y+B%Y
    ADD_POINTS%Z=A%Z+B%Z

  END FUNCTION ADD_POINTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  TYPE(POINT) FUNCTION SUBTRACT_POINTS(A,B)

    IMPLICIT NONE

    TYPE(POINT), INTENT(IN) :: A,B

    SUBTRACT_POINTS%X=A%X-B%X
    SUBTRACT_POINTS%Y=A%Y-B%Y
    SUBTRACT_POINTS%Z=A%Z-B%Z

  END FUNCTION SUBTRACT_POINTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  TYPE(POINT) FUNCTION SCALE_POINT_REAL(SCALE,A)

    use constants, only: DP
    IMPLICIT NONE

    real(kind=DP), INTENT(IN) :: SCALE

    TYPE(POINT), INTENT(IN) :: A

    SCALE_POINT_REAL%X=SCALE*A%X
    SCALE_POINT_REAL%Y=SCALE*A%Y
    SCALE_POINT_REAL%Z=SCALE*A%Z

  END FUNCTION SCALE_POINT_REAL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  TYPE(POINT) FUNCTION SCALE_POINT_INT(SCALE,A)

    use constants, only: DP
    IMPLICIT NONE

    integer, INTENT(IN) :: SCALE

    TYPE(POINT), INTENT(IN) :: A

    SCALE_POINT_INT%X=real(SCALE,kind=DP)*A%X
    SCALE_POINT_INT%Y=real(SCALE,kind=DP)*A%Y
    SCALE_POINT_INT%Z=real(SCALE,kind=DP)*A%Z

  END FUNCTION SCALE_POINT_INT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  TYPE(POINT) FUNCTION CROSS_PRODUCT(A,B)

    IMPLICIT NONE

    TYPE(POINT), INTENT(IN) :: A,B

    CROSS_PRODUCT%X=A%Y*B%Z-A%Z*B%Y
    CROSS_PRODUCT%Y=A%Z*B%X-A%X*B%Z
    CROSS_PRODUCT%Z=A%X*B%Y-A%Y*B%X

  END FUNCTION CROSS_PRODUCT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(kind=DP) FUNCTION DOT_PROD(A,B)

    IMPLICIT NONE

    TYPE(POINT), INTENT(IN) :: A,B

    DOT_PROD=A%X*B%X+A%Y*B%Y+A%Z*B%Z

  END FUNCTION DOT_PROD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE GEOMETRY












