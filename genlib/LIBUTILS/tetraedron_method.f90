module tetraedron

use tetraedron_method_tool


contains



 SUBROUTINE tetrahedron( z, e, weights) ! former int_tet1

 !!!--------------------------------------------------!!!
 !!!  These subroutines compute the weight factors for!!!
 !!!  the analytical tetrahedron method as described  !!!
 !!!  in the paper by Lambin and Vigneron, Physical   !!!
 !!!  Review B, Vol. 29,#6, page 3430.  The formulas  !!!
 !!!  used here are written slightly differently from !!!
 !!!  the ones given in the article and hopefully in  !!!
 !!!  a more stable form.  In particular we try to    !!!
 !!!  express the answer as much as possible in a new !!!
 !!!  function, the socalled Lambin-Vigneron function,!!!
 !!!  lv, which is easy to compute and has quite nice !!!
 !!!  asymptotic limits.                              !!!
 !!!--------------------------------------------------!!!

  complex(8), intent(out)  :: weights(1:4)
  complex(8), INTENT(in)   :: z           ! ( In ) Frequency
  complex(8), INTENT(in)   :: e(4)        ! ( In ) Energy eigenvalues
  INTEGER                  :: flg         ! Corner degeneracy index
  INTEGER                  :: isrt(1:4)   ! Corner reordering index
  complex(8)               :: ecp(1:4)    ! Local copy of energy values
  complex(8)               :: res(1:4)    ! Linear weight factors 



  ! Make local copy of corner energies
  ecp(:) = e(:)
  ! Reorder the corner energies
  CALL creorder( ecp,isrt,flg )
  ! Select the appropriate routine to compute the integral
  SELECT CASE ( flg )
  CASE(1)
     CALL arbek_optimized( z,ecp,res )     !  Arbitrary eigenvalues
  CASE(2)
     CALL ek2idn( z,ecp,res )    !  Two eigenvalues are identical
  CASE(3)
     CALL ek22idn( z,ecp,res )   !  Two pairs of energy are identical
  CASE(4)
     CALL ek3idn( z,ecp,res )    !  Three eigenvalues are identical
  CASE(5)
     CALL ek4idn( z,ecp,res )    !  Four (all) eigenvalues are identical
  END SELECT
  ! Put integration weights in correct order
  weights(isrt(1)) = res(1)
  weights(isrt(2)) = res(2)
  weights(isrt(3)) = res(3)
  weights(isrt(4)) = res(4)
END SUBROUTINE







  SUBROUTINE tetrahedron_int(z, e, weights) ! former int_tet1_int

!!!--------------------------------------------------!!!
!!!  These subroutines compute the weight factors for!!!
!!!  the analytical tetrahedron method as described  !!!
!!!  in the paper by Lambin and Vigneron, Physical   !!!
!!!  Review B, Vol. 29,#6, page 3430.  The formulas  !!!
!!!  used here are written slightly differently from !!!
!!!  the ones given in the article and hopefully in  !!!
!!!  a more stable form.  In particular we try to    !!!
!!!  express the answer as much as possible in a new !!!
!!!  function, the socalled Lambin-Vigneron function,!!!
!!!  lv, which is easy to compute and has quite nice !!!
!!!  asymptotic limits.                              !!!
!!!--------------------------------------------------!!!

    IMPLICIT NONE
    complex(8), intent(out) :: weights(1:4)
    complex(8), INTENT(in)  :: z         ! ( In ) Frequency
    complex(8), INTENT(in)  :: e(4)      ! ( In ) Energy eigenvalues
    real(8)                 :: MINIMUM, MAXIMUM
    INTEGER                 :: flg         ! Corner degeneracy index
    INTEGER                 :: isrt(1:4)   ! Corner reordering index
    complex(8)              :: ecp(1:4)    ! Local copy of energy values
    complex(8)              :: res(1:4)    ! Linear weight factors 
    real(8)                 :: emin, emax
    real(8),PARAMETER       :: pi = 3.14159266!358979323844
    complex(8)              :: prefact

    ! Trivial case that happens most of the time
    ! weight is either unity or zero

    prefact = -pi*cmplx(0,1.,8)
    emin = MINIMUM(e)
    emax = MAXIMUM(e)
    IF (real(z).GT.emax) THEN
       weights = 0.25*prefact
       RETURN
    ELSEIF (real(z).LT.emin) THEN
       weights = 0
       RETURN
    ENDIF
    ecp(:) = e(:)
    ! Reorder the corner energies
    CALL creorder( ecp,isrt,flg )
    ! Select the appropriate routine to compute the integral
    SELECT CASE ( flg )
    CASE(1)
       CALL arbek_int( z,ecp,res )     !  Arbitrary eigenvalues
    CASE(2)
       CALL ek2idn_int( z,ecp,res )    !  Two eigenvalues are identical
    CASE(3)
       CALL ek22idn_int( z,ecp,res )   !  Two pairs of energy are identical
    CASE(4)
       CALL ek3idn_int( z,ecp,res )    !  Three eigenvalues are identical
    CASE(5)
       CALL ek4idn_int( z,ecp,res )    !  Four (all) eigenvalues are identical
    END SELECT
    ! Put integration weights in correct order
    weights(isrt(1)) = res(1)*prefact
    weights(isrt(2)) = res(2)*prefact
    weights(isrt(3)) = res(3)*prefact
    weights(isrt(4)) = res(4)*prefact
  END SUBROUTINE








  SUBROUTINE ztetrahedron_int( om1, om2, e, weights)

!!!--------------------------------------------------!!!
!!!  These subroutines compute the weight factors for!!!
!!!  the analytical tetrahedron method -  integral   !!!
!!!  over momentum and frequency for complex         !!!
!!!  energies.                                       !!!
!!!  integral[1/(om-ek),{k in tetra},{om,om1,om2}]   !!!
!!!--------------------------------------------------!!!

    IMPLICIT NONE
    complex(8), intent(out) :: weights(1:4)
    complex(8), INTENT(in)  :: om1,om2   ! Frequencies
    complex(8), INTENT(in)  :: e(4)      ! Energy eigenvalues
    real(8)                 :: jsign
    INTEGER                 :: flg         ! Corner degeneracy index
    INTEGER                 :: isrt(1:4)   ! Corner reordering index
    complex(8)              :: ecp(1:4)    ! Local copy of energy values
    complex(8)              :: res(1:4)    ! Linear weight factors
    real(8)                 :: emin, emax, maxi, maxi1
    real(8),PARAMETER       :: pi = 3.14159266!358979323844
    complex(8), PARAMETER   :: small = cmplx(1e-6,0,8)
    INTEGER                 :: i

    ecp(:) = e(:)
    ! Reorder the corner energies

    CALL creorder( ecp,isrt,flg )

    ! Select the appropriate routine to compute the integral
    SELECT CASE ( flg )
    CASE(1)
       CALL zarbek_int( om1, om2, ecp, res )     !  Arbitrary eigenvalues
    CASE(2)
       ecp(2)=ecp(1)+small
       CALL zarbek_int( om1, om2, ecp, res )    !  Two eigenvalues are identical
    CASE(3)
       ecp(2)=ecp(1)+small
       ecp(4)=ecp(3)+small
       CALL zarbek_int( om1, om2, ecp, res )   !  Two pairs of energy are identical
    CASE(4)
       ecp(2)=ecp(1)+small
       ecp(3)=ecp(1)+2*small
       CALL zarbek_int( om1, om2, ecp, res )    !  Three eigenvalues are identical
    CASE(5)
       res(1) = 0.25*(log(om2-ecp(1))-log(om1-ecp(1))) !  Four (all) eigenvalues are identical
       res(2) = res(1)
       res(3) = res(1)
       res(4) = res(1)
    END SELECT
    ! Checking weather results are not exciding their upper limit
    maxi = 0.25*pi ! weight can not be more than that. Due to numerics, it can happen that it is more
    maxi1 = maxi+0.1
    
    IF (aimag(res(1)).GT.0) THEN
       !if (aimag(res(1)).GT.0.05) then
       !   print *, 'cs5', res(1), ecp(:), om2, flg
       !   CALL zarbek_int_debug(om1, om2, ecp, res)
       !endif
       res(1) = 0
    ENDIF
    IF (aimag(res(2)).GT.0) THEN
       !if (aimag(res(1)).GT.0.05) then
       !   print *, 'cs6', res(2), ecp(:), om2, flg
       !   CALL zarbek_int_debug(om1, om2, ecp, res)
       !endif
       res(2) = 0
    ENDIF
    IF (aimag(res(3)).GT.0) THEN
       !if (aimag(res(1)).GT.0.05) then
       !   print *, 'cs7', res(3), ecp(:), om2, flg
       !   CALL zarbek_int_debug(om1, om2, ecp, res)
       !endif
       res(3) = 0
    ENDIF
    IF (aimag(res(4)).GT.0) THEN
       !if (aimag(res(1)).GT.0.05) then
       !   print *, 'cs8', res(4), ecp(:), om2, flg
       !   CALL zarbek_int_debug(om1, om2, ecp, res)
       !endif
       res(4) = 0
    ENDIF

    IF ((-aimag(res(1))).GT.maxi1) THEN 
       !print *, 'cs1', res(1), ecp(:), flg
       !CALL zarbek_int_debug(om1, om2, ecp, res)
       res(1) = cmplx(0, -maxi,8)
    ENDIF
    IF ((-aimag(res(2))).GT.maxi1) THEN
       !print *, 'cs2', res(2), ecp(:), flg
       !CALL zarbek_int_debug(om1, om2, ecp, res)
       res(2) = cmplx(0, -maxi,8)
    ENDIF
    IF ((-aimag(res(3))).GT.maxi1) THEN
       !print *, 'cs3', res(3), ecp(:), flg
       !CALL zarbek_int_debug(om1, om2, ecp, res)
       res(3) = cmplx(0, -maxi,8)
    ENDIF
    IF ((-aimag(res(4))).GT.maxi1) THEN
       !print *, 'cs4', res(4), ecp(:), flg
       !CALL zarbek_int_debug(om1, om2, ecp, res)
       res(4) = cmplx(0, -maxi,8)
    ENDIF
    ! Put integration weights in correct order
    weights(isrt(1)) = res(1)
    weights(isrt(2)) = res(2)
    weights(isrt(3)) = res(3)
    weights(isrt(4)) = res(4)
    RETURN
  END SUBROUTINE







end module

