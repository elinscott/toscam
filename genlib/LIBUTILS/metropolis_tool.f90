MODULE com_metrop

use random

!BUG September 5th, next line is wrong:
! PRIVATE :: 
PRIVATE

  INTEGER                  :: stranica(2,6), istranica(3,4)
  real(8),     ALLOCATABLE  :: dst(:,:,:), dst0(:,:)
  INTEGER,    ALLOCATABLE  :: indx(:,:), indx0(:)
! INTEGER,    SAVE         :: indim
! PUBLIC  :: Distance, TDistance, 
  PUBLIC  :: Metropolis, TMetropolis, Init_VMetrop,Init_TMetrop, Destroy_Metrop
! PRIVATE :: ran, indim

contains

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

  SUBROUTINE Init_Vmetrop(ndim)
    IMPLICIT NONE
    INTEGER :: ndim
!    indim = ndim
    ALLOCATE(dst0(ndim,ndim), indx0(ndim))
  END SUBROUTINE Init_Vmetrop

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

  SUBROUTINE Init_Tmetrop(ndim)
    IMPLICIT NONE
    INTEGER :: ndim
!    indim = ndim
    stranica(1,1)=1
    stranica(2,1)=2
    stranica(1,2)=1
    stranica(2,2)=3
    stranica(1,3)=1
    stranica(2,3)=4
    stranica(1,4)=2
    stranica(2,4)=3
    stranica(1,5)=2
    stranica(2,5)=4
    stranica(1,6)=3
    stranica(2,6)=4
    istranica(1,1)=1
    istranica(2,1)=2
    istranica(3,1)=3
    istranica(1,2)=1
    istranica(2,2)=4
    istranica(3,2)=5
    istranica(1,3)=2
    istranica(2,3)=4
    istranica(3,3)=6
    istranica(1,4)=3
    istranica(2,4)=5
    istranica(3,4)=6
    ALLOCATE(dst(ndim,ndim,6), indx(ndim,4))
  END SUBROUTINE Init_Tmetrop

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

  SUBROUTINE Destroy_metrop
    IMPLICIT NONE
    IF (ALLOCATED (dst)  ) DEALLOCATE (dst)
    IF (ALLOCATED (indx) ) DEALLOCATE (indx)
    IF (ALLOCATED (dst0) ) DEALLOCATE (dst0)
    IF (ALLOCATED (indx0)) DEALLOCATE (indx0)
  END SUBROUTINE Destroy_metrop

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

  real(8) FUNCTION Metropolis(v0, v1, max_steps, indim)
    IMPLICIT NONE
    ! Input variables
    complex(8), intent(in) :: v0(indim), v1(indim)
    INTEGER, intent(in)    :: max_steps
    INTEGER, intent(in)    :: indim
    ! Local variables
    real(8), PARAMETER :: small = 1e-6
    INTEGER    :: i, j
    real(8)     :: distanc, distanc0, dDistance
    INTEGER    :: a, b, ia, ib, da, direction
    INTEGER    :: idum
    !  print *, max_steps
    DO i=1,indim
       DO j=1,indim
          dst0(i,j) = abs(v0(i)-v1(j))
       ENDDO
    ENDDO
    DO i=1,indim
       indx0(i) = i
    ENDDO
    distanc = 0.0
    DO i=1,indim
       distanc = distanc + dst0(i,i)
    ENDDO
    distanc0 = distanc
    idum = -1
    DO i=1,max_steps
       a = int(ran(idum)*indim) + 1

       da = int(-5*log(1-ran(idum)))+1
       direction = int(2*ran(idum))
       direction = 2*direction-1
       b = a+direction*da
       IF (b.LT.1) b=1
       IF (b.GT.indim) b=indim

       IF (a.EQ.b) THEN
          IF (a.GT.1) THEN 
             b=a-1
          ELSE
             b=a+1
          ENDIF
       ENDIF
       ia = indx0(a)
       ib = indx0(b)
       dDistance = -dst0(a,ia)-dst0(b,ib)+dst0(a,ib)+dst0(b,ia)
       IF (dDistance .LT. -small) THEN
          distanc = distanc + dDistance
          indx0(a) = ib
          indx0(b) = ia
       ENDIF
    ENDDO
    Metropolis = (distanc0-distanc)/indim
  END FUNCTION Metropolis

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

  real(8) FUNCTION TMetropolis(ve, max_steps, indim)
    IMPLICIT NONE
    complex(8), intent(in)  :: ve(indim,4)
    INTEGER, intent(in)     :: max_steps
    INTEGER, intent(in)     :: indim
    real(8), PARAMETER :: small = 1e-6
    INTEGER    :: i, j, k, l
    real(8)     :: distanc, distanc0, dDistance
    INTEGER    :: a, b, ia, ib, ii, da, direction, ika, ikb
    INTEGER    :: idum

    IF (indim.LT.2) THEN 
       TMetropolis=0
       return
    ENDIF

    DO l=1,6
       DO i=1,indim
          DO j=1,indim
             dst(j,i,l) = abs(ve(j,stranica(1,l))-ve(i,stranica(2,l)))
          ENDDO
       ENDDO
    ENDDO
    DO i=1,4
       DO j=1,indim
          indx(j,i) = j
       END DO
    END DO
    distanc=0
    DO l=1,6
       DO i=1,indim
          distanc = distanc + dst(i,i,l)
       END DO
    END DO
    distanc0 = distanc

    idum=-1
    DO i=1,max_steps
       ii = int(ran(idum)*4)+1
       a  = int(ran(idum)*indim)+1
       da = int(-7*log(1-ran(idum)))+1
       direction = int(2*ran(idum))
       direction = 2*direction-1
       b = a+direction*da

       IF (b.LT.1) b=1
       IF (b.GT.indim) b=indim
       IF (a.EQ.b) THEN
          IF (a.GT.1) THEN 
             b=a-1
          ELSE
             b=a+1
          ENDIF
       ENDIF
       ia = indx(a,ii)
       ib = indx(b,ii)
       dDistance=0
       DO j=1,3
          l = istranica(j,ii)
          IF (ii==stranica(1,l)) THEN
             k = stranica(2,l)
             ika = indx(a,k)
             ikb = indx(b,k)
             dDistance = dDistance - dst(ia,ika,l) - dst(ib,ikb,l) + dst(ib,ika,l) + dst(ia,ikb,l)
          ELSE 
             k = stranica(1,l)
             ika = indx(a,k)
             ikb = indx(b,k)
             dDistance = dDistance - dst(ika,ia,l) - dst(ikb,ib,l) + dst(ika,ib,l) + dst(ikb,ia,l)
          ENDIF
       END DO
       IF (dDistance<-small) THEN
          distanc = distanc + dDistance
          indx(a,ii) = ib
          indx(b,ii) = ia
       ENDIF
    END DO
    TMetropolis = (distanc0-distanc)/(6*indim)
  END FUNCTION TMetropolis

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
    
  real(8) FUNCTION TDistance(ve, indim)
    IMPLICIT NONE
    complex(8), intent(in) :: ve(indim,4)
    INTEGER, intent(in)    :: indim
    INTEGER    :: i, j, l
    real(8)     :: ds
    ds = 0.0
    DO l=1,6
       DO i=1,indim
          ds = ds + abs(ve(i,stranica(1,l))-ve(i,stranica(2,l)))
       END DO
    END DO
    TDistance = ds/(6*indim)
  END FUNCTION

!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************
!**********************************************************************************

END MODULE
