module symmetrize_tool

use rtm 

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

SUBROUTINE symmetrize_fast(nhh,ndim)
  IMPLICIT NONE
  complex(8), intent(inout) :: nhh(ndim,ndim)
  INTEGER, intent(in)       :: ndim
  complex(8)                :: temp1(ndim,ndim), temp2(ndim,ndim), Uk(ndim,ndim)
  INTEGER                   :: ig
  temp1=0
  DO ig=1,ngroup
     CALL cmp_Ukn(Uk, ig)
     temp2 =  nhh
     CALL ZTransform(temp2,Uk,'N')
     temp1 = temp1 + temp2
  ENDDO
  nhh = temp1/ngroup
END SUBROUTINE

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

SUBROUTINE symmetrize_spheric(nhh,isize)
  IMPLICIT NONE
  INTEGER, intent(in)       :: isize
  complex(8), intent(inout) :: nhh(isize,isize)
  complex(8), ALLOCATABLE   :: A(:,:,:,:)
  complex(8), ALLOCATABLE   :: B(:,:,:,:)
  complex(8), ALLOCATABLE   :: C(:,:,:,:)
  complex(8), ALLOCATABLE   :: temp(:,:)
  complex(8)                :: gs(2,2)
  INTEGER                   :: ig, s, s1, sp, s1p
  INTEGER                   :: i,j

  ALLOCATE(A(norb,norb,nspin,nspin),B(norb,norb,nspin,nspin),C(norb,norb,nspin,nspin), temp(norb,norb))

  IF (.NOT. allocated(Ug) .OR. .NOT. allocated(SJ)) THEN
     WRITE(0,*) 'In module "rtm" memory not allocated. Call rtm.init() first '
     return
  ENDIF

  IF (isize.EQ.ndim .AND. ndim.EQ.2*norb) THEN
     gs(1,1) = 1
     gs(2,2) = 1
     gs(1,2) = 1 !sqrt(2.)
     gs(2,1) = 1 !-sqrt(2.)

     A(:,:,1,1) = nhh(1:norb,     1:norb)
     A(:,:,2,2) = nhh(norb+1:ndim,norb+1:ndim)
     A(:,:,1,2) = nhh(1:norb,     norb+1:ndim)
     A(:,:,2,1) = nhh(norb+1:ndim,1:norb)

     C = 0
     DO ig=1,ngroup
        DO sp=1,2
           DO s1p=1,2
              B(:,:,sp,s1p) = matmul(matmul(transpose(Ug(:,:,ig)),A(:,:,sp,s1p)),conjg(Ug(:,:,ig)))
           ENDDO
        ENDDO
        B(:,:,1,2) = B(:,:,1,2)/gs(1,2)
        B(:,:,2,1) = B(:,:,2,1)/gs(2,1)

        DO s=1,2
           DO s1=1,2
              DO sp=1,2
                 DO s1p=1,2
                    C(:,:,s,s1) = C(:,:,s,s1) + (conjg(SJ(s,sp,ig))*SJ(s1,s1p,ig)*gs(s,s1))*B(:,:,sp,s1p)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     nhh(1:norb,     1:norb)      = C(:,:,1,1)/ngroup
     nhh(norb+1:ndim,norb+1:ndim) = C(:,:,2,2)/ngroup
     nhh(1:norb,     norb+1:ndim) = C(:,:,1,2)/ngroup
     nhh(norb+1:ndim,1:norb)      = C(:,:,2,1)/ngroup
  ELSEIF(isize.EQ.norb) THEN
     A(:,:,1,1) = nhh(1:norb,1:norb)
     C(:,:,1,1) = 0

     DO ig=1,ngroup
        temp(:,:) = matmul(transpose(Ug(:,:,ig)),A(:,:,1,1))
        B(:,:,1,1) = matmul(temp(:,:),conjg(Ug(:,:,ig)))
        C(:,:,1,1) = C(:,:,1,1) + B(:,:,1,1)
     ENDDO
     nhh(1:norb,1:norb) = C(:,:,1,1)/ngroup
  ELSE
     WRITE(0,*) "MOD_SYMT: Error when symmetrizing - size not ndim nor norb"
  ENDIF
  DEALLOCATE(A,B,C,temp)

END SUBROUTINE


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

SUBROUTINE symmetrize_cubic(nhh,T2CJ, ndim)
  IMPLICIT NONE
  complex(8), intent(inout) :: nhh(ndim,ndim)
  complex(8), intent(in)    :: T2CJ(ndim,ndim)
  INTEGER, intent(in)       :: ndim
  ! symmetrization works in spheric harmonics, therefore we transform
  ! from CUB->SPH
  CALL ZTransform(nhh, T2CJ, 'C',ndim)
  ! symmetrize in spheric
  CALL symmetrize_spheric(nhh,ndim)
  ! and back from SPH->CUB
  CALL ZTransform(nhh, T2CJ, 'N',ndim)
END SUBROUTINE

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


SUBROUTINE symmetrize_jj(njj, lndim)
  IMPLICIT NONE
  complex(8), intent(inout) :: njj(lndim,lndim)
  INTEGER, intent(in)       :: lndim
  complex(8), ALLOCATABLE   :: A(:,:,:,:)
  complex(8), ALLOCATABLE   :: B(:,:,:,:)
  complex(8), ALLOCATABLE   :: C(:,:,:,:)
  complex(8)                :: gs(2,2)
  INTEGER                   :: ig, s, s1, sp, s1p

  IF (.NOT. allocated(Ugll) .OR. .NOT. allocated(SJ)) THEN
     WRITE(0,*) 'In module "rtm" memory not allocated. Call rtm.init() first '
     return
  ENDIF
  ALLOCATE(A(lnorb,lnorb,nspin,nspin),B(lnorb,lnorb,nspin,nspin),C(lnorb,lnorb,nspin,nspin))
  IF (ndim.EQ.2*norb) THEN
     gs(1,1) = 1
     gs(2,2) = 1
     gs(1,2) = 1!sqrt(2.)
     gs(2,1) = 1!-sqrt(2.)

     A(:,:,1,1) = njj(1:lnorb,     1:lnorb)
     A(:,:,2,2) = njj(lnorb+1:lndim,lnorb+1:lndim)
     A(:,:,1,2) = njj(1:lnorb,     lnorb+1:lndim)
     A(:,:,2,1) = njj(lnorb+1:lndim,1:lnorb)

     C = 0
     DO ig=1,ngroup
        DO sp=1,2
           DO s1p=1,2
              B(:,:,sp,s1p) = matmul(matmul(transpose(Ugll(:,:,ig)),A(:,:,sp,s1p)),conjg(Ugll(:,:,ig)))
           ENDDO
        ENDDO
        B(:,:,1,2) = B(:,:,1,2)/gs(1,2)
        B(:,:,2,1) = B(:,:,2,1)/gs(2,1)

        DO s=1,2
           DO s1=1,2
              DO sp=1,2
                 DO s1p=1,2
                    C(:,:,s,s1) = C(:,:,s,s1) + (conjg(SJ(s,sp,ig))*SJ(s1,s1p,ig)*gs(s,s1))*B(:,:,sp,s1p)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     njj(1:lnorb,     1:lnorb)        = C(:,:,1,1)/ngroup
     njj(lnorb+1:lndim,lnorb+1:lndim) = C(:,:,2,2)/ngroup
     njj(1:lnorb,     lnorb+1:lndim)  = C(:,:,1,2)/ngroup
     njj(lnorb+1:lndim,1:lnorb)       = C(:,:,2,1)/ngroup
  ELSEIF (ndim.EQ.norb) THEN
     A(:,:,1,1) = njj(1:lnorb,1:lnorb)
     C(:,:,1,1) = 0
     DO ig=1,ngroup
        B(:,:,1,1) = matmul(matmul(transpose(Ugll(:,:,ig)),A(:,:,1,1)),conjg(Ugll(:,:,ig)))
        C(:,:,1,1) = C(:,:,1,1) + B(:,:,1,1)
     ENDDO
     njj(1:lnorb,1:lnorb) = C(:,:,1,1)/ngroup
  ELSE
     WRITE(0,*) "MOD_SYMT: Error when symmetrizing - size not ndim nor norb"
  ENDIF
  DEALLOCATE(A,B,C)
END SUBROUTINE

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

SUBROUTINE symmetrize_hj(nhj, ndim, lndim)
  IMPLICIT NONE
  complex(8), intent(inout) :: nhj(ndim,lndim)
  INTEGER, intent(in)       :: ndim, lndim
  complex(8), ALLOCATABLE   :: A(:,:,:,:)
  complex(8), ALLOCATABLE   :: B(:,:,:,:)
  complex(8), ALLOCATABLE   :: C(:,:,:,:)
  complex(8)                :: gs(2,2)
  INTEGER                   :: ig, s, s1, sp, s1p

  IF (.NOT. allocated(Ug) .OR. .NOT. allocated(SJ)) THEN
     WRITE(0,*) 'In module "rtm" memory not allocated. Call rtm.init() first '
     return
  ENDIF

  ALLOCATE(A(norb,lnorb,nspin,nspin),B(norb,lnorb,nspin,nspin),C(norb,lnorb,nspin,nspin))

  IF (ndim.EQ.2*norb) THEN
     gs(1,1) = 1
     gs(2,2) = 1
     gs(1,2) = 1!sqrt(2.)
     gs(2,1) = 1!-sqrt(2.)

     A(:,:,1,1) = nhj(1:norb,     1:lnorb)
     A(:,:,2,2) = nhj(norb+1:ndim,lnorb+1:lndim)
     A(:,:,1,2) = nhj(1:norb,     lnorb+1:lndim)
     A(:,:,2,1) = nhj(norb+1:ndim,1:lnorb)

     C = 0
     DO ig=1,ngroup
        DO sp=1,2
           DO s1p=1,2
              B(:,:,sp,s1p) = matmul(matmul(transpose(Ug(:,:,ig)),A(:,:,sp,s1p)),conjg(Ugll(:,:,ig)))
           ENDDO
        ENDDO
        B(:,:,1,2) = B(:,:,1,2)/gs(1,2)
        B(:,:,2,1) = B(:,:,2,1)/gs(2,1)

        DO s=1,2
           DO s1=1,2
              DO sp=1,2
                 DO s1p=1,2
                    C(:,:,s,s1) = C(:,:,s,s1) + (conjg(SJ(s,sp,ig))*SJ(s1,s1p,ig)*gs(s,s1))*B(:,:,sp,s1p)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     nhj(1:norb,     1:lnorb)        = C(:,:,1,1)/ngroup
     nhj(norb+1:ndim,lnorb+1:lndim) = C(:,:,2,2)/ngroup
     nhj(1:norb,     lnorb+1:lndim)  = C(:,:,1,2)/ngroup
     nhj(norb+1:ndim,1:lnorb)       = C(:,:,2,1)/ngroup
  ELSEIF (ndim.EQ.norb) THEN
     A(:,:,1,1) = nhj(1:norb,1:lnorb)
     C(:,:,1,1) = 0
     DO ig=1,ngroup
        B(:,:,1,1) = matmul(matmul(transpose(Ug(:,:,ig)),A(:,:,1,1)),conjg(Ugll(:,:,ig)))
        C(:,:,1,1) = C(:,:,1,1) + B(:,:,1,1)
     ENDDO
     nhj(1:norb,1:lnorb) = C(:,:,1,1)/ngroup
  ELSE
     WRITE(0,*) "MOD_SYMT: Error when symmetrizing - size not ndim nor norb"
  ENDIF
  DEALLOCATE(A,B,C)

END SUBROUTINE

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

SUBROUTINE symmetrize_jh(njh, lndim, ndim)
  IMPLICIT NONE
  complex(8), intent(inout) :: njh(lndim,ndim)
  INTEGER, intent(in)       :: lndim, ndim
  complex(8), ALLOCATABLE   :: A(:,:,:,:)
  complex(8), ALLOCATABLE   :: B(:,:,:,:)
  complex(8), ALLOCATABLE   :: C(:,:,:,:)
  complex(8)                :: gs(2,2)
  INTEGER                   :: ig, s, s1, sp, s1p

  IF (.NOT. allocated(Ug) .OR. .NOT. allocated(SJ)) THEN
     WRITE(0,*) 'In module "rtm" memory not allocated. Call rtm.init() first '
     return
  ENDIF

  ALLOCATE(A(lnorb,norb,nspin,nspin),B(lnorb,norb,nspin,nspin),C(lnorb,norb,nspin,nspin))

  IF (ndim.EQ.2*norb) THEN
     gs(1,1) = 1
     gs(2,2) = 1
     gs(1,2) = 1!sqrt(2.)
     gs(2,1) = 1!-sqrt(2.)

     A(:,:,1,1) = njh(1:lnorb,     1:norb)
     A(:,:,2,2) = njh(lnorb+1:lndim,norb+1:ndim)
     A(:,:,1,2) = njh(1:lnorb,      norb+1:ndim)
     A(:,:,2,1) = njh(lnorb+1:lndim,1:norb)

     C = 0
     DO ig=1,ngroup
        DO sp=1,2
           DO s1p=1,2
              B(:,:,sp,s1p) = matmul(matmul(transpose(Ugll(:,:,ig)),A(:,:,sp,s1p)),conjg(Ug(:,:,ig)))
           ENDDO
        ENDDO
        B(:,:,1,2) = B(:,:,1,2)/gs(1,2)
        B(:,:,2,1) = B(:,:,2,1)/gs(2,1)

        DO s=1,2
           DO s1=1,2
              DO sp=1,2
                 DO s1p=1,2
                    C(:,:,s,s1) = C(:,:,s,s1) + (conjg(SJ(s,sp,ig))*SJ(s1,s1p,ig)*gs(s,s1))*B(:,:,sp,s1p)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     njh(1:lnorb,      1:norb)       = C(:,:,1,1)/ngroup
     njh(lnorb+1:lndim,norb+1:ndim)  = C(:,:,2,2)/ngroup
     njh(1:lnorb,      norb+1:ndim)  = C(:,:,1,2)/ngroup
     njh(lnorb+1:lndim,1:norb)       = C(:,:,2,1)/ngroup
  ELSEIF (ndim.EQ.norb) THEN
     A(:,:,1,1) = njh(1:lnorb,1:norb)
     C(:,:,1,1) = 0
     DO ig=1,ngroup
        B(:,:,1,1) = matmul(matmul(transpose(Ugll(:,:,ig)),A(:,:,1,1)),conjg(Ug(:,:,ig)))
        C(:,:,1,1) = C(:,:,1,1) + B(:,:,1,1)
     ENDDO
     njh(1:lnorb,1:norb) = C(:,:,1,1)/ngroup
  ELSE
     WRITE(0,*) "MOD_SYMT: Error when symmetrizing - size not ndim nor norb"
  ENDIF
  DEALLOCATE(A,B,C)
END SUBROUTINE 

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

end module
