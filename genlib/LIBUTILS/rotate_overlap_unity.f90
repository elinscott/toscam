module rotate_overlap_unity

  USE msk

contains

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   !------------------------------------------------------------------------------!
   ! This file  contains a transformation of base such that the overlap           !
   ! in the new base is close to unity == the new base is alomst orthogonal       !
   ! This is convenient for DMFT applications since the impurity levels are than  !
   ! less sensitive to the Sigma_{infinity} and at the same time the impurity     !
   ! occupation is equal to the lattice occupation.                               !
   !------------------------------------------------------------------------------!

  SUBROUTINE get_overlap_numbers(hh,jj,moo,inpdir,bndind,is,ltsort,lmax,nspin,inf_buf,nsort,ndim,natom,lmaxu)
  IMPLICIT NONE
  complex(8), intent(out)         :: hh(ndim)
  complex(8), intent(out)         :: jj(ndim)
  complex(8), intent(out)         :: moo(4,0:lmax,nsort,nspin)
  CHARACTER*1000, intent(in)      :: inpdir
  INTEGER, intent(in)             :: bndind(ndim,4)
  INTEGER, intent(in)             :: is(natom)
  INTEGER, intent(in)             :: ltsort(0:lmaxu,nsort)
  INTEGER, intent(in)             :: lmax, nspin
  CHARACTER(10000), intent(inout) :: inf_buf
  integer                         :: ndim
  integer                         :: natom
  integer                         :: nsort
  integer                         :: lmaxu
  INTEGER                         :: fh_opar    ! Temporary file handle
  INTEGER                         :: nsort_tmp  ! Number of sorts from overlap numbers
  INTEGER                         :: nspin_tmp  ! Number of spins from overlap numbers
  INTEGER                         :: i, j, il, jl, ispin, jspin, isort, iatom, jatom
  INTEGER                         :: l, l_max
  CHARACTER(256)                  :: fh_info
  complex(8), ALLOCATABLE         :: hhjj0(:,:,:,:)
  real(8)                         :: mx(4), mg(4), mf, mw(22), meps, macc
  INTEGER                         :: mifun, miter, mflag
  complex(8)                      :: co1, co2, co3, co4

  fh_opar = 998
  OPEN(fh_opar,file=TRIM(inpdir)//"/opar.dat" ,form='formatted',status='old')

  !---------- Read dimension data ----------
  READ(fh_opar,*) l_max, nsort_tmp, nspin_tmp

  IF (nsort_tmp.NE.nsort) WRITE(0,*) "Fort-get_overlap: Error in reading opar.dat, nsort not consistent"
  IF (nspin_tmp.NE.nspin) WRITE(0,*) "Fort-get_overlap:: Error in reading opar.dat, nspin not consistent"
  IF (l_max.NE.lmax) WRITE(0,*) "Fort-get_overlap: Error in reading opar.dat, lmax not consistent"

  !---------- Reads overlap numbers into moo ----

  DO i=1,4
     READ(fh_opar,*) (((moo(i,l,isort,ispin),l=0,l_max),isort=1,nsort),ispin=1,nspin_tmp)
     if (i.lt.4) READ(fh_opar,*) ! empty line
!     WRITE(fh_info,*) (((moo(i,l,isort,ispin),l=0,l_max),isort=1,nsort),ispin=1,nspin_tmp)
  ENDDO
  CLOSE(fh_opar)

  ALLOCATE(hhjj0(0:l_max,nsort_tmp,nspin_tmp,2))    
  hhjj0 = 0

! Find the best coefficient h and j such that 
! o(1)=h^{*}h  ==  OHH
! o(2)=j^{*}h  ==  OJH
! o(3)=h^{*}j  ==  OHJ
! o(4)=j^{*}j  ==  OJJ
! O_k = (j . S_k - h)^{+} . (j . S_k - h) = o(1) - S_k^{+} o(2) - o(3) S_k + S_k^{+} o(4) S_k

  macc=10.D-20
  meps=.00001
  
  write(fh_info,*) 'Fort-get_overlap: Minimization for overlap coefficients gives the following results'
  inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)

  DO l=0,l_max
     DO isort=1,nsort
        if (.NOT. ltsort(l,isort)==1) CYCLE
        DO ispin=1,nspin
           ! Set initial estimates to the minimizer
           mx(1)=sqrt(abs(moo(1,l,isort,ispin)))
           mx(2)=0.
           mx(3)=sqrt(abs(moo(4,l,isort,ispin)))
           if (dble(moo(2,l,isort,ispin)).lt.0) mx(3) = -mx(3) 
           mx(4)=0.

           CALL CONMIN(4,mx,mf,mg,mifun,miter,meps,mflag,300,mw,0,22,6,macc,0,moo(:,l,isort,ispin))

           if(mflag.GT.0) THEN
             write(fh_info,*) 'Fort-get_overlap: Minimizer did not succeeded at ', l, isort, ispin
             inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
             write(fh_info,*) 'Fort-get_overlap: return code is', mflag, mf, miter , mifun
             inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
           endif

           write (fh_info,'(3x,A,i1,1x,A,i1,1x,A,i1,1x,A,f14.4,2x,A,f14.4,2x,A,2x,f14.4)') &
           & 'l=',l,'sort=',isort,'spin=',ispin,'h=',mx(1),'j=',mx(3), ' difference:', mf
           inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
           co1 = mx(1)**2+mx(2)**2
           co2 = cmplx(mx(1), mx(2))*cmplx(mx(3),-mx(4),8)
           co3 = cmplx(mx(1),-mx(2))*cmplx(mx(3), mx(4),8)
           co4 = mx(3)**2+mx(4)**2
           write (fh_info, '(3x,A,1x,f14.4,2x,f14.4,2x,f14.4,2x,f14.4)') 'minimized :', dble(co1), dble(co2), dble(co3), dble(co4)
           inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
           write (fh_info, '(3x,A,1x,f14.4,2x,f14.4,2x,f14.4,2x,f14.4)') 'original  :', dble(moo(1,l,isort,ispin)), &
                dble(moo(2,l,isort,ispin)), dble(moo(3,l,isort,ispin)), dble(moo(4,l,isort,ispin))
           inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)

           hhjj0(l,isort,ispin,1) = cmplx(mx(1), mx(2),8);
           hhjj0(l,isort,ispin,2) = cmplx(mx(3), mx(4),8);
        ENDDO
     ENDDO
  ENDDO


  !------ Store overlap numbers into normal vector: hh(1:ndim) and jj(1:ndim)
  !    noo=0

  DO i=1,ndim
     iatom = bndind(i,1)
     il    = bndind(i,2)
     ispin = bndind(i,4)
     isort = is(iatom)
     ! make a vector of hh,jj values
     hh(i) = hhjj0(il,isort,ispin,1)
     jj(i) = hhjj0(il,isort,ispin,2)
     !       DO j=1,ndim
     !          jatom = bndind(j,1)
     !          jl    = bndind(j,2)
     !          jspin = bndind(j,4)
     !         ! makes a matrix of overlap numbers. Somethimes it is convenient to have matrix instead of vector
     !          IF (iatom.EQ.jatom .AND. il.EQ.jl .AND. ispin.EQ.jspin) THEN
     !             noo(1,i,j) = moo(1,il,isort,ispin)
     !             noo(2,i,j) = moo(2,il,isort,ispin)
     !             noo(3,i,j) = moo(3,il,isort,ispin)
     !             noo(4,i,j) = moo(4,il,isort,ispin)
     !          ENDIF
     !       ENDDO
  ENDDO

  DEALLOCATE(hhjj0)

return
END SUBROUTINE

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************


SUBROUTINE cmp_Uloc(Uloc, symmetrize_spheric, Utk, transform, olap, wk, ndim, nkp)

  ! Calculates transformation matrix Uloc, which is equal to
  ! srqt(Oloc) = sqrt((\sum_k Ok^{-1})^{-1})

  IMPLICIT NONE

  complex(8), intent(out)   :: Uloc(ndim,ndim)
  complex(8), intent(in)    :: Utk(ndim,ndim,nkp)
  LOGICAL,    intent(in)    :: transform
  complex(8), intent(in)    :: olap(ndim,ndim,nkp)
  real(8),     intent(in)   :: wk(nkp)
  INTEGER,    intent(in)    :: ndim, nkp
  EXTERNAL                  :: symmetrize_spheric
  complex(8)                :: dummy_argm(ndim,ndim)
  complex(8), ALLOCATABLE   :: olap1(:,:), oloc(:,:)
  INTEGER                   :: LWORK
  complex(8),ALLOCATABLE    :: WORK(:), tmp(:,:)
  real(8),ALLOCATABLE       :: RWORK(:), W(:)
  INTEGER                   :: ikp, INFO, i, j

  CALL symmetrize_spheric(dummy_argm,ndim)
 
  LWORK = 3*ndim
  ALLOCATE(WORK(LWORK))
  ALLOCATE(olap1(ndim,ndim), oloc(ndim,ndim))

  oloc=0
  DO ikp=1,nkp
     olap1 = olap(:,:,ikp)
! Calculates Ok^-1
     CALL CINV(ndim,olap1,ndim,WORK,INFO)
! and transforms Ok^-1 to new base via Uk Ok^-1 Uk^+
     IF (transform) CALL ZTransform(olap1, Utk(:,:,ikp), 'C', ndim)
! sums up to get local O_loc
     oloc = oloc + olap1*wk(ikp)
  ENDDO

! To simulate summation over all k-points, symmetrize the quantity
  CALL symmetrize_spheric(oloc,ndim)

! O_loc^{-1} = (\sum_k O_k^{-1})
  olap1 = oloc
  
  ALLOCATE(W(ndim),RWORK(LWORK))

! eigensystem is calculated for O_loc^{-1}

  CALL zheev('V','U',ndim,olap1,ndim,W,WORK,LWORK,RWORK,INFO)
  IF(INFO.NE.0) THEN
     WRITE(0,*) "LAPACK_ZHEEV Error solving eigenvalue problem"
     return 
  ENDIF

  Uloc=0
  DO i=1,ndim
     Uloc(i,i) = 1/sqrt(w(i))
  ENDDO

  ALLOCATE(tmp(ndim,ndim))

  tmp = matmul(olap1,Uloc)
  Uloc = matmul(tmp,transpose(conjg(olap1)))

  DEALLOCATE(olap1, oloc)
  DEALLOCATE(tmp)
  DEALLOCATE(WORK,RWORK,W)

END SUBROUTINE

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

SUBROUTINE cmp_transformation1(symmetrize_spheric_, inpdir, olap, wk, is, ilsb, lmtind, ltpan, &
                             & bndind, norb, nspin, lmax, hh, jj, inf_buf, &
                             &   ndim, nkp, nsort, natom, ilsbmax2, lmaxu, Utk0, dioloc)
  IMPLICIT NONE

  complex(8), intent(out)       :: Utk0(ndim,ndim,nkp)
  real(8),     intent(out)      :: dioloc(ndim)
  complex(8), intent(in)        :: olap(ndim,ndim,nkp)
  real(8),     intent(in)       :: wk(nkp)
  INTEGER,    intent(in)        :: is(natom)
  INTEGER,    intent(in)        :: ilsb(nsort)
  INTEGER,    intent(in)        :: lmtind(ilsbmax2,natom)
  INTEGER,    intent(in)        :: ltpan(0:lmaxu,natom)
  INTEGER,    intent(in)        :: bndind(ndim,4)
  complex(8), intent(in)        :: hh(ndim), jj(ndim)
  CHARACTER*1000,intent(in)     :: inpdir
  INTEGER,    intent(in)        :: norb, nspin, lmax, ndim, nkp, nsort, natom, ilsbmax2, lmaxu 
  CHARACTER*1000, intent(inout) :: inf_buf
  EXTERNAL                      :: symmetrize_spheric_
  complex(8)                    :: dummy_argm(ndim,ndim)
  INTEGER                       :: i, j, ikp, s, iatom, isort, l, m, lm, orb
  INTEGER                       :: numsrt(nsort)
  complex(8)                    :: olap1(norb,norb), oloc(norb,norb)
  complex(8)                    :: Uloc(ndim,ndim), Sk(ndim,ndim)
  real(8)                       :: olav(nsort,0:lmaxu)
  real(8)                       :: rtmp
  CHARACTER*256                 :: fh_info

  CALL symmetrize_spheric_(dummy_argm,ndim)

  CALL open_file(inpdir, natom) ! small version


  !---- Create transformation matrix Ut_k = H - J S_k
  DO ikp=1,nkp
     CALL readS(Sk, ikp, bndind, bndind, ndim, ndim)
     DO i=1,ndim
        Utk0(i,:,ikp) = Sk(i,:)*(-jj(i));
        Utk0(i,i,ikp) = Utk0(i,i,ikp) + hh(i)
     ENDDO
  ENDDO
  CALL close_file


  !---------- Initialize the local overlap matrix ----------
  WRITE(fh_info,fmt='(1x,A)') "Fort-cmp_transform: Computes transformation Utk, &
               & which transforms to alomst orhogonal base"
  inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)


  !----- Improves tranformation, O_local should be closer to unity
! First calculates the transformation that makes O_local equal unity
! Uloc = sqrt((\sum_k Ok^{-1})^{-1})


  CALL cmp_Uloc(Uloc, symmetrize_spheric_, Utk0, .TRUE., olap, wk, ndim, nkp)

! Uloc, however, can not be used as it is, since
! it does not preserve L! For that, one needs to take an identity in each
! subspace of definite l.

  numsrt=0
  DO iatom = 1,natom
     isort = is(iatom)
     numsrt(isort) = numsrt(isort)+1
  ENDDO


  olav=0


  !--- Diagonal average for each sort and angular momentum ---
  DO iatom = 1,natom
     isort = is(iatom)
     DO l = 0,ilsb(isort)
        if (.NOT.ltpan(l,iatom)==1) cycle
        rtmp = 0.
        DO m = -l,l
           lm = l*(l+1)+m+1
           orb = lmtind(lm,iatom)
           rtmp = rtmp+dble(Uloc(orb,orb))
        ENDDO
        olav(isort,l) = olav(isort,l)+rtmp/(2.*l+1)/numsrt(isort)
     ENDDO
  ENDDO


  !--------- Writing out olav -----------------------
  WRITE(fh_info,fmt='(1x,A)') "Fort-cmp_transform: Additional local transformation is applied to get orthogonal DMFT impurity problem"
  inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
  WRITE(fh_info,fmt='(3x,A)') "Due to overlap in interstitial region, this transformation is sometimes relatively far from unity"
  inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
  WRITE(fh_info,fmt='(3x,A,1x)') 'dioloc: ' 
  inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)//'\n'
  DO isort=1,nsort
     WRITE(fh_info,fmt='("sort=",I2,":",1x)') isort
     inf_buf = TRIM(inf_buf)//"  "//TRIM(fh_info)//"  "
     DO l=0,lmax
        WRITE(fh_info,fmt='("l=",I1,":",1x)') l
        inf_buf = TRIM(inf_buf)//"  "//TRIM(fh_info)
        WRITE(fh_info,fmt='(f15.8)') olav(isort,l)
        inf_buf = TRIM(inf_buf)//" "//TRIM(ADJUSTL(fh_info))
     ENDDO
     inf_buf = TRIM(inf_buf)//"\n"
  ENDDO


  !---------- Make diagonal average vector ----------
  dioloc=0
  DO iatom = 1,natom
     isort = is(iatom)
     DO l = 0,ilsb(isort)
        if (.NOT.ltpan(l,iatom)==1) cycle
        DO m = -l,l
           DO s=1,nspin
              lm = l*(l+1)+m+1 
              orb = lmtind(lm,iatom)+ (s-1)*norb
              dioloc(orb) = olav(isort,l)
           ENDDO
        ENDDO
     ENDDO
  ENDDO


  !------- Finally update transformation -------------
  DO ikp=1,nkp
     DO i=1,ndim
        Utk0(i,:,ikp) = dioloc(i)*Utk0(i,:,ikp)
     ENDDO
  ENDDO

END SUBROUTINE 


!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

SUBROUTINE transform(ham, Utk0, ndim, nkp)
  IMPLICIT NONE
  complex(8), intent(inout)     :: ham(ndim,ndim,nkp)
  complex(8), intent(in)        :: Utk0(ndim,ndim,nkp)
  INTEGER, intent(in)           :: ndim, nkp
  INTEGER                       :: i, j, ikp
  complex(8), ALLOCATABLE       :: Stk(:,:)
  INTEGER                       :: INFO
  complex(8), ALLOCATABLE       :: WORK(:)

  ALLOCATE(WORK(ndim) )    ! For inversion
  ALLOCATE(Stk(ndim,ndim)) ! St_k=Utk^-1

  DO ikp=1,nkp
! Create transformation matrix St_k = Utk^{-1}, where Utk=H-J*Sk
     Stk(:,:) = Utk0(:,:,ikp)
     CALL CINV(ndim,Stk,ndim,WORK,INFO)
! Transform Hamiltonian and overlap into new (almost) orthonormal base
! H_new = St_k^{+} H_old St_k
     CALL ZTransform(ham(:,:,ikp),Stk,'N',ndim)
  ENDDO
  DEALLOCATE(WORK)
  DEALLOCATE(Stk)

END SUBROUTINE 


!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

SUBROUTINE transform2(ham1, ham2, Utk0, ndim, nkp)
  IMPLICIT NONE
  complex(8), intent(inout)     :: ham1(ndim,ndim, nkp)
  complex(8), intent(inout)     :: ham2(ndim,ndim, nkp)
  complex(8), intent(in)        :: Utk0(ndim,ndim, nkp)
  INTEGER, intent(in)           :: ndim, nkp
  INTEGER                       :: i, j, ikp
  complex(8), ALLOCATABLE       :: Stk(:,:)
  INTEGER                       :: INFO
  complex(8), ALLOCATABLE       :: WORK(:)

  ALLOCATE(WORK(ndim) )    ! For inversion
  ALLOCATE(Stk(ndim,ndim)) ! St_k=Utk^-1

  DO ikp=1,nkp
! Create transformation matrix St_k = Utk^{-1}, where Utk=H-J*Sk
     Stk(:,:) = Utk0(:,:,ikp)
     CALL CINV(ndim,Stk,ndim,WORK,INFO)
! Transform Hamiltonian and overlap into new (almost) orthonormal base
! H_new = St_k^{+} H_old St_k
     CALL ZTransform(ham1(:,:,ikp),Stk,'N',ndim)
     CALL ZTransform(ham2(:,:,ikp),Stk,'N',ndim)
  ENDDO

  DEALLOCATE(WORK)
  DEALLOCATE(Stk)

END SUBROUTINE

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

SUBROUTINE transform3(ve1,ve2,ve3,Utk0,ndim,nkp)
  IMPLICIT NONE
  complex(8), intent(inout)     :: ve1(ndim,ndim, nkp)
  complex(8), intent(inout)     :: ve2(ndim,ndim, nkp)
  complex(8), intent(inout)     :: ve3(ndim,ndim, nkp)
  complex(8), intent(in)        :: Utk0(ndim,ndim, nkp)
  INTEGER, intent(in)           :: ndim, nkp
  INTEGER                       :: i, j, ikp
  complex(8), ALLOCATABLE       :: Stk(:,:)
  INTEGER                       :: INFO
  complex(8), ALLOCATABLE       :: WORK(:)

  ALLOCATE(WORK(ndim) )    ! For inversion
  ALLOCATE(Stk(ndim,ndim)) ! St_k=Utk^-1


  DO ikp=1,nkp
!  Create transformation matrix St_k = Utk^{-1}, where Utk=H-J*Sk
     Stk(:,:) = Utk0(:,:,ikp)
     CALL CINV(ndim,Stk,ndim,WORK,INFO)
!  Transform Hamiltonian and overlap into new (almost) orthonormal base
!  H_new = St_k^{+} H_old St_k
     CALL ZTransform(ve1(:,:,ikp),Stk,'N',ndim)
     CALL ZTransform(ve2(:,:,ikp),Stk,'N',ndim)
     CALL ZTransform(ve3(:,:,ikp),Stk,'N',ndim)
  ENDDO

  DEALLOCATE(WORK)
  DEALLOCATE(Stk)

END SUBROUTINE 

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

end module
