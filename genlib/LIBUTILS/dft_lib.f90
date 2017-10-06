module dft_lib

  USE msk
  use genvar
  use symmetrize_tool

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

  SUBROUTINE cmp_olocef(olap, Utk, wk, T2CJ, inf_buf, ndim, nkp)
  IMPLICIT NONE
  !---------- Passed variables ---------!
  complex(8),     intent(in)      :: olap(ndim,ndim,nkp)
  complex(8),     intent(in)      :: Utk(ndim,ndim,nkp)
  real(8),         intent(in)     :: wk(nkp)
  complex(8),     intent(in)      :: T2CJ(ndim,ndim)
  INTEGER,        intent(in)      :: ndim, nkp
  CHARACTER(10000), intent(inout) :: inf_buf

  !---------- Local variables ----------!
  CHARACTER(256)                  :: fh_info
  INTEGER                         :: i, j
  INTEGER                         :: ikp
  complex(8)                      :: ctmp
  complex(8), ALLOCATABLE         :: olminv(:,:)

  ! For inverse and blas multiplication only

  complex(8), ALLOCATABLE         :: oloceff(:,:)
  complex(8), ALLOCATABLE         :: work(:)
  INTEGER                         :: info

  !------- Allocate local memory -------!
  ALLOCATE( olminv(ndim,ndim), work(ndim), oloceff(ndim,ndim))
  oloceff = 0

  !---- Compute the rotated average ----!
  ! The matrix oloceff is the inverse of local overlap. 
  ! To symmetrize oloceff, one needs to be in spheric
  ! harmonics, later we go to cubic harmonics

  DO ikp=1,nkp
     olminv = olap(:,:,ikp)
     CALL CINV (ndim,olminv,ndim,work,info)
     CALL ZTransform(olminv, Utk(:,:,ikp),'C',ndim)
     oloceff = oloceff + olminv*wk(ikp)
  ENDDO

!  Symmetrizes the result

!  CALL symmetrize_spheric(oloceff,ndim)
   CALL symmetrize_cubic(oloceff,T2CJ,ndim)
!  CALL ZTransform(oloceff,T2CJ,'N',ndim) ! To cubic harmonics

  inf_buf = ""
  WRITE(fh_info,*) "CMP_OLOCEF: The effective local overlap matrix"
  inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
  WRITE(fh_info,*) '          |      |       diagonal    off-diagonal |'
  inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
  DO i=1,ndim
     ctmp=0
     DO j=1,ndim
        IF (i.NE.j) ctmp = ctmp + abs(oloceff(i,j))
     ENDDO
     ctmp = ctmp/ndim
     WRITE(fh_info,fmt='(11x,"|",2x,I3,2x,"|",3x,2f12.6,"|",3x,2f12.6,5x,"|")') i, dble(oloceff(i,i)), dble(ctmp)
     inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
  ENDDO
  inf_buf = TRIM(inf_buf)//'\n'

  DEALLOCATE(olminv, work, oloceff)

  RETURN
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

  SUBROUTINE cmp_rho_coefficients(ghh, ghj, gjh, gjj, inpdir, weight, Al, Ar, bndind, lbndind, lndim, ndim, nkp)

  !******************************************************************************!
  !*  To calculate DOS, rho(r) or Gf, we need                                   *!
  !*   1) ghh = A^r.w.A^l                                                       *!
  !*   1) ghj = A^r.w.A^l.Sk^+                                                  *!
  !*   2) gjh = Sk.A^r.w.A^l                                                    *!
  !*   3) gjj = Sk.A^r.w.A^l.Sk^+                                               *!
  !******************************************************************************!

  IMPLICIT NONE

  !---------- Passed variables ---------!
  complex(8), intent(out)     :: ghh(ndim,ndim), ghj(ndim,lndim), gjh(lndim,ndim), gjj(lndim,lndim)
  CHARACTER(1000), intent(in) :: inpdir
  complex(8), intent(in)      :: weight(ndim,nkp)
  complex(8), intent(in)      ::  Al(ndim,ndim,nkp)
  complex(8), intent(in)      :: Ar(ndim,ndim,nkp)
  INTEGER, intent(in)         :: bndind(ndim,4)
  INTEGER, intent(in)         :: lbndind(lndim,4)
  INTEGER, intent(in)         :: lndim, ndim, nkp

  !---------- Local variables ----------!

  INTEGER                     :: ikp                  ! Loop index for k-point loop
  complex(8)                  :: Sk(lndim,ndim)       ! temporary memory for structure constant
  complex(8)                  :: cnhh(ndim,ndim), cnhj(ndim,lndim) ! temporary
  complex(8)                  :: temp(ndim,ndim)
  complex(8)                  :: ena
  INTEGER                     :: iatom, natom, p

  ena = cmplx(1.0,0.0,8)

  natom=0 ! will become max(bndind(#,1))
  DO p=1,ndim
     iatom = bndind(p,1)
     IF (iatom.GT.natom) natom = iatom
  ENDDO

  ghh=0; ghj=0 ; gjh=0; gjj=0

  CALL open_file(inpdir, natom)

  DO ikp=1,nkp
     CALL readS(Sk, ikp, lbndind, bndind, lndim, ndim)
     CALL ZProduct_ADAt(cnhh, Ar(:,:,ikp), weight(:,ikp), Al(:,:,ikp), ndim, ndim, ndim, temp) !cnhh=A^R.w.A^L
     CALL ZProduct_MM(cnhj, cnhh,Sk,'N','C', ndim, ndim, lndim, ndim, ndim, lndim) ! cnhj = cnhh.Sk^+
     CALL ZProduct_sum_MM(gjh, Sk, cnhh, 'N', 'N', lndim, ndim, ndim, ndim, lndim, ndim) ! gjh = gjh + Sk.cnhh
     CALL ZProduct_sum_MM(gjj, Sk, cnhj, 'N', 'N', lndim, ndim, ndim, lndim, lndim, lndim) ! gjj = gjj + Sk.cnhj
     ghh = ghh + cnhh
     ghj = ghj + cnhj
  ENDDO

  CALL close_file

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

SUBROUTINE give_density(nxy, gxy, gyx, nspin, ndim1, ndim2)

  !*********************************************************!
  !* From G computes density by n=(G-G^+)/(2*pi*i)         *!
  !*********************************************************!

  IMPLICIT NONE

  complex(8), intent(out) :: nxy(ndim1,ndim2)
  complex(8), intent(in)  :: gxy(ndim1,ndim2)
  complex(8), intent(in)  :: gyx(ndim2,ndim1)
  INTEGER, intent(in)     :: nspin, ndim1, ndim2
  complex(8)              :: cxy(ndim1,ndim2)
  
  cxy = transpose(conjg(gyx))

  ! additional conjugation for using advaced rather than retarded functions
  nxy = conjg(gxy-cxy)/(2*pi*cmplx(0.,1.,8)*dble(nspin)/2.)

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

SUBROUTINE cmp_occupancy(mocc, nhh,nhj,njh,njj, moo, bndind, is, llmtind, ndim, lndim, lmax, nsort, nspin, illsbmax2, natom)

  !***************************************************************************!
  !*  From the matrix elements of Density in the LMTO base computed above    *!
  !*  (nhh,nhj,njh,njj), computes a full density matrix                      *!
  !***************************************************************************!

  IMPLICIT NONE

  complex(8), intent(out) :: mocc(ndim,ndim)
  complex(8), intent(in)  :: nhh(ndim,ndim), nhj(ndim,lndim), njh(lndim,ndim), njj(lndim,lndim)
  complex(8), intent(in)  :: moo(4,0:lmax,nsort,nspin)
  INTEGER, intent(in)     :: bndind(ndim,4)
  INTEGER, intent(in)     :: is(natom)
  INTEGER, intent(in)     :: llmtind(illsbmax2,natom)
  INTEGER, intent(in)     :: ndim, lndim, lmax, nsort, nspin, illsbmax2, natom
  INTEGER                 :: iatom,  l,  m,  ispin,  lm,  idimh, idimj
  INTEGER                 :: iatom1, l1, m1, ispin1, lm1, jdimh, jdimj
  INTEGER                 :: isort


!  CALL give_density(nhh, ghh, ghh, nspin,  ndim,  ndim)
!  CALL give_density(nhj, ghj, gjh, nspin,  ndim, lndim)
!  CALL give_density(njh, gjh, ghj, nspin, lndim,  ndim)
!  CALL give_density(njj, gjj, gjj, nspin, lndim, lndim)

  mocc=0
  DO idimh=1,ndim
     iatom = bndind(idimh,1)
     l     = bndind(idimh,2)
     m     = bndind(idimh,3)
     ispin = bndind(idimh,4)
     lm = l*(l+1)+m+1
     idimj = llmtind(lm,iatom) + (ispin-1)*lndim/2
     isort = is(iatom)
     DO jdimh=1,ndim
        iatom1 = bndind(jdimh,1)
        l1     = bndind(jdimh,2)
        ispin1 = bndind(jdimh,4)
        IF (.NOT.(iatom.EQ.iatom1 .AND. l.EQ.l1 .AND. ispin.EQ.ispin1)) CYCLE
        m1     = bndind(jdimh,3)
        lm1 = l*(l+1)+m1+1
        jdimj = llmtind(lm1,iatom) + (ispin-1)*lndim/2
        mocc(idimh,jdimh) = conjg(nhh(idimh,jdimh))*moo(1,l,isort,ispin) &
             -conjg(njh(idimj,jdimh))*moo(2,l,isort,ispin)&
             -conjg(nhj(idimh,jdimj))*moo(3,l,isort,ispin)&
             +conjg(njj(idimj,jdimj))*moo(4,l,isort,ispin)
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

SUBROUTINE cmp_rho_valence(rhoval, isym_max, norbs, nhh,nhj,njh,njj, LmSym_a, &
           &         ltpan_a, lmtind_a, llmtind_a, ilsv_s, Cf, gam, Fih, Fij, &
           &         ilst_s, ndim, lndim, nrad_s, nspin, limvmem, lmax, ilsbmax2, illsbmax2, lmax0, lmax1)
  IMPLICIT NONE

  complex(8), intent(out)  :: rhoval(0:nrad_s,isym_max,nspin,norbs)
  complex(8), intent(in)   :: nhh(ndim,ndim), nhj(ndim,lndim), njh(lndim,ndim), njj(lndim,lndim)
  INTEGER, intent(in)      :: LmSym_a(limvmem)
  INTEGER, intent(in)      :: ltpan_a(0:lmax)
  INTEGER, intent(in)      :: lmtind_a(ilsbmax2)
  INTEGER, intent(in)      :: llmtind_a(illsbmax2)
  INTEGER, intent(in)      :: ilsv_s               ! maximum l for this sort
  real(8), intent(in)      :: Cf(0:lmax0-1,lmax1,lmax1), gam(0:ilst_s+1)
  complex(8), intent(in)   :: Fih(0:nrad_s,0:ilst_s,nspin)
  complex(8), intent(in)   :: Fij(0:nrad_s,0:ilst_s,nspin)
  INTEGER, intent(in)      :: isym_max, nspin, norbs
  INTEGER, intent(in)      :: ndim, lndim
  INTEGER, intent(in)      :: ilst_s               ! maximum l for this sort
  INTEGER, intent(in)      :: limvmem, lmax
  INTEGER, intent(in)      :: nrad_s
  INTEGER, intent(in)      :: ilsbmax2, illsbmax2
  INTEGER, intent(in)      :: lmax0, lmax1
  INTEGER                  :: l, l1, l2, isym
  INTEGER                  :: iupdn, iorbs, ispin, ispin1, m, m1, m2, lm, lm1, lm2
  INTEGER                  :: idim1h, idim1j, idimh, idimj, irad
  LOGICAL                  :: BLOG, BLOG1
  complex(8)               :: CGNT

  rhoval = 0

  DO iorbs=1,norbs !!! over spins
     DO iupdn=1,nspin !!! over spins
        IF(iorbs.EQ.1)THEN
           ispin1=iupdn
           ispin =iupdn
        ELSEIF(IORBS.EQ.2)THEN
           ispin1=iupdn
           ispin =nspin+1-iupdn
        ENDIF
        DO l1=0,ilst_s !!! over l1
           BLOG1=.FALSE.
           IF(ltpan_a(l1)==1) BLOG1=.TRUE.
           DO m1=-l1,l1 !!! over m1
              lm1=l1*(l1+1)+m1+1
              IF (BLOG1) idim1h = lmtind_a(lm1)+(ispin1-1)*ndim/2
              idim1j = llmtind_a(lm1)+(ispin1-1)*lndim/2
              DO l=0,ilst_s !!! over l
                 BLOG=.FALSE.
                 IF(ltpan_a(l)==1) BLOG=.TRUE.
                 DO m=-l,l !!! over m
                    lm=l*(l+1)+m+1
                    IF (BLOG) idimh = lmtind_a(lm)+(ispin-1)*ndim/2
                    idimj = llmtind_a(lm)+(ispin-1)*lndim/2
                    DO l2=ABS(l1-l),l1+l,2 !!! over l2 (expansion of ro)
                       IF(l2.LE.ilsv_s)THEN
                          m2=m-m1
                          lm2=l2*(l2+1)+m2+1
                          IF(ABS(m2).LE.l2)THEN !!! over m2 (expansion of ro)
                             isym=LmSym_a(lm2)
                             IF(isym.GT.0) THEN
                                CGNT=(0.D0,1.D0)**(L-L1-L2)*Cf(L2/2,LM,LM1)
                                DO irad=0,nrad_s
                                   IF(BLOG1.AND.BLOG)RHOVAL(irad,isym,iupdn,iorbs)=&
 & RHOVAL(irad,isym,iupdn,iorbs)+(nhh(idim1h,idimh)*CONJG(FIH(irad,l1,ispin1))*FIH(irad,l,ispin))*CGNT
                                   IF(BLOG1)         RHOVAL(irad,isym,iupdn,iorbs)=&
 & RHOVAL(irad,isym,iupdn,iorbs)-(nhj(idim1h,idimj)*CONJG(FIH(irad,l1,ispin1))*FIJ(irad,l,ispin)*gam(L))*CGNT
                                   IF(BLOG)          RHOVAL(irad,isym,iupdn,iorbs)=&
 & RHOVAL(irad,isym,iupdn,iorbs)-(njh(idim1j,idimh)*CONJG(FIJ(irad,l1,ispin1))*FIH(irad,l,ispin)*gam(L1))*CGNT
                                   RHOVAL(irad,isym,iupdn,iorbs)= &
 & RHOVAL(irad,isym,iupdn,iorbs)+(njj(idim1j,idimj)*CONJG(FIJ(irad,l1,ispin1))*FIJ(irad,l,ispin)*gam(L1)*gam(L))*CGNT
                                ENDDO
                             ENDIF
                          ENDIF !!! over m2
                       ENDIF !!!      l2
                    ENDDO !!! over l2 
                 ENDDO !!! over m
              ENDDO !!! over l
           ENDDO !!! over m1
        ENDDO !!! over l1
     ENDDO
  ENDDO

RETURN
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

SUBROUTINE cmp_tgf(ghh, ghj, gjh, gjj, hh, jj, dioloc, norb, gf, ndim, lndim)

  !******************************************************************************!
  !*  To calculate Gloc in DMFT base from coefficients ghh,ghj,gjh,gjj          *!
  !*   1) ghh = A^r.w.A^l                                                       *!
  !*   1) ghj = A^r.w.A^l.Sk^+                                                  *!
  !*   2) gjh = Sk.A^r.w.A^l                                                    *!
  !*   3) gjj = Sk.A^r.w.A^l.Sk^+                                               *!
  !******************************************************************************!

  IMPLICIT NONE

  complex(8), intent(out) :: gf(ndim,ndim)
  complex(8), intent(in)  :: ghh(ndim,ndim), ghj(ndim,lndim), gjh(lndim,ndim), gjj(lndim,lndim)
  complex(8), intent(in)  :: hh(ndim), jj(ndim)
  real(8), intent(in)     :: dioloc(ndim)
  INTEGER, intent(in)     :: norb
  INTEGER, intent(in)     :: ndim, lndim
  INTEGER                 :: ph, pj, qh, qj
  complex(8)              :: gf0

  do ph=1,ndim
     pj = ph + (lndim/2-ndim/2)*int(ph/(norb+1))
     do qh=1,ndim
        qj = qh + (lndim/2-ndim/2)*int(qh/(norb+1))
        gf0 = hh(ph)*ghh(ph,qh)*conjg(hh(qh))-hh(ph)*ghj(ph,qj)*conjg(jj(qh))- &
            & jj(ph)*gjh(pj,qh)*conjg(hh(qh))+jj(ph)*gjj(pj,qj)*conjg(jj(qh))
        gf(ph,qh) = dioloc(ph)*gf0*dioloc(qh)
     enddo
  enddo

RETURN
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


SUBROUTINE cmp_oho(oho, olap, hamf, wk, ndim, nkp)

  !******************************************************************************!
  ! *  Computes  sum_k O^{-1}*Hk*O^{-1} in DMFT base
  !******************************************************************************!

  IMPLICIT NONE

  complex(8), intent(out) :: oho(ndim,ndim)
  complex(8), intent(in)  :: olap(ndim,ndim,nkp)
  complex(8), intent(in)  :: hamf(ndim,ndim,nkp)
  real(8),     intent(in) :: wk(nkp)
  INTEGER, intent(in)     :: ndim, nkp
  complex(8)              :: ham(ndim,ndim), oinv(ndim,ndim), oh(ndim,ndim), work(ndim)
  INTEGER                 :: ikp, info

  oho = 0.0
  DO ikp=1,nkp
    ham = hamf(:,:,ikp)*wk(ikp)
    oinv = olap(:,:,ikp)
    CALL CINV (ndim, oinv, ndim, work, info)
    CALL ZProduct_MM(oh, oinv, ham, 'N','N', ndim, ndim, ndim, ndim, ndim, ndim)       ! oh = O^{-1}*Hk
    CALL ZProduct_sum_MM(oho, oh, oinv, 'N', 'N', ndim, ndim, ndim, ndim, ndim, ndim)  ! oho += oh*oinv 
  ENDDO

RETURN
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

SUBROUTINE cmp_ohol(oho, olap, hamf, Utk, wk, ndim, nkp)

  !******************************************************************************!
  ! *  Computes  sum_k O^{-1}*Hk*O^{-1} in LMTO base
  !******************************************************************************!

  IMPLICIT NONE

  complex(8), intent(out) :: oho(ndim,ndim)
  complex(8), intent(in)  :: olap(ndim,ndim,nkp)
  complex(8), intent(in)  :: hamf(ndim,ndim,nkp)
  complex(8), intent(in)  :: Utk(ndim,ndim,nkp)
  real(8),    intent(in)  :: wk(nkp)
  INTEGER, intent(in)     :: ndim, nkp
  complex(8)              :: ham(ndim,ndim), oinv(ndim,ndim), oh(ndim,ndim), woho(ndim,ndim), work(ndim)
  INTEGER                 :: ikp, info

  oho = 0.0
  DO ikp=1,nkp
     ham = hamf(:,:,ikp)*wk(ikp)
     oinv = olap(:,:,ikp)
     CALL CINV (ndim, oinv, ndim, work, info)
     CALL ZProduct_MM(oh, oinv, ham, 'N','N', ndim, ndim, ndim, ndim, ndim, ndim)       ! oh = O^{-1}*Hk
     CALL ZProduct_MM(woho, oh, oinv, 'N', 'N', ndim, ndim, ndim, ndim, ndim, ndim)     ! oho = oh*oinv 
     CALL ZTransform(woho, Utk(:,:,ikp),'C', ndim)
     oho(:,:) = oho(:,:) + woho(:,:)
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

SUBROUTINE cmp_oso(oso, olap, sinf, wk, ndim, nkp)

  !***************************************************!
  ! *  Computes  sum_k O^{-1}*Sigma*O^{-1}
  !***************************************************!

  IMPLICIT NONE
  complex(8), intent(out) :: oso(ndim,ndim)
  complex(8), intent(in)  :: olap(ndim,ndim,nkp)
  complex(8), intent(in)  :: sinf(ndim,ndim)
  real(8),     intent(in) :: wk(nkp)
  INTEGER, intent(in)     :: ndim, nkp
  complex(8)              :: ham(ndim,ndim), oinv(ndim,ndim), os(ndim,ndim), work(ndim)
  INTEGER                 :: ikp, info

  oso = 0.0
  DO ikp=1,nkp
     ham = sinf*wk(ikp)
     oinv = olap(:,:,ikp)
     CALL CINV (ndim, oinv, ndim, work, info)
     CALL ZProduct_MM(os, oinv, ham, 'N','N', ndim, ndim, ndim, ndim, ndim, ndim)     ! os = O^{-1}*Sigma
     CALL ZProduct_sum_MM(oso, os, oinv, 'N', 'N', ndim, ndim, ndim, ndim, ndim, ndim) ! oso += os*oinv 
  ENDDO
  oso(:,:) = oso(:,:)

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
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************


  SUBROUTINE cmp_osol(oso, olap, sinf, Utk, wk, ndim, nkp)

  !******************************************************************************!
  ! *  Computes  sum_k O^{-1}*Sigma*O^{-1}
  !******************************************************************************!

  IMPLICIT NONE

  complex(8), intent(out) :: oso(ndim,ndim)
  complex(8), intent(in)  :: olap(ndim,ndim,nkp)
  complex(8), intent(in)  :: sinf(ndim,ndim)
  complex(8), intent(in)  :: Utk(ndim,ndim,nkp)
  real(8),    intent(in)  :: wk(nkp)
  INTEGER, intent(in)     :: ndim, nkp
  complex(8)              :: ham(ndim,ndim), oinv(ndim,ndim), os(ndim,ndim), work(ndim)
  INTEGER                 :: ikp, info

  oso = 0.0
  DO ikp=1,nkp
     ham = sinf*wk(ikp) ! Maybe a bug! I guess sinf is still in DMFT!
     oinv = olap(:,:,ikp)
     CALL CINV (ndim, oinv, ndim, work, info)
     CALL ZProduct_MM(os, oinv, ham, 'N', 'N', ndim, ndim, ndim, ndim, ndim, ndim) ! os = O^{-1}*Sigma
     CALL ZProduct_MM(ham, os, oinv, 'N', 'N', ndim, ndim, ndim, ndim, ndim, ndim) ! ham = os*oinv 
     CALL ZTransform(ham, Utk(:,:,ikp),'C', ndim)
     oso(:,:) = oso(:,:) + ham(:,:)
  ENDDO
  oso(:,:) = oso(:,:)
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

SUBROUTINE read_c(filename, n1, n2, n3, dataout, beg_skip, rep_skip)
  IMPLICIT NONE
  CHARACTER*(*), intent(in) :: filename
  INTEGER, intent(in)       :: n1
  INTEGER, intent(in)       :: n2
  INTEGER, intent(in)       :: n3
  INTEGER, intent(in)       :: beg_skip
  INTEGER, intent(in)       :: rep_skip
  complex(8), intent(out)   :: dataout(n1,n2,n3)
  INTEGER                   :: fh_tmp, i, j, k

  fh_tmp = 113
  OPEN(fh_tmp,file=filename ,form='formatted',status='old')
  DO i=1,beg_skip
     READ(fh_tmp,*)
  ENDDO
  IF (rep_skip.NE.0) THEN
     DO k=1,n3
        DO i=1,rep_skip
           READ(fh_tmp,*)
        ENDDO
        READ(fh_tmp,*) ((dataout(i,j,k),i=1,n1),j=1,n2)
     ENDDO
  ELSE
     READ(fh_tmp,*) (((dataout(i,j,k),i=1,n1),j=1,n2),k=1,n3)
  ENDIF
  CLOSE(fh_tmp)
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


SUBROUTINE read_d(filename, n1, n2, n3, dataout, beg_skip, rep_skip)
  IMPLICIT NONE
  CHARACTER*(*), intent(in) :: filename
  INTEGER, intent(in)       :: n1
  INTEGER, intent(in)       :: n2
  INTEGER, intent(in)       :: n3
  INTEGER, intent(in)       :: beg_skip
  INTEGER, intent(in)       :: rep_skip
  real(8),  intent(out)     :: dataout(n1,n2,n3)
  INTEGER                   :: fh_tmp, i, j, k

  fh_tmp = 113
  OPEN(fh_tmp,file=filename ,form='formatted',status='old')
  DO i=1,beg_skip
     READ(fh_tmp,*)
  ENDDO
  IF (rep_skip.NE.0) THEN
     DO k=1,n3
        DO i=1,rep_skip
           READ(fh_tmp,*)
        ENDDO
        READ(fh_tmp,*) ((dataout(i,j,k),i=1,n1),j=1,n2)
     ENDDO
  ELSE
     READ(fh_tmp,*) (((dataout(i,j,k),i=1,n1),j=1,n2),k=1,n3)
  ENDIF
  CLOSE(fh_tmp)

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

SUBROUTINE read_sparse(filename, n1, n2, dataout, beg_skip)
  IMPLICIT NONE

  CHARACTER*(*), intent(in)  :: filename
  INTEGER, intent(in)        :: n1
  INTEGER, intent(in)        :: n2
  INTEGER, intent(in)        :: beg_skip
  complex(8),intent(out)     :: dataout(n1,n2)
  complex(8)                 :: ctmp
  INTEGER                    :: fh_tmp, i, j, i0, j0 

  fh_tmp = 113
  OPEN(fh_tmp,file=filename ,form='formatted',status='old')
  DO i=1,beg_skip
     READ(fh_tmp,*)
  ENDDO
  DO j=1,n2
     READ(fh_tmp,*) j0
     DO i = 1,n1
        READ(fh_tmp,*) i0, ctmp
        dataout(i0,j0) = ctmp
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

SUBROUTINE cmp_mtoverlap(Olap, inpdir, bndind, isrt, nkp, natom, ndim)
  IMPLICIT NONE

  complex(8), intent(out)    :: Olap(ndim,ndim,nkp)
  CHARACTER*1000, intent(in) :: inpdir
  INTEGER, intent(in)        :: bndind(ndim,4)
  INTEGER, intent(in)        :: isrt(natom)
  INTEGER, intent(in)        :: natom, ndim, nkp
  INTEGER                    :: ikp     
  complex(8)                 :: Sk(ndim,ndim), Skp(ndim,ndim), Ojj(ndim,ndim)   ! temporary memory for structure constant
  complex(8)                 :: temp(ndim,ndim)
  INTEGER                    :: iatom, p, q, i, l, isort, ispin, il, is
  INTEGER                    :: l_max, nsort_tmp, nspin_tmp, fh_opar
  complex(8)                 :: moox(ndim,4)
  complex(8), ALLOCATABLE    :: moo(:,:,:,:)

  fh_opar = 998
  OPEN(fh_opar,file=TRIM(inpdir)//"/opar.dat" ,form='formatted',status='old')
  !---------- Read dimension data ----------
  READ(fh_opar,*) l_max, nsort_tmp, nspin_tmp

  ALLOCATE(moo(4,0:l_max,nsort_tmp,nspin_tmp))
  !---------- Reads overlap numbers into moo ----
  DO i=1,4
     READ(fh_opar,*) (((moo(i,l,isort,ispin),l=0,l_max),isort=1,nsort_tmp),ispin=1,nspin_tmp)
     if (i.lt.4) READ(fh_opar,*) ! empty line
  ENDDO
  CLOSE(fh_opar)

  DO p=1,ndim
     iatom = bndind(p,1)
     isort = isrt(iatom)
     il    = bndind(p,2)
     ispin = bndind(p,4)
     DO i=1,4
        moox(p,i) = moo(i,il,isort,ispin)
     ENDDO
  ENDDO

  Olap=0
  ! Reads structure constant
  CALL open_file(inpdir, natom)
  DO ikp=1,nkp
     CALL readS(Sk, ikp, bndind, bndind, ndim, ndim)

     Skp = conjg(transpose(Sk))

     CALL ZProduct_ADAt(Ojj, Skp, moox(:,4), Sk, ndim, ndim, ndim, temp)

     DO p=1,ndim
        Olap(p,p,ikp) = Olap(p,p,ikp) + moox(p,1)

        DO q=1,ndim
           Olap(p,q,ikp) = Olap(p,q,ikp) - Skp(p,q)*moox(q,2)
           Olap(p,q,ikp) = Olap(p,q,ikp) - moox(p,3)*Sk(p,q)
           Olap(p,q,ikp) = Olap(p,q,ikp) + Ojj(p,q)
        ENDDO
     ENDDO
  ENDDO

  CALL close_file

  DEALLOCATE(moo)

END SUBROUTINE cmp_mtoverlap

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

end module
