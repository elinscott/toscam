module hartree_fock

 use specialFunction , only : SlaterF, cmp_all_Gaunt

contains

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

SUBROUTINE Hartree(Sinfty, Edc, U, Jc, occ, bndind, corind, da, ndim, subtractDC)

  !***************************************************************************************!
  !* This subroutine computes Hartree Self-energy for given charge density matrix mocc   *!
  !* given Coulomb repulsion in the form of (U,J). Here we use Slater integrals          *!
  !* F0, F2, F4, F6 rather than 4-dimensional array of Coulomb repulsion. The ration     *!
  !* between F2, F4, and F6 is fixed to default values (see for example Wien code)       *!
  !* so that only U==F0 and J is necessary in the input.                                 *!
  !***************************************************************************************

  IMPLICIT NONE

  complex(8), intent(out) :: Sinfty(ndim,ndim)                 ! Output LDA+U self-energy with double-counting subtracted
  complex(8), intent(out) :: Edc(ndim)                         ! Output LDA+U self-energy with double-counting subtracted
  real(8), intent(in)     :: U, Jc
  complex(8), intent(in)  :: occ(ndim,ndim)                    ! occupations (n)
  INTEGER, intent(in)     :: bndind(ndim,4)                    ! index to (atom,l,m,s)
  INTEGER, intent(in)     :: corind(ndim)                      ! index of correlated orbitals
  INTEGER, intent(in)     :: ndim                              ! dimension of occupation matrix
  real(8), intent(in)     :: da
  INTEGER, intent(in)     :: subtractDC                        ! subtract double-counting or not
  real(8)                 :: gck(0:3,-3:3,-3:3,0:3)            ! Gaunt coefficients are precomputed for speed
  real(8)                 :: Fn(0:3,0:3)
  real(8)                 :: ff, nf_s, nf
  INTEGER                 :: i, j, k, a, b
  INTEGER                 :: a1, l1, m1, s1, a2, l2, m2, s2, a3, l3, m3, s3, a4, l4, m4, s4
  INTEGER                 :: natom, lmax, nspin
  real(8), ALLOCATABLE    :: tnf(:,:,:) !(natom, lmax, nspin)
   
 
  ! Ratio between F2,F4,F6 and J! At the end of the day, we want to use U and J only!
  CALL SlaterF(Fn, U, Jc)
  ! Precomputes gaunt coefficients for speed
  CALL cmp_all_Gaunt(gck)
  ! need natom, lmax, and nspin for allocating tnf
  natom=0; lmax=0; nspin=0
  DO i=1,ndim
     IF (bndind(i,1).GT.natom) natom = bndind(i,1)
     IF (bndind(i,2).GT.lmax)  lmax  = bndind(i,2)
     IF (bndind(i,4).GT.nspin) nspin = bndind(i,4)
  ENDDO

  ALLOCATE(tnf(natom, 0:lmax, nspin))

  Sinfty=0
  DO i=1,ndim
     IF (corind(i).EQ.0) CYCLE ! This index is not correlated
     a1 = bndind(i,1) ! atom
     l1 = bndind(i,2) ! l
     m1 = bndind(i,3) ! m
     s1 = bndind(i,4) ! s
     DO j=1,ndim
        IF (corind(j).EQ.0) CYCLE ! This index is not correlated
        a2 = bndind(j,1)
        l2 = bndind(j,2)
        m2 = bndind(j,3)
        s2 = bndind(j,4)
        IF (a1.NE.a2 .OR. l1.NE.l2) CYCLE
        DO a=1,ndim
           IF (corind(a).EQ.0) CYCLE ! This index is not correlated
           a3 = bndind(a,1)
           l3 = bndind(a,2)
           m3 = bndind(a,3)
           s3 = bndind(a,4)
           IF (a1.NE.a3 .OR. l1.NE.l3) CYCLE
           DO b=1,ndim
              IF (corind(b).EQ.0) CYCLE ! This index is not correlated
              a4 = bndind(b,1)
              l4 = bndind(b,2)
              m4 = bndind(b,3)
              s4 = bndind(b,4)
              IF (a1.NE.a4 .OR. l1.NE.l4) CYCLE
              
              IF (s1.NE.s4 .OR. s2.NE.s3) CYCLE ! Hartree term vanishes
              IF (m1-m4 .NE. m3-m2) CYCLE       ! Hartree term vanishes
              ff = 0
              DO k=0,l1
                 ff = ff + gck(l1,m4,m1,k)*gck(l1,m2,m3,k)*Fn(k,l1)
              ENDDO
              Sinfty(i,b) = Sinfty(i,b) + occ(j,a)*ff ! Hartree
              Sinfty(i,a) = Sinfty(i,a) - occ(j,b)*ff ! Fock
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  
  tnf=0
  DO i=1,ndim
     IF (corind(i).EQ.0) CYCLE ! This index is not correlated
     a1 = bndind(i,1) ! atom
     l1 = bndind(i,2) ! l
     s1 = bndind(i,4) ! s
     tnf(a1,l1,s1) = tnf(a1,l1,s1) + real(occ(i,i))

     !print *, i, real(occ(i,i))
  ENDDO
  
  DO i=1,ndim
     IF (corind(i).EQ.0) CYCLE ! This index is not correlated
     a1 = bndind(i,1) ! atom
     l1 = bndind(i,2) ! l
     s1 = bndind(i,4) ! s

     nf = (tnf(a1,l1,1) + tnf(a1,l1,nspin))*nspin/2.

     nf_s = tnf(a1,l1,s1)
     if (nspin.eq.1) nf_s = 0.5

     Edc(i) = U*(nf-0.5-da) - Jc*(nf_s-0.5)
     
     IF (subtractDC.EQ.1) Sinfty(i,i) = Sinfty(i,i) - Edc(i)
     
  ENDDO
  

  DEALLOCATE(tnf)
END SUBROUTINE HARTREE

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

end module

