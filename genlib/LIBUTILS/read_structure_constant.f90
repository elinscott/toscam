MODULE msk

  ! This module read the structure constant, eack k-point separately

  PRIVATE 

  INTEGER                 :: natom
  INTEGER                 :: fh_Sk                  ! filehandle
  complex(8), ALLOCATABLE :: Sk_temp(:,:,:,:)       ! Temporary array for reading
  INTEGER                 :: ilm1_max, ilm2_max     ! Largest lm index - temporary

CONTAINS

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE open_file(inpdir, natom_)
    IMPLICIT NONE
    !---------- Passed variables ----------
    CHARACTER*1000, intent(in) :: inpdir
    INTEGER, intent(in)        :: natom_
    !---------- Local variables ----------
    INTEGER :: natom_tmp              ! Assumed number of atoms
    INTEGER :: iatom                  ! Atom index
    real(8)  :: latcon_tmp             ! Lattice constant
    
    fh_Sk = 996
    OPEN(fh_Sk,file=TRIM(inpdir)//"/Sk.dat" ,form='formatted',status='old')
    !---------- Read dimension data ----------
    READ(fh_Sk,*) ilm1_max, ilm2_max, natom_tmp, latcon_tmp
    natom = natom_
    IF (natom_tmp .NE. natom) WRITE(0,*) "Troubles reading Sk. natom!=natom_tmp"
    ! Temporary variable for reading
    ALLOCATE(Sk_temp(ilm1_max,natom,ilm2_max,natom))
  END SUBROUTINE open_file

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE readS(Sk, ikp, lbndind, bndind, lndim, ndim)
    IMPLICIT NONE
    !---------- Passed variables ----------
    INTEGER, intent(in) :: ikp
    INTEGER, intent(in) :: ndim, lndim
    INTEGER, intent(in) :: lbndind(lndim,4), bndind(ndim,4)
    complex(8), intent(out) :: Sk(lndim,ndim)
    !---------- Local variables ----------
    INTEGER :: iatom,iatom1,iatom2    ! Atom index
    INTEGER :: itmp                   ! Integer dummy variable
    INTEGER :: p, q, s, s1            ! Orbital index
    INTEGER :: l,m,lm,l1,m1,lm1,lm2   ! Angular momentum indices
    READ(fh_Sk,*) itmp
    IF(itmp.NE.ikp) WRITE(0,*) "mod_Sk.f90: Error in reading structure constant Sk.dat"
    ! Actual reading in tmp array
    READ(fh_Sk,*) ((((Sk_temp(lm1,iatom1,lm2,iatom2),lm1=1,ilm1_max),iatom1=1,natom),lm2=1,ilm2_max),iatom2=1,natom)

    ! Writing into output array Sk that is a normal matrix
    Sk=0
    DO p=1,lndim
       iatom = lbndind(p,1)
       l = lbndind(p,2)
       m = lbndind(p,3)
       s = lbndind(p,4)
       lm = l*(l+1)+m+1
       DO q=1,ndim
          iatom1 = bndind(q,1)
          l1 = bndind(q,2)
          m1 = bndind(q,3)
          s1 = bndind(q,4)
          lm1 = l1*(l1+1)+m1+1
          if (s.NE.s1) cycle
          Sk(p,q) = Sk_temp(lm,iatom,lm1,iatom1)
       ENDDO
    ENDDO
  END SUBROUTINE readS

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  SUBROUTINE close_file()
    IMPLICIT NONE
    IF (ALLOCATED(Sk_temp)) DEALLOCATE (Sk_temp)
    CLOSE(fh_Sk)
  END SUBROUTINE close_file

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

END MODULE msk
