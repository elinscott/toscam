
MODULE com_scs

!BUG : CW 5 September 2012 was before:
!PRIVATE ::  

PRIVATE 

  INTEGER   :: natom0, nsort0, nspin0, norbs0
  INTEGER, ALLOCATABLE   :: is0(:)
  CHARACTER :: text*64
  real(8)    :: rbas0(3,3), gbas0(3,3), parr0, orth0(2:3), strain0(3,3,2)
  INTEGER   :: nfftg0(3), nplw0
  LOGICAL   :: scfmade
  real(8),     ALLOCATABLE :: tao0(:,:), tau0(:,:)
  real(8) ,    ALLOCATABLE :: znuc0(:), smt0(:), axmag0(:,:)
  INTEGER,    ALLOCATABLE :: nrad0(:), nsym0(:)
  real(8),     ALLOCATABLE :: R0(:,:) !, RHOCOR0(:,:,:), POTFRZ0(:,:,:)
!  real(8),     ALLOCATABLE :: ENY0(:,:,:), DNY0(:,:,:)
  INTEGER,    ALLOCATABLE :: MNY0(:,:), NOD0(:,:)
  INTEGER,    ALLOCATABLE :: lmval(:,:)

!  real(8),     ALLOCATABLE :: RHOCOR1(:,:,:), POTFRZ1(:,:,:)
!  real(8),     ALLOCATABLE :: ENY1(:,:,:), DNY1(:,:,:)

  real(8)    :: rbas1(3,3), gbas1(3,3), parr1, orth1(2:3), strain1(3,3,2)
  INTEGER   :: nfftg1(3), nplw1
  real(8),     ALLOCATABLE :: tao1(:,:), tau1(:,:)

contains

!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************
!***************************************************

  subroutine ReadSCSFile1(filename, inf_buf, &
      &  efmu0, task, RHO0, nrad_max, nsym_max, admix,&
      &  is, natom, nsort, nspin, norbs,ilst_max)
    IMPLICIT NONE
    !Passed varibables
    CHARACTER*500, intent(in)      :: filename
    CHARACTER*10000, intent(inout) :: inf_buf
    integer    :: task,nsort,nspin,norbs,is,admix,ilst_max
    real(8)     :: RHOCOR0(0:nrad_max,nsort,nspin)
    real(8)     :: POTFRZ0(0:nrad_max,nsort,nspin)
    real(8)     :: ENY0(0:ilst_max,nsort,nspin)
    real(8)     :: DNY0(0:ilst_max,nsort,nspin),Ry2ev
    complex(8)  :: RHO0(0:nrad_max,nsym_max,natom,nspin,norbs)
    real(8), intent(out)    :: efmu0
    INTEGER, intent(in)    :: natom
    INTEGER, intent(in)    :: nrad_max, nsym_max
    integer,save :: fh_scs=1


    !Local variables
    INTEGER  :: iatom0, isort0, iorbs0, ispin0, isym0, lm0, irad0, l0, iupdn0
    INTEGER  :: iatom,  isort,  iorbs,  ispin,  isym,  lm,  irad,  l
    INTEGER  :: i, m
    INTEGER  :: nrad_max0, nsym_max0
    INTEGER  :: ilsv0(nsort), nkap0, ilst0(nsort), ilsb0(nsort)
    CHARACTER*256 :: fh_info   

    ! Open the self-consistent file

    fh_scs = fh_scs+1

    OPEN(fh_scs,file=TRIM(filename),status='old')

    WRITE(fh_info,fmt='(1x,A,1x,A1)') 'READ_SCS: Start reading SCFFILE in mode ', task
    inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)

    READ(fh_scs,fmt='(A)') text
    IF(text(2:2).NE.'<') THEN 
       WRITE(0,*) "READ_SCS The format of scs file seems to be wrong"
       WRITE(fh_info,*) "READ_SCS The format of scs file seems to be wrong"
       inf_buf = TRIM(inf_buf)//'\n'//TRIM(fh_info)
       return
    ENDIF

    READ(fh_scs,*) natom0, nsort0, nspin0, norbs0

    IF (.NOT. ALLOCATED(is0)) ALLOCATE(is0(natom0))
    READ(fh_scs,*)(is0(iatom0),iatom0=1,natom0)

    IF (.NOT.ALLOCATED(tao0))   ALLOCATE(tao0(3,natom0))
    IF (.NOT.ALLOCATED(tau0))   ALLOCATE(tau0(3,natom0))
    IF (.NOT.ALLOCATED(znuc0))  ALLOCATE(znuc0(nsort0))
    IF (.NOT.ALLOCATED(smt0))   ALLOCATE(smt0(nsort0))
    IF (.NOT.ALLOCATED(axmag0)) ALLOCATE(axmag0(3,nsort0))
    IF (.NOT.ALLOCATED(nrad0))  ALLOCATE(nrad0(nsort0))
    IF (.NOT.ALLOCATED(nsym0))  ALLOCATE(nsym0(natom0))

    READ(fh_scs,*) RBAS0, GBAS0, PARR0, ORTH0, STRAIN0, ((TAO0(I,IATOM0),I=1,3),IATOM0=1,NATOM0), ((TAU0(I,IATOM0),I=1,3),IATOM0=1,NATOM0),&
         (ZNUC0(ISORT0),ISORT0=1,NSORT0), (SMT0(ISORT0),ISORT0=1,NSORT0), ((AXMAG0(I,ISORT0),I=1,3),ISORT0=1,NSORT0)

    READ(fh_scs,*) (NRAD0(ISORT0),ISORT0=1,NSORT0), (NSYM0(IATOM0),IATOM0=1,NATOM0),(ILSV0(ISORT0),ISORT0=1,NSORT0)
    
    nrad_max0 = maxim(nsort0, nrad0)
    nsym_max0 = maxim(natom0,nsym0)

    IF (.NOT.ALLOCATED(R0)) ALLOCATE(R0(0:nrad_max0,nsort0))
    IF (.NOT.ALLOCATED(lmval)) ALLOCATE(lmval(nsym_max0,nsort0))

    READ(fh_scs,*)((R0(IRAD0,ISORT0),IRAD0=0,NRAD0(ISORT0)),ISORT0=1,NSORT0)

    ! Actual reading of the largest array of total density
    DO IATOM0=1,NATOM0        !!! over old atoms
       ISORT0=IS0(IATOM0)
       DO IORBS0=1,NORBS0     !!! over old iorbs
          DO IUPDN0=1,NSPIN0  !!! over old spins
             DO ISYM0=1,NSYM0(IATOM0) !!! over old symmetry
                READ(fh_scs,*) lmval(isym0,isort0)!LM0
                READ(fh_scs,*)(RHO0(irad0,isym0,iatom0,iupdn0,iorbs0),irad0=0,nrad0(isort0))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO ISORT0=1,NSORT0        !!! over old sorts
       DO ISPIN0=1,NSPIN0     !!! over old spins
          READ(fh_scs,*)(RHOCOR0(IRAD0,isort0,ispin0),IRAD0=0,NRAD0(ISORT0))
          READ(fh_scs,*)(POTFRZ0(IRAD0,isort0,ispin0),IRAD0=0,NRAD0(ISORT0))
       ENDDO
    ENDDO

    IF (.NOT.ALLOCATED(tao1))   ALLOCATE(tao1(3,natom0))
    IF (.NOT.ALLOCATED(tau1))   ALLOCATE(tau1(3,natom0))

    READ(fh_scs,*) RBAS1,GBAS1,ORTH1,PARR1,STRAIN1,NPLW1,NFFTG1, &
     & ((TAO1(I,IATOM0),I=1,3),IATOM0=1,NATOM0),((TAU1(I,IATOM0),I=1,3),IATOM0=1,NATOM0)

    ! Gap here!
    READ(fh_scs,*) NKAP0,(ILST0(ISORT0),ISORT0=1,NSORT0),(ILSB0(ISORT0),ISORT0=1,NSORT0)

    ilst_max = maxim(nsort,ilst0)

    READ(fh_scs,*) &
         (((ENY0(L0,ISORT0,ISPIN0),L0=0,ILST0(ISORT0)),ISORT0=1,NSORT0),ISPIN0=1,NSPIN0),&
         (((DNY0(L0,ISORT0,ISPIN0),L0=0,ILST0(ISORT0)),ISORT0=1,NSORT0),ISPIN0=1,NSPIN0),&
         ((MNY0(L0,ISORT0),L0=0,ILST0(ISORT0)),ISORT0=1,NSORT0),&
         ((NOD0(L0,ISORT0),L0=0,ILST0(ISORT0)),ISORT0=1,NSORT0)

    READ(fh_scs,*) efmu0, scfmade

    Ry2ev=1.2**2
    efmu0 = efmu0 * Ry2eV

    CLOSE(fh_scs)

  end subroutine ReadSCSFile1

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

FUNCTION maxim(nsort, nrad)
  IMPLICIT NONE
  INTEGER :: maxim
  ! Passed variables
  INTEGER :: nsort
  INTEGER :: nrad(nsort)
  ! Local variables
  INTEGER :: nrad_max, isort
  nrad_max = nrad(1)
  DO isort=2,nsort
     IF (nrad(isort).GT.nrad_max) nrad_max = nrad(isort)
  ENDDO
  maxim = nrad_max
END FUNCTION maxim

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

END MODULE com_scs
