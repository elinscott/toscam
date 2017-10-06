MODULE rtm

 
  complex(8), ALLOCATABLE :: Ug(:,:,:)       ! Transformation matrix for rot-mat withought phase factor: Ug(:,:,group)
  complex(8), ALLOCATABLE :: Ugll(:,:,:)     ! Transformation matrix for rot-mat withought phase factor: Ugll(:,:,group)
  complex(8), ALLOCATABLE :: Ph(:,:,:,:)     ! Phase factor vector for rot-mat: Ph(:,k1,k2,k3)
  complex(8), ALLOCATABLE :: SJ(:,:,:)
  INTEGER :: norb, lnorb, ndim, ngroup, nspin ! local copy of dimansions
  INTEGER :: lndim


CONTAINS

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

  subroutine clean_up()
    if (allocated(Ug)) deallocate (Ug)
    if (allocated(Ugll)) deallocate (Ugll)
    if (allocated(Ph)) deallocate (Ph) 
    if (allocated(SJ)) deallocate (SJ)
  end subroutine clean_up

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

  subroutine init(rdc, bk, wigmat, ltpan, ilsb, ilst, igatom, lmtind, llmtind, is, ikov, SJ_, &
       norb_, lnorb_, ndim_, lndim_, nspin_, kirrptr, gropptr, &
       ngroup_, natom, nkp, ng, imax, lmax, nsort, ilsbmax2, illsbmax2, ndiv1, ndiv2, ndiv3)  

!!!--------------------------------------------------------------------------------------------------!!!
!!! Allocates memory for global transformation matrices Ug and Ph and initializes them.              !!!
!!! The matrix Ug depends only on group operation and Ph is a vector (equivalent to diagonal matrix) !!! 
!!! of phase factors and depends only on k point.                                                    !!!
!!! The rotation matrix is finaly : Ug . Ph' , where . is a matrix product and Ph'                   !!!
!!! is assumed to be digonal matrix build from the vector Ph.                                        !!!
!!!--------------------------------------------------------------------------------------------------!!!

    IMPLICIT NONE

    real(8),    INTENT(in) :: rdc(3,natom,ngroup_) ! coordinate of the atom after rotation
    real(8),    INTENT(in) :: bk(3,nkp)           ! k-point coordinate iside 1BZ
    complex(8), INTENT(in) :: wigmat(imax,ng)     ! Wigner matrices
    INTEGER,    INTENT(in) :: ltpan(0:lmax,natom) ! l-presence array
    INTEGER,    INTENT(in) :: ilsb(nsort)         ! maximum l in base
    INTEGER,    INTENT(in) :: ilst(nsort)         ! maximum l in base
    INTEGER,    INTENT(in) :: igatom(natom,ngroup_)! index of the atom after rotation
    INTEGER,    INTENT(in) :: lmtind(ilsbmax2,natom)
    INTEGER,    INTENT(in) :: llmtind(illsbmax2,natom)
    INTEGER,    INTENT(in) :: is(natom)
    INTEGER,    INTENT(in) :: ikov(ngroup_)
    INTEGER,    intent(in) :: kirrptr(ndiv1,ndiv2,ndiv3)
    INTEGER,    intent(in) :: gropptr(ndiv1,ndiv2,ndiv3)
    complex(8), intent(in) :: SJ_(2,2,ngroup_)
    INTEGER,    intent(in) :: norb_, lnorb_, ndim_, lndim_, nspin_, ngroup_, natom, nkp, ng, imax, lmax, nsort, ilsbmax2, illsbmax2
    INTEGER,    intent(in) :: ndiv1, ndiv2, ndiv3
    INTEGER                :: i, ikp, ir, i1, i2, i3
    complex(8),ALLOCATABLE :: unity(:,:)

    ngroup = ngroup_
    norb   = norb_
    lnorb  = lnorb_
    ndim   = ndim_
    lndim  = lndim_
    nspin  = nspin_
    
! Allocate module memory - will be shared among members of the module
    allocate( SJ(2,2,ngroup))
    allocate( Ug(norb,norb,ngroup))
    allocate( Ugll(lnorb,lnorb,ngroup))
    allocate( Ph(norb,ndiv1,ndiv2,ndiv3))

! Initialize the rest of module members
    SJ = SJ_

    allocate( unity(norb,norb))

    unity=0
    do i=1,norb
       unity(i,i)=1.0
    enddo

    do i=1,ngroup
       call rot_mat_right_nophase_(Ug(:,:,i), unity, .FALSE., i, &
            wigmat, ltpan, ilsb, igatom, lmtind, is, ikov, norb, natom, ngroup, ng, imax, lmax, nsort, ilsbmax2)
       ! call rot_mat_right_nophase(unity,Ug(:,:,i),i)
    enddo
    
    deallocate(unity)
    allocate(unity(1:lnorb,1:lnorb))
    unity=0
    do i=1,lnorb
       unity(i,i)=1.0
    enddo

    do i=1,ngroup
       call rot_mat_right_nophase_(Ugll(:,:,i), unity, .TRUE., i, &
            wigmat, ltpan, ilst, igatom, llmtind, is, ikov, lnorb, natom, ngroup, ng, imax, lmax, nsort, illsbmax2)
       ! call rot_mat_right_nophase_long_long(unity,Ugll(:,:,i),i)
    enddo

    do i3=1,ndiv3
       do i2=1,ndiv2
          do i1=1,ndiv1
             ikp = kirrptr(i1,i2,i3)
             ir  = gropptr(i1,i2,i3)
             call phase_right(Ph(:,i1,i2,i3),ikp,ir, &
                  rdc,bk,ltpan,ilsb,lmtind,is,ikov,norb,natom,nkp,ngroup,lmax,nsort,ilsbmax2)
             ! call phase_right(Ph(:,i1,i2,i3),ikp,ir)
          enddo
       enddo
    enddo
    deallocate(unity)
  endsubroutine init


!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

  subroutine cmp_Ukn_phase(Uk, ir, i1,i2,i3)

    !------------------------------------------------------------------------------------------------
    ! Returns the rotaion matrix Uk that can be used to rotate Hamiltonian as follows: Uk^+ . H_k . Uk
    ! The rotation matrix is obtained by : Ug . Ph' , where . is a matrix product and Ph'
    ! is assumed to be digonal matrix build from the vector Ph.
    !------------------------------------------------------------------------------------------------

    IMPLICIT NONE

    complex(8), intent(out) :: Uk(ndim,ndim)        ! output transformation matrix 
    INTEGER, intent(in)     :: i1, i2, i3           ! k-point
    INTEGER, intent(in)     :: ir                   ! group operation
    INTEGER                 :: i, j
    complex(8)              :: Ukw(1:norb,1:norb)

    IF (.NOT. allocated(Ug) .OR. .NOT. allocated(Ph)) THEN
       WRITE(0,*) 'In module "rtm" memory not allocated. Call rtm.init() first '
       return
    ENDIF

    Uk=0
    do j=1,norb
       do i=1,norb
          Ukw(i,j) = Ug(i,j,ir)*Ph(j,i1,i2,i3)
       enddo
    enddo
    if (ndim.eq.norb) then
       Uk = conjg(Ukw)
    else
       Uk(1:norb,1:norb) = conjg(Ukw)*SJ(1,1,ir)
       Uk(norb+1:ndim,norb+1:ndim) = conjg(Ukw)*SJ(2,2,ir)
       Uk(1:norb,norb+1:ndim) = conjg(Ukw)*SJ(2,1,ir)
       Uk(norb+1:ndim,1:norb) = conjg(Ukw)*SJ(1,2,ir)
    endif

  return
  end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

  subroutine cmp_Ukn(Uk, ir)

    !------------------------------------------------------------------------------------------------
    ! Returns the rotaion matrix Uk that can be used to rotate Hamiltonian as follows: Uk^+ . H_k . Uk
    ! The rotation matrix is obtained by : Ug . Ph' , where . is a matrix product and Ph'
    ! is assumed to be digonal matrix build from the vector Ph.
    !------------------------------------------------------------------------------------------------

    IMPLICIT NONE
    complex(8), intent(out) :: Uk(1:ndim,1:ndim)    ! output transformation matrix 
    INTEGER, intent(in)     :: ir                   ! group operation
    INTEGER                 :: i, j

    IF (.NOT. allocated(Ug) .OR. .NOT. allocated(SJ)) THEN
       WRITE(0,*) 'In module "rtm" memory not allocated. Call rtm.init() first '
       return
    ENDIF

    Uk=0
    if (ndim.eq.norb) then
       Uk = conjg(Ug(:,:,ir))
    else
       Uk(1:norb,1:norb) = conjg(Ug(:,:,ir))*SJ(1,1,ir)
       Uk(norb+1:ndim,norb+1:ndim) = conjg(Ug(:,:,ir))*SJ(2,2,ir)
       Uk(1:norb,norb+1:ndim) = conjg(Ug(:,:,ir))*SJ(2,1,ir)
       Uk(norb+1:ndim,1:norb) = conjg(Ug(:,:,ir))*SJ(1,2,ir)
    endif
    return
  end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************


  subroutine cmp_Tall(Tall, T2CJ, gropptr, ndiv1, ndiv2, ndiv3)

    !---------------------------------------------------------------------
    ! Calculates the large Transformation matrix for each k-point
    ! that rotates the matrix to particular k-point and at the same time 
    ! transform to the cubic harmonics : Tall = Ug . Ph' . T2c
    !---------------------------------------------------------------------

    IMPLICIT NONE
    complex(8), intent(out) :: Tall(ndim,ndim,ndiv1,ndiv2,ndiv3)
    complex(8), intent(in)  :: T2CJ(ndim,ndim)
    integer, intent(in)     :: gropptr(ndiv1,ndiv2,ndiv3)
    INTEGER, intent(in)     :: ndiv1, ndiv2, ndiv3
    complex(8)              :: UtoJ(ndim,ndim)
    complex(8)              :: R(1:ndim,1:ndim)
    INTEGER                 :: i3, i2, i1, ir, i, j, iat
    real(8)                 :: j0, mj, l, ml, ms, s0
    do i3=1,ndiv3
       do i2=1,ndiv2
          do i1=1,ndiv1
             ir  = gropptr(i1,i2,i3)
             CALL cmp_Ukn_phase(R,ir,i1,i2,i3)
             ! Add transformation to cubic harmonics
             Tall(:,:,i1,i2,i3) = matmul(R,T2CJ)
          enddo
       enddo
    enddo
  end subroutine

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

subroutine rot_mat(matout, matin, ikp, ir, rdc, bk, wigmat, ltpan, ilsb, igatom, lmtind, is, ikov, &
                   norb, natom, nkp, ngroup, ng, imax, lmax, nsort, ilsbmax2)

!!!--------------------------------------------------------------------------!!!
!!!  This routine can be used to rotate a matrix, such as the Hamiltonian or !!!
!!!  the overlap matrix.  The original matrix which is assumed to correspond !!!
!!!  to k point ikp, is stored in matin and the rotated matrix is stored in  !!!
!!!  matout.  The matrices are stored according to the convention used in    !!!
!!!  Sergei Savrasov's program and we use the Wigner matrices from his code. !!!
!!!--------------------------------------------------------------------------!!!

  IMPLICIT NONE

  INTEGER,    INTENT(in)  :: ikp                 ! Index to irreducible k-point
  INTEGER,    INTENT(in)  :: ir                  ! Index to group operation
  complex(8), INTENT(out) :: matout(norb,norb)   ! input matrix
  complex(8), INTENT(in)  :: matin(norb,norb)    ! output rotated matrix
  real(8),     INTENT(in) :: rdc(3,natom,ngroup) ! coordinate of the atom after rotation
  real(8),     INTENT(in) :: bk(3,nkp)           ! k-point coordinate iside 1BZ
  complex(8), INTENT(in)  :: wigmat(imax,ng)     ! Wigner matrices
  INTEGER,    INTENT(in)  :: ltpan(0:lmax,natom) ! l-presence array
  INTEGER,    INTENT(in)  :: ilsb(nsort)         ! maximum l in base
  INTEGER,    INTENT(in)  :: igatom(natom,ngroup)! index of the atom after rotation
  INTEGER,    INTENT(in)  :: lmtind(ilsbmax2,natom)
  INTEGER,    INTENT(in)  :: is(natom)
  INTEGER,    INTENT(in)  :: ikov(ngroup)
  INTEGER                 :: norb, natom, nkp, ngroup, ng, imax, lmax, nsort, ilsbmax2
  INTEGER                 :: i                 ! Loop index
  INTEGER                 :: lm                ! Combined angular momentum index
  INTEGER                 :: ig                ! Group operation index for Wigner matrices 
  INTEGER                 :: iatom1,iatom2     ! Pointers to atoms
  INTEGER                 :: ratom1,ratom2     ! Pointers to rotated atoms
  INTEGER                 :: isort1,isort2     ! Atomic sorts
  INTEGER                 :: l1,l2             ! Total angular momentum indices
  INTEGER                 :: m1,m2             ! Angular momentum component indices
  INTEGER                 :: m1s,m2s           ! Angular momentum component indices for sum
  INTEGER                 :: orb1,orb2         ! Matrix orbit indices
  INTEGER                 :: sorb1,sorb2       ! Matrix orbit indices for sum 
  INTEGER                 :: cidx1,cidx2       ! Complex indices for Wigner matrices 
  real(8)                 :: r1(1:3)           ! Coordinates of shift vector for atom1
  real(8)                 :: r2(1:3)           ! Coordinates of shift vector for atom2
  real(8)                 :: kirr(1:3)         ! Coordinates of irreducible k-point
  real(8)                 :: rtmp              ! Real dummy variable
  real(8)                 :: pi                ! Good old pi
  complex(8)              :: arg            ! Argument for exponential factor
  complex(8)              :: expfac         ! Exponential factor 
  complex(8)              :: ctmp1          ! Complex dummy variable
  complex(8)              :: ctmp2          ! Complex dummy variable

  pi = 4*ATAN(1.d0)

  !---------- Change from index ir to ig ----------
  ig = ikov(ir)

  !---------- Create matrix element of rotated matrix ----------
  matout = 0.d0

  !---------- Loop over atoms ----------
  DO iatom1 = 1,natom
     DO iatom2 = 1,natom
        !---------- Compute exponential factor ----------
        r1 = rdc(:,iatom1,ir)
        r2 = rdc(:,iatom2,ir)
        kirr = bk(:,ikp)
        rtmp = 0.d0
        DO i = 1,3
           rtmp = rtmp+kirr(i)*(r1(i)-r2(i))
        ENDDO
        arg = cmplx(0.d0,2*pi*rtmp,8)
        expfac = EXP(arg)
        !---------- Find rotated atoms ----------
        ratom1 = igatom(iatom1,ir)
        ratom2 = igatom(iatom2,ir)
        !---------- Find what sort the atoms are ----------
        isort1 = is(iatom1)
        isort2 = is(iatom2)
        !---------- Loop over total angular momentum ----------
        DO l1 = 0,ilsb(isort1)
           !---------- Cycle if this l not present in matrix ----------
           IF(.NOT.ltpan(l1,iatom1)==1) CYCLE
           DO l2 = 0,ilsb(isort2)
              !---------- Cycle if this l not present in matrix ----------
              IF(.NOT.ltpan(l2,iatom2)==1) CYCLE
              !---------- Loop over angular momentum component ----------
              DO m1 = -l1,l1
                 DO m2 = -l2,l2
                    !---------- Find corresponding orbits of matrix ----------
                    lm = l1*(l1+1)+m1+1
                    orb1 = lmtind(lm,iatom1)
                    lm = l2*(l2+1)+m2+1
                    orb2 = lmtind(lm,iatom2)
                    !---------- Sum over angular momentum components ----------
                    ctmp1 = 0.d0
                    DO m1s = -l1,l1
                       !---------- Find orbits of matrix in sum ----------
                       lm = l1*(l1+1)+m1s+1
                       sorb1 = lmtind(lm,ratom1)
                       !---------- Find complex index for Wigner matrices ----------
                       cidx1 = l1*(2*l1-1)*(2*l1+1)/3+(2*l1+1)*(l1+m1s)+l1+m1+1
                       ctmp2 = 0.d0
                       DO m2s = -l2,l2
                          !---------- Find orbits of matrix in sum ----------
                          lm = l2*(l2+1)+m2s+1
                          sorb2 = lmtind(lm,ratom2)
                          !---------- Find complex index for Wigner matrices ----------
                          cidx2 = l2*(2*l2-1)*(2*l2+1)/3+(2*l2+1)*(l2+m2s)+l2+m2+1
                          !---------- Do the  sum ----------
                          ctmp2 = ctmp2+matin(sorb1,sorb2)*CONJG(wigmat(cidx2,ig))
                       ENDDO
                       ctmp1 = ctmp1+ctmp2*wigmat(cidx1,ig)
                    ENDDO
                    !---------- Sum is finished, store results ----------
                    matout(orb1,orb2) = ctmp1*expfac
                 ENDDO
              ENDDO
              !---------- End of angular momentum component loop ----------
           ENDDO
        ENDDO
        !---------- That was the end of the total angular momentum loops ----------
     ENDDO
  ENDDO
  !---------- And this was the end the atom loops ----------

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


subroutine phase_right(matout,ikp,ir,rdc,bk,ltpan,ilsb,lmtind,is,ikov,norb,natom,nkp,ngroup,lmax,nsort,ilsbmax2)

  IMPLICIT NONE

  INTEGER,    INTENT(in)  :: ikp                 ! Index to irreducible k-point
  INTEGER,    INTENT(in)  :: ir                  ! Index to group operation
  complex(8), INTENT(out) :: matout(norb)        ! input matrix
  real(8),     INTENT(in) :: rdc(3,natom,ngroup) ! coordinate of the atom after rotation !
  real(8),     INTENT(in) :: bk(3,nkp)           ! k-point coordinate iside 1BZ          !
  INTEGER,    INTENT(in)  :: ltpan(0:lmax,natom) ! l-presence array                      !
  INTEGER,    INTENT(in)  :: ilsb(nsort)       ! maximum l in base                     !
  INTEGER,    INTENT(in)  :: lmtind(ilsbmax2,natom)                                    !
  INTEGER,    INTENT(in)  :: is(natom)         !                                       !
  INTEGER,    INTENT(in)  :: ikov(ngroup)      !                                       !
  INTEGER                 :: norb, natom, nkp, ngroup, lmax, nsort, ilsbmax2
  INTEGER                 :: i                 ! Loop index
  INTEGER                 :: lm                ! Combined angular momentum index
  INTEGER                 :: ig                ! Group operation index for Wigner matrices 
  INTEGER                 :: iatom1            ! Pointers to atoms
  INTEGER                 :: isort1            ! Atomic sorts
  INTEGER                 :: l1                ! Total angular momentum indices
  INTEGER                 :: m1                ! Angular momentum component indices
  INTEGER                 :: orb1              ! Matrix orbit indices
  real(8)                 :: r1(1:3)           ! Coordinates of shift vector for atom1
  real(8)                 :: kirr(1:3)         ! Coordinates of irreducible k-point
  real(8)                 :: rtmp              ! Real dummy variable
  real(8)                 :: pi                ! Good old pi
  complex(8)              :: arg               ! Argument for exponential factor
  complex(8)              :: expfac            ! Exponential factor 

  pi = 4*ATAN(1.d0)
  !---------- Change from index ir to ig ----------
  ig = ikov(ir)
  !---------- Create matrix element of rotated matrix ----------
  matout = 1.d0
  !---------- Loop over atoms ----------
  DO iatom1 = 1,natom
     !---------- Compute exponential factor ----------
     r1 = rdc(:,iatom1,ir)
     kirr = bk(:,ikp)
     rtmp = 0.d0
     DO i = 1,3
        rtmp = rtmp+kirr(i)*(-r1(i))
     ENDDO
     arg = cmplx(0.d0,2*pi*rtmp,8)
     expfac = EXP(arg)
     !---------- Find what sort the atoms are ----------
     isort1 = is(iatom1)
     !---------- Loop over total angular momentum ----------
     DO l1 = 0,ilsb(isort1)
        !---------- Cycle if this l not present in matrix ----------
        IF(.NOT.ltpan(l1,iatom1)==1) CYCLE
        DO m1 = -l1,l1
           !---------- Find corresponding orbits of matrix ----------
           lm = l1*(l1+1)+m1+1
           orb1 = lmtind(lm,iatom1)
           !---------- Sum over angular momentum components ----------
           matout(orb1) = expfac
        ENDDO
     ENDDO
     !---------- End of angular momentum component loop ----------
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

subroutine rot_mat_right_nophase(matout, matin, ir, wigmat, ltpan, ilsb, igatom, lmtind, is, ikov, &
     norb, natom, ngroup, ng, imax, lmax, nsort, ilsbmax2)

!--------------------------------------------------------------------------!
!  This routine can be used to rotate a matrix, such as the Hamiltonian or !
!  the overlap matrix.  The original matrix which is assumed to correspond !
!  to k point ikp, is stored in matin and the rotated matrix is stored in  !
!  matout.  The matrices are stored according to the convention used in    !
!  Sergei Savrasov's program and we use the Wigner matrices from his code. !
!--------------------------------------------------------------------------!

  IMPLICIT NONE

  INTEGER,    INTENT(in) :: ir                  ! Index to group operation
  complex(8), INTENT(out):: matout(norb,norb)   ! input matrix
  complex(8), INTENT(in) :: matin(norb,norb)    ! output rotated matrix
  complex(8), INTENT(in) :: wigmat(imax,ng)     ! Wigner matrices
  INTEGER,    INTENT(in) :: ltpan(0:lmax,natom) ! l-presence array
  INTEGER,    INTENT(in) :: ilsb(nsort)         ! maximum l in base
  INTEGER,    INTENT(in) :: igatom(natom,ngroup)! index of the atom after rotation
  INTEGER,    INTENT(in) :: lmtind(ilsbmax2,natom)
  INTEGER,    INTENT(in) :: is(natom)
  INTEGER,    INTENT(in) :: ikov(ngroup)
  INTEGER                :: norb, natom, ngroup, ng, imax, lmax, nsort, ilsbmax2
  INTEGER                :: i                 ! Loop index
  INTEGER                :: lm                ! Combined angular momentum index
  INTEGER                :: ig                ! Group operation index for Wigner matrices 
  INTEGER                :: iatom1,iatom2     ! Pointers to atoms
  INTEGER                :: ratom2            ! Pointers to rotated atoms
  INTEGER                :: isort1,isort2     ! Atomic sorts
  INTEGER                :: l1,l2             ! Total angular momentum indices
  INTEGER                :: m1,m2             ! Angular momentum component indices
  INTEGER                :: m1s,m2s           ! Angular momentum component indices for sum
  INTEGER                :: orb1,orb2         ! Matrix orbit indices
  INTEGER                :: sorb1,sorb2       ! Matrix orbit indices for sum 
  INTEGER                :: cidx1,cidx2       ! Complex indices for Wigner matrices 
  real(8)                :: r1(1:3)           ! Coordinates of shift vector for atom1
  real(8)                :: r2(1:3)           ! Coordinates of shift vector for atom2
  real(8)                :: kirr(1:3)         ! Coordinates of irreducible k-point
  real(8)                :: rtmp              ! Real dummy variable
  real(8)                :: pi                ! Good old pi
  complex(8)             :: arg               ! Argument for exponential factor
  complex(8)             :: expfac            ! Exponential factor 
  complex(8)             :: ctmp1             ! Complex dummy variable
  complex(8)             :: ctmp2             ! Complex dummy variable

  pi = 4*ATAN(1.d0)
  !---------- Change from index ir to ig ----------
  ig = ikov(ir)
  !---------- Create matrix element of rotated matrix ----------
  matout = 0.d0
  !---------- Loop over atoms ----------
  DO iatom1 = 1,natom
     DO iatom2 = 1,natom
        !---------- Find rotated atoms ----------
        ratom2 = igatom(iatom2,ir)
        !---------- Find what sort the atoms are ----------
        isort1 = is(iatom1)
        isort2 = is(iatom2)
        !---------- Loop over total angular momentum ----------
        DO l1 = 0,ilsb(isort1)
           !---------- Cycle if this l not present in matrix ----------
           IF(.NOT.ltpan(l1,iatom1)==1) CYCLE
           DO l2 = 0,ilsb(isort2)
              !---------- Cycle if this l not present in matrix ----------
              IF(.NOT.ltpan(l2,iatom2)==1) CYCLE
              !---------- Loop over angular momentum component ----------
              DO m1 = -l1,l1
                 DO m2 = -l2,l2
                    !---------- Find corresponding orbits of matrix ----------
                    lm = l1*(l1+1)+m1+1
                    orb1 = lmtind(lm,iatom1)
                    lm = l2*(l2+1)+m2+1
                    orb2 = lmtind(lm,iatom2)
                    !---------- Sum over angular momentum components ----------
                    ctmp2 = 0.d0
                    DO m2s = -l2,l2
                       !---------- Find orbits of matrix in sum ----------
                       lm = l2*(l2+1)+m2s+1
                       sorb2 = lmtind(lm,ratom2)
                       !---------- Find complex index for Wigner matrices ----------
                       cidx2 = l2*(2*l2-1)*(2*l2+1)/3+(2*l2+1)*(l2+m2s)+l2+m2+1
                       ctmp2 = ctmp2+matin(orb1,sorb2)*CONJG(wigmat(cidx2,ig))
                    ENDDO
                    !---------- Sum is finished, store results ----------
                    matout(orb1,orb2) = ctmp2
                 ENDDO
              ENDDO
              !---------- End of angular momentum component loop ----------
           ENDDO
        ENDDO
        !---------- That was the end of the total angular momentum loops ----------
     ENDDO
  ENDDO
  !---------- And this was the end the atom loops ----------

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


! How to call for small or long version:
! long   -> (false/true)
! norb   -> (norb/lnorb)
! ilsx   -> (ilsb/ilst)
! lmtind -> (lmtind/llmtind)


subroutine rot_mat_right_nophase_(matout, matin, long, ir, wigmat, ltpan, ilsx, igatom, lmtind, is, ikov, &
     norb, natom, ngroup, ng, imax, lmax, nsort, ilsbmax2)

!--------------------------------------------------------------------------!
!  This routine can be used to rotate a matrix, such as the Hamiltonian or !
!  the overlap matrix.  The original matrix which is assumed to correspond !
!  to k point ikp, is stored in matin and the rotated matrix is stored in  !
!  matout.  The matrices are stored according to the convention used in    !
!  Sergei Savrasov's program and we use the Wigner matrices from his code. !
!--------------------------------------------------------------------------!

  IMPLICIT NONE
  INTEGER,    INTENT(in) :: ir                  ! Index to group operation
  complex(8), INTENT(out):: matout(norb,norb)   ! input matrix
  complex(8), INTENT(in) :: matin(norb,norb)    ! output rotated matrix
  LOGICAL,    INTENT(in) :: long                ! does it go over l in base (norb) or l in expansion(norb)
  complex(8), INTENT(in) :: wigmat(imax,ng)     ! Wigner matrices
  INTEGER,    INTENT(in) :: ltpan(0:lmax,natom) ! l-presence array
  INTEGER,    INTENT(in) :: ilsx(nsort)         ! maximum l in base
  INTEGER,    INTENT(in) :: igatom(natom,ngroup)! index of the atom after rotation
  INTEGER,    INTENT(in) :: lmtind(ilsbmax2,natom)
  INTEGER,    INTENT(in) :: is(natom)
  INTEGER,    INTENT(in) :: ikov(ngroup)
  INTEGER                :: norb, natom, ngroup, ng, imax, lmax, nsort, ilsbmax2
  INTEGER                :: i                 ! Loop index
  INTEGER                :: lm                ! Combined angular momentum index
  INTEGER                :: ig                ! Group operation index for Wigner matrices 
  INTEGER                :: iatom1,iatom2     ! Pointers to atoms
  INTEGER                :: ratom2            ! Pointers to rotated atoms
  INTEGER                :: isort1,isort2     ! Atomic sorts
  INTEGER                :: l1,l2             ! Total angular momentum indices
  INTEGER                :: m1,m2             ! Angular momentum component indices
  INTEGER                :: m1s,m2s           ! Angular momentum component indices for sum
  INTEGER                :: orb1,orb2         ! Matrix orbit indices
  INTEGER                :: sorb1,sorb2       ! Matrix orbit indices for sum 
  INTEGER                :: cidx1,cidx2       ! Complex indices for Wigner matrices 
  real(8)                :: r1(1:3)           ! Coordinates of shift vector for atom1
  real(8)                :: r2(1:3)           ! Coordinates of shift vector for atom2
  real(8)                :: kirr(1:3)         ! Coordinates of irreducible k-point
  real(8)                :: rtmp              ! Real dummy variable
  real(8)                :: pi                ! Good old pi
  complex(8)             :: arg               ! Argument for exponential factor
  complex(8)             :: expfac            ! Exponential factor 
  complex(8)             :: ctmp1             ! Complex dummy variable
  complex(8)             :: ctmp2             ! Complex dummy variable

  pi = 4*ATAN(1.d0)
  !---------- Change from index ir to ig ----------
  ig = ikov(ir)
  !---------- Create matrix element of rotated matrix ----------
  matout = 0.d0
  !---------- Loop over atoms ----------
  DO iatom1 = 1,natom
     DO iatom2 = 1,natom
        !---------- Find rotated atoms ----------
        ratom2 = igatom(iatom2,ir)
        !---------- Find what sort the atoms are ----------
        isort1 = is(iatom1)
        isort2 = is(iatom2)
        !---------- Loop over total angular momentum ----------
        DO l1 = 0,ilsx(isort1)
           !---------- Cycle if this l not present in matrix ----------
           IF(.NOT. long .AND. .NOT.ltpan(l1,iatom1)==1) CYCLE
           DO l2 = 0,ilsx(isort2)
              !---------- Cycle if this l not present in matrix ----------
              IF(.NOT. long .AND. .NOT.ltpan(l2,iatom2)==1) CYCLE
              !---------- Loop over angular momentum component ----------
              DO m1 = -l1,l1
                 DO m2 = -l2,l2
                    !---------- Find corresponding orbits of matrix ----------
                    lm = l1*(l1+1)+m1+1
                    orb1 = lmtind(lm,iatom1)
                    lm = l2*(l2+1)+m2+1
                    orb2 = lmtind(lm,iatom2)
                    !---------- Sum over angular momentum components ----------
                    ctmp2 = 0.d0
                    DO m2s = -l2,l2
                       !---------- Find orbits of matrix in sum ----------
                       lm = l2*(l2+1)+m2s+1
                       sorb2 = lmtind(lm,ratom2)
                       !---------- Find complex index for Wigner matrices ----------
                       cidx2 = l2*(2*l2-1)*(2*l2+1)/3+(2*l2+1)*(l2+m2s)+l2+m2+1
                       ctmp2 = ctmp2+matin(orb1,sorb2)*CONJG(wigmat(cidx2,ig))
                    ENDDO
                    !---------- Sum is finished, store results ----------
                    matout(orb1,orb2) = ctmp2
                 ENDDO
              ENDDO
              !---------- End of angular momentum component loop ----------
           ENDDO
        ENDDO
        !---------- That was the end of the total angular momentum loops ----------
     ENDDO
  ENDDO
  !---------- And this was the end the atom loops ----------
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


end module
