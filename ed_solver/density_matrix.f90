MODULE density_matrix

  use HAIMupdo_class


  IMPLICIT NONE
 
  private

  REAL(DBL), PARAMETER, PRIVATE :: zero=0.0_DBL,one=1.0_DBL
  LOGICAL,   PARAMETER, PRIVATE :: F=.FALSE.,T=.TRUE.

  !----------------------------------------------------!
  ! COMPUTE THE REDUCED DENSITY MATRIX OF THE IMPURITY !
  !----------------------------------------------------!

CONTAINS

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE compute_density_matrix(dmat,AIM,beta,GS)

    !-----------------------------------------------------!
    ! COMPUTE THE REDUCED DENSITY MATRIX OF THE IMPURITY  !
    !-----------------------------------------------------!

    TYPE(rcmatrix_type),        INTENT(INOUT) :: dmat
    TYPE(AIM_type),             INTENT(IN)    :: AIM
    REAL(DBL),                  INTENT(IN)    :: beta
    TYPE(eigensectorlist_type), INTENT(IN)    :: GS      ! list of lowest eigenstates
    REAL(DBL)                                 :: boltz,Zpart,E0
    INTEGER                                   :: AIMrank1,AIMrank2,IMPrank1,IMPrank2,BATHrank
    INTEGER                                   :: AIMstate1,AIMstate2
    INTEGER                                   :: nIMPstates,nBATHstates,Nb,Nc
    INTEGER                                   :: start_compute
    INTEGER                                   :: IMPchunk(nproc),IMPstatemin(nproc),IMPstatemax(nproc)
#ifdef _complex 
    COMPLEX(DBL)                              :: coeff1,coeff2 
    COMPLEX(DBL), ALLOCATABLE                 :: dmat_vec(:),dmat_vec_tot(:) 
#else
    REAL(DBL)                                 :: coeff1,coeff2 
    REAL(DBL),    ALLOCATABLE                 :: dmat_vec(:),dmat_vec_tot(:) 
#endif
    INTEGER,      ALLOCATABLE                 :: rankmin(:),rankmax(:),rankchunk(:) 
    INTEGER                                   :: thisrank,isector,ieigen 
    TYPE(eigensector_type), POINTER           :: es    => NULL()
    TYPE(eigen_type),       POINTER           :: eigen => NULL()
    
    CALL dump_message(TEXT="# START COMPUTING THE REDUCED DENSITY MATRIX... ")
    CALL reset_timer(start_compute)
    
    nIMPstates  = AIM%impurity%nstates
    Nc          = AIM%impurity%Nc
    nBATHstates = AIM%bath    %nstates
    Nb          = AIM%bath    %Nb

    CALL split(nIMPstates,IMPstatemin,IMPstatemax,IMPchunk)

    Zpart       = partition(beta,GS)
    E0          = GSenergy(GS)
    dmat%rc     = zero

    !=========================================================================================!

    DO isector=1,GS%nsector

      es => GS%es(isector)

      !----------------------------------------------!
      ! PARSE THE LIST OF EIGENSTATES IN THIS SECTOR !
      !----------------------------------------------!

      DO ieigen=1,es%lowest%neigen

        eigen => es%lowest%eigen(ieigen)

        boltz   = DEXPc(-beta*(eigen%val-E0)) ! Boltzman factor

        !----------------------------------------------!
        ! LOOP OVER MATRIX ELEMENTS WE WANT TO COMPUTE !
        !----------------------------------------------!

!$OMP PARALLEL PRIVATE(IMPrank1,IMPrank2,BATHrank,AIMstate1,AIMrank1,AIMstate2,AIMrank2,coeff1,coeff2)
!$OMP DO
        DO IMPrank1=IMPstatemin(iproc),IMPstatemax(iproc)
          DO IMPrank2=1,IMPrank1

            !---------------------------------------!
            ! TRACE OUT THE BATH DEGREES OF FREEDOM !
            !---------------------------------------!

        !----------------------------------------------------------------------------------------------!
            DO BATHrank=1,nBATHstates
              CALL IMPBATH2AIMstate(AIMstate1,IMPrank1-1,BATHrank-1,Nc,Nb) 
              IF(is_in_sector(AIMstate1,es%sector))THEN
               AIMrank1 = rank_func(AIMstate1,es%sector)
               coeff1   = eigen%vec%rc(AIMrank1)
                IF(coeff1/=zero)THEN
                  CALL IMPBATH2AIMstate(AIMstate2,IMPrank2-1,BATHrank-1,Nc,Nb)
                  IF(is_in_sector(AIMstate2,es%sector))THEN
                   AIMrank2 = rank_func(AIMstate2,es%sector)
                   coeff2   = eigen%vec%rc(AIMrank2)
                    IF(coeff2/=zero)THEN
                      dmat%rc(IMPrank1,IMPrank2) = dmat%rc(IMPrank1,IMPrank2) +      coeff1  * conj(coeff2) * boltz
                      IF(IMPrank1/=IMPrank2) &
                      dmat%rc(IMPrank2,IMPrank1) = dmat%rc(IMPrank2,IMPrank1) + conj(coeff1) *       coeff2  * boltz
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO ! end loop over bath states
        !----------------------------------------------------------------------------------------------!

          ENDDO 
        ENDDO ! loop over impurity states
!$OMP END DO
!$OMP END PARALLEL 
      ENDDO ! loop over eigenstates

    ENDDO ! loop over eigensectors
    !=========================================================================================!


    write(*,*) '............ end main loops density matrix ............'

    if(size2>1.and..not.no_mpi) call mpi_collect_

 
    !-------------------------------------------------------!
    ! RENORMALIZE USING PARTITION FUNCTION                  !
    ! 1 if T=0, 1 also at T>0 if we had all the eigenstates !
    !-------------------------------------------------------!

    dmat%rc = dmat%rc / Zpart 

    CALL timer_fortran(start_compute,"# ... TOOK ")



  contains

  !--------------!
  !--------------!
  !--------------!
  !--------------!

   subroutine mpi_collect_

    write(*,*) 'start collecting MPI chunks'

    ALLOCATE(rankmin(nproc),rankmax(nproc),rankchunk(nproc))
    rankmin   = IMPstatemin * ( IMPstatemin - 1 ) / 2 + 1
    rankmax   = IMPstatemax * ( IMPstatemax + 1 ) / 2
    rankchunk = rankmax - rankmin + 1

    ALLOCATE(dmat_vec(rankchunk(iproc)),dmat_vec_tot(rankmax(nproc)))

    DO IMPrank1=IMPstatemin(iproc),IMPstatemax(iproc)
      thisrank = IMPrank1*(IMPrank1-1)/2
      DO IMPrank2=1,IMPrank1
        dmat_vec_tot(thisrank+IMPrank2)                  = dmat%rc(IMPrank1,IMPrank2)
        dmat_vec    (thisrank+IMPrank2-rankmin(iproc)+1) = dmat%rc(IMPrank1,IMPrank2)
      ENDDO
    ENDDO

write(*,*) 'COLLECTING DATA IN DENSITY MATRIX'

if(size2>1.and..not.no_mpi)then
#ifdef _complex
      CALL MPI_ALLGATHERV( &
      & dmat_vec,rankchunk(iproc),MPI_DOUBLE_COMPLEX,dmat_vec_tot,rankchunk,rankmin-1,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
#else
      CALL MPI_ALLGATHERV( &
      & dmat_vec,rankchunk(iproc),MPI_DOUBLE_PRECISION,dmat_vec_tot,rankchunk,rankmin-1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
#endif
endif

    IF(ierr/=0.and.size2>1.and..not.no_mpi) CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)

    DO IMPrank1=1,nIMPstates
 IF(IMPrank1<IMPstatemin(iproc).OR.IMPrank1>IMPstatemax(iproc))THEN
      thisrank = IMPrank1*(IMPrank1-1)/2
      dmat%rc(IMPrank1,IMPrank1)   = dmat_vec_tot(thisrank+IMPrank1)
      DO IMPrank2=1,IMPrank1-1
        dmat%rc(IMPrank1,IMPrank2) = dmat_vec_tot(thisrank+IMPrank2)
        dmat%rc(IMPrank2,IMPrank1) = conj(dmat%rc(IMPrank1,IMPrank2))
      ENDDO
    ENDIF
 ENDDO

    if(allocated(dmat_vec)    ) deallocate(dmat_vec)
    if(allocated(dmat_vec_tot)) deallocate(dmat_vec_tot)
    if(allocated(rankmin))      deallocate(rankmin)
    if(allocated(rankmax))      deallocate(rankmax)
    if(allocated(rankchunk))    deallocate(rankchunk)

  end subroutine

  !--------------!
  !--------------!
  !--------------!
  !--------------!

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE analyze_density_matrix(dmat,IMPiorb,NAMBU) 

    !--------------------------------------------------------------!
    ! COMPUTE & DISPLAY THE FIRST EIGENPAIRS OF THE DENSITY MATRIX !
    ! COMPUTE THE ENTANGLEMENT ENTROPY:                            !
    ! S = - SUM_i ( lambda_i * log(lambda_i) )                     !
    ! WHERE lambda_i ARE THE EIGENVALUES OF THE DENSITY MATRIX     !
    !--------------------------------------------------------------!

    INTEGER,      INTENT(IN) :: IMPiorb(:,:)
    LOGICAL,      INTENT(IN) :: NAMBU
#ifdef _complex
    COMPLEX(DBL), INTENT(IN) :: dmat(:,:)
    COMPLEX(DBL)             :: VECP(SIZE(dmat,1),SIZE(dmat,1))
#else
    REAL(DBL),    INTENT(IN) :: dmat(:,:)
    REAL(DBL)                :: VECP(SIZE(dmat,1),SIZE(dmat,1))
#endif
    REAL(DBL)                :: VALP(SIZE(dmat,1)),absvec(SIZE(dmat,1)),states(SIZE(dmat,1)) 
    REAL(DBL)                :: ENTROPY,TRACE
    INTEGER                  :: ivp,nvp,icomp
    CHARACTER(LEN=1000)      :: cvec,prefix
    INTEGER                  :: npairs,ncomp,ncompi
    character(100)           :: intwri

    !-------------------------------------------------------------------------------!
    ! NUMBER OF EIGENVALUES TO DISPLAY, NUMBER OF EIGENVECTOR COMPONENTS TO DISPLAY !
    !-------------------------------------------------------------------------------!

    ncomp=min(size(dmat,1),16) 
 npairs=ncomp 

    if(size(dmat,1)<=16) call write_array(dmat," DENSITY MATRIX ", unit=log_unit, short=.true.)

    write(*,*) 'analyse density matrix'
 
    nvp = SIZE(dmat,1)
    CALL diagonalize(dmat,VALP,VECP,EIGENVAL_ONLY=F)
    CALL dump_message(TEXT="# "//c2s(i2c(npairs))//" LARGEST EIGENVALUES OF THE DENSITY MATRIX:") 

    !===============================================================================================!
    DO ivp=nvp,nvp-(npairs-1),-1 ! loop over the npairs largest eigenvalues
      absvec = ABS(VECP(:,ivp)) 
      states = DBLE((/(icomp,icomp=1,nvp)/))
      CALL dsort(absvec,states,nvp,-2) ! sort components of eigenvector in decreasing order and sort state accordingly

      ncompi=ncomp
      do icomp=ncomp,1,-1
       if(abs(VECP(INT(states(icomp)),ivp))>1.d-10)then
        ncompi=icomp
        exit
       endif
      enddo

      write(log_unit,*) '======================================================='
#ifdef _complex
      intwri = '('//c2s(i2c(ncompi))//'(2f9.6,a))'
      write(log_unit,*) 'non zerp elements in state vector  : ', ncompi
      write(log_unit,*)  (VECP(INT(states(icomp)),ivp)," "//cket_from_state(INT(states(icomp))-1,IMPiorb,NAMBU)//" ",icomp=1,ncompi)
      WRITE(cvec,ADJUSTL(TRIM(intwri))) &
      (real(VECP(INT(states(icomp)),ivp)),aimag(VECP(INT(states(icomp)),ivp))," "//cket_from_state(INT(states(icomp))-1,IMPiorb,NAMBU)//" ",icomp=1,ncompi)
#else
      intwri = '('//c2s(i2c(ncompi))//'(f9.6,a))'
      write(log_unit,*) 'non zero elements in state vector  : ', ncompi
      write(log_unit,*)  (VECP(INT(states(icomp)),ivp)," "//cket_from_state(INT(states(icomp))-1,IMPiorb,NAMBU)//" ",icomp=1,ncompi)
      WRITE(cvec,ADJUSTL(TRIM(intwri))) &
      (VECP(INT(states(icomp)),ivp)," "//cket_from_state(INT(states(icomp))-1,IMPiorb,NAMBU)//" ",icomp=1,ncompi)
#endif
      write(log_unit,*) '-------------------------------------------------------'
      prefix = "# "//c2s(i2c(nvp-ivp+1))//"   "
      WRITE(log_unit,'(a'//c2s(i2c(LEN_TRIM(prefix)))//',f9.6,a)') TRIM(prefix),VALP(ivp)," -> "//TRIM(cvec)
      write(log_unit,*) '======================================================='

    ENDDO
    !===============================================================================================!

    write(*,*) '.... start entropy ....'

    write(log_unit,*) '-----------------------------------------'
    write(log_unit,*) 'EIGENVALUES ARE : '
    write(log_unit,'(1000f6.3)') VALP
    write(log_unit,*) '-----------------------------------------'

    ENTROPY = zero
 TRACE = zero
    DO ivp = nvp,1,-1
      IF(VALP(ivp)>1.d-16) ENTROPY = ENTROPY - VALP(ivp) * LOG(VALP(ivp))
      TRACE = TRACE + VALP(ivp)
    ENDDO
    write(*,*) '.... done ....'
    WRITE(log_unit,'(a,f9.6)') "# CHECK TRACE          = ",TRACE 
    WRITE(log_unit,'(a,f9.6)') "# ENTANGLEMENT ENTROPY = ",ENTROPY 
    CALL flush(log_unit)

  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
 
END MODULE 
