MODULE Lanczos_Cullum_wrapper

#ifdef _complex
  USE Lanczos_Cullum_hlevec
#else
  USE Lanczos_Cullum_levec
#endif

  IMPLICIT NONE

  CHARACTER(LEN=100), PRIVATE :: VALFILE1  = 'ED_out/cullum_eval_lancmat.out'
  CHARACTER(LEN=100), PRIVATE :: VALFILE3  = 'ED_out/cullum_eval.out'
  CHARACTER(LEN=100), PRIVATE :: VALFILE4  = 'ED_out/cullum_eval_err.out'
  CHARACTER(LEN=100), PRIVATE :: VALFILE5  = 'ED_out/cullum_eval.in'
  CHARACTER(LEN=100), PRIVATE :: VECFILE5  = 'ED_out/cullum_evec.in'
  CHARACTER(LEN=100), PRIVATE :: VECFILE10 = 'ED_out/cullum_evec_bnds.out'

  PRIVATE

  public ::   init_Lanczos_Cullum, Lanczos_Cullum_diagonalize
  
CONTAINS

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE init_Lanczos_Cullum
#ifdef _complex
    CALL init_Lanczos_Cullum_leval
    CALL init_Lanczos_Cullum_hlevec
#else
    CALL init_Lanczos_Cullum_leval
    CALL init_Lanczos_Cullum_levec
#endif
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE Lanczos_Cullum_diagonalize(lowest,eigenvalue_eigenvector)
    TYPE(eigenlist_type), INTENT(INOUT) :: lowest
    INTEGER                             :: start_diagH,eigenvalue_eigenvector

    CALL reset_timer(start_diagH)

#ifdef _complex
    call mpibarrier
    write(*,*) ' ------ Lanczos Cullum compute eigenvalues ------ '

  if(eigenvalue_eigenvector==1)then
    ! COMPUTE EIGENVALUES 
    CALL find_eigenval_c(VALFILE1,VALFILE3,VALFILE4,VALFILE5,lowest)
    write(*,*) '------- Lanczos Cullum compute eigenvectors -------'
  elseif(eigenvalue_eigenvector==2)then
    ! COMPUTE EIGENVECTORS
    call mpibarrier
    CALL find_eigenvec_c(VALFILE1,VALFILE3,VALFILE4,VECFILE5,VECFILE10,lowest)
  endif

    write(*,*) ' ---- done -----'
    CALL timer_fortran(start_diagH,"# DIAGONALIZATION OF "//TRIM(ADJUSTL(title_H_()))//" USING HLEVAL/HLEVEC TOOK ")
    call mpibarrier
#else

  if(eigenvalue_eigenvector==1)then
    ! COMPUTE EIGENVALUES 
    call mpibarrier
    CALL find_eigenval_r(VALFILE1,VALFILE3,VALFILE4,VALFILE5,lowest)
  elseif(eigenvalue_eigenvector==2)then
    ! COMPUTE EIGENVECTORS
    call mpibarrier
    CALL find_eigenvec_r(VALFILE1,VALFILE3,VALFILE4,VECFILE5,VECFILE10,lowest)
  endif

    CALL timer_fortran(start_diagH,"# DIAGONALIZATION OF "//TRIM(ADJUSTL(title_H_()))//" USING LEVAL/LEVEC TOOK ")
    call mpibarrier

#endif

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************



END MODULE
