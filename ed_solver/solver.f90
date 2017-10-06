MODULE solver

  use correlations 
  use namelistmod 
  use Lanczos_Cullum_wrapper
  use vertex

  IMPLICIT NONE

  TYPE(eigensectorlist_type),  PRIVATE,SAVE :: sector2diagH      ! sectors to diagonalize 
  TYPE(eigensectorlist_type),  PRIVATE,SAVE :: oldGS             ! GS from previous run

  REAL(DBL),    PARAMETER,     PRIVATE      :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL
  LOGICAL,      PARAMETER,     PRIVATE      :: F=.FALSE.,T=.TRUE.
  LOGICAL,                     PRIVATE      :: start_from_old_gs = F

CONTAINS

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

#include "solver_init.h"
#include "solver_write.h"

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

  SUBROUTINE solve_AIM(GS,AIM,OLDGSFILE,COMPUTE_ALL_CORREL,skip_vector,skip_filter)

    !---------------------!    
    ! AIM DIAGONALIZATION !
    !---------------------!

    TYPE(eigensectorlist_type), INTENT(INOUT) :: GS   
    TYPE(AIM_type),             INTENT(IN)    :: AIM
    INTEGER                                   :: isector
    REAL(DBL)                                 :: E0
    LOGICAL, SAVE                             :: first_call = T 
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL    :: OLDGSFILE
    INTEGER                                   :: uup,ddn,nstu,nstd,itot,ssz
    LOGICAL                                   :: NOT_COMMENSURATE
    LOGICAL,OPTIONAL                          :: COMPUTE_ALL_CORREL
    LOGICAL,OPTIONAL                          :: skip_vector,skip_filter

    itot=AIM%bath%Nb+AIM%impurity%Nc

     if(size2>1.and..not.no_mpi)then
      write(*,*) 'AIM SOLVER : start SOLVE AIM, RANK = ', rank
      call mpibarrier
     endif

    IF(first_call.AND.start_from_old_gs)THEN
      !--------------------------------------------------!
      ! READ GS FROM PREVIOUS RUN (FIRST ITERATION ONLY) !
      !--------------------------------------------------!
       write(*,*) ' FIRST_CALL AND STARTING FROM OLD GS - RANK ', rank
       write(log_unit,*) '- READ GS FROM PREVIOUS RUN -'
       call dumpmess1
       CALL copy_eigensectorlist(GS,oldGS)
       CALL delete_eigensectorlist (oldGS)
    ELSE
      !-------------!
      ! DIAGONALIZE !
      !-------------!
       write(*,*) ' AIM SOLVER : NOT STARTING FROM OLD GS - RANK ', rank
       call dumpmess2
       write(log_unit,*) '- INITIALIZE GS TO LIST OF EMPTY SECTORS -'
       CALL copy_eigensectorlist(GS,sector2diagH) 
       CALL get_eigen_energies
       if(AIM%bath%SUPER) call nambu_energy_shift(GS)
       E0 = GSenergy(GS)
       if(.not.FLAG_FULL_ED_GREEN.and..not.present(skip_filter)) CALL filter_eigensector( GS , [E0,E0+dEmax] )
       CALL print_eigensectorlist(GS)

       if(AIM%bath%SUPER) call nambu_energy_shift(GS,back=.true.)

       if(.not.(which_lanczos=='FULL_ED'.or.which_lanczos=='ARPACK').and..not.present(skip_vector)) call get_GS

    ENDIF

    call get_new_bounds

    IF(present(OLDGSFILE).and.iproc==1.AND.dump_ground_state &
    & .AND.((.NOT.first_call).OR.(.NOT.start_from_old_gs))) then
         write(*,*) 'WRITING RAW EIGENSECTOR LIST'
         CALL write_raw_eigensectorlist(GS,OLDGSFILE)
    endif
    IF(first_call.and..not.ALL_FIRST_CALL) first_call = F 

  contains

#include "solver.h"

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

END MODULE 



