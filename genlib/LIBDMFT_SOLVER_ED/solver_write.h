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


  SUBROUTINE write_info_solver(UNIT)
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    INTEGER                       :: unit_

    CALL dump_message(UNIT=UNIT,TEXT="############################")
    CALL dump_message(UNIT=UNIT,TEXT="### ED-SOLVER PARAMETERS ###")
    CALL dump_message(UNIT=UNIT,TEXT="############################")
                      unit_ = log_unit
    IF(PRESENT(UNIT)) unit_ = UNIT
    CALL dump_message(UNIT=UNIT,TEXT="# NUMBER OF LANCZOS ITERATIONS = "//c2s(i2c(Nitermax)))

    IF(Block_size==0)THEN
      CALL dump_message(UNIT=UNIT,TEXT=  "# USING PLAIN LANCZOS")
    ELSE
#ifdef _complex
      ! THERE IS NO AVAILABLE BLOCK LANCZOS PACKAGE FOR HERMITIAN MATRICES!
      CALL dump_message(UNIT=UNIT,TEXT=  "# USING CULLUM LANCZOS FOR COMPLEX HAMILTONIAN")
#else
      IF(Block_size<0)THEN
        CALL dump_message(UNIT=UNIT,TEXT="# USING CULLUM LANCZOS FOR REAL HAMILTONIAN")
      ELSEif(Block_size>0)THEN
        CALL dump_message(UNIT=UNIT,TEXT="# USING BLOCK LANCZOS FOR REAL HAMILTONIAN WITH BLOCK SIZE="//c2s(i2c(Block_size)))
      ENDIF
#endif
    ENDIF

#ifndef _complex
    IF(Neigen>1.and.Block_size>0) CALL dump_message(UNIT=UNIT,TEXT="# BLOCK SIZE = "//c2s(i2c(Block_size)))
    WRITE(unit_,'(a,E14.7)') "# LANCZOS TOLERANCE = ",tolerance
#endif

    WRITE(unit_,'(a,E14.7,a)') "# RETAIN THE "//c2s(i2c(Neigen))//" LOWEST EIGENSTATES WITHIN ",dEmax," ABOVE GROUND STATE" 
    CALL dump_message(UNIT=unit_,TEXT="# IN ...")
    CALL print_eigensectorlist(sector2diagH,UNIT=UNIT)

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
