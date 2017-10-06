MODULE drive_class

   use genvar
   use common_def

  IMPLICIT NONE


  REAL(DBL),    PARAMETER, PRIVATE                 :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,      PARAMETER, PRIVATE                 :: F=.FALSE.,T=.TRUE.

    
  ! DRIVE QUANTITY TOWARDS TARGET VALUE BY DICHOTOMY
  ! SUBTLETY: QUANTITY IS AN OBSERVABLE (LIKE DENSITY)
  ! BUT DICHOTOMY IS ON IMPURITY HAMILTONIAN CONJUGATE PARAMETER (LIKE CHEMICAL POTENTIAL)
  ! WARNING: OBSERVABLE IS ASSUMED TO BE AN INCREASING FCTN OF PARAM. IS ASSUMED!


  TYPE drive_type
    CHARACTER(LEN=100) :: obsname='\0',paramname='\0'    ! name of the obs./param.
    LOGICAL   :: to_target  = F                          ! T[F] TO DRIVE DENSITY TOWARDS TARGET VALUE
    REAL(DBL) :: target_obs = zero                       ! target value we want to approach
    REAL(DBL) :: param      = zero                       ! current value of the impurity Hamiltonian parameter
    REAL(DBL) :: obs        = zero                       ! current value of the resulting observable
    REAL(DBL) :: obs_error_max=zero,param_error_max=zero ! obs./param. tolerance
    REAL(DBL) :: obs_bnds(2)=zero,param_bnds(2)=zero     ! bounds of obs./param. dichotomy intervals
    INTEGER   :: bnds2find  = -1                         ! # of bounds of dichotomy interval yet to be determined 
  END TYPE 


CONTAINS


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION read_drive(target,targetobs,obserror,paramerror,OBS,PARAM) RESULT(drive)
    TYPE(drive_type)                :: drive
    logical                         :: target 
    REAL(DBL)                       :: targetobs,obserror,paramerror 
    REAL(DBL), OPTIONAL, INTENT(IN) :: OBS,PARAM
    IF(PRESENT(PARAM)) drive%param = PARAM ! starting value for impurity Hamiltonian parameter
    IF(PRESENT(OBS))   drive%obs   = OBS   ! starting value for the resulting impurity observable
    drive%to_target=target
    drive%target_obs=targetobs    
    drive%obs_error_max=obserror
    drive%param_error_max=paramerror
    drive%param_error_max = drive%param_error_max * drive%obs_error_max
    drive%bnds2find = 2
  END FUNCTION

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE write_info_drive(drive,UNIT)
    TYPE(drive_type),  INTENT(IN) :: drive
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    INTEGER :: unit_
    IF(drive%to_target)THEN
      unit_ = log_unit
      IF(PRESENT(UNIT)) unit_ = UNIT
      WRITE(unit_,'(a,f0.6)')   "### ADJUSTING "//TRIM(ADJUSTL(drive%paramname))//" WITH MAX. AMPLITUDE ",drive%param_error_max
      WRITE(unit_,'(2(a,f0.6))') "### TO  DRIVE "//TRIM(ADJUSTL(drive%obsname))//" WITHIN ",drive%obs_error_max," OF ",drive%target_obs
      CALL flush(unit_)
    ENDIF
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION read_raw_drive(UNIT) RESULT(drive)
    TYPE(drive_type)    :: drive
    INTEGER, INTENT(IN) :: UNIT
    READ(UNIT,*) drive
  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
  
  SUBROUTINE write_raw_drive(drive,UNIT)
    TYPE(drive_type), INTENT(IN) :: drive
    INTEGER,          INTENT(IN) :: UNIT
    WRITE(UNIT,*) drive
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  
  FUNCTION new_param(drive)
 
    ! CHANGE THE HAMILTONIAN PARAMETER FROM 'drive%param' TO 'new_param' 
    ! TO DRIVE THE OBSERVABLE 'drive%obs' TOWARDS 'drive%target_obs'
    ! WITHIN 'drive%obs_error_max'
 
    REAL(DBL)                       :: new_param
    TYPE(drive_type), INTENT(INOUT) :: drive

    ! WRITE HEADER
    CALL write_header_drive(drive)

    IF(ABS(drive%obs-drive%target_obs)<=drive%obs_error_max)THEN
      ! ALREADY IN THE DESIRED WINDOW
      new_param = drive%param
    ELSE
      ! FIND PARAMETER DICHOTOMY INTERVAL
      new_param = find_dichotomy_interval(drive)
    ENDIF

    drive%param = new_param ! UPDATE IMPURITY HAMILTONIAN PARAMETER

  END FUNCTION

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  RECURSIVE FUNCTION find_dichotomy_interval(drive) RESULT(new_param)
    REAL(DBL)                       :: new_param
    TYPE(drive_type), INTENT(INOUT) :: drive

    SELECT CASE(drive%bnds2find)
      CASE(2)

        ! NO BOUND WAS FOUND YET: FIND 1 OF 2 obs_bnds 

        new_param = slide_dichotomy_interval(drive)
        drive%bnds2find  = 1
      CASE(1)

        ! 1 BOUND WAS ALREADY FOUND: ADJUST PARAMETER INTERVAL 

        IF(drive%param==drive%param_bnds(2))THEN 
          ! PARAMETER IS AT UPPER BOUND
          IF(drive%obs<drive%target_obs)THEN 
            ! ... BUT WE'RE STILL BELOW TARGET => SLIDE PARAMETER TO RIGHT
            new_param = slide_dichotomy_interval(drive)
          ELSE                           
            ! WE RE ABOVE TARGET => WE HAVE A DICHOTOMY INTERVAL FOR THE OBSERVABLE
            drive%bnds2find = 0
            new_param = find_dichotomy_interval(drive) 
          ENDIF
        ENDIF
        IF(drive%param==drive%param_bnds(1))THEN
          ! PARAMETER IS AT LOWER BOUND
          IF(drive%obs>drive%target_obs)THEN
            ! ... BUT WE'RE STILL ABOVE TARGET => SLIDE PARAMETER INTERVAL TO LEFT
            new_param = slide_dichotomy_interval(drive)
          ELSE
            ! WE'RE BELOW TARGET => WE HAVE A DICHOTOMY INTERVAL FOR THE OBSERVABLE
            drive%bnds2find  = 0
            new_param = find_dichotomy_interval(drive) 
          ENDIF
        ENDIF
      CASE(0)

        ! WE ALREADY HAVE A DICHOTOMY INTERVAL 

        IF(drive%obs<drive%target_obs)THEN
          drive%param_bnds(1) = drive%param
          drive%obs_bnds(1)   = drive%obs
          new_param           = SUM(drive%param_bnds) / 2
        ELSE
          drive%param_bnds(2) = drive%param
          drive%obs_bnds(2)   = drive%obs
          new_param           = SUM(drive%param_bnds) / 2
        ENDIF
    END SELECT
    
    CALL write_footer_drive(new_param,drive)

  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE write_header_drive(drive,UNIT)
    TYPE(drive_type),  INTENT(IN) :: drive
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    INTEGER :: unit_
    unit_ = log_unit
    IF(PRESENT(UNIT)) unit_ = UNIT
    CALL dump_message(UNIT=UNIT,TEXT="######################################")
    CALL dump_message(UNIT=UNIT,TEXT="### START DRIVE "//TRIM(ADJUSTL(drive%obsname))//" ###")
    CALL dump_message(UNIT=UNIT,TEXT="######################################")
    WRITE(unit_,'(3(a,f0.6),a)') &
    "# "//TRIM(ADJUSTL(drive%paramname))//" = ",drive%param," "//TRIM(ADJUSTL(drive%obsname))//" = ",drive%obs," (target = ",drive%target_obs,")"
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE write_footer_drive(new_param,drive,UNIT)
    REAL(DBL),         INTENT(IN) :: new_param
    TYPE(drive_type),  INTENT(IN) :: drive
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    INTEGER :: unit_
    unit_ = log_unit
    IF(PRESENT(UNIT)) unit_ = UNIT
    WRITE(unit_,'(a,f0.6,a,I0,a)') "# new "//TRIM(ADJUSTL(drive%paramname))//" = ",new_param," (",drive%bnds2find," remaining bounds to find)"
    WRITE(unit_,'(2(a,f0.6),a)')   "# "//TRIM(ADJUSTL(drive%paramname))//" dichotomy interval = (",drive%param_bnds(1),",",drive%param_bnds(2),")"
    WRITE(unit_,'(2(a,f0.6),a)')   "# "//TRIM(ADJUSTL(drive%obsname))//" dichotomy interval = (",drive%obs_bnds(1),",",drive%obs_bnds(2),")"
    CALL dump_message(UNIT=UNIT,TEXT="####################################")
    CALL dump_message(UNIT=UNIT,TEXT="### END DRIVE "//TRIM(ADJUSTL(drive%obsname))//" ###")
    CALL dump_message(UNIT=UNIT,TEXT="####################################")
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  FUNCTION slide_dichotomy_interval(drive) RESULT(new_param)
    REAL(DBL)                       :: new_param
    TYPE(drive_type), INTENT(INOUT) :: drive
    IF(drive%obs<drive%target_obs)THEN 
      ! BELOW TARGET => NEW PARAMETER INTERVAL EXTENDS TO THE RIGHT 
      drive%param_bnds  = (/drive%param,drive%param + drive%param_error_max/)
      drive%obs_bnds(1) = drive%obs
      new_param         = drive%param_bnds(2)
    ELSE 
      ! ABOVE TARGET => NEW PARAMETER INTERVAL EXTENDS TO THE LEFT
      drive%param_bnds  = (/drive%param - drive%param_error_max,drive%param/)
      drive%obs_bnds(2) = drive%obs
      new_param         = drive%param_bnds(1)
    ENDIF
  END FUNCTION

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

END MODULE 
