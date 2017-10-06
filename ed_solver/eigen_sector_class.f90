MODULE eigen_sector_class

  USE sector_class
  USE eigen_class

  IMPLICIT NONE

  REAL(DBL),    PARAMETER, PRIVATE       :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,      PARAMETER, PRIVATE       :: F=.FALSE.,T=.TRUE.

! PAIR (SECTOR,LOWEST EIGENSTATES)
  
  TYPE eigensector_type
     TYPE(sector_type)    :: sector
     TYPE(eigenlist_type) :: lowest
  END TYPE


  INTERFACE new_eigensector
     MODULE PROCEDURE copy_eigensector
     MODULE PROCEDURE new_eigensector_from_scratch
  END INTERFACE


  TYPE eigensectorlist_type
!CEDRIC
     INTEGER                         :: nsector=0
     TYPE(eigensector_type), POINTER :: es(:) => NULL()
  END TYPE


  INTERFACE new_rcvec_from_file
     MODULE PROCEDURE new_rcvec_from_file_c,new_rcvec_from_file_r
  END INTERFACE


CONTAINS


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine new_rcvec_from_file_r(unit,vec,sector)
  integer           :: unit,istate,jstate,stateup,statedo
  TYPE(sector_type) :: sector
  real(8)           :: coeff,vec(:)
    DO istate=1,size(vec)
      READ(unit,*)stateup,statedo,coeff
      jstate = stateup + IBITS(NOT(statedo),0,12)*(2**12)
      vec(rank_func(jstate,sector)) = coeff
    ENDDO
  end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine new_rcvec_from_file_c(unit,vec,sector)
  integer           :: unit,istate,jstate,stateup,statedo
  TYPE(sector_type) :: sector
  complex(8)        :: coeff,vec(:)
    DO istate=1,size(vec)
      READ(unit,*)stateup,statedo,coeff
      jstate = stateup + IBITS(NOT(statedo),0,12)*(2**12)
      vec(rank_func(jstate,sector)) = coeff
    ENDDO
  end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_eigensectorlist(list,nsector)
    TYPE(eigensectorlist_type), INTENT(INOUT) :: list
    INTEGER,                    INTENT(IN)    :: nsector

    if(messages3) write(log_unit,*) 'delete eigensector'
    CALL delete_eigensectorlist(list)
    list%nsector = nsector
    if(messages3) write(log_unit,*) 'new eigensector list nsector ---> ' , nsector

    IF(list%nsector>0) then
      if(associated(list%es)) deallocate(list%es)
      if(messages3) write(log_unit,*) 'allocate new list....'
      ALLOCATE(list%es(nsector))
    ELSE
      if(messages3) write(log_unit,*) 'create list with 0 eigensector'
      if(associated(list%es)) deallocate(list%es)
      list%es=>NULL() 
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

  SUBROUTINE remove_eigensector(sector,list)
    TYPE(eigensectorlist_type) :: list
    TYPE(sector_type)          :: sector
    INTEGER                    :: isector,jsector,truerank
    TYPE(eigensectorlist_type) :: tmp

    truerank = rank_sector_in_list(sector,list)
    IF(truerank==0) then  
      write(*,*) "ERROR IN remove_eigensector: SECTOR "//TRIM(title_sector(sector))//" ISNT IN LIST!" 
      stop
    endif
    CALL copy_eigensectorlist(tmp,list)
    CALL  new_eigensectorlist(list,list%nsector-1)
    jsector = 0
    DO isector=1,tmp%nsector
      IF(isector/=truerank)THEN
        jsector = jsector + 1
        CALL copy_eigensector(list%es(jsector),tmp%es(isector))
      ENDIF
    ENDDO
    CALL delete_eigensectorlist(tmp)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE add_eigensector(es,list)
    TYPE(eigensectorlist_type), INTENT(INOUT) :: list
    TYPE(eigensector_type),     INTENT(IN)    :: es
    INTEGER                                   :: isector
    TYPE(eigensectorlist_type)                :: tmp

    CALL copy_eigensectorlist(tmp,list)
    CALL  new_eigensectorlist(list,list%nsector+1)
    DO isector=1,tmp%nsector
      CALL copy_eigensector(list%es(isector),tmp%es(isector))
    ENDDO
    CALL copy_eigensector(list%es(list%nsector),es)
    CALL delete_eigensectorlist(tmp)

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_eigensector_from_scratch(es,sector_in)
    !------------------------!
    ! NEW EMPTY EIGENSECTOR  !
    !------------------------!
    TYPE(eigensector_type), INTENT(INOUT) :: es
    TYPE(sector_type),      INTENT(IN)    :: sector_in
    CALL delete_eigensector(es)
    CALL copy_sector(es%sector,sector_in)
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE copy_eigensector(esout,esin) 
    TYPE(eigensector_type), INTENT(INOUT) :: esout
    TYPE(eigensector_type), INTENT(IN)    :: esin
    CALL delete_eigensector(esout)
    CALL copy_sector(esout%sector,esin%sector)
    CALL copy_eigenlist(esout%lowest,esin%lowest)
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE copy_eigensectorlist(listout,listin)
    TYPE(eigensectorlist_type), INTENT(INOUT) :: listout
    TYPE(eigensectorlist_type), INTENT(IN)    :: listin
    INTEGER                                   :: isector

    CALL new_eigensectorlist(listout,listin%nsector)
    if(listin%nsector==0) return
    DO isector=1,listin%nsector
      CALL copy_eigensector(listout%es(isector),listin%es(isector))
    ENDDO

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_eigensector(es)
    TYPE(eigensector_type), INTENT(INOUT) :: es
    CALL delete_sector(es%sector) 
    CALL delete_eigenlist(es%lowest)
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_eigensectorlist(list)
    TYPE(eigensectorlist_type), INTENT(INOUT) :: list
    INTEGER                                   :: isector

    DO isector=1,list%nsector
      CALL delete_eigensector(list%es(isector))
    ENDDO
    if(messages3) write(log_unit,*) 'delete list'
    IF(ASSOCIATED(list%es)) DEALLOCATE(list%es)
    list%nsector = 0
    list%es=>NULL()
    if(messages3) write(log_unit,*) '---done---'

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE filter_eigensector(list,window) 
    TYPE(eigensectorlist_type),INTENT(IN) :: list
    REAL(DBL),                 INTENT(IN) :: window(2)
    INTEGER                               :: isector,nn

  !-----------------------------------------------------------------------!
  ! FILTER OUT ALL EIGENSTATES OUTSIDE window IN ALL EIGENSECTORS OF list !
  !-----------------------------------------------------------------------!

    isector = 1
    nn      = size(list%es(:))

    DO 
      if(messages4)then
       write(log_unit,*) '.....start to filter sector    : ', isector,window
       write(log_unit,*) '.....ISECTOR, NN, SIZE OF LIST : ', isector,nn,size(list%es(:)),list%nsector
      endif
      if(isector>nn) exit
      if(list%es(isector)%lowest%neigen>0) CALL filter_eigen(list%es(isector)%lowest,window)
      IF(list%es(isector)%lowest%neigen==0) then
        if(messages4) write(log_unit,*) '----remove empty sector----', isector
        CALL remove_eigensector(list%es(isector)%sector,list) 
        if(messages4) write(log_unit,*) 'done, now remains .. sectors : ' , nn-1
        nn=nn-1
        if(nn/=list%nsector) stop 'error in filter sector, nn different from list nsector'
      else
        isector=isector+1
      endif
    ENDDO

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE nambu_energy_shift(list,back)
   TYPE(eigensectorlist_type) :: list
   INTEGER                    :: isector,nn,k,kk
   logical,optional           :: back

   nn = size(list%es(:))

   DO isector=1,nn

    if(ASSOCIATED(list%es(isector)%lowest%eigen)) then

     kk = size( list%es(isector)%lowest%eigen(:) )

     do k=1,kk
       if(.not.present(back))then
        list%es(isector)%lowest%eigen(k)%val = list%es(isector)%lowest%eigen(k)%val + energy_global_shift + energy_global_shift2
       else
        list%es(isector)%lowest%eigen(k)%val = list%es(isector)%lowest%eigen(k)%val - energy_global_shift - energy_global_shift2
       endif
     enddo

    else

     write(*,*) 'SECTOR IN NAMBU NOT ASSOCIATED' 
     write(*,*) 'SECTOR : ' , isector , '  /  ' , nn
      
     if(strongstop) stop 'critical'

    endif

   ENDDO

   write(log_unit,*) ' --- the ground state energy was shifted by --- '
   write(log_unit,*) ' ---           shift due to mu              --- ' 
   write(log_unit,*)               energy_global_shift 
   write(log_unit,*) ' ---       shift due to epsilons bath       --- ' 
   write(log_unit,*)               energy_global_shift2

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION GSenergy(list) 
 
    ! FIND THE LOWEST EIGENVALUE IN list
 
    REAL(DBL)                              :: GSenergy
    TYPE(eigensectorlist_type), INTENT(IN) :: list
    INTEGER                                :: isector
    REAL(DBL)                              :: E0

    GSenergy = huge_
    DO isector=1,list%nsector
      E0 = min_eigen(list%es(isector)%lowest) ! GS ENERGY OF THE LOCAL SECTOR
      IF(E0<GSenergy) GSenergy = E0
    ENDDO

  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  logical function not_commensurate_sector(list,isector)
    TYPE(eigensectorlist_type), INTENT(IN) :: list
    INTEGER                                :: isector

     if(.not.associated(list%es(isector)%sector%updo))then
      not_commensurate_sector=.true.
      return
     endif

      not_commensurate_sector=ANY(list%es(isector)%sector%updo%up%chunk/= &
                                  list%es(isector)%sector%updo%up%chunk(1)) .or. &
                              ANY(list%es(isector)%sector%updo%down%chunk/=&
                                  list%es(isector)%sector%updo%down%chunk(1)) .or. &
                              mod(list%es(isector)%sector%updo%up%chunk(1),size2)/=0
  end function

                   !-----------------------------------------!

  logical function not_commensurate_sector_(sector)
    TYPE(sector_type), INTENT(IN) :: sector

     if(.not.associated(sector%updo))then
      not_commensurate_sector_=.true.
      return
     endif

      not_commensurate_sector_=ANY(sector%updo%up%chunk/=sector%updo%up%chunk(1)) .or.  &
                               ANY(sector%updo%down%chunk/=sector%updo%down%chunk(1)) .or. &
                               mod(sector%updo%up%chunk(1),size2)/=0
  end function


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  FUNCTION partition(beta,list) RESULT(Zpart)
    REAL(DBL)                              :: Zpart
    REAL(DBL),                  INTENT(IN) :: beta
    TYPE(eigensectorlist_type), INTENT(IN) :: list 
    INTEGER                                :: isector 
    REAL(DBL)                              :: E0 

    E0    =  GSenergy(list)  
    Zpart =  zero
    DO isector=1,list%nsector
      Zpart = Zpart + sum_boltz(list%es(isector)%lowest,beta,E0) 
    ENDDO

  END FUNCTION

          !---------------------------------------------!

  FUNCTION average_energy(beta,list) RESULT(EE)
    REAL(DBL)                              :: Zpart,EE
    REAL(DBL),                  INTENT(IN) :: beta
    TYPE(eigensectorlist_type), INTENT(IN) :: list
    INTEGER                                :: isector
    REAL(DBL)                              :: E0

    E0    =  GSenergy(list)
    Zpart =  zero
    EE    =  zero
    DO isector=1,list%nsector
     Zpart = Zpart +         sum_boltz(list%es(isector)%lowest,beta,E0)
     EE    = EE    + sum_boltz_times_E(list%es(isector)%lowest,beta,E0)
    ENDDO
    EE=EE/Zpart

  END FUNCTION

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION nlowest(list) 
    INTEGER                                :: nlowest
    TYPE(eigensectorlist_type), INTENT(IN) :: list 
    INTEGER                                :: isector 
    nlowest = 0
    DO isector=1,list%nsector
      nlowest = nlowest + list%es(isector)%lowest%neigen
    ENDDO 
  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE print_eigensectorlist(list,UNIT) 
    TYPE(eigensectorlist_type), INTENT(IN) :: list
    INTEGER, OPTIONAL,          INTENT(IN) :: UNIT
    INTEGER                                :: isector
    DO isector=1,list%nsector
      CALL dump_message(UNIT=UNIT,TEXT=&
      "# SECTOR "//TRIM(title_sector(list%es(isector)%sector))//" [dim="//c2s(i2c(dimen_func(list%es(isector)%sector)))//"]")
      CALL print_eigenlist(list%es(isector)%lowest,UNIT=UNIT)
    ENDDO
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE write_raw_eigensectorlist(list,FILEOUT) 
    TYPE(eigensectorlist_type), INTENT(IN) :: list
    CHARACTER(LEN=*),           INTENT(IN) :: FILEOUT
    INTEGER                                :: isector
    INTEGER                                :: UNIT 
    CALL open_safe(UNIT,FILEOUT,'UNKNOWN','WRITE',get_unit=.true.)
    DO isector=1,list%nsector
      WRITE(UNIT,*) T ! there is a new item
      CALL write_raw_sector   (list%es(isector)%sector,UNIT)
      CALL write_raw_eigenlist(list%es(isector)%lowest,UNIT)
    ENDDO
    WRITE(UNIT,*) F ! last item
    CALL close_safe(UNIT)
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE read_raw_eigensectorlist(list,FILEIN) 
   TYPE(eigensectorlist_type), INTENT(INOUT) :: list
   CHARACTER(LEN=*),           INTENT(IN)    :: FILEIN
   TYPE(eigensector_type)                    :: es
   INTEGER                                   :: UNIT
   LOGICAL                                   :: new_item

    CALL open_safe(UNIT,FILEIN,'UNKNOWN','READ',get_unit=.true.)

    DO
      READ(UNIT,*) new_item
      IF(.NOT.new_item) EXIT
      CALL read_raw_sector(es%sector,UNIT)
      CALL read_raw_eigenlist(es%lowest,UNIT)
      CALL add_eigensector(es,list)
    ENDDO
    CALL close_safe(UNIT)

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION rank_sector_in_list(sector,list)
    INTEGER                                :: rank_sector_in_list
    TYPE(sector_type),          INTENT(IN) :: sector
    TYPE(eigensectorlist_type), INTENT(IN) :: list
    INTEGER :: isector
    rank_sector_in_list = 0
    DO isector=1,list%nsector
      IF(equal_sector(sector,list%es(isector)%sector))THEN
        rank_sector_in_list = isector
        EXIT 
      ENDIF
    ENDDO
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
