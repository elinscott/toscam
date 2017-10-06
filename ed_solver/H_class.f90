MODULE H_class

  use apply_P 

  IMPLICIT NONE

  !------------------------------------------!  
  ! sector WILL SELECT APPROPRIATE ROUTINES  !
  !------------------------------------------!

  TYPE(sector_type), POINTER :: sector_h => NULL()

  INTERFACE Hmult__
    MODULE PROCEDURE Hmultr,Hmultc
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

  SUBROUTINE new_H(AIM,sector_in)
    TYPE(AIM_type),    INTENT(IN)         :: AIM
    TYPE(sector_type), INTENT(IN), TARGET :: sector_in
    sector_h => sector_in
    IF     (ASSOCIATED(sector_h%updo))THEN
      CALL new_HAIMupdo(AIM,sector_h%updo)
    ELSE IF(ASSOCIATED(sector_h%sz))  THEN
      CALL new_HAIMsz(  AIM,sector_h%sz)
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

  SUBROUTINE delete_H
    ! IF sector=>NULL THEN H WAS NEVER CREATED!
    IF(ASSOCIATED(sector_h))THEN
      IF     (ASSOCIATED(sector_h%updo))THEN
       CALL delete_HAIMupdo()
      ELSE IF(ASSOCIATED(sector_h%sz))  THEN
       CALL delete_HAIMsz()
      ENDIF
      NULLIFY(sector_h)
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

  SUBROUTINE Hmultc(n,vec_out,vec_in)
    INTEGER                     ::  n
    COMPLEX(DBL), INTENT(INOUT) ::  vec_out(n)
    COMPLEX(DBL), INTENT(IN)    ::  vec_in(n)

    IF(ASSOCIATED(sector_h%updo))THEN
     if(.not.ON_FLY)then
      if(.not.USE_TRANSPOSE_TRICK_MPI)then
       CALL HAIMupdo_mult(vec_out,vec_in)
      else
       CALL HAIMupdo_mult_split(vec_out,vec_in)
      endif
     else
      stop 'ON FLY and up dn basis'
     endif
    ELSE IF(ASSOCIATED(sector_h%sz))THEN
     if(use_cuda_lanczos)then
      CALL Hmult_sz_complex_cuda(vec_out,vec_in)
     else
      if(.not.ON_FLY)then
        CALL HAIMsz_mult(vec_out,vec_in) 
      else
        CALL HAIMsz_mult_fly(vec_out,vec_in) 
      endif
     endif
    ELSEIF(.true.) then
      STOP 'Hmultc, no array are associated!'
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

  SUBROUTINE Hmultr(n,vec_out,vec_in)  
    INTEGER                     :: n
    REAL(DBL),    INTENT(INOUT) :: vec_out(n)
    REAL(DBL),    INTENT(IN)    :: vec_in(n)
    IF(ASSOCIATED(sector_h%updo))THEN
     if(.not.ON_FLY)then
      if(.not.USE_TRANSPOSE_TRICK_MPI)then
       CALL HAIMupdo_mult(vec_out,vec_in)
      else
       CALL HAIMupdo_mult_split(vec_out,vec_in)
      endif
     else
      stop 'ON FLY and up dn basis'
     endif
    ELSE IF(ASSOCIATED(sector_h%sz))THEN
     if(use_cuda_lanczos)then
       CALL Hmult_sz_real_cuda(vec_out,vec_in)
      else
       if(.not.ON_FLY)then
         CALL HAIMsz_mult(vec_out,vec_in)
       else
         CALL HAIMsz_mult_fly(vec_out,vec_in)
       endif
      endif
    ELSEIF(.true.) then
      STOP 'Hmultr, no array are associated!'
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

     !-------------------------------!

  pure FUNCTION dimen_H() 
    INTEGER :: dimen_H
    dimen_H = dimen_func(sector_h)
  END FUNCTION 

     !-------------------------------!

  FUNCTION chunk_H() 
    INTEGER :: chunk_H
    chunk_H = chunk_func(sector_h)
  END FUNCTION 

     !-------------------------------!
  
  FUNCTION istatemin_H() 
    INTEGER :: istatemin_H
    istatemin_H = istatemin__(sector_h)
  END FUNCTION 

     !-------------------------------!
  
  FUNCTION istatemax_H() 
    INTEGER :: istatemax_H
    istatemax_H = istatemax__(sector_h)
  END FUNCTION 

     !-------------------------------!
  
  FUNCTION title_H_() 
    CHARACTER(LEN=100) :: title_H_
    title_H_ = 'HAMILTONIAN IN SECTOR '//TRIM(ADJUSTL(title_sector(sector_h)))
  END FUNCTION 

     !-------------------------------!

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

END MODULE
