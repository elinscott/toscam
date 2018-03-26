MODULE H_class

   use genvar,       only: DBL
   use sector_class, only: sector_type

   IMPLICIT NONE

   private

   !------------------------------------------!
   ! sector WILL SELECT APPROPRIATE ROUTINES  !
   !------------------------------------------!

   TYPE(sector_type), POINTER :: sector_h => NULL()

   INTERFACE Hmult__
      MODULE PROCEDURE Hmultr, Hmultc
   END INTERFACE

   public :: dimen_H
   public :: sector_H
   public :: title_H_
   public :: hmultc
   public :: hmultr
   public :: Hmult__

contains

   subroutine new_H(AIM, sector_in)

      use sector_class,   only: sector_type
      use haimsz_class,   only: new_haimsz
      use haimupdo_class, only: new_haimupdo
      use aim_class,      only: aim_type

      implicit none

      TYPE(AIM_type), INTENT(IN)            :: AIM
      TYPE(sector_type), INTENT(IN), TARGET :: sector_in

      sector_h => sector_in
      IF     (ASSOCIATED(sector_h%updo))THEN
         CALL new_HAIMupdo(AIM, sector_h%updo)
      ELSE IF(ASSOCIATED(sector_h%sz))  THEN
         CALL new_HAIMsz(  AIM, sector_h%sz)
      ENDIF
   end subroutine

   subroutine delete_H()

      ! IF sector => NULL THEN H WAS NEVER CREATED

      use haimupdo_class, only: delete_haimupdo
      use haimsz_class,   only: delete_haimsz

      implicit none


      IF(ASSOCIATED(sector_h))THEN
         IF     (ASSOCIATED(sector_h%updo))THEN
            CALL delete_HAIMupdo()
         ELSE IF(ASSOCIATED(sector_h%sz))  THEN
            CALL delete_HAIMsz()
         ENDIF
         NULLIFY(sector_h)
      ENDIF
   end subroutine

   subroutine Hmultc(n, vec_out, vec_in)

      use haimsz_class,        only: haimsz_mult, haimsz_mult_fly
      use haimupdo_class,      only: haimupdo_mult, haimupdo_mult_split
      use genvar,              only: dbl
      use globalvar_ed_solver, only: ON_FLY, USE_TRANSPOSE_TRICK_MPI

      implicit none

      COMPLEX(DBL), INTENT(INOUT) :: vec_out(n)
      COMPLEX(DBL), INTENT(IN)    :: vec_in(n)
      INTEGER :: n

      IF(ASSOCIATED(sector_h%updo))THEN
         if(.not.ON_FLY)then
            if(.not.USE_TRANSPOSE_TRICK_MPI)then
               CALL HAIMupdo_mult(vec_out, vec_in)
            else
               CALL HAIMupdo_mult_split(vec_out, vec_in)
            endif
         else
            stop 'ON FLY and up dn basis'
         endif
      ELSE IF(ASSOCIATED(sector_h%sz))THEN
         ! ebl: removing GPU functionality
         ! if(use_cuda_lanczos)then
         !  CALL Hmult_sz_complex_cuda(vec_out, vec_in)
         ! else
         if(.not.ON_FLY)then
            CALL HAIMsz_mult(vec_out, vec_in)
         else
            CALL HAIMsz_mult_fly(vec_out, vec_in)
         endif
         ! endif
      ELSEIF(.true.) then
         STOP 'Hmultc, no array are associated!'
      ENDIF

   end subroutine

   subroutine Hmultr(n, vec_out, vec_in)

      use genvar,              only: dbl
      use globalvar_ed_solver, only: ON_FLY, USE_TRANSPOSE_TRICK_MPI
      use haimsz_class,        only: haimsz_mult, haimsz_mult_fly
      use haimupdo_class,      only: haimupdo_mult, haimupdo_mult_split

      implicit none

      REAL(DBL), INTENT(INOUT) :: vec_out(n)
      REAL(DBL), INTENT(IN)    :: vec_in(n)
      INTEGER :: n

      IF(ASSOCIATED(sector_h%updo))THEN
         if(.not.ON_FLY)then
            if(.not.USE_TRANSPOSE_TRICK_MPI)then
               CALL HAIMupdo_mult(vec_out, vec_in)
            else
               CALL HAIMupdo_mult_split(vec_out, vec_in)
            endif
         else
            stop 'ON FLY and up dn basis'
         endif
      ELSE IF(ASSOCIATED(sector_h%sz))THEN
         ! ebl: remove GPU functionality
         ! if(use_cuda_lanczos)then
         !    CALL Hmult_sz_real_cuda(vec_out, vec_in)
         ! else
         if(.not.ON_FLY)then
            CALL HAIMsz_mult(vec_out, vec_in)
         else
            CALL HAIMsz_mult_fly(vec_out, vec_in)
         endif
         ! endif
      ELSEIF(.true.) then
         STOP 'Hmultr, no array are associated!'
      ENDIF
   end subroutine

   pure function dimen_H()

      use sector_class, only: dimen_func

      implicit none

      INTEGER :: dimen_H

      dimen_H = dimen_func(sector_h)
   end function

   function chunk_H()

      use sector_class, only: chunk_func

      implicit none

      INTEGER :: chunk_H

      chunk_H = chunk_func(sector_h)
   end function

   function istatemin_H()

      use sector_class, only: istatemin__

      implicit none

      INTEGER :: istatemin_H

      istatemin_H = istatemin__(sector_h)
   end function

   function istatemax_H()

      use sector_class, only: istatemax__

      implicit none

      INTEGER :: istatemax_H

      istatemax_H = istatemax__(sector_h)
   end function

   function title_H_()

      use sector_class, only: title_sector

      implicit none

      CHARACTER(LEN = 100) :: title_H_

      title_H_ = 'HAMILTONIAN IN SECTOR ' // &
           TRIM(ADJUSTL(title_sector(sector_h)))
   end function

end module
