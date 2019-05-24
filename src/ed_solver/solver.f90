MODULE solver

   use eigen_sector_class, only: eigensectorlist_type

   IMPLICIT NONE

   private

   TYPE(eigensectorlist_type), PRIVATE, SAVE :: sector2diagH ! sectors to
   ! diagonalize
   TYPE(eigensectorlist_type), PRIVATE, SAVE :: oldGS ! GS from previous run

   LOGICAL, PRIVATE      :: start_from_old_gs = .false.

   public :: init_solver
   public :: remove_filef
   public :: solve_aim
   public :: write_info_solver

contains

#include "solver_init.h"
#include "solver_write.h"

   subroutine solve_AIM(GS, AIM, OLDGSFILE, COMPUTE_ALL_CORREL, skip_vector, &
                        skip_filter)

      ! AIM DIAGONALIZATION

      use eigen_sector_class, only: copy_eigensectorlist, &
         delete_eigensectorlist, eigensectorlist_type, filter_eigensector, &
         gsenergy, nambu_energy_shift, print_eigensectorlist, &
         write_raw_eigensectorlist
      use genvar, only: dbl, iproc, log_unit, no_mpi, rank, size2
      use globalvar_ed_solver, only: ALL_FIRST_CALL, dump_ground_state, &
         dEmax, FLAG_FULL_ED_GREEN, which_lanczos
      use mpi_mod, only: mpibarrier
      use aim_class, only: aim_type

      implicit none

      TYPE(eigensectorlist_type), INTENT(INOUT) :: GS
      TYPE(AIM_type), INTENT(IN)                :: AIM
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: OLDGSFILE
      INTEGER           :: isector
      REAL(DBL)         :: E0
      LOGICAL, SAVE     :: first_call = .true.
      INTEGER           :: uup, ddn, nstu, nstd, itot, ssz
      LOGICAL           :: NOT_COMMENSURATE
      LOGICAL, OPTIONAL :: COMPUTE_ALL_CORREL
      LOGICAL, OPTIONAL :: skip_vector, skip_filter

      itot = AIM%bath%Nb + AIM%impurity%Nc

      if (size2 > 1 .and. .not. no_mpi) then
         write (*, *) 'AIM SOLVER : start SOLVE AIM, RANK = ', rank
         call mpibarrier
      endif

      IF (first_call .AND. start_from_old_gs) THEN
         !--------------------------------------------------!
         ! READ GS FROM PREVIOUS RUN (FIRST ITERATION ONLY) !
         !--------------------------------------------------!
         write (*, *) ' FIRST_CALL AND STARTING FROM OLD GS - RANK ', rank
         write (log_unit, *) '- READ GS FROM PREVIOUS RUN -'
         call dumpmess1
         CALL copy_eigensectorlist(GS, oldGS)
         CALL delete_eigensectorlist(oldGS)
      ELSE
         !-------------!
         ! DIAGONALIZE !
         !-------------!
         write (*, *) ' AIM SOLVER : NOT STARTING FROM OLD GS - RANK ', rank
         call dumpmess2
         write (log_unit, *) '- INITIALIZE GS TO LIST OF EMPTY SECTORS -'
         CALL copy_eigensectorlist(GS, sector2diagH)
         CALL get_eigen_energies
         if (AIM%bath%SUPER) call nambu_energy_shift(GS)
         E0 = GSenergy(GS)
         if (.not. FLAG_FULL_ED_GREEN .and. .not. present(skip_filter)) CALL &
            filter_eigensector(GS, [E0, E0 + dEmax])
         CALL print_eigensectorlist(GS)

         if (AIM%bath%SUPER) call nambu_energy_shift(GS, back=.true.)

         if (.not. (which_lanczos == 'FULL_ED' .or. which_lanczos == 'ARPACK') &
             .and. .not. present(skip_vector)) call get_GS

      ENDIF

      call get_new_bounds

      IF (present(OLDGSFILE) .and. iproc == 1 .AND. dump_ground_state .AND. &
          ((.NOT. first_call) .OR. (.NOT. start_from_old_gs))) then
         write (*, *) 'WRITING RAW EIGENSECTOR LIST'
         CALL write_raw_eigensectorlist(GS, OLDGSFILE)
      endif
      IF (first_call .and. .not. ALL_FIRST_CALL) first_call = .false.
   contains

#include "solver.h"

   end subroutine

   subroutine remove_filef(filename1)

      use genvar, only: status

      implicit none

      character*(*) :: filename1
      integer       :: unit_, ios
      logical       :: is_it_opened, checkfile

      write (*, *) 'ERASING FILE : ', trim(adjustl(filename1))
      inquire (file=trim(adjustl(filename1)), exist=checkfile)
      if (.not. checkfile) return
      unit_ = 20
      do
         unit_ = unit_ + 1
         INQUIRE (unit=unit_, OPENED=is_it_opened, iostat=ios)
         if (.not. is_it_opened .and. ios == 0) exit
      enddo
      open (unit=unit_, file=trim(adjustl(filename1)))
      close (unit_, status='delete')
   end subroutine

end module

