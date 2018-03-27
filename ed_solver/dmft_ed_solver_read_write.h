#ifdef NO_SYS_CALL
   subroutine remove_filef(filename1)

      implicit none

      character*(*) :: filename1
      integer       :: unit_, ios
      logical       :: is_it_opened, checkfile

      write(*, *) 'ERASING FILE : ', trim(adjustl(filename1))
      inquire(file = trim(adjustl(filename1)), exist = checkfile)
      if(.not.checkfile) return
      unit_ = 20
      do
         unit_ = unit_ + 1
         INQUIRE(unit = unit_, OPENED = is_it_opened, iostat = ios)
         if(.not.is_it_opened .and. ios == 0)exit
      enddo
      open(unit = unit_, file = trim(adjustl(filename1)))
      close(unit_, status = 'delete')
   end subroutine
#endif

   subroutine write_header(AIM, UNIT)

      use aim_class,      only: aim_type
      use common_def,     only: c2s, dump_message, i2c
      use correlations,   only: write_info_correlations
      use genvar,         only: iproc, nproc, procname
      use impurity_class, only: write_impurity
      use solver,         only: write_info_solver

      implicit none

      TYPE(AIM_type), INTENT(IN)    :: AIM
      INTEGER, OPTIONAL, INTENT(IN) :: UNIT

      CALL dump_message(UNIT = UNIT, TEXT = "############")
      CALL dump_message(UNIT = UNIT, TEXT = "############")
      CALL dump_message(UNIT = UNIT, TEXT = "### DMFT ###")
      CALL dump_message(UNIT = UNIT, TEXT = "############")
      CALL dump_message(UNIT = UNIT, TEXT = "############")

      CALL dump_message(UNIT = UNIT, TEXT = "# THIS IS PROCESS " // &
           c2s(i2c(iproc)) // " OF " // c2s(i2c(nproc)) // " RUNNING ON " // &
           TRIM(ADJUSTL(procname)))

#ifdef _complex
      CALL dump_message(UNIT = UNIT, TEXT = "# FLAG: COMPLEX HAMILTONIAN")
#else
      CALL dump_message(UNIT = UNIT, TEXT = "# FLAG: REAL HAMILTONIAN")
#endif

      CALL write_info_DMFT_loop(UNIT = UNIT)
      CALL write_impurity(AIM%impurity, UNIT = UNIT)
      CALL write_info_correlations(UNIT = UNIT)
      CALL write_info_solver(UNIT = UNIT)

   end subroutine

   subroutine write_begin_iter(UNIT)

      use genvar,              only: log_unit
      use globalvar_ed_solver, only: iterdmft
      use common_def,          only: dump_message

      implicit none

      INTEGER, INTENT(IN), OPTIONAL :: UNIT
      INTEGER :: unit_

      unit_ = log_unit
      IF(PRESENT(UNIT))  unit_ = UNIT

      CALL dump_message(UNIT = UNIT, TEXT = "#")
      CALL dump_message(UNIT = UNIT, TEXT = "#########################")
      CALL dump_message(UNIT = UNIT, TEXT = "#########################")
      WRITE(unit_, '(a, I2, a)') "### DMFT iteration ", iterdmft, " ###"
      CALL dump_message(UNIT = UNIT, TEXT = "#########################")
      CALL dump_message(UNIT = UNIT, TEXT = "#########################")
      CALL dump_message(UNIT = UNIT, TEXT = "#")
   end subroutine

   subroutine read_DMFT_parameters(nmatsu_frequ, FLAGIMP, nstep, Nc)

      use fermion_Hilbert_class, only: FLAG_FORBID_SPLITTING
      use frequency_class,       only: init_frequency_module
      use genvar
      use globalvar_ed_solver
      use lockmod,               only: MAXT
      use mpirout,               only: mpibarrier
      use namelistmod,           only: look_for_command_line_argument, &
           look_for_namelist_in_file, namelist_init, namelist_set
      use stringmanip,           only: tostring

      implicit none

      integer            :: ff, ii, nmatsu_frequ, FLAGIMP, nstep, Nc
      type(namelist_set) :: nm
      logical            :: target
      real(8)            :: targetobs, obserror, paramerror
      character(20)      :: fileb, filec

      ii         =   FLAGIMP
      niterdmft  =   nstep
      fileb      =  'ed_hybrid' // TRIM(ADJUSTL(toString(ii)))
      filec      =  'ed_correl' // TRIM(ADJUSTL(toString(ii)))

      call namelist_init(nm, 400, name_of_namelist = 'DMFT_ED_SOLVER main &
           &input')
#include "dmft_ed_solver_variables.h"
      call look_for_command_line_argument(nm) ! to define input files
      call look_for_namelist_in_file(nm, EDfile)
      call look_for_command_line_argument(nm)
#include "dmft_ed_solver_flags.h"

      if(size2 > 1 .and. .not.no_mpi) then
         write(*, *) 'DONE WITH READING INPUT FILES, RANK = ', rank
         call mpibarrier
      endif
   end subroutine

   subroutine write_info_DMFT_loop(UNIT)

      use genvar,              only: log_unit
      use globalvar_ed_solver, only: niterdmft

      implicit none

      INTEGER, OPTIONAL, INTENT(IN) :: UNIT
      INTEGER              :: unit_, rank
      CHARACTER(LEN = 100) :: fmt_archive


      unit_ = log_unit
      IF(PRESENT(UNIT)) unit_ = UNIT
      WRITE(unit_, '(a, I0, a)') "# WILL PERFORM UP TO ", niterdmft, " DMFT &
           &ITERATIONS"
      CALL flush(unit_)

   end subroutine

   subroutine write_footer(AIM, GS, UNIT)

      use aim_class,          only: aim_type
      use common_def,         only: dump_message
      use correlations,       only: write_static
      use bath_class,         only: write_bath
      use eigen_sector_class, only: eigensectorlist_type
      use genvar,             only: log_unit

      implicit none

      TYPE(AIM_type), INTENT(IN)             :: AIM
      TYPE(eigensectorlist_type), INTENT(IN) :: GS
      INTEGER, OPTIONAL, INTENT(IN)          :: UNIT
      INTEGER :: unit_

      unit_ = log_unit
      IF(PRESENT(UNIT)) unit_ = UNIT
      CALL dump_message(UNIT = unit_, TEXT = "##############################")
      CALL dump_message(UNIT = unit_, TEXT = "##############################")
      CALL dump_message(UNIT = unit_, TEXT = "### END OF DMFT ITERATIONS ###")
      CALL dump_message(UNIT = unit_, TEXT = "##############################")
      CALL dump_message(UNIT = unit_, TEXT = "##############################")
      CALL write_bath(AIM%bath, UNIT = unit_)
      CALL write_static(AIM, GS, UNIT = unit_)

   end subroutine
