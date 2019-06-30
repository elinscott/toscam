module timer_mod
   ! Basic module for timing how long we spend in given routines
   ! Not yet compatible with MPI or OpenMP

   use genvar, only: DP

   implicit none

   private

   integer, parameter :: max_number_of_timed_routines = 20

   type timing_info
      character(len=:), allocatable :: name
      integer                       :: num_calls
      real(kind = DP)               :: cum_time
      real(kind = DP)               :: start_time
      logical                       :: running
   end type timing_info

   integer :: num_entries = 0
   logical :: timer_mod_initialized = .false.
   integer :: program_start_time
   integer :: clock_rate

   integer :: timer_log_unit

   type(timing_info) :: routines(max_number_of_timed_routines)

   public :: start_timer
   public :: stop_timer
   public :: initialize_timing
   public :: finalize_timing

contains

   subroutine initialize_timing()
      ! Call this at the start of a program
   
      use common_def, only: utils_unit
      
      implicit none

      timer_log_unit = utils_unit()
      open(unit = timer_log_unit, file = "timing_calls.log", status = "replace")
      
      call system_clock(count=program_start_time, count_rate=clock_rate)
      timer_mod_initialized = .true.

   end subroutine initialize_timing

   function percentage_accounted_for()
      implicit none

      real(kind=DP) :: percentage_accounted_for

      integer :: i

      percentage_accounted_for = 0

      do i = 1, num_entries
         percentage_accounted_for = percentage_accounted_for + routines(i)%cum_time
      end do

      percentage_accounted_for = percentage_accounted_for/elapsed_time()*100

   end function

   integer function num_running()
      implicit none

      integer :: i

      num_running = 0

      do i = 1, num_entries
         if (routines(i)%running) num_running = num_running + 1
      end do

      return

   end function num_running

   subroutine start_timer(name)
      ! Starts the timer corresponding to the routine "name"

      use common_def, only: utils_assert

      implicit none

      character(len=*), intent(in) :: name
      integer                      :: i

      i = routine_index(name)

      if (i == 0) then
         ! We have not called this routine before; initialise it
         i = num_entries + 1

         if (num_entries == 0) then
            ! This is the first timer we've ever called; initialise the routines object
            routines(:)%num_calls = 0
            routines(:)%cum_time = 0.0_DP
         end if
         
         num_entries = i
         routines(i)%name = trim(adjustl(name))
      else
         call utils_assert(.not. routines(i)%running, 'Error in timer_mod: '&
               // trim(adjustl(routines(i)%name)) // ' timer is still running when &
               &attempting to start timer')
      end if

      routines(i)%num_calls = routines(i)%num_calls + 1
      routines(i)%start_time = elapsed_time()
      routines(i)%running = .true.

      write(timer_log_unit, '(3a,f9.2)') repeat(" ", num_running()), "Starting ", name, percentage_accounted_for()

   end subroutine start_timer

   subroutine stop_timer(name)
      ! Stops the timer corresponding to the routine "name"

      use common_def, only: utils_assert

      implicit none

      character(len=*), intent(in) :: name

      real(kind = DP)         :: final_time
      integer                 :: i

      ! Find index corresponding to timer
      i = routine_index(name)
      call utils_assert(i > 0, 'Error in timer_mod: ' // trim(adjustl(name)) &
            // ' not found in list of timed routines when trying to stop timer')

      ! ebl: stop the timer and record the elapsed time
      final_time = elapsed_time()
      routines(i)%cum_time = routines(i)%cum_time + (final_time - routines(i)%start_time)
      routines(i)%running = .false.

   end subroutine stop_timer

   function routine_index(name)
      ! Finds the index of the "name" routine. Returns zero if the routine is
      ! not in the list of stored routines.

      implicit none

      character(len=*), intent(in) :: name
      integer                      :: routine_index
      integer                      :: i

      routine_index = 0

      do i = 1, num_entries
         if (trim(adjustl(routines(i)%name)) == trim(adjustl(name))) routine_index = i
      end do

      return

   end function routine_index

   subroutine finalize_timing(logfile)

      use common_def, only: utils_assert

      implicit none

      integer, optional :: logfile

      integer :: i
      integer :: loc_logfile
      character(len=42) :: name_string
      character(len=14) :: calls_string
      character(len=14) :: time_string

      ! Make sure timing was initialised
      call utils_assert(timer_mod_initialized, "Error in timer_mod: initialize_timing &
            &needs to be called at the start of the program")

      name_string = "Routine"
      calls_string = "Calls"
      time_string = "Walltime"

      ! By default, write to screen
      loc_logfile = 6
      if (present(logfile)) loc_logfile = logfile

      ! ebl: header
      write(loc_logfile, '(a)') " " // repeat("=", 31) // " TIMING SUMMARY " &
            // repeat("=", 31) // " "
      write(loc_logfile, '(a2,a42,a3,a14,a3,a14)') "  ", adjustl(name_string), " | ", &
            adjustl(calls_string), " | ", adjustl(time_string)

      ! ebl: loop over routines
      do i = 1, num_entries

         ! ebl: check all timers have been stopped
         call utils_assert(.not. routines(i)%running, 'Error in timer_mod: '&
               // trim(adjustl(routines(i)%name)) // ' timer is still running when &
               &attempting to print summary')

         ! ebl: print summary; by using "name_string", the names will be left-aligned
         name_string = routines(i)%name
         write(loc_logfile, '(a2,a42,a3,i14,a3,f14.2)') "  ", adjustl(name_string), &
               " | ", routines(i)%num_calls, " | ", routines(i)%cum_time

      end do

      ! ebl: footer
      write(loc_logfile, '(a61,a3,f14.2)') 'TOTAL', " | ", elapsed_time()
      write(loc_logfile, '(a,a78)') " ", repeat("=", 78)

      close(timer_log_unit)

   end subroutine finalize_timing

   real(kind = DP) function elapsed_time()

      ! Fetch the time elapsed since the start of the calculation
      ! Currently written as a function so that we will be able to wrap MPI_WTIME etc
      ! when parallelism is implemented

      implicit none

      integer :: current_time
      integer :: count_rate

      call system_clock(count = current_time, count_rate = count_rate)

      elapsed_time = real(current_time - program_start_time, DP)/real(count_rate, DP)

   end function elapsed_time

end module timer_mod
