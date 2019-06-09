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

   type(timing_info) :: routines(max_number_of_timed_routines)

   public :: start_timer
   public :: stop_timer
   public :: timing_summary

contains

   subroutine start_timer(name)
      ! Starts the timer corresponding to the routine "name"

      use common_def, only: utils_assert

      implicit none

      character(len=*), intent(in) :: name
      integer                      :: i

      i = routine_index(name)

      if (i == 0) then
         ! We have not called this routine before; initialise it
         i = free_routine_index()

         if (i == 1) then
            ! This is the first timer we've ever called; initialise the routines object
            routines(:)%num_calls = 0
            routines(:)%cum_time = 0.0_DP
         end if
            
         routines(i)%name = trim(adjustl(name))
      else
         call utils_assert(.not. routines(i)%running, 'Error in timer_mod: '&
               // trim(adjustl(routines(i)%name)) // ' timer is still running when &
               &attempting to start timer')
      end if

      routines(i)%num_calls = routines(i)%num_calls + 1
      routines(i)%start_time = elapsed_time()
      routines(i)%running = .true.

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
      routines(i)%cum_time = final_time - routines(i)%start_time
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

      do i = 1, max_number_of_timed_routines
         if (trim(adjustl(routines(i)%name)) == trim(adjustl(name))) routine_index = i
      end do

      return

   end function routine_index

   function free_routine_index()
      ! Finds the first routine index that has not been assigned

      use common_def, only: utils_abort

      implicit none
      
      integer           :: free_routine_index
      integer           :: i

      do free_routine_index = 1, max_number_of_timed_routines
         if (routines(i)%num_calls == 0) return
      end do

      call utils_abort("Error in timer_mod: no free routine indices available")

   end function free_routine_index

   subroutine timing_summary(logfile)

      use common_def, only: utils_assert

      implicit none

      integer, optional :: logfile

      integer :: i

      ! ebl: header
      if (present(logfile)) then
         write(logfile, '(a)') repeat("=", 80)
         write(logfile, '(a50,a10,a10)') "Routine", "Calls", "Walltime"
      else
         write(*, '(a)') repeat("=", 80)
         write(*, '(a50,a10,a10)') "Routine", "Calls", "Walltime"
      end if

      ! ebl: loop over routines
      do i = 1, max_number_of_timed_routines

         ! ebl: break out of loop if we reach the end of the stored number of
         ! routines
         if (routines(i)%num_calls > 1) exit

         ! ebl: check all timers have been stopped
         call utils_assert(.not. routines(i)%running, 'Error in timer_mod: '&
               // trim(adjustl(routines(i)%name)) // ' timer is still running when &
               &attempting to print summary')

         ! ebl: print summary
         if (present(logfile)) then
            write(logfile, '(a50,i10,f10.3)') trim(adjustl(routines(i)%name)), &
                  routines(i)%num_calls, routines(i)%cum_time
         else
            write(*, '(a50,i10,f10.3)') trim(adjustl(routines(i)%name)), &
                  routines(i)%num_calls, routines(i)%cum_time
         end if

      end do

      ! ebl: footer
      if (present(logfile)) then
         write(logfile, '(a60,f10.3)') 'TOTAL', elapsed_time()
         write(logfile, '(a)') repeat("=", 80)
      else
         write(*, '(a60,f10.3)') 'TOTAL', elapsed_time()
         write(*, '(a)') repeat("=", 80)
      end if

   end subroutine timing_summary

   function elapsed_time()
      ! Fetch the time elapsed since the start of the calculation
      integer, parameter :: SP = kind(1.0)
      real(kind=SP)      :: delta, time_array(2), etime
      real(kind=DP)      :: elapsed_time
      
      delta = etime(time_array)

      elapsed_time = real(delta, DP)

   end function elapsed_time

end module timer_mod
