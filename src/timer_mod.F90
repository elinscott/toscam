! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
#ifdef ACCELRYS
!=============================================================================!
! This module contains a dummy interface to a timer routine which could be    !
! used in future to provide a breakdown of the time spent in various routines !
! in the code.                                                                !
!-----------------------------------------------------------------------------!
! This version by Peter Haynes, 2 March 2005                                  !
!=============================================================================!

module timer

  use constants
  use rundat, only : timings_level

  implicit none

  private

  public :: timer_clock

  real(kind=DP) :: start_time

contains

  subroutine timer_clock(name,option,work)

    use comms

    implicit none

    ! Arguments
    character(len=*), intent(in) :: name
    integer, intent(in) :: option
    real(kind=DP), optional, intent(in) :: work

    ! Local variables
    real(kind=DP) :: tot_time

    if (timings_level == 0 .and. (option == 1 .or. option == 2)) return

    select case (option)

    case (0)

       start_time = wrappers_etime()

    case (1:2)

       ! Do nothing

    case (3)

       tot_time = wrappers_etime() - start_time
       call comms_reduce('MAX',tot_time)

       if (pub_on_root) write(stdout,'(/a,f12.3,a,i6,a)') 'TOTAL TIME:', tot_time, &
            's on ', pub_total_num_nodes, ' node(s)'

    case default

       if (pub_on_root) write(stdout,'(a,i6)') &
            'Error in timer_clock: unknown option ',option
       call comms_abort

    end select

  end subroutine timer_clock


  real(kind=DP) function wrappers_etime()

    implicit none

    ! Local variables
#ifdef MPI
    real(kind=DP), external :: MPI_WTIME
#else
    real(kind=kind(0.0)) :: tmp_cpu_time   ! Intrinsic call is single-precision
#endif

#ifdef MPI
    wrappers_etime = MPI_WTIME()
#else
    call cpu_time(tmp_cpu_time)
    wrappers_etime = real(tmp_cpu_time,DP)
#endif

  end function wrappers_etime

end module timer

#else
!######################################################################!
!#                                                                    #!
!#                                                                    #!
!#              This module belongs to the package "ONES"             #!
!#                                                                    #!
!#                                                                    #!
!#                          Oswaldo Dieguez                           #!
!#                 The Theory of Condensed Matter Group               #!
!#            Cavendish Laboratory (University of Cambridge)          #!
!#                           Madingley Road                           #!
!#                         Cambridge  CB3 0HE                         #!
!#                           United Kingdom                           #!
!#                                                                    #!
!#                        odl22@phy.cam.ac.uk                         #!
!#                                                                    #!
!#                                                                    #!
!#                          27 November 2001                          #!
!#                                                                    #!
!#                                                                    #!
!######################################################################!

! USAGE:
! This module can be used to measure the computation time for different
! parts of a program.
! Follow these steps:
! 1) Include the line "USE timer" in the subroutines where timings are
!    to be performed.
! 2) Include the line "CALL clock('TOTAL TIME:',0)" before the first
!    executable order in the main program. Include the line
!    "CALL clock('TOTAL TIME:',3)" after the last executable order in
!    the main program.
! 3) Include the lines "CALL clock(string,1)" and "CALL clock(string,2)"
!    delimiting the beginning and end of a part of the program to be
!    timed. "string" is a string of characters chosen by the user that
!    will identify on output the part of the program timed.

!=====================================================================!
! Adapted for the parallel version of ONETEP by Chris-Kriton Skylaris !
! in November 2003.                                                   !
! Revised for MPI by Peter Haynes, July 2004                          !
! Now also counts amount of work (floating points operations or bytes !
! transferred).                                                       !
!=====================================================================!

module timer

  use constants, only: DP

  implicit none

  private

  type timing_info
     integer :: num_calls           ! no. of calls
     real(kind=DP) :: cum_work      ! cumulative amount of work
     real(kind=DP) :: cum_time      ! cumulative time
     real(kind=DP) :: self_time     ! cum. time - cum. time of timed children
     real(kind=DP) :: last_time     ! last time
     character(len=50) :: name      ! name label
  end type timing_info

  integer, parameter :: max_names = 100  ! no. of different timings possible
  integer, parameter :: max_call_stack_depth = 1000 ! max. depth of call stack
  integer :: num_names = 0
  type(timing_info) :: routines(max_names)
  integer :: call_stack_depth = 0
  integer :: call_stack(max_call_stack_depth)

  public :: timer_clock

contains

  subroutine timer_clock(name,option,work)

    use comms, only: comms_barrier, comms_abort, pub_on_root, comms_reduce, &
         pub_my_node_id, pub_total_num_nodes
    use constants, only: DP, stdout
    use rundat, only: timings_level

    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in) :: option
    real(kind=DP), optional, intent(in) :: work

    ! cks: internal declarations
    integer :: iname
    real(kind=DP) :: tot_time, av_tot_time
    real(kind=DP) :: percent
    real(kind=DP) :: sum_time
    real(kind=DP) :: time_passed
    real(kind=DP) :: sum_work
    integer :: node_count
    integer :: sum_calls
    integer :: calls_per_node      ! no. of calls per node
    integer :: min_global_num_names,max_global_num_names
    integer :: parent
    character(len=8) :: rate_string
    character(len=50) :: local_name

    if (timings_level == 0 .and. (option == 1 .or. option == 2)) return

#ifdef BARRIER_IN_TIMER
    call comms_barrier()
#endif

    select case (option)

    case (0)
       ! cks: initialise variables corresponding to the start of time
       num_names = 1
       routines(:)%cum_work = 0.0_DP
       routines(:)%cum_time = 0.0_DP
       routines(:)%self_time = 0.0_DP
       routines(:)%num_calls = 0
       routines(1)%num_calls = 1
       routines(1)%name = name
       routines(1)%last_time = wrappers_etime()
       ! jd: Manually put 'total_time' as the first entry in the call stack
       call_stack(1) = 1
       call_stack_depth = 1
       return

    case (1)
       ! cks: Loop over all names used so far, looking for this one
       do iname=1,num_names
          ! cks: break out if current time "tag" has been used before
          if (name == routines(iname)%name) exit
       end do

       ! jd: If the loop ran to completion, iname = num_names + 1 and we
       !     don't have this "tag" yet. In both cases, update this tag,
       !     but first make sure we don't go past the end of the array

       ! cks: make sure you don't exceed allocated memory
       if (iname > max_names) then
          if (pub_on_root) write(stdout,'(a,i3,a)') 'Error in timer_clock: &
               &maximum number of timing tags (',max_names, &
               ' exceeded'
          call comms_abort
       end if

       routines(iname)%name = name
       routines(iname)%last_time = wrappers_etime()
       routines(iname)%num_calls = routines(iname)%num_calls + 1

       ! jd: Put the index of this "tag" on top of the call stack
       call_stack_depth = call_stack_depth + 1
       if (call_stack_depth > max_call_stack_depth) then
          if (pub_on_root) write(stdout,*) 'Error in timer_clock: maximum &
               &call stack depth exceeded. Check for unbounded recursion.'
          call comms_abort
       end if
       call_stack(call_stack_depth) = iname

       ! jd: If this was a new tag, increase count
       if(iname == num_names + 1) num_names = num_names + 1

    case (2)

       if (routines(call_stack(call_stack_depth))%name /= name) then
          if (pub_on_root) write(stdout,*) 'Error in timer_clock: &
               &call stack inconsistent. Check for mismatched calls to &
               &timer_clock. Timer indicates exit from ',name, &
               'while the most recently entered routine was ', &
               routines(call_stack(call_stack_depth))%name
          call comms_abort
       end if

       ! cks: Loop over all names used so far, looking for this one
       do iname=1,num_names
          if (name == routines(iname)%name) then

             ! cks: add time to the time already stored for the current "tag"
             time_passed = wrappers_etime() - routines(iname)%last_time
             routines(iname)%cum_time = routines(iname)%cum_time + time_passed
             routines(iname)%self_time = routines(iname)%self_time + &
                  time_passed

             ! jd: Add the work for this "tag"
             if (present(work)) routines(iname)%cum_work = &
                  routines(iname)%cum_work + work

             ! jd: Pop this finished "tag" from the top of the call stack
             call_stack(call_stack_depth) = -1
             call_stack_depth = call_stack_depth - 1
             if (call_stack_depth < 0) then
                if (pub_on_root) write(stdout,*) 'Error in timer_clock: &
                     &call stack inconsistent. Check for mismatched calls to &
                     &timer_clock.'
                call comms_abort
             end if

             ! jd: Subtract the timing of this routine from the self time
             !     of the routine which is now the top of the call stack
             !     (and which is the parent)
             parent = call_stack(call_stack_depth)
             routines(parent)%self_time = routines(parent)%self_time - &
                  time_passed

             return

          end if
       end do

       if (pub_on_root) write(stdout,'(a)') 'WARNING in timer_clock: &
            &matching tag "',trim(name),'" not found for option 2'

    case (3)
       ! cks: total time the program was running
       tot_time = wrappers_etime() - routines(1)%last_time
       av_tot_time = tot_time
       call comms_reduce('SUM',av_tot_time)
       av_tot_time = av_tot_time / pub_total_num_nodes
       if (pub_on_root) then
          write(stdout,'(/a)') '------------------------------ TIMING INFORMATI&
               &ON ------------------------------'
          write(stdout,'(a,f12.3,a,i6,a)') 'TOTAL TIME:',av_tot_time,'s on', &
               pub_total_num_nodes,' node(s)'
       end if
       routines(1)%cum_time = tot_time
       routines(1)%self_time = routines(1)%self_time + tot_time

       if (timings_level == 0) return

       ! ********************************
       ! *** Timing breakdown by node ***
       ! ********************************
       if (iand(timings_level,4) /= 0) then

          ! cks: print name of each "tag", number of calls, times and
          !      percent of total time
          node_loop: do node_count=0,pub_total_num_nodes-1

             call comms_barrier

             if (pub_my_node_id == node_count) then

                ! cks: print timings regarding pub_my_node_id
                write(stdout,'(a)') ' TAG                                     &
                     &n calls  cpu time  %total   node'
                do iname=1,num_names
                   percent = 100.0_DP * routines(iname)%cum_time / tot_time

                   write(stdout,'(a40,a1,i7,f10.3,a1,f8.3,a1,i5,a2)') &
                        routines(iname)%name,':',routines(iname)%num_calls, &
                        routines(iname)%cum_time,'s',percent,'%', &
                        pub_my_node_id,' |'

                end do

             end if

          end do node_loop

       end if

       ! cks: print maximum timings across all nodes

       call comms_barrier

       ! cks: see if every node has timings for the same tasks
       min_global_num_names = num_names
       call comms_reduce('MIN', min_global_num_names)
       max_global_num_names = num_names
       call comms_reduce('MAX', max_global_num_names)
       if (min_global_num_names /= max_global_num_names .and. pub_on_root) &
            write(stdout,'(a)') 'WARNING in timer_clock: inconsistent timing &
                 &data'

       ! **********************************
       ! *** Usual (cumulative) timings ***
       ! **********************************
       if (iand(timings_level,1) /= 0) then

          if (pub_on_root) then

             write(stdout,'(a)') '================== AVERAGE TIMINGS FROM ALL &
                  &NODES (CUMULATIVE) ================='
             write(stdout,'(a)') '|| TAG                                   &
                  &n calls    cpu time   %total  Gflops ||'

          end if

          do iname=1,min_global_num_names

             sum_calls = routines(iname)%num_calls
             call comms_reduce('SUM',sum_calls)

             sum_time = routines(iname)%cum_time
             call comms_reduce('SUM',sum_time)

             sum_work = routines(iname)%cum_work
             call comms_reduce('SUM',sum_work)

             percent = 100.0_DP * sum_time / (av_tot_time * pub_total_num_nodes)

             if (sum_time > epsilon(1.0_DP) .and. sum_work > 0.0_DP) then
                write(rate_string,'(f8.3)') 1.0e-9 * sum_work * &
                     pub_total_num_nodes / sum_time
             else
                rate_string = '  ------'
             end if

             if (pub_on_root) then
                calls_per_node = sum_calls/pub_total_num_nodes
                if(calls_per_node < 100000000) then
                   write(stdout,&
                        '(a3,a36,a1,i8,a1,f10.2,a1,f8.3,a1,a8,a3)') '|| ', &
                        routines(iname)%name,':',calls_per_node,&
                        ' ',sum_time/pub_total_num_nodes,'s',percent,'%',&
                        rate_string, ' ||'
                else if(calls_per_node < 200000000) then
                   write(stdout,&
                        '(a3,a36,a1,i7,a2,f10.2,a1,f8.3,a1,a8,a3)') '|| ', &
                        routines(iname)%name,':',calls_per_node/1000,&
                        'K ',sum_time/pub_total_num_nodes,'s',percent,'%',&
                        rate_string, ' ||'
                else
                   write(stdout,&
                        '(a3,a36,a1,i7,a2,f10.2,a1,f8.3,a1,a8,a3)') '|| ', &
                        routines(iname)%name,':',calls_per_node/1000000,&
                        'M ',sum_time/pub_total_num_nodes,'s',percent,'%',&
                        rate_string, ' ||'
                end if
             end if

          end do

          if (pub_on_root) write(stdout,'(a)') '===============================&
               &================================================='

       end if

       ! ******************
       ! *** Self-times ***
       ! ******************
       if (iand(timings_level,2) /= 0) then

          call comms_barrier

          if (pub_on_root) then

             write(stdout,'(a)') '==================== AVERAGE TIMINGS FROM ALL&
                  & NODES (SELF) ====================='
             write(stdout,'(a)') '|| TAG                                   &
                  &n calls    cpu time   %total  Gflops ||'

          end if

          do iname=1,min_global_num_names

             sum_calls = routines(iname)%num_calls
             call comms_reduce('SUM',sum_calls)

             sum_time = routines(iname)%self_time
             call comms_reduce('SUM',sum_time)

             sum_work = routines(iname)%cum_work
             call comms_reduce('SUM',sum_work)

             percent = 100.0_DP * sum_time / (av_tot_time * pub_total_num_nodes)

             if (sum_time > epsilon(1.0_DP) .and. sum_work > 0.0_DP) then
                write(rate_string,'(f8.3)') 1.0e-9 * sum_work * &
                     pub_total_num_nodes / sum_time
             else
                rate_string = '  ------'
             end if

             local_name = routines(iname)%name
             if(iname == 1) local_name = 'main program (onetep.F90)'

             if (pub_on_root) then
                calls_per_node = sum_calls/pub_total_num_nodes
                if(calls_per_node < 100000000) then
                   write(stdout,&
                        '(a3,a36,a1,i8,a1,f10.2,a1,f8.3,a1,a8,a3)') '|| ', &
                        routines(iname)%name,':',calls_per_node,&
                        ' ',sum_time/pub_total_num_nodes,'s',percent,'%',&
                        rate_string, ' ||'
                else if(calls_per_node < 200000000) then
                   write(stdout,&
                        '(a3,a36,a1,i7,a2,f10.2,a1,f8.3,a1,a8,a3)') '|| ', &
                        routines(iname)%name,':',calls_per_node/1000,&
                        'K ',sum_time/pub_total_num_nodes,'s',percent,'%',&
                        rate_string, ' ||'
                else
                   write(stdout,&
                        '(a3,a36,a1,i7,a2,f10.2,a1,f8.3,a1,a8,a3)') '|| ', &
                        routines(iname)%name,':',calls_per_node/1000000,&
                        'M ',sum_time/pub_total_num_nodes,'s',percent,'%',&
                        rate_string, ' ||'
                end if
             end if

          end do

          if (pub_on_root) write(stdout,'(a)') '===============================&
               &================================================='

       end if


    case default

       if (pub_on_root) write(stdout,'(a,i6)') &
            'Error in timer_clock: unknown option ',option
       call comms_abort

    end select

  end subroutine timer_clock


  real(kind=DP) function wrappers_etime()

    implicit none

    ! internal variables
#ifdef MPI
    real(kind=DP), external :: MPI_WTIME

    wrappers_etime = MPI_WTIME()
#else
    integer, parameter :: SP  = KIND(1.0)
    real(kind=SP) :: delta,tarray(2),etime

    delta = etime(tarray)
    wrappers_etime = real(delta,DP)
#endif

  end function wrappers_etime

end module timer
#endif




