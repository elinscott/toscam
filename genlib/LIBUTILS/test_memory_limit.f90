module test_memory_limit



contains

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


subroutine test_mem

!    MEMORY_TEST declares larger and larger arrays of various types,
!    to see when an error occurs because of memory limitations or other
!    limits.

  implicit none

  integer n
  integer n_log
  integer n_log_max

  n_log_max = 27

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MEMORY_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Try to see how big a single vector can be.'

  n = 1
  n_log = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4VEC Memory Test'
  write ( *, '(a)' ) ' '  
  write ( *, '(a)' ) &
    '        Log2(N)            N     Average         CPU         Real'
  write ( *, '(a)' ) ' '

  do while ( n_log <= n_log_max )

    call i4vec_memory_test ( n_log, n )

    n = n * 2
    n_log = n_log + 1

  end do

  n = 1
  n_log = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R4VEC Memory Test'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '        Log2(N)            N     Average         CPU         Real'
  write ( *, '(a)' ) ' '

  do while ( n_log <= n_log_max )

    call r4vec_memory_test ( n_log, n )

    n = n * 2
    n_log = n_log + 1

  end do

  n = 1
  n_log = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC Memory Test'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '        Log2(N)            N     Average         CPU         Real'
  write ( *, '(a)' ) ' '

  do while ( n_log <= n_log_max )

    call r8vec_memory_test ( n_log, n )

    n = n * 2
    n_log = n_log + 1

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MEMORY_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************



subroutine i4vec_memory_test ( n_log, n )
  implicit none

  real ( kind = 4 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer i
  integer ( kind = 4 ), allocatable, dimension ( : ) :: i4vec
  integer n
  integer n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2
  real ( kind = 4 ) x

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( i4vec(1:n) )

  do i = 1, n
    call random_number ( harvest = x )
    i4vec(i) = int ( 3.0E+00 * x )
  end do

  average = real ( sum ( i4vec(1:n) ), kind = 4 ) / real ( n, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i12,2x,i12,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    n_log, n, average, cpu_diff, real_diff

  deallocate ( i4vec )

  return
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************



subroutine r4vec_memory_test ( n_log, n )
  implicit none

  real ( kind = 4 ) average
  real cpu_diff
  real cpu_time1
  real cpu_time2
  integer i
  real ( kind = 4 ), allocatable, dimension ( : ) :: r4vec
  integer n
  integer n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2
  real ( kind = 4 ) x

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( r4vec(1:n) )

  do i = 1, n
    call random_number ( harvest = x )
    r4vec(i) = 2.0E+00 * x
  end do

  average = sum ( r4vec(1:n) ) / real ( n, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i12,2x,i12,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    n_log, n, average, cpu_diff, real_diff

  deallocate ( r4vec )

  return
end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************



subroutine r8vec_memory_test ( n_log, n )
  implicit none

  real ( kind = 8 ) average
  real cpu_diff
  real cpu_time1
  real cpu_time2
  integer i
  real ( kind = 8 ), allocatable, dimension ( : ) :: r8vec
  integer n
  integer n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2
  real ( kind = 8 ) x

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( r8vec(1:n) )

  do i = 1, n
    call random_number ( harvest = x )
    r8vec(i) = 2.0D+00 * x
  end do

  average = sum ( r8vec(1:n) ) / real ( n, kind = 8 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i12,2x,i12,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    n_log, n, average, cpu_diff, real_diff

  deallocate ( r8vec )

  return
end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************



subroutine real_time ( seconds )
  implicit none

  integer clock_count
  integer clock_max
  integer clock_rate
  real ( kind = 8 ) seconds

  call system_clock ( clock_count, clock_rate, clock_max )

  seconds = real ( clock_count, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  return
end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


subroutine timestamp ( )
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


end module

