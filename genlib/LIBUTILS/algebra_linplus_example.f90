program main

!*****************************************************************************80
!
!! MAIN is the main program for LINPLUS_PRB.
!
!  Discussion:
!
!    LINPLUS_PRB calls the LINPLUS test routines.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPLUS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test the routines in the LINPLUS library.'

  call test01
  call test014
  call test015
  call test016
  call test017
  call test019
  call test0192
  call test0193
  call test0194
  call test0195
  call test0196
  call test02
  call test03
  call test035
  call test037
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09
  call test10
  call test105

  call test11
  call test115
  call test12
  call test13
  call test14
  call test15
  call test151
  call test152
  call test153
  call test154
  call test155
  call test156
  call test1565
  call test1566
  call test157
  call test16
  call test17
  call test18
  call test19
  call test193
  call test195
  call test197
  call test20

  call test21
  call test22
  call test23
  call test235
  call test24
  call test25
  call test26
  call test265
  call test2655
  call test27
  call test275
  call test28
  call test285
  call test29
  call test295
 
  call test30
  call test31
  call test315
  call test317
  call test32
  call test33
  call test34
  call test35
  call test36
  call test37
  call test38
  call test385
  call test39

  call test40
  call test41
  call test42
  call test422
  call test423
  call test425
  call test426
  call test428
  call test43
  call test44
  call test443
  call test445
  call test45
  call test46
  call test47
  call test48
  call test485
  call test49

  call test50
  call test505
  call test51
  call test515
  call test517
  call test52
  call test525
  call test527
  call test53
  call test531
  call test532
  call test533
  call test534
  call test535
  call test54
  call test55
  call test555
  call test56
  call test57
  call test5705
  call test571
  call test572
  call test5722
  call test5724
  call test5725
  call test573
  call test574
  call test5745
  call test575
  call test577
  call test58
  call test581
  call test583
  call test585
  call test587
  call test589
  call test59

  call test60
  call test605
  call test61
  call test62
  call test63

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LINPLUS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 tests C83_CR_FA, C83_CR_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) diag(2*n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_10 = 10
  complex ( kind = 8 ) rhs(0:2*n)
  complex ( kind = 8 ) subd(0:2*n)
  complex ( kind = 8 ) supd(0:2*n)
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  C83_CR_FA factors a complex tridiagonal matrix;'
  write ( *, '(a)' ) '  C83_CR_SL solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  do i = 1, n
    diag(i) = cmplx ( 2.0D+00, 2.0D+00 * real ( i, kind = 8 ) )
  end do

  do i = 1, n-1
    subd(i) = cmplx ( -1.0D+00, - real ( i, kind = 8 ) )
  end do

  do i = 1, n-1
    supd(i) = cmplx ( -1.0D+00, real ( -i-1, kind = 8 ) )
  end do
!
!  Set the desired solution.
!
  do i = 1, n
    x(i) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  i = 1
  rhs(1) = cmplx ( 2.0D+00, 2.0D+00 * real ( i, kind = 8 ) ) * x(i) &
         + cmplx ( -1.0D+00, real ( -i-1, kind = 8 ) ) * x(i+1)

  do i = 2, n-1
    rhs(i) = cmplx ( -1.0D+00, real ( -i+1, kind = 8 ) ) * x(i-1) &
           + cmplx (  2.0D+00, real ( 2*i, kind = 8  ) ) * x(i) &
           + cmplx ( -1.0D+00, real ( -i-1, kind = 8 ) ) * x(i+1)
  end do

  rhs(n) = cmplx ( -1.0D+00, real ( -i+1, kind = 8 ) ) * x(i-1) &
         + cmplx (  2.0D+00, real ( 2*i, kind = 8 ) ) * x(i)
!
!  Factor the matrix.
!
  call c83_cr_fa ( n, subd, diag, supd )
!
!  Solve the linear system.
!
  call c83_cr_sl ( n, subd, diag, supd, rhs )

  call c8vec_print_some ( n, rhs, i4_10, '  Solution:' )

  return
end
subroutine test014

!*****************************************************************************80
!
!! TEST014 tests C83_NP_FA, C83_NP_SL, C83_NP_ML.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) a1(2:n)
  complex ( kind = 8 ) a2(1:n)
  complex ( kind = 8 ) a3(1:n-1)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  For a complex tridiagonal matrix that can be'
  write ( *, '(a)' ) '    factored with no pivoting,'
  write ( *, '(a)' ) '  C83_NP_FA factors;'
  write ( *, '(a)' ) '  C83_NP_SL solves a factored system.'
  write ( *, '(a)' ) '  C83_NP_ML multiplies A*X when A has been factored.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c83_random ( n, seed, a1, a2, a3 )

  call c83_print ( n, a1, a2, a3, '  The tridiagonal matrix' )
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  call c83_mxv ( n, a1, a2, a3, x, b )

  call c8vec_print ( n, b, '  The right hand side' )
!
!  Factor the matrix.
!
  call c83_np_fa ( n, a1, a2, a3, info )
!
!  Solve the linear system.
!
  job = 0
  call c83_np_sl ( n, a1, a2, a3, b, job )

  call c8vec_print ( n, b, '  The solution' )
!
!  Now set a SECOND desired solution.
!
  do i = 1, n
    x(i) = cmplx ( real ( 10 * i, kind = 8 ), real ( i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side, using the FACTORED matrix.
!
  call c83_np_ml ( n, a1, a2, a3, x, b, job )

  call c8vec_print ( n, b, '  The second right hand side' )
!
!  Solve the linear system.
!
  call c83_np_sl ( n, a1, a2, a3, b, job )

  call c8vec_print ( n, b, '  The second solution' )

  return
end
subroutine test015

!*****************************************************************************80
!
!! TEST015 tests C83_NP_FA, C83_NP_ML and C83_NP_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) a1(2:n)
  complex ( kind = 8 ) a2(1:n)
  complex ( kind = 8 ) a3(1:n-1)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  For a complex tridiagonal matrix that can be'
  write ( *, '(a)' ) '    factored with no pivoting,'
  write ( *, '(a)' ) '  C83_NP_FA factors;'
  write ( *, '(a)' ) '  C83_NP_SL solves a factored system.'
  write ( *, '(a)' ) '  C83_NP_ML multiplies A*X when A has been factored.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will look at the TRANSPOSED linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c83_random ( n, seed, a1, a2, a3 )

  call c83_print ( n, a1, a2, a3, '  The tridiagonal matrix' )
!
!  Set the desired solution
!
  do i = 1, n
    x(i) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  call c83_vxm ( n, a1, a2, a3, x, b )

  call c8vec_print ( n, b, '  The right hand side B1' )
!
!  Factor the matrix.
!
  call c83_np_fa ( n, a1, a2, a3, info )
!
!  Solve the linear system.
!
  job = 1
  call c83_np_sl ( n, a1, a2, a3, b, job )
 
  call c8vec_print ( n, b, '  The solution to At * X1 = B1' )
!
!  Set the second solution.
!
  do i = 1, n
    x(i) = cmplx ( real ( 10 * i, kind = 8 ), real ( i, kind = 8 ) )
  end do
!
!  Compute the corresponding right hand side.
!
  call c83_np_ml ( n, a1, a2, a3, x, b, job )

  call c8vec_print ( n, b, '  The second right hand side B2' )
!
!  Solve the linear system.
!
  call c83_np_sl ( n, a1, a2, a3, b, job )
 
  call c8vec_print ( n, b, '  Solution to At * X2 = B2' )
 
  return
end
subroutine test016

!*****************************************************************************80
!
!! TEST016 tests C8CI_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  C8CI_SL solves a complex circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c8ci_random ( n, seed, a )

  call c8ci_print ( n, a, '  The circulant matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    do i = 1, n
      x(i) = cmplx ( real ( i, kind = 8 ), real ( 10 * i, kind = 8 ) )
    end do
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call c8ci_mxv ( n, a, x, b )
    else
      call c8ci_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call c8ci_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call c8vec_print ( n, x, '  Solution:' )
    else
      call c8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine test017

!*****************************************************************************80
!
!! TEST017 tests C8TO_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  complex ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  C8TO_SL solves a complex Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call c8to_random ( n, seed, a )

  call c8to_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call c8vec_indicator ( n, x )

    call c8vec_print_some ( n, x, i4_10, '  Desired solution:' )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call c8to_mxv ( n, a, x, b )
    else
      call c8to_vxm ( n, a, x, b )
    end if

    call c8vec_print_some ( n, b, i4_10, '  Right hand side:' )
!
!  Solve the linear system.
!
    call c8to_sl ( n, a, b, x, job )

    call c8vec_print_some ( n, x, i4_10, '  Computed solution:' )

  end do

  return
end
subroutine test019

!*****************************************************************************80
!
!! TEST019 tests C8VEC_UNITY.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 5

  integer ( kind = 4 ) n
  complex ( kind = 8 ) x(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  C8VEC_UNITY returns the N complex roots of unity.'
  write ( *, '(a)' ) ' '

  do n = 1, n_max

    call c8vec_unity ( n, x )

    call c8vec_print ( n, x, '  Roots of unity:' )

  end do

  return
end
subroutine test0192

!*****************************************************************************80
!
!! TEST0192 tests R8CC_WRITE.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real    ( kind = 8 ) a(nz_num)
  character ( len = 127 ) :: a_file = 'R8CC_a.txt'
  integer ( kind = 4 ), dimension (n+1) :: col = (/ 1, 4, 6, 8, 10, 13 /)
  character ( len = 127 ) :: col_file = 'R8CC_col.txt'
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  character ( len = 127 ) :: row_file = 'R8CC_row.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0192'
  write ( *, '(a)' ) '  For a matrix in the R8CC format,'
  write ( *, '(a)' ) '  (double precision compressed column sparse)'
  write ( *, '(a)' ) '  R8CC_WRITE writes the matrix to 3 files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + i4_1, col, '  The COL vector:' )

  call i4vec_print ( nz_num, row, '  The ROW vector:' )

  call R8CC_indicator ( m, n, nz_num, col, row, a )

  call R8CC_print ( m, n, nz_num, col, row, a, '  The R8CC matrix:' )

  call R8CC_write ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

  return
end
subroutine test0193

!*****************************************************************************80
!
!! TEST0193 tests R8CC_READ_SIZE and R8CC_READ.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ), allocatable, dimension ( : ) :: a
  character ( len = 127 ) :: a_file = 'R8CC_a.txt'
  integer ( kind = 4 ) base
  integer ( kind = 4 ), allocatable, dimension ( : ) :: col
  character ( len = 127 ) :: col_file = 'R8CC_col.txt'
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: row
  character ( len = 127 ) :: row_file = 'R8CC_row.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0193'
  write ( *, '(a)' ) '  For a matrix in the R8CC format,'
  write ( *, '(a)' ) '  (double precision compressed column sparse)'
  write ( *, '(a)' ) '  R8CC_READ_SIZE reads the sizes of the data;'
  write ( *, '(a)' ) '  R8CC_READ reads the data.'

  call R8CC_read_size ( col_file, row_file, m, n, nz_num, base )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
  write ( *, '(a,i8)' ) '  Index base (0/1)  = ', base

  allocate ( a(1:nz_num) )
  allocate ( col(1:n+1) )
  allocate ( row(1:nz_num) )

  call R8CC_read ( col_file, row_file, a_file, m, n, nz_num, col, row, a )

  call i4vec_print ( n + i4_1, col, '  The COL vector:' )

  call i4vec_print ( nz_num, row, '  The ROW vector:' )

  call R8CC_print ( m, n, nz_num, col, row, a, '  The R8CC matrix:' )

  deallocate ( a )
  deallocate ( col )
  deallocate ( row )
!
!  Delete the files.  
!  If you want to look at the files, comment out these lines
!  and rerun the problem!
!
  call file_delete ( a_file )
  call file_delete ( col_file )
  call file_delete ( row_file )

  return
end
subroutine test0194

!*****************************************************************************80
!
!! TEST0194 tests R8VEC_TO_R8CB, R8CB_TO_R8VEC.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real    ( kind = 8 ) a(ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) x((ml+mu+1)*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0194'
  write ( *, '(a)' ) '  For a compressed banded matrix,'
  write ( *, '(a)' ) '  R8VEC_TO_R8CB converts a real vector to an R8CB matrix.'
  write ( *, '(a)' ) '  R8CB_TO_R8VEC converts an R8CB matrix to a real vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M      = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N   = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB indicator matrix:' )

  call r8cb_to_r8vec ( m, n, ml, mu, a, x )

  k = 0
  do j = 1, n
    do i = 1, ml+mu+1
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_r8cb ( m, n, ml, mu, x, a )

  call r8cb_print ( m, n, ml, mu, a, '  The recovered R8CB indicator matrix:' )

  return
end
subroutine test0195
!
!*****************************************************************************80
!
!! TEST0195 tests R8VEC_TO_R8GB, R8GB_TO_R8VEC.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) x((2*ml+mu+1)*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0195'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8VEC_TO_R8GB converts a real vector to an R8GB matrix.'
  write ( *, '(a)' ) '  R8GB_TO_R8VEC converts an R8GB matrix to a real vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB indicator matrix:' )

  call r8gb_to_r8vec ( m, n, ml, mu, a, x )

  k = 0
  do j = 1, n
    do i = 1, 2*ml+mu+1
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_r8gb ( m, n, ml, mu, x, a )

  call r8gb_print ( m, n, ml, mu, a, '  The recovered R8GB indicator matrix:' )

  return
end
subroutine test0196

!*****************************************************************************80
!
!! TEST0196 tests R8VEC_TO_DGE, DGE_TO_R8VEC.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 6

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) x(m*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0196'
  write ( *, '(a)' ) '  For a general matrix,'
  write ( *, '(a)' ) '  R8VEC_TO_SGE converts a real vector to an DGE matrix.'
  write ( *, '(a)' ) '  DGE_TO_R8VEC converts an DGE matrix to a real vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  call dge_indicator ( m, n, a )

  call dge_print ( m, n, a, '  The DGE indicator matrix:' )

  call dge_to_r8vec ( m, n, a, x )

  k = 0
  do j = 1, n
    do i = 1, m
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_dge ( m, n, x, a )

  call dge_print ( m, n, a, '  The recovered DGE indicator matrix:' )

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 tests R83_CR_FA, R83_CR_SL.
!
!  Modified:
!
!    11 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) a_cr(3,0:2*n)
  real    ( kind = 8 ) b(n)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n2
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R83_CR_FA factors a real tridiagonal matrix;'
  write ( *, '(a)' ) '  R83_CR_SL solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  Demonstrate multiple system solution method.'
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  a(1,1) = 0.0D+00
  a(1,2:n) = -1.0D+00

  a(2,1:n) = 2.0D+00

  a(3,1:n-1) = -1.0D+00
  a(3,n) = 0.0D+00
!
!  Print the matrix.
!
  if ( debug ) then
    call r83_print ( n, a, '  Input matrix:' )
  end if
!
!  Factor the matrix once.
!
  call r83_cr_fa ( n, a, a_cr )
!
!  Print the factor information.
!
  if ( debug ) then
    n2 = 2 * n + 1
    call r83_print ( n2, a_cr, '  Cyclic reduction factor information:' )
  end if

  do j = 1, 2

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) '  Solve linear system number #', j, '.'

    if ( j == 1 ) then
      b(1:n-1) = 0.0D+00
      b(n) = real ( n + 1, kind = 8 )
    else
      b(1) = 1.0D+00
      b(2:n-1) = 0.0D+00
      b(n) = 1.0D+00
    end if
!
!  Solve the linear system.
!
    call r83_cr_sl ( n, a_cr, b, x )

    call r8vec_print_some ( n, x, i4_10, '  Solution:' )

  end do

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 tests R83_CR_FA, R83_CR_SL.
!
!  Modified:
!
!    11 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) a_cr(3,0:2*n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For a real tridiagonal matrix,'
  write ( *, '(a)' ) '  using CYCLIC REDUCTION,'
  write ( *, '(a)' ) '  R83_CR_FA factors;'
  write ( *, '(a)' ) '  R83_CR_SL solves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  The matrix is NOT symmetric.'
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  a(1,1) = 0.0D+00
  do j = 2, n
    a(1,j) = real ( j, kind = 8 )
  end do

  do j = 1, n
    a(2,j) = 4.0D+00 * real ( j, kind = 8 )
  end do

  do j = 1, n - 1
    a(3,j) = real ( j, kind = 8 )
  end do
  a(3,n) = 0.0D+00

  call r83_print ( n, a, '  The matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83_mxv ( n, a, x, b )
  x(1:n) = 0.0D+00
!
!  Factor the matrix.
!
  call r83_cr_fa ( n, a, a_cr )
!
!  Solve the linear system.
!
  call r83_cr_sl ( n, a_cr, b, x )

  call r8vec_print_some ( n, x, i4_10, '  Solution:' )

  return
end
subroutine test035

!*****************************************************************************80
!
!! TEST035 tests R83_GS_SL.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) it
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: it_max = 1000
  real    ( kind = 8 ) :: tol = 0.000001D+00
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  For a real tridiagonal system,'
  write ( *, '(a)' ) '  R83_GS_SL solves a linear system using'
  write ( *, '(a)' ) '    Gauss-Seidel iteration'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  a(1,1)     =  0.0D+00
  a(1,2:n)   = -1.0D+00
  a(2,1:n)   =  2.0D+00
  a(3,1:n-1) = -1.0D+00
  a(3,n)     =  0.0D+00

  do job = 0, 1

    if ( job == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solving A * x = b.'
      write ( *, '(a)' ) ' '
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solving A'' * x = b.'
      write ( *, '(a)' ) ' '
    end if
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r83_mxv ( n, a, x, b )
    else
      call r83_vxm ( n, a, x, b )
    end if
!
!  Set the starting solution.
!
    x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
    do i = 1, 3

      call r83_gs_sl ( n, a, b, x, tol, it_max, job, it, diff )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
      write ( *, '(a,g14.6)' ) '  Maximum solution change on last step = ', diff

      call r8vec_print_some ( n, x, i4_10, '  Current solution estimate:' )

    end do

  end do

  return
end
subroutine test037

!*****************************************************************************80
!
!! TEST037 tests R83_INDICATOR.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
 implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037'
  write ( *, '(a)' ) '  R83_INDICATOR sets up an R83 indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83_indicator ( n, a )

  call r83_print ( n, a, '  The R83 indicator matrix:' )

  return
end
subroutine test04

!*****************************************************************************80
!
!! TEST04 tests R83_JAC_SL.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) it
  integer ( kind = 4 ) job
  integer ( kind = 4 ) it_max
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) x(n)

  it_max = 1000
  tol = 0.0000001D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For a real tridiagonal system,'
  write ( *, '(a)' ) '  R83_JAC_SL solves a linear system using'
  write ( *, '(a)' ) '    Jacobi iteration'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  a(1,1)     =  0.0D+00
  a(1,2:n)   = -1.0D+00
  a(2,1:n)   =  2.0D+00
  a(3,1:n-1) = -1.0D+00
  a(3,n)     =  0.0D+00

  do job = 0, 1

    if ( job == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solving A * x = b.'
      write ( *, '(a)' ) ' '
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solving A'' * x = b.'
      write ( *, '(a)' ) ' '
    end if
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r83_mxv ( n, a, x, b )
    else
      call r83_vxm ( n, a, x, b )
    end if

    call r8vec_print_some ( n, b, i4_10, '  The right hand side:' )
!
!  Set the starting solution.
!
    x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
    do i = 1, 3

      call r83_jac_sl ( n, a, b, x, tol, it_max, job, it, diff )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
      write ( *, '(a,g14.6)' ) '  Maximum solution change on last step = ', diff

      call r8vec_print_some ( n, x, i4_10, '  Current solution estimate:' )

    end do

  end do

  return
end
subroutine test05

!*****************************************************************************80
!
!! TEST05 tests R83_NP_DET, R83_NP_FA.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For a tridiagonal matrix that can be factored'
  write ( *, '(a)' ) '    with no pivoting,'
  write ( *, '(a)' ) '  R83_NP_FA factors,'
  write ( *, '(a)' ) '  R83_NP_DET computes the determinant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83_random ( n, seed, a )
!
!  Copy the matrix into general storage.
!
  call r83_to_dge ( n, a, b )
!
!  Factor the matrix.
!
  call r83_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Warning!'
    write ( *, '(a)' ) '  R83_NP_FA returns INFO = ', info
  end if

  call r83_print ( n, a, '  The factored R83 matrix:' )
!
!  Compute the determinant.
!
  call r83_np_det ( n, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R83_NP_DET computes determinant = ', det
!
!  Factor the matrix in DGE storage.
!
  call dge_np_fa ( n, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Warning!'
    write ( *, '(a)' ) '  DGE_NP_FA returns INFO = ', info
  end if
!
!  Compute the determinant of the matrix in DGE storage.
!
  call dge_np_det ( n, b, det )

  write ( *, '(a,g14.6)' ) '  DGE_NP_DET computes determinant = ', det

  return
end
subroutine test06

!*****************************************************************************80
!
!! TEST06 tests R83_NP_FA, R83_NP_SL.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For a tridiagonal matrix that can be factored'
  write ( *, '(a)' ) '    with no pivoting,'
  write ( *, '(a)' ) '  R83_NP_FA factors;'
  write ( *, '(a)' ) '  R83_NP_SL solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83_random ( n, seed, a )

  call r83_print ( n, a, '  The tridiagonal matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r83_mxv ( n, a, x, b )
  x(1:n) = 0.0D+00
!
!  Factor the matrix.
!
  call r83_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06 - Fatal error!'
    write ( *, '(a)' ) '  The test matrix is singular.'
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r83_np_sl ( n, a, b, job )
 
  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side, using the factored matrix.
!
  job = 1
  call r83_np_ml ( n, a, x, b, job )
!
!  Solve the linear system.
!
  job = 1
  call r83_np_sl ( n, a, b, job )
 
  call r8vec_print ( n, b, '  Solution to tranposed system:' )
 
  return
end
subroutine test07

!*****************************************************************************80
!
!! TEST07 tests R83_NP_FS.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  R83_NP_FS factors and solves a tridiagonal'
  write ( *, '(a)' ) '    linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix elements.
!
  call r83_random ( n, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute b = A * x.
!
  call r83_mxv ( n, a, x, b )
!
!  Wipe out the solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the system.
!
  call r83_np_fs ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution:' )

  return
end
subroutine test08

!*****************************************************************************80
!
!! TEST08 tests R83_NP_ML.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  R83_NP_ML computes A*x or A''*x'
  write ( *, '(a)' ) '    where A has been factored by R83_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call r83_random ( n, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r83_mxv ( n, a, x, b )
    else
      call r83_vxm ( n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r83_np_fa ( n, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08 - Fatal error!'
      write ( *, '(a)' ) '  R83_NP_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r83_np_ml ( n, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, i4_10, '  A*x and PLU*x:' )
    else
      call r8vec2_print_some ( n, b, b2, i4_10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine test09

!*****************************************************************************80
!
!! TEST09 tests R83P_DET.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) work2(n-1)
  real    ( kind = 8 ) work3(n-1)
  real    ( kind = 8 ) work4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  R83P_DET, determinant of a tridiagonal'
  write ( *, '(a)' ) '    periodic matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_random ( n, seed, a )

  call r83p_print ( n, a, '  The periodic tridiagonal matrix:' )
!
!  Copy the matrix into a general array.
!
  call r83p_to_dge ( n, a, b )
!
!  Factor the matrix.
!
  call r83p_fa ( n, a, info, work2, work3, work4 )
!
!  Compute the determinant.
!
  call r83p_det ( n, a, work4, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R83P_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call dge_fa ( n, b, pivot, info )
!
!  Compute the determinant.
!
  call dge_det ( n, b, pivot, det )

  write ( *, '(a,g14.6)' ) '  DGE_DET computes the determinant = ', det

  return
end
subroutine test10

!*****************************************************************************80
!
!! TEST10 tests R83P_FA, R83P_SL.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) work2(n-1)
  real    ( kind = 8 ) work3(n-1)
  real    ( kind = 8 ) work4
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  R83P_FA factors a tridiagonal periodic system.'
  write ( *, '(a)' ) '  R83P_SL solves a tridiagonal periodic system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_random ( n, seed, a )
!
!  Factor the matrix.
!
  call r83p_fa ( n, a, info, work2, work3, work4 )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a,i8)' ) '  R83P_FA returns INFO = ', info
    return
  end if

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    call r83p_ml ( n, a, x, b, job )
!
!  Solve the linear system.
!
    call r83p_sl ( n, a, b, x, job, work2, work3, work4 )

    if ( job == 0 ) then
      call r8vec_print ( n, x, '  Solution:' )
    else
      call r8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine test105

!*****************************************************************************80
!
!! TEST105 tests R83P_INDICATOR.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) '  R83P_INDICATOR sets up an R83P indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r83p_indicator ( n, a )

  call r83p_print ( n, a, '  The R83P indicator matrix:' )

  return
end
subroutine test11

!*****************************************************************************80
!
!! TEST11 tests R83P_ML.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) work2(n-1)
  real    ( kind = 8 ) work3(n-1)
  real    ( kind = 8 ) work4
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  R83P_ML computes A*x or A''*X'
  write ( *, '(a)' ) '    where A has been factored by R83P_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call r83p_random ( n, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r83p_mxv ( n, a, x, b )
    else
      call r83p_vxm ( n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r83p_fa ( n, a, info, work2, work3, work4 )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11 - Fatal error!'
      write ( *, '(a)' ) '  R83P_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r83p_ml ( n, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, i4_10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, i4_10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine test115

!*****************************************************************************80
!
!! TEST115 tests R85_INDICATOR.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(5,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115'
  write ( *, '(a)' ) '  R85_INDICATOR sets up an R85 indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r85_indicator ( n, a )

  call r85_print ( n, a, '  The R85 indicator matrix:' )

  return
end
subroutine test12

!*****************************************************************************80
!
!! TEST12 tests R85_NP_FS.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(5,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  R85_NP_FS factors and solves an R85'
  write ( *, '(a)' ) '    linear system, with no pivoting.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix to a random value.
!
  call r85_random ( n, seed, a )

  call r85_print ( n, a, '  The pentadiagonal matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute b = A * x.
!
  call r85_mxv ( n, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Wipe out the solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the system.
!
  call r85_np_fs ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution:' )

  return
end
subroutine test13

!*****************************************************************************80
!
!! TEST13 tests R8BB_FA, R8BB_PRINT, R8BB_RANDOM, R8BB_SL.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 8
  integer ( kind = 4 ), parameter :: n2 = 2
  integer ( kind = 4 ), parameter :: n = n1 + n2
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real    ( kind = 8 ) a(na)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  For a border banded matrix:'
  write ( *, '(a)' ) '  R8BB_FA factors;'
  write ( *, '(a)' ) '  R8BB_PRINT prints;'
  write ( *, '(a)' ) '  R8BB_RANDOM randomizes;'
  write ( *, '(a)' ) '  R8BB_SL solves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_random ( n1, n2, ml, mu, seed, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The border-banded matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n1+n2, x )
!
!  Compute the corresponding right hand side.
!
  call r8bb_mxv ( n1, n2, ml, mu, a, x, b )

  call r8vec_print ( n1+n2, b, '  The right hand side vector:' )
!
!  Factor the matrix.
!
  call r8bb_fa ( n1, n2, ml, mu, a, pivot, info )

  call r8bb_print ( n1, n2, ml, mu, a, '  The FACTORED border-banded matrix:' )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST13 - Fatal error!'
    write ( *, '(a)' ) '  R8BB_FA claims the matrix is singular.'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8bb_sl ( n1, n2, ml, mu, a, pivot, b )

  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine test14

!*****************************************************************************80
!
!! TEST14 tests R8BB_FA, R8BB_SL.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 0
  integer ( kind = 4 ), parameter :: n2 = 10
  integer ( kind = 4 ), parameter :: n = n1 + n2
  integer ( kind = 4 ), parameter :: ml = 0
  integer ( kind = 4 ), parameter :: mu = 0
  integer ( kind = 4 ), parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real    ( kind = 8 ) a(na)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  For a border banded matrix:'
  write ( *, '(a)' ) '  R8BB_FA factors;'
  write ( *, '(a)' ) '  R8BB_SL solves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_random ( n1, n2, ml, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8bb_mxv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8bb_fa ( n1, n2, ml, mu, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST14 - Fatal error!'
    write ( *, '(a)' ) '  R8BB_FA claims the matrix is singular.'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8bb_sl ( n1, n2, ml, mu, a, pivot, b )

  call r8vec_print ( n, b, '  Solution:' )
 
  return
end
subroutine test15

!*****************************************************************************80
!
!! TEST15 tests R8BB_FA, R8BB_SL.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 10
  integer ( kind = 4 ), parameter :: n2 = 0
  integer ( kind = 4 ), parameter :: n = n1 + n2
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real    ( kind = 8 ) a(na)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For a border banded matrix:'
  write ( *, '(a)' ) '  R8BB_FA factors;'
  write ( *, '(a)' ) '  R8BB_SL solves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_random ( n1, n2, ml, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8bb_mxv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix
!
  call r8bb_fa ( n1, n2, ml, mu, a, pivot, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Fatal error!'
    write ( *, '(a)' ) '  R8BB_FA claims the matrix is singular.'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8bb_sl ( n1, n2, ml, mu, a, pivot, b )

  call r8vec_print ( n, b, '  Solution:' )
 
  return
end
subroutine test151

!*****************************************************************************80
!
!! TEST151 tests R8BB_INDICATOR.
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 6
  integer ( kind = 4 ), parameter :: n2 = 2
  integer ( kind = 4 ), parameter :: n = n1 + n2
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real    ( kind = 8 ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST151'
  write ( *, '(a)' ) '  R8BB_INDICATOR sets up an R8BB indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB indicator matrix:' )

  return
end
subroutine test152

!*****************************************************************************80
!
!! TEST152 tests R8BLT_INDICATOR.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6
  integer ( kind = 4 ), parameter :: ml = 2

  real    ( kind = 8 ) a(ml+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST152'
  write ( *, '(a)' ) '  R8BLT_INDICATOR sets up an R8BLT indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
!
!  Set the matrix.
!
  call r8blt_indicator ( n, ml, a )

  call r8blt_print ( n, ml, a, '  The R8BLT indicator matrix:' )

  return
end
subroutine test153

!*****************************************************************************80
!
!! TEST153 tests R8BLT_MXV, R8BLT_PRINT, R8BLT_RANDOM, R8BLT_SL, R8BLT_VXM.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(ml+1,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST153'
  write ( *, '(a)' ) '  For a band matrix in lower triangular storage,'
  write ( *, '(a)' ) '  R8BLT_RANDOM sets a random value;'
  write ( *, '(a)' ) '  R8BLT_SL solves systems;'
  write ( *, '(a)' ) '  R8BLT_MXV computes matrix-vector products;'
  write ( *, '(a)' ) '  R8BLT_VXM computes vector-matrix products;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml

  call r8blt_random ( n, ml, seed, a )

  call r8blt_print ( n, ml, a, '  The R8BLT matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8blt_mxv ( n, ml, a, x, b )
    else
      call r8blt_vxm ( n, ml, a, x, b )
    end if

    call r8vec_print ( n, b, '  The right hand side:' )
!
!  Solve the linear system.
!
    call r8blt_sl ( n, ml, a, b, job )
 
    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution to the untransposed system:' )
    else
      call r8vec_print ( n, b, '  Solution to the transposed system:' )
    end if

  end do

  return
end
subroutine test154

!***********************************************************************
!
!! TEST154 tests R8BTO_INDICATOR.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  real    ( kind = 8 ), dimension ( m, m, 2*l-1 ) ::  a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST154'
  write ( *, '(a)' ) '  For a real block Toeplitz matrix,'
  write ( *, '(a)' ) '  R8BTO_INDICATOR sets up an indicator matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', m * l

  call r8bto_indicator ( m, l, a )

  call r8bto_print ( m, l, a, '  The block Toeplitz matrix:' )

  return
end
subroutine test155

!***********************************************************************
!
!! TEST155 tests R8BTO_MXV, R8BTO_VXM.
!
!  Modified:
!
!    24 March 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 3
  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ), parameter :: n = m * l
  integer ( kind = 4 ), parameter :: p = 2 * l - 1

  real    ( kind = 8 ), dimension ( m, m, p ) ::  a = reshape ( (/ &
    1.0D+00, 5.0D+00, 2.0D+00, 5.0D+00, &
    3.0D+00, 6.0D+00, 4.0D+00, 6.0D+00, &
    5.0D+00, 7.0D+00, 6.0D+00, 7.0D+00, &
    7.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
    9.0D+00, 9.0D+00, 0.0D+00, 9.0D+00 /), (/ m, m, p /) )
  real    ( kind = 8 ) b(m,l)
  real    ( kind = 8 ) x(m,l)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST155'
  write ( *, '(a)' ) '  For a real block Toeplitz matrix,'
  write ( *, '(a)' ) '  R8BTO_MXV computes A * x.'
  write ( *, '(a)' ) '  R8BTO_VXM computes A''* x.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8bto_print ( m, l, a, '  The block Toeplitz matrix:' )

  call dge_indicator ( m, l, x )

  call dge_print ( m, l, x, '  The "vector" x:' )

  call r8bto_mxv ( m, l, a, x, b )

  call dge_print ( m, l, b, '  The product A*x:' )

  call r8bto_vxm ( m, l, a, x, b )

  call dge_print ( m, l, b, '  The product A''*x:' )

  return
end
subroutine test156

!*****************************************************************************80
!
!! TEST156 tests R8BTO_SL.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: l = 3

  integer ( kind = 4 ), parameter :: n = m * l
  integer ( kind = 4 ), parameter :: p = 2 * l - 1

  real    ( kind = 8 ), dimension ( m, m, p ) ::  a = reshape ( (/ &
    9.0D+00, 2.0D+00, 1.0D+00, 8.0D+00, &
    3.0D+00, 6.0D+00, 4.0D+00, 6.0D+00, &
    5.0D+00, 7.0D+00, 6.0D+00, 7.0D+00, &
    7.0D+00, 8.0D+00, 8.0D+00, 8.0D+00, &
    9.0D+00, 9.0D+00, 0.0D+00, 9.0D+00 /), (/ m, m, p /) )
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST156'
  write ( *, '(a)' ) '  For a real block Toeplitz matrix,'
  write ( *, '(a)' ) '  R8BTO_SL solves a block Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Block order M =  ', m
  write ( *, '(a,i8)' ) '  Block number L = ', l
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8bto_print ( m, l, a, '  The block Toeplitz matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the right hand side.
!
  call r8bto_mxv ( m, l, a, x, b )

  call r8vec_print ( n, b, '  The right hand side B:' )

  call r8bto_sl ( m, l, a, b, x )

  call r8vec_print ( n, x, '  The computed solution X:' )

  return
end
subroutine test1565

!*****************************************************************************80
!
!! TEST1565 tests R8BUT_INDICATOR.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6
  integer ( kind = 4 ), parameter :: mu = 2

  real    ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1565'
  write ( *, '(a)' ) '  R8BUT_INDICATOR sets up an R8BUT indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8but_indicator ( n, mu, a )

  call r8but_print ( n, mu, a, '  The R8BUT indicator matrix:' )

  return
end
subroutine test1566

!*****************************************************************************80
!
!! TEST1566 tests R8BUT_MXV, R8BUT_PRINT, R8BUT_RANDOM, R8BUT_SL, R8BUT_VXM.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: mu = 3
  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(mu+1,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1566'
  write ( *, '(a)' ) '  For a band matrix in upper triangular storage,'
  write ( *, '(a)' ) '  R8BUT_RANDOM sets a random value;'
  write ( *, '(a)' ) '  R8BUT_SL solves systems;'
  write ( *, '(a)' ) '  R8BUT_MXV computes matrix-vector products;'
  write ( *, '(a)' ) '  R8BUT_VXM computes vector-matrix products;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8but_mxv ( n, mu, a, x, b )
    else
      call r8but_vxm ( n, mu, a, x, b )
    end if

    call r8vec_print ( n, b, '  The right hand side:' )
!
!  Solve the linear system.
!
    call r8but_sl ( n, mu, a, b, job )
 
    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution to the untransposed system:' )
    else
      call r8vec_print ( n, b, '  Solution to the transposed system:' )
    end if

  end do

  return
end
subroutine test157

!*****************************************************************************80
!
!! TEST157 tests R8CB_INDICATOR.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real    ( kind = 8 ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST157'
  write ( *, '(a)' ) '  For a compact band matrix:'
  write ( *, '(a)' ) '  R8CB_INDICATOR computes the indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB indicator matrix:' )

  return
end
subroutine test16

!*****************************************************************************80
!
!! TEST16 tests R8CB_NP_FA, R8CB_DET.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 3

  real    ( kind = 8 ) a(ml+mu+1,n)
  real    ( kind = 8 ) a2(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  For a compact band matrix, no pivoting:'
  write ( *, '(a)' ) '  R8CB_NP_FA factors;'
  write ( *, '(a)' ) '  R8CB_DET computes the determinant;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cb_random ( n, ml, mu, seed, a )

  call r8cb_print ( n, n, ml, mu, a, '  The compact band matrix:' )
!
!  Copy the matrix into a general array.
!
  call r8cb_to_dge ( n, ml, mu, a, a2 )
!
!  Factor the matrix.
!
  call r8cb_np_fa ( n, ml, mu, a, info )
!
!  Compute the determinant.
!
  call r8cb_det ( n, ml, mu, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8CB_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call dge_fa ( n, a2, pivot, info )
!
!  Compute the determinant.
!
  call dge_det ( n, a2, pivot, det )

  write ( *, '(a,g14.6)' ) '  DGE_DET computes the determinant = ', det

  return
end
subroutine test17

!*****************************************************************************80
!
!! TEST17 tests R8CB_NP_FA, R8CB_NP_SL.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 2

  real    ( kind = 8 ) a(ml+mu+1,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  For a compact band matrix, no pivoting:'
  write ( *, '(a)' ) '  R8CB_NP_FA factors;'
  write ( *, '(a)' ) '  R8CB_NP_SL solves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8cb_random ( n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the right hand side.
!
    if ( job == 0 ) then
      call r8cb_mxv ( n, ml, mu, a, x, b )
    else
      call r8cb_vxm ( n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8cb_np_fa ( n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Fatal error!'
      write ( *, '(a)' ) '  R8CB_NP_FA claims the matrix is singular.'
      write ( *, '(a,i8)' ) '  The value of info is ', info
      return
    end if
!
!  Solve the system.
!
    call r8cb_np_sl ( n, ml, mu, a, b, job )

    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution:' )
    else
      call r8vec_print ( n, b, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine test18

!*****************************************************************************80
!
!! TEST18 tests R8CB_ML.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 2

  real    ( kind = 8 ) a(ml+mu+1,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  For a compact band matrix:'
  write ( *, '(a)' ) '  R8CB_ML computes A*x or A''*X'
  write ( *, '(a)' ) '    where A has been factored by R8CB_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8cb_random ( n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8cb_mxv ( n, ml, mu, a, x, b )
    else
      call r8cb_vxm ( n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8cb_np_fa ( n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18 - Fatal error!'
      write ( *, '(a)' ) '  R8CB_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8cb_ml ( n, ml, mu, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, i4_10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, i4_10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine test19

!*****************************************************************************80
!
!! TEST19 tests R8CBB_FA, R8CBB_PRINT, R8CBB_RANDOM, R8CBB_SL.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 8
  integer ( kind = 4 ), parameter :: n2 = 2
  integer ( kind = 4 ), parameter :: n = n1 + n2
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real    ( kind = 8 ) a(na)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  For a compressed border banded matrix:'
  write ( *, '(a)' ) '  R8CBB_RANDOM randomly generates;'
  write ( *, '(a)' ) '  R8CBB_PRINT prints;'
  write ( *, '(a)' ) '  R8CBB_FA factors (no pivoting);'
  write ( *, '(a)' ) '  R8CBB_SL solves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cbb_random ( n1, n2, ml, mu, seed, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8cbb_mxv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix
!
  call r8cbb_fa ( n1, n2, ml, mu, a, info )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The factored R8CBB matrix:' )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST19 - Fatal error!'
    write ( *, '(a)' ) '  R8CBB_FA claims the matrix is singular.'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8vec_print ( n, b, '  The right hand side vector:' )

  call r8cbb_sl ( n1, n2, ml, mu, a, b )

  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine test193

!*****************************************************************************80
!
!! TEST193 tests R8CBB_INDICATOR.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 8
  integer ( kind = 4 ), parameter :: n2 = 2
  integer ( kind = 4 ), parameter :: n = n1 + n2
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1
  integer ( kind = 4 ), parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real    ( kind = 8 ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST193'
  write ( *, '(a)' ) '  For a compressed border banded matrix:'
  write ( *, '(a)' ) '  R8CBB_INDICATOR sets an indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cbb_indicator ( n1, n2, ml, mu, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The compact border-banded matrix:' )

  return
end
subroutine test195

!*****************************************************************************80
!
!! TEST195 tests R8CI_EVAL.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n)
  complex ( kind = 8 ) lambda(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST195'
  write ( *, '(a)' ) '  R8CI_EVAL finds the eigenvalues of '
  write ( *, '(a)' ) '  a real circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The circulant matrix:' )

  call r8ci_eval ( n, a, lambda )

  call c8vec_print ( n, lambda, '  The eigenvalues:' )

  return
end
subroutine test197

!*****************************************************************************80
!
!! TEST197 tests R8CI_INDICATOR.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST197'
  write ( *, '(a)' ) '  R8CI_INDICATOR sets up an R8CI indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print ( n, a, '  The circulant matrix:' )

  return
end
subroutine test20

!*****************************************************************************80
!
!! TEST20 tests R8CI_SL.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  R8CI_SL solves a circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The circulant matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8ci_mxv ( n, a, x, b )
    else
      call r8ci_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call r8ci_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call r8vec_print ( n, x, '  Solution:' )
    else
      call r8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine test21

!*****************************************************************************80
!
!! TEST21 tests R8GB_DET.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: mu = 2

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  real    ( kind = 8 ) a2(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8GB_DET computes the determinant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =       ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =    ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML  = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU  = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )

  call r8gb_print ( m, n, ml, mu, a, '  A random R8GB matrix:' )
!
!  Copy the matrix into a general array.
!
  call r8gb_to_dge ( m, n, ml, mu, a, a2 )
!
!  Factor the matrix.
!
  call r8gb_fa ( n, ml, mu, a, pivot, info )
!
!  Compute the determinant.
!
  call r8gb_det ( n, ml, mu, a, pivot, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8GB_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call dge_fa ( n, a2, pivot, info )
!
!  Compute the determinant.
!
  call dge_det ( n, a2, pivot, det )

  write ( *, '(a,g14.6)' ) '  DGE_DET computes the determinant = ', det

  return
end
subroutine test22

!*****************************************************************************80
!
!! TEST22 tests R8GB_FA, R8GB_SL.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: ml = 1 
  integer ( kind = 4 ), parameter :: mu = 2 

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8GB_FA computes the PLU factors.'
  write ( *, '(a)' ) '  R8GB_SL solves a factored linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )

  call r8gb_print ( m, n, ml, mu, a, '  The banded matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8gb_mxv ( m, n, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8gb_fa ( n, ml, mu, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST22 - Fatal error!'
    write ( *, '(a)' ) '  R8GB_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r8gb_sl ( n, ml, mu, a, pivot, b, job )
 
  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8gb_ml ( n, ml, mu, a, pivot, x, b, job )

  call r8vec_print ( n, b, '  Right hand side of transposed system:' )
!
!  Solve the linear system.
!
  job = 1
  call r8gb_sl ( n, ml, mu, a, pivot, b, job )
 
  call r8vec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine test23

!*****************************************************************************80
!
!! TEST23 tests R8GB_FA, R8GB_TRF.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 1

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8GB_FA factors, using LINPACK conventions;'
  write ( *, '(a)' ) '  R8GB_TRF factors, using LAPACK conventions;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  seed = 123456789
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Factor the matrix.
!
  call r8gb_fa ( n, ml, mu, a, pivot, info )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB_FA factors:' )
!
!  Set the matrix.
!
  seed = 123456789
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Factor the matrix.
!
  call r8gb_trf ( m, n, ml, mu, a, pivot, info )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB_TRF factors:')

  return
end
subroutine test235

!*****************************************************************************80
!
!! TEST235 tests R8GB_INDICATOR.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 3
  integer ( kind = 4 ), parameter :: mu = 2

  integer ( kind = 4 ), parameter :: row_num = 2 * ml + mu + 1

  real    ( kind = 8 ) a(row_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST235'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8GB_INDICATOR computes the indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =       ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =    ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML  = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU  = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call dge_print ( row_num, n, a, '  The banded matrix in DGE format:' )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB indicator matrix:' )

  return
end
subroutine test24

!*****************************************************************************80
!
!! TEST24 tests R8GB_ML.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 2

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8GB_ML computes A*x or A''*X'
  write ( *, '(a)' ) '    where A has been factored by R8GB_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8gb_mxv ( m, n, ml, mu, a, x, b )
    else
      call r8gb_vxm ( m, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8gb_fa ( n, ml, mu, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24 - Fatal error!'
      write ( *, '(a)' ) '  R8GB_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8gb_ml ( n, ml, mu, a, pivot, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, i4_10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, i4_10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine test25

!*****************************************************************************80
!
!! TEST25 tests R8GB_PRINT.
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 3

  real    ( kind = 8 ) a(2*ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8GB_PRINT prints the matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB matrix:' )

  return
end
subroutine test26

!*****************************************************************************80
!
!! TEST26 tests R8GB_NZ_NUM.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 2

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) diag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8GB_R8GB_NZ_NUM counts the nonzero entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Make some zero entries.
!
  do j = 1, n
    do diag = 1, 2*ml+mu+1
      if ( a(diag,j) < 0.3D+00 ) then
        a(diag,j) = 0.0D+00
      end if
    end do
  end do

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB matrix:' )

  call r8gb_nz_num ( m, n, ml, mu, a, nz_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Nonzero entries = ', nz_num

  return
end
subroutine test265

!*****************************************************************************80
!
!! TEST265 tests R8GB_TO_DGE, DGE_TO_R8GB.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  real    ( kind = 8 ) b(m,n)
  real    ( kind = 8 ) c(2*ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST265'
  write ( *, '(a)' ) '  R8GB_TO_DGE copies a R8GB matrix to a DGE matrix.'
  write ( *, '(a)' ) '  DGE_TO_R8GB copies a DGE matrix to a R8GB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB matrix:' )

  call r8gb_to_dge ( m, n, ml, mu, a, b )

  call dge_print ( m, n, b, '  The DGE matrix:' )

  call dge_to_r8gb ( m, n, ml, mu, b, c )

  call r8gb_print ( m, n, ml, mu, c, '  The recovered R8GB matrix:' )

  return
end
subroutine test2655

!*****************************************************************************80
!
!! TEST2655 tests R8GB_TO_DS3.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: ml = 2
  integer ( kind = 4 ), parameter :: mu = 1

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  real    ( kind = 8 ), allocatable, dimension ( : ) :: b
  integer ( kind = 4 ), allocatable, dimension ( : ) :: col
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: row

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2655'
  write ( *, '(a)' ) '  R8GB_TO_DS3 copies a R8GB matrix to a DS3 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB matrix:' )

  call r8gb_nz_num ( m, n, ml, mu, a, nz_num )

  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM =    ', nz_num

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( b(1:nz_num) )

  call r8gb_to_ds3 ( m, n, ml, mu, a, nz_num, isym, row, col, b )

  call ds3_print ( m, n, nz_num, isym, row, col, b, '  The DS3 matrix:' )

  deallocate ( row )
  deallocate ( col )
  deallocate ( b )

  return
end
subroutine test27

!*****************************************************************************80
!
!! TEST27 tests R8GB_TRF, R8GB_TRS.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: m = n
  integer ( kind = 4 ), parameter :: ml = 1
  integer ( kind = 4 ), parameter :: mu = 2
  integer ( kind = 4 ), parameter :: nrhs = 1

  real    ( kind = 8 ) a(2*ml+mu+1,n)
  real    ( kind = 8 ) b(n,nrhs)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  For a general banded matrix,'
  write ( *, '(a)' ) '  R8GB_TRF computes the PLU factors.'
  write ( *, '(a)' ) '  R8GB_TRS solves a factored linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8gb_mxv ( m, n, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8gb_trf ( m, n, ml, mu, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST27 - Fatal error!'
    write ( *, '(a)' ) '  R8GB_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8gb_trs ( n, ml, mu, nrhs, 'N', a, pivot, b, info )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8gb_mu ( n, ml, mu, a, pivot, x, b, job )
!
!  Solve the linear system.
!
  call r8gb_trs ( n, ml, mu, nrhs, 'T', a, pivot, b, info )

  call r8vec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine test275

!*****************************************************************************80
!
!! TEST275 tests DGD_INDICATOR.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real    ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ) offset(ndiag)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST275'
  write ( *, '(a)' ) '  For a general diagonal matrix:'
  write ( *, '(a)' ) '  DGD_INDICATOR sets up an indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call dgd_indicator ( n, ndiag, offset, a )

  call dgd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )

  return
end
subroutine test28

!*****************************************************************************80
!
!! TEST28 tests DGD_MXV, DGD_PRINT, DGD_RANDOM, DGD_VXM.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 4

  real    ( kind = 8 ) a(n,ndiag)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) offset(ndiag)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  For a general diagonal matrix:'
  write ( *, '(a)' ) '  DGD_MXV computes A * x;'
  write ( *, '(a)' ) '  DGD_PRINT prints it;'
  write ( *, '(a)' ) '  DGD_RANDOM randomly generates one;'
  write ( *, '(a)' ) '  DGD_VXM computes A''*x;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N            = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals NDIAG = ', ndiag
!
!  Set the matrix.
!
  offset(1:4) = (/ -2, 0, 1, n - 1 /)

  call i4vec_print ( ndiag, offset, '  The offset vector:' )

  call dgd_random ( n, ndiag, offset, seed, a )

  call dge_print ( n, ndiag, a, '  The raw matrix: ' )

  call dgd_print ( n, ndiag, offset, a, '  The general diagonal matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call dgd_mxv ( n, ndiag, offset, a, x, b )

  call r8vec_print ( n, b, '  A * x:' )
!
!  Compute the corresponding right hand side.
!
  call dgd_vxm ( n, ndiag, offset, a, x, b )

  call r8vec_print ( n, b, '  A'' * x:' )

  return
end
subroutine test285

!*****************************************************************************80
!
!! TEST285 tests DGE_CO.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) a_inverse(n,n)
  real    ( kind = 8 ) a_inverse_norm_l1
  real    ( kind = 8 ) a_norm_l1
  real    ( kind = 8 ) cond_l1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) pivot(n)
  real    ( kind = 8 ) rcond
  real    ( kind = 8 ) :: x = 2.0D+00
  real    ( kind = 8 ) :: y = 3.0D+00
  real    ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST285'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_CO estimates the condition number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = x + y
      else
        a(i,j) = y
      end if
    end do
  end do

  a_norm_l1 = 0.0D+00
  do j = 1, n
    a_norm_l1 = max ( a_norm_l1, sum ( abs ( a(1:n,j) ) ) )
  end do

  a_inverse(1:n,1:n) = a(1:n,1:n)
  call dge_fa ( n, a_inverse, pivot, info )
  call dge_inverse ( n, a_inverse, pivot )

  a_inverse_norm_l1 = 0.0D+00
  do j = 1, n
    a_inverse_norm_l1 = max ( a_inverse_norm_l1, &
      sum ( abs ( a_inverse(1:n,j) ) ) )
  end do

  cond_l1 = a_norm_l1 * a_inverse_norm_l1

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The L1 condition number is ', cond_l1
!
!  Factor the matrix.
!
  call dge_co ( n, a, pivot, rcond, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The DGE_CO estimate is     ', 1.0D+00 / rcond

  return
end
subroutine test29

!*****************************************************************************80
!
!! TEST29 tests DGE_DET.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) det
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) pivot(n)
  real    ( kind = 8 ) :: x = 2.0D+00
  real    ( kind = 8 ) :: y = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_DET computes the determinant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = x + y
      else
        a(i,j) = y
      end if
    end do
  end do
!
!  Factor the matrix.
!
  call dge_fa ( n, a, pivot, info )
!
!  Compute the determinant.
!
  call dge_det ( n, a, pivot, det )

  exact = x**(n-1) * ( x + real ( n ) * y )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  DGE_DET computes the determinant = ', det
  write ( *, '(a,g14.6)' ) '  Exact determinant =                ', exact

  return
end
subroutine test295

!*****************************************************************************80
!
!! TEST295 tests DGE_DILU.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 3
  integer ( kind = 4 ), parameter :: nrow = 3
  integer ( kind = 4 ), parameter :: n = nrow * ncol
  integer ( kind = 4 ), parameter :: m = n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) d(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST295'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_DILU returns the DILU factors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  do i = 1, nrow * ncol
    do j = 1, nrow * ncol

      if ( i == j ) then
        a(i,j) = 4.0D+00
      else if ( i == j + 1 .or. i == j - 1 .or. &
                i == j + nrow .or. i == j - nrow ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  call dge_print ( m, n, a, '  Matrix A:' )
!
!  Compute the incomplete LU factorization.
!
  call dge_dilu ( m, n, a, d )

  call r8vec_print ( m, d, '  DILU factor:' )

  return
end
subroutine test30

!*****************************************************************************80
!
!! TEST30 tests DGE_FA, DGE_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_FA computes the LU factors,'
  write ( *, '(a)' ) '  DGE_SL solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call dge_mxv ( n, n, a, x, b )
!
!  Factor the matrix.
!
  call dge_fa ( n, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST30 - Fatal error!'
    write ( *, '(a)' ) '  DGE_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call dge_sl ( n, a, pivot, b, job )
 
  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  job = 0
  call dge_ml ( n, a, pivot, x, b, job )
!
!  Solve the system
!
  job = 0
  call dge_sl ( n, a, pivot, b, job )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call dge_ml ( n, a, pivot, x, b, job )
!
!  Solve the system
!
  job = 1
  call dge_sl ( n, a, pivot, b, job )

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine test31

!*****************************************************************************80
!
!! TEST31 tests DGE_FA, DGE_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_FA computes the LU factors,'
  write ( *, '(a)' ) '  DGE_SL solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dge_random ( n, n, seed, a )

  call dge_print ( n, n, a, '  The matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call dge_mxv ( n, n, a, x, b )
!
!  Factor the matrix.
!
  call dge_fa ( n, a, pivot, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST31 - Fatal error!'
    write ( *, '(a)' ) '  DGE_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Display the gory details.
!
  call dmat_print ( n, n, a, '  The compressed LU factors:' )

  call i4vec_print ( n, pivot, '  The pivot vector P:' )
!
!  Solve the linear system.
!
  job = 0
  call dge_sl ( n, a, pivot, b, job )

  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine test315

!*****************************************************************************80
!
!! TEST315 tests DGE_ILU.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ncol = 3
  integer ( kind = 4 ), parameter :: nrow = 3
  integer ( kind = 4 ), parameter :: n = nrow * ncol
  integer ( kind = 4 ), parameter :: m = n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) l(m,m)
  real    ( kind = 8 ) lu(m,n)
  real    ( kind = 8 ) u(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST315'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_ILU returns the ILU factors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  do i = 1, nrow * ncol
    do j = 1, nrow * ncol

      if ( i == j ) then
        a(i,j) = 4.0D+00
      else if ( i == j + 1 .or. i == j - 1 .or. &
                i == j + nrow .or. i == j - nrow ) then
        a(i,j) = -1.0D+00
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  call dge_print ( m, n, a, '  Matrix A:' )
!
!  Compute the incomplete LU factorization.
!
  call dge_ilu ( m, n, a, l, u )

  call dge_print ( m, m, l, '  Factor L:' )

  call dge_print ( m, n, u, '  Factor U:' )

  lu(1:m,1:n) = matmul ( l(1:m,1:m), u(1:m,1:n) )

  call dge_print ( m, n, lu, '  Product L*U:' )

  return
end
subroutine test317

!*****************************************************************************80
!
!! TEST317 tests DGE_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST317'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_INDICATOR sets up the indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  call dge_indicator ( m, n, a )

  call dge_print ( m, n, a, '  The DGE indicator matrix:' )

  return
end
subroutine test32

!*****************************************************************************80
!
!! TEST32 tests DGE_NP_FA, DGE_NP_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_NP_FA computes the LU factors without pivoting,'
  write ( *, '(a)' ) '  DGE_NP_SL solves factored systems.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call dge_mxv ( n, n, a, x, b )
!
!  Factor the matrix.
!
  call dge_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST32 - Fatal error!'
    write ( *, '(a)' ) '  DGE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call dge_np_sl ( n, a, b, job )
 
  call r8vec_print_some ( n, b, i4_10, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call dge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 0
  call dge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call dge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 1
  call dge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine test33

!*****************************************************************************80
!
!! TEST33 tests DGE_NP_FA, DGE_NP_INVERSE.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_NP_FA computes LU factors without pivoting,'
  write ( *, '(a)' ) '  DGE_NP_INVERSE computes the inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dge_random ( n, n, seed, a )

  call dge_print ( n, n, a, '  The random matrix:' )
!
!  Factor and invert the matrix.
!
  b(1:n,1:n) = a(1:n,1:n)

  call dge_np_fa ( n, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST33 - Fatal error!'
    write ( *, '(a)' ) '  DGE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if

  call dge_np_inverse ( n, b )

  call dge_print ( n, n, b, '  The inverse matrix:' )
!
!  Compute A * B = C.
!
  call dge_mxm ( n, a, b, c )

  call dge_print ( n, n, c, '  The product:' )

  return
end
subroutine test34

!*****************************************************************************80
!
!! TEST34 tests DGE_FS.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_FS factors and solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call dge_mxv ( n, n, a, x, b )
!
!  Factor and solve the system.
!
  call dge_fs ( n, a, b, info )
  
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST34 - Fatal error!'
    write ( *, '(a)' ) '  DGE_FS reports the matrix is singular.'
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine test35

!*****************************************************************************80
!
!! TEST35 tests DGE_INVERSE.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) pivot(n)
  real    ( kind = 8 ), parameter :: x = 2.0D+00
  real    ( kind = 8 ), parameter :: y = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_INVERSE computes the inverse matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  do j = 1, n
    do i = 1, n
      if ( i == j ) then
        a(i,i) = x + y
      else
        a(i,j) = y
      end if
    end do
  end do

  call dge_print ( n, n, a, '  Matrix A:' )
!
!  Factor and invert the matrix.
!
  b(1:n,1:n) = a(1:n,1:n)

  call dge_fa ( n, b, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST35 - Fatal error!'
    write ( *, '(a)' ) '  DGE_FA reports the matrix is singular.'
    return
  end if

  call dge_inverse ( n, b, pivot )

  call dge_print ( n, n, b, '  Inverse matrix B:' )
!
!  Check.
!
  call dge_mxm ( n, a, b, c )

  call dge_print ( n, n, c, '  Product matrix:' )

  return
end
subroutine test36

!*****************************************************************************80
!
!! TEST36 tests DGE_ML.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_ML computes A*x or A''*X'
  write ( *, '(a)' ) '    where A has been factored by DGE_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call dge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call dge_mxv ( n, n, a, x, b )
    else
      call dge_vxm ( n, n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call dge_fa ( n, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST36 - Fatal error!'
      write ( *, '(a)' ) '  DGE_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call dge_ml ( n, a, pivot, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, i4_10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, i4_10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine test37

!*****************************************************************************80
!
!! TEST37 tests DGE_NP_ML.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_NP_ML computes A*x or A''*X'
  write ( *, '(a)' ) '    where A has been factored by DGE_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call dge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call dge_mxv ( n, n, a, x, b )
    else
      call dge_vxm ( n, n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call dge_np_fa ( n, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST37 - Fatal error!'
      write ( *, '(a)' ) '  DGE_NP_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      cycle
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call dge_np_ml ( n, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, i4_10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, i4_10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine test38

!*****************************************************************************80
!
!! TEST38 tests DGE_MU.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 3

  real    ( kind = 8 ) amn(m,n)
  real    ( kind = 8 ) anm(n,m)
  real    ( kind = 8 ) bm(m)
  real    ( kind = 8 ) bn(n)
  real    ( kind = 8 ) cm(m)
  real    ( kind = 8 ) cn(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(m+n)
  integer ( kind = 4 ) :: seed = 123456789
  character trans
  real    ( kind = 8 ) xm(m)
  real    ( kind = 8 ) xn(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST38'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_MU computes A*x or A''*X'
  write ( *, '(a)' ) '    where A has been factored by DGE_TRF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  do job = 0, 1

    if ( job == 0 ) then 
      trans = 'N'
    else
      trans = 'T'
    end if
!
!  Set the matrix.
!
    call dge_random ( m, n, seed, amn )

    if ( job == 0 ) then

      call r8vec_indicator ( n, xn )

      call dge_mxv ( m, n, amn, xn, cm )

    else

      call r8vec_indicator ( m, xm )

      call dge_vxm ( m, n, amn, xm, cn )

    end if
!
!  Factor the matrix.
!
    call dge_trf ( m, n, amn, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST38 - Fatal error!'
      write ( *, '(a)' ) '  DGE_TRF declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      cycle
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    if ( job == 0 ) then

      call dge_mu ( m, n, amn, trans, pivot, xn, bm )

      call r8vec2_print_some ( m, cm, bm, i4_10, '  A*x and PLU*x' )

    else

      call dge_mu ( m, n, amn, trans, pivot, xm, bn )

      call r8vec2_print_some ( n, cn, bn, i4_10, '  A''*x and (PLU)''*x' )

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Matrix is ', n, ' by ', m

  do job = 0, 1

    if ( job == 0 ) then 
      trans = 'N'
    else
      trans = 'T'
    end if
!
!  Set the matrix.
!
    call dge_random ( n, m, seed, anm )

    if ( job == 0 ) then

      call r8vec_indicator ( m, xm )

      call dge_mxv ( n, m, anm, xm, cn )

    else

      call r8vec_indicator ( n, xn )

      call dge_vxm ( n, m, anm, xn, cm )

    end if
!
!  Factor the matrix.
!
    call dge_trf ( n, m, anm, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST38 - Fatal error!'
      write ( *, '(a)' ) '  DGE_TRF declares the matrix is singular!'
      write ( *, '(a)' ) '  The value of INFO is ', info
      cycle
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    if ( job == 0 ) then

      call dge_mu ( n, m, anm, trans, pivot, xm, bn )

      call r8vec2_print_some ( n, cn, bn, i4_10, '  A*x and PLU*x' )

    else

      call dge_mu ( n, m, anm, trans, pivot, xn, bm )

      call r8vec2_print_some ( m, cm, bm, i4_10, '  A''*x and (PLU)''*x' )

    end if

  end do

  return
end
subroutine test385

!*****************************************************************************80
!
!! TEST385 tests DGE_PLU.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) l(m,m)
  real    ( kind = 8 ) p(m,m)
  real    ( kind = 8 ) plu(m,n)
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) u(m,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST385'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_PLU returns the PLU factors of a matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  call dge_random ( m, n, seed, a )

  call dge_print ( m, n, a, '  Matrix A:' )
!
!  Compute the PLU factors.
!
  call dge_plu ( m, n, a, p, l, u )

  call dge_print ( m, m, p, '  Factor P:' )

  call dge_print ( m, m, l, '  Factor L:' )

  call dge_print ( m, n, u, '  Factor U:' )

  plu(1:m,1:n) = matmul ( p(1:m,1:m), &
                 matmul ( l(1:m,1:m), u(1:m,1:n) ) )
        
  call dge_print ( m, n, plu, '  Product P*L*U:')

  return
end
subroutine test39

!*****************************************************************************80
!
!! TEST39 tests DGE_POLY.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) j
  real    ( kind = 8 ) p(0:n)
  real    ( kind = 8 ) p_true(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST39'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_POLY computes the characteristic polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  p_true(0:12) = &   
       (/ 1.0D+00,    - 23.0D+00,     231.0D+00,  - 1330.0D+00,    4845.0D+00, &
    - 11628.0D+00,   18564.0D+00, - 19448.0D+00,   12870.0D+00,  - 5005.0D+00, &
       1001.0D+00,    - 78.0D+00,       1.0D+00 /)
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n
      a(i,j) = real ( min ( i, j ), kind = 8 )
    end do
  end do
!
!  Get the characteristic polynomial.
!
  call dge_poly ( n, a, p )
!
!  Compare.
!
  call r8vec2_print_some ( n + i4_1, p, p_true, i4_10, 'I, P(I), True P(I)' )

  return
end
subroutine test40

!*****************************************************************************80
!
!! TEST40 tests DGE_SL_IT.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) alu(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  real    ( kind = 8 ) r(n)
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST40'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_SL_IT applies one step of iterative '
  write ( *, '(a)' ) '    refinement to an DGE_SL solution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the coefficient matrix.
!
  call hilbert_inverse ( n, a )
!
!  Set the right hand side b.
!
  b(1:n-1) = 0.0D+00
  b(n) = 1.0D+00
!
!  It is necessary to keep both an unfactored and factored copy
!  of the coefficient matrix.
!
  alu(1:n,1:n) = a(1:n,1:n)
!
!  Compute the factored coefficient matrix.
!
  call dge_fa ( n, alu, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST40 - Fatal error!'
    write ( *, '(a)' ) '  DGE_FA declares the matrix is singular!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!  (Careful!  DGE_SL overwrites the right hand side with the solution!)
!
  x(1:n) = b(1:n)

  call dge_sl ( n, alu, pivot, x, job )
!
!  Compute and print the residual.
!
  call dge_res ( n, n, a, x, b, r )

  call r8vec2_print_some ( n, x, r, i4_10, '  i, x, b-A*x' )
!
!  Take a few steps of iterative refinement.
!
  do j = 1, 5

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'Iterative refinement step ', j
    write ( *, '(a)' ) ' '
!
!  Improve the solution.
!
    call dge_sl_it ( n, a, alu, pivot, b, job, x, r )

    call r8vec_print_some ( n, r, i4_10, '  I, DX:' )
!
!  Compute and print the residual.
!
    call dge_res ( n, n, a, x, b, r )

    call r8vec2_print_some ( n, x, r, i4_10, '  i, x, b-A*x' )

  end do

  return
end
subroutine test41

!*****************************************************************************80
!
!! TEST41 tests DGE_TRF, DGE_TRS.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: m = n
  integer ( kind = 4 ), parameter :: nrhs = 1

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,nrhs)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) pivot(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST41'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_TRF computes the LU factors,'
  write ( *, '(a)' ) '  DGE_TRS solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Number of matrix columns N = ', n

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0D+00
      else if ( i == j - 1 ) then
        a(i,j) = - 1.0D+00
      else if ( i == j + 1 ) then
        a(i,j) = - 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call dge_trf ( m, n, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST41 - Fatal error!'
    write ( *, '(a)' ) '  DGE_TRF declares the matrix is singular!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  b(1:n-1,1) = 0.0D+00
  b(n,1) = n + 1

  call dge_trs ( n, nrhs, 'N', a, pivot, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST41 - Fatal error!'
    write ( *, '(a)' ) '  DGE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )

  b(1:n-1,1) = 0.0D+00
  b(n,1) = n + 1

  call dge_trs ( n, nrhs, 'T', a, pivot, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST41 - Fatal error!'
    write ( *, '(a)' ) '  DGE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine test42

!*****************************************************************************80
!
!! TEST42 tests DGE_NP_TRF, DGE_NP_TRM, DGE_NP_TRS.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = m
  integer ( kind = 4 ), parameter :: nrhs = 1

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) b(m,nrhs)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST42'
  write ( *, '(a)' ) '  For a matrix in general storage,'
  write ( *, '(a)' ) '  DGE_NP_TRF factors without pivoting,'
  write ( *, '(a)' ) '  DGE_NP_TRS solves factored systems.'
  write ( *, '(a)' ) '  DGE_NP_TRM computes A*X for factored A.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call dge_random ( m, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call dge_mxv ( m, n, a, x, b )
!
!  Factor the matrix.
!
  call dge_np_trf ( m, n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST42 - Fatal error!'
    write ( *, '(a)' ) '  DGE_NP_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call dge_np_trs ( n, nrhs, 'N', a, b, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST42 - Fatal error!'
    write ( *, '(a)' ) '  DGE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call dge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system
!
  call dge_np_trs ( n, nrhs, 'N', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST42 - Fatal error!'
    write ( *, '(a)' ) '  DGE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call dge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system.
!
  call dge_np_trs ( n, nrhs, 'T', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST42 - Fatal error!'
    write ( *, '(a)' ) '  DGE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine test422

!*****************************************************************************80
!
!! TEST422 tests R8CC_GET, R8CC_IJK, R8CC_INC, R8CC_KIJ, R8CC_SET.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real    ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 20
  real    ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST422'
  write ( *, '(a)' ) '  For a matrix in the R8CC format,'
  write ( *, '(a)' ) '  (double precision Harwell-Boeing Unsymmetric Assembled)'
  write ( *, '(a)' ) '  R8CC_GET gets an entry;'
  write ( *, '(a)' ) '  R8CC_IJK gets K from (I,J)'
  write ( *, '(a)' ) '  R8CC_INC increments an entry;'
  write ( *, '(a)' ) '  R8CC_KIJ gets (I,J) from K;'
  write ( *, '(a)' ) '  R8CC_SET sets an entry;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call i4vec_print ( n + i4_1, colptr, '  The COLPTR vector:' )

  call i4vec_print ( nz_num, rowind, '  The ROWIND vector:' )

  a(1:nz_num) = 0.0D+00
!
!  Initialize the matrix to random values.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The initial R8CC matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_IJK locates some (I,J) entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K'
  write ( *, '(a)' ) ' '

  do test = 1, test_num  
    i = i4_uniform ( 1, m, seed )
    j = i4_uniform ( 1, n, seed )
    call r8cc_ijk ( m, n, nz_num, colptr, rowind, i, j, k )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_KIJ locates some K entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         K         I         J'
  write ( *, '(a)' ) ' '

  do test = 1, test_num  
    k = i4_uniform ( 1, nz_num, seed )
    call r8cc_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) k, i, j
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_SET sets 10 entries at random.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K      NEW_VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform ( 1, nz_num, seed )
    call r8cc_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    value = 100.0D+00 + real ( test, kind = 8 )
    call r8cc_set ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_INC increments 10 entries at random.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K       NEW_VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform ( 1, nz_num, seed )
    call r8cc_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    value = 20.0D+00 + real ( test, kind = 8 )
    call r8cc_inc ( m, n, nz_num, colptr, rowind, a, i, j, value )
    call r8cc_get ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8CC_GET retrieves 10 entries.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K      VALUE'
  write ( *, '(a)' ) ' '

  do test = 1, 10
    k = i4_uniform ( 1, nz_num, seed )
    call r8cc_kij ( m, n, nz_num, colptr, rowind, k, i, j )
    call r8cc_get ( m, n, nz_num, colptr, rowind, a, i, j, value )
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) i, j, k, value
  end do

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The final R8CC matrix:' )

  return
end
subroutine test423

!*****************************************************************************80
!
!! TEST423 tests R8CC_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real    ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST423'
  write ( *, '(a)' ) '  For a matrix in the R8CC format,'
  write ( *, '(a)' ) '  (double precision Harwell-Boeing Unsymmetric Assembled)'
  write ( *, '(a)' ) '  R8CC_INDICATOR sets up the indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num

  call r8cc_indicator ( m, n, nz_num, colptr, rowind, a )

  call r8cc_print ( m, n, nz_num, colptr, rowind, a, &
    '  The R8CC indicator matrix:' )

  return
end
subroutine test425

!*****************************************************************************80
!
!! TEST425 tests R8CC_MXV, R8CC_VXM.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real    ( kind = 8 ) a(nz_num)
  real    ( kind = 8 ) bm(m)
  real    ( kind = 8 ) bn(n)
  real    ( kind = 8 ) c(m,n)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) xm(m)
  real    ( kind = 8 ) xn(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST425'
  write ( *, '(a)' ) '  For a matrix in the R8CC format,'
  write ( *, '(a)' ) '  (double precision Harwell-Boeing Unsymmetric Assembled)'
  write ( *, '(a)' ) '  R8CC_MXV multiplies an R8CC matrix by a vector;'
  write ( *, '(a)' ) '  R8CC_VXM multiplies a vector by an R8CC matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
!
!  Set the matrix.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )
!
!  Make an DGE copy.
!
  call r8cc_to_dge ( m, n, nz_num, colptr, rowind, a, c )
!
!  Print the DGE copy.
!
  call dge_print ( m, n, c, '  The R8CC matrix, in DGE form:' )
!
!  Compute BM = A(MxN) * XN
!
  xn(1) = 1.0D+00
  xn(2:n-1) = 0.0D+00
  xn(n) = -1.0D+00

  call r8vec_print ( n, xn, '  The vector xn:' )

  call r8cc_mxv ( m, n, nz_num, colptr, rowind, a, xn, bm )

  call r8vec_print ( m, bm, '  The product A * xn:' )
!
!  Compute BN = XM * A(MxN)
!
  xm(1) = 2.0D+00
  xm(2:m-1) = 0.0D+00
  xm(m) = -3.0D+00

  call R8CC_vxm ( m, n, nz_num, colptr, rowind, a, xm, bn )

  call r8vec_print ( n, bn, '  The product A'' * xm:' )

  return
end
subroutine test426

!*****************************************************************************80
!
!! TEST426 tests R8CC_PRINT.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 12

  real    ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension (n+1) :: colptr = (/ 1, 4, 6, 8, 10, 13 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: rowind = (/ &
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 /)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST426'
  write ( *, '(a)' ) '  For a matrix in the R8CC format,'
  write ( *, '(a)' ) '  (double precision Harwell-Boeing Unsymmetric Assembled)'
  write ( *, '(a)' ) '  R8CC_PRINT prints the matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M     = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N  = ', n
  write ( *, '(a,i8)' ) '  Nonzeros NZ_NUM   = ', nz_num
!
!  Set the matrix.
!
  call r8cc_random ( m, n, nz_num, colptr, rowind, a, seed )
!
!  Print the matrix.
!
  call r8cc_print ( m, n, nz_num, colptr, rowind, a, '  The R8CC matrix:' )

  return
end
subroutine test428

!*****************************************************************************80
!
!! TEST428 tests DLT_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST428'
  write ( *, '(a)' ) '  For a matrix in lower triangular storage,'
  write ( *, '(a)' ) '  DLT_INDICATOR sets up an indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  call dlt_indicator ( m, n, a )

  call dlt_print ( m, n, a, '  The SLT indicator matrix:' )

  return
end
subroutine test43

!*****************************************************************************80
!
!! TEST43 tests DLT_MXV, DLT_SL, DLT_VXM.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST43'
  write ( *, '(a)' ) '  For a matrix in lower triangular storage,'
  write ( *, '(a)' ) '  DLT_SL solves systems;'
  write ( *, '(a)' ) '  DLT_MXV computes matrix-vector products;'
  write ( *, '(a)' ) '  DLT_VXM computes vector-matrix products;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do i = 1, n
    do j = 1, n
      if ( j <= i ) then
        a(i,j) = j
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call dlt_print ( n, n, a, '  The lower triangular matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call dlt_mxv ( n, n, a, x, b )
    else
      call dlt_vxm ( n, n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call dlt_sl ( n, a, b, job )
 
    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution:' )
    else
      call r8vec_print ( n, b, '  Solution to the transposed system:' )
    end if

  end do

  return
end
subroutine test44

!*****************************************************************************80
!
!! TEST44 tests DLT_DET, DLT_INVERSE, DLT_MXM, DLT_RANDOM.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) c(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST44'
  write ( *, '(a)' ) '  For a matrix in lower triangular storage,'
  write ( *, '(a)' ) '  DLT_DET computes the determinant.'
  write ( *, '(a)' ) '  DLT_INVERSE computes the inverse.'
  write ( *, '(a)' ) '  DLT_MXM computes matrix products.'
  write ( *, '(a)' ) '  DLT_RANDOM sets a random value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dlt_random ( n, n, seed, a )

  b(1:n,1:n) = a(1:n,1:n)

  call dlt_print ( n, n, a, '  Matrix A:' )
!
!  Compute the determinant.
!
  call dlt_det ( n, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant is ', det
!
!  Compute the inverse matrix.
!
  call dlt_inverse ( n, b )

  call dlt_print ( n, n, b, '  Inverse matrix B:' )
!
!  Check
!
  call dlt_mxm ( n, a, b, c )

  call dlt_print ( n, n, c, '  Product A * B:' )

  return
end
subroutine test443

!*****************************************************************************80
!
!! TEST443 tests DNCF_INDICATOR.
!
!  Modified:
!
!    11 July 2007
!
!  Author:
!
!   John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_2 = 2
  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 15

  real    ( kind = 8 ), dimension ( nz_num ) :: a
  integer ( kind = 4 ), dimension ( i4_2, nz_num ) :: rowcol = reshape ( &
    (/ &
    1, 1, &
    2, 2, &
    3, 3, &
    4, 4, &
    5, 5, &
    2, 1, &
    5, 1, &
    1, 2, &
    5, 2, &
    1, 4, &
    2, 4, &
    3, 4, &
    4, 5, &
    4, 6, &
    1, 7 /), (/ i4_2, nz_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST443'
  write ( *, '(a)' ) '  DNCF_INDICATOR sets up a DNCF indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call dncf_indicator ( m, n, nz_num, rowcol, a )

  call dncf_print ( m, n, nz_num, rowcol, a, '  The DNCF indicator matrix:' )

  return
end
subroutine test445

!*****************************************************************************80
!
!! TEST445 tests DPBL_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: mu = 3

  real    ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST445'
  write ( *, '(a)' ) '  DPBL_INDICATOR sets up an DPBL indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth MU =   ', mu

  call dpbl_indicator ( n, mu, a )

  call dpbl_print ( n, mu, a, '  The DPBL indicator matrix:' )

  return
end
subroutine test45

!*****************************************************************************80
!
!! TEST45 tests DPBU_CG.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50
  integer ( kind = 4 ), parameter :: mu = 1

  real    ( kind = 8 ) a(mu+1,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) err
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_10 = 10
  real    ( kind = 8 ) r(n)
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST45'
  write ( *, '(a)' ) '  DPBU_CG applies the conjugate gradient method'
  write ( *, '(a)' ) '    to a symmetric positive definite banded '
  write ( *, '(a)' ) '    linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix values.
!
  a(2,1:n) = 2.0D+00
  a(1,2:n) = -1.0D+00

  call dpbu_print_some ( n, mu, a, i4_1, i4_1, i4_10, i4_10, &
    'The symmetric banded matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the right hand side.
!
  call dpbu_mxv ( n, mu, a, x, b )
!
!  Set the approximate solution.
!
  x(1:n) = 1.0D+00
!
!  Call the conjugate gradient method.
!
  call dpbu_cg ( n, mu, a, b, x )
!
!  Compute the residual, A*x-b
!
  call dpbu_mxv ( n, mu, a, x, r )
 
  err = maxval ( abs ( r(1:n) - b(1:n) ) )
 
  call r8vec_print_some ( n, x, i4_10, '  Solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum residual = ', err
 
  return
end
subroutine test46

!*****************************************************************************80
!
!! TEST46 tests DPBU_DET.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: mu = 3

  real    ( kind = 8 ) a(mu+1,n)
  real    ( kind = 8 ) a2(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST46'
  write ( *, '(a)' ) '  DPBU_DET, determinant of a positive definite'
  write ( *, '(a)' ) '    symmetric banded matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call dpbu_random ( n, mu, seed, a )

  call dpbu_print ( n, mu, a, '  The DPBU matrix:' )
!
!  Copy the matrix into a general array.
!
  call dpbu_to_dge ( n, mu, a, a2 )
!
!  Factor the matrix.
!
  call dpbu_fa ( n, mu, a, info )

  call dpbu_print ( n, mu, a, '  The DPBU factored matrix:' )
!
!  Compute the determinant.
!
  call dpbu_det ( n, mu, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  DPBU_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call dge_fa ( n, a2, pivot, info )
!
!  Compute the determinant.
!
  call dge_det ( n, a2, pivot, det )

  write ( *, '(a,g14.6)' ) '  DGE_DET computes the determinant =  ', det

  return
end
subroutine test47

!*****************************************************************************80
!
!! TEST47 tests DPBU_FA, DPBU_RANDOM, DPBU_SL.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50
  integer ( kind = 4 ), parameter :: mu = 1

  real    ( kind = 8 ) a(mu+1,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST47'
  write ( *, '(a)' ) '  For a banded positive definite symmetric matrix,'
  write ( *, '(a)' ) '  DPBU_FA computes the LU factors.'
  write ( *, '(a)' ) '  DPBU_RANDOM sets a random value.'
  write ( *, '(a)' ) '  DPBU_SL solves a linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
  write ( *, '(a)' ) ' '
!
!  Set the matrix values.
!
  call dpbu_random ( n, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the right hand side.
!
  call dpbu_mxv ( n, mu, a, x, b )
!
!  Factor the matrix.
!
  call dpbu_fa ( n, mu, a, info )
!
!  Solve the linear system.
!
  call dpbu_sl ( n, mu, a, b )
 
  call r8vec_print_some ( n, b, i4_10, '  Solution:' )
 
  return
end
subroutine test48

!*****************************************************************************80
!
!! TEST48 tests DPBU_ML.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: mu = 3

  real    ( kind = 8 ) a(mu+1,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST48'
  write ( *, '(a)' ) '  DPBU_ML computes A*x '
  write ( *, '(a)' ) '    where A has been factored by DPBU_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call dpbu_random ( n, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call dpbu_mxv ( n, mu, a, x, b )
!
!  Factor the matrix.
!
  call dpbu_fa ( n, mu, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a)' ) '  DPBU_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
  call dpbu_ml ( n, mu, a, x, b2 )

  call r8vec2_print_some ( n, b, b2, i4_10, '  A*x and PLU*x' )

  return
end
subroutine test485

!*****************************************************************************80
!
!! TEST485 tests DPBU_INDICATOR.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: mu = 3

  real    ( kind = 8 ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST485'
  write ( *, '(a)' ) '  DPBU_INDICATOR sets up an DPBU indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Bandwidth MU =   ', mu

  call dpbu_indicator ( n, mu, a )

  call dpbu_print ( n, mu, a, '  The DPBU indicator matrix:' )

  return
end
subroutine test49

!*****************************************************************************80
!
!! TEST49 tests DPBU_SOR.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50
  integer ( kind = 4 ), parameter :: mu = 1

  real    ( kind = 8 ) a(mu+1,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  real    ( kind = 8 ) eps
  real    ( kind = 8 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) itchk
  integer ( kind = 4 ) itknt
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) k
  real    ( kind = 8 ) omega
  real    ( kind = 8 ), parameter :: pi = 3.14159265D+00
  real    ( kind = 8 ) t
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST49'
  write ( *, '(a)' ) '  DPBU_SOR, SOR routine for iterative'
  write ( *, '(a)' ) '    solution of A*x=b.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do k = 1, 3
 
    if ( k == 1 ) then
      omega = 0.25D+00
    else if ( k == 2 ) then
      omega = 0.75D+00
    else
      omega = 1.00D+00
    end if
!
!  Set matrix values.
!
    a(2,1:n) = 2.0D+00

    a(1,1) = 0.0D+00
    a(1,2:n) = -1.0D+00
!
!  Set the desired solution.
!
    do i = 1, n
      t = pi * real ( i - 1 ) / real ( n - 1 )
      x(i) = sin ( t )
    end do
!
!  Compute the right hand side.
!
    call dpbu_mxv ( n, mu, a, x, b ) 
!
!  Set the initial solution estimate.
!
    x(1:n) = 1.0D+00
 
    itchk = 1
    itmax = 8000
    eps = 0.0001D+00

    call dpbu_sor ( n, mu, a, b, eps, itchk, itknt, itmax, omega, x )
!
!  Compute residual, A*x-b
!
    call dpbu_mxv ( n, mu, a, x, b2 )
 
    err = 0.0D+00
    do i = 1, n
      err = max ( err, abs ( b2(i) - b(i) ) )
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SOR iteration.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Relaxation factor OMEGA = ', omega
    write ( *, '(a,i8)' ) '  Iterations taken = ', itknt

    call r8vec_print_some ( n, x, i4_10, '  Solution:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Maximum error = ', err
 
  end do
 
  return
end
subroutine test50

!*****************************************************************************80
!
!! TEST50 tests DPO_FA, DPO_SL.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST50'
  write ( *, '(a)' ) '  DPO_FA factors a positive definite symmetric'
  write ( *, '(a)' ) '    linear system,'
  write ( *, '(a)' ) '  DPO_SL solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( min ( i, j ), kind = 8 )
    end do
  end do
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
  call dpo_mxv ( n, a, x, b )
!
!  Factor the matrix.
!
  call dpo_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a)' ) '  DPO_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call dpo_sl ( n, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  x(1:n) = 1
!
!  Compute the corresponding right hand side, using the factored matrix.
!
  call dpo_ml ( n, a, x, b )
!
!  Solve the linear system.
!
  call dpo_sl ( n, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine test505

!*****************************************************************************80
!
!! TEST505 tests DPO_FA;
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST505'
  write ( *, '(a)' ) '  DPO_FA factors a positive definite symmetric'
  write ( *, '(a)' ) '    linear system,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( min ( i, j ), kind = 8 )
    end do
  end do

  call dpo_print ( n, a, '  The matrix A:' )
!
!  Factor the matrix.
!
  call dpo_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a)' ) '  DPO_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if

  call dut_print ( n, n, a, '  The factor R (an DUT matrix):' )
!
!  Compute the product R' * R.
!
  a(1:n,1:n) = matmul ( transpose ( a(1:n,1:n) ), a(1:n,1:n) )

  call dge_print ( n, n, a, '  The product R'' * R:' )

  return
end
subroutine test51

!*****************************************************************************80
!
!! TEST51 tests DPO_DET, DPO_INVERSE.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) c(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST51'
  write ( *, '(a)' ) '  For a symmetric positive definite matrix'
  write ( *, '(a)' ) '    factored by DPO_FA,'
  write ( *, '(a)' ) '  DPO_DET computes the determinant;'
  write ( *, '(a)' ) '  DPO_INVERSE computes the inverse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  do i = 1, n
    do j = 1, n
      a(i,j) = real ( min ( i, j ), kind = 8 )
    end do
  end do

  b(1:n,1:n) = a(1:n,1:n)

  call dpo_print ( n, a, '  Matrix A:' )
!
!  Factor the matrix.
!
  call dpo_fa ( n, b, info )
!
!  Compute the determinant.
!
  call dpo_det ( n, b, det )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Matrix determinant = ', det
!
!  Compute the inverse.
!
  call dpo_inverse ( n, b )

  call dpo_print ( n, b, '  Inverse matrix B:' )
!
!  Check.
!
  call dpo_mxm ( n, a, b, c )

  call dpo_print ( n, c, '  Product A * B:' )

  return
end
subroutine test515

!*****************************************************************************80
!
!! TEST515 tests DPO_FA, DPO_SL, DPO_ML, DPO_MXV.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST515'
  write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
  write ( *, '(a)' ) '  DPO_FA computes the Cholesky factor.'
  write ( *, '(a)' ) '  DPO_SL solves a factored linear system.'
  write ( *, '(a)' ) '  DPO_MXV multiplies unfactored A * x'
  write ( *, '(a)' ) '  DPO_ML multiplies factored A * x'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dpo_random ( n, seed, a )

  call dpo_print ( n, a, '  The matrix A:' )
!
!  Compute the desired right hand side.
!
  call r8vec_indicator ( n, x )

  call dpo_mxv ( n, a, x, b )

  call r8vec_print ( n, b, '  Right hand side, computed by DPO_MXV' )
!
!  Factor the matrix.
!
  call dpo_fa ( n, a, info )
!
!  Solve the linear system.
!
  call dpo_sl ( n, a, b )

  call r8vec_print ( n, b, '  Solution (should be 1,2,3...)' )
!
!  Recompute the desired right hand side.
!  Since A has been factored, we have to call DPO_ML now.
!
  call r8vec_indicator ( n, x )

  call dpo_ml ( n, a, x, b )

  call r8vec_print ( n, b, '  Right hand side, computed by DPO_ML' )
!
!  Solve the linear system.
!
  call dpo_sl ( n, a, b )

  call r8vec_print ( n, b, '  Solution (should be 1,2,3...)' )

  return
end
subroutine test517

!*****************************************************************************80
!
!! TEST517 tests DPO_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST517'
  write ( *, '(a)' ) '  DPO_INDICATOR sets up an DPO indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dpo_indicator ( n, a )

  call dpo_print ( n, a, '  The DPO indicator matrix:' )
 
  return
end
subroutine test52

!*****************************************************************************80
!
!! TEST52 tests DPO_RANDOM, DPO_TO_DGE.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST52'
  write ( *, '(a)' ) '  DPO_RANDOM computes a random positive definite'
  write ( *, '(a)' ) '    symmetric matrix.'
  write ( *, '(a)' ) '  DPO_TO_DGE converts an DPO matrix to DGE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dpo_random ( n, seed, a )

  call dpo_print ( n, a, '  The random DPO matrix:' )
 
  call dge_print ( n, n, a, '  The random DPO matrix (printed by DGE_PRINT):' )

  call dpo_to_dge ( n, a, b )

  call dge_print ( n, n, b, '  The random DGE matrix (printed by DGE_PRINT):' )

  return
end
subroutine test525

!*****************************************************************************80
!
!! TEST525 tests DPP_FA, DPP_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a((n*(n+1))/2)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST525'
  write ( *, '(a)' ) '  DPP_FA factors an DPP system,'
  write ( *, '(a)' ) '  DPP_SL solves an DPP system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dpp_random ( n, seed, a )

  call dpp_print ( n, a, '  The DPP matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )

  call r8vec_print ( n, x, '  The desired solution:' )
!
!  Compute the corresponding right hand side.
!
  call dpp_mxv ( n, a, x, b )

  call r8vec_print ( n, b, '  The right hand side:' )
!
!  Factor the matrix.
!
  call dpp_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a)' ) '  DPP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The DPP matrix has been factored.'
  end if
!
!  Solve the linear system.
!
  call dpp_sl ( n, a, b )
 
  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine test527

!*****************************************************************************80
!
!! TEST527 tests DPP_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a((n*(n+1))/2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST527'
  write ( *, '(a)' ) '  DPP_INDICATOR sets up an DPP indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dpp_indicator ( n, a )

  call dpp_print ( n, a, '  The DPP indicator matrix:' )
 
  return
end
subroutine test53

!*****************************************************************************80
!
!! TEST53 tests DPP_RANDOM.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a((n*(n+1))/2)
  real    ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST53'
  write ( *, '(a)' ) '  DPP_RANDOM, compute a random positive definite'
  write ( *, '(a)' ) '    symmetric packed matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dpp_random ( n, seed, a )

  call dpp_print ( n, a, '  The matrix (printed by DPP_PRINT):' )
 
  call dpp_to_dge ( n, a, b )

  call dge_print ( n, n, b, '  The random DPP matrix (printed by DGE_PRINT):' )

  return
end
subroutine test531

!*****************************************************************************80
!
!! TEST531 tests DS3_INDICATOR.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 6
  integer ( kind = 4 ), parameter :: nz_num = 20

  real    ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    5, 6, 2, 2, 3, 4, 4, 5, 1, 6, &
    4, 6, 5, 1, 6, 3, 1, 2, 1, 3 /)
  integer ( kind = 4 ) isym
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 3, 4, 6, 5, 2, 6, 3, 1, 2, &
    4, 6, 5, 4, 4, 3, 6, 2, 3, 4 /)

  isym = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST531'
  write ( *, '(a)' ) '  DS3_INDICATOR sets up an DS3 indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =  ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros = ', nz_num

  call ds3_indicator ( n, nz_num, isym, row, col, a )

  call ds3_print ( m, n, nz_num, isym, row, col, a, &
    '  The DS3 indicator matrix:' )

  return
end
subroutine test532

!*****************************************************************************80
!
!! TEST532 tests DS3_DIAGONAL.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 6
  integer ( kind = 4 ), parameter :: nz_num = 20

  real    ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    5, 6, 2, 2, 3, 4, 4, 5, 1, 6, &
    4, 6, 5, 1, 6, 3, 1, 2, 1, 3 /)
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 3, 4, 6, 5, 2, 6, 3, 1, 2, &
    4, 6, 5, 4, 4, 3, 6, 2, 3, 4 /)

  isym = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST532'
  write ( *, '(a)' ) '  DS3_DIAGONAL rearranges a square DS3 matrix'
  write ( *, '(a)' ) '  structure so that the diagonal is listed first.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =  ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros = ', nz_num

  call ds3_indicator ( n, nz_num, isym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Before rearrangement:'
  write ( *, '(a)' ) '       K  ROW(K)  COL(K)      A(K)'
  write ( *, '(a)' ) ' '
  do k = 1, nz_num
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) k, row(k), col(k), a(k)
  end do

  call ds3_diagonal ( n, nz_num, isym, row, col, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After rearrangement:'
  write ( *, '(a)' ) '       K  ROW(K)  COL(K)      A(K)'
  write ( *, '(a)' ) ' '
  do k = 1, nz_num
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) k, row(k), col(k), a(k)
  end do

  return
end
subroutine test533

!*****************************************************************************80
!
!! TEST533 tests DS3_JAC_SL.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: nz_num = 3 * n - 2

  real    ( kind = 8 ) a(nz_num)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) col(nz_num)
  real    ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nz_num)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) x(n)

  it_max = 1000
  tol = 0.000001D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST533'
  write ( *, '(a)' ) '  For a DS3 system,'
  write ( *, '(a)' ) '  DS3_JAC_SL solves a linear system using'
  write ( *, '(a)' ) '    Jacobi iteration'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =      ', n
  write ( *, '(a,i8)' ) '  Iterations per call = ', it_max

  isym = 0
!
!  Set the matrix values.
!
  k = 0
  do i = 1, n
    k = k + 1
    row(k) = i
    col(k) = i
    a(k) = 2.0D+00    
  end do

  do i = 2, n
    k = k + 1
    row(k) = i
    col(k) = i-1
    a(k) = -1.0D+00    
  end do

  do i = 1, n-1
    k = k + 1
    row(k) = i
    col(k) = i+1
    a(k) = -1.0D+00    
  end do

  do job = 0, 1

    if ( job == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solving A * x = b.'
      write ( *, '(a)' ) ' '
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solving A'' * x = b.'
      write ( *, '(a)' ) ' '
    end if
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call ds3_mxv ( n, n, nz_num, isym, row, col, a, x, b )
    else
      call ds3_vxm ( n, n, nz_num, isym, row, col, a, x, b )
    end if

    call r8vec_print_some ( n, b, i4_10, '  The right hand side:' )
!
!  Set the starting solution.
!
    x(1:n) = 0.0D+00
!
!  Solve the linear system.
!
    do i = 1, 3

      call ds3_jac_sl ( n, nz_num, isym, row, col, a, b, x, tol, it_max, &
        job, it, diff )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
      write ( *, '(a,g14.6)' ) '  Maximum solution change on last step = ', diff

      call r8vec_print_some ( n, x, i4_10, '  Current solution estimate:' )

    end do

  end do

  return
end
subroutine test534

!*****************************************************************************80
!
!! TEST534 tests DS3_WRITE.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: nz_num = 3 * n - 2

  real    ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 80 ) output_file
  integer ( kind = 4 ) row(nz_num)

  output_file = 'ds3_matrix.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST534'
  write ( *, '(a)' ) '  For a DS3 system,'
  write ( *, '(a)' ) '  DS3_WRITE writes the matrix to a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
  write ( *, '(a)' ) ' '

  isym = 0
!
!  Set the matrix values.
!
  k = 0
  do i = 1, n
    k = k + 1
    j = i
    row(k) = i
    col(k) = j
    a(k) = real ( 100 * i + j, kind = 8 )
  end do

  do i = 2, n
    j = i - 1
    k = k + 1
    row(k) = i
    col(k) = j
    a(k) = real ( 100 * i + j, kind = 8 )
  end do

  do i = 1, n-1
    j = i + 1
    k = k + 1
    row(k) = i
    col(k) = j
    a(k) = real ( 100 * i + j, kind = 8 )
  end do

  call ds3_print_some ( n, n, nz_num, isym, row, col, a, i4_1, i4_1, &
    i4_10, i4_10, '  Initial 10x10 block of DS3 matrix:' )

  call ds3_write ( n, nz_num, isym, row, col, a, output_file )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DS3_WRITE wrote the matrix data to "' &
    // trim ( output_file ) // '".'

  return
end
subroutine test535

!*****************************************************************************80
!
!! TEST535 tests DS3_READ, DS3_READ_SIZE.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ), allocatable, dimension ( : ) :: a
  integer ( kind = 4 ), allocatable, dimension ( : ) :: col
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_10 = 10
  character ( len = 80 ) input_file
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: row

  input_file = 'ds3_matrix.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST535'
  write ( *, '(a)' ) '  For a DS3 system,'
  write ( *, '(a)' ) '  DS3_READ_SIZE reads the size of the matrix.'
  write ( *, '(a)' ) '  DS3_READ reads the matrix.'

  call ds3_read_size ( input_file, n, nz_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DS3_READ_SIZE reports matrix size data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call ds3_read ( input_file, n, nz_num, row, col, a )

  isym = 0

  call ds3_print_some ( n, n, nz_num, isym, row, col, a, i4_1, i4_1, &
    i4_10, i4_10, '  Initial 10x10 block of recovered DS3 matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Deleting the matrix data file "' &
    // trim ( input_file ) // '".'

  call file_delete ( input_file )

  deallocate ( row )
  deallocate ( col )
  deallocate ( a )

  return
end
subroutine test54

!*****************************************************************************80
!
!! TEST54 tests DSD_CG.
!
!  Discussion:
!
!    NX and NY are the number of grid points in the X and Y directions.
!    N is the number of unknowns.
!    NDIAG is the number of nonzero diagonals we will store.  We only
!      store the main diagonal, and the superdiagonals.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_0 = 0
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ), parameter :: ndiag = 3
  integer ( kind = 4 ), parameter :: nx = 10
  integer ( kind = 4 ), parameter :: ny = 10
  integer ( kind = 4 ), parameter :: n = nx * ny

  real    ( kind = 8 ) a(n,ndiag)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  real    ( kind = 8 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ i4_0, i4_1, nx /)
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST54'
  write ( *, '(a)' ) '  DSD_CG applies the conjugate gradient method'
  write ( *, '(a)' ) '    to a symmetric positive definite linear'
  write ( *, '(a)' ) '    system stored by diagonals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a,i8)' ) '  Number of diagonals is ', ndiag
  write ( *, '(a)' ) ' '
!
!  Now we compute the numbers that go into the diagonals.  For this
!  problem, we could simply store a column of 4's, and two columns of
!  -1's, but I wanted to go through the motions of thinking about the
!  value of each entry.  "K" counts the row of the original matrix
!  that we are working on.
!
  k = 0
  do j = 1, ny
    do i = 1, nx

      k = k + 1
!
!  Central
!
      a(k,1) = 4.0D+00
!
!  East ( = West )
!
      if ( i == nx ) then
        a(k,2) = 0.0D+00
      else
        a(k,2) = -1.0D+00
      end if
!
!  North ( = South )
!
      if ( j == ny ) then
        a(k,3) = 0.0D+00
      else
        a(k,3) = -1.0D+00
      end if

    end do
  end do
!
!  Print some of the matrix.
!
  call dsd_print_some ( n, ndiag, offset, a, i4_1, i4_1, i4_10, i4_10, &
    '  First 10 rows and columns of matrix.' )
!
!  Set the desired solution.
!
  k = 0
  do j = 1, ny
    do i = 1, nx
      k = k + 1
      x(k) = real ( 10 * i + j, kind = 8 )
    end do
  end do
!
!  Compute the corresponding right hand side.
!
  call dsd_mxv ( n, ndiag, offset, a, x, b )

  call r8vec_print_some ( n, b, i4_10, '  Right hand side:' )
!
!  Set X to zero so no one accuses us of cheating.
!
  x(1:n) = 0.0D+00
!
!  Call the conjugate gradient method.
!
  call dsd_cg ( n, ndiag, offset, a, b, x )
!
!  Compute the residual, A*x-b
!
  call dsd_mxv ( n, ndiag, offset, a, x, b2 )
 
  err = maxval ( abs ( b2(1:n) - b(1:n) ) )
 
  call r8vec_print_some ( n, x, i4_10, '  Solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum residual = ', err
!
!  Note that if we're not satisfied with the solution, we can
!  call again, using the computed X as our starting estimate.
!
!
!  Call the conjugate gradient method AGAIN.
!
  call dsd_cg ( n, ndiag, offset, a, b, x )
!
!  Compute the residual, A*x-b
!
  call dsd_mxv ( n, ndiag, offset, a, x, b2 )
 
  err = maxval ( abs ( b2(1:n) - b(1:n) ) )
 
  call r8vec_print_some ( n, x, i4_10, '  Second attempt at solution:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum residual of second attempt = ', err

  return
end
subroutine test55

!*****************************************************************************80
!
!! TEST55 tests DSD_CG.
!
!  Discussion:
!
!    This is a sample demonstration of how to compute some eigenvalues
!    and corresponding eigenvectors of a matrix.  The matrix is the
!    discretized Laplacian operator, which can be stored by diagonals,
!    and handled by the conjugate gradient method.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxvec = 3
  integer ( kind = 4 ), parameter :: ndiag = 3
  integer ( kind = 4 ), parameter :: nx = 10
  integer ( kind = 4 ), parameter :: ny = 10
  integer ( kind = 4 ), parameter :: n = nx * ny
  real    ( kind = 8 ), parameter :: pi = 3.141592653589D+00

  real    ( kind = 8 ) a(n,ndiag)
  real    ( kind = 8 ) del
  real    ( kind = 8 ) dot
  real    ( kind = 8 ) eval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) ivec
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) lambda
  real    ( kind = 8 ) lambda_old
  real    ( kind = 8 ) lamvec(maxvec)
  real    ( kind = 8 ) norm
  integer ( kind = 4 ) nvec
  integer ( kind = 4 ) offset(ndiag)
  real    ( kind = 8 ) vec(n,maxvec)
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xnew(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST55'
  write ( *, '(a)' ) '  DSD_CG is used for linear equation solving'
  write ( *, '(a)' ) '    in a demonstration of inverse iteration to'
  write ( *, '(a)' ) '    compute eigenvalues and eigenvectors of a'
  write ( *, '(a)' ) '    symmetric matrix stored by diagonals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here are 25 of the smallest eigenvalues:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, J, eigenvalue(I,J):'
  write ( *, '(a)' ) ' '

  do i = 1, min ( 5, nx )
    do j = 1, min ( 5, ny )
      eval = 4.0D+00 - 2.0D+00 * cos ( real ( i ) * pi / real ( nx + 1 ) ) &
                     - 2.0D+00 * cos ( real ( j ) * pi / real ( ny + 1 ) )
      write ( *, '(2i8,g14.6)' ) i, j, eval
    end do
  end do
!
!  OFFSET tells us where the nonzero diagonals are.  It does this
!  by recording how "high" or to the right the diagonals are from
!  the main diagonal.
!
  offset(1) =   0
  offset(2) =   1
  offset(3) =  nx
!
!  Now we compute the numbers that go into the diagonals.  For this
!  problem, we could simply store a column of 4's, and two columns of
!  -1's, but I wanted to go through the motions of thinking about the
!  value of each entry.  "K" counts the row of the original matrix
!  that we are working on.
!
  k = 0
  do j = 1, ny
    do i = 1, nx

      k = k + 1
!
!  Central
!
      a(k,1) = 4.0D+00
!
!  East ( = West )
!
      if ( i == nx ) then
        a(k,2) = 0.0D+00
      else
        a(k,2) = -1.0D+00
      end if
!
!  North ( = South )
!
      if ( j == ny ) then
        a(k,3) = 0.0D+00
      else
        a(k,3) = -1.0D+00
      end if

    end do
  end do

  nvec = 0
!
!  Set the starting eigenvector and eigenvalue estimates.
!
  do 

    write ( *, '(a)' ) ' '

    lambda = 0.0D+00

    k = 0
    do j = 1, ny
      do i = 1, nx
        k = k + 1
        x(k) = 1.0D+00
      end do
    end do
!
!  Remove any components of previous eigenvectors.
!
    do ivec = 1, nvec
      dot = dot_product ( x(1:n), vec(1:n,ivec) )
      x(1:n) = x(1:n) - dot * vec(1:n,ivec)
    end do

    xnew(1:n) = x(1:n)
!
!  Iterate
!
    do iter = 1, 40

      norm = sqrt ( sum ( xnew(1:n)**2 ) )

      xnew(1:n) = xnew(1:n) / norm

      lambda_old = lambda
      lambda = 1.0D+00 / norm
!
!  Check for convergence.
!
      if ( 1 < iter ) then
        del = abs ( lambda - lambda_old )
        if ( del < 0.000001D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,g14.6)' ) 'Lambda estimate = ', lambda
          write ( *, '(a,i8)' ) 'Converged on step ', iter
          exit
        end if
      end if
!
!  Call the conjugate gradient method, solving
!    A * XNEW = X.
!
      x(1:n) = xnew(1:n)

      call dsd_cg ( n, ndiag, offset, a, x, xnew  )

      do ivec = 1, nvec
        dot = dot_product ( xnew(1:n), vec(1:n,ivec) )
        xnew(1:n) = xnew(1:n) - dot * vec(1:n,ivec)
      end do

    end do

    nvec = nvec + 1
    lamvec(nvec) = lambda
    vec(1:n,nvec) = xnew(1:n)

    if ( maxvec <= nvec ) then
      exit
    end if

  end do

  return
end
subroutine test555

!*****************************************************************************80
!
!! TEST555 tests DSD_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: ndiag = 3

  real    ( kind = 8 ) a(n,ndiag)
  integer ( kind = 4 ), dimension ( ndiag ) :: offset = (/ 0, 1, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST555'
  write ( *, '(a)' ) '  DSD_INDICATOR sets up an DSD indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =         ', n
  write ( *, '(a,i8)' ) '  Matrix diagonals NDIAG = ', ndiag

  call dsd_indicator ( n, ndiag, offset, a )

  call dsd_print ( n, ndiag, offset, a, '  The DSD indicator matrix:' )

  return
end
subroutine test56

!*****************************************************************************80
!
!! TEST56 tests DSM_ML.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = m

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) u(m)
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST56'
  write ( *, '(a)' ) '  DSM_ML computes A*x or A''*X'
  write ( *, '(a)' ) '    where A is a Sherman Morrison matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call dsm_random ( m, n, seed, a, u, v )

    call dsm_print ( m, n, a, u, v, '  The Sherman Morrison matrix:' )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call dsm_mxv ( m, n, a, u, v, x, b )
    else
      call dsm_vxm ( m, n, a, u, v, x, b )
    end if
!
!  Factor the matrix.
!
    call dge_fa ( n, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Fatal error!'
      write ( *, '(a)' ) '  DGE_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call dsm_ml ( n, a, u, v, pivot, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, i4_10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, i4_10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine test57

!*****************************************************************************80
!
!! TEST57 tests DSM_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = m

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) u(m)
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST57'
  write ( *, '(a)' ) '  DSM_SL implements the Sherman-Morrison method '
  write ( *, '(a)' ) '    for solving a perturbed linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  do job = 1, 0, -1
!
!  Set the matrix.
!
    call dsm_random ( m, n, seed, a, u, v )

    call dsm_print ( m, n, a, u, v, '  The Sherman-Morrison matrix A:' )
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call dsm_mxv ( m, n, a, u, v, x, b )
    else
      call dsm_vxm ( m, n, a, u, v, x, b )
    end if

    call r8vec_print ( n, b, '  The right hand side vector B:' )
!
!  Factor the matrix.
!
    call dge_fa ( n, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Fatal error!'
      write ( *, '(a)' ) '  DGE_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      cycle
    end if
!
!  Solve the linear system.
!
    call dsm_sl ( n, a, u, v, b, ierror, pivot, job )
 
    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution to A * X = B:' )
    else
      call r8vec_print ( n, b, '  Solution to A'' * X = B:' )
    end if
 
  end do

  return
end
subroutine test5705

!*****************************************************************************80
!
!! TEST5705 tests DSP_IJ_TO_K.
!
!  Modified:
!
!    11 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 10

  logical check
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    2, 5, 1, 5, 1, 2, 3, 4, 4, 1 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 1, 2, 2, 4, 4, 4, 5, 6, 7 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5705'
  write ( *, '(a)' ) '  DSP_IJ_TO_K returns the DSP index of (I,J).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call dsp_check ( m, n, nz_num, row, col, check )

  if ( .not. check ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DSP_CHECK - Error!'
    write ( *, '(a)' ) '  The matrix is not in the proper sorted format.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J         K'
  write ( *, '(a)' ) ' '
  do i = 1, m
    do j = 1, n
      call dsp_ij_to_k ( nz_num, row, col, i, j, k )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
    end do
  end do

  return
end
subroutine test571

!*****************************************************************************80
!
!! TEST571 tests DSP_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 10

  real    ( kind = 8 ), dimension ( nz_num ) :: a
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    2, 5, 1, 5, 1, 2, 3, 4, 4, 1 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 1, 2, 2, 4, 4, 4, 5, 6, 7 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST571'
  write ( *, '(a)' ) '  DSP_INDICATOR sets up a DSP indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call dsp_indicator ( m, n, nz_num, row, col, a )

  call dsp_print ( m, n, nz_num, row, col, a, '  The DSP indicator matrix:' )

  return
end
subroutine test572

!*****************************************************************************80
!
!! TEST572 tests DSP_MXV, DSP_VXM.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 7
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz_num = 10

  real    ( kind = 8 ), dimension ( nz_num ) :: a
  real    ( kind = 8 ) b(m)
  real    ( kind = 8 ) c(m,n)
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    2, 5, 1, 5, 1, 2, 3, 4, 4, 1 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    1, 1, 2, 2, 4, 4, 4, 5, 6, 7 /)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST572'
  write ( *, '(a)' ) '  DSP_MXV multiplies an DSP matrix by a vector;'
  write ( *, '(a)' ) '  DSP_VXM multiplies a vector by an DSP matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num
!
!  Set the matrix.
!
  call dsp_random ( m, n, nz_num, row, col, seed, a )
!
!  Make an DGE copy.
!
  call dsp_to_dge ( m, n, nz_num, row, col, a, c )
!
!  Print the DGE copy.
!
  call dge_print ( m, n, c, '  The DSP matrix, in DGE form:' )

  x(1) = 1.0D+00
  x(2:n-1) = 0.0D+00
  x(n) = -1.0D+00

  call r8vec_print ( n, x, '  The vector x:' )

  call dsp_mxv ( m, n, nz_num, row, col, a, x, b )

  call r8vec_print ( m, b, '  The product A * x:' )

  b(1) = 1.0D+00
  b(2:m-1) = 0.0D+00
  b(m) = -1.0D+00

  call r8vec_print ( m, b, '  The vector x:' )

  call dsp_vxm ( m, n, nz_num, row, col, a, b, x )

  call r8vec_print ( n, x, '  The product A'' * x:' )

  return
end
subroutine test5722

!*****************************************************************************80
!
!! TEST5722 tests DSP_PRINT.
!
!  Discussion:
!
!    Because MATLAB seems to allow a DSP matrix to store the same index
!    several times, presumably with the matrix entry being the SUM of
!    these occurrences, I modified DSP_PRINT to handle this situation
!    (I hope).
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 7
  integer ( kind = 4 ), parameter :: nz_num = 12

  real    ( kind = 8 ), dimension ( nz_num ) :: a = (/ &
    21.0D+00,  51.0D+00, 12.0D+00, 52.0D+00, 14.0D+00, &
    24.0D+00,  34.0D+00, 45.0D+00, 46.0D+00, 17.0D+00, &
   100.0D+00, 200.0D+00 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: col = (/ &
    1, 1, 2, 2, 4, 4, 4, 5, 6, 7, 2, 4 /)
  integer ( kind = 4 ), dimension ( nz_num ) :: row = (/ &
    2, 5, 1, 5, 1, 2, 3, 4, 4, 1, 1, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5722'
  write ( *, '(a)' ) '  DSP_PRINT prints a DSP matrix;'
  write ( *, '(a)' ) '  In this example, we have listed several matrix'
  write ( *, '(a)' ) '  locations TWICE.  DSP_PRINT should compute the'
  write ( *, '(a)' ) '  sum of these values.' 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In particular, we want A(1,2) = 112 and A(3,4) = 234.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros =  ', nz_num

  call dsp_print ( m, n, nz_num, row, col, a, '  The DSP matrix:' )

  return
end
subroutine test5724

!*****************************************************************************80
!
!! TEST5724 tests DSP_WRITE.
!
!  Modified:
!
!    03 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 100
  integer ( kind = 4 ), parameter :: n = 100
  integer ( kind = 4 ), parameter :: nz_num = 3 * n - 2

  real    ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 80 ) output_file
  integer ( kind = 4 ) row(nz_num)

  output_file = 'dsp_matrix.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5724'
  write ( *, '(a)' ) '  For a DSP system,'
  write ( *, '(a)' ) '  DSP_WRITE writes the matrix to a file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =          ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =       ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num
!
!  Set the matrix values.
!
  k = 0
  do i = 1, n
    k = k + 1
    j = i
    row(k) = i
    col(k) = j
    a(k) = real ( 100 * i + j, kind = 8 )
  end do

  do i = 2, n
    j = i - 1
    k = k + 1
    row(k) = i
    col(k) = j
    a(k) = real ( 100 * i + j, kind = 8 )
  end do

  do i = 1, n-1
    j = i + 1
    k = k + 1
    row(k) = i
    col(k) = j
    a(k) = real ( 100 * i + j, kind = 8 )
  end do

  call dsp_print_some ( m, n, nz_num, row, col, a, i4_1, i4_1, &
    i4_10, i4_10, '  Initial 10x10 block of DS3 matrix:' )

  call dsp_write ( m, n, nz_num, row, col, a, output_file )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DSP_WRITE wrote the matrix data to "' &
    // trim ( output_file ) // '".'

  return
end
subroutine test5725

!*****************************************************************************80
!
!! TEST5725 tests DSP_READ, DSP_READ_SIZE.
!
!  Modified:
!
!    03 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ), allocatable, dimension ( : ) :: a
  integer ( kind = 4 ), allocatable, dimension ( : ) :: col
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_10 = 10
  character ( len = 80 ) input_file
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: row

  input_file = 'dsp_matrix.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5725'
  write ( *, '(a)' ) '  For a DSP system,'
  write ( *, '(a)' ) '  DSP_READ_SIZE reads the size of the matrix.'
  write ( *, '(a)' ) '  DSP_READ reads the matrix.'

  call dsp_read_size ( input_file, m, n, nz_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DSP_READ_SIZE reports matrix size data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =          ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =       ', n
  write ( *, '(a,i8)' ) '  Matrix nonzeros NZ_NUM = ', nz_num

  allocate ( row(1:nz_num) )
  allocate ( col(1:nz_num) )
  allocate ( a(1:nz_num) )

  call dsp_read ( input_file, m, n, nz_num, row, col, a )

  call dsp_print_some ( m, n, nz_num, row, col, a, i4_1, i4_1, &
    i4_10, i4_10, '  Initial 10x10 block of recovered DSP matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Deleting the matrix data file "' &
    // trim ( input_file ) // '".'

  call file_delete ( input_file )

  deallocate ( row )
  deallocate ( col )
  deallocate ( a )

  return
end
subroutine test573

!*****************************************************************************80
!
!! TEST573 tests DSR_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz = 7

  integer ( kind = 4 ), dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) off(nz)
  integer ( kind = 4 ), dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST573'
  write ( *, '(a)' ) '  DSR_INDICATOR sets up an DSR indicator matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dsr_indicator ( n, nz, row, col, diag, off )

  call dsr_print ( n, nz, row, col, diag, off, '  The DSR indicator matrix:' )

  return
end
subroutine test574

!*****************************************************************************80
!
!! TEST574 tests DSR_MXV, DSR_VXM.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz = 7

  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) c(n,n)
  integer ( kind = 4 ), dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) off(nz)
  integer ( kind = 4 ), dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST574'
  write ( *, '(a)' ) '  DSR_MXV multiplies an DSR matrix by a vector;'
  write ( *, '(a)' ) '  DSR_VXM multiplies a vector by an DSR matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dsr_random ( n, nz, row, col, diag, off, seed )
!
!  Make an DGE copy.
!
  call dsr_to_dge ( n, nz, row, col, diag, off, c )
!
!  Print the DGE copy.
!
  call dge_print ( n, n, c, '  The DSR matrix, in DGE form:' )

  x(1) = 1.0D+00
  x(2:n-1) = 0.0D+00
  x(n) = -1.0D+00

  call r8vec_print ( n, x, '  The vector x:' )

  call dsr_mxv ( n, nz, row, col, diag, off, x, b )

  call r8vec_print ( n, b, '  The product A * x:' )

  call dsr_vxm ( n, nz, row, col, diag, off, x, b )

  call r8vec_print ( n, b, '  The product A'' * x:' )

  return
end
subroutine test5745

!*****************************************************************************80
!
!! TEST5745 tests DSR_PRINT.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz = 7

  integer ( kind = 4 ), dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) off(nz)
  integer ( kind = 4 ), dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5745'
  write ( *, '(a)' ) '  DSR_PRINT prints an DSR matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dsr_random ( n, nz, row, col, diag, off, seed )
!
!  Print the matrix.
!
  call dsr_print ( n, nz, row, col, diag, off, '  The DSR matrix:' )

  return
end
subroutine test575

!*****************************************************************************80
!
!! TEST575 tests DSR_RANDOM.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nz = 7

  real    ( kind = 8 ) b(n,n)
  integer ( kind = 4 ), dimension ( nz ) :: col = (/ 2, 5, 5, 1, 1, 2, 3 /)
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) off(nz)
  integer ( kind = 4 ), dimension (n+1) :: row = (/ 1, 3, 4, 5, 6, 8 /)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST575'
  write ( *, '(a)' ) '  DSR_RANDOM randomizes an DSR matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dsr_random ( n, nz, row, col, diag, off, seed )
!
!  Make an DGE copy.
!
  call dsr_to_dge ( n, nz, row, col, diag, off, b )
!
!  Print the DGE copy.
!
  call dge_print ( n, n, b, '  The DSR matrix, in DGE form:' )

  return
end
subroutine test577

!*****************************************************************************80
!
!! TEST577 tests DSS_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  real    ( kind = 8 ) a((n*(n+1))/2)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) na

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST577'
  write ( *, '(a)' ) '  For a symmetric skyline storage matrix,'
  write ( *, '(a)' ) '  DSS_INDICATOR computes an indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dss_indicator ( n, na, diag, a )

  call dss_print ( n, na, diag, a, '  The DSS indicator matrix:' )

  return
end
subroutine test58

!*****************************************************************************80
!
!! TEST58 tests DSS_MXV, DSS_PRINT.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  real    ( kind = 8 ) a((n*(n+1))/2)
  real    ( kind = 8 ) a2(n,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) b2(n)
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) na
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST58'
  write ( *, '(a)' ) '  For a symmetric skyline storage matrix,'
  write ( *, '(a)' ) '  DSS_MXV computes A*x,'
  write ( *, '(a)' ) '  DSS_PRINT prints it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dss_random ( n, na, diag, a, seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nonzero entries stored is ', na
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Diagonal storage indices:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i8)' ) i, diag(i)
  end do
!
!  Replace the random entries by marker values.
!
  ij = 0
  do j = 1, n

    if ( j == 1 ) then
      ilo = 1
    else
      ilo = diag(j-1) - diag(j) + j + 1
    end if

    do i = ilo, j
      ij = ij + 1
      a(ij) = real ( 10 * i + j, kind = 8 )
    end do

  end do

  call dss_print ( n, na, diag, a, '  The DSS matrix:' )
!
!  Copy the matrix into a general matrix.
!
  call dss_to_dge ( n, na, diag, a, a2 )
!
!  Set the vector X.
!
  call r8vec_indicator ( n, x )
!
!  Compute the product.
!
  call dss_mxv ( n, na, diag, a, x, b )
!
!  Compute the product using the general matrix.
!
  call dge_mxv ( n, n, a2, x, b2 )
!
!  Compare the results.
!
  call r8vec2_print_some ( n, b, b2, i4_10, '  DSS_MXV verse DGE_MXV' )

  return
end
subroutine test581

!*****************************************************************************80
!
!! TEST581 tests DSTO_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real    ( kind = 8 ), dimension ( n ) :: a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST581'
  write ( *, '(a)' ) '  DSTO_INDICATOR sets up an DSTO indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dsto_indicator ( n, a )

  call dsto_print ( n, a, '  The DSTO indicator matrix:' )

  return
end
subroutine test583

!*****************************************************************************80
!
!! TEST583 tests DSTO_INVERSE.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real    ( kind = 8 ), dimension ( n ) :: a = (/ 4.0D+00, 2.0D+00, 0.8D+00 /)
  real    ( kind = 8 ) a2(n,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) c(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST583'
  write ( *, '(a)' ) '  DSTO_INVERSE computes the inverse of a positive'
  write ( *, '(a)' ) '  definite symmetric Toeplitz matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dsto_print ( n, a, '  The symmetric Toeplitz matrix A:' )

  call dsto_inverse ( n, a, b )

  call dge_print ( n, n, b, '  The inverse matrix B:' )

  call dsto_to_dge ( n, a, a2 )

  c(1:n,1:n) = matmul ( a2(1:n,1:n), b(1:n,1:n) )

  call dge_print ( n, n, c, '  The product C = A * B:' )

  return
end
subroutine test585

!*****************************************************************************80
!
!! TEST585 tests DSTO_MXV, DSTO_YW_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ), dimension ( 0:n ) :: r = (/ &
    1.0D+00, 0.5D+00, 0.2D+00, 0.1D+00 /)
  real    ( kind = 8 ) x(n)

  a(1:n) = r(0:n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST585'
  write ( *, '(a)' ) '  DSTO_YW_SL solves the Yule-Walker equations for a'
  write ( *, '(a)' ) '  symmetric Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dsto_print ( n, a, '  The symmetric Toeplitz matrix:' )

  b(1:n) = -r(1:n)
  call r8vec_print ( n, b, '  The right hand side, B:' )

  b(1:n) = -b(1:n)
  call dsto_yw_sl ( n, b, x )

  call r8vec_print ( n, x, '  The computed solution, X:' )

  call dsto_mxv ( n, a, x, b )

  call r8vec_print ( n, b, '  The product A*X:' )

  return
end
subroutine test587

!*****************************************************************************80
!
!! TEST587 tests DSTO_SL.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real    ( kind = 8 ), dimension ( 0:n-1 ) :: a = (/ 1.0D+00, 0.5D+00, 0.2D+00 /)
  real    ( kind = 8 ), dimension ( n ) :: b = (/ 4.0D+00, -1.0D+00, 3.0D+00 /)
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST587'
  write ( *, '(a)' ) '  DSTO_SL solves a positive definite symmetric '
  write ( *, '(a)' ) '    Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dsto_print ( n, a, '  The symmetric Toeplitz matrix A:' )

  call r8vec_print ( n, b, '  The right hand side vector B:' )

  call dsto_sl ( n, a, b, x )

  call r8vec_print ( n, x, '  The solution X:' )

  call dsto_mxv ( n, a, x, b )

  call r8vec_print ( n, b, '  The product vector  B = A * X:' )

  return
end
subroutine test589

!*****************************************************************************80
!
!! TEST589 tests DTO_INDICATOR.
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(2*n-1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST589'
  write ( *, '(a)' ) '  DTO_INDICATOR sets up an DTO indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call dto_indicator ( n, a )

  call dto_print ( n, a, '  The DTO indicator matrix:' )

  return
end
subroutine test59

!*****************************************************************************80
!
!! TEST59 tests DTO_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(2*n-1)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST59'
  write ( *, '(a)' ) '  DTO_SL solves a Toeplitz system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dto_random ( n, seed, a )

  call dto_print ( n, a, '  The Toeplitz matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call dto_mxv ( n, a, x, b )
    else
      call dto_vxm ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call dto_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call r8vec_print_some ( n, x, i4_10, '  Solution:' )
    else
      call r8vec_print_some ( n, x, i4_10, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine test60

!*****************************************************************************80
!
!! TEST60 tests DUT_DET, DUT_INVERSE, DUT_MXM, DUT_RANDOM.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) c(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST60'
  write ( *, '(a)' ) '  For an upper triangular matrix,'
  write ( *, '(a)' ) '  DUT_DET computes the determinant.'
  write ( *, '(a)' ) '  DUT_INVERSE computes the inverse.'
  write ( *, '(a)' ) '  DUT_MXM computes matrix products.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call dut_random ( n, n, seed, a )

  call dut_print ( n, n, a, '  The matrix A:' )
!
!  Compute the determinant.
!
  call dut_det ( n, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant is ', det
!
!  Compute the inverse matrix B.
!
  b(1:n,1:n) = a(1:n,1:n)

  call dut_inverse ( n, b )

  call dut_print ( n, n, b, '  The inverse matrix B:' )
!
!  Check
!
  call dut_mxm ( n, a, b, c )

  call dut_print ( n, n, c, '  The product A * B:' )

  return
end
subroutine test605

!*****************************************************************************80
!
!! TEST605 tests DUT_INDICATOR.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 8
  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST605'
  write ( *, '(a)' ) '  For an upper triangular matrix,'
  write ( *, '(a)' ) '  DUT_INDICATOR sets up an indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  call dut_indicator ( m, n, a )

  call dut_print ( m, n, a, '  The DUT indicator matrix:' )

  return
end
subroutine test61

!*****************************************************************************80
!
!! TEST61 tests DUT_MXV, DUT_SL, DUT_VXM.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST61'
  write ( *, '(a)' ) '  For an upper triangular matrix,'
  write ( *, '(a)' ) '  DUT_SL solves systems;'
  write ( *, '(a)' ) '  DUT_MXV computes matrix-vector products.'
  write ( *, '(a)' ) '  DUT_VXM computes vector-matrix products.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do i = 1, n
    do j = 1, n
      if ( i <= j ) then
        a(i,j) = real ( j, kind = 8 )
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do

  call dut_print ( n, n, a, '  The upper triangular matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call dut_mxv ( n, n, a, x, b )
    else
      call dut_vxm ( n, n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call dut_sl ( n, a, b, job )

    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution:' )
    else
      call r8vec_print ( n, b, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine test62

!*****************************************************************************80
!
!! TEST62 tests R8VM_DET.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) a2(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST62'
  write ( *, '(a)' ) '  R8VM_DET, determinant of a Vandermonde matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8vm_random ( n, n, seed, a )

  call r8vm_print ( n, n, a, '  The Vandermonde matrix:' )
!
!  Copy the matrix into a general array.
!
  call r8vm_to_dge ( n, n, a, a2 )
!
!  Compute the determinant.
!
  call r8vm_det ( n, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8VM_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call dge_fa ( n, a2, pivot, info )
!
!  Compute the determinant.
!
  call dge_det ( n, a2, pivot, det )

  write ( *, '(a,g14.6)' ) '  DGE_DET computes the determinant = ', det

  return
end
subroutine test63

!*****************************************************************************80
!
!! TEST63 tests R8VM_SL.
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ), parameter :: i4_10 = 10
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST63'
  write ( *, '(a)' ) '  R8VM_SL solves a Vandermonde system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8vm_random ( n, n, seed, a )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8vm_mxv ( n, n, a, x, b )
    else
      call r8vm_vxm ( n, n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call r8vm_sl ( n, a, b, x, job, info )

    if ( job == 0 ) then
      call r8vec_print_some ( n, x, i4_10, '  Solution:' )
    else
      call r8vec_print_some ( n, x, i4_10, '  Solution to transposed system:' )
    end if

  end do
 
  return
end

