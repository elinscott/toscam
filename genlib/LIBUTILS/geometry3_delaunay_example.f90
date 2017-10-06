program main

!*****************************************************************************80
!
!! MAIN is the main program for GEOMPACK_PRB.
!
!  Discussion:
!
!    GEOMPACK_PRB calls a set of problems for GEOMPACK.
!
!  Modified:
!
!    02 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Tests for the routines in the GEOMPACK library.'
 
  call test01
  call test02
  call test03
  call test05
  call test06
  call test07
  call test08

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
   
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 tests POINTS_DELAUNAY_NAIVE_2D.
!
!  Diagram:
!
!   !....3&11....
!   !............
!   !............
!   X..9.........
!   !.....5......
!   !...........6
!   !.4.2...10...
!   !.....8...12.
!   V............
!   !..7.........
!   !......1.....
!   !............
!   !............
!   !----V----X--
! 
!
  implicit none

  integer, parameter :: node_num = 12
  integer, parameter :: dim_num = 2

  integer, parameter :: maxtri = 2 * node_num - 3

  integer, dimension ( 3, maxtri ) :: triangle_node
  integer triangle_num
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: node_xy = reshape ( (/ &
     7.0D+00,  3.0D+00, &
     4.0D+00,  7.0D+00, &
     5.0D+00, 13.0D+00, &
     2.0D+00,  7.0D+00, &
     6.0D+00,  9.0D+00, &
    12.0D+00, 10.0D+00, &
     3.0D+00,  4.0D+00, &
     6.0D+00,  6.0D+00, &
     3.0D+00, 10.0D+00, &
     8.0D+00,  7.0D+00, &
     5.0D+00, 13.0D+00, &
    10.0D+00,  6.0D+00 /), (/ dim_num, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  POINTS_DELAUNAY_NAIVE_2D computes the Delaunay'
  write ( *, '(a)' ) '  triangulation of a set of points.'

  call r8mat_transpose_print ( dim_num, node_num, node_xy, '  The points:' )

  call points_delaunay_naive_2d ( node_num, node_xy, maxtri, triangle_num, &
    triangle_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of triangles is TRIANGLE_NUM = ', &
    triangle_num

  call i4mat_transpose_print ( 3, triangle_num, triangle_node, &
    '  The triangles:' )

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 tests D2VEC_PART_QUICK_A.
!
  implicit none

  integer, parameter :: node_num = 12
  integer, parameter :: dim_num = 2

  real ( kind = 8 ), dimension ( dim_num, node_num ) :: a
  real ( kind = 8 ) :: b = 0.0D+00
  real ( kind = 8 ) :: c = 10.0D+00
  integer l
  integer r
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  D2VEC_PART_QUICK_A reorders a D2 vector'
  write ( *, '(a)' ) '    as part of a quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r8mat_uniform ( dim_num, node_num, b, c, seed, a )

  call r8mat_transpose_print ( dim_num, node_num, a, '  Before rearrangment:' )

  call r82vec_part_quick_a ( node_num, a, l, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rearranged array'
  write ( *, '(a,i6)' ) '  Left index =  ', l
  write ( *, '(a,i6)' ) '  Key index =   ', l+1
  write ( *, '(a,i6)' ) '  Right index = ', r

  call r8mat_transpose_print ( dim_num, l,     a(1:dim_num,1:l),   &
    '  Left half:' )

  call r8mat_transpose_print ( dim_num, 1,     a(1:dim_num,l+1),   '  Key:' )

  call r8mat_transpose_print ( dim_num, node_num-l-1, &
    a(1:dim_num,l+2:node_num), &
    '  Right half:' )

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 tests D2VEC_SORT_QUICK_A.
!
  implicit none

  integer, parameter :: node_num = 12
  integer, parameter :: dim_num = 2

  real ( kind = 8 ) a(dim_num,node_num)
  real ( kind = 8 ) :: b = 0.0D+00
  real ( kind = 8 ) :: c = 10.0D+00
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  R82VEC_SORT_QUICK_A sorts a D2 vector'
  write ( *, '(a)' ) '    using quick sort.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed

  call r8mat_uniform ( dim_num, node_num, b, c, seed, a )
!
!  Give a few elements the same first component.
!
  a(1,3) = a(1,5)
  a(1,4) = a(1,12)
!
!  Give a few elements the same second component.
!
  a(2,6) = a(2,1)
  a(2,2) = a(2,9)
!
!  Make two entries equal.
!
  a(1:2,7) = a(1:2,11)

  call r8mat_transpose_print ( dim_num, node_num, a, '  Before rearrangement:' )

  call r82vec_sort_quick_a ( node_num, a )

  call r8mat_transpose_print ( dim_num, node_num, a, '  Sorted array:' )

  return
end
subroutine test05

!*****************************************************************************80
!
!! TEST05 tests DTRIS2.
!
  implicit none

  integer, parameter :: node_num = 9
  integer, parameter :: dim_num = 2

  integer triangle_neighbor(3,2*node_num)
  integer triangle_node(3,2*node_num)
  integer triangle_num
  real ( kind = 8 ), dimension (dim_num,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       0.2D+00, 0.5D+00, &
       0.3D+00, 0.6D+00, &
       0.4D+00, 0.5D+00, &
       0.6D+00, 0.4D+00, &
       0.6D+00, 0.5D+00, &
       1.0D+00, 0.0D+00, &
       1.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  DTRIS2 computes the Delaunay triangulation of'
  write ( *, '(a)' ) '    a pointset in 2D.'
!
!  Set up the Delaunay triangulation.
!
  call dtris2 ( node_num, node_xy, triangle_num, triangle_node, &
    triangle_neighbor )

  call triangulation_order3_print ( node_num, triangle_num, node_xy, &
    triangle_node, triangle_neighbor )

  return
end
subroutine test06

!*****************************************************************************80
!
!! TEST06 tests TRIANGLE_CIRCUMCENTER_2D;
!
  implicit none

  integer, parameter :: dim_num = 2
  integer, parameter :: test_num = 4

  real ( kind = 8 ) center(dim_num)
  integer i
  real ( kind = 8 ) t(dim_num,3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For a triangle in 2D:'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCENTER_2D computes the circumcenter.'

  do i = 1, test_num

    if ( i == 1 ) then
      t(1:dim_num,1:3) = reshape ( (/ &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.0D+00,  1.0D+00 /), (/ dim_num, 3 /) )
    else if ( i == 2 ) then
      t(1:dim_num,1:3) = reshape ( (/ &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00,  0.86602539D+00 /), (/ dim_num, 3 /) )
    else if ( i == 3 ) then
      t(1:dim_num,1:3) = reshape ( (/ &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00, 10.0D+00 /), (/ dim_num, 3 /) )
    else if ( i == 4 ) then
      t(1:dim_num,1:3) = reshape ( (/ &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
        10.0D+00,  2.0D+00 /), (/ dim_num, 3 /) )
    end if

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices' )

    call triangle_circumcenter_2d ( t, center )

    call r8vec_print ( dim_num, center, '  Circumcenter :' )

 end do

  return
end
subroutine test07

!*****************************************************************************80
!
!! TEST07 tests TRIANGULATION_ORDER3_PLOT.
!
  implicit none

  integer, parameter :: dim_num = 2
  integer, parameter :: node_num = 9
  integer, parameter :: triangle_num = 12

  character ( len = 80 ) :: file_name = 'triangulation_plot.eps'
  integer :: node_show = 2
  real ( kind = 8 ), dimension (dim_num,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       0.2D+00, 0.5D+00, &
       0.3D+00, 0.6D+00, &
       0.4D+00, 0.5D+00, &
       0.6D+00, 0.4D+00, &
       0.6D+00, 0.5D+00, &
       1.0D+00, 0.0D+00, &
       1.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  integer, dimension ( 3, triangle_num ) :: triangle_node = reshape ( (/ &
       2, 1, 3, &
       3, 1, 6, &
       2, 3, 4, &
       4, 3, 5, &
       7, 4, 5, &
       5, 3, 6, &
       7, 5, 6, &
       9, 4, 7, &
       6, 1, 8, &
       7, 6, 8, &
       7, 8, 9, &
       2, 4, 9 /), (/ 3, triangle_num /) )
  integer :: triangle_show = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PLOT can plot a triangulation.'

  call triangulation_order3_plot ( file_name, node_num, node_xy, &
    triangle_num, triangle_node, node_show, triangle_show )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PLOT has created an'
  write ( *, '(a)' ) '  Encapsulated PostScript file (EPS) containing'
  write ( *, '(a)' ) '  an aimage of the triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This file is called "' // trim ( file_name ) //'".'

  return
end
subroutine test08

!*****************************************************************************80
!
!! TEST08 tests TRIANGULATION_ORDER3_PRINT.
!
  implicit none

  integer, parameter :: node_num = 9
  integer, parameter :: triangle_num = 12

  real ( kind = 8 ), dimension (2,node_num) :: node_xy = reshape ( (/ &
       0.0D+00, 0.0D+00, &
       0.0D+00, 1.0D+00, &
       0.2D+00, 0.5D+00, &
       0.3D+00, 0.6D+00, &
       0.4D+00, 0.5D+00, &
       0.6D+00, 0.4D+00, &
       0.6D+00, 0.5D+00, &
       1.0D+00, 0.0D+00, &
       1.0D+00, 1.0D+00 /), (/ 2, node_num /) )
  integer, dimension (3,triangle_num) :: triangle_node = reshape ( (/ &
       2, 1, 3, &
       3, 1, 6, &
       2, 3, 4, &
       4, 3, 5, &
       7, 4, 5, &
       5, 3, 6, &
       7, 5, 6, &
       9, 4, 7, &
       6, 1, 8, &
       7, 6, 8, &
       7, 8, 9, &
       2, 4, 9 /), (/ 3, triangle_num /) )
  integer, dimension (3,triangle_num) :: triangle_neighbor = reshape ( (/ &
       -28,   2,  3, &
         1,   9,  6, &
         1,   4, 12, &
         3,   6,  5, &
         8,   4,  7, &
         4,   2,  7, &
         5,   6, 10, &
        12,   5, 11, &
         2, -34, 10, &
         7,   9, 11, &
        10, -38,  8, &
         3,   8, -3 /), (/ 3, triangle_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  TRIANGULATION_ORDER3_PRINT prints out a triangulation.'

  call triangulation_order3_print ( node_num, triangle_num, node_xy, &
    triangle_node, triangle_neighbor )

  return
end
  
