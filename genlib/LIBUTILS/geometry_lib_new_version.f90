subroutine angle_contains_ray_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )
!
!*******************************************************************************
!
!! ANGLE_CONTAINS_RAY_2D determines if an angle contains a ray, in 2D.
!
!
!  Discussion:
!
!    The angle is defined by the sequence of points (X1,Y1), (X2,Y2)
!    and (X3,Y3).
!
!    The ray is defined by the sequence of points (X2,Y2), (X,Y).
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates of
!    the angle.
!
!    Input, real X, Y, the end point of the ray to be checked.
!    The ray is assumed to have an origin at (X2,Y2).
!
!    Output, logical INSIDE, is .TRUE. if the ray is inside
!    the angle or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  real a1
  real a2
  real angle_deg_2d
  logical inside
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  a1 = angle_deg_2d ( x1, y1, x2, y2, x, y )
  a2 = angle_deg_2d ( x1, y1, x2, y2, x3, y3 )

  if ( a1 <= a2 ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
function angle_deg_2d ( x1, y1, x2, y2, x3, y3 )
!
!*******************************************************************************
!
!! ANGLE_DEG_2D returns the angle swept out between two rays in 2D.
!
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_DEG_2D(X1,Y1,X2,Y2,X3,Y3)
!    + ANGLE_DEG_2D(X3,Y3,X2,Y2,X1,Y1) = 360.0
!
!  Modified:
!
!    14 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, define the rays
!    ( X1-X2, Y1-Y2 ) and ( X3-X2, Y3-Y2 ) which in turn define the
!    angle, counterclockwise from ( X1-X2, Y1-Y2 ).
!
!    Output, real ANGLE_DEG_2D, the angle swept out by the rays, measured
!    in degrees.  0 <= ANGLE_DEG_2D < 360.  If either ray has zero length,
!    then ANGLE_DEG_2D is set to 0.
!
  implicit none
!
  real angle_deg_2d
  real angle_rad_2d
  real r_pi
  real radians_to_degrees
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
  y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )

  if ( x == 0.0E+00 .and. y == 0.0E+00 ) then

    angle_deg_2d = 0.0E+00

  else

    angle_rad_2d = atan2 ( y, x )

    if ( angle_rad_2d < 0.0E+00 ) then
      angle_rad_2d = angle_rad_2d + 2.0E+00 * r_pi ( )
    end if

    angle_deg_2d = radians_to_degrees ( angle_rad_2d )

  end if

  return
end
function angle_rad_2d ( x1, y1, x2, y2, x3, y3 )
!
!*******************************************************************************
!
!! ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
!
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D(X1,Y1,X2,Y2,X3,Y3)
!    + ANGLE_RAD_2D(X3,Y3,X2,Y2,X1,Y1) = 2 * PI
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, define the rays
!    ( X1-X2, Y1-Y2 ) and ( X3-X2, Y3-Y2 ) which in turn define the
!    angle, counterclockwise from ( X1-X2, Y1-Y2 ).
!
!    Output, real ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
!
  implicit none
!
  real angle_rad_2d
  real r_pi
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
  y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )

  if ( x == 0.0E+00 .and. y == 0.0E+00 ) then

    angle_rad_2d = 0.0E+00

  else

    angle_rad_2d = atan2 ( y, x )

    if ( angle_rad_2d < 0.0E+00 ) then
      angle_rad_2d = angle_rad_2d + 2.0E+00 * r_pi ( )
    end if

  end if

  return
end
function angle_rad_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! ANGLE_RAD_3D returns the angle in radians between two rays in 3D.
!
!
!  Discussion:
!
!    The routine always computes the SMALLER of the two angles between
!    two rays.  Thus, if the rays make an (exterior) angle of
!    1.5 radians, the (interior) angle of 0.5 radians will be reported.
!
!  Formula:
!
!    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are three points
!    which define the rays.  The rays are:
!    ( X1-X2, Y1-Y2, Z1-Z2 ) and ( X3-X2, Y3-Y2, Z3-Z2 ).
!
!    Output, real ANGLE_RAD_3D, the angle between the two rays, in radians.
!    This value will always be between 0 and PI.  If either ray has
!    zero length, then the angle is returned as zero.
!
  implicit none
!
  real angle_rad_3d
  real arc_cosine
  real dot
  real dot0_3d
  real enorm0_3d
  real theta
  real v1norm
  real v2norm
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  v1norm = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
  v2norm = enorm0_3d ( x3, y3, z3, x2, y2, z2 )

  if ( v1norm == 0.0E+00 .or. v2norm == 0.0E+00 ) then
    angle_rad_3d = 0.0E+00
  else
    dot = dot0_3d ( x2, y2, y2, x1, y1, z1, x3, y3, z3 )
    angle_rad_3d = arc_cosine ( dot / ( v1norm * v2norm ) )
  end if

  return
end
function angle_rad_nd ( n, vec1, vec2 )
!
!*******************************************************************************
!
!! ANGLE_RAD_ND returns the angle in radians between two rays in ND.
!
!
!  Discussion:
!
!    This routine always computes the SMALLER of the two angles between
!    two rays.  Thus, if the rays make an (exterior) angle of 1.5 PI,
!    then the (interior) angle of 0.5 PI is reported.
!
!  Formula:
!
!    X dot Y = Norm(X) * Norm(Y) * Cos( Angle(X,Y) )
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the rays.
!
!    Input, real VEC1(N), VEC2(N), the two rays to be considered.
!
!    Output, real ANGLE_RAD_ND, the angle between the rays, in radians.
!    This value will always be between 0 and PI.
!
  implicit none
!
  integer n
!
  real angle_rad_nd
  real arc_cosine
  real dot
  real dot_nd
  real enorm_nd
  real theta
  real v1norm
  real v2norm
  real vec1(n)
  real vec2(n)
!
  dot = dot_nd ( n, vec1, vec2 )

  v1norm = enorm_nd ( n, vec1 )

  v2norm = enorm_nd ( n, vec2 )

  if ( v1norm == 0.0E+00 .or. v2norm == 0.0E+00 ) then
    angle_rad_nd = 0.0E+00
  else
    angle_rad_nd = arc_cosine ( dot / ( v1norm * v2norm ) )
  end if

  return
end
function anglei_deg_2d ( x1, y1, x2, y2, x3, y3 )
!
!*******************************************************************************
!
!! ANGLEI_DEG_2D returns the interior angle in degrees between two rays in 2D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, define the rays
!    (X1-X2,Y1-Y2) and (X3-X2,Y3-Y2) which in turn define the angle.
!
!    Output, real ANGLEI_DEG_2D, the angle swept out by the rays, measured
!    in degrees.  This value satisfies 0 <= ANGLEI_DEG_2D < 180.0.  If either
!    ray is of zero length, then ANGLEI_DEG_2D is returned as 0.
!
  implicit none
!
  real anglei_deg_2d
  real anglei_rad_2d
  real r_pi
  real radians_to_degrees
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
  y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )

  if ( x == 0.0E+00 .and. y == 0.0E+00 ) then

    anglei_deg_2d = 0.0E+00

  else

    anglei_rad_2d = atan2 ( y, x )

    if ( anglei_rad_2d < 0.0E+00 ) then
      anglei_rad_2d = anglei_rad_2d + 2.0E+00 * r_pi ( )
    end if

    anglei_deg_2d = radians_to_degrees ( anglei_rad_2d )

    if ( anglei_deg_2d > 180.0E+00 ) then
      anglei_deg_2d = 360.0E+00 - anglei_deg_2d
    end if

  end if

  return
end
function anglei_rad_2d ( x1, y1, x2, y2, x3, y3 )
!
!*******************************************************************************
!
!! ANGLEI_RAD_2D returns the interior angle in radians between two rays in 2D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, define the rays
!    (X1-X2,Y1-Y2) and (X3-X2,Y3-Y2) which in turn define the angle.
!
!    Output, real ANGLEI_RAD_2D, the angle swept out by the rays, measured
!    in radians.  This value satisfies 0 <= ANGLEI_RAD_2D < PI.  If either
!    ray is of zero length, then ANGLEI_RAD_2D is returned as 0.
!
  implicit none
!
  real anglei_rad_2d
  real r_pi
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
  y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )

  if ( x == 0.0E+00 .and. y == 0.0E+00 ) then

    anglei_rad_2d = 0.0E+00

  else

    anglei_rad_2d = atan2 ( y, x )

    if ( anglei_rad_2d < 0.0E+00 ) then
      anglei_rad_2d = anglei_rad_2d + 2.0E+00 * r_pi ( )
    end if

  end if

  return
end
function arc_cosine ( c )
!
!*******************************************************************************
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
!    This routine truncates arguments outside the range.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real C, the argument.
!
!    Output, real ARC_COSINE, an angle whose cosine is C.
!
  implicit none
!
  real arc_cosine
  real c
  real c2
!
  c2 = c
  c2 = max ( c2, -1.0E+00 )
  c2 = min ( c2, +1.0E+00 )

  arc_cosine = acos ( c2 )

  return
end
function atan4 ( y, x )
!
!*******************************************************************************
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
!
!
!  Discussion:
!
!    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
!    the built in functions ATAN and ATAN2 already do.
!
!    However:
!
!    * ATAN4 always returns a positive angle, between 0 and 2 PI,
!      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
!      and [-PI,+PI] respectively;
!
!    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
!     function by contrast always returns an angle in the first or fourth
!     quadrants.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Y, X, two quantities which represent the tangent of
!    an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ATAN4, an angle between 0 and 2 * PI, whose tangent is
!    (Y/X), and which lies in the appropriate quadrant so that the signs
!    of its cosine and sine match those of X and Y.
!
  implicit none
!
  real abs_x
  real abs_y
  real atan4
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real theta
  real theta_0
  real x
  real y
!
!  Special cases:
!
  if ( x == 0.0E+00 ) then

    if ( y > 0.0E+00 ) then
      theta = pi / 2.0E+00
    else if ( y < 0.0E+00 ) then
      theta = 3.0E+00 * pi / 2.0E+00
    else if ( y == 0.0E+00 ) then
      theta = 0.0E+00
    end if

  else if ( y == 0.0E+00 ) then

    if ( x > 0.0E+00 ) then
      theta = 0.0E+00
    else if ( x < 0.0E+00 ) then
      theta = PI
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( x > 0.0E+00 .and. y > 0.0E+00 ) then
      theta = theta_0
    else if ( x < 0.0E+00 .and. y > 0.0E+00 ) then
      theta = pi - theta_0
    else if ( x < 0.0E+00 .and. y < 0.0E+00 ) then
      theta = pi + theta_0
    else if ( x > 0.0E+00 .and. y < 0.0E+00 ) then
      theta = 2.0E+00 * pi - theta_0
    end if

  end if

  atan4 = theta

  return
end
subroutine ball_unit_sample_2d ( x )
!
!*******************************************************************************
!
!! BALL_UNIT_SAMPLE_2D picks a random point in the unit ball in 2D.
!
!
!  Discussion:
!
!    The unit ball is the set of points (X,Y) such that
!
!      X(1)**2 + X(2)**2 <= 1.
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real X(2), a random point in the unit ball.
!
  implicit none
!
  real r
  real r_pi
  real theta
  real u(2)
  real x(2)
!
  call random_number ( harvest = u(1:2) )

  r = sqrt ( u(1) )
  theta = 2.0E+00 * r_pi ( ) * u(2)

  x(1) = r * cos ( theta )
  x(2) = r * sin ( theta )

  return
end
subroutine ball_unit_sample_3d ( x )
!
!*******************************************************************************
!
!! BALL_UNIT_SAMPLE_3D picks a random point in the unit ball in 3D.
!
!
!  Modified:
!
!    19 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real X(3), the sample point.
!
  implicit none
!
  real arc_cosine
  integer i
  real phi
  real r
  real r_pi
  real theta
  real u(3)
  real vdot
  real x(3)
!
  call random_number ( harvest = u(1:3) )
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  vdot = 2.0E+00 * u(1) - 1.0E+00

  phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  theta = 2.0E+00 * r_pi ( ) * u(2)
!
!  Pick a random radius R.
!
  r = u(3)**( 1.0E+00 / 3.0E+00 )

  x(1) = r * cos ( theta ) * sin ( phi )
  x(2) = r * sin ( theta ) * sin ( phi )
  x(3) = r * cos ( phi )

  return
end
subroutine ball_unit_sample_nd ( n, x )
!
!*******************************************************************************
!
!! BALL_UNIT_SAMPLE_ND picks a random point in the unit ball in ND.
!
!
!  Discussion:
!
!    N-1 random Givens rotations are applied to the point ( 1, 0, 0, ..., 0 ).
!
!    The I-th Givens rotation is in the plane of coordinate axes I and I+1,
!    and has the form:
!
!     [ cos ( theta )  - sin ( theta ) ] * x(i)      = x'(i)
!     [ sin ( theta )    cos ( theta ) ]   x(i+1)      x'(i+1)
!
!    Finally, a scaling is applied to set the point at a distance R
!    from the center, in a way that results in a uniform distribution.
!
!  Modified:
!
!    19 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Output, real X(N), the random point.
!
  implicit none
!
  integer n
!
  integer i
  real r
  real random_cosine(1:n-1)
  real random_sign(1:n-1)
  real random_sine(1:n-1)
  real x(n)
  real xi
!
  x(1) = 1.0E+00
  x(2:n) = 0.0E+00

  call random_number ( harvest = random_cosine(1:n-1) )
  random_cosine(1:n-1) = 2.0E+00 * random_cosine(1:n-1) - 1.0E+00

  call random_number ( harvest = random_sign(1:n-1) )
  random_sign(1:n-1) = real ( 2 * int ( 2.0E+00 * random_sign(1:n-1) ) - 1 )

  random_sine(1:n-1) = &
    random_sign(1:n-1) * sqrt ( 1.0E+00 - random_cosine(1:n-1)**2 )

  do i = 1, n-1
    xi = x(i)
    x(i  ) = random_cosine(i) * xi
    x(i+1) = random_sine(i)   * xi
  end do

  call random_number ( harvest = r )
  r = r**( 1.0E+00 / real ( n ) )

  x(1:n) = r * x(1:n)

  return
end
subroutine basis_map_3d ( u1, u2, u3, v1, v2, v3, a, ierror )
!
!*******************************************************************************
!
!! BASIS_MAP_3D computes the matrix which maps one basis to another.
!
!
!  Discussion:
!
!    As long as the vectors U1, U2 and U3 are linearly independent,
!    a matrix A will be computed that maps U1 to V1, U2 to V2, and
!    U3 to V3.
!
!    Depending on the values of the vectors, A may represent a
!    rotation, reflection, dilation, project, or a combination of these
!    basic linear transformations.
!
!  Modified:
!
!    20 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real U1(3), U2(3), U3(3), the three "domain" or "preimage"
!    vectors, which should be linearly independent.
!
!    Input, real V1(3), V2(3), V3(3), the three "range" or "image" vectors.
!
!    Output, real A(3,3), a matrix with the property that A * U1 = V1,
!    A * U2 = V2 and A * U3 = V3.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    nonzero, the matrix [ U1 | U2 | U3 ] is exactly singular.
!
  implicit none
!
  real a(3,3)
  real b(3,3)
  real c(3,3)
  real det
  integer i
  integer ierror
  integer j
  integer k
  real u1(3)
  real u2(3)
  real u3(3)
  real v1(3)
  real v2(3)
  real v3(3)
!
  ierror = 0
!
!  Set B = [ U1 | U2 | U3 ].
!
  b(1:3,1) = u1(1:3)
  b(1:3,2) = u2(1:3)
  b(1:3,3) = u3(1:3)
!
!  Compute C = the inverse of [ U1 | U2 | U3 ].
!
  call rmat3_inverse ( b, c, det )

  if ( det == 0.0E+00 ) then
    ierror = 1
    return
  end if
!
!  Set B = [ V1 | V2 | V3 ].
!
  b(1:3,1) = v1(1:3)
  b(1:3,2) = v2(1:3)
  b(1:3,3) = v3(1:3)
!
!  A = [ V1 | V2 | V3 ] * inverse [ U1 | U2 | U3 ].
!
  a(1:3,1:3) = matmul ( b(1:3,1:3), c(1:3,1:3) )

  return
end
subroutine bin_search_one_2d ( bin, nset, pset, nbin, bin_start, bin_next, &
  ptest, found_a_neighbor, i_min, d_min_sq, compares )
!
!*******************************************************************************
!
!! BIN_SEARCH_ONE_2D searches one cell in a 2D array of bins.
!
!
!  Discussion:
!
!    The bins are presumed to have been set up by successive calls to:
!
!      R2VEC_BIN_EVEN2,
!      R2VEC_BINNED_REORDER, and
!      R2VEC_BINNED_SORT_A.
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer BIN(2), the indices of the cell to be examined.
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN(2), the number of cells in the horizontal and
!    vertical directions.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    indicates the index of the first and last element in the bin, or -1 
!    if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, real PTEST(2), the coordinates of the test point.
!
!    Input/output, logical FOUND_A_NEIGHBOR, is set to TRUE if at least
!    one point of PSET is found in the current bin.  Otherwise, it retains its
!    input value.
!
!    Input/output, integer I_MIN, the index of the nearest neighbor in
!    PSET to PTEST, if at least one neighbor has been found.
!
!    Input/output, real D_MIN_SQ, the square of the distance from the nearest
!    neighbor in PSET to PTEST, if at least one neighbor has been found.
!
!    Input/output, integer COMPARES, the number of elements of PSET whose
!    distance to PTEST has been computed.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  integer nbin(ndim)
  integer nset
!
  integer bin(ndim)
  integer bin_next(nset)
  integer bin_start(nbin(1),nbin(2))
  integer compares
  real d_min_sq
  real d_sq
  logical found_a_neighbor
  integer i_min
  integer node
  real pset(ndim,nset)
  real ptest(ndim)
!
  node = bin_start(bin(1),bin(2))

  do while ( node > 0 )

    found_a_neighbor = .true.

    d_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )
    compares = compares + 1

    if ( d_sq < d_min_sq ) then
      d_min_sq = d_sq
      i_min = node
    end if

    node = bin_next(node)

  end do

  return
end
subroutine bin_to_r_even2 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R_EVEN2 returns the limits for a given "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A = 10, B = 20
!
!    BIN      CMIN  CMAX
!
!    1         10    12
!    2         12    14
!    3         14    16
!    4         16    18
!    5         18    20
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, integer BIN, the index of the bin to be considered.
!    If BIN is less than 1, or greater than NBIN, the user will get what
!    the user deserves.
!
!    Input, real A, B, the lower and upper limits of the bin interval.
!    While A is expected to be less than B, the code should return useful
!    results if A is actually greater than B.
!
!    Output, real CMIN, CMAX, the minimum and maximum limits on the bin.
!
  implicit none
!
  real a
  real b
  integer bin
  real cmax
  real cmin
  integer nbin
!
!  Compute the bin limits.
!
  if ( bin < 1 ) then
    cmin = - huge ( cmin )
    cmax = a
  else if ( bin <= nbin ) then
    cmin = ( real ( nbin - bin + 1 ) * a + real ( bin - 1 ) * b ) &
      / real ( nbin )
    cmax = ( real ( nbin - bin ) * a + real ( bin ) * b ) &
      / real ( nbin )
  else if ( bin > nbin ) then
    cmin = b
    cmax = huge ( cmax )
  end if

  return
end
subroutine bin_to_r2_even2 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R2_EVEN2 returns the limits for a given R2 "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A(1) = 5, B(1) = 15
!              A(2) = 0, B(2) = 20
!
!     BIN         CMIN      CMAX
!    ------   -----------  --------
!    1, 1       5     0     7     4
!    2, 2       7     4     9     8
!    3, 3       9     8    11    12
!    4, 4      11    12    13    16
!    5, 5      13    16    15    20
!
!  Modified:
!
!    07 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, integer BIN(2), the index of the bin to be considered.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Output, real CMIN(2), CMAX(2), the minimum and maximum limits on the bin.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real cmax(ndim)
  real cmin(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call bin_to_r_even2 ( nbin, bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r2_even3 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R2_EVEN3 returns the limits for a given R2 "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A(I) to B(I) is divided into NBIN(I) equal
!    subintervals or bins.
!
!  Example:
!
!    NBIN = (/ 4, 5, /)
!
!    A(1) = 5, B(1) = 15
!    A(2) = 0, B(2) = 20
!
!     BIN         CMIN      CMAX
!    ------   -----------  --------
!    1, 1       5     0     7     4
!    2, 2       7     4     9     8
!    3, 3       9     8    11    12
!    4, 4      11    12    13    16
!    5, 5      13    16    15    20
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!
!    Input, integer BIN(2), the index of the bin to be considered.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Output, real CMIN(2), CMAX(2), the minimum and maximum limits on the bin.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real cmax(ndim)
  real cmin(ndim)
  integer i
  integer nbin(ndim)
!
  do i = 1, ndim
    call bin_to_r_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r3_even2 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R3_EVEN2 returns the limits for a given R3 "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, integer BIN(3), the index of the bin to be considered.
!
!    Input, real A(3), B(3), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Output, real CMIN(3), CMAX(3), the minimum and maximum limits on the bin.
!
  implicit none
!
  integer, parameter :: ndim = 3
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real cmax(ndim)
  real cmin(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call bin_to_r_even2 ( nbin, bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r3_even3 ( nbin, bin, a, b, cmin, cmax )
!
!*******************************************************************************
!
!! BIN_TO_R3_EVEN3 returns the limits for a given R3 "bin" in [A,B].
!
!
!  Discussion:
!
!    The interval from A(I) to B(I) is divided into NBIN(I) equal
!    subintervals or bins.
!
!  Modified:
!
!    30 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN(3), the number of bins in each dimension.
!
!    Input, integer BIN(3), the index of the bin to be considered.
!
!    Input, real A(3), B(3), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Output, real CMIN(3), CMAX(3), the minimum and maximum limits on the bin.
!
  implicit none
!
  integer, parameter :: ndim = 3
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real cmax(ndim)
  real cmin(ndim)
  integer i
  integer nbin(ndim)
!
  do i = 1, ndim
    call bin_to_r_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine box_clip_line_2d ( xmin, ymin, xmax, ymax, x1, y1, x2, y2, x3, y3, &
  x4, y4, ival )
!
!*******************************************************************************
!
!! BOX_CLIP_LINE_2D uses a box to clip a line segment in 2D.
!
!
!  Discussion:
!
!    The box is assumed to be a rectangle with sides aligned on coordinate
!    axes.
!
!  Modified:
!
!    18 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XMIN, YMIN, XMAX, YMAX, the minimum and maximum X and Y
!    values, which define the box.
!
!    Input, real X1, Y1, X2, Y2, the coordinates of the endpoints of the
!    line segment.
!
!    Output, real X3, Y3, X4, Y4, the clipped coordinates.
!
!    Output, integer IVAL:
!    -1, no part of the line segment is within the box.
!     0, no clipping was necessary.  The line segment is entirely within the box.
!     1, (X1,Y1) was clipped.
!     2, (X2,Y2) was clipped.
!     3, (X1,Y1) and (X2,Y2) were clipped.
!
  implicit none
!
  integer ival
  logical l1
  logical l2
  real x
  real x1
  real x2
  real x3
  real x4
  real xmax
  real xmin
  real y
  real y1
  real y2
  real y3
  real y4
  real ymax
  real ymin
!
  l1 = .false.
  l2 = .false.

  x3 = x1
  y3 = y1
  x4 = x2
  y4 = y2
!
!  Require that XMIN <= X.
!
  if ( x3 < xmin .and. x4 < xmin ) then
    ival = -1
    return
  end if

  if ( x3 < xmin .and. xmin <= x4 ) then
    x = xmin
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x3 = x
    y3 = y
    l1 = .true.
  else if ( xmin <= x3 .and. x4 < xmin ) then
    x = xmin
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x4 = x
    y4 = y
    l2 = .true.
  end if
!
!  Require that X <= XMAX.
!
  if ( xmax < x3 .and. xmax < x4 ) then
    ival = -1
    return
  end if

  if ( xmax < x3 .and. x4 <= xmax ) then
    x = xmax
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x3 = x
    y3 = y
    l1 = .true.
  else if ( x3 <= xmax .and. xmax < x4 ) then
    x = xmax
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x4 = x
    y4 = y
    l2 = .true.
  end if
!
!  Require that YMIN <= Y.
!
  if ( y3 < ymin .and. y4 < ymin ) then
    ival = -1
    return
  end if

  if ( y3 < ymin .and. ymin <= y4 ) then
    y = ymin
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y3 = y
    x3 = x
    l1 = .true.
  else if ( ymin <= y3 .and. y4 < ymin ) then
    y = ymin
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y4 = y
    x4 = x
    l2 = .true.
  end if
!
!  Require that Y <= YMAX.
!
  if ( ymax < y3 .and. ymax < y4 ) then
    ival = -1
    return
  end if

  if ( ymax < y3 .and. y4 <= ymax ) then
    y = ymax
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y3 = y
    x3 = x
    l1 = .true.
  else if ( y3 <= ymax .and. ymax < y4 ) then
    y = ymax
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y4 = y
    x4 = x
    l2 = .true.
  end if

  ival = 0

  if ( l1 ) then
    ival = ival + 1
  end if

  if ( l2 ) then
    ival = ival + 2
  end if

  return
end
function box_contains_line_seg_2d ( xmin, ymin, xmax, ymax, x1, y1, x2, y2  )
!
!*******************************************************************************
!
!! BOX_CONTAINS_LINE_SEG_2D reports if a box contains a line segment in 2D.
!
!
!  Discussion:
!
!    The box is assumed to be a rectangle with sides aligned on coordinate
!    axes.
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    18 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XMIN, YMIN, XMAX, YMAX, the minimum and maximum X and Y
!    values, which define the box.
!
!    Input, real X1, Y1, X2, Y2, the coordinates of the endpoints of the
!    line segment.
!
!    Output, logical BOX_CONTAINS_LINE_SEG_2D, is TRUE if the box contains
!    the line segment.
!
  implicit none
!
  logical box_contains_line_seg_2d
  logical box_contains_point_2d
  real x1
  real x2
  real xmax
  real xmin
  real y1
  real y2
  real ymax
  real ymin
!
  if ( box_contains_point_2d ( xmin, ymin, xmax, ymax, x1, y1 ) .and. &
       box_contains_point_2d ( xmin, ymin, xmax, ymax, x2, y2 ) ) then
    box_contains_line_seg_2d = .true.
  else
    box_contains_line_seg_2d = .false.
  end if

  return
end
function box_contains_point_2d ( xmin, ymin, xmax, ymax, x, y )
!
!*******************************************************************************
!
!! BOX_CONTAINS_POINT_2D reports if a box contains a point in 2D.
!
!
!  Discussion:
!
!    The box is assumed to be a rectangle with sides aligned on coordinate
!    axes.
!
!  Modified:
!
!    18 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XMIN, YMIN, XMAX, YMAX, the minimum and maximum X and Y
!    values, which define the box.
!
!    Input, real X, Y, the coordinates of the point.
!
!    Output, logical BOX_CONTAINS_POINT_2D, is TRUE if the box contains
!    the point.
!
  implicit none
!
  logical box_contains_point_2d
  real x
  real xmax
  real xmin
  real y
  real ymax
  real ymin
!
  if ( xmin <= x .and. x <= xmax .and. ymin <= y .and. y <= ymax ) then
    box_contains_point_2d = .true.
  else
    box_contains_point_2d = .false.
  end if

  return
end
subroutine box_ray_int_2d ( xmin, ymin, xmax, ymax, xa, ya, xb, yb, xi, yi )
!
!*******************************************************************************
!
!! BOX_RAY_INT_2D: intersection ( box, ray ) in 2D.
!
!
!  Discussion:
!
!    The box is assumed to be a rectangle with sides aligned on coordinate
!    axes.
!
!    The origin of the ray is assumed to be inside the box.  This
!    guarantees that the ray will intersect the box in exactly one point.
!
!  Modified:
!
!    16 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XMIN, YMIN, the lower left corner of the box.
!
!    Input, real XMAX, YMAX, the upper right corner of the box.
!
!    Input, real XA, YA, the origin of the ray, which should be
!    inside the box.
!
!    Input, real XB, YB, a second point on the ray.
!
!    Output, real XI, YI, the point on the box intersected by the ray.
!
  implicit none
!
  logical inside
  integer ival
  integer side
  real xa
  real xb
  real xc
  real xd
  real xi
  real xmax
  real xmin
  real ya
  real yb
  real yc
  real yd
  real yi
  real ymax
  real ymin
!
  do side = 1, 4

    if ( side == 1 ) then
      xc = xmin
      yc = ymin
      xd = xmax
      yd = ymin
    else if ( side == 2 ) then
      xc = xmax
      yc = ymin
      xd = xmax
      yd = ymax
    else if ( side == 3 ) then
      xc = xmax
      yc = ymax
      xd = xmin
      yd = ymax
    else if ( side == 4 ) then
      xc = xmin
      yc = ymax
      xd = xmin
      yd = ymin
    end if

    call angle_contains_ray_2d ( xc, yc, xa, ya, xd, yd, xb, yb, inside )

    if ( inside ) then
      exit
    end if

    if ( side == 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOX_RAY_INT_2D - Fatal error!'
      write ( *, '(a)' ) '  No intersection could be found.'
      stop
    end if

  end do

  call lines_exp_int_2d ( xa, ya, xb, yb, xc, yc, xd, yd, ival, xi, yi )

  return
end
subroutine ch_cap ( c )
!
!*******************************************************************************
!
!! CH_CAP capitalizes a single character.
!
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none
!
  character c
  integer itemp
!
  itemp = ichar ( c )
 
  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if
 
  return
end
function ch_eqi ( c1, c2 )
!
!*******************************************************************************
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!
!  Examples:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none
!
  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap
!
  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine circle_area_2d ( r, area )
!
!*******************************************************************************
!
!! CIRCLE_AREA_2D computes the area of a circle in 2D.
!
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Output, real AREA, the area of the circle.
!
  implicit none
!
  real area
  real r_pi
  real r
!
  area = r_pi ( ) * r * r

  return
end
subroutine circle_dia2imp_2d ( x1, y1, x2, y2, r, xc, yc )
!
!*******************************************************************************
!
!! CIRCLE_DIA2IMP_2D converts a diameter to an implicit circle in 2D.
!
!
!  Discussion:
!
!    The diameter form of a circle is:
!
!      (X1,Y1) and (X2,Y2) are endpoints of a diameter.
!
!    The implicit form of a circle in 2D is:
!
!      (X-XC)**2 + (Y-YC)**2 = R**2
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, are the X and Y coordinates
!    of two points which form a diameter of the circle.
!
!    Output, real R, the computed radius of the circle.
!
!    Output, real XC, YC, the computed center of the circle.
!
  implicit none
!
  real enorm0_2d
  real r
  real x1
  real x2
  real xc
  real y1
  real y2
  real yc
!
  r = 0.5E+00 * enorm0_2d ( x1, y1, x2, y2 )

  xc = 0.5E+00 * ( x1 + x2 )
  yc = 0.5E+00 * ( y1 + y2 )

  return
end
subroutine circle_exp2imp_2d ( x1, y1, x2, y2, x3, y3, r, xc, yc )
!
!*******************************************************************************
!
!! CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
!
!
!  Formula:
!
!    The explicit form of a circle in 2D is:
!
!      The circle passing through (X1,Y1), (X2,Y2), (X3,Y3).
!
!    The implicit form of a circle in 2D is:
!
!      (X-XC)**2 + (Y-YC)**2 = R**2
!
!  Discussion:
!
!    Any three points define a circle, as long as they don't lie on a straight
!    line.  (If the points do lie on a straight line, we could stretch the
!    definition of a circle to allow an infinite radius and a center at
!    some infinite point.)
!
!    Instead of the formulas used here, you can use the linear system
!    approach in the routine TRIANGLE_OUTCIRCLE_2D.
!
!    The diameter of the circle can be found by solving a 2 by 2 linear system.
!    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
!    and each forms a right triangle with the diameter.  Hence, the dot product
!    of P2 - P1 with the diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
!    diameter vector originating at P1.
!
!  Reference:
!
!    Joseph O'Rourke,
!    Computational Geometry,
!    Cambridge University Press,
!    Second Edition, 1998, page 187.
!
!  Modified:
!
!    12 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, are the X and Y coordinates
!    of three points that lie on the circle.  These points should be
!    distinct, and not collinear.
!
!    Output, real R, the radius of the circle.  Normally, R will be positive.
!    R is returned as -1 in the unlikely event that the points are
!    numerically collinear.
!
!    Output, real XC, YC, the center of the circle.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real e
  real f
  real g
  real r
  real x1
  real x2
  real x3
  real xc
  real y1
  real y2
  real y3
  real yc
!
  a = x2 - x1
  b = y2 - y1
  c = x3 - x1
  d = y3 - y1

  e = a * ( x1 + x2 ) + b * ( y1 + y2 )
  f = c * ( x1 + x3 ) + d * ( y1 + y3 )
!
!  Our formula is:
!
!    G = a * ( d - b ) - b * ( c - a )
!
!  but we get slightly better results using the original data.
!
  g = a * ( y3 - y2 ) - b * ( x3 - x2 )
!
!  We check for collinearity.  A more useful check would compare the
!  absolute value of G to a small quantity.
!
  if ( g == 0.0E+00 ) then
    xc = 0.0E+00
    yc = 0.0E+00
    r = -1.0E+00
    return
  end if
!
!  The center is halfway along the diameter vector from (X1,Y1).
!
  xc = 0.5E+00 * ( d * e - b * f ) / g
  yc = 0.5E+00 * ( a * f - c * e ) / g
!
!  Knowing the center, the radius is now easy to compute.
!
  r = sqrt ( ( x1 - xc )**2 + ( y1 - yc )**2 )

  return
end
subroutine circle_exp_contains_point_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )
!
!*******************************************************************************
!
!! CIRCLE_EXP_CONTAINS_POINT_2D determines if an explicit circle contains a point in 2D.
!
!
!  Formula:
!
!    The explicit form of a circle in 2D is:
!
!      The circle passing through (X1,Y1), (X2,Y2), (X3,Y3).
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the (X,Y) coordinates of three
!    points that lie on a circle.
!
!    Input, real X, Y, the (X,Y) coordinates of a point, whose position
!    relative to the circle is desired.
!
!    Output, integer INSIDE:
!
!   -1, the three points are distinct and noncolinear,
!    and (x,y) lies inside the circle.
!    0, the three points are distinct and noncolinear,
!    and (x,y) lies on the circle.
!    1, the three points are distinct and noncolinear,
!    and (x,y) lies outside the circle.
!
!    2, the three points are distinct and colinear,
!    and (x,y) lies on the line.
!    3, the three points are distinct and colinear,
!    and (x,y) does not lie on the line.
!
!    4, two points are distinct, and (x,y) lies on the line.
!    5, two points are distinct, and (x,y) does not lie on the line.
!
!    6, all three points are equal, and (x,y) is equal to them,
!    7, all three points are equal, and (x,y) is not equal to them.
!
  implicit none
!
  real a(4,4)
  real det
  real rmat4_det
  integer inside
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
!  P1 = P2?
!
  if ( x1 == x2 .and. y1 == y2 ) then

    if ( x1 == x3 .and. y1 == y3 ) then

      if ( x1 == x .and. y1 == y ) then
        inside = 6
      else
        inside = 7
      end if

    else

      det = ( x1 - x3 ) * ( y - y3 ) - ( x - x3 ) * ( y1 - y3 )

      if ( det == 0.0E+00 ) then
        inside = 4
      else
        inside = 5
      end if
    end if

    return

  end if
!
!  P1 does not equal P2.  Does P1 = P3?
!
  if ( x1 == x3 .and. y1 == y3 ) then

    det = ( x1 - x2 ) * ( y - y2 ) - ( x - x2 ) * ( y1 - y2 )

    if ( det == 0.0E+00 ) then
      inside = 4
    else
      inside = 5
    end if

    return

  end if
!
!  The points are distinct.  Are they colinear?
!
  det = ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 )

  if ( det == 0.0E+00 ) then

    det = ( x1 - x2 ) * ( y - y2 ) - ( x - x2 ) * ( y1 - y2 )

    if ( det == 0.0E+00 ) then
      inside = 2
    else
      inside = 3
    end if

    return

  end if
!
!  The points are distinct and non-colinear.
!
!  Compute the determinant
!
  a(1,1) = x1
  a(1,2) = y1
  a(1,3) = x1 * x1 + y1 * y1
  a(1,4) = 1.0E+00

  a(2,1) = x2
  a(2,2) = y2
  a(2,3) = x2 * x2 + y2 * y2
  a(2,4) = 1.0E+00

  a(3,1) = x3
  a(3,2) = y3
  a(3,3) = x3 * x3 + y3 * y3
  a(3,4) = 1.0E+00

  a(4,1) = x
  a(4,2) = y
  a(4,3) = x * x + y * y
  a(4,4) = 1.0E+00

  det = rmat4_det ( a )

  if ( det < 0.0E+00 ) then
    inside = 1
  else if ( det == 0.0E+00 ) then
    inside = 0
  else
    inside = -1
  end if

  return
end
subroutine circle_imp_contains_point_2d ( r, xc, yc, x, y, inside )
!
!*******************************************************************************
!
!! CIRCLE_IMP_CONTAINS_POINT_2D determines if an implicit circle contains a point in 2D.
!
!
!  Formula:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 = R**2
!
!  Modified:
!
!    21 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real XC, YC, the coordinates of the center of the circle.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is inside or on the circle,
!    FALSE otherwise.
!
  implicit none
!
  real enormsq0_2d
  logical inside
  real r
  real x
  real xc
  real y
  real yc
!
  if ( enormsq0_2d ( x, y, xc, yc ) <= r * r ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
subroutine circle_imp_line_par_int_2d ( r, xc, yc, x0, y0, f, g, num_int, x, y )
!
!*******************************************************************************
!
!! CIRCLE_IMP_LINE_PAR_INT_2D: ( implicit circle, parametric line ) intersection in 2D.
!
!
!  Formula:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 = R**2
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real XC, YC, the coordinates of the center of the circle.
!
!    Input, real F, G, X0, Y0, the parametric parameters of the line.
!
!    Output, integer NUM_INT, the number of intersecting points found.
!    NUM_INT will be 0, 1 or 2.
!
!    Output, real X(2), Y(2), contains the X and Y coordinates of
!    the intersecting points.
!
  implicit none
!
  real f
  real g
  integer num_int
  real r
  real root
  real t
  real x(2)
  real x0
  real xc
  real y(2)
  real y0
  real yc
!
  root = r * r * ( f * f + g * g ) - ( f * ( yc - y0 ) - g * ( xc - x0 ) )**2

  if ( root < 0.0E+00 ) then

    num_int = 0

  else if ( root == 0.0E+00 )then

    num_int = 1

    t = ( f * ( xc - x0 ) + g * ( yc - y0 ) ) / ( f * f + g * g )
    x(1) = x0 + f * t
    y(1) = y0 + g * t

  else if ( root > 0.0E+00 ) then

    num_int = 2

    t = ( ( f * ( xc - x0 ) + g * ( yc - y0 ) ) - sqrt ( root ) ) &
      / ( f * f + g * g )

    x(1) = x0 + f * t
    y(1) = y0 + g * t

    t = ( ( f * ( xc - x0 ) + g * ( yc - y0 ) ) + sqrt ( root ) ) &
      / ( f * f + g * g )

    x(2) = x0 + f * t
    y(2) = y0 + g * t

  end if

  return
end
subroutine circle_imp_point_dist_2d ( r, xc, yc, x, y, dist )
!
!*******************************************************************************
!
!! CIRCLE_IMP_POINT_DIST_2D: distance ( implicit circle, point ) in 2D.
!
!
!  Discussion:
!
!    The distance is zero if the point is on the circle.
!
!  Formula:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 = R**2
!
!  Modified:
!
!    03 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real XC, YC, the coordinates of the center of the circle.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, real DIST, the distance of the point to the circle.
!
  implicit none
!
  real dist
  real r
  real x
  real xc
  real y
  real yc
!
  dist = sqrt ( abs ( ( x - xc )**2 + ( y - yc )**2  - r**2 ) )

  return
end
subroutine circle_imp_point_dist_signed_2d ( r, xc, yc, x, y, dist )
!
!*******************************************************************************
!
!! CIRCLE_IMP_POINT_DIST_SIGNED_2D: signed distance ( implicit circle, point ) in 2D.
!
!
!  Discussion:
!
!    The signed distance is zero if the point is on the circle.
!    The signed distance is negative if the point is inside the circle.
!
!  Formula:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 = R**2
!
!  Modified:
!
!    03 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real XC, YC, the coordinates of the center of the circle.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, real DIST, the signed distance of the point to the circle.
!    If the point is inside the circle, the signed distance is negative.
!
  implicit none
!
  real dist
  real r
  real t
  real x
  real xc
  real y
  real yc
!
  t = ( x - xc )**2 + ( y - yc )**2  - r**2

  dist = sign ( 1.0E+00, t ) * sqrt ( abs ( t ) )

  return
end
subroutine circle_imp_points_2d ( r, xc, yc, n, x, y )
!
!*******************************************************************************
!
!! CIRCLE_IMP_POINTS_2D returns N equally spaced points on an implicit circle in 2D.
!
!
!  Note:
!
!    The first point is always ( XC + R, YC ), and subsequent points
!    proceed counterclockwise around the circle.
!
!  Definition:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 = R**2
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real XC, YC, the coordinates of the center of the circle.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real X(N), Y(N), the coordinates of points on the circle.
!
  implicit none
!
  integer n
!
  integer i
  real r_pi
  real r
  real theta
  real x(n)
  real xc
  real y(n)
  real yc
!
  do i = 1, n
    theta = ( 2.0E+00 * r_pi ( ) * real ( i - 1 ) ) / real ( n )
    x(i) = xc + r * cos ( theta )
    y(i) = yc + r * sin ( theta )
  end do

  return
end
subroutine circle_imp_points_arc_2d ( r, xc, yc, theta1, theta2, n, x, y )
!
!*******************************************************************************
!
!! CIRCLE_IMP_POINTS_ARC_2D returns N points on an arc of an implicit circle in 2D.
!
!
!  Discussion:
!
!    The first point is ( XC + R * COS ( THETA1 ), YC + R * SIN ( THETA1 ) );
!    The last point is  ( XC + R * COS ( THETA2 ), YC + R * SIN ( THETA2 ) );
!    and the intermediate points are evenly spaced in angle between these,
!    and in counterclockwise order.
!
!  Definition:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 = R**2
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real XC, YC, the coordinates of the center of the circle.
!
!    Input, real THETA1, THETA2, the angular coordinates of the first
!    and last points to be drawn, in radians.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real X(N), Y(N), the coordinates of points on the circle.
!
  implicit none
!
  integer n
!
  integer i
  real r_pi
  real r
  real r_modp
  real theta
  real theta1
  real theta2
  real theta3
  real x(n)
  real xc
  real y(n)
  real yc
!
!  THETA3 is the smallest angle, no less than THETA1, which
!  coincides with THETA2.
!
  theta3 = theta1 + r_modp ( theta2 - theta1, 2.0E+00 * r_pi ( ) )

  do i = 1, n

    if ( n > 1 ) then
      theta = ( real ( n - i ) * theta1 + real ( i - 1 ) * theta3 ) &
        / real ( n - 1 )
    else
      theta = 0.5E+00 * ( theta1 + theta3 )
    end if

    x(i) = xc + r * cos ( theta )
    y(i) = yc + r * sin ( theta )

  end do

  return
end
subroutine circle_lune_area_2d ( r, theta, area )
!
!*******************************************************************************
!
!! CIRCLE_LUNE_AREA_2D returns the area of a circular lune in 2D.
!
!
!  Discussion:
!
!    A lune is formed by drawing a circular arc, and joining its endpoints.
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real THETA, the angle subtended by the arc.
!
!    Output, real AREA, the area of the lune.
!
  implicit none
!
  real area
  real area_sector
  real area_triangle
  real r
  real theta
!
  call circle_sector_area_2d ( r, theta, area_sector )
  call circle_triangle_area_2d ( r, theta, area_triangle )

  area = area_sector - area_triangle

  return
end
subroutine circle_lune_centroid_2d ( r, xc, yc, theta1, theta2, x, y )
!
!*******************************************************************************
!
!! CIRCLE_LUNE_CENTROID_2D returns the centroid of a circular lune in 2D.
!
!
!  Discussion:
!
!    A lune is formed by drawing a circular arc, and joining its endpoints.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real XC, YC, the coordinates of the center of the circle.
!
!    Input, real THETA1, THETA2, the angles of the first and last points
!    on the circular arc.
!
!    Output, real X, Y, the coordinates of the centroid of the lune.
!
  implicit none
!
  real area
  real d
  real r
  real theta
  real theta1
  real theta2
  real x
  real xc
  real y
  real yc
!
  theta = theta2 - theta1

  if ( theta == 0.0E+00 ) then
    d = r
  else
    d = 4.0E+00 * r * ( sin ( 0.5E+00 * theta ) )**3 / &
      ( 3.0E+00 * ( theta - sin ( theta ) ) )
  end if

  x = xc + d * cos ( theta )
  y = yc + d * sin ( theta )

  return
end
subroutine circle_sector_area_2d ( r, theta, area )
!
!*******************************************************************************
!
!! CIRCLE_SECTOR_AREA_2D computes the area of a circular sector in 2D.
!
!
!  Discussion:
!
!    A  circular sector is formed by a circular arc, and the two straight line 
!    segments that join its ends to the center of the circle.
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real THETA, the angle defining the size of the sector, in radians.
!
!    Output, real AREA, the area of the circle.
!
  implicit none
!
  real area
  real r
  real theta
!
  area = 0.5E+00 * r * r * theta

  return
end
subroutine circle_sector_centroid_2d ( r, xc, yc, theta1, theta2, x, y )
!
!*******************************************************************************
!
!! CIRCLE_SECTOR_CENTROID_2D returns the centroid of a circular sector in 2D.
!
!
!  Discussion:
!
!    A  circular sector is formed by a circular arc, and the two straight line 
!    segments that join its ends to the center of the circle.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real XC, YC, the coordinates of the center of the circle.
!
!    Input, real THETA1, THETA2, the angles of the first and last points
!    on the circular arc.
!
!    Output, real X, Y, the coordinates of the centroid of the sector.
!
  implicit none
!
  real area
  real d
  real r
  real theta
  real theta1
  real theta2
  real x
  real xc
  real y
  real yc
!
  theta = theta2 - theta1

  if ( theta == 0.0E+00 ) then
    d = 2.0E+00 * r / 3.0E+00
  else
    d = 4.0E+00 * r * sin ( 0.5E+00 * theta ) / &
      ( 3.0E+00 * theta )
  end if

  x = xc + d * cos ( theta )
  y = yc + d * sin ( theta )

  return
end
subroutine circle_triangle_area_2d ( r, theta, area )
!
!*******************************************************************************
!
!! CIRCLE_TRIANGLE_AREA_2D returns the area of a circle triangle in 2D.
!
!
!  Discussion:
!
!    A circle triangle is formed by drawing a circular arc, and considering
!    the triangle formed by the endpoints of the arc plus the center of
!    the circle.
!
!    Note that for angles greater than PI, the triangle will actually
!    have NEGATIVE area.
!
!  Modified:
!
!    12 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle.
!
!    Input, real THETA, the angle subtended by the arc.
!
!    Output, real AREA, the (signed) area of the triangle.
!
  implicit none
!
  real area
  real r
  real theta
!
  area = 0.5E+00 * r**2 * sin ( theta )

  return
end
subroutine circles_imp_int_2d ( r1, xc1, yc1, r2, xc2, yc2, num_int, x, y )
!
!*******************************************************************************
!
!! CIRCLES_IMP_INT_2D: finds the intersection of two implicit circles in 2D.
!
!
!  Discussion:
!
!    Two circles can intersect in 0, 1, 2 or infinitely many points.
!
!    The 0 and 2 intersection cases are numerically robust; the 1 and
!    infinite intersection cases are numerically fragile.  The routine
!    uses a tolerance to try to detect the 1 and infinite cases.
!
!  Formula:
!
!    An implicit circle in 2D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 = R**2
!
!  Modified:
!
!    02 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R1, the radius of the first circle.
!
!    Input, real XC1, YC1, the coordinates of the center of the first circle.
!
!    Input, real R2, the radius of the second circle.
!
!    Input, real XC2, YC2, the coordinates of the center of the second circle.
!
!    Output, integer NUM_INT, the number of intersecting points found.
!    NUM_INT will be 0, 1, 2 or 3.  3 indicates that there are an infinite
!    number of intersection points.
!
!    Output, real X(2), Y(2), if NUM_INT is 1 or 2, contains the X and Y
!    coordinates of the intersecting points.
!
  implicit none
!
  real distsq
  integer num_int
  real r1
  real r2
  real root
  real sc1
  real sc2
  real t1
  real t2
  real tol
  real x(2)
  real xc1
  real xc2
  real y(2)
  real yc1
  real yc2
!
  tol = epsilon ( tol )

  x(1:2) = 0.0E+00
  y(1:2) = 0.0E+00
!
!  Take care of the case in which the circles have the same center.
!
  t1 = ( abs ( xc1 - xc2 ) + abs ( yc1 - yc2 ) ) / 2.0E+00
  t2 = ( abs ( xc1 ) + abs ( xc2 ) + abs ( yc1 ) + abs ( yc2 ) + 1.0E+00 ) &
    / 5.0E+00

  if ( t1 <= tol * t2 ) then

    t1 = abs ( r1 - r2 )
    t2 = ( abs ( r1 ) + abs ( r2 ) + 1.0E+00 ) / 3.0E+00

    if ( t1 <= tol * t2 ) then
      num_int = 3
    else
      num_int = 0
    end if

    return

  end if

  distsq = ( xc1 - xc2 )**2 + ( yc1 - yc2 )**2

  root = 2.0E+00 * ( r1**2 + r2**2 )* distsq - distsq**2 &
    - ( r1 - r2 )**2 * ( r1 + r2 )**2

  if ( root < -tol ) then
    num_int = 0
    return
  end if

  sc1 = ( distsq - ( r2**2 - r1**2 ) ) / distsq

  if ( root < tol ) then
    num_int = 1
    x(1) = xc1 + 0.5E+00 * sc1 * ( xc2 - xc1 )
    y(1) = yc1 + 0.5E+00 * sc1 * ( yc2 - yc1 )
    return
  end if

  sc2 = sqrt ( root ) / distsq

  num_int = 2

  x(1) = xc1 + 0.5E+00 * sc1 * ( xc2 - xc1 ) - 0.5E+00 * sc2 * ( yc2 - yc1 )
  y(1) = yc1 + 0.5E+00 * sc1 * ( yc2 - yc1 ) + 0.5E+00 * sc2 * ( xc2 - xc1 )

  x(2) = xc1 + 0.5E+00 * sc1 * ( xc2 - xc1 ) + 0.5E+00 * sc2 * ( yc2 - yc1 )
  y(2) = yc1 + 0.5E+00 * sc1 * ( yc2 - yc1 ) - 0.5E+00 * sc2 * ( xc2 - xc1 )

  return
end
subroutine cone_area_3d ( h, r, area )
!
!*******************************************************************************
!
!! CONE_AREA_3D computes the surface area of a right circular cone in 3D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real H, R, the height of the cone, and the radius of the
!    circle that forms the base of the cone.
!
!    Output, real AREA, the surface area of the cone.
!
  implicit none
!
  real area
  real h
  real r_pi
  real r
!
  area = r_pi ( ) * r * sqrt ( h * h + r * r )

  return
end
subroutine cone_centroid_3d ( r, xc, yc, zc, xh, yh, zh, x, y, z )
!
!*******************************************************************************
!
!! CONE_CENTROID_2D returns the centroid of a cone in 3D.
!
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the circle at the base of the cone.
!
!    Input, real XC, YC, ZC, the coordinates of the center of the circle.
!
!    Input, real XH, YH, ZH, the coordinates of the tip of the cone.
!
!    Output, real X, Y, Z, the coordinates of the centroid of the cone.
!
  implicit none
!
  real r
  real x
  real xc
  real xh
  real y
  real yc
  real yh
  real z
  real zc
  real zh
!
  x = 0.75E+00 * xc + 0.25E+00 * xh
  y = 0.75E+00 * yc + 0.25E+00 * yh
  z = 0.75E+00 * zc + 0.25E+00 * zh

  return
end
subroutine cone_volume_3d ( h, r, volume )
!
!*******************************************************************************
!
!! CONE_VOLUME_3D computes the volume of a right circular cone in 3D.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real H, R, the height of the cone, and the radius of the
!    circle that forms the base of the cone.
!
!    Output, real VOLUME, the volume of the cone.
!
  implicit none
!
  real h
  real r_pi
  real r
  real volume
!
  volume = r_pi ( ) * r * r * h / 3.0E+00

  return
end
subroutine conv3d ( cor, cor2, cor3, maxcor2, maxcor3, ncor2, theta )
!
!*******************************************************************************
!
!! CONV3D converts 3D data to a 2D projection.
!
!
!  Discussion:
!
!    A "presentation angle" THETA is used to project the 3D point
!    (X3D, Y3D, Z3D) to the 2D projection (XVAL,YVAL).
!
!  Formula:
!
!    If COR = 'X':
!
!      X2D = Y3D - sin ( THETA ) * X3D
!      Y2D = Z3D - sin ( THETA ) * X3D
!
!  Modified:
!
!    01 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character COR, the coordinate to be projected.
!    COR should be 'X', 'Y', or 'Z'.
!
!    Output, real COR2(2,NCOR), the 2D projections.
!
!    Input, real COR3(3,MAXCOR3), the X, Y, and Z components of points.
!
!    Input, integer MAXCOR2, MAXCOR3, the maximum number of 2D and
!    3D points allowed.
!
!    Input, integer NCOR2, the number of 2D values to be computed.
!
!    Input, real THETA, the presentation angle in degrees.
!
  implicit none
!
  integer maxcor2
  integer maxcor3
!
  logical ch_eqi
  character cor
  real cor2(2,maxcor2)
  real cor3(3,maxcor3)
  real degrees_to_radians
  integer i
  integer ncor2
  real stheta
  real theta
!
  stheta = sin ( degrees_to_radians ( theta ) )

  if ( ch_eqi ( cor, 'X' ) ) then

    cor2(1,1:ncor2) = cor3(2,1:ncor2) - stheta * cor3(1,1:ncor2)
    cor2(2,1:ncor2) = cor3(3,1:ncor2) - stheta * cor3(1,1:ncor2)

  else if ( ch_eqi ( cor, 'Y' ) ) then

    cor2(1,1:ncor2) = cor3(1,1:ncor2) - stheta * cor3(2,1:ncor2)
    cor2(2,1:ncor2) = cor3(3,1:ncor2) - stheta * cor3(2,1:ncor2)

  else if ( ch_eqi ( cor, 'Z' ) ) then

    cor2(1,1:ncor2) = cor3(1,1:ncor2) - stheta * cor3(3,1:ncor2)
    cor2(2,1:ncor2) = cor3(2,1:ncor2) - stheta * cor3(3,1:ncor2)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CONV3D - Fatal error!'
    write ( *, '(a)' ) '  Illegal coordinate index = ' // cor
    stop

  end if

  return
end
subroutine corpl_2d ( dist, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5 )
!
!*******************************************************************************
!
!! CORPL_2D "boxes" an angle defined by three points in 2D.
!
!
!  Discussion:
!
!    CORPL2 is given
!      (X1,Y1), (X2,Y2), and (X3,Y3)
!    determining the two lines:
!      (X1,Y1), (X2,Y2)
!    and
!      (X2,Y2), (X3,Y3)
!    and a nonnegative distance
!      DIST.
!
!    CORPL2 returns a pair of "corner" points
!      (X4,Y4) and (X5,Y5)
!    both of which are a distance DIST from both lines, and in fact,
!    both of which are a distance DIST from (X2,Y2).
!
!                         /   3
!                        /   /   /
!     - - - - - - - - - 4 - / - 6 - - -
!                      /   /   /
!     1---------------/---2-----------------
!                    /   /   /
!     - - - - - - - 7 - / - 5 - - - - -
!                  /   /   /
!
!    In the illustration, the numbers "1" "2" and "3" represent
!    the points defining the lines.
!
!    The numbers "4" and "5" represent the desired "corner points", which
!    are on the positive or negative sides of both lines.
!
!    The numbers "6" and "7" represent the undesired points, which
!    are on the positive side of one line and the negative of the other.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real DIST, the nonnegative distance from (X1,Y1) to the computed
!    points (X3,Y3) and (X4,Y4).
!
!    Input, real X1, Y1, X2, Y2, X3, Y3.
!    (X1,Y1) and (X2,Y2) are distinct points that define a line.
!    (X2,Y2) and (X3,Y3) are distinct points that define a line.
!
!    Special cases:
!
!    if ( X1,Y1) = (X2,Y2), this is the same as extending the line from
!    (X3,Y3) through (X2,Y2) without a bend.
!
!    if ( X3,Y3) = (X2,Y2), this is the same as extending the line from
!    (X1,Y1) through (X2,Y2) without a bend.
!
!    if ( X1,Y1) = (X2,Y2) = (X3,Y3) this is an error.
!
!    Output, real X4, Y4, X5, Y5.
!    (X4,Y4) and (X5,Y5) are points which lie DIST units from
!    the line between (X1,Y1) and (X2,Y2), and DIST units from
!    the line between (X2,Y2) and (X3,Y3).
!
  implicit none
!
  real dist
  real stheta
  real temp
  real ux
  real ux1
  real ux2
  real uy
  real uy1
  real uy2
  real x1
  real x1copy
  real x2
  real x3
  real x3copy
  real x4
  real x5
  real y1
  real y1copy
  real y2
  real y3
  real y3copy
  real y4
  real y5
!
!  Fail if all three points are equal.
!
  if ( x1 == x2 .and. x2 == x3 .and. y1 == y2 .and. y2 == y3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CORPL2 - Fatal error!'
    write ( *, '(a)' ) '  Input points (X1,Y1) = (X2,Y2) = (X3,Y3).'
    write ( *, '(a,2g14.6)' ) '  (X1,Y1)= ', x1, y1
    stop
  end if
!
!  If P1 = P2 or P2 = P3, extend the line through the doubled point.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    x1copy = 2.0E+00 * x2 - x3
    y1copy = 2.0E+00 * y2 - y3
  else
    x1copy = x1
    y1copy = y1
  end if

  if ( x2 == x3 .and. y2 == y3 ) then
    x3copy = 2.0E+00 * x2 - x1
    y3copy = 2.0E+00 * y2 - y1
  else
    x3copy = x3
    y3copy = y3
  end if
!
!  Now compute the unit normal vectors to each line.
!  We choose the sign so that the unit normal to line 1 has
!  a positive dot product with line 2.
!
  ux1 = y1copy - y2
  uy1 = x2 - x1copy
  temp = sqrt ( ux1 * ux1 + uy1 * uy1 )
  ux1 = ux1 / temp
  uy1 = uy1 / temp

  if ( ux1 * ( x3copy - x2 ) + uy1 * ( y3copy - y2 ) < 0.0E+00 ) then
    ux1 = - ux1
    uy1 = - uy1
  end if

  ux2 = y3copy - y2
  uy2 = x2 - x3copy
  temp = sqrt ( ux2 * ux2 + uy2 * uy2)
  ux2 = ux2 / temp
  uy2 = uy2 / temp

  if ( ux2 * ( x1copy - x2 ) + uy2 * ( y1copy - y2 ) < 0.0E+00 ) then
    ux2 = - ux2
    uy2 = - uy2
  end if
!
!  Try to catch the case where we can't determine the
!  sign of U1, because both U1 and -U1 are perpendicular
!  to (P3-P2)...and similarly for U2 and (P1-P2).
!
  if ( ux1 * ( x3copy - x2) + uy1 * ( y3copy - y2 ) == 0.0E+00 .or. &
       ux2 * ( x1copy - x2) + uy2 * ( y1copy - y2 ) == 0.0E+00 ) then

    if ( ux1 * ux2 + uy1 * uy2 < 0.0E+00 ) then
      ux1 = -ux1
      uy1 = -uy1
    end if

  end if
!
!  Try to catch a line turning back on itself, evidenced by
!    Cos(theta) = (P3-P2) dot (P2-P1) / ( norm(P3-P2) * norm(P2-P1) )
!  being -1, or very close to -1.
!
  temp = ( ( x3copy - x2 ) * ( x2 - x1copy ) &
         + ( y3copy - y2 ) * ( y2 - y1copy ) )

  temp = temp / &
        ( sqrt ( ( x3copy - x2 )**2 + ( y3copy - y2 )**2 ) &
        * sqrt ( ( x2 - x1copy )**2 + ( y2 - y1copy )**2 ) )

  if ( temp < -0.99E+00 ) then
    temp = sqrt ( ( x2 - x1copy)**2 + ( y2 - y1copy )**2 )
    x4 = x2 + dist * ( x2 - x1copy ) / temp + dist * ux1
    y4 = y2 + dist * ( y2 - y1copy ) / temp + dist * uy1
    x5 = x2 + dist * ( x2 - x1copy ) / temp - dist * ux1
    y5 = y2 + dist * ( y2 - y1copy ) / temp - dist * uy1
    return
  end if
!
!  Compute the "average" unit normal vector.
!
!  The average of the unit normals could be zero, but only when
!  the second line has the same direction and opposite sense
!  of the first, and we've already checked for that case.
!
  ux = 0.5E+00 * ( ux1 + ux2 )
  uy = 0.5E+00 * ( uy1 + uy2 )
  temp = sqrt ( ux * ux + uy * uy )
  ux = ux / temp
  uy = uy / temp
!
!  You must go DIST/STHETA units along this unit normal to
!  result in a distance DIST from line1 (and line2).
!
  stheta = ux * ux1 + uy * uy1

  x4 = x2 + dist * ux / stheta
  y4 = y2 + dist * uy / stheta

  x5 = x2 - dist * ux / stheta
  y5 = y2 - dist * uy / stheta

  return
end
function cot ( angle )
!
!*******************************************************************************
!
!! COT returns the cotangent of an angle.
!
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle, in radians.
!
!    Output, real COT, the cotangent of the angle.
!
  implicit none
!
  real angle
  real cot
!
  cot  = cos ( angle ) / sin ( angle )

  return
end
function cot_degrees ( angle )
!
!*******************************************************************************
!
!! COT_DEGREES returns the cotangent of an angle given in degrees.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle, in degrees.
!
!    Output, real COT_DEGREES, the cotangent of the angle.
!
  implicit none
!
  real angle
  real angle_rad
  real cot_degrees
  real degrees_to_radians
!
  angle_rad = degrees_to_radians ( angle )

  cot_degrees  = cos ( angle_rad ) / sin ( angle_rad )

  return
end
function cross0_2d ( x0, y0, x1, y1, x2, y2 )
!
!*******************************************************************************
!
!! CROSS0_2D finds the cross product of (P1-P0) and (P2-P0) in 2D.
!
!
!  Discussion:
!
!    Strictly speaking, the vectors lie in the (X,Y) plane, and
!    the cross product here is a vector in the Z direction.
!
!    The vectors are specified with respect to a basis point P0.
!    We are computing the normal to the triangle (P0,P1,P2).
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real X3, Y3, Z3, the cross product (P1-P0) x (P2-P0), or
!    (X1-X0,Y1-Y0,Z1-Z0) x (X2-X0,Y2-Y0,Z2-Z0).
!
!    Input, real X0, Y0, X1, Y1, X2, Y2, the coordinates of the three
!    points.  The basis point P0 is (X0,Y0).
!
!    Output, real CROSS0_2D, the Z component of the cross product
!    (X1-X0,Y1-Y0,0) x (X2-X0,Y2-Y0,0).
!
  implicit none
!
  real cross0_2d
  real x0
  real x1
  real x2
  real y0
  real y1
  real y2
!
  cross0_2d = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )

  return
end
subroutine cross0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! CROSS0_3D computes the cross product of (P1-P0) and (P2-P0) in 3D.
!
!
!  Discussion:
!
!    The vectors are specified with respect to a basis point P0.
!    We are computing the normal to the triangle (P0,P1,P2).
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, the coordinates of
!    three points.  The basis point is (X0,Y0,Z0).
!
!    Output, real X3, Y3, Z3, the cross product (P1-P0) x (P2-P0), or
!    (X1-X0,Y1-Y0,Z1-Z0) x (X2-X0,Y2-Y0,Z2-Z0).
!
  implicit none
!
  real x0
  real x1
  real x2
  real x3
  real y0
  real y1
  real y2
  real y3
  real z0
  real z1
  real z2
  real z3
!
  x3 = ( y1 - y0 ) * ( z2 - z0 ) - ( z1 - z0 ) * ( y2 - y0 )
  y3 = ( z1 - z0 ) * ( x2 - x0 ) - ( x1 - x0 ) * ( z2 - z0 )
  z3 = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )

  return
end
function cross_2d ( x1, y1, x2, y2 )
!
!*******************************************************************************
!
!! CROSS_2D finds the cross product of a pair of vectors in 2D.
!
!
!  Discussion:
!
!    Strictly speaking, the vectors lie in the (X,Y) plane, and
!    the cross product here is a vector in the Z direction.
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the X and Y coordinates of the vectors.
!
!    Output, real CROSS_2D, the Z component of the cross product
!    of (X1,Y1,0) and (X2,Y2,0).
!
  implicit none
!
  real cross_2d
  real x1
  real x2
  real y1
  real y2
!
  cross_2d = x1 * y2 - y1 * x2

  return
end
subroutine cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! CROSS_3D computes the cross product of two vectors in 3D.
!
!
!  Definition:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
!
!    Output, real X3, Y3, Z3, the cross product vector.
!
  implicit none
!
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  x3 = y1 * z2 - z1 * y2
  y3 = z1 * x2 - x1 * z2
  z3 = x1 * y2 - y1 * x2

  return
end
subroutine cube_shape_3d ( max_num, max_order, point_num, face_num, &
  face_order, point_coord, face_point )
!
!*******************************************************************************
!
!! CUBE_SHAPE_3D describes a cube in 3D.
!
!
!  Discussion:
!
!    The vertices of the cube lie on the unit sphere.
!
!  Modified:
!
!    11 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAX_NUM, the maximum number of faces, and of
!    points.  This quantity is used to dimension arrays.
!    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
!    not be set.
!
!    Input, integer MAX_ORDER, the maximum order of any face.
!    This quantity is used to dimension arrays.  If FACE_ORDER
!    exceeds MAX_ORDER, the arrays will not be set.
!
!    Output, integer POINT_NUM, the number of points in the shape.
!
!    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
!    the X, Y and Z coordinates of point J.
!
!    Output, integer FACE_NUM, the number of faces in the shape.
!
!    Output, integer FACE_ORDER, the number of vertices per face.
!
!    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
!    contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter-clockwise direction defined
!    by the outward normal at the face.
!
  implicit none
!
  integer max_num
  integer max_order
!
  real a
  integer face_num
  integer face_order
  integer face_point(max_order,max_num)
  integer point_num
  real point_coord(3,max_num)
!
  point_num = 8
  face_num = 6
  face_order = 4
!
!  Check.
!
  if ( point_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBE_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of vertices exceeds MAX_NUM.'
    return
  end if

  if ( face_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBE_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of faces exceeds MAX_NUM.'
    return
  end if

  if ( face_order > max_order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBE_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Face order exceeds MAX_ORDER.'
    return
  end if
!
!  Set point coordinates.
!
  a = sqrt ( 1.0E+00 / 3.0E+00 )

  point_coord(1,1) = - a
  point_coord(2,1) = - a
  point_coord(3,1) = - a

  point_coord(1,2) =   a
  point_coord(2,2) = - a
  point_coord(3,2) = - a

  point_coord(1,3) =   a
  point_coord(2,3) =   a
  point_coord(3,3) = - a

  point_coord(1,4) = - a
  point_coord(2,4) =   a
  point_coord(3,4) = - a

  point_coord(1,5) = - a
  point_coord(2,5) = - a
  point_coord(3,5) =   a

  point_coord(1,6) =   a
  point_coord(2,6) = - a
  point_coord(3,6) =   a

  point_coord(1,7) =   a
  point_coord(2,7) =   a
  point_coord(3,7) =   a

  point_coord(1,8) = - a
  point_coord(2,8) =   a
  point_coord(3,8) =   a
!
!  Set faces.
!
  face_point(1,1) = 1
  face_point(2,1) = 4
  face_point(3,1) = 3
  face_point(4,1) = 2

  face_point(1,2) = 1
  face_point(2,2) = 2
  face_point(3,2) = 6
  face_point(4,2) = 5

  face_point(1,3) = 2
  face_point(2,3) = 3
  face_point(3,3) = 7
  face_point(4,3) = 6

  face_point(1,4) = 3
  face_point(2,4) = 4
  face_point(3,4) = 8
  face_point(4,4) = 7

  face_point(1,5) = 1
  face_point(2,5) = 5
  face_point(3,5) = 8
  face_point(4,5) = 4

  face_point(1,6) = 5
  face_point(2,6) = 6
  face_point(3,6) = 7
  face_point(4,6) = 8

  return
end
function degrees_to_radians ( angle )
!
!*******************************************************************************
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, an angle in degrees.
!
!    Output, real DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  implicit none
!
  real angle
  real degrees_to_radians
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
!
  degrees_to_radians = ( angle / 180.0E+00 ) * pi

  return
end
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )
!
!*******************************************************************************
!
!! DIAEDG chooses a diagonal edge.
!
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge
!    that should be chosen, based on the circumcircle criterion, where
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!    quadrilateral in counterclockwise order.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, real X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none
!
  real ca
  real cb
  integer diaedg
  real dx10
  real dx12
  real dx30
  real dx32
  real dy10
  real dy12
  real dy30
  real dy32
  real s
  real tol
  real tola
  real tolb
  real x0
  real x1
  real x2
  real x3
  real y0
  real y1
  real y2
  real y3
!
  tol = 100.0E+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( ca > tola .and. cb > tolb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( s > tola ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end
subroutine direction_pert_3d ( sigma, vbase, vran )
!
!*******************************************************************************
!
!! DIRECTION_PERT_3D randomly perturbs a direction vector in 3D.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real SIGMA, determines the strength of the perturbation.
!    SIGMA <= 0 results in a completely random direction.
!    1 <= SIGMA results in VBASE.
!    0 < SIGMA < 1 results in a perturbation from VBASE, which is
!    large when SIGMA is near 0, and small when SIGMA is near 1.
!
!    Input, real VBASE(3), the base direction vector, which should have
!    unit norm.
!
!    Output, real VRAN(3), the perturbed vector, which will have unit norm.
!
  implicit none
!
  real arc_cosine
  real dphi
  integer i
  real phi
  real r_pi
  real psi
  real r
  real sigma
  real theta
  real v(3)
  real vbase(3)
  real vdot
  real vran(3)
  real x
!
!  SIGMA >= 0, just use the base vector.
!
  if ( sigma >= 1.0E+00 ) then

    vran(1:3) = vbase(1:3)

  else if ( sigma <= 0.0E+00 ) then

    call random_number ( harvest = vdot )
    vdot = 2.0E+00 * vdot - 1.0E+00

    phi = arc_cosine ( vdot )

    call random_number ( harvest = theta )
    theta = 2.0E+00 * r_pi ( ) * theta

    vran(1) = cos ( theta ) * sin ( phi )
    vran(2) = sin ( theta ) * sin ( phi )
    vran(3) = cos ( phi )

  else

    phi = arc_cosine ( vbase(3) )
    theta = atan2 ( vbase(2), vbase(1) )
!
!  Pick VDOT, which must be between -1 and 1.  This represents
!  the dot product of the perturbed vector with the base vector.
!
!  RANDOM_NUMBER returns a uniformly random value between 0 and 1.
!  The operations we perform on this quantity tend to bias it
!  out towards 1, as SIGMA grows from 0 to 1.
!
!  VDOT, in turn, is a value between -1 and 1, which, for large
!  SIGMA, we want biased towards 1.
!
    call random_number ( harvest = r )
    x = exp ( ( 1.0E+00 - sigma ) * log ( r ) )
    dphi = arc_cosine ( 2.0E+00 * x - 1.0E+00 )
!
!  Now we know enough to write down a vector that is rotated DPHI
!  from the base vector.
!
    v(1) = cos ( theta ) * sin ( phi + dphi )
    v(2) = sin ( theta ) * sin ( phi + dphi )
    v(3) = cos ( phi + dphi )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the base vector.
!
    call random_number ( harvest = psi )
    psi = 2.0E+00 * r_pi ( ) * psi
!
!  Carry out the rotation.
!
    call rotation_axis_vector_3d ( vbase, psi, v, vran )

  end if

  return
end
subroutine direction_random_2d ( vran )
!
!*******************************************************************************
!
!! DIRECTION_RANDOM_2D picks a random direction vector in 2D.
!
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VRAN(2), the random direction vector, with unit norm.
!
  implicit none
!
  real r_pi
  real theta
  real vran(2)
!
  call random_number ( harvest = theta )
  theta = 2.0E+00 * r_pi ( ) * theta

  vran(1) = cos ( theta )
  vran(2) = sin ( theta )

  return
end
subroutine direction_random_3d ( vran )
!
!*******************************************************************************
!
!! DIRECTION_RANDOM_3D picks a random direction vector in 3D.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real VRAN(3), the random direction vector, with unit norm.
!
  implicit none
!
  real arc_cosine
  real phi
  real r_pi
  real r
  real theta
  real vdot
  real vran(3)
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  call random_number ( harvest = vdot )
  vdot = 2.0E+00 * vdot - 1.0E+00

  phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  call random_number ( harvest = theta )
  theta = 2.0E+00 * r_pi ( ) * theta

  vran(1) = cos ( theta ) * sin ( phi )
  vran(2) = sin ( theta ) * sin ( phi )
  vran(3) = cos ( phi )

  return
end
subroutine direction_random_nd ( n, w )
!
!*******************************************************************************
!
!! DIRECTION_RANDOM_ND generates a random direction vector in ND.
!
!
!  Modified:
!
!    13 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Output, real W(N), a random direction vector, with unit norm.
!
  implicit none
!
  integer n
!
  integer i
  real norm
  real w(n)
!
!  Get N values from a standard normal distribution.
!
  call normal_01_vector ( n, w )
!
!  Compute the length of the vector.
!
  norm = sqrt ( sum ( w(1:n)**2 ) )
!
!  Normalize the vector.
!
  w(1:n) = w(1:n) / norm

  return
end
subroutine dms_to_radians ( degrees, minutes, seconds, radians )
!
!*******************************************************************************
!
!! DMS_TO_RADIANS converts an angle from degrees/minutes/seconds to radians.
!
!
!  Modified:
!
!    05 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DEGREES, MINUTES, SECONDS, an angle in degrees, minutes,
!    and seconds.
!
!    Output, real RADIANS, the equivalent angle in radians.
!
  implicit none
!
  real angle
  integer degrees
  integer minutes
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real radians
  integer seconds
!
  angle =   real ( degrees ) &
        + ( real ( minutes ) &
        + ( real ( seconds ) / 60.0E+00 ) ) / 60.0E+00

  radians = ( angle / 180.0E+00 ) * pi

  return
end
subroutine dodec_shape_3d ( max_num, max_order, point_num, face_num, &
  face_order, point_coord, face_point )
!
!*******************************************************************************
!
!! DODEC_SHAPE_3D describes an dodecahedron in 3D.
!
!
!  Discussion:
!
!    The vertices of the dodecahedron lie on the unit sphere.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAX_NUM, the maximum number of faces, and of
!    points.  This quantity is used to dimension arrays.
!    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
!    not be set.
!
!    Input, integer MAX_ORDER, the maximum order of any face.
!    This quantity is used to dimension arrays.  If FACE_ORDER
!    exceeds MAX_ORDER, the arrays will not be set.
!
!    Output, integer POINT_NUM, the number of points in the shape.
!
!    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
!    the X, Y and Z coordinates of point J.
!
!    Output, integer FACE_NUM, the number of faces in the shape.
!
!    Output, integer FACE_ORDER, the number of vertices per face.
!
!    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
!    contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter-clockwise direction defined
!    by the outward normal at the face.
!
  implicit none
!
  integer max_num
  integer max_order
!
  real a
  real b
  real c
  integer face_num
  integer face_order
  integer face_point(max_order,max_num)
  real phi
  integer point_num
  real point_coord(3,max_num)
!
  point_num = 20
  face_num = 12
  face_order = 5
!
!  Check.
!
  if ( point_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DODEC_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of vertices exceeds MAX_NUM.'
    return
  end if

  if ( face_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DODEC_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of faces exceeds MAX_NUM.'
    return
  end if

  if ( face_order > max_order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DODEC_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Face order exceeds MAX_ORDER.'
    return
  end if
!
!  Set point coordinates.
!
  phi = 0.5E+00 * ( sqrt ( 5.0E+00 ) + 1.0E+00 )

  a = 1.0E+00 / sqrt ( 3.0E+00 )
  b = phi / sqrt ( 3.0E+00 )
  c = ( phi - 1.0E+00 ) / sqrt ( 3.0E+00 )

  point_coord(1,1) =   a
  point_coord(2,1) =   a
  point_coord(3,1) =   a

  point_coord(1,2) =   a
  point_coord(2,2) =   a
  point_coord(3,2) = - a

  point_coord(1,3) =   a
  point_coord(2,3) = - a
  point_coord(3,3) =   a

  point_coord(1,4) =   a
  point_coord(2,4) = - a
  point_coord(3,4) = - a

  point_coord(1,5) = - a
  point_coord(2,5) =   a
  point_coord(3,5) =   a

  point_coord(1,6) = - a
  point_coord(2,6) =   a
  point_coord(3,6) = - a

  point_coord(1,7) = - a
  point_coord(2,7) = - a
  point_coord(3,7) =   a

  point_coord(1,8) = - a
  point_coord(2,8) = - a
  point_coord(3,8) = - a

  point_coord(1,9) =    c
  point_coord(2,9) =    b
  point_coord(3,9) =    0.0E+00

  point_coord(1,10) = - c
  point_coord(2,10) =   b
  point_coord(3,10) =   0.0E+00

  point_coord(1,11) =   c
  point_coord(2,11) = - b
  point_coord(3,11) =   0.0E+00

  point_coord(1,12) = - c
  point_coord(2,12) = - b
  point_coord(3,12) =   0.0E+00

  point_coord(1,13) =   b
  point_coord(2,13) =   0.0E+00
  point_coord(3,13) =   c

  point_coord(1,14) =   b
  point_coord(2,14) =   0.0E+00
  point_coord(3,14) = - c

  point_coord(1,15) = - b
  point_coord(2,15) =   0.0E+00
  point_coord(3,15) =   c

  point_coord(1,16) = - b
  point_coord(2,16) =   0.0E+00
  point_coord(3,16) = - c

  point_coord(1,17) =   0.0E+00
  point_coord(2,17) =   c
  point_coord(3,17) =   b

  point_coord(1,18) =   0.0E+00
  point_coord(2,18) = - c
  point_coord(3,18) =   b

  point_coord(1,19) =   0.0E+00
  point_coord(2,19) =   c
  point_coord(3,19) = - b

  point_coord(1,20) =   0.0E+00
  point_coord(2,20) = - c
  point_coord(3,20) = - b
!
!  Set faces.
!
  face_point(1,1) = 2
  face_point(2,1) = 9
  face_point(3,1) = 1
  face_point(4,1) = 13
  face_point(5,1) = 14

  face_point(1,2) = 5
  face_point(2,2) = 10
  face_point(3,2) = 6
  face_point(4,2) = 16
  face_point(5,2) = 15

  face_point(1,3) = 3
  face_point(2,3) = 11
  face_point(3,3) = 4
  face_point(4,3) = 14
  face_point(5,3) = 13

  face_point(1,4) = 8
  face_point(2,4) = 12
  face_point(3,4) = 7
  face_point(4,4) = 15
  face_point(5,4) = 16

  face_point(1,5) = 3
  face_point(2,5) = 13
  face_point(3,5) = 1
  face_point(4,5) = 17
  face_point(5,5) = 18

  face_point(1,6) = 2
  face_point(2,6) = 14
  face_point(3,6) = 4
  face_point(4,6) = 20
  face_point(5,6) = 19

  face_point(1,7) = 5
  face_point(2,7) = 15
  face_point(3,7) = 7
  face_point(4,7) = 18
  face_point(5,7) = 17

  face_point(1,8) = 8
  face_point(2,8) = 16
  face_point(3,8) = 6
  face_point(4,8) = 19
  face_point(5,8) = 20

  face_point(1,9) = 5
  face_point(2,9) = 17
  face_point(3,9) = 1
  face_point(4,9) = 9
  face_point(5,9) = 10

  face_point(1,10) = 3
  face_point(2,10) = 18
  face_point(3,10) = 7
  face_point(4,10) = 12
  face_point(5,10) = 11

  face_point(1,11) = 2
  face_point(2,11) = 19
  face_point(3,11) = 6
  face_point(4,11) = 10
  face_point(5,11) = 9

  face_point(1,12) = 8
  face_point(2,12) = 20
  face_point(3,12) = 4
  face_point(4,12) = 11
  face_point(5,12) = 12

  return
end
function dot0_2d ( x0, y0, x1, y1, x2, y2 )
!
!*******************************************************************************
!
!! DOT0_2D computes the dot product of (P1-P0) and (P2-P0) in 2D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, X1, Y1, X2, Y2, the coordinates of
!    the points P0, P1 and P1.
!
!    Output, real DOT0_2D, the dot product of (P1-P0) and (P2-P0).
!
  implicit none
!
  real dot0_2d
  real x0
  real x1
  real x2
  real y0
  real y1
  real y2
!
  dot0_2d = ( x1 - x0 ) * ( x2 - x0 ) + ( y1 - y0 ) * ( y2 - y0 )

  return
end
function dot0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2 )
!
!*******************************************************************************
!
!! DOT0_3D computes the dot product of (P1-P0) and (P2-P0) in 3D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, the coordinates of the point P0.
!
!    Input, real X1, Y1, Z1, the coordinates of the point P1.
!
!    Input, real X2, Y2, Z2, the coordinates of the point P2.
!
!    Output, real DOT0_3D, the dot product of (P1-P0) and (P2-P0).
!
  implicit none
!
  real dot0_3d
  real x0
  real x1
  real x2
  real y0
  real y1
  real y2
  real z0
  real z1
  real z2
!
  dot0_3d = ( x1 - x0 ) * ( x2 - x0 ) + ( y1 - y0 ) * ( y2 - y0 ) + &
            ( z1 - z0 ) * ( z2 - z0 )

  return
end
function dot_2d ( x1, y1, x2, y2 )
!
!*******************************************************************************
!
!! DOT_2D computes the dot product of a pair of vectors in 2D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the coordinates of the vectors.
!
!    Output, real DOT_2D, the dot product.
!
  real dot_2d
  real x1
  real x2
  real y1
  real y2
!
  dot_2d = x1 * x2 + y1 * y2

  return
end
function dot_3d ( x1, y1, z1, x2, y2, z2 )
!
!*******************************************************************************
!
!! DOT_3D computes the dot product of a pair of vectors in 3D.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
!
!    Output, real DOT_3D, the dot product.
!
  implicit none
!
  real dot_3d
  real x1
  real x2
  real y1
  real y2
  real z1
  real z2
!
  dot_3d = x1 * x2 + y1 * y2 + z1 * z2

  return
end
function dot_nd ( n, v1, v2 )
!
!*******************************************************************************
!
!! DOT_ND computes the dot product of a pair of vectors in ND.
!
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real V1(N), V2(N), the two vectors.
!
!    Output, real DOT_ND, the dot product.
!
  implicit none
!
  integer n
!
  real dot_nd
  integer i
  real v1(n)
  real v2(n)
!
  dot_nd = dot_product ( v1(1:n), v2(1:n) )

  return
end
subroutine dual_shape_3d ( max_num, max_order, point_num1, face_num1, &
  face_order1, point_coord1, face_point1, point_num2, face_num2, face_order2, &
  point_coord2, face_point2 )
!
!*******************************************************************************
!
!! DUAL_SHAPE_3D constructs the dual of a shape in 3D.
!
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAX_NUM, the maximum number of faces, and of
!    points.  This quantity is used to dimension arrays.
!    If POINT_NUM1 or FACE_NUM1 exceeds MAX_NUM, the arrays will
!    not be set.
!
!    Input, integer MAX_ORDER, the maximum order of any face.
!    This quantity is used to dimension arrays.  If FACE_ORDER1
!    or FACE_ORDER2 exceeds MAX_ORDER, the arrays will not be set.
!
!    Input, integer POINT_NUM1, the number of points in the shape.
!
!    Input, integer FACE_NUM1, the number of faces in the shape.
!
!    Input, integer FACE_ORDER1, the number of vertices per face.
!
!    Input, real POINT_COORD1(3,MAX_NUM); POINT_COORD1(*,J) contains
!    the X, Y and Z coordinates of point J.
!
!    Input, integer FACE_POINT1(MAX_ORDER,MAX_NUM); FACE_POINT1(I,J)
!    contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter-clockwise direction defined
!    by the outward normal at the face.
!
!    Output, integer POINT_NUM2, the number of points in the dual.
!
!    Output, integer FACE_NUM2, the number of faces in the dual.
!
!    Output, integer FACE_ORDER2, the number of vertices per face
!    in the dual.
!
!    Output, real POINT_COORD2(3,MAX_NUM), the point coordinates
!    of the dual.
!
!    Output, integer FACE_POINT2(MAX_ORDER,MAX_NUM), the vertices
!    of each face in the dual.
!
  implicit none
!
  integer max_num
  integer max_order
!
  integer face_num1
  integer face_num2
  integer face_order1
  integer face_order2
  integer face_point1(max_order,max_num)
  integer face_point2(max_order,max_num)
  integer i
  integer icol
  integer iface
  integer inext
  integer iprev
  integer irow
  integer istop
  integer j
  integer k
  real norm
  integer point_num1
  integer point_num2
  real point_coord1(3,max_num)
  real point_coord2(3,max_num)
  real x
  real y
  real z
!
  point_num2 = face_num1
  face_num2 = point_num1
!
!  This computation should really compute the center of gravity
!  of the face, in the general case.
!
!  We'll also assume the vertices of the original and the dual
!  are to lie on the unit sphere, so we can normalize the
!  position vector of the vertex.
!
  do i = 1, face_num1

    x = 0.0E+00
    y = 0.0E+00
    z = 0.0E+00
    do j = 1, face_order1
      k = face_point1(j,i)
      x = x + point_coord1(1,k)
      y = y + point_coord1(2,k)
      z = z + point_coord1(3,k)
    end do

    norm = sqrt ( x * x + y * y + z * z )

    point_coord2(1,i) = x / norm
    point_coord2(2,i) = y / norm
    point_coord2(3,i) = z / norm

  end do
!
!  Now build the face in the dual associated with each node IFACE.
!
  do iface = 1, face_num2
!
!  Initialize the order.
!
    face_order2 = 0
!
!  Find the first occurrence of IFACE in an edge of polyhedron 1.
!
    call icol_find_item ( face_point1, max_order, face_order1, face_num1, &
      iface, irow, icol )

    if ( irow == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DUAL_SHAPE_3D - Fatal error!'
      write ( *, '(a,i6)' ) '  Could not find an edge using node ', iface
      return
    end if
!
!  Save the following node as ISTOP.
!  When we encounter ISTOP again, this will mark the end of our search.
!
    i = irow + 1
    if ( i > face_order1 ) then
      i = 1
    end if

    istop = face_point1(i,icol)
!
!  Save the previous node as INEXT.
!
    do

      i = irow - 1
      if ( i < 1 ) then
        i = i + face_order1
      end if

      inext = face_point1(i,icol)

      face_order2 = face_order2 + 1

      face_point2(face_order2,iface) = icol
!
!  If INEXT =/= ISTOP, continue.
!
      if ( inext == istop ) then
        exit
      end if
!
!  Set IPREV:= INEXT.
!
      iprev = inext
!
!  Search for the occurrence of the edge IFACE-IPREV.
!
      call icol_find_pair_wrap ( face_point1, max_order, face_order1, &
        face_num1, iface, iprev, irow, icol )

      if ( irow == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DUAL_SHAPE_3D - Fatal error!'
        write ( *, '(a,i6)' ) '  No edge from node ', iprev
        write ( *, '(a,i6)' ) '  to node ', iface
        return
      end if

    end do

  end do

  return
end
subroutine ellipse_area_2d ( r1, r2, area )
!
!*******************************************************************************
!
!! ELLIPSE_AREA_2D returns the area of an ellipse in 2D.
!
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    Output, real AREA, the area of the ellipse.
!
  implicit none
!
  real area
  real r_pi
  real r1
  real r2
!
  area = r_pi ( ) * r1 * r2

  return
end
subroutine ellipse_points_2d ( x0, y0, r1, r2, psi, n, x, y )
!
!*******************************************************************************
!
!! ELLIPSE_POINTS_2D returns N points on an tilted ellipse in 2D.
!
!
!  Discussion:
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter of the ellipse.
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, the coordinates of the center of the ellipse.
!
!    Input, real R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    Input, real PSI, the angle that the major axis of the ellipse
!    makes with the X axis.  A value of 0.0 means that the major and
!    minor axes of the ellipse will be the X and Y coordinate axes.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real X(N), Y(N), the coordinates of points on the ellipse.
!
  implicit none
!
  integer n
!
  integer i
  real r_pi
  real psi
  real r1
  real r2
  real theta
  real x0
  real x(n)
  real y0
  real y(n)
!
  do i = 1, n

    theta = ( 2.0E+00 * r_pi ( ) * real ( i - 1 ) ) / real ( n )

    x(i) = x0 + r1 * cos ( psi ) * cos ( theta ) &
              - r2 * sin ( psi ) * sin ( theta )

    y(i) = y0 + r1 * sin ( psi ) * cos ( theta ) &
              + r2 * cos ( psi ) * sin ( theta )

  end do

  return
end
subroutine ellipse_points_arc_2d ( x0, y0, r1, r2, psi, theta1, theta2, n, &
  x, y )
!
!*******************************************************************************
!
!! ELLIPSE_POINTS_ARC_2D returns N points on a tilted elliptical arc in 2D.
!
!
!  Discussion:
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter of the ellipse.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, the coordinates of the center of the ellipse.
!
!    Input, real R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    Input, real PSI, the angle that the major axis of the ellipse
!    makes with the X axis.  A value of 0.0 means that the major and
!    minor axes of the ellipse will be the X and Y coordinate axes.
!
!    Input, real THETA1, THETA2, the angular coordinates of the first
!    and last points to be drawn, in radians.  This angle is measured
!    with respect to the (possibly tilted) major axis.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real X(N), Y(N), the coordinates of points on the ellipse.
!
  implicit none
!
  integer n
!
  integer i
  real r_pi
  real psi
  real r1
  real r2
  real r_modp
  real theta
  real theta1
  real theta2
  real theta3
  real x0
  real x(n)
  real y0
  real y(n)
!
!  THETA3 is the smallest angle, no less than THETA1, which
!  coincides with THETA2.
!
  theta3 = theta1 + r_modp ( theta2 - theta1, 2.0E+00 * r_pi ( ) )

  do i = 1, n

    if ( n > 1 ) then
      theta = ( real ( n - i ) * theta1 + real ( i - 1 ) * theta3 ) &
        / real ( n - 1 )
    else
      theta = 0.5E+00 * ( theta1 + theta3 )
    end if

    x(i) = x0 + r1 * cos ( psi ) * cos ( theta ) &
              - r2 * sin ( psi ) * sin ( theta )

    y(i) = y0 + r1 * sin ( psi ) * cos ( theta ) &
              + r2 * cos ( psi ) * sin ( theta )

  end do

  return
end
function enorm0_2d ( x0, y0, x1, y1 )
!
!*******************************************************************************
!
!! ENORM0_2D computes the Euclidean norm of (P1-P0) in 2D.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, X1, Y1, the coordinates of the points P0 and P1.
!
!    Output, real ENORM0_2D, the Euclidean norm of (P1-P0).
!
  implicit none
!
  real enorm0_2d
  real x0
  real x1
  real y0
  real y1
!
  enorm0_2d = sqrt ( ( x1 - x0 )**2 + ( y1 - y0 )**2 )

  return
end
function enorm0_3d ( x0, y0, z0, x1, y1, z1 )
!
!*******************************************************************************
!
!! ENORM0_3D computes the Euclidean norm of (P1-P0) in 3D.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, X1, Y1, Z1, the coordinates of the points
!    P0 and P1.
!
!    Output, real ENORM0_3D, the Euclidean norm of (P1-P0).
!
  implicit none
!
  real enorm0_3d
  real x0
  real x1
  real y0
  real y1
  real z0
  real z1
!
  enorm0_3d = sqrt ( ( x1 - x0 )**2 + ( y1 - y0 )**2 + ( z1 - z0 )**2 )

  return
end
function enorm0_nd ( n, x0, x1 )
!
!*******************************************************************************
!
!! ENORM0_ND computes the Euclidean norm of (P1-P0) in ND.
!
!
!  Modified:
!
!    18 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Input, real X0(N), the base vector.
!
!    Input, real X1(N), the displacement vector.
!
!    Output, real ENORM0_ND, the Euclidean norm of the vector
!    ( X1 - X0 ).
!
  integer n
!
  real enorm0_nd
  real x0(n)
  real x1(n)
!
  enorm0_nd = sqrt ( sum ( ( x1(1:n) - x0(1:n) )**2 ) )

  return
end
function enorm_2d ( x1, y1 )
!
!*******************************************************************************
!
!! ENORM_2D computes the Euclidean norm of a vector in 2D.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, the coordinates of the vector.
!
!    Output, real ENORM_2D, the Euclidean norm of the vector.
!
  implicit none
!
  real enorm_2d
  real x1
  real y1
!
  enorm_2d = sqrt ( x1 * x1 + y1 * y1 )

  return
end
function enorm_3d ( x1, y1, z1 )
!
!*******************************************************************************
!
!! ENORM_3D computes the Euclidean norm of a vector in 3D.
!
!
!  Modified:
!
!    27 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, the coordinates of the vector.
!
!    Output, real ENORM_3D, the Euclidean norm of the vector.
!
  implicit none
!
  real enorm_3d
  real x1
  real y1
  real z1
!
  enorm_3d = sqrt ( x1 * x1 + y1 * y1 + z1 * z1 )

  return
end
function enorm_nd ( n, x )
!
!*******************************************************************************
!
!! ENORM_ND computes the Euclidean norm of a vector in ND.
!
!
!  Modified:
!
!    27 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Input, real X(N), the coordinates of the vector.
!
!    Output, real ENORM_ND, the Euclidean norm of the vector.
!
  implicit none
!
  integer n
!
  real enorm_nd
  integer i
  real x(n)
!
  enorm_nd = sqrt ( sum ( x(1:n)**2 ) )

  return
end
function enormsq0_2d ( x0, y0, x1, y1 )
!
!*******************************************************************************
!
!! ENORMSQ0_2D computes the square of the Euclidean norm of (P1-P0) in 2D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, X1, Y1, the coordinates of the points
!    P0 and P1.
!
!    Output, real ENORMSQ0_2D, the square of the Euclidean norm of (P1-P0).
!
  implicit none
!
  real enormsq0_2d
  real x0
  real x1
  real y0
  real y1
!
  enormsq0_2d = ( x1 - x0 )**2 + ( y1 - y0 )**2

  return
end
function enormsq0_3d ( x0, y0, z0, x1, y1, z1 )
!
!*******************************************************************************
!
!! ENORMSQ0_3D computes the square of the Euclidean norm of (P1-P0) in 3D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, Z0, X1, Y1, Z1, the coordinates of the points
!    P0 and P1.
!
!    Output, real ENORMSQ0_3D, the square of the Euclidean norm of (P1-P0).
!
  implicit none
!
  real enormsq0_3d
  real x0
  real x1
  real y0
  real y1
  real z0
  real z1
!
  enormsq0_3d =( x1 - x0 )**2 + ( y1 - y0 )**2 + ( z1 - z0 )**2

  return
end
function enormsq0_nd ( n, x0, x1 )
!
!*******************************************************************************
!
!! ENORMSQ0_ND computes the squared Euclidean norm of (P1-P0) in ND.
!
!
!  Modified:
!
!    18 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Input, real X0(N), the coordinates of the base vector.
!
!    Input, real X1(N), the coordinates of the displacement vector.
!
!    Output, real ENORMSQ0_ND, the squared Euclidean norm of the vector
!    ( X1 - X0 ).
!
  implicit none
!
  integer n
!
  real enormsq0_nd
  real x0(n)
  real x1(n)
!
  enormsq0_nd = sum ( ( x1(1:n) - x0(1:n) )**2 )

  return
end
subroutine glob2loc_3d ( cospitch, cosroll, cosyaw, locpts, globas, glopts, &
  sinpitch, sinroll, sinyaw )
!
!*******************************************************************************
!
!! GLOB2LOC_3D converts from a global to a local coordinate system in 3D.
!
!
!  Discussion:
!
!    A global coordinate system is given.
!
!    A local coordinate system has been translated to the point with
!    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
!    a roll.
!
!    A point has global coordinates GLOPTS, and it is desired to know
!    the point's local coordinates LOCPTS.
!
!    The transformation may be written as
!
!      LOC = M_ROLL * M_PITCH * M_YAW * ( GLOB - GLOBAS )
!
!    where
!
!               (       1            0            0      )
!    M_ROLL =   (       0        cos(Roll)    sin(Roll)  )
!               (       0      - sin(Roll)    cos(Roll)  )
!
!               (   cos(Pitch)       0      - sin(Pitch) )
!    M_PITCH =  (       0            1            0      )
!               (   sin(Pitch)       0        cos(Pitch) )
!
!               (   cos(Yaw)     sin(Yaw)         0      )
!    M_YAW    = ( - sin(Yaw)     cos(Yaw)         0      )
!               (       0            0            1      )
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COSPITCH, COSROLL, COSYAW, the cosines of the pitch,
!    roll and yaw angles.
!
!    Input, real GLOBAS(3), the X, Y, and Z coordinates of the
!    global base vector.
!
!    Input, real GLOPTS(3), the global ( X, Y and Z ) coordinates of
!    the point whose coordinates are to be transformed.
!
!    Output, real LOCPTS(3), the local coordinates of the point
!    whose global coordinates were given in GLOPTS.
!
!    Input, real SINPITCH, SINROLL, SINYAW, the sines of the pitch,
!    roll and yaw angles.
!
  implicit none
!
  real cospitch
  real cosroll
  real cosyaw
  real globas(3)
  real glopts(3)
  real locpts(3)
  real sinpitch
  real sinroll
  real sinyaw
!
  locpts(1) = ( cosyaw * cospitch ) * ( glopts(1) - globas(1) ) &
            + ( sinyaw * cospitch ) * ( glopts(2) - globas(2) ) &
            -   sinpitch * ( glopts(3) - globas(3) )

  locpts(2) = ( cosyaw * sinpitch * sinroll - sinyaw * cosroll ) &
    * ( glopts(1) - globas(1) ) &
    + ( sinyaw * sinpitch * sinroll + cosyaw * cosroll ) &
    * ( glopts(2) - globas(2) ) &
    +   cospitch * sinroll * ( glopts(3) - globas(3) )

  locpts(3) = ( cosyaw * sinpitch * cosroll + sinyaw * sinroll ) &
    * ( glopts(1) - globas(1) ) &
    + ( sinyaw * sinpitch * cosroll - cosyaw * sinroll  ) &
    * ( glopts(2) - globas(2) ) &
    + ( cospitch * cosroll ) * ( glopts(3) - globas(3) )

  return
end
function halfplane_contains_point_2d ( xa, ya, xb, yb, x, y )
!
!*******************************************************************************
!
!! HALFPLANE_CONTAINS_POINT_2D reports if a half-plane contains a point in 2d.
!
!
!  Discussion:
!
!    The halfplane is assumed to be all the points "to the left" of the
!    line segment from PA = (XA,YA) to PB = (XB,YB).  Thus, one way to
!    understand where the point P = (X,Y) is, is to compute the signed
!    area of the triangle ( PA, PB, P ).
!
!    If this area is
!      positive, the point is strictly inside the halfplane,
!      zero, the point is on the boundary of the halfplane,
!      negative, the point is strictly outside the halfplane.
!
!  Modified:
!
!    09 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XA, YA, XB, YB, two distinct points (XA,YA) and (XB,YB)
!    that lie on the line defining the half plane.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, logical HALFPLANE_CONTAINS_POINT_2D, is TRUE if the halfplane
!    contains the point, and FALSE otherwise.
!
  implicit none
!
  real area_signed
  logical halfplane_contains_point_2d
  real x
  real xa
  real xb
  real y
  real ya
  real yb
!
  area_signed = 0.5E+00 * &
    ( xa * ( yb - y ) + xb * ( y - ya ) + x * ( ya - yb ) )

  halfplane_contains_point_2d = ( area_signed >= 0.0E+00 )

  return
end
subroutine halfspace_imp_triangle_int_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  a, b, c, d, num_int, x, y, z )
!
!*******************************************************************************
!
!! HALFSPACE_IMP_TRIANGLE_INT_3D: intersection ( implicit halfspace, triangle ) in 3D.
!
!
!  Definition:
!
!    The implicit form of a half-space in 3D may be described as the set
!    of points (X,Y,Z) on or "above" an implicit plane:
!
!      A * X + B * Y + C * Z + D >= 0
!
!    The triangle is specified by listing its three vertices.
!
!  Discussion:
!
!    The intersection may be described by the number of vertices of the
!    triangle that are included in the halfspace, and by the location of
!    points between vertices that separate a side of the triangle into
!    an included part and an unincluded part.
!
!    0 vertices, 0 separators    (no intersection)
!    1 vertex, 0 separators      (point intersection)
!    2 vertices, 0 separators    (line intersection)
!    3 vertices, 0 separators    (triangle intersection)
!
!    1 vertex, 2 separators,     (intersection is a triangle)
!    2 vertices, 2 separators,   (intersection is a quadrilateral).
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of the vertices of the triangle.
!
!    Input, real A, B, C, D, the parameters that define the implicit plane,
!    which in turn define the implicit halfspace.
!
!    Output, integer NUM_INT, the number of intersection points returned,
!    which will always be between 0 and 4.
!
!    Output, real X(4), Y(4), Z(4), the coordinates of the NUM_INT
!    intersection points.  The points will lie in sequence on the triangle.
!    Some points will be vertices, and some may be separators.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real dist1
  real dist2
  real dist3
  integer num_int
  real x(4)
  real x1
  real x2
  real x3
  real y(4)
  real y1
  real y2
  real y3
  real z(4)
  real z1
  real z2
  real z3
!
!  Compute the signed distances between the vertices and the plane.
!
  dist1 = a * x1 + b * y1 + c * z1 + d
  dist2 = a * x2 + b * y2 + c * z2 + d
  dist3 = a * x3 + b * y3 + c * z3 + d
!
!  Now we can find the intersections.
!
  call halfspace_triangle_int_3d ( x1, y1, z1, x2, y2, &
    z2, x3, y3, z3, dist1, dist2, dist3, num_int, x, y, z )

  return
end
subroutine halfspace_norm_triangle_int_3d ( x1, y1, z1, x2, y2, z2, &
  x3, y3, z3, xp, yp, zp, xn, yn, zn, num_int, x, y, z )
!
!*******************************************************************************
!
!! HALFSPACE_NORM_TRIANGLE_INT_3D: intersection ( normal halfspace, triangle ) in 3D.
!
!
!  Definition:
!
!    The normal form of a halfspace in 3D may be described as the set
!    of points (X,Y,Z) on or "above" a plane described in normal form:
!
!      ( Xp, Yp, Zp ) is a point on the plane,
!      ( Xn, Yn, Zn ) is the unit normal vector, pointing "out" of the
!      halfspace.
!
!    The triangle is specified by listing its three vertices.
!
!  Discussion:
!
!    The intersection may be described by the number of vertices of the
!    triangle that are included in the halfspace, and by the location of
!    points between vertices that separate a side of the triangle into
!    an included part and an unincluded part.
!
!    0 vertices, 0 separators    (no intersection)
!    1 vertex, 0 separators      (point intersection)
!    2 vertices, 0 separators    (line intersection)
!    3 vertices, 0 separators    (triangle intersection)
!
!    1 vertex, 2 separators,     (intersection is a triangle)
!    2 vertices, 2 separators,   (intersection is a quadrilateral).
!
!  Modified:
!
!    03 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of the vertices of the triangle.
!
!    Input, real XP, YP, ZP, a point on the bounding plane that defines
!    the halfspace.
!
!    Input, real XN, YN, ZN, the components of the normal vector to the
!    bounding plane that defines the halfspace.  By convention, the
!    normal vector points "outwards" from the halfspace.
!
!    Output, integer NUM_INT, the number of intersection points returned,
!    which will always be between 0 and 4.
!
!    Output, real X(4), Y(4), Z(4), the coordinates of the NUM_INT
!    intersection points.  The points will lie in sequence on the triangle.
!    Some points will be vertices, and some may be separators.
!
  implicit none
!
  real d
  real dist1
  real dist2
  real dist3
  integer num_int
  real x(4)
  real x1
  real x2
  real x3
  real xn
  real xp
  real y(4)
  real y1
  real y2
  real y3
  real yn
  real yp
  real z(4)
  real z1
  real z2
  real z3
  real zn
  real zp
!
!  Compute the signed distances between the vertices and the plane.
!
  d = - xn * xp - yn * yp - zn * zp

  dist1 = xn * x1 + yn * y1 + zn * z1 + d
  dist2 = xn * x2 + yn * y2 + zn * z2 + d
  dist3 = xn * x3 + yn * y3 + zn * z3 + d
!
!  Now we can find the intersections.
!
  call halfspace_triangle_int_3d ( x1, y1, z1, x2, y2, &
    z2, x3, y3, z3, dist1, dist2, dist3, num_int, x, y, z )

  return
end
subroutine halfspace_triangle_int_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  dist1, dist2, dist3, num_int, x, y, z )
!
!*******************************************************************************
!
!! HALFSPACE_TRIANGLE_INT_3D: intersection ( halfspace, triangle ) in 3D.
!
!
!  Definition:
!
!    The triangle is specified by listing its three vertices.
!
!  Discussion:
!
!    The halfspace is not described in the input data.  Rather, the
!    distances from the triangle vertices to the halfspace are given.
!
!    The intersection may be described by the number of vertices of the
!    triangle that are included in the halfspace, and by the location of
!    points between vertices that separate a side of the triangle into
!    an included part and an unincluded part.
!
!    0 vertices, 0 separators    (no intersection)
!    1 vertex, 0 separators      (point intersection)
!    2 vertices, 0 separators    (line intersection)
!    3 vertices, 0 separators    (triangle intersection)
!
!    1 vertex, 2 separators,     (intersection is a triangle)
!    2 vertices, 2 separators,   (intersection is a quadrilateral).
!
!  Modified:
!
!    03 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of the vertices of the triangle.
!
!    Input, real DIST1, DIST2, DIST3, the distances from each of the
!    three vertices of the triangle to the halfspace.  The distance is
!    zero if a vertex lies within the halfspace, or on the plane that
!    defines the boundary of the halfspace.  Otherwise, it is the
!    distance from that vertex to the bounding plane.
!
!    Output, integer NUM_INT, the number of intersection points returned,
!    which will always be between 0 and 4.
!
!    Output, real X(4), Y(4), Z(4), the coordinates of the NUM_INT
!    intersection points.  The points will lie in sequence on the triangle.
!    Some points will be vertices, and some may be separators.
!
  implicit none
!
  real dist1
  real dist2
  real dist3
  integer num_int
  real x(4)
  real x1
  real x2
  real x3
  real y(4)
  real y1
  real y2
  real y3
  real z(4)
  real z1
  real z2
  real z3
!
!  Walk around the triangle, looking for vertices that are included,
!  and points of separation.
!
  num_int = 0

  if ( dist1 <= 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x1
    y(num_int) = y1
    z(num_int) = z1

  end if

  if ( dist1 * dist2 < 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = ( dist1 * x2 - dist2 * x1 ) / ( dist1 - dist2 )
    y(num_int) = ( dist1 * y2 - dist2 * y1 ) / ( dist1 - dist2 )
    z(num_int) = ( dist1 * z2 - dist2 * z1 ) / ( dist1 - dist2 )

  end if

  if ( dist2 <= 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x2
    y(num_int) = y2
    z(num_int) = z2

  end if

  if ( dist2 * dist3 < 0.0E+00 ) then

    num_int = num_int + 1

    x(num_int) = ( dist2 * x3 - dist3 * x2 ) / ( dist2 - dist3 )
    y(num_int) = ( dist2 * y3 - dist3 * y2 ) / ( dist2 - dist3 )
    z(num_int) = ( dist2 * z3 - dist3 * z2 ) / ( dist2 - dist3 )

  end if

  if ( dist3 <= 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x3
    y(num_int) = y3
    z(num_int) = z3

  end if

  if ( dist3 * dist1 < 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = ( dist3 * x1 - dist1 * x3 ) / ( dist3 - dist1 )
    y(num_int) = ( dist3 * y1 - dist1 * y3 ) / ( dist3 - dist1 )
    z(num_int) = ( dist3 * z1 - dist1 * z3 ) / ( dist3 - dist1 )

  end if

  return
end
function haversine ( a )
!
!*******************************************************************************
!
!! HAVERSINE computes the haversine of an angle.
!
!
!  Discussion:
!
!    haversine(A) = ( 1 - cos ( A ) ) / 2
!
!    The haversine is useful in spherical trigonometry.
!
!  Modified:
!
!    02 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the angle.
!
!    Output, real HAVERSINE, the haversine of the angle.
!
  implicit none
!
  real a
  real haversine
!
  haversine = ( 1.0E+00 - cos ( a ) ) / 2.0E+00

  return
end
subroutine helix_shape_3d ( a, n, r, theta1, theta2, x, y, z )
!
!*******************************************************************************
!
!! HELIX_SHAPE_3D computes points on a helix in 3D.
!
!
!  Discussion:
!
!    The user specifies the parameters A and R, the first and last
!    THETA values, and the number of equally spaced THETA values
!    at which (X,Y,Z) values are to be computed.
!
!  Formula:
!
!    X = R * COS ( THETA )
!    Y = R * SIN ( THETA )
!    Z = A * THETA
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the rate at which Z advances with THETA.
!
!    Input, integer N, the number of points to compute on the helix.
!
!    Input, real R, the radius of the helix.
!
!    Input, real THETA1, THETA2, the first and last THETA values at
!    which to compute points on the helix.  THETA is measured in
!    radians.
!
!    Output, real X(N), Y(N), Z(N), the X, Y and Z coordinates
!    of points on the helix.
!
  implicit none
!
  integer n
!
  real a
  integer i
  real r
  real theta
  real theta1
  real theta2
  real x(n)
  real y(n)
  real z(n)
!
  do i = 1, n

    if ( n == 1 ) then
      theta = 0.5E+00 * ( theta1 + theta2 )
    else
      theta = ( real ( n - i ) * theta1 + real ( i - 1 ) * theta2 ) &
        / real ( n - 1 )
    end if

    x(i) = r * cos ( theta )
    y(i) = r * sin ( theta )
    z(i) = a * theta

  end do

  return
end
function hexagon_area_2d ( rad )
!
!*******************************************************************************
!
!! HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
!
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RAD, the radius of the hexagon.
!
!    Output, real HEXAGON_AREA_2D, the area of the hexagon.
!
  implicit none
!
  real hexagon_area_2d
  real hexagon_unit_area_2d
  real rad
!
  hexagon_area_2d = rad**2 * hexagon_unit_area_2d ( )

  return
end
subroutine hexagon_shape_2d ( angle, x, y )
!
!*******************************************************************************
!
!! HEXAGON_SHAPE_2D returns points on the unit hexagon in 2D.
!
!
!  Diagram:
!
!      120_____60
!        /     \
!    180/       \0
!       \       /
!        \_____/
!      240     300
!
!  Discussion:
!
!    The unit hexagon has maximum radius 1, and is the hull of the points
!
!      (   1,              0 ),
!      (   0.5,   sqrt (3)/2 ),
!      ( - 0.5,   sqrt (3)/2 ),
!      ( - 1,              0 ),
!      ( - 0.5, - sqrt (3)/2 ),
!      (   0.5, - sqrt (3)/2 ).
!
!  Modified:
!
!    12 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle, in degrees, of the point.
!
!    Output, real X, Y, the coordinates of the point.
!
  implicit none
!
  real angle
  real angle2
  real cot_degrees
  real tan_degrees
  real x
  real y
!
  angle2 = mod ( angle, 360.0E+00 )
  if ( angle2 < 0.0E+00 ) then
    angle2 = angle2 + 360.0E+00
  end if
!
!  y = - sqrt(3) * x + sqrt(3)
!
  if ( 0.0E+00 <= angle2 .and. angle2 <= 60.0E+00 ) then

    x = sqrt ( 3.0E+00 ) / ( tan_degrees ( angle2 ) + sqrt ( 3.0E+00 ) )
    y = tan_degrees ( angle2 ) * x
!
!  y = sqrt(3) / 2
!
  else if ( angle2 <= 120.0E+00 ) then

    y = sqrt ( 3.0E+00 ) / 2.0E+00
    x = cot_degrees ( angle2 ) * y
!
!  y = sqrt(3) * x + sqrt(3)
!
  else if ( angle2 <= 180.0E+00 ) then

    x = sqrt ( 3.0E+00 ) / ( tan_degrees ( angle2 ) - sqrt ( 3.0E+00 ) )
    y = tan_degrees ( angle2 ) * x
!
!  y = - sqrt(3) * x - sqrt(3)
!
  else if ( angle2 <= 240.0E+00 ) then

    x = - sqrt ( 3.0E+00 ) / ( tan_degrees ( angle2 ) + sqrt ( 3.0E+00 ) )
    y = tan_degrees ( angle2 ) * x
!
!  y = - sqrt(3) / 2
!
  else if ( angle2 <= 300.0E+00 ) then

    y = - sqrt ( 3.0E+00 ) / 2.0E+00
    x = y * cot_degrees ( angle2 )
!
!  y = sqrt(3) * x - sqrt(3)
!
  else if ( angle2 <= 360.0E+00 ) then

    x = - sqrt ( 3.0E+00 ) / ( tan_degrees ( angle2 ) - sqrt ( 3.0E+00 ) )
    y = tan_degrees ( angle2 ) * x

  end if

  return
end
function hexagon_unit_area_2d ( )
!
!*******************************************************************************
!
!! HEXAGON_UNIT_AREA_2D returns the area of a unit regular hexagon in 2D.
!
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real HEXAGON_UNIT_AREA_2D, the area of the hexagon.
!
  implicit none
!
  real hexagon_unit_area_2d
  real rad
!
  hexagon_unit_area_2d = 3.0D+00 * sqrt ( 3.0D+00 ) / 2.0D+00

  return
end
subroutine hexagon_vertices_2d ( x, y )
!
!*******************************************************************************
!
!! HEXAGON_VERTICES_2D returns the vertices of the unit hexagon in 2D.
!
!
!  Diagram:
!
!      120_____60
!        /     \
!    180/       \0
!       \       /
!        \_____/
!      240     300
!
!  Discussion:
!
!    The unit hexagon has maximum radius 1, and is the hull of the points
!
!      (   1,              0 ),
!      (   0.5,   sqrt (3)/2 ),
!      ( - 0.5,   sqrt (3)/2 ),
!      ( - 1,              0 ),
!      ( - 0.5, - sqrt (3)/2 ),
!      (   0.5, - sqrt (3)/2 ).
!
!  Modified:
!
!    21 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real X(6), Y(6), the coordinates of the vertices.
!
  implicit none
!
  real, parameter :: a = 0.8660254037844386E+00
  real x(6)
  real y(6)
!
  x(1:6) = (/ 1.0E+00, 0.5E+00, -0.5E+00, -1.0E+00, -0.5E+00,  0.5E+00 /)
  y(1:6) = (/ 0.0E+00, a,        a,        0.0E+00, -a,       -a /)

  return
end
function i_modp ( i, j )
!
!*******************************************************************************
!
!! I_MODP returns the nonnegative remainder of integer division.
!
!
!  Formula:
!
!    If
!      NREM = I_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I_MODP(A,360) is between 0 and 360, always.
!
!  Examples:
!
!        I     J     MOD  I_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none
!
  integer i
  integer i_modp
  integer j
!
  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I_MODP ( I, J ) called with J = ', j
    stop
  end if

  i_modp = mod ( i, j )

  if ( i_modp < 0 ) then
    i_modp = i_modp + abs ( j )
  end if

  return
end
subroutine i_random ( ilo, ihi, i )
!
!*******************************************************************************
!
!! I_RANDOM returns a random integer in a given range.
!
!
!  Modified:
!
!    23 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ILO, IHI, the minimum and maximum acceptable values.
!
!    Output, integer I, the randomly chosen integer.
!
  implicit none
!
  logical, save :: seed = .false.
  integer i
  integer ihi
  integer ilo
  real r
  real rhi
  real rlo
!
  if ( .not. seed ) then
    call random_seed ( )
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = r )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo ) - 0.5E+00
  rhi = real ( ihi ) + 0.5E+00
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0E+00 - r ) * rlo + r * rhi )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
subroutine i_swap ( i, j )
!
!*******************************************************************************
!
!! I_SWAP switches two integer values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none
!
  integer i
  integer j
  integer k
!
  k = i
  i = j
  j = k

  return
end
function i_wrap ( ival, ilo, ihi )
!
!*******************************************************************************
!
!! I_WRAP forces an integer to lie between given limits by wrapping.
!
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, an integer value.
!
!    Input, integer ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer I_WRAP, a "wrapped" version of IVAL.
!
  implicit none
!
  integer i_modp
  integer i_wrap
  integer ihi
  integer ilo
  integer ival
  integer wide
!
  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i_wrap = ilo
  else
    i_wrap = ilo + i_modp ( ival-ilo, wide )
  end if

  return
end
subroutine icol_find_item ( itable, maxrow, nrow, ncol, item, irow, icol )
!
!*******************************************************************************
!
!! ICOL_FIND_ITEM searches a table by columns for a given value.
!
!
!  Modified:
!
!    11 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ITABLE(MAXROW,NCOL), the table to search.
!
!    Input, integer MAXROW, the leading dimension of ITABLE.
!
!    Input, integer NROW, NCOL, the number of rows and columns
!    in the table.
!
!    Input, integer ITEM, the value to search for.
!
!    Output, integer IROW, ICOL, the row and column indices
!    of the first occurrence of the value ITEM.  The search
!    is conducted by rows.  If the item is not found, then
!    IROW = ICOL = 0.
!
  implicit none
!
  integer maxrow
  integer ncol
!
  integer i
  integer icol
  integer irow
  integer itable(maxrow,ncol)
  integer item
  integer nrow
  integer j
!
  do j = 1, ncol
    do i = 1, nrow
      if ( itable(i,j) == item ) then
        irow = i
        icol = j
        return
      end if
    end do
  end do

  irow = 0
  icol = 0

  return
end
subroutine icol_find_pair_wrap ( itable, maxrow, nrow, ncol, item1, item2, &
  irow, icol )
!
!*******************************************************************************
!
!! ICOL_FIND_ITEM wrap searches a table by columns for a pair of items.
!
!
!  Discussion:
!
!    The items must occur consecutively, with ITEM1 occurring
!    first.  However, wrapping is allowed.  That is, if ITEM1
!    occurs in the last row, and ITEM2 in the first, this
!    is also regarded as a match.
!
!  Modified:
!
!    11 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ITABLE(MAXROW,NCOL), the table to search.
!
!    Input, integer MAXROW, the leading dimension of ITABLE.
!
!    Input, integer NROW, NCOL, the number of rows and columns
!    in the table.
!
!    Input, integer ITEM1, ITEM2, the values to search for.
!
!    Output, integer IROW, ICOL, the row and column indices
!    of the first occurrence of the value ITEM1 followed immediately
!    by ITEM2.  The search is conducted by columns.  If the pair of
!    items is not found, then IROW = ICOL = 0.  If IROW = NROW,
!    the ITEM1 occurs in row NROW and ITEM2 occurs in row 1.
!
  implicit none
!
  integer maxrow
  integer ncol
!
  integer i
  integer icol
  integer ip1
  integer irow
  integer itable(maxrow,ncol)
  integer item1
  integer item2
  integer j
  integer nrow
!
  do j = 1, ncol
    do i = 1, nrow

      if ( itable(i,j) == item1 ) then

        if ( i < nrow ) then
          ip1 = i + 1
        else
          ip1 = 1
        end if

        if ( itable(ip1,j) == item2 ) then
          irow = i
          icol = j
          return
        end if

      end if

    end do
  end do

  irow = 0
  icol = 0

  return
end
subroutine icos_shape_3d ( max_num, max_order, point_num, face_num, &
  face_order, point_coord, face_point )
!
!*******************************************************************************
!
!! ICOS_SHAPE_3D describes an icosahedron in 3D.
!
!
!  Discussion:
!
!    The vertices of the icosahedron lie on the unit sphere.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAX_NUM, the maximum number of faces, and of
!    points.  This quantity is used to dimension arrays.
!    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
!    not be set.
!
!    Input, integer MAX_ORDER, the maximum order of any face.
!    This quantity is used to dimension arrays.  If FACE_ORDER
!    exceeds MAX_ORDER, the arrays will not be set.
!
!    Output, integer POINT_NUM, the number of points in the shape.
!
!    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
!    the X, Y and Z coordinates of point J.
!
!    Output, integer FACE_NUM, the number of faces in the shape.
!
!    Output, integer FACE_ORDER, the number of vertices per face.
!
!    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
!    contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter-clockwise direction defined
!    by the outward normal at the face.
!
  implicit none
!
  integer max_num
  integer max_order
!
  real a
  real b
  integer face_num
  integer face_order
  integer face_point(max_order,max_num)
  real phi
  integer point_num
  real point_coord(3,max_num)
!
  point_num = 12
  face_num = 20
  face_order = 3
!
!  Check.
!
  if ( point_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ICOS_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of vertices exceeds MAX_NUM.'
    return
  end if

  if ( face_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ICOS_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of faces exceeds MAX_NUM.'
    return
  end if

  if ( face_order > max_order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ICOS_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Face order exceeds MAX_ORDER.'
    return
  end if
!
!  Set point coordinates.
!
  phi = 0.5E+00 * ( sqrt ( 5.0E+00 ) + 1.0E+00 )
  b = 1.0E+00 / sqrt ( 1.0E+00 + phi * phi )
  a = phi * b

  point_coord(1,1) =   a
  point_coord(2,1) =   b
  point_coord(3,1) = 0.0E+00

  point_coord(1,2) = - a
  point_coord(2,2) =   b
  point_coord(3,2) = 0.0E+00

  point_coord(1,3) =   a
  point_coord(2,3) = - b
  point_coord(3,3) = 0.0E+00

  point_coord(1,4) = - a
  point_coord(2,4) = - b
  point_coord(3,4) = 0.0E+00

  point_coord(1,5) =   b
  point_coord(2,5) = 0.0E+00
  point_coord(3,5) =   a

  point_coord(1,6) =   b
  point_coord(2,6) = 0.0E+00
  point_coord(3,6) = - a

  point_coord(1,7) = - b
  point_coord(2,7) = 0.0E+00
  point_coord(3,7) =   a

  point_coord(1,8) = - b
  point_coord(2,8) = 0.0E+00
  point_coord(3,8) = - a

  point_coord(1,9) =  0.0E+00
  point_coord(2,9) =    a
  point_coord(3,9) =    b

  point_coord(1,10) = 0.0E+00
  point_coord(2,10) = - a
  point_coord(3,10) =   b

  point_coord(1,11) = 0.0E+00
  point_coord(2,11) =   a
  point_coord(3,11) = - b

  point_coord(1,12) = 0.0E+00
  point_coord(2,12) = - a
  point_coord(3,12) = - b
!
!  Set faces.
!
  face_point(1,1) = 1
  face_point(2,1) = 9
  face_point(3,1) = 5

  face_point(1,2) = 1
  face_point(2,2) = 6
  face_point(3,2) = 11

  face_point(1,3) = 3
  face_point(2,3) = 5
  face_point(3,3) = 10

  face_point(1,4) = 3
  face_point(2,4) = 12
  face_point(3,4) = 6

  face_point(1,5) = 2
  face_point(2,5) = 7
  face_point(3,5) = 9

  face_point(1,6) = 2
  face_point(2,6) = 11
  face_point(3,6) = 8

  face_point(1,7) = 4
  face_point(2,7) = 10
  face_point(3,7) = 7

  face_point(1,8) = 4
  face_point(2,8) = 8
  face_point(3,8) = 12

  face_point(1,9) = 1
  face_point(2,9) = 11
  face_point(3,9) = 9

  face_point(1,10) = 2
  face_point(2,10) = 9
  face_point(3,10) = 11

  face_point(1,11) = 3
  face_point(2,11) = 10
  face_point(3,11) = 12

  face_point(1,12) = 4
  face_point(2,12) = 12
  face_point(3,12) = 10

  face_point(1,13) = 5
  face_point(2,13) = 3
  face_point(3,13) = 1

  face_point(1,14) = 6
  face_point(2,14) = 1
  face_point(3,14) = 3

  face_point(1,15) = 7
  face_point(2,15) = 2
  face_point(3,15) = 4

  face_point(1,16) = 8
  face_point(2,16) = 4
  face_point(3,16) = 2

  face_point(1,17) = 9
  face_point(2,17) = 7
  face_point(3,17) = 5

  face_point(1,18) = 10
  face_point(2,18) = 5
  face_point(3,18) = 7

  face_point(1,19) = 11
  face_point(2,19) = 6
  face_point(3,19) = 8

  face_point(1,20) = 12
  face_point(2,20) = 8
  face_point(3,20) = 6

  return
end
subroutine imat_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! IMAT_PRINT prints an integer matrix.
!
!
!  Modified:
!
!    31 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, integer A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer, parameter :: inc = 10
  integer j
  integer jhi
  integer jlo
  integer m
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, inc
    jhi = min ( jlo + inc - 1, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,10(i7))' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )
!
!*******************************************************************************
!
!! INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
!
!
!  Discussion:
!
!    The box has center at (IC,JC), and has half-widths N1 and N2.
!    The indices are exactly those which are between (IC-N1,JC-N2) and
!    (IC+N1,JC+N2) with the property that at least one of I and J
!    is an "extreme" value.
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the half-widths of the box, that is, the
!    maximum distance allowed between (IC,JC) and (I,J).
!
!    Input, integer IC, JC, the central cell of the box.
!
!    Input/output, integer I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none
!
  integer i
  integer ic
  integer j
  integer jc
  logical more
  integer n1
  integer n2
!
  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 ) then
    more = .false.
    return
  end if
!
!  Increment J.
!
  j = j + 1
!
!  Check J.
!
  if ( j > jc + n2 ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. ( i == ic - n1 .or. i == ic + n1 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, more )
!
!*******************************************************************************
!
!! INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
!
!
!  Discussion:
!
!    The box has a central cell of (IC,JC,KC), with a half widths of
!    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3) and
!    (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J, and K
!    is an "extreme" value.
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, the "half widths" of the box, that is, the
!    maximum distances from the central cell allowed for I, J and K.
!
!    Input, integer IC, JC, KC, the central cell of the box.
!
!    Input/output, integer I, J, K.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none
!
  integer i
  integer ic
  integer j
  integer jc
  integer k
  integer kc
  logical more
  integer n1
  integer n2
  integer n3
!
  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    k = kc - n3
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 .and. k == kc + n3 ) then
    more = .false.
    return
  end if
!
!  Increment K.
!
  k = k + 1
!
!  Check K.
!
  if ( k > kc + n3 ) then
    k = kc - n3
    j = j + 1
  else if ( k < kc + n3 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. j == jc - n2 .or. j == jc + n2 ) ) then
    return
  else
    k = kc + n3
    return
  end if
!
!  Check J.
!
  if ( j > jc + n2 ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. k == kc - n3 .or. k == kc + n3 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine irow_compare ( lda, m, n, a, i, j, isgn )
!
!*******************************************************************************
!
!! IROW_COMPARE compares two rows of a integer array.
!
!
!  Example:
!
!    Input:
!
!  M = 3, N = 4, I = 2, J = 3
!
!  A = (
!    1  2  3  4
!    5  6  7  8
!    9 10 11 12 )
!
!    Output:
!
!  ISGN = -1
!
!  Modified:
!
!    27 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array, which must
!    be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(LDA,N), an array of M rows of vectors of length N.
!
!    Input, integer I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row I > row J.
!
  implicit none
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer isgn
  integer j
  integer k
  integer m
!
!  Check that I and J are legal.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is less than 1.'
    write ( *, '(a,i6)' ) '  I = ', i
    stop
  else if ( i > m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i6)' ) '  I = ', i
    write ( *, '(a,i6)' ) '  Maximum legal value is M = ', m
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is less than 1.'
    write ( *, '(a,i6)' ) '  J = ', j
    stop
  else if ( j > m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i6)' ) '  J = ', j 
    write ( *, '(a,i6)' ) '  Maximum legal value is M = ', m
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = - 1
      return
    else if ( a(i,k) > a(j,k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine irow_sort_a ( lda, m, n, a )
!
!*******************************************************************************
!
!! IROW_SORT_A ascending sorts the rows of an integer array.
!
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, X is less than Y if, at the first index where they
!    differ, the X value is less than the Y value.
!
!  Example:
!
!    Input:
!
!      LDA = 5, M = 5, N = 3
!
!      A =
!        3  2  1
!        2  4  3
!        3  1  8
!        2  4  2
!        1  9  9
!
!    Output:
!
!      A =
!        1  9  9
!        2  4  2
!        2  4  3
!        3  1  8
!        3  2  1
!
!  Modified:
!
!    16 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, the number of rows of A.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(LDA,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in ascending
!    lexicographic order.
!
  implicit none
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer i
  integer indx
  integer isgn
  integer j
  integer m
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call irow_swap ( lda, m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call irow_compare ( lda, m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine irow_swap ( lda, m, n, a, irow1, irow2 )
!
!*******************************************************************************
!
!! IROW_SWAP swaps two rows of an integer array.
!
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, integer A(LDA,N), an array of data.
!
!    Input, integer IROW1, IROW2, the two rows to swap.
!
  implicit none
!
  integer lda
  integer n
!
  integer a(lda,n)
  integer m
  integer irow1
  integer irow2
  integer row(n)
!
!  Check.
!
  if ( irow1 < 1 .or. irow1 > m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  IROW1 is out of range.'
    stop
  end if

  if ( irow2 < 1 .or. irow2 > m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  IROW2 is out of range.'
    stop
  end if

  if ( irow1 == irow2 ) then
    return
  end if

  row(1:n) = a(irow1,1:n)
  a(irow1,1:n) = a(irow2,1:n)
  a(irow2,1:n) = row(1:n)

  return
end
subroutine ivec_heap_d ( n, a )
!
!*******************************************************************************
!
!! IVEC_HEAP_D reorders an array of integers into an descending heap.
!
!
!  Definition:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the input array.
!
!    Input/output, integer A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
  integer ifree
  integer key
  integer m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( m > n ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m+1) > a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine ivec_print ( n, a, title )
!
!*******************************************************************************
!
!! IVEC_PRINT prints an integer vector.
!
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer big
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(i6,1x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(i6,1x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i6,i11)' ) i, a(i)
    end do
  end if

  return
end
subroutine ivec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! IVEC_RANDOM returns a random integer vector in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ALO, AHI, the range allowed for the entries.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, integer A(N), the vector of randomly chosen integers.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer ahi
  integer alo
  integer i
!
  do i = 1, n

    call i_random ( alo, ahi, a(i) )

  end do

  return
end
subroutine ivec_sort_heap_a ( n, a )
!
!*******************************************************************************
!
!! IVEC_SORT_HEAP_A ascending sorts an integer array using heap sort.
!
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer n1
!
  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call ivec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call ivec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i_swap ( a(1), a(n1) )

  end do

  return
end
subroutine ivec_unique ( n, a, nuniq )
!
!*******************************************************************************
!
!! IVEC_UNIQUE finds the number of unique elements in a sorted integer array.
!
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in A.
!
!    Input/output, integer A(N).  On input, the sorted
!    integer array.  On output, the unique elements in A.
!
!    Output, integer NUNIQ, the number of unique elements in A.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer itest
  integer nuniq
!
  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a(itest) /= a(nuniq) ) then
      nuniq = nuniq + 1
      a(nuniq) = a(itest)
    end if

  end do

  return
end
subroutine ivec2_compare ( n, a1, a2, i, j, isgn )
!
!*******************************************************************************
!
!! IVEC2_COMPARE compares pairs of integers stored in two vectors.
!
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data items.
!
!    Input, integer A1(N), A2(N), contain the two components of each item.
!
!    Input, integer I, J, the items to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  implicit none
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer i
  integer isgn
  integer j
!
  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(i) > a2(j) ) then
      isgn = +1
    end if

  else if ( a1(i) > a1(j) ) then

    isgn = +1

  end if

  return
end
subroutine ivec2_sort_a ( n, a1, a2 )
!
!*******************************************************************************
!
!! IVEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, integer A1(N), A2(N), the data to be sorted..
!
  implicit none
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer i
  integer indx
  integer isgn
  integer j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call i_swap ( a1(i), a1(j) )
      call i_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call ivec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine ivec2_unique ( n, a1, a2, nuniq )
!
!*******************************************************************************
!
!! IVEC2_UNIQUE keeps the unique elements in a array of pairs of integers.
!
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items.
!
!    Input/output, integer A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of NUNIQ unique items.
!
!    Output, integer NUNIQ, the number of unique items.
!
  implicit none
!
  integer n
!
  integer a1(n)
  integer a2(n)
  integer itest
  integer nuniq
!
  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a1(itest) /= a1(nuniq) .or. a2(itest) /= a2(nuniq) ) then

      nuniq = nuniq + 1

      a1(nuniq) = a1(itest)
      a2(nuniq) = a2(itest)

    end if

  end do

  return
end
subroutine line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )
!
!*******************************************************************************
!
!! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
!
!
!  Formula:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are
!    two points on the line. (X1,Y1) must be different
!    from (X2,Y2).
!
!    Output, real A, B, C, three coefficients which describe
!    the line that passes through (X1,Y1) and (X2,Y2).
!
  implicit none
!
  real a
  real b
  real c
  real x1
  real x2
  real y1
  real y2
!
!  Take care of degenerate cases.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Fatal error!'
    write ( *, '(a)' ) '  (X1,Y1) = (X2,Y2)'
    write ( *, '(a,2g14.6)' ) '  (X1,Y1) = ', x1, y1
    write ( *, '(a,2g14.6)' ) '  (X2,Y2) = ', x2, y2
    stop
  end if

  a = y2 - y1
  b = x1 - x2
  c = x2 * y1 - x1 * y2

  return
end
subroutine line_exp2par_2d ( x1, y1, x2, y2, f, g, x0, y0 )
!
!*******************************************************************************
!
!! LINE_EXP2PAR_2D converts a line from explicit to parametric form in 2D.
!
!
!  Formula:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, two points on the line.
!
!    Output, real F, G, X0, Y0, the parametric parameters of the line.
!
  implicit none
!
  real f
  real g
  real norm
  real x0
  real x1
  real x2
  real y0
  real y1
  real y2
!
  x0 = x1
  y0 = y1

  f = x2 - x1
  g = y2 - y1

  norm = sqrt ( f * f + g * g )

  if ( norm /= 0.0E+00 ) then
    f = f / norm
    g = g / norm
  end if

  return
end
subroutine line_exp2par_3d ( x1, y1, z1, x2, y2, z2, f, g, h, x0, y0, z0 )
!
!*******************************************************************************
!
!! LINE_EXP2PAR_3D converts a line from explicit to parametric form in 3D.
!
!
!  Formula:
!
!    The explicit form of a line in 3D is:
!
!      (X1,Y1,Z1), (X2,Y2,Z2).
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z3, two points on the line.
!
!    Output, real F, G, H, X0, Y0, Z0, the parametric parameters of the line.
!
  implicit none
!
  real f
  real g
  real h
  real norm
  real x0
  real x1
  real x2
  real y0
  real y1
  real y2
  real z0
  real z1
  real z2
!
  x0 = x1
  y0 = y1
  z0 = z1

  f = x2 - x1
  g = y2 - y1
  h = z2 - z1

  norm = sqrt ( f * f + g * g + h * h )

  if ( norm /= 0.0E+00 ) then
    f = f / norm
    g = g / norm
    h = h / norm
  end if

  return
end
subroutine line_exp_normal_2d ( x1, y1, x2, y2, n1, n2 )
!
!*******************************************************************************
!
!! LINE_EXP_NORMAL_2D computes the unit normal vector to a line in 2D.
!
!
!  Formula:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Modified:
!
!    17 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, two points on the line.
!
!    Output, real N1, N2, the components of a unit normal vector to the line.
!    If the two points are equal, then N1 = N2 = 0.
!
  implicit none
!
  real n1
  real n2
  real norm
  real x1
  real x2
  real y1
  real y2
!
  norm = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

  if ( norm == 0.0E+00 ) then
    n1 = 0.0E+00
    n2 = 0.0E+00
    return
  end if

  n1 =   ( y2 - y1 ) / norm
  n2 = - ( x2 - x1 ) / norm

  return
end
subroutine line_exp_perp_2d ( x1, y1, x2, y2, x3, y3, x4, y4 )
!
!*******************************************************************************
!
!! LINE_EXP_PERP_2D computes a line perpendicular to a line and through a point.
!
!
!  Formula:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Modified:
!
!    02 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, two points on the given line.
!
!    Input, real X3, Y3, a point not on the given line, through which the
!    perpendicular must pass.
!
!    Output, real X4, Y4, a point on the given line, such that the line
!    ((X3,Y3),(X4,Y4) is perpendicular to the given line.
!
  implicit none
!
  real bot
  real enormsq0_2d
  real t
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
!
  bot = enormsq0_2d ( x1, y1, x2, y2 )

  if ( bot == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP_PERP_2D - Fatal error!'
    write ( *, '(a)' ) '  The points (X1,Y1) and (X2,Y2) are identical.'
    stop
  end if
!
!  (P3-P1) dot (P2-P1) = Norm(P3-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P3-P1) dot (P2-P1) / Norm(P3-P1)**2 = normalized coordinate T
!  of the projection of (P3-P1) onto (P2-P1).
!
  t = ( ( x1 - x3 ) * ( x1 - x2 ) + ( y1 - y3 ) * ( y1 - y2 ) ) / bot

  x4 = x1 + t * ( x2 - x1 )
  y4 = y1 + t * ( y2 - y1 )

  return
end
subroutine line_exp_point_dist_2d ( x1, y1, x2, y2, x, y, dist )
!
!*******************************************************************************
!
!! LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
!
!
!  Formula:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are two points on
!    the line.
!
!    Input, real X, Y, the point whose distance from the line is
!    to be measured.
!
!    Output, real DIST, the distance from the point to the line.
!
  implicit none
!
  real bot
  real dist
  real dot
  real enorm0_2d
  real enormsq0_2d
  real t
  real x
  real xn
  real x1
  real x2
  real y
  real yn
  real y1
  real y2
!
  bot = enormsq0_2d ( x1, y1, x2, y2 )

  if ( bot == 0.0E+00 ) then

    xn = x1
    yn = y1
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  else

    dot = ( x - x1 ) * ( x2 - x1 ) + ( y - y1 ) * ( y2 - y1 )

    t = dot / bot

    xn = x1 + t * ( x2 - x1 )
    yn = y1 + t * ( y2 - y1 )

  end if

  dist = enorm0_2d ( xn, yn, x, y )

  return
end
subroutine line_exp_point_dist_3d ( x1, y1, z1, x2, y2, z2, x, y, z, dist )
!
!*******************************************************************************
!
!! LINE_EXP_POINT_DIST_3D: distance ( explicit line, point ) in 3D.
!
!
!  Formula:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1,Z1), (X2,Y2,Z2).
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, (X1,Y1,Z1) and (X2,Y2,Z2) are
!    two points on the line.  If the points are identical, then
!    the line will be treated as a single point.
!
!    Input, real X, Y, Z, the point whose distance from the line is
!    to be measured.
!
!    Output, real DIST, the distance from the point to the line.
!
  implicit none
!
  real bot
  real dist
  real enorm0_3d
  real enormsq0_3d
  real t
  real x
  real xn
  real x1
  real x2
  real y
  real yn
  real y1
  real y2
  real z
  real zn
  real z1
  real z2
!
  bot = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )

  if ( bot == 0.0E+00 ) then

    xn = x1
    yn = y1
    zn = z1
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  else

    t = ( ( x - x1 ) * ( x2 - x1 ) + &
          ( y - y1 ) * ( y2 - y1 ) + &
          ( z - z1 ) * ( z2 - z1 ) ) / bot

    xn = x1 + t * ( x2 - x1 )
    yn = y1 + t * ( y2 - y1 )
    zn = z1 + t * ( z2 - z1 )

  end if
!
!  Now compute the distance between the projection point and P.
!
  dist = enorm0_3d ( x, y, z, xn, yn, zn )

  return
end
subroutine line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dist_signed )
!
!*******************************************************************************
!
!! LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( explicit line, point ) in 2D.
!
!
!  Discussion:
!
!    The signed distance has two interesting properties:
!
!    *  The absolute value of the signed distance is the
!        usual (Euclidean) distance.
!
!    *  Points with signed distance 0 lie on the line,
!       points with a negative signed distance lie on one side
!         of the line,
!       points with a positive signed distance lie on the
!         other side of the line.
!
!    Assuming that C is nonnegative, then if a point is a positive
!    distance away from the line, it is on the same side of the
!    line as the point (0,0), and if it is a negative distance
!    from the line, it is on the opposite side from (0,0).
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, define the two points
!    (X1,Y1) and (X2,Y2) that determine the line.
!
!    Input, real X, Y, the point (X,Y) whose signed distance is desired.
!
!    Output, real DIST_SIGNED, the signed distance from the point to the line.
!
  implicit none
!
  real a
  real b
  real c
  real dist_signed
  real x
  real x1
  real x2
  real y
  real y1
  real y2
!
!  Convert the line to implicit form.
!
  call line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )
!
!  Compute the signed distance from the point to the line.
!
  dist_signed = ( a * x + b * y + c ) / sqrt ( a * a + b * b )

  return
end
subroutine line_exp_point_near_2d ( x1, y1, x2, y2, x, y, xn, yn, dist, t )
!
!*******************************************************************************
!
!! LINE_EXP_POINT_NEAR_2D computes the point on an explicit line nearest a point in 2D.
!
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are
!    two points on the line.  (X1,Y1) must be different from (X2,Y2).
!
!    Input, real X, Y, the point whose nearest neighbor on the line is
!    to be determined.
!
!    Output, real XN, YN, the nearest point on the line to (X,Y).
!
!    Output, real DIST, the distance from the point to the line.
!
!    Output, real T, the relative position of the point
!    (XN,YN) to the points (X1,Y1) and (X2,Y2).
!
!    (XN,YN) = (1-T)*(X1,Y1) + T*(X2,Y2).
!
!    Less than 0, (XN,YN) is furthest away from (X2,Y2).
!    Between 0 and 1, (XN,YN) is between (X1,Y1) and (X2,Y2).
!    Greater than 1, (XN,YN) is furthest away from (X1,Y1).
!
  implicit none
!
  real bot
  real dist
  real enorm0_2d
  real enormsq0_2d
  real t
  real x
  real xn
  real x1
  real x2
  real y
  real yn
  real y1
  real y2
!
  bot = enormsq0_2d ( x1, y1, x2, y2 )

  if ( bot == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_POINT_NEAR_2D - Fatal error!'
    write ( *, '(a)' ) '  The points (X1,Y1) and (X2,Y2) are identical.'
    stop
  end if
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  t = ( ( x1 - x ) * ( x1 - x2 ) + ( y1 - y ) * ( y1 - y2 ) ) / bot

  xn = x1 + t * ( x2 - x1 )
  yn = y1 + t * ( y2 - y1 )

  dist = enorm0_2d ( xn, yn, x, y )

  return
end
subroutine line_exp_point_near_3d ( x1, y1, z1, x2, y2, z2, x, y, z, &
  xn, yn, zn, dist, t )
!
!*******************************************************************************
!
!! LINE_EXP_POINT_NEAR_3D: nearest point on explicit line to point in 3D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2.
!    (X1,Y1,Z1) and (X2,Y2,Z2) are two points on the line.  They
!    must be distinct.
!
!    Input, real X, Y, Z, the point whose nearest neighbor on the line
!    is to be determined.
!
!    Output, real XN, YN, ZN, the point which is the nearest point on
!    the line to (X,Y,Z).
!
!    Output, real DIST, the distance from the point to the nearest point
!    on the line.
!
!    Output, real T, the relative position of the point
!    (XN,YN,ZN) to the points (X1,Y1,Z1) and (X2,Y2,Z2).
!
!      (XN,YN,ZN) = (1-T)*(X1,Y1,Z1) + T*(X2,Y2,Z2).
!
!    Less than 0, (XN,YN,ZN) is furthest away from (X2,Y2,Z2).
!    Between 0 and 1, (XN,YN,ZN) is between (X1,Y1,Z1) and (X2,Y2,Z2).
!    Greater than 1, (XN,YN,ZN) is furthest away from (X1,Y1,Z1).
!
  implicit none
!
  real bot
  real dist
  real enorm0_3d
  real enormsq0_3d
  real t
  real x
  real xn
  real x1
  real x2
  real y
  real yn
  real y1
  real y2
  real z
  real zn
  real z1
  real z2
!
  bot = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )

  if ( bot == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP_POINT_NEAR_3D - Fatal error!'
    write ( *, '(a)' ) '  The points (X1,Y1,Z1) and (X2,Y2,Z2) '
    write ( *, '(a)' ) '  are identical.'
    stop
  end if
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  t = (   ( x - x1 ) * ( x2 - x1 ) &
        + ( y - y1 ) * ( y2 - y1 ) &
        + ( z - z1 ) * ( z2 - z1 ) ) / bot
!
!  Now compute the location of the projection point.
!
  xn = x1 + t * ( x2 - x1 )
  yn = y1 + t * ( y2 - y1 )
  zn = z1 + t * ( z2 - z1 )
!
!  Now compute the distance between the projection point and P.
!
  dist = enorm0_3d ( xn, yn, zn, x, y, z )

  return
end
subroutine line_imp2par_2d ( a, b, c, f, g, x0, y0 )
!
!*******************************************************************************
!
!! LINE_IMP2PAR_2D converts an implicit line to parametric form in 2D.
!
!
!  Formula:
!
!    The implicit form of line in 2D is:
!
!      A * X + B * Y + C = 0
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the implicit parameters of the line.
!
!    Output, real F, G, X0, Y0, the parametric parameters of the line.
!
  implicit none
!
  real a
  real b
  real c
  real f
  real g
  real x0
  real y0
!
  if ( a * a + b * b == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_IMP2PAR_2D - Fatal error!'
    write ( *, '(a)' ) '  A * A + B * B = 0.'
    stop
  end if

  x0 = - a * c / ( a * a + b * b )
  y0 = - b * c / ( a * a + b * b )
  f =   b / sqrt ( a * a + b * b )
  g = - a / sqrt ( a * a + b * b )

  return
end
subroutine line_imp_point_dist_2d ( a, b, c, x, y, dist )
!
!*******************************************************************************
!
!! LINE_IMP_POINT_DIST_2D: distance ( implicit line, point ) in 2D.
!
!
!  Formula:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the implicit line parameters.
!
!    Input, real X, Y, the point whose distance from the line is
!    to be measured.
!
!    Output, real DIST, the distance from the point to the line.
!
  implicit none
!
  real a
  real b
  real c
  real dist
  real x
  real y
!
  if ( a * a + b * b == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_IMP_POINT_DIST_2D - Fatal error!'
    write ( *, '(a)' ) '  A * A + B * B = 0.'
    stop
  end if

  dist = abs ( a * x + b * y + c ) / sqrt ( a * a + b * b )

  return
end
subroutine line_imp_point_dist_3d ( a, b, c, d, x, y, z, dist )
!
!*******************************************************************************
!
!! LINE_IMP_POINT_DIST_3D: distance ( implicit line, point ) in 3D.
!
!
!  Formula:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, the implicit line parameters.
!
!    Input, real X, Y, Z, the point whose distance from the line is
!    to be measured.
!
!    Output, real DIST, the distance from the point to the line.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real dist
  real x
  real y
  real z
!
  if ( a * a + b * b + c * c == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_IMP_POINT_DIST_3D - Fatal error!'
    write ( *, '(a)' ) '  A * A + B * B + C * C = 0.'
    stop
  end if

  dist = abs ( a * x + b * y + c * z + d ) / sqrt ( a * a + b * b + c * c )

  return
end
subroutine line_imp_point_dist_signed_2d ( a, b, c, x, y, dist_signed )
!
!*******************************************************************************
!
!! LINE_IMP_POINT_DIST_SIGNED_2D: signed distance ( implicit line, point ) in 2D.
!
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the equation of the line is A*X + B*Y + C = 0.
!
!    Input, real X, Y, the coordinates of the point.
!
!    Output, real DIST_SIGNED, the signed distance from the point to
!    the line.
!
  implicit none
!
  real a
  real b
  real c
  real dist_signed
  real x
  real y
!
  if ( a * a + b * b == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_IMP_POINT_DIST_SIGNED_2D - Fatal error!'
    write ( *, '(a)' ) '  A * A + B * B = 0.'
    stop
  end if

  dist_signed = - sign ( 1.0E+00, c ) * ( a * x + b * y + c ) / &
    sqrt ( a * a + b * b )

  return
end
subroutine line_par2imp_2d ( f, g, x0, y0, a, b, c )
!
!*******************************************************************************
!
!! LINE_PAR2IMP_2D converts a parametric line to implicit form in 2D.
!
!
!  Formula:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real F, G, X0, Y0, the parametric parameters of the line.
!
!    Output, real A, B, C, the implicit parameters of the line.
!
  implicit none
!
  real a
  real b
  real c
  real f
  real g
  real x0
  real y0
!
  a = - g
  b = f
  c = g * x0 - f * y0

  return
end
subroutine line_par_point_dist_2d ( f, g, x0, y0, x, y, dist )
!
!*******************************************************************************
!
!! LINE_PAR_POINT_DIST_2D: distance ( parametric line, point ) in 2D.
!
!
!  Formula:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real F, G, X0, Y0, the parametric line parameters.
!
!    Input, real X, Y, the point whose distance from the line is
!    to be measured.
!
!    Output, real DIST, the distance from the point to the line.
!
  implicit none
!
  real dist
  real dx
  real dy
  real f
  real g
  real x
  real x0
  real y
  real y0
!
  dx =   g * g * ( x - x0 ) - f * g * ( y - y0 )
  dy = - f * g * ( x - x0 ) + f * f * ( y - y0 )

  dist = sqrt ( dx * dx + dy * dy ) / ( f * f + g * g )

  return
end
subroutine line_par_point_dist_3d ( f, g, h, x0, y0, z0, x, y, z, dist )
!
!*******************************************************************************
!
!! LINE_PAR_POINT_DIST_3D: distance ( parametric line, point ) in 3D.
!
!
!  Formula:
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real F, G, H, X0, Y0, Z0, the parametric line parameters.
!
!    Input, real X, Y, Z, the point whose distance from the line is
!    to be measured.
!
!    Output, real DIST, the distance from the point to the line.
!
  implicit none
!
  real dist
  real dx
  real dy
  real dz
  real f
  real g
  real h
  real x
  real x0
  real y
  real y0
  real z
  real z0
!
  dx =   g * ( f * ( y - y0 ) - g * ( x - x0 ) ) &
       + h * ( f * ( z - z0 ) - h * ( x - x0 ) )

  dy =   h * ( g * ( z - z0 ) - h * ( y - y0 ) )  &
       - f * ( f * ( y - y0 ) - g * ( x - x0 ) )

  dz = - f * ( f * ( z - z0 ) - h * ( x - x0 ) ) &
       - g * ( g * ( z - z0 ) - h * ( y - y0 ) )

  dist = sqrt ( dx * dx + dy * dy + dz * dz ) &
    / ( f * f + g * g + h * h )

  return
end
subroutine line_seg_contains_point_1d ( x1, x2, x3, u )
!
!*******************************************************************************
!
!! LINE_SEG_CONTAINS_POINT_1D reports if a line segment contains a point in 1D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, X2, two points defining a line segment.
!    The line segment has origin at X1, and unit at X2.
!
!    Input, real X3, a point to be tested.
!
!    Output, real U, the coordinate of X3 in units of (X2-X1).
!    The point X3 is contained in the line segment if 0 <= U <= 1.
!
  implicit none
!
  real u
  real unit
  real x1
  real x2
  real x3
!
  unit = x2 - x1

  if ( unit == 0.0E+00 ) then

    if ( x3 == x1 ) then
      u = 0.5E+00
    else if ( x3 < x1 ) then
      u = - huge ( u )
    else if ( x3 > x1 ) then
      u = huge ( u )
    end if

  else

    u = ( x3 - x1 ) / unit

  end if

  return
end
subroutine line_seg_contains_point_2d ( x1, y1, x2, y2, x3, y3, u, v )
!
!*******************************************************************************
!
!! LINE_SEG_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!    In exact arithmetic, point P3 = (X3,Y3) is on the line segment between
!    P1=(X1,Y1) and P2=(X2,Y2) if and only if 0 <= U <= 1 and V = 0.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the endpoints P1 and P2 of a line segment.
!
!    Input, real X3, Y3, a point P3 to be tested.
!
!    Output, real U, the coordinate of (X3,Y3) along the axis from
!    with origin at P1 and unit at P2.
!
!    Output, real V, the magnitude of the off-axis portion of the
!    vector P3-P1, measured in units of (P2-P1).
!
  implicit none
!
  real u
  real unit
  real v
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
!
  unit = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

  if ( unit == 0.0E+00 ) then

    if ( x3 == x1 .and. y3 == y1 ) then
      u = 0.5E+00
      v = 0.0E+00
    else
      u = 0.5E+00
      v = huge ( v )
    end if

  else

    u = ( ( x3 - x1 ) * ( x2 - x1 ) + ( y3 - y1 ) * ( y2 - y1 ) ) / unit**2

    v = sqrt ( ( ( u - 1.0E+00 ) * x1 - u * x2 + x3 )**2 &
             + ( ( u - 1.0E+00 ) * y1 - u * y2 + y3 )**2 ) / unit

  end if

  return
end
subroutine line_seg_point_dist_2d ( x1, y1, x2, y2, x, y, dist )
!
!*******************************************************************************
!
!! LINE_SEG_POINT_DIST_2D: distance ( line segment, point ) in 2D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the endpoints of the line segment.
!
!    Input, real X, Y, the point whose nearest neighbor on the line
!    segment is to be determined.
!
!    Output, real DIST, the distance from the point to the line segment.
!
  implicit none
!
  real bot
  real dist
  real enorm0_2d
  real enormsq0_2d
  real t
  real x
  real xn
  real x1
  real x2
  real y
  real yn
  real y1
  real y2
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( x1 == x2 .and. y1 == y2 ) then

    xn = x1
    yn = y1

  else

    bot = enormsq0_2d ( x1, y1, x2, y2 )

    t = ( ( x1 - x ) * ( x1 - x2 ) + ( y1 - y ) * ( y1 - y2 ) ) / bot

    t = max ( t, 0.0E+00 )
    t = min ( t, 1.0E+00 )

    xn = x1 + t * ( x2 - x1 )
    yn = y1 + t * ( y2 - y1 )

  end if

  dist = enorm0_2d ( xn, yn, x, y )

  return
end
subroutine line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, x, y, z, dist )
!
!*******************************************************************************
!
!! LINE_SEG_POINT_DIST_3D: distance ( line segment, point ) in 3D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the endpoints of the line
!    segment.
!
!    Input, real X, Y, Z, the point whose nearest neighbor on the line
!    segment is to be determined.
!
!    Output, real DIST, the distance from the point to the line segment.
!
  implicit none
!
  real bot
  real dist
  real enorm0_3d
  real enormsq0_3d
  real t
  real x
  real xn
  real x1
  real x2
  real y
  real yn
  real y1
  real y2
  real z
  real zn
  real z1
  real z2
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( x1 == x2 .and. y1 == y2 .and. z1 == z2 ) then

    xn = x1
    yn = y1
    zn = z1

  else

    bot = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )

    t = (   ( x1 - x ) * ( x1 - x2 ) &
          + ( y1 - y ) * ( y1 - y2 ) &
          + ( z1 - z ) * ( z1 - z2 ) ) / bot

    t = max ( t, 0.0E+00 )
    t = min ( t, 1.0E+00 )

    xn = x1 + t * ( x2 - x1 )
    yn = y1 + t * ( y2 - y1 )
    zn = z1 + t * ( z2 - z1 )

  end if

  dist = enorm0_3d ( x, y, z, xn, yn, zn )

  return
end
subroutine line_seg_point_near_2d ( x1, y1, x2, y2, x, y, xn, yn, dist, t )
!
!*******************************************************************************
!
!! LINE_SEG_POINT_NEAR_2D finds the point on a line segment nearest a point in 2D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the two endpoints of the line segment.
!
!    (X1,Y1) should generally be different from (X2,Y2), but
!    if they are equal, the program will still compute a
!    meaningful result.
!
!    Input, real X, Y, the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real XN, YN, the point on the line segment which is
!    nearest the point (X,Y).
!
!    Output, real DIST, the distance from the point to the nearest point
!    on the line segment.
!
!    Output, real T, the relative position of the point (XN,YN) to the
!    points (X1,Y1) and (X2,Y2).
!
!      (XN,YN) = (1-T) * (X1,Y1) + T * (X2,Y2).
!
!    T will always be between 0 and 1.
!
  implicit none
!
  real bot
  real dist
  real enorm0_2d
  real enormsq0_2d
  real t
  real x
  real xn
  real x1
  real x2
  real y
  real yn
  real y1
  real y2
!
  if ( x1 == x2 .and. y1 == y2 ) then

    t = 0.0E+00
    xn = x1
    yn = y1

  else

    bot = enormsq0_2d ( x1, y1, x2, y2 )

    t = (   ( x1 - x ) * ( x1 - x2 ) &
          + ( y1 - y ) * ( y1 - y2 ) ) / bot

    if ( t < 0.0E+00 ) then
      t = 0.0E+00
      xn = x1
      yn = y1
    else if ( t > 1.0E+00 ) then
      t = 1.0E+00
      xn = x2
      yn = y2
    else
      xn = x1 + t * ( x2 - x1 )
      yn = y1 + t * ( y2 - y1 )
    end if

  end if

  dist = enorm0_2d ( x, y, xn, yn )

  return
end
subroutine line_seg_point_near_3d ( x1, y1, z1, x2, y2, z2, x, y, z, &
  xn, yn, zn, dist, t )
!
!*******************************************************************************
!
!! LINE_SEG_POINT_NEAR_3D finds the point on a line segment nearest a point in 3D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the two endpoints of the line segment.
!
!    (X1,Y1,Z1) should generally be different from (X2,Y2,Z2), but
!    if they are equal, the program will still compute a meaningful
!    result.
!
!    Input, real X, Y, Z, the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real XN, YN, ZN, the point on the line segment which is
!    nearest the point (X,Y,Z).
!
!    Output, real DIST, the distance from the point to the nearest point
!    on the line segment.
!
!    Output, real T, the relative position of the nearest point
!    (XN,YN,ZN) to the defining points (X1,Y1,Z1) and (X2,Y2,Z2).
!
!      (XN,YN,ZN) = (1-T)*(X1,Y1,Z1) + T*(X2,Y2,Z2).
!
!    T will always be between 0 and 1.
!
  implicit none
!
  real bot
  real dist
  real enorm0_3d
  real enormsq0_3d
  real t
  real x
  real xn
  real x1
  real x2
  real y
  real yn
  real y1
  real y2
  real z
  real zn
  real z1
  real z2
!
  if ( x1 == x2 .and. y1 == y2 .and. z1 == z2 ) then

    t = 0.0E+00
    xn = x1
    yn = y1
    zn = z1

  else

    bot = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )

    t = (   ( x1 - x ) * ( x1 - x2 ) &
          + ( y1 - y ) * ( y1 - y2 ) &
          + ( z1 - z ) * ( z1 - z2 ) ) / bot

    if ( t < 0.0E+00 ) then
      t = 0.0E+00
      xn = x1
      yn = y1
      zn = z1
    else if ( t > 1.0E+00 ) then
      t = 1.0E+00
      xn = x2
      yn = y2
      zn = z2
    else
      xn = x1 + t * ( x2 - x1 )
      yn = y1 + t * ( y2 - y1 )
      zn = z1 + t * ( z2 - z1 )
    end if

  end if

  dist = enorm0_3d ( xn, yn, zn, x, y, z )

  return
end
subroutine lines_exp_angle_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, angle )
!
!*******************************************************************************
!
!! LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
!
!
!  Formula:
!
!    The explicit form of a line in 3D is:
!
!      (X1,Y1,Z1), (X2,Y2,Z2).
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, two distince points on the first line.
!
!    Input, real X3, Y3, Z3, X4, Y4, Z4, two distinct points on the second line.
!
!    Output, real ANGLE, the angle in radians between the two lines.
!    ANGLE is computed using the ACOS function, and so lies between 0 and PI.
!    But if one of the lines is degenerate, ANGLE is returned as -1.0.
!
  implicit none
!
  real angle
  real arc_cosine
  real ctheta
  real enorm0_3d
  real pdotq
  real pnorm
  real qnorm
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
  real z1
  real z2
  real z3
  real z4
!
  pnorm = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
  qnorm = enorm0_3d ( x3, y3, z3, x4, y4, z4 )

  pdotq =    ( x2 - x1 ) * ( x4 - x3 ) &
           + ( y2 - y1 ) * ( y4 - y3 ) &
           + ( z2 - z1 ) * ( z4 - z3 )

  if ( pnorm == 0.0E+00 .or. qnorm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Fatal error!'
    write ( *, '(a)' ) '  One of the lines is degenerate!'
    angle = - 1.0E+00

  else

    ctheta = pdotq / ( pnorm * qnorm )
    angle = arc_cosine ( ctheta )

  end if

  return
end
subroutine lines_exp_angle_nd ( p1, p2, q1, q2, n, angle )
!
!*******************************************************************************
!
!! LINES_EXP_ANGLE_ND returns the angle between two explicit lines in ND.
!
!
!  Modified:
!
!    27 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real P1(N), P2(N), two points on the first line.
!
!    Input, real Q1(N), Q2(N), two points on the second line.
!
!    Input, integer N, the dimension of the space.
!
!    Output, real ANGLE, the angle in radians between the two lines.
!    ANGLE is computed using the ACOS function, and so lies between 0 and PI.
!    But if one of the lines is degenerate, ANGLE is returned as -1.0.
!
  implicit none
!
  integer n
!
  real angle
  real arc_cosine
  real ctheta
  real enorm0_nd
  integer i
  real p1(n)
  real p2(n)
  real pdotq
  real pnorm
  real q1(n)
  real q2(n)
  real qnorm
!
  pnorm = enorm0_nd ( n, p1, p2 )
  qnorm = enorm0_nd ( n, q1, q2 )

  pdotq = 0.0E+00
  do i = 1, n
    pdotq = pdotq + ( p2(i) - p1(i) ) * ( q2(i) - q1(i) )
  end do

  if ( pnorm == 0.0E+00 .or. qnorm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINES_EXP_ANGLE_ND - Fatal error!'
    write ( *, '(a)' ) '  One of the lines is degenerate!'
    angle = - 1.0E+00
  else
    ctheta = pdotq / ( pnorm * qnorm )
    angle = arc_cosine ( ctheta )

  end if

  return
end
subroutine lines_exp_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
  dist )
!
!*******************************************************************************
!
!! LINES_EXP_DIST_3D computes the distance between two explicit lines in 3D.
!
!
!  Formula:
!
!    The explicit form of a line in 3D is:
!
!      (X1,Y1,Z1), (X2,Y2,Z2).
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, (X1,Y1,Z1) and (X2,Y2,Z2) are
!    two points on the first line.  They must be distinct.
!
!    Input, real X3, Y3, Z3, X4, Y4, Z4, (X3,Y3,Z3) and (X4,Y4,Z4) are
!    two points on the second line.  They must be distinct.
!
!    Output, real DIST, the distance between the lines.
!
  implicit none
!
  real a11
  real a12
  real a13
  real a21
  real a22
  real a23
  real a31
  real a32
  real a33
  real bot
  real cr1
  real cr2
  real cr3
  real dist
  real top
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
  real z1
  real z2
  real z3
  real z4
!
!  The distance is found by computing the volume of a parallelipiped,
!  and dividing by the area of its base.
!
!  But if the lines are parallel, we compute the distance by
!  finding the distance between the first line and any point
!  on the second line.
!
  a11 = x3 - x1
  a12 = y3 - y1
  a13 = z3 - z1

  a21 = x2 - x1
  a22 = y2 - y1
  a23 = z2 - z1

  a31 = x4 - x3
  a32 = y4 - y3
  a33 = z4 - z3

  call cross_3d ( a21, a22, a23, a31, a32, a33, cr1, cr2, cr3 )

  bot = sqrt ( cr1 * cr1 + cr2 * cr2 + cr3 * cr3 )

  if ( bot == 0.0E+00 ) then

    call line_exp_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, dist )

  else

    top = abs (   a11 * ( a22 * a33 - a23 * a32 ) &
                - a12 * ( a21 * a33 - a23 * a31 ) &
                + a13 * ( a21 * a32 - a22 * a31 ) )

    dist = top / bot

  end if

  return
end
subroutine lines_exp_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, ival, x, y )
!
!*******************************************************************************
!
!! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
!
!
!  Formula:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Modified:
!
!    16 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, define the first line.
!
!    Input, real X3, Y3, X4, Y4, define the second line.
!
!    Output, integer IVAL, reports on the intersection:
!    0, no intersection, the lines may be parallel or degenerate.
!    1, one intersection point, returned in X, Y.
!    2, infinitely many intersections, the lines are identical.
!
!    Output, real X, Y, if IVAl = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
  implicit none
!
  real a1
  real a2
  real b1
  real b2
  real c1
  real c2
  integer ival
  logical point_1
  logical point_2
  real x
  real x1
  real x2
  real x3
  real x4
  real y
  real y1
  real y2
  real y3
  real y4
!
  ival = 0
  x = 0.0E+00
  y = 0.0E+00
!
!  Check whether either line is a point.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    point_1 = .true.
  else
    point_1 = .false.
  end if

  if ( x3 == x4 .and. y3 == y4 ) then
    point_2 = .true.
  else
    point_2 = .false.
  end if
!
!  Convert the lines to ABC format.
!
  if ( .not. point_1 ) then
    call line_exp2imp_2d ( x1, y1, x2, y2, a1, b1, c1 )
  end if

  if ( .not. point_2 ) then
    call line_exp2imp_2d ( x3, y3, x4, y4, a2, b2, c2 )
  end if
!
!  Search for intersection of the lines.
!
  if ( point_1 .and. point_2 ) then
    if ( x1 == x3 .and. y1 == y3 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_1 ) then
    if ( a2 * x1 + b2 * y1 == c2 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_2 ) then
    if ( a1 * x3 + b1 * y3 == c1 ) then
      ival = 1
      x = x3
      y = y3
    end if
  else
    call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )
  end if

  return
end
function lines_exp_parallel_2d ( x1, y1, x2, y2, x3, y3, x4, y4 )
!
!*******************************************************************************
!
!! LINES_EXP_PARALLEL_2D determines if two lines are parallel in 2D.
!
!
!  Discussion:
!
!    The test is essentially a comparison of slopes, but should be
!    more accurate than an explicit slope comparison, and unfazed
!    by degenerate cases.
!
!    If the lines are determined to be parallel, then you can
!    determine whether they are identical or distinct by evaluating:
!
!      lines_exp_parallel_2d ( x1, y1, x4, y4, x3, y3, x2, y2 )
!
!  Modified:
!
!    09 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, define the first line.
!
!    Input, real X3, Y3, X4, Y4, define the second line.
!
!    Output, logical LINES_EXP_PARALLEL_2D is TRUE if the lines are parallel.
!
  implicit none
!
  logical lines_exp_parallel_2d
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
!
  lines_exp_parallel_2d = &
    ( ( y2 - y1 ) * ( x4 - x3 ) == ( y4 - y3 ) * ( x2 - x1 ) )

  return
end
subroutine lines_imp_angle_2d ( a1, b1, c1, a2, b2, c2, theta )
!
!*******************************************************************************
!
!! LINES_IMP_ANGLE_2D finds the angle between two implicit lines in 2D.
!
!
!  Formula:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A1, B1, C1, the implicit parameters of the first line.
!
!    Input, real A2, B2, C2, the implicit parameters of the second line.
!
!    Output, real THETA, the angle between the two lines.
!
  implicit none
!
  real a1
  real a2
  real arc_cosine
  real b1
  real b2
  real c1
  real c2
  real pdotq
  real pnorm
  real qnorm
  real theta
!
  pdotq = a1 * a2 + b1 * b2
  pnorm = sqrt ( a1 * a1 + b1 * b1 )
  qnorm = sqrt ( a2 * a2 + b2 * b2 )

  theta = arc_cosine ( pdotq / ( pnorm * qnorm ) )

  return
end
subroutine lines_imp_dist_2d ( a1, b1, c1, a2, b2, c2, dist )
!
!*******************************************************************************
!
!! LINES_IMP_DIST_2D determines the distance between two implicit lines in 2D.
!
!
!  Discussion:
!
!    If the lines intersect, then their distance is zero.
!    If the two lines are parallel, then they have a nonzero distance.
!
!  Formula:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Modified:
!
!    12 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, real DIST, the distance between the two lines.
!
  implicit none
!
  real a1
  real a2
  real b1
  real b2
  real c1
  real c2
  real dist
!
!  Refuse to handle degenerate lines.
!
  if ( a1 == 0.0E+00 .and. b1 == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINES_IMP_DIST_2D - Fatal error!'
    write ( *, '(a)' ) '  Line 1 is degenerate.'
    stop
  else if ( a2 == 0.0E+00 .and. b2 == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINES_IMP_DIST_2D - Fatal error!'
    write ( *, '(a)' ) '  Line 2 is degenerate.'
    stop
  end if
!
!  Determine if the lines intersect.
!
  if ( a1 * b2 /= a2 * b1 ) then
    dist = 0.0E+00
    return
  end if
!
!  Determine the distance between the parallel lines.
!
  dist = abs ( c2 / sqrt ( a2**2 + b2**2 ) - c1 / sqrt ( a1**2 + b1**2 ) )

  return
end
subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )
!
!*******************************************************************************
!
!! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
!
!
!  Formula:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, integer IVAL, reports on the intersection.
!
!    -1, both A1 and B1 were zero.
!    -2, both A2 and B2 were zero.
!     0, no intersection, the lines are parallel.
!     1, one intersection point, returned in X, Y.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, real X, Y, if IVAL = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
  implicit none
!
  real a(2,2)
  real a1
  real a2
  real b(2,2)
  real b1
  real b2
  real c1
  real c2
  real det
  integer ival
  real x
  real y
!
  x = 0.0E+00
  y = 0.0E+00
!
!  Refuse to handle degenerate lines.
!
  if ( a1 == 0.0E+00 .and. b1 == 0.0E+00 ) then
    ival = -1
    return
  else if ( a2 == 0.0E+00 .and. b2 == 0.0E+00 ) then
    ival = -2
    return
  end if
!
!  Set up a linear system, and compute its inverse.
!
  a(1,1) = a1
  a(1,2) = b1
  a(2,1) = a2
  a(2,2) = b2

  call rmat2_inverse ( a, b, det )
!
!  If the inverse exists, then the lines intersect.
!  Multiply the inverse times -C to get the intersection point.
!
  if ( det /= 0.0E+00 ) then

    ival = 1
    x = - b(1,1) * c1 - b(1,2) * c2
    y = - b(2,1) * c1 - b(2,2) * c2
!
!  If the inverse does not exist, then the lines are parallel
!  or coincident.  Check for parallelism by seeing if the
!  C entries are in the same ratio as the A or B entries.
!
  else

    ival = 0

    if ( a1 == 0.0E+00 ) then
      if ( b2 * c1 == c2 * b1 ) then
        ival = 2
      end if
    else
      if ( a2 * c1 == c2 * a1 ) then
        ival = 2
      end if
    end if

  end if

  return
end
subroutine lines_par_angle_2d ( f1, g1, x01, y01, f2, g2, x02, y02, theta )
!
!*******************************************************************************
!
!! LINES_PAR_ANGLE_2D finds the angle between two parametric lines in 2D.
!
!
!  Formula:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real F1, G1, X01, Y01, the parametric parameters of the
!    first line.
!
!    Input, real F2, G2, X02, Y02, the parametric parameters of the
!    second line.
!
!    Output, real THETA, the angle between the two lines.
!
  implicit none
!
  real arc_cosine
  real f1
  real f2
  real g1
  real g2
  real pdotq
  real pnorm
  real qnorm
  real theta
  real x01
  real x02
  real y01
  real y02
!
  pdotq = f1 * f2 + g1 * g2
  pnorm = sqrt ( f1 * f1 + g1 * g1 )
  qnorm = sqrt ( f2 * f2 + g2 * g2 )

  theta = arc_cosine ( pdotq / ( pnorm * qnorm ) )

  return
end
subroutine lines_par_angle_3d ( f1, g1, h1, x01, y01, z01, f2, g2, h2, &
  x02, y02, z02, theta )
!
!*******************************************************************************
!
!! LINES_PAR_ANGLE_3D finds the angle between two parametric lines in 3D.
!
!
!  Formula:
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real F1, G1, H1, X01, Y01, Z01, the parametric parameters
!    of the first line.
!
!    Input, real F2, G2, H2, X02, Y02, Z02, the parametric parameters
!    of the second line.
!
!    Output, real THETA, the angle between the two lines.
!
  implicit none
!
  real arc_cosine
  real f1
  real f2
  real g1
  real g2
  real h1
  real h2
  real pdotq
  real pnorm
  real qnorm
  real theta
  real x01
  real x02
  real y01
  real y02
  real z01
  real z02
!
  pdotq = f1 * f2 + g1 * g2 + h1 * h2
  pnorm = sqrt ( f1 * f1 + g1 * g1 + h1 * h1 )
  qnorm = sqrt ( f2 * f2 + g2 * g2 + h2 * h2 )

  theta = arc_cosine ( pdotq / ( pnorm * qnorm ) )

  return
end
subroutine lines_par_dist_3d ( f1, g1, h1, x01, y01, z01, f2, g2, h2, &
  x02, y02, z02, dist )
!
!*******************************************************************************
!
!! LINES_PAR_DIST_3D finds the distance between two parametric lines in 3D.
!
!
!  Formula:
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!  Warning:
!
!    This code does not work for parallel or near parallel lines.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real F1, G1, H1, X01, Y01, Z01, the parametric parameters
!    of the first line.
!
!    Input, real F2, G2, H2, X02, Y02, Z02, the parametric parameters
!    of the second line.
!
!    Output, real DIST, the distance between the two lines.
!
  implicit none
!
  real dist
  real f1
  real f2
  real g1
  real g2
  real h1
  real h2
  real x01
  real x02
  real y01
  real y02
  real z01
  real z02
!
  dist = abs ( ( x02 - x01 ) * ( g1 * h2 - g2 * h1 ) &
             + ( y02 - y01 ) * ( h1 * f2 - h2 * f1 ) &
             + ( z02 - z01 ) * ( f1 * g2 - f2 * g1 ) )  / &
             ( ( f1 * g2 - f2 * g1 )**2 &
             + ( g1 * h2 - g2 * h1 )**2 &
             + ( h1 * f2 - h2 * f1 )**2 )

  return
end
subroutine lines_par_int_2d ( f1, g1, x1, y1, f2, g2, x2, y2, t1, t2, &
  xint, yint )
!
!*******************************************************************************
!
!! LINES_PAR_INT_2D determines where two parametric lines intersect in 2D.
!
!
!  Formula:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real F1, G1, X1, Y1, define the first parametric line.
!
!    Input, real F2, G2, X2, Y2, define the second parametric line.
!
!    Output, real DET, the determinant of the system.  If DET = 0,
!    then the lines are parallel and no intersection was computed.
!
!    Output, real T1, T2, the T parameters on the first and second
!    lines of the intersection point.
!
!    Output, real XINT, YINT, the intersection point.
!
  implicit none
!
  real det
  real f1
  real f2
  real g1
  real g2
  real t1
  real t2
  real x1
  real x2
  real xint
  real y1
  real y2
  real yint
!
  det = f2 * g1 - f1 * g2

  if ( det == 0.0E+00 ) then
    t1 = 0.0E+00
    t2 = 0.0E+00
    xint = 0.0E+00
    yint = 0.0E+00
  else
    t1 = ( f2 * ( y2 - y1 ) - g2 * ( x2 - x1 ) ) / det
    t2 = ( f1 * ( y2 - y1 ) - g1 * ( x2 - x1 ) ) / det
    xint = x1 + f1 * t1
    yint = y1 + g1 * t1
  end if

  return
end
subroutine lines_seg_dist_2d ( x1, y1, x2, y2, x3, y3, x4, y4, dist )
!
!*******************************************************************************
!
!! LINES_SEG_DIST_2D computes the distance between two line segments in 2D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    29 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the endpoints of the first segment.
!
!    Input, real X3, Y3, X4, Y4, the endpoints of the second segment.
!
!    Output, real DIST, the distance between the line segments.
!
  implicit none
!
  real dist
  real dl
  real dm
  real dr
  real tl
  real tm
  real tmin
  real tr
  real x1
  real x2
  real x3
  real x4
  real xm
  real xmin
  real y1
  real y2
  real y3
  real y4
  real ym
  real ymin
!
!  Label the left, middle and right of line 1 as L1, M1, R1.
!
!  Compute the distance from each point to the opposite line segment.
!
  call line_seg_point_dist_2d ( x3, y3, x4, y4, x1, y1, dl )

  xm = 0.5E+00 * ( x1 + x2 )
  ym = 0.5E+00 * ( y1 + y2 )

  call line_seg_point_dist_2d ( x3, y3, x4, y4, xm, ym, dm )

  call line_seg_point_dist_2d ( x3, y3, x4, y4, x2, y2, dr )
!
!  Now find the "theoretical" minimum of the distance function.
!
  tl = 0.0E+00
  tm = 0.5E+00
  tr = 1.0E+00
  call minabs ( tl, dl, tm, dm, tr, dr, tmin, dist )
!
!  Evaluate the distance at the minimum, to account for
!  a flattening of the absolute value function caused by parallel
!  or coincident portions of the lines.
!
  xmin = x1 + tmin * ( x2 - x1 )
  ymin = y1 + tmin * ( y2 - y1 )
  call line_seg_point_dist_2d ( x3, y3, x4, y4, xmin, ymin, dist )

  return
end
subroutine lines_seg_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, dist )
!
!*******************************************************************************
!
!! LINES_SEG_DIST_3D computes the distance between two line segments in 3D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    03 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the endpoints of the first segment.
!
!    Input, real X3, Y3, Z3, X4, Y4, Z4, the endpoints of the second segment.
!
!    Output, real DIST, the distance between the line segments.
!
  implicit none
!
  real d1
  real d2
  real dist
  real dl
  real dm
  real dr
  real t1
  real t2
  real tl
  real tm
  real tmin
  real tr
  real x1
  real x2
  real x3
  real x4
  real xn1
  real xn2
  real y1
  real y2
  real y3
  real y4
  real yn1
  real yn2
  real z1
  real z2
  real z3
  real z4
  real zn1
  real zn2
!
!  Find the nearest points on line 2 to the endpoints of line 1.
!
  call line_seg_point_near_3d ( x3, y3, z3, x4, y4, z4, x1, y1, z1, &
    xn1, yn1, zn1, d1, t1 )

  call line_seg_point_near_3d ( x3, y3, z3, x4, y4, z4, x2, y2, z2, &
     xn2, yn2, zn2, d2, t2 )

  if ( t1 == t2 ) then
    call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, dist )
    return
  end if
!
!  On line 2, over the interval between the points nearest to line 1,
!  the square of the distance of any point to line 1 is a quadratic function.
!  Evaluate it at three points, and seek its local minimum.
!
  call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, dl )

  call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, &
    0.5E+00 * ( xn1 + xn2 ), 0.5E+00 * ( yn1 + yn2 ), &
    0.5E+00 * ( zn1 + zn2 ), dm )

  call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, xn2, yn2, zn2, dr )

  tl = 0.0E+00
  tm = 0.5E+00
  tr = 1.0E+00

  dl = dl * dl
  dm = dm * dm
  dr = dr * dr

  call minquad ( tl, dl, tm, dm, tr, dr, tmin, dist )

  dist = sqrt ( dist )

  return
end
subroutine lines_seg_int_1d ( x1, x2, x3, x4, flag, x5, x6 )
!
!*******************************************************************************
!
!! LINES_SEG_INT_1D computes the intersection of two line segments in 1D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!    In 1D, two line segments "intersect" if they overlap.
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, X2, the endpoints of the first segment.
!
!    Input, real X3, X4, the endpoints of the second segment.
!
!    Output, integer FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real X5, X6, the endpoints of the intersection segment.
!    If FLAG = 0, X5 = X6 = 0.
!
  implicit none
!
  integer flag
  real x1
  real x2
  real x3
  real x4
  real x5
  real x6
  real y1
  real y2
  real y3
  real y4
!
  y1 = min ( x1, x2 )
  y2 = max ( x1, x2 )
  y3 = min ( x3, x4 )
  y4 = max ( x3, x4 )

  flag = 0
  x5 = 0.0E+00
  x6 = 0.0E+00

  if ( y4 < y1 ) then
    return
  else if ( y3 > y2 ) then
    return
  end if

  flag = 1
  x5 = max ( y1, y3 )
  x6 = min ( y2, y4 )

  return
end
subroutine lines_seg_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5 )
!
!*******************************************************************************
!
!! LINES_SEG_INT_2D computes the intersection of two line segments in 2D.
!
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!    In 2D, two line segments might not intersect, even though the
!    lines, of which they are portions, intersect.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the endpoints of the first segment.
!
!    Input, real X3, Y3, X4, Y4, the endpoints of the second segment.
!
!    Output, integer FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real X5, Y5.
!    If FLAG = 0, X5 = Y5 = 0.
!    If FLAG = 1, then (X5,Y5) is a point of intersection.
!
  implicit none
!
  integer flag
  integer ival
  real u
  real v
  real x1
  real x2
  real x3
  real x4
  real x5
  real y1
  real y2
  real y3
  real y4
  real y5
!
  x5 = 0.0E+00
  y5 = 0.0E+00
!
!  Find the intersection of the two lines.
!
  call lines_exp_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, ival, x5, y5 )

  if ( ival == 0 ) then
    flag = 0
    return
  end if
!
!  Is the point on the first segment?
!
  call line_seg_contains_point_2d ( x1, y1, x2, y2, x5, y5, u, v )

  if ( u < 0.0E+00 .or. 1.0E+00 < u .or. v > 0.001E+00 ) then
    flag = 0
    return
  end if
!
!  Is the point on the second segment?
!
  call line_seg_contains_point_2d ( x3, y3, x4, y4, x5, y5, u, v )

  if ( u < 0.0E+00 .or. 1.0E+00 < u .or. v > 0.001E+00 ) then
    flag = 0
    return
  end if

  flag = 1

  return
end
subroutine loc2glob_3d ( cospitch, cosroll, cosyaw, locpts, globas, glopts, &
  sinpitch, sinroll, sinyaw )
!
!*******************************************************************************
!
!! LOC2GLOB_3D converts from a local to global coordinate system in 3D.
!
!
!  Discussion:
!
!    A global coordinate system is given.
!
!    A local coordinate system has been translated to the point with
!    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
!    a roll.
!
!    A point has local coordinates LOCPTS, and it is desired to know
!    the point's global coordinates GLOPTS.
!
!    The transformation may be written as
!
!      GLOB = GLOBAS + N_YAW * N_PITCH * N_ROLL * LOC
!
!    where
!
!               (  cos(Yaw)   -sin(Yaw)        0      )
!    N_YAW    = (  sin(Yaw)    cos(Yaw)        0      )
!               (      0           0           1      )
!
!               (  cos(Pitch)      0       sin(Pitch) )
!    N_PITCH =  (      0           1           0      )
!               ( -sin(Pitch)      0       cos(Pitch) )
!
!               (      1           0           0      )
!    N_ROLL =   (      0       cos(Roll)  -sin(Roll)  )
!               (      0       sin(Roll)   cos(Roll)  )
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COSPITCH, COSROLL, COSYAW, the cosines of the pitch,
!    roll and yaw angles.
!
!    Input, real GLOBAS(3), the global coordinates of the base vector.
!
!    Output, real GLOPTS(3), the global coordinates of the point.
!
!    Input, real LOCPTS(3), the local coordinates of the point.
!
!    Input, real SINPITCH, SINROLL, SINYAW, the sines of the pitch,
!    roll and yaw angles.
!
  implicit none
!
  real cospitch
  real cosroll
  real cosyaw
  real globas(3)
  real glopts(3)
  real locpts(3)
  real sinpitch
  real sinroll
  real sinyaw
!
  glopts(1) = globas(1) + (  cosyaw * cospitch ) * locpts(1) &
    + (  cosyaw * sinpitch * sinroll - sinyaw * cosroll ) * locpts(2) &
    + (  cosyaw * sinpitch * cosroll + sinyaw * sinroll ) * locpts(3)

  glopts(2) = globas(2) + (  sinyaw * cospitch ) * locpts(1) &
    + (  sinyaw * sinpitch * sinroll + cosyaw * cosroll ) * locpts(2) &
    + (  sinyaw * sinpitch * cosroll - cosyaw * sinroll ) * locpts(3)

  glopts(3) = globas(3) + ( -sinpitch ) * locpts(1) &
    + (  cospitch * sinroll ) * locpts(2) &
    + (  cospitch * cosroll ) * locpts(3)

  return
end
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )
!
!*******************************************************************************
!
!! LRLINE determines where a point lies in relation to a directed line.
!
!
!  Discussion:
!
!    LRLINE determines whether a point is to the left of, right of,
!    or on a directed line parallel to a line through given points.
!
!  Modified:
!
!    18 June 2001
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, real XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
!    directed line is parallel to and at signed distance DV to the left of
!    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
!    which the position relative to the directed line is to be determined.
!
!    Input, real DV, the signed distance, positive for left.
!
!    Output, integer LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
!    to the right of, on, or left of the directed line.  LRLINE is 0 if
!    the line degenerates to a point.
!
  implicit none
!
  real dv
  real dx
  real dxu
  real dy
  real dyu
  integer lrline
  real t
  real, parameter :: tol = 0.0000001E+00
  real tolabs
  real xu
  real xv1
  real xv2
  real yu
  real yv1
  real yv2
!
  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), abs ( dyu ), &
    abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx**2 + dy**2 )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else if ( t < -tolabs ) then
    lrline = -1
  end if

  return
end
subroutine minabs ( x1, y1, x2, y2, x3, y3, xmin, ymin )
!
!*******************************************************************************
!
!! MINABS finds a local minimum of F(X) = A * abs ( X ) + B.
!
!
!  Modified:
!
!    28 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X1, Y1, X2, Y2, X3, Y3, are three sets of data
!    of the form ( X, F(X) ).  The three X values must be distinct.
!    On output, the data has been sorted so that X1 < X2 < X3,
!    and the Y values have been rearranged accordingly.
!
!    Output, real XMIN, YMIN.  XMIN is a point within the interval
!    spanned by X1, X2 and X3, at which F takes its local minimum
!    value YMIN.
!
  implicit none
!
  real slope
  real slope12
  real slope13
  real slope23
  real temp
  real x1
  real x2
  real x3
  real xmin
  real y1
  real y2
  real y3
  real ymin
!
!  Refuse to deal with coincident data.
!
  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MINABS - Fatal error!'
    write ( *, '(a)' ) '  X values are equal.'
    return
  end if
!
!  Sort the data.
!
  if ( x2 < x1 ) then
    call r_swap ( x1, x2 )
    call r_swap ( y1, y2 )
  end if

  if ( x3 < x1 ) then
    call r_swap ( x1, x3 )
    call r_swap ( y1, y3 )
  end if

  if ( x3 < x2 ) then
    call r_swap ( x2, x3 )
    call r_swap ( y2, y3 )
  end if
!
!  Now determine the slopes.
!
  slope12 = ( y2 - y1 ) / ( x2 - x1 )
  slope23 = ( y3 - y2 ) / ( x3 - x2 )
  slope13 = ( y3 - y1 ) / ( x3 - x1 )
!
!  Case 1: Minimum must be at an endpoint.
!
  if ( slope12 >= slope13 .or. slope12 >= 0.0E+00 ) then

    if ( y1 < y3 ) then
      xmin = x1
      ymin = y1
    else
      xmin = x3
      ymin = y3
    end if
!
!  Case 2: The curve decreases, and decreases faster than the line
!  joining the endpoints.
!
!  Whichever of SLOPE12 and SLOPE23 is the greater in magnitude
!  represents the actual slope of the underlying function.
!  Find where two lines of that slope, passing through the
!  endpoint data, intersect.
!
  else

    slope = max ( abs ( slope12 ), slope23 )

    xmin = 0.5E+00 * ( x1 + x3 + ( y1 - y3 ) / slope )
    ymin = y1 - slope * ( xmin - x1 )

  end if

  return
end
subroutine minquad ( x1, y1, x2, y2, x3, y3, xmin, ymin )
!
!*******************************************************************************
!
!! MINQUAD finds a local minimum of F(X) = A * X * X + B * X + C.
!
!
!  Note:
!
!    MINQUAD is primarily intended as a utility routine.
!    The square of the distance function between a point
!    and a line segment has the form of F(X).  Hence, we can seek
!    the line on the second segment which minimizes the square of
!    the distance to the other line segment.
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X1, Y1, X2, Y2, X3, Y3, are three sets of data
!    of the form ( X, F(X) ).  The three X values must be distinct.
!    On output, the data has been sorted so that X1 < X2 < X3,
!    and the Y values have been rearranged accordingly.
!
!    Output, real XMIN, YMIN.  XMIN is a point within the interval
!    spanned by X1, X2 and X3, at which F takes its local minimum
!    value YMIN.
!
  implicit none
!
  integer ierror
  real x
  real x1
  real x2
  real x3
  real xleft
  real xmin
  real xrite
  real y
  real y1
  real y2
  real y3
  real ymin
!
!  Refuse to deal with coincident data.
!
  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MINQUAD - Fatal error!'
    write ( *, '(a)' ) '  X values are equal.'
    return
  end if
!
!  Find the interval endpoints.
!
  xleft = min ( x1, x2, x3 )
  xrite = max ( x1, x2, x3 )
!
!  Find the minimizer and its function value, over the three input points.
!
  if ( y1 <= y2 .and. y1 <= y3 ) then
    xmin = x1
    ymin = y1
  else if ( y2 <= y1 .and. y2 <= y3 ) then
    xmin = x2
    ymin = y2
  else
    xmin = x3
    ymin = y3
  end if
!
!  Find the minimizer and its function value over the real line.
!
  call parabola_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )
!
!  If F is linear, then take the already computed min.
!
  if ( ierror == 2 ) then
!
!  If F has a maximum, then take the already computed min.
!
  else if ( y > ymin ) then
!
!  If the minimizer is to the left, take the already computed min.
!
  else if ( x < xleft ) then
!
!  If the minimizer is to the right, take the already computed min.
!
  else if ( x > xrite ) then

  else

    xmin = x
    ymin = y

  end if

  return
end
subroutine normal_01_sample ( x )
!
!*******************************************************************************
!
!! NORMAL_01_SAMPLE samples the standard Normal PDF.
!
!
!  Discussion:
!
!    The standard normal distribution has mean 0 and standard
!    deviation 1.
!
!  Method:
!
!    The Box-Muller method is used.
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real X, a sample of the PDF.
!
  implicit none
!
  integer, save :: iset = -1
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real v1
  real v2
  real x
  real, save :: xsave = 0.0E+00
!
  if ( iset == -1 ) then
    call random_seed ( )
    iset = 0
  end if

  if ( iset == 0 ) then

    call random_number ( harvest = v1 )

    if ( v1 <= 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  V1 <= 0.'
      write ( *, '(a,g14.6)' ) '  V1 = ', v1
      stop
    end if

    call random_number ( harvest = v2 )

    if ( v2 <= 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NORMAL_01_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  V2 <= 0.'
      write ( *, '(a,g14.6)' ) '  V2 = ', v2
      stop
    end if

    x = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * PI * v2 )

    xsave = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * PI * v2 )

    iset = 1

  else

    x = xsave
    iset = 0

  end if

  return
end
subroutine normal_01_vector ( n, x )
!
!*******************************************************************************
!
!! NORMAL_01_VECTOR samples the standard normal probability distribution.
!
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!  Method:
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Modified:
!
!    03 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values desired.  If N is negative,
!    then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Output, real X(N), a sample of the standard normal PDF.
!
  implicit none
!
  integer n
!
  integer m
  integer, save :: made = 0
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real r(n+1)
  integer, save :: saved = 0
  real x(n)
  integer xhi
  integer xlo
  real, save :: y = 0.0E+00
!
!  I'd like to allow the user to reset the random number seed.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0E+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  xlo = 1
  xhi = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    xlo = 2
  end if
!
!  If we don't need any more values, return.
!
  if ( xhi - xlo + 1 == 0 ) then

    return
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( xhi - xlo + 1 == 1 ) then

    call random_number ( harvest = r(1:2) )

    x(xhi) = sqrt ( -2.0E+00 * log ( r(1) ) ) * cos ( 2.0E+00 * pi * r(2) )
    y =      sqrt ( -2.0E+00 * log ( r(1) ) ) * sin ( 2.0E+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( xhi-xlo+1, 2 ) == 0 ) then

    m = ( xhi-xlo+1 ) / 2

    call random_number ( harvest = r(1:2*m) )

    x(xlo:xhi-1:2) = sqrt ( -2.0E+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0E+00 * pi * r(2:2*m:2) )

    x(xlo+1:xhi:2) = sqrt ( -2.0E+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0E+00 * pi * r(2:2*m:2) )

    made = made + xhi - xlo + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    xhi = xhi - 1

    m = ( xhi-xlo+1 ) / 2 + 1

    call random_number ( harvest = r(1:2*m) )

    x(xlo:xhi-1:2) = sqrt ( -2.0E+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0E+00 * pi * r(2:2*m-2:2) )

    x(xlo+1:xhi:2) = sqrt ( -2.0E+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0E+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0E+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0E+00 * pi * r(2*m) )

    y = sqrt ( -2.0E+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0E+00 * pi * r(2*m) )

    saved = 1

    made = made + xhi - xlo + 2

  end if

  return
end
subroutine octahedron_shape_3d ( max_num, max_order, point_num, face_num, &
  face_order, point_coord, face_point )
!
!*******************************************************************************
!
!! OCTAHEDRON_SHAPE_3D describes an octahedron in 3D.
!
!
!  Discussion:
!
!    The vertices of the octahedron lie on the unit sphere.
!
!  Modified:
!
!    11 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAX_NUM, the maximum number of faces, and of
!    points.  This quantity is used to dimension arrays.
!    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
!    not be set.
!
!    Input, integer MAX_ORDER, the maximum order of any face.
!    This quantity is used to dimension arrays.  If FACE_ORDER
!    exceeds MAX_ORDER, the arrays will not be set.
!
!    Output, integer POINT_NUM, the number of points in the shape.
!
!    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
!    the X, Y and Z coordinates of point J.
!
!    Output, integer FACE_NUM, the number of faces in the shape.
!
!    Output, integer FACE_ORDER, the number of vertices per face.
!
!    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
!    contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter-clockwise direction defined
!    by the outward normal at the face.
!
  implicit none
!
  integer max_num
  integer max_order
!
  integer face_num
  integer face_order
  integer face_point(max_order,max_num)
  integer point_num
  real point_coord(3,max_num)
!
  point_num = 6
  face_num = 8
  face_order = 3
!
!  Check.
!
  if ( point_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OCTAHEDRON_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of vertices exceeds MAX_NUM.'
    return
  end if

  if ( face_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OCTAHEDRON_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of faces exceeds MAX_NUM.'
    return
  end if

  if ( face_order > max_order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OCTAHEDRON_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Face order exceeds MAX_ORDER.'
    return
  end if
!
!  Set point coordinates.
!
  point_coord(1,1) =   0.0E+00
  point_coord(2,1) =   0.0E+00
  point_coord(3,1) = - 1.0E+00

  point_coord(1,2) =   0.0E+00
  point_coord(2,2) = - 1.0E+00
  point_coord(3,2) =   0.0E+00

  point_coord(1,3) =   1.0E+00
  point_coord(2,3) =   0.0E+00
  point_coord(3,3) =   0.0E+00

  point_coord(1,4) =   0.0E+00
  point_coord(2,4) =   1.0E+00
  point_coord(3,4) =   0.0E+00

  point_coord(1,5) =   0.0E+00
  point_coord(2,5) = - 1.0E+00
  point_coord(3,5) =   0.0E+00

  point_coord(1,6) =   0.0E+00
  point_coord(2,6) =   0.0E+00
  point_coord(3,6) =   1.0E+00
!
!  Set faces.
!
  face_point(1,1) = 1
  face_point(2,1) = 3
  face_point(3,1) = 2

  face_point(1,2) = 1
  face_point(2,2) = 4
  face_point(3,2) = 3

  face_point(1,3) = 1
  face_point(2,3) = 5
  face_point(3,3) = 4

  face_point(1,4) = 1
  face_point(2,4) = 2
  face_point(3,4) = 5

  face_point(1,5) = 2
  face_point(2,5) = 3
  face_point(3,5) = 6

  face_point(1,6) = 3
  face_point(2,6) = 4
  face_point(3,6) = 6

  face_point(1,7) = 4
  face_point(2,7) = 5
  face_point(3,7) = 6

  face_point(1,8) = 5
  face_point(2,8) = 2
  face_point(3,8) = 6

  return
end
subroutine para_contains_point_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )
!
!*******************************************************************************
!
!! PARA_CONTAINS_POINT_2D determines if a point is inside a parallelogram in 2D.
!
!
!  Diagram:
!
!         (X3,Y3).............
!            /              .
!           /              .
!          /              .
!    (X1,Y1)--------->(X2,Y2)
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of three
!    corners of the parallelogram.  (X1,Y1) should be directly connected
!    to (X2,Y2) and (X3,Y3).
!
!    Input, real X, Y, the point to be checked.
!
!    Output, logical INSIDE, is .TRUE. if (X,Y) is inside the
!    parallelogram, or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  real dot0_2d
  real enormsq0_2d
  logical inside
  real p21dot
  real p21normsq
  real p31dot
  real p31normsq
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  p21normsq = enormsq0_2d ( x1, y1, x2, y2 )
  p31normsq = enormsq0_2d ( x1, y1, x3, y3 )

  p21dot = dot0_2d ( x1, y1, x2, y2, x, y )
  p31dot = dot0_2d ( x1, y1, x3, y3, x, y )

  if ( 0.0E+00 <= p21dot .and. p21dot <= p21normsq .and. &
       0.0E+00 <= p31dot .and. p31dot <= p31normsq ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
subroutine para_contains_point_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x, y, z, inside )
!
!*******************************************************************************
!
!! PARA_CONTAINS_POINT_3D determines if a point is inside a parallelogram in 3D.
!
!
!  Diagram:
!
!         (X2,Y2,Z2).............
!            /                 .
!           /                 .
!          /                 .
!    (X1,Y1,Z1)--------->(X3,Y3,Z3)
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
!    the corners of the parallelogram.
!
!    Input, real X, Y, Z, the point to be checked.
!
!    Output, logical INSIDE, is .TRUE. if (X,Y,Z) is inside the
!    parallelogram, or on its boundary, and .FALSE. otherwise.
!    A slight amount of leeway is allowed for error, since a three
!    dimensional point may lie exactly in the plane of the parallelogram,
!    and yet be computationally slightly outside it.
!
  implicit none
!
  real dot
  real dotb
  real dott
  real enorm0_3d
  logical inside
  real, parameter :: tol = 0.00001E+00
  real v
  real x
  real x1
  real x2
  real x3
  real xn12
  real xn23
  real xn31
  real y
  real y1
  real y2
  real y3
  real yn12
  real yn23
  real yn31
  real z
  real z1
  real z2
  real z3
  real zn12
  real zn23
  real zn31
!
!  Compute V3, the vector normal to V1 = P2-P1 and V2 = P3-P1.
!
  call cross0_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xn12, yn12, zn12 )
!
!  If the component of V = P-P1 in the V3 direction is too large,
!  then it does not lie in the parallelogram.
!
  dot = ( x - x1 ) * xn12 + ( y - y1 ) * yn12 +  ( z - z1 ) * zn12

  v = enorm0_3d ( x, y, z, x2, y2, z2 )

  if ( abs ( dot ) > tol * ( 1.0E+00 + v ) ) then
    inside = .false.
    return
  end if
!
!  Compute V23, the vector normal to V2 and V3.
!
  call cross_3d ( x3-x1, y3-y1, z3-z1, xn12, yn12, zn12, xn23, yn23, zn23 )
!
!  Compute ALPHA = ( V dot V23 ) / ( V1 dot V23 )
!
  dott = ( x - x1 ) * xn23 + ( y - y1 ) * yn23 + ( z - z1 ) * zn23

  dotb = ( x2 - x1 ) * xn23 + ( y2 - y1 ) * yn23 + ( z2 - z1 ) * zn23

  if ( dotb < 0.0E+00 ) then
    dott = - dott
    dotb = - dotb
  end if

  if ( dott < 0.0E+00 .or. dott > dotb ) then
    inside = .false.
    return
  end if
!
!  Compute V31, the vector normal to V3 and V1.
!
  call cross_3d ( xn12, yn12, zn12, x2-x1, y2-y1, z2-z1, xn31, yn31, zn31 )
!
!  Compute BETA = ( V dot V31 ) / ( V2 dot V31 )
!
  dott = ( x - x1 ) * xn31 + ( y - y1 ) * yn31 + ( z - z1 ) * zn31

  dotb = ( x3 - x1 ) * xn31 + ( y3 - y1 ) * yn31 + ( z3 - z1 ) * zn31

  if ( dotb < 0.0E+00 ) then
    dott = - dott
    dotb = - dotb
  end if

  if ( dott < 0.0E+00 .or. dott > dotb ) then
    inside = .false.
    return
  end if
!
!  V = ALPHA * V1 + BETA * V2, where both ALPHA and BETA are between
!  0 and 1.
!
  inside = .true.

  return
end
subroutine para_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z, &
  dist )
!
!*******************************************************************************
!
!! PARA_POINT_DIST_3D: distance ( parallelogram, point ) in 3D.
!
!
!  Diagram:
!
!         (X2,Y2,Z2).............
!            /                 .
!           /                 .
!          /                 .
!    (X1,Y1,Z1)--------->(X3,Y3,Z3)
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, determine the
!    parallelogram, generated by the vectors from (X1,Y1) to (X2,Y2)
!    and from (X1,Y1) to (X3,Y3).
!
!    Input, real X, Y, Z, the point which is to be checked.
!
!    Output, real DIST, the distance from the point to the
!    parallelogram.  DIST is zero if the point lies exactly on the
!    parallelogram.
!
  implicit none
!
  real dis13
  real dis21
  real dis34
  real dis42
  real dist
  real enorm0_3d
  logical inside
  real t
  real temp
  real x
  real x1
  real x2
  real x3
  real x4
  real xn
  real xp
  real y
  real y1
  real y2
  real y3
  real y4
  real yn
  real yp
  real z
  real z1
  real z2
  real z3
  real z4
  real zn
  real zp
!
!  Compute P, the unit normal to X2-X1 and X3-X1:
!
  call cross0_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xp, yp, zp )

  temp = sqrt ( xp * xp + yp * yp + zp * zp )

  if ( temp == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARA_POINT_DIST_3D - Fatal error!'
    write ( *, '(a)' ) '  The normal vector is zero.'
    stop
  end if

  xp = xp / temp
  yp = yp / temp
  zp = zp / temp
!
!  Find ( XN, YN, ZN), the nearest point to ( X, Y, Z ) in the plane.
!
  t = xp * ( x - x1 ) + yp * ( y - y1 ) + zp * ( z - z1 )

  xn = x - xp * t
  yn = y - yp * t
  zn = z - zp * t
!
!  If ( XN, YN, ZN ) lies WITHIN the parallelogram, we're done.
!
  call para_contains_point_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
    x, y, z, inside )

  if ( inside ) then
    dist = enorm0_3d ( x, y, z, xn, yn, zn )
    return
  end if
!
!  Otherwise, find the distance between ( X, Y, Z ) and each of the
!  four line segments that make up the boundary of the parallelogram.
!
  x4 = x2 + x3 - x1
  y4 = y2 + y3 - y1
  z4 = z2 + z3 - z1

  call line_seg_point_dist_3d ( x1, y1, z1, x3, y3, z3, x, y, z, dis13 )

  call line_seg_point_dist_3d ( x3, y3, z3, x4, y4, z4, x, y, z, dis34 )

  call line_seg_point_dist_3d ( x4, y4, z4, x2, y2, z2, x, y, z, dis42 )

  call line_seg_point_dist_3d ( x2, y2, z2, x1, y1, z1, x, y, z, dis21 )

  dist = min ( dis13, dis34, dis42, dis21 )

  return
end
subroutine parabola_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )
!
!*******************************************************************************
!
!! PARABOLA_EX finds the extremal point of a parabola determined by three points.
!
!
!  Modified:
!
!    02 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
!    on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real X, Y, the X coordinate of the extremal point of the
!    parabola, and the value of the parabola at that point.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal point.
!
  implicit none
!
  real bot
  integer ierror
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  ierror = 0

  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    ierror = 1
    return
  end if

  if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
    x = x1
    y = y1
    return
  end if

  bot = ( x2 - x3 ) * y1 - ( x1 - x3 ) * y2 + ( x1 - x2 ) * y3

  if ( bot == 0.0E+00 ) then
    ierror = 2
    return
  end if

  x = 0.5E+00 * (    x1 * x1 * ( y3 - y2 ) &
               + x2 * x2 * ( y1 - y3 ) &
               + x3 * x3 * ( y2 - y1 ) ) / bot

  y = (    ( x - x2 ) * ( x - x3 ) * ( x2 - x3 ) * y1 &
         - ( x - x1 ) * ( x - x3 ) * ( x1 - x3 ) * y2 &
          + ( x - x1 ) * ( x - x2 ) * ( x1 - x2 ) * y3 ) / &
          ( ( x1 - x2 ) * ( x2 - x3 ) * ( x1 - x3 ) )

  return
end
subroutine parabola_ex2 ( x1, y1, x2, y2, x3, y3, x, y, a, b, c, ierror )
!
!*******************************************************************************
!
!! PARABOLA_EX2 finds the extremal point of a parabola determined by three points.
!
!
!  Modified:
!
!    29 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
!    on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real X, Y, the X coordinate of the extremal point of the
!    parabola, and the value of the parabola at that point.
!
!    Output, real A, B, C, the coefficients that define the parabola:
!    P(X) = A * X * X + B * X + C.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal
!    point.
!
  implicit none
!
  real a
  real b
  real c
  real det
  integer ierror
  real v(3,3)
  real w(3,3)
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  ierror = 0

  if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
    ierror = 1
    return
  end if

  if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
    x = x1
    y = y1
    return
  end if
!
!  Set up the Vandermonde matrix.
!
  v(1,1) = 1.0E+00
  v(1,2) = x1
  v(1,3) = x1 * x1

  v(2,1) = 1.0E+00
  v(2,2) = x2
  v(2,3) = x2 * x2

  v(3,1) = 1.0E+00
  v(3,2) = x3
  v(3,3) = x3 * x3
!
!  Get the inverse.
!
  call rmat3_inverse ( v, w, det )
!
!  Compute the parabolic coefficients.
!
  c = w(1,1) * y1 + w(1,2) * y2 + w(1,3) * y3
  b = w(2,1) * y1 + w(2,2) * y2 + w(2,3) * y3
  a = w(3,1) * y1 + w(3,2) * y2 + w(3,3) * y3
!
!  Determine the extremal point.
!
  if ( a == 0.0E+00 ) then
    ierror = 2
    return
  end if

  x = - b / ( 2.0E+00 * a )
  y = a * x * x + b * x + c

  return
end
subroutine parapp_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, x, y, z, dist )
!
!*******************************************************************************
!
!! PARAPP_POINT_DIST_3D: distance ( parallelepiped, point ) in 3D.
!
!
!  Definition:
!
!    A parallelepiped is a "slanted box", that is, opposite
!    sides are parallel planes.
!
!  Diagram:
!
!           7----8
!          /|   /|
!         / 3--/-5
!        4----6 /
!        |/   |/
!        1----2
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, half of
!    the corners of the box, from which the other corners can be
!    deduced.  The corners should be chosen so that the first corner
!    is directly connected to the other three.  The locations of
!    corners 5, 6, 7 and 8 will be computed by the parallelogram
!    relation.
!
!    Input, real X, Y, Z, the point which is to be checked.
!
!    Output, real DIST, the distance from the point to the box.  DIST is
!    zero if the point lies exactly on the box.
!
  implicit none
!
  real dis
  real dist
  real x
  real x1
  real x2
  real x3
  real x4
  real x5
  real x6
  real x7
  real x8
  real y
  real y1
  real y2
  real y3
  real y4
  real y5
  real y6
  real y7
  real y8
  real z
  real z1
  real z2
  real z3
  real z4
  real z5
  real z6
  real z7
  real z8
!
!  Fill in the other corners
!
  x5 = x2 + x3 - x1
  y5 = y2 + y3 - y1
  z5 = z2 + z3 - z1

  x6 = x2 + x4 - x1
  y6 = y2 + y4 - y1
  z6 = z2 + z4 - z1

  x7 = x3 + x4 - x1
  y7 = y3 + y4 - y1
  z7 = z3 + z4 - z1

  x8 = x2 + x3 + x4 - 2.0E+00 * x1
  y8 = y2 + y3 + y4 - 2.0E+00 * y1
  z8 = z2 + z3 + z4 - 2.0E+00 * z1
!
!  Compute the distance from the point ( X, Y, Z ) to each of the six
!  paralleogram faces.
!
  call para_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z, dis )

  dist = dis

  call para_point_dist_3d ( x1, y1, z1, x2, y2, z2, x4, y4, z4, x, y, z, dis )

  dist = min ( dist, dis )

  call para_point_dist_3d ( x1, y1, z1, x3, y3, z3, x4, y4, z4, x, y, z, dis )

  dist = min ( dist, dis )

  call para_point_dist_3d ( x8, y8, z8, x5, y5, z5, x6, y6, z6, x, y, z, dis )

  dist = min ( dist, dis )

  call para_point_dist_3d ( x8, y8, z8, x5, y5, z5, x7, y7, z7, x, y, z, dis )

  dist = min ( dist, dis )

  call para_point_dist_3d ( x8, y8, z8, x6, y6, z6, x7, y7, z7, x, y, z, dis )

  dist = min ( dist, dis )

  return
end
subroutine plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!*******************************************************************************
!
!! PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
!
!
!  Definition:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, X2, X3, Y3, Z3, are three points
!    on the plane, which must be distinct, and not collinear.
!
!    Output, real A, B, C, D, coefficients which describe the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3

  a = ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 )
  b = ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 )
  c = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )
  d = - x2 * a - y2 * b - z2 * c

  return
end
subroutine plane_exp2norm_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xp, yp, zp, &
  xn, yn, zn )
!
!*******************************************************************************
!
!! PLANE_EXP2NORM_3D converts an explicit plane to normal form in 3D.
!
!
!  Definition:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!    The normal form of a plane in 3D is
!
!      (Xp,Yp,Zp), a point on the plane, and
!      (Xn,Yn,Zn), the unit normal to the plane.
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, X2, X3, Y3, Z3, are three points
!    on the plane, which must be distinct, and not collinear.
!
!    Output, real XP, YP, ZP, a point on the plane.
!
!    Output, real XN, YN, ZN, the unit normal vector to the plane.
!
  implicit none
!
  real norm
  real x1
  real x2
  real x3
  real xn
  real xp
  real y1
  real y2
  real y3
  real yn
  real yp
  real z1
  real z2
  real z3
  real zn
  real zp

  xp = x1
  yp = y1
  zp = z1

  xn = ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 )
  yn = ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 )
  zn = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )

  norm = sqrt ( xn * xn + yn * yn + zn * zn )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_EXP2NORM_3D - Fatal error!'
    write ( *, '(a)' ) '  The normal vector is null.'
    write ( *, '(a)' ) '  Two points coincide, or nearly so.'
    stop
  end if

  xn = xn / norm
  yn = yn / norm
  zn = zn / norm

  return
end
subroutine plane_exp_normal_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xn, yn, &
  zn )
!
!*******************************************************************************
!
!! PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
!
!
!  Definition:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of three points that constitute a line.  These points should not
!    be identical, nor colinear.
!
!    Output, real XN, YN, ZN, the coordinates of the unit normal
!    vector to the plane containing the three points.
!
  implicit none
!
  real temp
  real x1
  real x2
  real x3
  real xn
  real y1
  real y2
  real y3
  real yn
  real z1
  real z2
  real z3
  real zn
!
!  The cross product (P2-P1) x (P3-P1) is a vector normal to
!  (P2-P1) and (P3-P1).
!
  call cross0_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xn, yn, zn )

  temp = sqrt ( xn * xn + yn * yn + zn * zn )

  if ( temp == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_EXP_NORMAL_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane is poorly defined.'
    stop
  else
    xn = xn / temp
    yn = yn / temp
    zn = zn / temp
  end if

  return
end
subroutine plane_exp_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x, y, z, dist )
!
!*******************************************************************************
!
!! PLANE_EXP_POINT_DIST_3D: distance ( explicit plane, point ) in 3D.
!
!
!  Definition:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!  Modified:
!
!    15 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
!    coordinates of three points on the plane.
!
!    Input, real X, Y, Z, the coordinates of the point.
!
!    Output, real DIST, the distance from the point to the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real dist
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
  real z
  real z1
  real z2
  real z3
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )

  call plane_imp_point_dist_3d ( a, b, c, d, x, y, z, dist )

  return
end
subroutine plane_exp_project_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  xf, yf, zf, npnt, xo, yo, zo, xp, yp, zp, ivis )
!
!*******************************************************************************
!
!! PLANE_EXP_PROJECT_3D projects points through a point onto a plane in 3D.
!
!
!  Definition:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
!    coordinates of three points on the plane.
!
!    Input, real XF, YF, ZF, are the coordinates of the focus point.
!
!    Input, integer NPNT, the number of points to project.
!
!    Input, real XO(NPNT), YO(NPNT), ZO(NPNT), are the coordinates of
!    the object points.
!
!    Output, real XP(NPNT), YP(NPNT), ZP(NPNT), are the coordinates of the
!    projections of the object points through the focus point onto
!    the plane.  XP, YP, and ZP may share the same memory as XO, YO,
!    and ZO, in which case the projections will overwrite the original data.
!
!    Output, integer IVIS(NPNT), visibility indicator:
!    3, the object was behind the plane;
!    2, the object was already on the plane;
!    1, the object was between the focus and the plane;
!    0, the line from the object to the focus is parallel to the plane,
!    so the object is "invisible".
!    -1, the focus is between the object and the plane.  The object
!    might be considered invisible.
!
  implicit none
!
  integer npnt
!
  real a
  real alpha
  real angle_rad_3d
  real b
  real beta
  real c
  real d
  real disfo
  real disfn
  integer i
  integer ivis(npnt)
  real x1
  real x2
  real x3
  real xf
  real xn
  real xo(npnt)
  real xp(npnt)
  real y1
  real y2
  real y3
  real yf
  real yn
  real yo(npnt)
  real yp(npnt)
  real z1
  real z2
  real z3
  real zf
  real zn
  real zo(npnt)
  real zp(npnt)
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!  Get the nearest point on the plane to the focus.
!
  call plane_imp_point_near_3d ( a, b, c, d, xf, yf, zf, xn, yn, zn )
!
!  Get the distance from the focus to the plane.
!
  call points_dist_3d ( xf, yf, zf, xn, yn, zn, disfn )
!
!  If the focus lies in the plane, this is bad.  We could still
!  project points that actually lie in the plane, but we'll
!  just bail out.
!
  if ( disfn == 0.0E+00 ) then
    ivis(1:npnt) = 0
    xp(1:npnt) = xf
    yp(1:npnt) = yf
    zp(1:npnt) = zf
    return
  end if
!
!  Process the points.
!
  do i = 1, npnt
!
!  Get the distance from the focus to the object.
!
    call points_dist_3d ( xf, yf, zf, xo(i), yo(i), zo(i), disfo )

    if ( disfo == 0.0E+00 ) then

      ivis(i) = 0
      xp(i) = xn
      yp(i) = yn
      zp(i) = zn

    else
!
!  Compute ALPHA, the angle between (OBJECT-FOCUS) and (NEAREST-FOCUS).
!
      alpha = angle_rad_3d ( xo(i), yo(i), zo(i), xf, yf, zf, xn, yn, zn )

      if ( cos ( alpha ) == 0.0E+00 ) then

        ivis(i) = 0
        xp(i) = xn
        yp(i) = yn
        zp(i) = zn

      else
!
!  BETA is Dist(NEAREST-FOCUS) / ( Cos(ALPHA)*Dist(OBJECT-FOCUS) )
!
        beta = disfn / ( cos ( alpha ) * disfo )

        if ( beta > 1.0E+00 ) then
          ivis(i) = 1
        else if ( beta == 1.0E+00 ) then
          ivis(i) = 2
        else if ( beta > 0.0E+00 ) then
          ivis(i) = 3
        else
          ivis(i) = -1
        end if
!
!  Set the projected point.
!
        xp(i) = xf + beta * ( xo(i) - xf )
        yp(i) = yf + beta * ( yo(i) - yf )
        zp(i) = zf + beta * ( zo(i) - zf )

      end if

    end if

  end do

  return
end
subroutine plane_grid_3d ( cor3, ierror, lines, maxcor3, maxline, &
  ncor3, nline, x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! PLANE_GRID_3D computes points and lines making up a planar grid in 3D.
!
!
!  Note:
!
!    The data format used is that of SGI Inventor.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real COR3(3,MAXCOR3), the X, Y, Z coordinates of points
!    used in the grid.
!
!    Output, integer IERROR, error indicator.
!    0, no error.
!    1, more space for point coordinates is needed.
!    2, more space for line data is needed.
!
!    Output, integer LINES(MAXLINE), the indices of points used in
!    the lines of the grid.  Successive entries of LINES are joined
!    by a line, unless an entry equals -1.  Note that indices begin
!    with 0.
!
!    Input, integer MAXCOR3, the maximum number of points that can be
!    stored.
!
!    Input, integer MAXLINE, the maximum number of line items that can
!    be stored.
!
!    Input/output, integer NCOR3, the number of points stored in COR3.
!
!    On input, if NCOR3 is zero, then the data computed by this routine
!    will be stored normally in COR3.  If NCOR3 is nonzero, then it
!    is assumed that COR3 already contains some useful data.  The
!    new data is appended to COR3.
!
!    On output, NCOR3 is increased by the number of points computed
!    by this routine.
!
!    Input/output, integer NLINE, the number of line data items.
!
!    On input, if NLINE is zero, then the data computed by this routine
!    will be stored normally in LINES.  If NLINE is nonzero, then it
!    is assumed that LINES already contains some useful data.  The
!    new data is appended to LINES.
!
!    On output, NLINE is increased by the number of points computed
!    by this routine.
!
!    Input, real X1, Y1, Z1, X2, Y2, X2, X3, Y3, Z3, are three points
!    on the plane, which must be distinct, and not collinear.
!
  implicit none
!
  integer maxcor3
  integer maxline
!
  real a
  real amax
  real amin
  real b
  real bmax
  real bmin
  real cor3(3,maxcor3)
  real dot
  integer i
  integer ierror
  integer j
  integer k
  integer lines(maxline)
  integer nbase
  integer ncor3
  integer nline
  integer, parameter :: nx = 5
  integer, parameter :: ny = 5
  real v1(3)
  real v2(3)
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  ierror = 0

  nbase = ncor3
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1) = x2 - x1
  v1(2) = y2 - y1
  v1(3) = z2 - z1

  call vector_unit_nd ( 3, v1 )

  v2(1) = x3 - x1
  v2(2) = y3 - y1
  v2(3) = z3 - z1

  dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

  v2(1:3) = v2(1:3) - dot * v1(1:3)

  call vector_unit_nd ( 3, v2 )
!
!  Compute the (V1,V2) coordinate range of the data.
!
  amax = 0.0E+00
  amin = 0.0E+00
  bmax = 0.0E+00
  bmin = 0.0E+00

  do i = 1, ncor3

    a = 0.0E+00
    b = 0.0E+00
    do j = 1, 3
      a = a + v1(j) * cor3(j,i)
      b = b + v2(j) * cor3(j,i)
    end do

    if ( i == 1 ) then
      amax = a
      amin = a
      bmax = b
      bmin = b
    else
      amax = max ( amax, a )
      amin = min ( amin, a )
      bmax = max ( bmax, b )
      bmin = min ( bmin, b )
    end if

  end do
!
!  Generate the points we will use.
!
  if ( ncor3 + nx * ny > maxcor3 ) then
    ierror = 1
    return
  end if

  do j = 1, ny

    b = ( real ( ny - j ) * bmin + real ( j - 1 ) * bmax ) / real ( ny - 1 )

    do i = 1, nx

      a = ( real ( nx - i ) * amin + real ( i - 1 ) * amax ) / real ( nx - 1 )

      ncor3 = ncor3 + 1
      cor3(1:3,ncor3) = a * v1(1:3) + b * v2(1:3)

    end do

  end do
!
!  Do the "horizontals".
!
  do i = 1, nx

    do j = 1, ny

      if ( nline >= maxline ) then
        ierror = 2
        return
      end if

      nline = nline + 1
      lines(nline) = nbase + ( j - 1 ) * nx + i

    end do

    if ( nline >= maxline ) then
      ierror = 2
      return
    end if

    nline = nline + 1
    lines(nline) = 0

  end do
!
!  Do the "verticals".
!
  do j = 1, ny

    do i = 1, nx

      if ( nline >= maxline ) then
        ierror = 2
        return
      end if

      nline = nline + 1
      lines(nline) = nbase + ( j - 1 ) * nx + i

    end do

    if ( nline >= maxline ) then
      ierror = 2
      return
    end if

    nline = nline + 1
    lines(nline) = 0

  end do

  return
end
subroutine plane_imp2exp_3d ( a, b, c, d, x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! PLANE_IMP2EXP_3D converts an implicit plane to explicit form in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0.
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!  Modified:
!
!    27 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, parameters that define the implicit plane.
!
!    Output, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
!    three points that lie on the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real x1
  real x2
  real x3
  real xn
  real xp
  real y1
  real y2
  real y3
  real yn
  real yp
  real z1
  real z2
  real z3
  real zn
  real zp
!
  call plane_imp2norm_3d ( a, b, c, d, xp, yp, zp, xn, yn, zn )

  call plane_norm2exp_3d ( xp, yp, zp, xn, yn, zn, x1, y1, z1, &
    x2, y2, z2, x3, y3, z3 )

  return
end
subroutine plane_imp2norm_3d ( a, b, c, d, xp, yp, zp, xn, yn, zn )
!
!*******************************************************************************
!
!! PLANE_IMP2NORM_3D converts an implicit plane to normal form in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0.
!
!    The normal form of a plane in 3D is
!
!      (Xp,Yp,Zp), a point on the plane, and
!      (Xn,Yn,Zn), the unit normal to the plane.
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, parameters that define the implicit plane.
!
!    Output, real XP, YP, ZP, a point on the plane.
!
!    Output, real XN, YN, ZN, the unit normal vector to the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real norm
  real xn
  real xp
  real yn
  real yp
  real zn
  real zp
!
  norm = sqrt ( a * a + b * b + c * c )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP2NORM_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane (A,B,C) has zero norm.'
    stop
  end if

  xn = a / norm
  yn = b / norm
  zn = c / norm

  if ( a /= 0.0E+00 ) then
    xp = - d / a
    yp = 0.0E+00
    zp = 0.0E+00
  else if ( b /= 0.0E+00 ) then
    xp = 0.0E+00
    yp = - d / b
    zp = 0.0E+00
  else if ( c /= 0.0E+00 ) then
    xp = 0.0E+00
    yp = 0.0E+00
    zp = - d / c
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP2NORM_3D - Fatal error!'
    write ( *, '(a)' ) '  The (A,B,C) vector is null.'
    stop
  end if

  return
end
function plane_imp_is_degenerate_3d ( a, b, c )
!
!*******************************************************************************
!
!! PLANE_IMP_IS_DEGENERATE_3D is TRUE if an implicit plane is degenerate.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    The implicit plane is degenerate if A = B = C = 0.
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, the coefficients of X, Y, and Z.
!
!    Output, logical PLANE_IMP_IS_DEGENERATE_3D, is TRUE if the plane 
!    is degenerate.
!
  implicit none
!
  real a
  real b
  real c
  logical plane_imp_is_degenerate_3d
!
  if ( a == 0.0E+00 .and. b == 0.0E+00 .and. c == 0.0E+00 ) then
    plane_imp_is_degenerate_3d = .true.
  else
    plane_imp_is_degenerate_3d = .false.
  end if

  return
end
subroutine plane_imp_line_par_int_3d ( a, b, c, d, x0, y0, z0, f, g, h, &
  intersect, x, y, z )
!
!*******************************************************************************
!
!! PLANE_IMP_LINE_PAR_INT_3D: intersection ( implicit plane, parametric line ) in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983, page 111.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, parameters that define the implicit plane.
!
!    Input, real X0, Y0, Z0, F, G, H, parameters that define the
!    parametric line.
!
!    Output, logical INTERSECT, is TRUE if the line and the plane
!    intersect, and false otherwise.
!
!    Output, real X, Y, Z, is a point of intersection of the line
!    and the plane, if INTERSECT is TRUE.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real denom
  real f
  real g
  real h
  logical intersect
  real norm1
  real norm2
  real t
  real, parameter :: tol = 0.00001E+00
  real x
  real x0
  real y
  real y0
  real z
  real z0
!
!  Check.
!
  norm1 = sqrt ( a * a + b * b + c * c )

  if ( norm1 == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    stop
  end if

  norm2 = sqrt ( f * f + g * g + h * h )

  if ( norm2 == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The line direction vector is null.'
    stop
  end if

  denom = a * f + b * g + c * h
!
!  The line and the plane may be parallel.
!
  if ( abs ( denom ) < TOL * norm1 * norm2 ) then

    if ( a * x0 + b * y0 + c * z0 + d == 0.0E+00 ) then
      intersect = .TRUE.
      x = x0
      y = y0
      z = z0
    else
      intersect = .FALSE.
      x = 0.0E+00
      y = 0.0E+00
      z = 0.0E+00
    end if
!
!  If they are not parallel, they must intersect.
!
  else

    intersect = .TRUE.
    t = - ( a * x0 + b * y0 + c * z0 + d ) / denom
    x = x0 + t * f
    y = y0 + t * g
    z = z0 + t * h

  end if

  return
end
subroutine plane_imp_line_seg_near_3d ( x1, y1, z1, x2, y2, z2, &
  a, b, c, d, dist, xp, yp, zp, xls, yls, zls )
!
!*******************************************************************************
!
!! PLANE_IMP_LINE_SEG_NEAR_3D: nearest ( implicit plane, line segment ) in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Modified:
!
!    17 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the endpoints of the line
!    segment.
!
!    Input, real A, B, C, D, the parameters that define the implicit
!    plane.
!
!    Output, real DIST, the distance between the line segment and
!    the plane.
!
!    Output, real XP, YP, ZP, the nearest point on the plane.
!
!    Output, real XLS, YLS, ZLS, the nearest point on the line segment
!    to the plane.  If DIST is zero, the (XLS,YLS,ZLS) is a point of
!    intersection of the plane and the line segment.
!
  implicit none
!
  real a
  real alpha
  real an
  real b
  real bn
  real c
  real cn
  real d
  real dist
  real dn
  real norm
  real p1
  real p2
  real x1
  real x2
  real xls
  real xp
  real y1
  real y2
  real yls
  real yp
  real z1
  real z2
  real zls
  real zp
!
  xls = 0.0E+00
  yls = 0.0E+00
  zls = 0.0E+00
  xp = 0.0E+00
  yp = 0.0E+00
  zp = 0.0E+00

  norm = sqrt ( a * a + b * b + c * c )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LINE_SEG_NEAR_3D - Fatal error!'
    write ( *, '(a)' ) '  Plane normal vector is null.'
    stop
  end if
!
!  The normalized coefficients allow us to compute the (signed) distance.
!
  an = a / norm
  bn = b / norm
  cn = c / norm
  dn = d / norm
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( x1 == x2 .and. y1 == y2 .and. z1 == z2 ) then

    p1 = an * x1 + bn * y1 + cn * z1 + dn
    dist = abs ( p1 )
    xls = x1
    yls = y1
    zls = z1
    xp = xls - an * p1
    yp = yls - bn * p1
    zp = zls - cn * p1
    return

  end if
!
!  Compute the projections of the two points onto the normal vector.
!
  p1 = an * x1 + bn * y1 + cn * z1 + dn
  p2 = an * x2 + bn * y2 + cn * z2 + dn
!
!  If these have the same sign, then the line segment does not
!  cross the plane, and one endpoint is the nearest point.
!
  if ( ( p1 > 0.0E+00 .and. p2 > 0.0E+00 ) .or. &
       ( p1 < 0.0E+00 .and. p2 < 0.0E+00 ) ) then

    p1 = abs ( p1 )
    p2 = abs ( p2 )

    if ( p1 < p2 ) then
      xls = x1
      yls = y1
      zls = z1
      xp = xls - an * p1
      yp = yls - bn * p1
      zp = zls - cn * p1
      dist = p1
    else
      xls = x2
      yls = y2
      zls = z2
      dist = p2
      xp = xls - an * p2
      yp = yls - bn * p2
      zp = zls - cn * p2
    end if
!
!  If the projections differ in sign, the line segment crosses the plane.
!
  else

    if ( p1 == 0.0E+00 ) then
      alpha = 0.0E+00
    else if ( p2 == 0.0E+00 ) then
      alpha = 1.0E+00
    else
      alpha = p2 / ( p2 - p1 )
    end if

    xls = alpha * x1 + ( 1.0E+00 - alpha ) * x2
    yls = alpha * y1 + ( 1.0E+00 - alpha ) * y2
    zls = alpha * z1 + ( 1.0E+00 - alpha ) * z2
    xp = xls
    yp = yls
    zp = zls

    dist = 0.0E+00

  end if

  return
end
subroutine plane_imp_point_dist_3d ( a, b, c, d, x, y, z, dist )
!
!*******************************************************************************
!
!! PLANE_IMP_POINT_DIST_3D: distance ( implicit plane, point ) in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    23 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, coefficients that define the plane as
!    the set of points for which A*X+B*Y+C*Z+D = 0.
!
!    Input, real X, Y, Z, the coordinates of the point.
!
!    Output, real DIST, the distance from the point to the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real dist
  real norm
  real x
  real y
  real z
!
  norm = sqrt ( a * a + b * b + c * c )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_POINT_DIST_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    stop
  end if

  dist = abs ( a * x + b * y + c * z + d ) / norm

  return
end
subroutine plane_imp_point_dist_signed_3d ( a, b, c, d, x, y, z, dist_signed )
!
!*******************************************************************************
!
!! PLANE_IMP_POINT_DIST_SIGNED_3D: signed distance ( implicit plane, point) in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, determine the equation of the
!    plane, which is:
!      A*X + B*Y + C*Z + D = 0.
!
!    Input, real X, Y, Z, the coordinates of the point.
!
!    Output, real DIST_SIGNED, the signed distance from the point to
!    the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real dist_signed
  real norm
  real x
  real y
  real z
!
  norm = sqrt ( a * a + b * b + c * c )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_POINT_DIST_SIGNED_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    stop
  end if

  dist_signed = - sign ( 1.0E+00, d ) * ( a * x + b * y + c * z + d ) / norm

  return
end
subroutine plane_imp_point_near_3d ( a, b, c, d, x, y, z, xn, yn, zn )
!
!*******************************************************************************
!
!! PLANE_IMP_POINT_NEAR_3D: nearest point on a implicit plane to a point in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, D, coefficients that define the plane as
!    the set of points for which A*X+B*Y+C*Z+D = 0.
!
!    Input, real X, Y, Z, the coordinates of the point.
!
!    Output, real XN, YN, ZN, the coordinates of the nearest point on
!    the plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  logical plane_imp_is_degenerate_3d
  real t
  real x
  real xn
  real y
  real yn
  real z
  real zn
!
  if ( plane_imp_is_degenerate_3d ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_POINT_NEAR_3D - Fatal error!'
    write ( *, '(a)' ) '  A = B = C = 0.'
    stop
  end if
!
!  The normal N to the plane is (A,B,C).
!
!  The line defined by (XN-X)/A = (YN-Y)/B = (ZN-Z)/C = T
!  goes through (X,Y,Z) and is parallel to N.
!
!  Solving for the point (XN,YN,ZN) we get
!
!    XN = A*T+X
!    YN = B*T+Y
!    ZN = C*T+Z
!
!  Now place these values in the equation for the plane:
!
!    A*(A*T+X) + B*(B*T+Y) + C*(C*T+Z) + D = 0
!
!  and solve for T:
!
!    T = (-A*X-B*Y-C*Z-D) / (A * A + B * B + C * C )
!
  t = - ( a * x + b * y + c * z + d ) / ( a * a + b * b + c * c )

  xn = x + a * t
  yn = y + b * t
  zn = z + c * t

  return
end
subroutine plane_imp_triangle_int_3d ( x1, y1, z1, x2, y2, z2, &
  x3, y3, z3, a, b, c, d, num_int, x, y, z )
!
!*******************************************************************************
!
!! PLANE_IMP_TRIANGLE_INT_3D: intersection ( implicit plane, triangle ) in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Discussion:
!
!    There may be 0, 1, 2 or 3 points of intersection returned.
!
!    If two intersection points are returned, then the entire line
!    between them comprises points of intersection.
!
!    If three intersection points are returned, then all points of
!    the triangle intersect the plane.
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of the vertices of the triangle.
!
!    Input, real A, B, C, D, the parameters that define the implicit plane.
!
!    Output, integer NUM_INT, the number of intersection points returned.
!
!    Output, real X(3), Y(3), Z(3), the coordinates of the intersection points.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real dist1
  real dist2
  real dist3
  integer num_int
  real x(3)
  real x1
  real x2
  real x3
  real y(3)
  real y1
  real y2
  real y3
  real z(3)
  real z1
  real z2
  real z3
!
  num_int = 0
!
!  Compute the signed distances between the vertices and the plane.
!
  dist1 = a * x1 + b * y1 + c * z1 + d
  dist2 = a * x2 + b * y2 + c * z2 + d
  dist3 = a * x3 + b * y3 + c * z3 + d
!
!  Consider any zero distances.
!
  if ( dist1 == 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x1
    y(num_int) = y1
    z(num_int) = z1

  end if

  if ( dist2 == 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x2
    y(num_int) = y2
    z(num_int) = z2

  end if

  if ( dist3 == 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x3
    y(num_int) = y3
    z(num_int) = z3

  end if
!
!  If 2 or 3 of the nodes intersect, we're already done.
!
  if ( num_int >= 2 ) then
    return
  end if
!
!  If one node intersects, then we're done unless the other two
!  are of opposite signs.
!
  if ( num_int == 1 ) then

    if ( dist1 == 0.0E+00 ) then

      call plane_triangle_int_add_3d ( x2, y2, z2, x3, y3, z3, &
        dist2, dist3, num_int, x, y, z )

    else if ( dist2 == 0.0E+00 ) then

      call plane_triangle_int_add_3d ( x1, y1, z1, x3, y3, z3, &
        dist1, dist3, num_int, x, y, z )

    else if ( dist3 == 0.0E+00 ) then

      call plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, &
        dist1, dist2, num_int, x, y, z )

    end if

    return

  end if
!
!  All nodal distances are nonzero, and there is at least one
!  positive and one negative.
!
  if ( dist1 * dist2 < 0.0E+00 .and. dist1 * dist3 < 0.0E+00 ) then

    call plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, &
      dist1, dist2, num_int, x, y, z )

    call plane_triangle_int_add_3d ( x1, y1, z1, x3, y3, z3, &
      dist1, dist3, num_int, x, y, z )

  else if ( dist2 * dist1 < 0.0E+00 .and.  dist2 * dist3 < 0.0E+00 ) then

    call plane_triangle_int_add_3d ( x2, y2, z2, x1, y1, z1, &
      dist2, dist1, num_int, x, y, z )

    call plane_triangle_int_add_3d ( x2, y2, z2, x3, y3, z3, &
      dist2, dist3, num_int, x, y, z )

  else if ( dist3 * dist1 < 0.0E+00 .and. dist3 * dist2 < 0.0E+00 ) then

    call plane_triangle_int_add_3d ( x3, y3, z3, x1, y1, z1, &
      dist3, dist1, num_int, x, y, z )

    call plane_triangle_int_add_3d ( x3, y3, z3, x2, y2, z2, &
      dist3, dist2, num_int, x, y, z )

  end if

  return
end
subroutine plane_imp_triangle_near_3d ( x1, y1, z1, x2, y2, z2, &
  x3, y3, z3, a, b, c, d, dist, num_near, x, y, z )
!
!*******************************************************************************
!
!! PLANE_IMP_TRIANGLE_NEAR_3D: nearest ( implicit plane, triangle ) in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Comments:
!
!    Please see to it that the underlying distance routine always returns
!    one of the endpoints if the entire line segment is at zero distance.
!
!  Modified:
!
!    18 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of the vertices of the triangle.
!
!    Input, real A, B, C, D, the parameters that define the implicit plane.
!
!    Output, real DIST, the distance between the triangle and the plane.
!
!    Output, integer NUM_NEAR, the number of nearest points returned.
!
!    Output, real X(6), Y(6), Z(6), a collection of nearest points.
!
!    If DIST = 0, then each point is a point of intersection, and there
!    will be at most 3 such points returned.
!
!    If DIST > 0, then the points are listed in pairs, with the first
!    being on the triangle, and the second on the plane.  Two points will
!    be listed in the most common case, but possibly 4 or 6.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real dist
  real dist12
  real dist23
  real dist31
  integer num_near
  real x(6)
  real x1
  real x2
  real x3
  real xp
  real xt
  real y(6)
  real y1
  real y2
  real y3
  real yp
  real yt
  real z(6)
  real z1
  real z2
  real z3
  real zp
  real zt
!
  num_near = 0
!
!  Consider the line segment P1 - P2.
!
  call plane_imp_line_seg_near_3d ( x1, y1, z1, x2, y2, z2, &
    a, b, c, d, dist12, xp, yp, zp, xt, yt, zt )

  dist = dist12

  num_near = num_near + 1
  x(num_near) = xt
  y(num_near) = yt
  z(num_near) = zt

  if ( dist12 > 0.0E+00 ) then
    num_near = num_near + 1
    x(num_near) = xp
    y(num_near) = yp
    z(num_near) = zp
  end if
!
!  Consider the line segment P2 - P3.
!
  call plane_imp_line_seg_near_3d ( x2, y2, z2, x3, y3, z3, &
    a, b, c, d, dist23, xp, yp, zp, xt, yt, zt )

  if ( dist23 < dist ) then

    num_near = 0
    dist = dist23

    num_near = num_near + 1
    x(num_near) = xt
    y(num_near) = yt
    z(num_near) = zt

    if ( dist23 > 0.0E+00 ) then
      num_near = num_near + 1
      x(num_near) = xp
      y(num_near) = yp
      z(num_near) = zp
    end if

  else if ( dist23 == dist ) then

    num_near = num_near + 1
    x(num_near) = xt
    y(num_near) = yt
    z(num_near) = zt

    if ( dist23 > 0.0E+00 ) then
      num_near = num_near + 1
      x(num_near) = xp
      y(num_near) = yp
      z(num_near) = zp
    end if

  end if
!
!  Consider the line segment P3 - P1.
!
  call plane_imp_line_seg_near_3d ( x3, y3, z3, x1, y1, z1, &
    a, b, c, d, dist31, xp, yp, zp, xt, yt, zt )

  if ( dist31 < dist ) then

    num_near = 0
    dist = dist31

    num_near = num_near + 1
    x(num_near) = xt
    y(num_near) = yt
    z(num_near) = zt

    if ( dist31 > 0.0E+00 ) then
      num_near = num_near + 1
      x(num_near) = xp
      y(num_near) = yp
      z(num_near) = zp
    end if

  else if ( dist31 == dist ) then

    num_near = num_near + 1
    x(num_near) = xt
    y(num_near) = yt
    z(num_near) = zt

    if ( dist31 > 0.0E+00 ) then
      num_near = num_near + 1
      x(num_near) = xp
      y(num_near) = yp
      z(num_near) = zp
    end if

  end if

  return
end
subroutine plane_norm2exp_3d ( xp, yp, zp, xn, yn, zn, x1, y1, z1, &
  x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! PLANE_NORM2EXP_3D converts a normal plane to explicit form in 3D.
!
!
!  Definition:
!
!    The normal form of a plane in 3D is
!
!      (Xp,Yp,Zp), a point on the plane, and
!      (Xn,Yn,Zn), the unit normal to the plane.
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!  Modified:
!
!    27 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XP, YP, ZP, a point on the plane.  (Actually, we never
!    need to know these values to do the calculation!)
!
!    Input, real XN, YN, ZN, a normal vector N to the plane.  The
!    vector must not have zero length, but it is not necessary for N
!    to have unit length.
!
!    Output, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
!    three points that lie on the plane.
!
  implicit none
!
  real x1
  real x2
  real x3
  real xn
  real xp
  real xq
  real xr
  real y1
  real y2
  real y3
  real yn
  real yp
  real yq
  real yr
  real z1
  real z2
  real z3
  real zn
  real zp
  real zq
  real zr

  call plane_norm_basis_3d ( xp, yp, zp, xn, yn, zn, xq, yq, zq, xr, yr, zr )

  x1 = xp
  y1 = yp
  z1 = zp

  x2 = xp + xq
  y2 = yp + yq
  z2 = zp + zq

  x3 = xp + xr
  y3 = yp + yr
  z3 = zp + zr

  return
end
subroutine plane_norm2imp_3d ( xp, yp, zp, xn, yn, zn, a, b, c, d )
!
!*******************************************************************************
!
!! PLANE_NORM2IMP_3D converts a normal form plane to implicit form in 3D.
!
!
!  Definition:
!
!    The normal form of a plane in 3D is
!
!      (Xp,Yp,Zp), a point on the plane, and
!      (Xn,Yn,Zn), the unit normal to the plane.
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0.
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XP, YP, ZP, a point on the plane.
!
!    Input, real XN, YN, ZN, the unit normal vector to the plane.
!
!    Output, real A, B, C, D, parameters that define the implicit plane.
!
  implicit none
!
  real a
  real b
  real c
  real d
  real xn
  real xp
  real yn
  real yp
  real zn
  real zp

  a = xn
  b = yn
  c = zn
  d = - a * xp - b * yp - c * zp

  return
end
subroutine plane_norm_basis_3d ( xp, yp, zp, xn, yn, zn, xq, yq, zq, &
  xr, yr, zr )
!
!*******************************************************************************
!
!! PLANE_NORM_BASIS_3D finds two perpendicular vectors in a plane in 3D.
!
!
!  Discussion:
!
!    Given a plane in point, normal form P = (XP,YP,ZP) and N = (XN,YN,ZN),
!    any point in that plane can be described in terms of the point P
!    plus a weighted sum of two vectors Q = (XQ,YQ,ZQ) and R = (XR,YR,ZR):
!
!      (X,Y,Z) = (XP,YP,ZP) + a * (XQ,YQ,ZQ) + b * (XR,YR,ZR).
!
!    The vector Q has unit length, and is perpendicular to P and R;
!    the vector R has unit length, and is perpendicular to P and Q.
!
!  Modified:
!
!    24 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XP, YP, ZP, a point on the plane.  (Actually, we never
!    need to know these values to do the calculation!)
!
!    Input, real XN, YN, ZN, a normal vector N to the plane.  The
!    vector must not have zero length, but it is not necessary for N
!    to have unit length.
!
!    Output, real XQ, YQ, ZQ, a vector of unit length, perpendicular
!    to the vector N and the vector R.
!
!    Output, real XR, YR, ZR, a vector of unit length, perpendicular
!    to the vector N and the vector Q.
!
  implicit none
!
  real dot
  real min_com
  real norm_n
  real norm_q
  real xn
  real xp
  real xq
  real xr
  real yn
  real yp
  real yq
  real yr
  real zn
  real zp
  real zq
  real zr
!
!  Compute the length of N = (XN,YN,ZN).
!
  norm_n = sqrt ( xn * xn + yn * yn + zn * zn )

  if ( norm_n == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_NORM_BASIS_3D - Fatal error!'
    write ( *, '(a)' ) '  The normal vector is 0.'
    stop
  end if
!
!  To find a vector distinct from N, find the minimum component
!  of N, and set the corresponding component of Q to 1, the
!  rest to zero.
!
  xq = 0.0E+00
  yq = 0.0E+00
  zq = 0.0E+00

  min_com = abs ( xn )

  if ( abs ( yn ) < min_com ) then
    min_com = abs ( yn )
  end if

  if ( abs ( zn ) < min_com ) then
    min_com = abs ( zn )
  end if

  if ( min_com == abs ( xn ) ) then
    xq = 1.0E+00

  else if ( min_com == abs ( yn ) ) then
    yq = 1.0E+00
  else if ( min_com == abs ( zn ) ) then
    zq = 1.0E+00
  end if
!
!  Now subtract off the component of Q in the direction of N,
!  computing Q = Q - Q dot N / || N ||,
!  and then normalize.
!
  dot = ( xq * xn + yq * yn + zq * zn ) / norm_n

  xq = xq - dot * xn / norm_n
  yq = yq - dot * yn / norm_n
  zq = zq - dot * zn / norm_n

  norm_q = sqrt ( xq * xq + yq * yq + zq * zq )

  xq = xq / norm_q
  yq = yq / norm_q
  zq = zq / norm_q
!
!  Now just take the cross product N x Q to get the R vector.
!  Plus, if we did things right, R will already have unit length.
!
  xr = ( yn * zq - zn * yq ) / norm_n
  yr = ( zn * xq - xn * zq ) / norm_n
  zr = ( xn * yq - yn * xq ) / norm_n

  return
end
subroutine plane_norm_triangle_int_3d ( x1, y1, z1, x2, y2, z2, &
  x3, y3, z3, xp, yp, zp, xn, yn, zn, num_int, x, y, z )
!
!*******************************************************************************
!
!! PLANE_NORM_TRIANGLE_INT_3D: intersection ( normal plane, triangle ) in 3D.
!
!
!  Definition:
!
!    The normal form of a plane in 3D is:
!
!      (Xp,Yp,Zp) is a point on the plane,
!      (Xn,Yn,Zn) is a normal vector to the plane.
!
!  Discussion:
!
!    There may be 0, 1, 2 or 3 points of intersection returned.
!
!    If two intersection points are returned, then the entire line
!    between them comprises points of intersection.
!
!    If three intersection points are returned, then all points of
!    the triangle intersect the plane.
!
!  Modified:
!
!    03 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of the vertices of the triangle.
!
!    Input, real XP, YP, ZP, a point on the plane.
!
!    Input, real XN, YN, ZN, a normal vector to the plane.
!
!    Output, integer NUM_INT, the number of intersection points returned.
!
!    Output, real X(3), Y(3), Z(3), the coordinates of the intersection points.
!
  implicit none
!
  real d
  real dist1
  real dist2
  real dist3
  integer num_int
  real x(3)
  real x1
  real x2
  real x3
  real xn
  real xp
  real y(3)
  real y1
  real y2
  real y3
  real yn
  real yp
  real z(3)
  real z1
  real z2
  real z3
  real zn
  real zp
!
  num_int = 0
!
!  Compute the signed distances between the vertices and the plane.
!
  d = - xn * xp - yn * yp - zn * zp
  dist1 = xn * x1 + yn * y1 + zn * z1 + d
  dist2 = xn * x2 + yn * y2 + zn * z2 + d
  dist3 = xn * x3 + yn * y3 + zn * z3 + d
!
!  Consider any zero distances.
!
  if ( dist1 == 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x1
    y(num_int) = y1
    z(num_int) = z1

  end if

  if ( dist2 == 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x2
    y(num_int) = y2
    z(num_int) = z2

  end if

  if ( dist3 == 0.0E+00 ) then

    num_int = num_int + 1
    x(num_int) = x3
    y(num_int) = y3
    z(num_int) = z3

  end if
!
!  If 2 or 3 of the nodes intersect, we're already done.
!
  if ( num_int >= 2 ) then
    return
  end if
!
!  If one node intersects, then we're done unless the other two
!  are of opposite signs.
!
  if ( num_int == 1 ) then

    if ( dist1 == 0.0E+00 ) then

      call plane_triangle_int_add_3d ( x2, y2, z2, x3, y3, z3, &
        dist2, dist3, num_int, x, y, z )

    else if ( dist2 == 0.0E+00 ) then

      call plane_triangle_int_add_3d ( x1, y1, z1, x3, y3, z3, &
        dist1, dist3, num_int, x, y, z )

    else if ( dist3 == 0.0E+00 ) then

      call plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, &
        dist1, dist2, num_int, x, y, z )

    end if

    return

  end if
!
!  All nodal distances are nonzero, and there is at least one
!  positive and one negative.
!
  if ( dist1 * dist2 < 0.0E+00 .and. dist1 * dist3 < 0.0E+00 ) then

    call plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, &
      dist1, dist2, num_int, x, y, z )

    call plane_triangle_int_add_3d ( x1, y1, z1, x3, y3, z3, &
      dist1, dist3, num_int, x, y, z )

  else if ( dist2 * dist1 < 0.0E+00 .and. dist2 * dist3 < 0.0E+00 ) then

    call plane_triangle_int_add_3d ( x2, y2, z2, x1, y1, z1, &
      dist2, dist1, num_int, x, y, z )

    call plane_triangle_int_add_3d ( x2, y2, z2, x3, y3, z3, &
      dist2, dist3, num_int, x, y, z )

  else if ( dist3 * dist1 < 0.0E+00 .and. dist3 * dist2 < 0.0E+00 ) then

    call plane_triangle_int_add_3d ( x3, y3, z3, x1, y1, z1, &
      dist3, dist1, num_int, x, y, z )

    call plane_triangle_int_add_3d ( x3, y3, z3, x2, y2, z2, &
      dist3, dist2, num_int, x, y, z )

  end if

  return
end
subroutine plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, dist1, dist2, &
  num_int, x, y, z )
!
!*******************************************************************************
!
!! PLANE_TRIANGLE_INT_ADD_3D is a utility for plane/triangle intersections.
!
!
!  Discussion:
!
!    This routine is called to consider the value of the signed distance
!    from a plane of two nodes of a triangle.  If the two values
!    have opposite signs, then there is a point of intersection between
!    them.  The routine computes this point and adds it to the list.
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
!    of two vertices of a triangle.
!
!    Input, real DIST1, DIST2, the signed distances of the two vertices
!    from a plane.
!
!    Input/output, integer NUM_INT, the number of intersection points.
!
!    Input/output, real X(NUM_INT), Y(NUM_INT), Z(NUM_INT), the coordinates
!    of the intersection points.
!
  implicit none
!
  real alpha
  real dist1
  real dist2
  integer num_int
  real x(*)
  real x1
  real x2
  real y(*)
  real y1
  real y2
  real z(*)
  real z1
  real z2
!
  if ( dist1 == 0.0E+00 ) then
    num_int = num_int + 1
    x(num_int) = x1
    y(num_int) = y1
    z(num_int) = z1
  else if ( dist2 == 0.0E+00 ) then
    num_int = num_int + 1
    x(num_int) = x2
    y(num_int) = y2
    z(num_int) = z2
  else if ( dist1 * dist2 < 0.0E+00 ) then
    alpha = dist2 / ( dist2 - dist1 )
    num_int = num_int + 1
    x(num_int) = alpha * x1 + ( 1.0E+00 - alpha ) * x2
    y(num_int) = alpha * y1 + ( 1.0E+00 - alpha ) * y2
    z(num_int) = alpha * z1 + ( 1.0E+00 - alpha ) * z2
  end if

  return
end
subroutine planes_imp_angle_3d ( a1, b1, c1, d1, a2, b2, c2, d2, angle )
!
!*******************************************************************************
!
!! PLANES_IMP_ANGLE_3D: dihedral angle between implicit planes in 3D.
!
!
!  Definition:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    If two planes P1 and P2 intersect in a nondegenerate way, then there is a
!    line of intersection L0.  Consider any plane perpendicular to L0.  The
!    dihedral angle of P1 and P2 is the angle between the lines L1 and L2, where
!    L1 is the intersection of P1 and P0, and L2 is the intersection of P2 
!    and P0.
!
!    The dihedral angle may also be calculated as the angle between the normal
!    vectors of the two planes.  Note that if the planes are parallel or
!    coincide, the normal vectors are identical, and the dihedral angle is 0.
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Math Tables and Formulae, 30th edition,
!    Section 4.13, "Planes",
!    CRC Press, 1996, pages 305-306.
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A1, B1, C1, D1, coefficients that define the first plane.
!
!    Input, real A2, B2, C2, D2, coefficients that define the second plane.
!
!    Output, real ANGLE, the dihedral angle, in radians, defined by the two
!    planes.  If either plane is degenerate, or they do not intersect, or
!    they coincide, then ANGLE is set to HUGE(1.0).
!    Otherwise, ANGLE is between 0 and PI.
!
  implicit none
!
  real a1
  real a2
  real angle
  real b1
  real b2
  real c1
  real c2
  real cosine
  real d1
  real d2
  real norm1
  real norm2
!
  angle = huge ( angle )

  norm1 = sqrt ( a1**2 + b1**2 + c1**2 )

  if ( norm1 == 0.0E+00 ) then
    return
  end if

  norm2 = sqrt ( a2**2 + b2**2 + c2**2 )

  if ( norm2 == 0.0E+00 ) then
    return
  end if

  cosine = ( a1 * a2 + b1 * b2 + c1 * c2 ) / ( norm1 * norm2 )

  angle = acos ( angle )

  return
end
function point_inside_box_2d ( x1, y1, x2, y2, x, y )
!
!*******************************************************************************
!
!! POINT_INSIDE_BOX_2D determines if a point is inside a box in 2D.
!
!
!  Definition:
!
!    A "box" is defined by its "left down" corner and its
!    "right up" corner, and all the points between.  It is
!    assumed that the sides of the box align with coordinate directions.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the two corners of the box.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, logical POINT_INSIDE_BOX_2D, is .TRUE. if (X,Y) is inside the
!    box, or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  logical point_inside_box_2d
  real x
  real x1
  real x2
  real y
  real y1
  real y2
!
  if ( x1 <= x .and. x <= x2 .and. &
       y1 <= y .and. y <= y2 ) then
    point_inside_box_2d = .true.
  else
    point_inside_box_2d = .false.
  end if

  return
end
function point_inside_box_3d ( x1, y1, z1, x2, y2, z2, x, y, z )
!
!*******************************************************************************
!
!! POINT_INSIDE_BOX_3D determines if a point is inside a box in 3D.
!
!
!  Definition:
!
!    A "box" is defined by its "left down front" corner and its
!    "right up back" corner, and all the points between.  It is
!    assumed that the sides of the box align with coordinate directions.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the two corners of the box.
!
!    Input, real X, Y, Z, the point to be checked.
!
!    Output, logical POINT_INSIDE_BOX_3D, is .TRUE. if (X,Y,Z) is inside the
!    box, or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  logical point_inside_box_3d
  real x
  real x1
  real x2
  real y
  real y1
  real y2
  real z
  real z1
  real z2
!
  if ( x1 <= x .and. x <= x2 .and. &
       y1 <= y .and. y <= y2 .and. &
       z1 <= z .and. z <= z2 ) then
    point_inside_box_3d = .true.
  else
    point_inside_box_3d = .false.
  end if

  return
end
function point_inside_parallelipiped_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, x, y, z )
!
!*******************************************************************************
!
!! POINT_INSIDE_PARALLELIPIPED_3D determines if a point is inside a parallelepiped in 3D.
!
!
!  Definition:
!
!    A parallelepiped is a "slanted box", that is, opposite
!    sides are parallel planes.
!
!  Diagram:
!
!                   *------------------*
!                  / \                / \
!                 /   \              /   \
!                /     \            /     \
!          (X4,Y4,Z4)--------------*       \
!                \        .         \       \
!                 \        .         \       \
!                  \        .         \       \
!                   \   (X2,Y2,Z2).....\-------\
!                    \     /            \     /
!                     \   /              \   /
!                      \ /                \ /
!                (X1,Y1,Z1)-----------(X3,Y3,Z3)
!
!  Modified:
!
!    04 February 199
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
!    coordinates of four corners of the parallelepiped, which we will
!    call P1, P2, P3 and P4.  It is assumed that P2, P3 and P4 are
!    immediate neighbors of P1.
!
!    Input, real X, Y, Z, the point to be checked.
!
!    Output, logical POINT_INSIDE_PARALLELIPIPED_3D, is .TRUE. if (X,Y,Z) 
!    is inside the parallelepiped, or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  real dot0_3d
  real enormsq0_3d
  real p21dot
  real p21normsq
  real p31dot
  real p31normsq
  real p41dot
  real p41normsq
  logical point_inside_parallelipiped_3d
  real x
  real x1
  real x2
  real x3
  real x4
  real y
  real y1
  real y2
  real y3
  real y4
  real z
  real z1
  real z2
  real z3
  real z4
!
  p21normsq = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )
  p31normsq = enormsq0_3d ( x1, y1, z1, x3, y3, z3 )
  p41normsq = enormsq0_3d ( x1, y1, z1, x4, y4, z4 )

  p21dot = dot0_3d ( x1, y1, z1, x2, y2, z2, x, y, z )
  p31dot = dot0_3d ( x1, y1, z1, x3, y3, z3, x, y, z )
  p41dot = dot0_3d ( x1, y1, z1, x4, y4, z4, x, y, z )

  if ( 0.0E+00 <= p21dot .and. p21dot <= p21normsq .and. &
       0.0E+00 <= p31dot .and. p31dot <= p31normsq .and. &
       0.0E+00 <= p41dot .and. p41dot <= p41normsq ) then
    point_inside_parallelipiped_3d = .true.
  else
    point_inside_parallelipiped_3d = .false.
  end if

  return
end
function points_avoid_point_naive_2d ( n, xy_set, xy_test )
!
!*******************************************************************************
!
!! POINTS_AVOID_POINT_NAIVE_2D determines if a point is "far enough" from a set of points in 2D.
!
!
!  Discussion:
!
!    The routine discards points that are too close to other points.
!    The method used to check this is quadratic in the number of points,
!    and may take an inordinate amount of time if there are a large
!    number of points.  But in that case, what do you want?  If you want
!    lots of points, you don't want to delete any because it won't matter.
!
!    The test point is "far enough" from an accepted point if
!    the Euclidean distance is at least 100 times EPSILON.
!
!  Modified:
!
!    24 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of accepted points.
!
!    Input, real XY_SET(2,N), the accepted points.
!
!    Input, real XY_TEST(2), a point to be tested.
!
!    Output, logical POINTS_AVOID_POINT_NAIVE_2D, is TRUE if XY_TEST is
!    "far enough" from all the accepted points.
!
  implicit none
!
  integer n
  integer, parameter :: ndim = 2
!
  integer j
  logical points_avoid_point_naive_2d
  real tol
  real xy_set(ndim,n)
  real xy_test(ndim)
!
  tol = 100.0E+00 * epsilon ( tol )

  points_avoid_point_naive_2d = .true.

  do j = 1, n

    if ( sqrt ( sum ( ( xy_set(1:ndim,j) - xy_test(1:ndim) )**2 ) ) < tol ) then
      points_avoid_point_naive_2d = .false.
      return
    end if

  end do

  return
end
subroutine points_bisect_line_imp_2d ( x1, y1, x2, y2, a, b, c )
!
!*******************************************************************************
!
!! POINTS_BISECT_LINE_IMP_2D finds the implicit line bisecting the line between two points in 2D.
!
!
!  Formula:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the coordinates of two points.
!
!    Output, real A, B, C, the parameters of the implicit line
!    equidistant from both points.
!
  implicit none
!
  real a
  real b
  real c
  real x1
  real x2
  real y1
  real y2
!
  a = x1 - x2
  b = y1 - y2
  c = - 0.5E+00 * ( ( x1 * x1 + y1 * y1 ) - ( x2 * x2 + y2 * y2 ) )

  return
end
subroutine points_bisect_line_par_2d ( x1, y1, x2, y2, f, g, x, y )
!
!*******************************************************************************
!
!! POINTS_BISECT_LINE_PAR_2D finds the parametric line bisecting the line between two points in 2D.
!
!
!  Formula:
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, the coordinates of two points.
!
!    Output, real F, G, X, Y, the parameters of the parametric line
!    equidistant from both points.
!
  implicit none
!
  real f
  real g
  real x
  real x1
  real x2
  real y
  real y1
  real y2
!
  f = 0.5E+00 * ( x1 + x2 )
  g = 0.5E+00 * ( y1 + y2 )
  x = - ( y2 - y1 )
  y = + ( x2 - x1 )

  return
end
subroutine points_centroid_2d ( n, x, y, cent )
!
!*******************************************************************************
!
!! POINTS_CENTROID_2D computes the discrete centroid of a point set in 2D.
!
!
!  Formula:
!
!    Given a discrete set of points S, the discrete centroid z is defined by
!
!                           Sum ( x in S ) ( x - z )**2
!        = min ( y in S ) { Sum ( x in S ) ( x - y )**2
!
!    In other words, the discrete centroid is a point in the set whose distance
!    to the other points is minimized.  The discrete centroid of a point set
!    need not be unique.  Consider a point set that comprises the
!    vertices of an equilateral triangle.
!
!    This discrete centroid may also be referred to as the K-means cluster.
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real X(N), Y(N), the coordinates of the points.
!
!    Output, integer CENT, the index of a discrete centroid of the set.
!
  implicit none
!
  integer n
!
  integer cent
  real dist
  real dist_min
  integer i
  integer j
  real x(n)
  real y(n)
!
  dist_min = 0.0E+00
  cent = 0

  do i = 1, n

    dist = 0.0E+00
    do j = 1, n
      dist = dist + ( x(i) - x(j) )**2 + ( y(i) - y(j) )**2
    end do

    if ( i == 1 ) then
      dist_min = dist
      cent = i
    else if ( dist < dist_min ) then
      dist_min = dist
      cent = i
    end if

  end do

  return
end
subroutine points_colin_2d ( x1, y1, x2, y2, x3, y3, colin )
!
!*******************************************************************************
!
!! POINTS_COLIN_2D estimates the colinearity of 3 points in 2D.
!
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of the points.
!
!    Output, real COLIN, an estimate of colinearity, namely, the ratio
!    of the area of the triangle spanned by the points to the area
!    of the equilateral triangle with the same perimeter.
!    COLIN is 1 if the points are maximally noncolinear, 0 if the
!    points are exactly colinear, and otherwise is closer to 1 or 0 depending
!    on whether the points are far or close to colinearity.
!
  implicit none
!
  real area
  real area2
  real colin
  real enorm0_2d
  real perim
  real side
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
!
  call triangle_area_2d ( x1, y1, x2, y2, x3, y3, area )

  if ( area == 0.0E+00 ) then

    colin = 0.0E+00

  else

    perim =   enorm0_2d ( x1, y1, x2, y2 ) &
            + enorm0_2d ( x2, y2, x3, y3 ) &
            + enorm0_2d ( x3, y3, x1, y1 )

    side = perim / 3.0E+00

    area2 = 0.25E+00 * sqrt ( 3.0E+00 ) * side * side

    colin = area / area2

  end if

  return
end
subroutine points_colin_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, colin )
!
!*******************************************************************************
!
!! POINTS_COLIN_3D estimates the colinearity of 3 points in 3D.
!
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
!    the points.
!
!    Output, real COLIN, an estimate of colinearity, namely, the ratio
!    of the area of the triangle spanned by the points to the area
!    of the equilateral triangle with the same perimeter.
!    COLIN is 1 if the points are maximally noncolinear, 0 if the
!    points are exactly colinear, and otherwise is closer to 1 or 0 depending
!    on whether the points are far or close to colinearity.
!
  implicit none
!
  real area
  real area2
  real colin
  real enorm0_3d
  real perim
  real side
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  call triangle_area_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, area )

  area = abs ( area )

  if ( area == 0.0E+00 ) then

    colin = 0.0E+00

  else

    perim = enorm0_3d ( x1, y1, z1, x2, y2, z2 ) &
          + enorm0_3d ( x2, y2, z2, x3, y3, z3 ) &
          + enorm0_3d ( x3, y3, z3, x1, y1, z1 )

    side = perim / 3.0E+00

    area2 = 0.25E+00 * sqrt ( 3.0E+00 ) * side * side

    colin = area / area2

  end if

  return
end
subroutine points_delaunay_naive_2d ( n, x, y, maxtri, ntri, tri )
!
!*******************************************************************************
!
!! POINTS_DELAUNAY_NAIVE_2D is a naive Delaunay triangulation scheme.
!
!
!  Discussion:
!
!    This routine is only suitable as a demonstration code for small
!    problems.  Its running time is of order N**4.  Much faster
!    algorithms are available.
!
!    Given a set of nodes in the plane, a triangulation is set of
!    triples of distinct nodes, forming triangles, so that every
!    point within the convex hull of the set of nodes is either
!    one of the nodes, or lies on an edge of one or more triangles,
!    or lies within exactly one triangle.
!
!    A Delaunay triangulation is a triangulation with additional
!    properties.
!
!  Reference:
!
!    Joseph O'Rourke,
!    Computational Geometry,
!    Cambridge University Press,
!    Second Edition, 1998, page 187.
!
!  Modified:
!
!    08 November 2000
!
!  Parameters:
!
!    Input, integer N, the number of nodes.  N must be at least 3.
!
!    Input, real X(N), Y(N), the coordinates of the nodes.
!
!    Input, integer MAXTRI, the maximum number of triangles.
!
!    Output, integer NTRI, the number of triangles in the triangulation.
!
!    Output, integer TRI(3,MAXTRI), contains in TRI(1,I), TRI(2,I)
!    and TRI(3,I) the indices of the nodes making the I-th triangle.
!
  implicit none
!
  integer maxtri
  integer n
!
  logical flag
  integer i
  integer j
  integer k
  integer m
  integer ntri
  integer tri(3,maxtri)
  real x(n)
  real xn
  real y(n)
  real yn
  real z(n)
  real zn
!
  ntri = 0

  if ( n < 3 ) then
    return
  end if
!
!  Compute Z = X*X + Y*Y.
!
  z(1:n) = x(1:n)**2 + y(1:n)**2
!
!  For each triple (I,J,K):
!
  do i = 1, n-2
    do j = i+1, n
      do k = i+1, n

        if ( j /= k ) then

          xn = ( y(j) - y(i) ) * ( z(k) - z(i) ) &
             - ( y(k) - y(i) ) * ( z(j) - z(i) )

          yn = ( x(k) - x(i) ) * ( z(j) - z(i) ) &
             - ( x(j) - x(i) ) * ( z(k) - z(i) )

          zn = ( x(j) - x(i) ) * ( y(k) - y(i) ) &
             - ( x(k) - x(i) ) * ( y(j) - y(i) )

          flag = ( zn < 0.0E+00 )

          if ( flag ) then
            do m = 1, n
              flag = flag .and. &
                ( ( x(m) - x(i) ) * xn + ( y(m) - y(i) ) * yn &
                  + ( z(m) - z(i) ) * zn <= 0.0E+00 )
            end do
          end if

          if ( flag ) then
            if ( ntri < maxtri ) then
              ntri = ntri + 1
              tri(1,ntri) = i
              tri(2,ntri) = j
              tri(3,ntri) = k
            end if
          end if

        end if

      end do
    end do
  end do

  return
end
subroutine points_dist_2d ( x1, y1, x2, y2, dist )
!
!*******************************************************************************
!
!! POINTS_DIST_2D finds the distance between two points in 2D.
!
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, determines the pair of
!    points (X1,Y1) and (X2,Y2) whose distance apart is be determined.
!
!    Output, real DIST, the distance between the points.
!
  implicit none
!
  real dist
  real x1
  real x2
  real y1
  real y2
!
  dist = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 )

  return
end
subroutine points_dist_3d ( x1, y1, z1, x2, y2, z2, dist )
!
!*******************************************************************************
!
!! POINTS_DIST_3D finds the distance between two points in 3D.
!
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, determines the pair of points
!    (X1,Y1,Z1) and (X2,Y2,Z2) whose distance apart is be determined.
!
!    Output, real DIST, the distance between the points.
!
  implicit none
!
  real dist
  real x1
  real x2
  real y1
  real y2
  real z1
  real z2
!
  dist = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 + ( z1 - z2 )**2 )

  return
end
subroutine points_dist_nd ( x1, x2, n, dist )
!
!*******************************************************************************
!
!! POINTS_DIST_ND finds the distance between two points in ND.
!
!
!  Modified:
!
!    31 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1(N), X2(N), the coordinates of two points.
!
!    Input, integer N, the dimension of the space.
!
!    Output, real DIST, the distance between the points.
!
  implicit none
!
  integer n
!
  real dist
  real x1(n)
  real x2(n)
!
  dist = sqrt ( sum ( ( x1(1:n) - x2(1:n) )**2 ) )

  return
end
subroutine points_dist_sphere ( lat1, long1, lat2, long2, radius, dist )
!
!*******************************************************************************
!
!! POINTS_DIST_SPHERE finds the distance between two points on a sphere.
!
!
!  Discussion:
!
!    The distance is measured on the surface of the sphere.
!
!  Modified:
!
!    05 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real LAT1, LONG1, LAT2, LONG2, the latitude and longitude of
!    the two points, in radians.
!
!    Input, real RADIUS, the radius of the sphere.
!
!    Output, real DIST, the distance between the points.
!
  implicit none
!
  real dist
  real lat1
  real lat2
  real long1
  real long2
  real radius
  real theta
!
  theta = acos ( sin ( lat1 ) * sin ( lat2 ) &
               + cos ( lat1 ) * cos ( lat2 ) * cos ( long1 - long2 ) )

  dist = radius * theta

  return
end
subroutine points_hull_2d ( ival, n, nval, x, y )
!
!*******************************************************************************
!
!! POINTS_HULL_2D computes the convex hull of a set of points in 2D.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IVAL(N).  Entries 1 through NVAL of IVAL contain
!    the indices of the nodes that form the convex hull, in order.
!
!    Input, integer N, the number of nodes.
!
!    Output, integer NVAL, the number of nodes that lie on the convex hull.
!
!    Input, real X(N), Y(N), the X and Y coordinates of the nodes.
!
  implicit none
!
  integer n
!
  real ai
  real angle_deg_2d
  real angmax
  real di
  real dr
  real enorm0_2d
  integer i
  integer iq
  integer ir
  integer istart
  integer ival(n)
  integer nval
  real x(n)
  real xp
  real xq
  real xr
  real y(n)
  real yp
  real yq
  real yr
!
  if ( n < 1 ) then
    nval = 0
    return
  end if
!
!  If N = 1, the hull is the point.
!
  if ( n == 1 ) then
    nval = 1
    ival(1) = 1
    return
  end if
!
!  If N = 2, then the convex hull is either the two distinct points,
!  or possibly a single (repeated) point.
!
  if ( n == 2 ) then

    if ( x(1) /= x(2) .or. y(1) /= y(2) ) then
      nval = 2
      ival(1) = 1
      ival(2) = 2
    else
      nval = 1
      ival(1) = 1
    end if

    return

  end if
!
!  Find the leftmost point, and take the bottom-most in a tie.
!  Call it "Q".
!
  iq = 1
  do i = 2, n
    if ( x(i) < x(iq) .or. ( x(i) == x(iq) .and. y(i) < y(iq) ) ) then
      iq = i
    end if
  end do

  xq = x(iq)
  yq = y(iq)
!
!  Remember the starting point.
!
  istart = iq
  nval = 1
  ival(1) = iq
!
!  For the first point, make a dummy previous point, 1 unit south,
!  and call it "P".
!
  xp = xq
  yp = yq - 1.0E+00
!
!  Now, having old point P, and current point Q, find the new point R
!  so the angle PQR is maximal.
!
!  Watch out for the possibility that the two nodes are identical.
!
  do

    ir = 0
    angmax = 0.0E+00

    do i = 1, n

      if ( i /= iq .and. ( x(i) /= xq .or. y(i) /= yq ) ) then

        ai = angle_deg_2d ( xp, yp, xq, yq, x(i), y(i) )

        if ( ir == 0 .or. ai > angmax ) then

          ir = i
          xr = x(ir)
          yr = y(ir)
          angmax = ai
!
!  In case of ties, choose the nearer point.
!
        else if ( ir /= 0 .and. ai == angmax ) then

          di = enorm0_2d ( xq, yq, x(i), y(i) )
          dr = enorm0_2d ( xq, yq, xr, yr )

          if ( di < dr ) then
            ir = i
            xr = x(ir)
            yr = y(ir)
            angmax = ai
          end if

        end if

      end if

    end do
!
!  If we've returned to our starting point, exit.
!
    if ( ir == istart ) then
      exit
    end if

    nval = nval + 1

    if ( nval > n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm failed.'
      stop
    end if
!
!  Set Q := P, P := R, and repeat.
!
    ival(nval) = ir
    xp = xq
    yp = yq
    iq = ir
    xq = xr
    yq = yr

  end do

  return
end
subroutine points_nearest_point_bins2_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_BINS2_2D finds the nearest point to a given point in 2D.
!
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINT_BINS_2D by calling
!    a subroutine to compute the next bin index.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    08 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, real PTEST(2), the coordinates of the test points.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
!
!    Output, real D_MIN, the distance between PTEST and PSET(*,I_MIN).
!
!    Output, integer COMPARES, the number of point-to-point comparisons.
!
  implicit none
!
  integer nbin
  integer, parameter :: ndim = 2
  integer nset
!
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin)
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares
  real d_min
  real d_min_sq
  real d_sq
  integer i
  integer i_min
  integer ic
  integer il
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim)
  real search_radius
!
  compares = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min = huge ( d_min )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    d_min = sqrt ( sum ( ( ptest(1:ndim) - pset(1:ndim,1) )**2 ) )
    compares = 1
    i_min = 1
    return
  end if

  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Search for each of the NTEST points.
!
  d_min_sq = huge ( d_min_sq )
  i_min = 0
  search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
  call r2_to_bin_even2 ( nbin, bin_min, bin_max, ptest, bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r2_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
  ic = bin(1)
  jc = bin(2)

  layer = 0
!
!  Search all legal bins in layer LAYER.
!
  do

    more_bins = .false.
    call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )
!
!  In layer LAYER, search each BIN I, J.
!
    do

      if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

        node = bin_start(i,j)

        do while ( node > 0 )

          d_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )
          compares = compares + 1

          if ( d_sq < d_min_sq ) then
            d_min_sq = d_sq
            i_min = node
          end if

          node = bin_next(node)

        end do

      end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
      do

        call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )

        if ( .not. more_bins ) then
          exit
        end if

        if ( 1 <= i .and. i <= nbin .and. &
             1 <= j .and. j <= nbin ) then
          exit
        end if

      end do

      if ( .not. more_bins ) then
        exit
      end if

    end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
    if ( layer == 0 ) then
      search_radius = min ( &
        minval ( abs ( ptest(1:ndim) - c_min(1:ndim) ) ), &
        minval ( abs ( ptest(1:ndim) - c_max(1:ndim) ) ) )
    else
      search_radius = search_radius + layer_width
    end if
!
!  We are done with PTEST if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
    if ( i_min /= 0 ) then
      d_min = sqrt ( d_min_sq )
      if ( search_radius >= d_min ) then
        exit
      end if
    end if

    layer = layer + 1

  end do
!
!  We are now done with all the layers.
!
  return
end
subroutine points_nearest_point_bins3_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_BINS3_2D finds the nearest point to a given point in 2D.
!
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by allowing the
!    user to specify the number of bins in each dimension.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN(1) by NBIN(2) regular grid of
!    cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = (/ 5, 4 /)
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt,
!    Mathematics Department,
!    Iowa State University.
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN(2), the number of cells in the horizontal and
!    vertical directions.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, real PTEST(2), the coordinates of the test point.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
!
!    Output, real D_MIN, the distance between PTEST and PSET(*,I_MIN).
!
!    Output, integer COMPARES, the number of point-to-point comparisons.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  integer nbin(ndim)
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2))
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin(1),nbin(2))
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares
  real d_min
  real d_min_sq
  real d_sq
  integer i
  integer i_min
  integer ic
  integer il
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim)
  real search_radius
!
  compares = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min = huge ( d_min )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    d_min = sqrt ( sum ( ( ptest(1:ndim) - pset(1:ndim,1) )**2 ) )
    i_min = 1
    compares = 1
    return
  end if
!
!  The efficiency of the code will suffer if the data in the vector
!
!    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) )
!
!  varies significantly.
!
  layer_width = minval ( &
    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) ) )

  d_min_sq = huge ( d_min_sq )
  i_min = 0
  search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
  call r2_to_bin_even3 ( nbin, bin_min, bin_max, ptest, bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r2_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
  ic = bin(1)
  jc = bin(2)

  layer = 0
!
!  Search all legal bins in layer LAYER.
!
  do

    more_bins = .false.
    call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )
!
!  In layer LAYER, search each BIN I, J.
!
    do

      if ( 1 <= i .and. i <= nbin(1) .and. 1 <= j .and. j <= nbin(2) ) then

        node = bin_start(i,j)

        do while ( node > 0 )

          d_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )
          compares = compares + 1

          if ( d_sq < d_min_sq ) then
            d_min_sq = d_sq
            i_min = node
          end if

          node = bin_next(node)

        end do

      end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
      do

        call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )

        if ( .not. more_bins ) then
          exit
        end if

        if ( 1 <= i .and. i <= nbin(1) .and. &
             1 <= j .and. j <= nbin(2) ) then
          exit
        end if

      end do

      if ( .not. more_bins ) then
        exit
      end if

    end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
    if ( layer == 0 ) then
      search_radius = min ( &
        minval ( abs ( ptest(1:ndim) - c_min(1:ndim) ) ), &
        minval ( abs ( ptest(1:ndim) - c_max(1:ndim) ) ) )
    else
      search_radius = search_radius + layer_width
    end if
!
!  We are done if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
    if ( i_min /= 0 ) then
      d_min = sqrt ( d_min_sq )
      if ( search_radius >= d_min ) then
        exit
      end if
    end if

    layer = layer + 1

  end do

  return
end
subroutine points_nearest_point_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, p, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_BINS_2D finds the nearest point to a given point in 2D.
!
!
!  Discussion:
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point P, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if P lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing P, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to P.  We now know that
!       we don't need to search any cell whose points will all be further
!       from P than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       P than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    15 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, real P(2), the point to be tested.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to P.
!
!    Output, real D_MIN, the distance between P and PSET(*,I_MIN).
!
!    Output, integer COMPARES, the number of point-to-point comparisons.
!
  implicit none
!
  integer nbin
  integer, parameter :: ndim = 2
  integer nset
!
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin)
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares
  real d_min
  real d_min_sq
  real d_sq
  integer i
  integer i_min
  integer ic
  integer il
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real p(ndim)
  real pset(ndim,nset)
  real search_radius
!
  compares = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min = huge ( d_min )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    d_min = sqrt ( sum ( ( p(1:ndim) - pset(1:ndim,1) )**2 ) )
    compares = 1
    i_min = 1
    return
  end if
!
!  Initialize.
!
  d_min_sq = huge ( d_min_sq )
  i_min = 0
  search_radius = 0.0E+00
  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Determine the bin coordinates of the point P.
!
  call r2_to_bin_even2 ( nbin, bin_min, bin_max, p, bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r2_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
  ic = bin(1)
  jc = bin(2)
!
!  Set
!  * the current layer,
!  * the starting bin of the current layer,
!  * the current bin
!
  layer = 0
  il = ic
  jl = jc
  i = il
  j = jl

  do
!
!  Search all legal bins in layer LAYER.
!
    do
!
!  Search BIN I, J.
!
      if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

        node = bin_start(i,j)

        do while ( node > 0 )

          d_sq = sum ( ( p(1:ndim) - pset(1:ndim,node) )**2 )
          compares = compares + 1

          if ( d_sq < d_min_sq ) then
            d_min_sq = d_sq
            i_min = node
          end if

          node = bin_next(node)

        end do

      end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
      more_bins = .true.

      do

        if ( i < ic + layer .and. j == jc - layer ) then
          i = i + 1
        else if ( i == ic + layer .and. j < jc + layer ) then
          j = j + 1
        else if ( ic - layer < i .and. j == jc + layer ) then
          i = i - 1
        else if ( i == ic - layer .and. jc - layer + 1 < j ) then
          j = j - 1
        else
          more_bins = .false.
          exit
        end if

        if ( 1 <= i .and. i <= nbin .and. &
             1 <= j .and. j <= nbin ) then
          exit
        end if

      end do

      if ( .not. more_bins ) then
        exit
      end if

    end do
!
!  We've completed this layer.
!  Update the radius of the searched area.
!
    if ( layer == 0 ) then
      search_radius = min ( &
        minval ( abs ( p(1:ndim) - c_min(1:ndim) ) ), &
        minval ( abs ( p(1:ndim) - c_max(1:ndim) ) ) )
    else
      search_radius = search_radius + layer_width
    end if
!
!  We can stop if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with P at the center and the nearest N on the circumference.
!
   if ( i_min /= 0 ) then
     d_min = sqrt ( d_min_sq )
     if ( search_radius >= d_min ) then
       exit
     end if
   end if
!
!  Prepare to search the next layer.
!
    layer = layer + 1

    il = ic - layer
    jl = jc - layer

    i = il
    j = jl

  end do

  return
end
subroutine points_nearest_point_del_2d ( point_num, xc, xd, nabes_first, &
  nabes_num, nabes_dim, nabes, nnear, dnear )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_DEL_2D searches a Delaunay triangulation for the nearest neighbor.
!
!
!  Discussion:
!
!    A set of points XC is given, along with its Delaunay triangulation.  
!    Now a new point XD is given, and we need to know the closest point in XC.
!
!    This algorithm starts at a random point in XC, and then repeatedly moves
!    to a neighboring point that is closer to XD.  This is guaranteed to be
!    possible because the triangulation is Delaunay.  Otherwise, it
!    would be possible to reach a vertex which was not the closest,
!    but for which all neighbors were further away.
!
!    This algorithm is able to handle the case where the point XD lies 
!    outside the convex hull.
!
!    The algorithm is very simple to code.  In the most likely
!    case, it should have an expected time complexity of O(sqrt(N)).
!
!    Overhead occurs in the development of the vertex adjacency data
!    structure.  The size of this array should be roughly 6*N on average.
!    Given the list of nodes that make up triangles, the vertex adjacency
!    data can be constructed by storing every pair of nodes (I,J) and (J,I),
!    and sorting the data into dictionary order.
!
!  Modified:
!
!    20 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, real XC(2,POINT_NUM), the coordinates of the set of points.
!
!    Input, real XD(2), a point whose nearest neighbor is to be found.
!
!    Input, integer NABES_FIRST(POINT_NUM), the index in NABES of the first
!    neighbor in the list for each node.
!
!    Input, integer NABES_NUM(POINT_NUM), the number of neighbors of each node.
!
!    Input, integer NABES_DIM, the dimension of NABES.
!
!    Input, integer NABES(NABES_DIM), a list of the neighbors of all the nodes.
!    Neighbors of node 1 are listed first, and so on.
!
!    Output, integer NNEAR, the nearest node to XD.
!
!    Output, real DNEAR, the distance of the nearest node to XD.
!
  implicit none
!
  integer nabes_dim
  integer point_num
!
  real dist
  real dnear
  integer i
  integer i1
  integer j
  integer nabes(nabes_dim)
  integer nabes_first(point_num)
  integer nabes_num(point_num)
  integer nnear
  integer nnear_old
  real x
  real x1
  real xc(2,point_num)
  real xd(2)
  real y
  real y1
!
  x = xd(1)
  y = xd(2)
!
!  Select a random vertex.
!
  nnear = 1
  x1 = xc(1,nnear)
  y1 = xc(2,nnear)
  dnear = sqrt ( ( x1 - x )**2 + ( y1 - y )**2 )
!
!  From the current vertex, consider all neighboring vertices.
!  For each neighbor, compute the distance to the point.
!  If no neighbor is closer, then the current vertex is the closest.
!  Otherwise, set the current vertex to the neighbor that was closest,
!  and repeat.
!
  do

    nnear_old = nnear

    j = nabes_first(nnear_old)

    do i = 1, nabes_num(nnear_old)

      i1 = nabes(j)
      x1 = xc(1,i1)
      y1 = xc(2,i1)
      dist = sqrt ( ( x1 - x )**2 + ( y1 - y )**2 )

      if ( dist < dnear ) then
        dnear = dist
        nnear = i1
      end if

      j = j + 1

    end do
!
!  If no neighbor was closer, we're done.
!
    if ( nnear == nnear_old ) then
      exit
    end if

  end do

  return
end
subroutine points_nearest_point_naive_2d ( nset, pset, ptest, i_min, d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_NAIVE_2D finds the nearest point to a given point in 2D.
!
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, real PTEST(2), the point whose nearest neighbor is sought.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to P.
!
!    Output, real D_MIN, the distance between P and PSET(*,I_MIN).
!
  implicit none
!
  integer nset
  integer, parameter :: ndim = 2
!
  real d
  real d_min
  integer i
  integer i_min
  real pset(ndim,nset)
  real ptest(ndim)
!
  d_min = huge ( d_min )
  i_min = 0

  do i = 1, nset
    d = sum ( ( ptest(1:ndim) - pset(1:ndim,i) )**2 )
    if ( d < d_min ) then
      d_min = d
      i_min = i
    end if
  end do

  d_min = sqrt ( d_min )

  return
end
subroutine points_nearest_point_naive_3d ( nset, pset, ptest, i_min, d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_NAIVE_3D finds the nearest point to a given point in 3D.
!
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(3,NSET), the coordinates of the points in the set.
!
!    Input, real PTEST(3), the point whose nearest neighbor is sought.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to P.
!
!    Output, real D_MIN, the distance between P and PSET(*,I_MIN).
!
  implicit none
!
  integer nset
  integer, parameter :: ndim = 3
!
  real d
  real d_min
  integer i
  integer i_min
  real pset(ndim,nset)
  real ptest(ndim)
!
  d_min = huge ( d_min )
  i_min = 0

  do i = 1, nset
    d = sum ( ( ptest(1:ndim) - pset(1:ndim,i) )**2 )
    if ( d < d_min ) then
      d_min = d
      i_min = i
    end if
  end do

  d_min = sqrt ( d_min )

  return
end
subroutine points_nearest_point_naive_nd ( ndim, n, pset, p, i_min, d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINT_NAIVE_ND finds the nearest point to a given point in ND.
!
!
!  Discussion:
!
!    A naive algorithm is used.  No attempt is made to optimize the
!    calculation, so there will be N distance calculations done.
!
!    For a large dataset, it would be better to group the points into
!    clusters, so that far fewer distance calculations need to be done.
!
!  Modified:
!
!    31 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of the points.
!
!    Input, integer N, the number of points in the set.
!
!    Input, real PSET(NDIM,N), the coordinates of the points in the set.
!
!    Input, real P(NDIM), the point to be tested.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to P.
!
!    Output, real D_MIN, the distance between P and PSET(*,I_MIN).
!
  implicit none
!
  integer n
  integer ndim
!
  real d
  real d_min
  integer i
  integer i_min
  real p(ndim)
  real pset(ndim,n)
!
  d_min = huge ( d_min )
  i_min = 0

  do i = 1, n

    d = sum ( ( p(1:ndim) - pset(1:ndim,i) )**2 )

    if ( d < d_min ) then
      d_min = d
      i_min = i
    end if

  end do
!
!  We save a little work by waiting til the end to take the square root.
!
  d_min = sqrt ( d_min )

  return
end
subroutine points_nearest_points_bins2_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS2_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by calling
!    a subroutine to compute the next bin index.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    08 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  implicit none
!
  integer nbin
  integer, parameter :: ndim = 2
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin)
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer il
  integer itest
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r2_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r2_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)

    layer = 0
!
!  Search all legal bins in layer LAYER.
!
    do

      more_bins = .false.
      call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )
!
!  In layer LAYER, search each BIN I, J.
!
      do

        if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

          node = bin_start(i,j)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

            node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
        do

          call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )

          if ( .not. more_bins ) then
            exit
          end if

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
      if ( i_min(itest) /= 0 ) then
        d_min(itest) = sqrt ( d_min_sq )
        if ( search_radius >= d_min(itest) ) then
          exit
        end if
      end if

      layer = layer + 1

    end do
!
!  We are now done with all the layers.
!
  end do
!
!  We are now done with all the test points.
!
  return
end
subroutine points_nearest_points_bins2_3d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS2_3D finds the nearest point to given points in 3D.
!
!
!  Discussion:
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J,K) and BIN_LAST(I,J,K) are given the coordinates
!    I, J, K of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J, K cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    15 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(3,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.  NBIN must be at least 3.
!
!    Input, real BIN_MIN(3), BIN_MAX(3), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(3,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  implicit none
!
  integer nbin
  integer, parameter :: ndim = 3
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin,nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin,nbin)
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer itest
  integer j
  integer jc
  integer k
  integer kc
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r3_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r3_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)
    kc = bin(3)

    layer = 0
!
!  Search all legal bins in layer LAYER.
!
    do

      more_bins = .false.

      call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, i, j, k, &
        more_bins )
!
!  In layer LAYER, search each BIN I, J, K.
!
      do

        if ( 1 <= i .and. i <= nbin .and. &
             1 <= j .and. j <= nbin .and. &
             1 <= k .and. k <= nbin ) then

          node = bin_start(i,j,k)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

             node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J, K.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
        do

          call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, &
            i, j, k, more_bins )

          if ( .not. more_bins ) then
            exit
          end if

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin .and. &
               1 <= k .and. k <= nbin ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
     if ( i_min(itest) /= 0 ) then
       d_min(itest) = sqrt ( d_min_sq )
       if ( search_radius >= d_min(itest) ) then
         exit
       end if
     end if
!
!  Prepare to search the next layer.
!
      layer = layer + 1

    end do

  end do

  return
end
subroutine points_nearest_points_bins3_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS3_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by allowing the
!    user to specify the number of bins in each dimension.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN(1) by NBIN(2) regular grid of
!    cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = (/ 5, 4 /)
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN(2), the number of cells in the horizontal and
!    vertical directions.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  integer nbin(ndim)
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2))
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin(1),nbin(2))
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer il
  integer itest
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min(1) )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if
!
!  The efficiency of the code will suffer if the data in the vector
!
!    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) )
!
!  varies significantly.
!
  layer_width = minval ( &
    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) ) )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r2_to_bin_even3 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r2_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)

    layer = 0
!
!  Search all legal bins in layer LAYER.
!
    do

      more_bins = .false.
      call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )
!
!  In layer LAYER, search each BIN I, J.
!
      do

        if ( 1 <= i .and. i <= nbin(1) .and. 1 <= j .and. j <= nbin(2) ) then

          node = bin_start(i,j)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

            node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
        do

          call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )

          if ( .not. more_bins ) then
            exit
          end if

          if ( 1 <= i .and. i <= nbin(1) .and. &
               1 <= j .and. j <= nbin(2) ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
      if ( i_min(itest) /= 0 ) then
        d_min(itest) = sqrt ( d_min_sq )
        if ( search_radius >= d_min(itest) ) then
          exit
        end if
      end if

      layer = layer + 1

    end do
!
!  We are now done with all the layers.
!
  end do
!
!  We are now done with all the test points.
!
  return
end
subroutine points_nearest_points_bins3_3d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS3_3D finds the nearest point to given points in 3D.
!
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_3D by allowing the
!    user to specify the number of bins in each dimension.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    box.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN(1) by NBIN(2) by NBIN(3)
!    regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = (/ 5, 4, 2 /)
!
!             Z LAYER 1                       Z LAYER 2
!     *----*----*----*----*----*     *----*----*----*----*----*
!     |    |    |    |    |    |     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |     | 36 | 37 | 38 | 39 | 40 |
!     |    |    |    |    |    |     |    |    |    |    |    |
!     *----*----*----*----*----*     *----*----*----*----*----*
!     |    |    |    |    |    |     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |     | 31 | 32 | 33 | 34 | 35 |
!     |    |    |    |    |    |     |    |    |    |    |    |
!     *----*----*----*----*----*     *----*----*----*----*----*
!     |    |    |    |    |    |     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |     | 26 | 27 | 28 | 29 | 30 |
!     |    |    |    |    |    |     |    |    |    |    |    |
!     *----*----*----*----*----*     *----*----*----*----*----*
!     |    |    |    |    |    |     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |     |    |    |    |    |    |
!     *----*----*----*----*----*     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!      * P1 is in the same cell as P2, P1.X = P2.X, P1.Y = P2.Y, but P1.Z < P2.Z
!
!    The arrays BIN_START(I,J,K) and BIN_LAST(I,J,K) are given the coordinates
!    I, J, K of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J, K cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    30 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(3,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN(3), the number of cells in the X, Y and Z directions.
!
!    Input, real BIN_MIN(3), BIN_MAX(3), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(3,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  implicit none
!
  integer, parameter :: ndim = 3
!
  integer nbin(ndim)
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2),nbin(3))
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin(1),nbin(2),nbin(3))
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer itest
  integer j
  integer jc
  integer k
  integer kc
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min(1) )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if
!
!  The efficiency of the code will suffer if the data in the vector
!
!    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) )
!
!  varies significantly.
!
  layer_width = minval ( &
    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) ) )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r3_to_bin_even3 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r3_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)
    kc = bin(3)

    layer = 0
!
!  Search all legal bins in layer LAYER.
!
    do

      more_bins = .false.
      call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, i, j, k, &
        more_bins )
!
!  In layer LAYER, search each BIN I, J, K.
!
      do

        if ( 1 <= i .and. i <= nbin(1) .and. &
             1 <= j .and. j <= nbin(2) .and. &
             1 <= k .and. k <= nbin(3) ) then

          node = bin_start(i,j,k)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

            node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J, K.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
        do

          call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, i, j, k, &
            more_bins )

          if ( .not. more_bins ) then
            exit
          end if

          if ( 1 <= i .and. i <= nbin(1) .and. &
               1 <= j .and. j <= nbin(2) .and. &
               1 <= k .and. k <= nbin(3) ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
      if ( i_min(itest) /= 0 ) then
        d_min(itest) = sqrt ( d_min_sq )
        if ( search_radius >= d_min(itest) ) then
          exit
        end if
      end if

      layer = layer + 1

    end do
!
!  We are now done with all the layers.
!
  end do
!
!  We are now done with all the test points.
!
  return
end
subroutine points_nearest_points_bins4_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS4_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by allowing the
!    user to specify the number of bins in each dimension.  The main reason
!    for doing this is to efficiently handle problems where the extent
!    of the region varies widely from one dimension to another.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN(1) by NBIN(2) regular grid of
!    cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = (/ 5, 4 /)
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN(2), the number of cells in the horizontal and
!    vertical directions.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  integer nbin(ndim)
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2))
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin(1),nbin(2))
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  real cell_width_i
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  real disc
  real first_dist
  logical found_a_neighbor
  integer i
  integer i_min(ntest)
  integer ic
  integer il
  integer itest
  integer j
  integer jc
  integer jl
  integer k
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  integer rank
  integer search1(ndim)
  integer search2(ndim)
  logical searched_everywhere
  integer searched_hi(ndim)
  integer searched_hi_new(ndim)
  integer searched_lo(ndim)
  integer searched_lo_new(ndim)
  integer searched_new(ndim)
  real w
  real wall
  real wall_dist
  real width
  real width_i
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min(1) )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest
!
!  PART ONE: Initialize data.
!
!  Determine the bin coordinates of the test point.
!  Determine the limits of the bin containing P.
!  Search the bin.
!  Set the indices of the searched region to the index of this bin.
!
    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    found_a_neighbor = .false.
    searched_everywhere = .false.

    call r2_to_bin_even3 ( nbin, bin_min, bin_max, ptest(1,itest), bin )

    call bin_to_r2_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )

    call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, bin_next, &
      ptest(1,itest), found_a_neighbor, i_min(itest), d_min_sq, compares(itest) )

    if ( found_a_neighbor) then
      d_min(itest) = sqrt ( d_min_sq )
    end if

    searched_lo(1:ndim) = bin(1:ndim)
    searched_hi(1:ndim) = bin(1:ndim)

    if ( all ( nbin(1:ndim) == 1 ) ) then
      if ( found_a_neighbor ) then
        d_min(itest) = sqrt ( d_min_sq )
      else
        d_min(itest) = huge ( d_min(1) )
      end if
      searched_everywhere = .true.
      cycle
    end if
!
!  PART TWO. Look for a neighbor.
!
!  Expand the search area in each dimension.
!  Consider push the search limits out in every direction by 1 index.
!  Determine the maximum width of the search region achieved in this way.
!  Where possible, move all smaller indices out to the maximum width.
!
!  Organize the search of the annexed region by dimension.
!  Decrease and increase the first dimension to its new limits.
!  Repeat for each dimension.
!
!  Jump to next section ALMOST as soon as you have found a neighbor;
!  just finish up the search in that coordinate direction, so that the
!  search can be picked up cleanly in the final section.
!
    do while ( .not. found_a_neighbor )

      do i = 1, ndim
        searched_hi_new(i) = min ( searched_hi(i) + 1, nbin(i) )
        searched_lo_new(i) = max ( searched_lo(i) - 1, 1 )
      end do

      width = 0.0E+00
      do i = 1, ndim
        width_i = ( searched_hi_new(i) + 1 - searched_lo_new(i) ) * &
            ( bin_max(i) - bin_min(i) ) / nbin(i)
        width = max ( width, width_i )
      end do

      do i = 1, ndim

        cell_width_i = ( bin_max(i) - bin_min(i) ) / nbin(i)

        width_i = ( searched_hi_new(i) + 1 - searched_lo_new(i) ) * cell_width_i

        disc = width - width_i

        if ( disc >= 2.0E+00 * cell_width_i ) then

          k = int ( disc / ( 2.0E+00 * cell_width_i ) )

          searched_hi_new(i) = searched_hi_new(i) + k
          searched_lo_new(i) = searched_lo_new(i) - k

          searched_hi_new(i) = min ( searched_hi_new(i), nbin(i) )
          searched_lo_new(i) = max ( searched_lo_new(i), 1 )

        end if

      end do

      searched_everywhere = .true.
      do i = 1, ndim
        if ( searched_hi_new(i) > searched_hi(i) ) then
          searched_everywhere = .false.
          exit
        end if
        if ( searched_lo_new(i) < searched_lo(i) ) then
          searched_everywhere = .false.
          exit
        end if
      end do

      if ( searched_everywhere ) then
        exit
      end if

      do i = 1, ndim

        if ( searched_lo_new(i) < searched_lo(i) ) then

          search1(1:ndim) = searched_lo(1:ndim)
          search1(i) = searched_lo_new(i)

          search2(1:i-1) = searched_hi(1:i-1)
          search2(i) = searched_lo(i) - 1
          search2(i+1:ndim) = searched_lo(i+1:ndim)

          rank = 0

          do

            call tuple_next2 ( ndim, search1, search2, bin, rank )

            if ( rank == 0 ) then
              exit
            end if

            call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, &
              bin_next, ptest(1,itest), found_a_neighbor, i_min(itest), &
              d_min_sq, compares(itest) )

          end do

          searched_lo(i) = searched_lo_new(i)

          if ( found_a_neighbor ) then
            exit
          end if

        end if

        if ( searched_hi_new(i) > searched_hi(i) ) then

          search1(1:i-1) = searched_lo(1:ndim)
          search1(i) = searched_hi(i) + 1
          search1(i+1:ndim) = searched_hi(1:ndim)

          search2(1:ndim) = searched_hi(1:ndim)
          search2(i) = searched_hi_new(i)

          rank = 0

          do

            call tuple_next2 ( ndim, search1, search2, bin, rank )

            if ( rank == 0 ) then
              exit
            end if

            call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, &
              bin_next, ptest(1,itest), found_a_neighbor, i_min(itest), &
              d_min_sq, compares(itest) )

          end do

          searched_hi(i) = searched_hi_new(i)

          if ( found_a_neighbor ) then
            exit
          end if

        end if

      end do

      d_min(itest) = sqrt ( d_min_sq )

    end do
!
!  PART THREE: Final search
!
!  You have found a neighbor in PSET to PTEST.
!  If the neighbor is closer than the nearest wall which might have
!    something on the other side, you're done.
!  Otherwise, expand the region enough in each direction so that, once
!    it is searched, we are sure to be done.
!
    wall_dist = huge ( wall_dist )

    do i = 1, ndim
      if ( searched_lo(i) > 1 ) then
        wall = ( real ( nbin(i) + 1 - searched_lo(i) ) * bin_min(i) &
                 + real ( searched_lo(i) - 1 ) * bin_max(i) ) / nbin(i)
        wall_dist = min ( wall_dist, ptest(i,itest) - wall )
      end if
      if ( searched_hi(i) < nbin(i) ) then
        wall = ( real ( nbin(i) - searched_hi(i) ) * bin_min(i) &
                 + real ( searched_hi(i) ) * bin_max(i) ) / nbin(i)
        wall_dist = min ( wall_dist, wall - ptest(i,itest) )
      end if

    end do

    first_dist = d_min(itest)

    if ( first_dist < wall_dist ) then
      cycle
    end if

    do i = 1, ndim

      if ( searched_lo(i) > 1 ) then
!
!  Solve for SEARCH_NEW(I) so that PTEST(I,ITEST) - WALL > FIRST_DIST.
!
        w = ( ptest(i,itest) - first_dist - bin_min(i) ) &
          / real ( bin_max(i) - bin_min(i) )
        k = int ( real ( nbin(i) ) * w )

        k = max ( k, 1 )

        search1(1:ndim) = searched_lo(1:ndim)
        search1(i) = k
        search2(1:ndim) = searched_hi(1:ndim)
        search2(i) = searched_lo(i) - 1

        if ( search1(i) <= search2(i) ) then

          rank = 0

          do

            call tuple_next2 ( ndim, search1, search2, bin, rank )

            if ( rank == 0 ) then
              exit
            end if

            call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, bin_next, &
              ptest(1,itest), found_a_neighbor, i_min(itest), d_min_sq, &
              compares(itest) )

          end do

          searched_lo(i) = k

        end if

      end if

      if ( searched_hi(i) < nbin(i) ) then
!
!  Solve for SEARCH_NEW(I) so that WALL - PTEST(I,ITEST) > FIRST_DIST.
!
        w = ( ptest(i,itest) + first_dist - bin_min(i) ) &
          / real ( bin_max(i) - bin_min(i) )
        k = 1 + int ( real ( nbin(i) ) * w )

        k = min ( k, nbin(i) )

        search1(1:ndim) = searched_lo(1:ndim)
        search1(i) = searched_hi(i)+1
        search2(1:ndim) = searched_hi(1:ndim)
        search2(i) = k

        if ( search1(i) <= search2(i) ) then

          rank = 0

          do

            call tuple_next2 ( ndim, search1, search2, bin, rank )

            if ( rank == 0 ) then
              exit
            end if

            call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, bin_next, &
              ptest(1,itest), found_a_neighbor, i_min(itest), d_min_sq, &
              compares(itest) )

          end do

          searched_hi(i) = k

        end if

      end if

    end do

    d_min(itest) = sqrt ( d_min_sq )

  end do

  return
end
subroutine points_nearest_points_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, d_min, compares )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_BINS_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Modified:
!
!    08 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NBIN, the number of cells.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the minimum and maximum bin values.
!
!    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), indicates
!    the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST) and
!    PSET(*,I_MIN).
!
!    Output, integer COMPARES(NTEST), the number of point-to-point comparisons.
!
  implicit none
!
  integer nbin
  integer, parameter :: ndim = 2
  integer nset
  integer ntest
!
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer bin_start(nbin,nbin)
  integer bin_next(nset)
  real c_max(ndim)
  real c_min(ndim)
  integer compares(ntest)
  real d_min(ntest)
  real d_min_sq
  real d_sq
  integer i
  integer i_min(ntest)
  integer ic
  integer il
  integer itest
  integer j
  integer jc
  integer jl
  integer layer
  real layer_width
  logical more_bins
  integer node
  real pset(ndim,nset)
  real ptest(ndim,ntest)
  real search_radius
!
  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min(1:ntest) = huge ( d_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      d_min(itest) = sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    d_min_sq = huge ( d_min_sq )
    i_min(itest) = 0
    search_radius = 0.0E+00
!
!  Determine the bin coordinates of the point P.
!
    call r2_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r2_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)
!
!  Set
!  * the current layer,
!  * the starting bin of the current layer,
!  * the current bin
!
    layer = 0
    il = ic
    jl = jc
    i = il
    j = jl

    do
!
!  Search all legal bins in layer LAYER.
!
      do
!
!  Search BIN I, J.
!
        if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

          node = bin_start(i,j)

          do while ( node > 0 )

            d_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( d_sq < d_min_sq ) then
              d_min_sq = d_sq
              i_min(itest) = node
            end if

            node = bin_next(node)

          end do

        end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
        more_bins = .true.

        do

          if ( i < ic + layer .and. j == jc - layer ) then
            i = i + 1
          else if ( i == ic + layer .and. j < jc + layer ) then
            j = j + 1
          else if ( ic - layer < i .and. j == jc + layer ) then
            i = i - 1
          else if ( i == ic - layer .and. jc - layer + 1 < j ) then
            j = j - 1
          else
            more_bins = .false.
            exit
          end if

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed this layer.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
     if ( i_min(itest) /= 0 ) then
       d_min(itest) = sqrt ( d_min_sq )
       if ( search_radius >= d_min(itest) ) then
         exit
       end if
     end if
!
!  Prepare to search the next layer.
!
      layer = layer + 1

      il = ic - layer
      jl = jc - layer

      i = il
      j = jl

    end do

  end do

  return
end
subroutine points_nearest_points_naive_2d ( nset, pset, ntest, ptest, i_min, &
  d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_NAIVE_2D finds the nearest point to given points in 2D.
!
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(2,NSET), the coordinates of the points in the set.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(2,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST)
!    and PSET(*,I_MIN).
!
  implicit none
!
  integer nset
  integer, parameter :: ndim = 2
  integer ntest
!
  real d
  real d_min(ntest)
  integer i
  integer i_min(ntest)
  integer itest
  real pset(ndim,nset)
  real ptest(ndim,ntest)
!
  do itest = 1, ntest

    d_min(itest) = huge ( d_min )
    i_min(itest) = 0

    do i = 1, nset
      d = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,i) )**2 )
      if ( d < d_min(itest) ) then
        d_min(itest) = d
        i_min(itest) = i
      end if
    end do

    d_min(itest) = sqrt ( d_min(itest) )

  end do

  return
end
subroutine points_nearest_points_naive_3d ( nset, pset, ntest, ptest, i_min, &
  d_min )
!
!*******************************************************************************
!
!! POINTS_NEAREST_POINTS_NAIVE_3D finds the nearest point to given points in 3D.
!
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real PSET(3,NSET), the coordinates of the points in the set.
!
!    Input, integer NTEST, the number of test points.
!
!    Input, real PTEST(3,NTEST), the coordinates of the test points.
!
!    Output, integer I_MIN(NTEST), the index of the nearest point in PSET
!    to PTEST(ITEST).
!
!    Output, real D_MIN(NTEST), the distance between PTEST(*,ITEST)
!    and PSET(*,I_MIN).
!
  implicit none
!
  integer nset
  integer, parameter :: ndim = 3
  integer ntest
!
  real d
  real d_min(ntest)
  integer i
  integer i_min(ntest)
  integer itest
  real pset(ndim,nset)
  real ptest(ndim,ntest)
!
  do itest = 1, ntest

    d_min(itest) = huge ( d_min )
    i_min(itest) = 0

    do i = 1, nset
      d = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,i) )**2 )
      if ( d < d_min(itest) ) then
        d_min(itest) = d
        i_min(itest) = i
      end if
    end do

    d_min(itest) = sqrt ( d_min(itest) )

  end do

  return
end
subroutine polygon_area_2_2d ( n, x, y, area )
!
!*******************************************************************************
!
!! POLYGON_AREA_2_2D computes the area of a polygon in 2D.
!
!
!  Formula:
!
!    The area is the sum of the areas of the triangles formed by
!    node N with consecutive pairs of nodes.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real X(N), Y(N), the coordinates of the vertices.
!
!    Output, real AREA, the absolute area of the polygon.
!
  implicit none
!
  integer n
!
  real area
  real areat
  integer i
  real x(n)
  real y(n)
!
  area = 0.0E+00

  do i = 1, n - 2

    call triangle_area_2d ( x(i), y(i), x(i+1), y(i+1), x(n), y(n), areat )

    area = area + areat

  end do

  return
end
subroutine polygon_area_2_3d ( n, x, y, z, area )
!
!*******************************************************************************
!
!! POLYGON_AREA_2_3D computes the area of a polygon in 3D.
!
!
!  Formula:
!
!    The area is the sum of the areas of the triangles formed by
!    node N with consecutive pairs of nodes.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real X(N), Y(N), Z(N), the coordinates of the vertices.
!
!    Output, real AREA, the absolute area of the polygon.
!
  implicit none
!
  integer n
!
  real area
  real areat
  integer i
  real x(n)
  real y(n)
  real z(n)
!
  area = 0.0E+00

  do i = 1, n - 2

    call triangle_area_3d ( x(i), y(i), z(i), x(i+1), y(i+1), z(i+1), &
      x(n), y(n), z(n), areat )

    area = area + areat

  end do

  return
end
subroutine polygon_area_2d ( n, x, y, area )
!
!*******************************************************************************
!
!! POLYGON_AREA_2D computes the area of a polygon in 2D.
!
!
!  Formula:
!
!    AREA = ABS ( 0.5 * SUM ( I = 1 to N ) X(I) * ( Y(I+1) - Y(I-1) ) )
!    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real X(N), Y(N), the coordinates of the vertices.
!
!    Output, real AREA, the absolute area of the polygon.
!
  implicit none
!
  integer n
!
  real area
  integer i
  integer im1
  integer ip1
  real x(n)
  real y(n)
!
  area = 0.0E+00

  do i = 1, n

    if ( i > 1 ) then
      im1 = i - 1
    else
      im1 = n
    end if

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    area = area + x(i) * ( y(ip1) - y(im1) )

  end do

  area = 0.5E+00 * abs ( area )

  return
end
subroutine polygon_area_3d ( n, x, y, z, area, normal )
!
!*******************************************************************************
!
!! POLYGON_AREA_3D computes the area of a polygon in 3D.
!
!
!  Restriction:
!
!    The computation is not valid unless the vertices really do lie
!    in a plane, so that the polygon that is defined is "flat".
!    The polygon does not have to be "regular", that is, neither its
!    sides nor its angles need to be equal.
!
!  Reference:
!
!    Allen Van Gelder,
!    Efficient Computation of Polygon Area and Polyhedron Volume,
!    Graphics Gems V, edited by Alan Paeth,
!    AP Professional, 1995.
!
!  Modified:
!
!    30 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of vertices.
!
!    Input, real X(N), Y(N), Z(N), the coordinates of the vertices.
!    The vertices should be listed in neighboring order.
!
!    Output, real AREA, the area of the polygon.
!
!    Output, real NORMAL(3), the unit normal vector to the polygon.
!
  implicit none
!
  integer n
!
  real area
  integer i
  integer ip1
  real normal(3)
  real x(n)
  real x1
  real x2
  real x3
  real y(n)
  real y1
  real y2
  real y3
  real z(n)
  real z1
  real z2
  real z3
!
  normal(1) = 0.0E+00
  normal(2) = 0.0E+00
  normal(3) = 0.0E+00

  do i = 1, n

    x1 = x(i)
    y1 = y(i)
    z1 = z(i)

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    x2 = x(ip1)
    y2 = y(ip1)
    z2 = z(ip1)

    call cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )

    normal(1) = normal(1) + x3
    normal(2) = normal(2) + y3
    normal(3) = normal(3) + z3

  end do

  area = sqrt ( normal(1)**2 + normal(2)**2 + normal(3)**2 )

  if ( area /= 0.0E+00 ) then
    normal(1) = normal(1) / area
    normal(2) = normal(2) / area
    normal(3) = normal(3) / area
  else
    normal(1) = 1.0E+00
    normal(2) = 0.0E+00
    normal(3) = 0.0E+00
  end if

  area = 0.5E+00 * area

  return
end
subroutine polygon_centroid_2_2d ( n, x, y, cx, cy )
!
!*******************************************************************************
!
!! POLYGON_CENTROID_2_2D computes the centroid of a polygon in 2D.
!
!
!  Method:
!
!    The centroid is the area-weighted sum of the centroids of
!    disjoint triangles that make up the polygon.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real X(N), Y(N), the coordinates of the vertices.
!
!    Output, real CX, CY, the coordinates of the centroid.
!
  implicit none
!
  integer n
!
  real area
  real areat
  real cx
  real cy
  integer i
  real x(n)
  real y(n)
!
  area = 0.0E+00
  cx = 0.0E+00
  cy = 0.0E+00

  do i = 1, n - 2

    call triangle_area_2d ( x(i), y(i), x(i+1), y(i+1), x(n), y(n), areat )

    area = area + areat
    cx = cx + areat * ( x(i) + x(i+1) + x(n) ) / 3.0E+00
    cy = cy + areat * ( y(i) + y(i+1) + y(n) ) / 3.0E+00

  end do

  cx = cx / area
  cy = cy / area

  return
end
subroutine polygon_centroid_2d ( n, x, y, cx, cy )
!
!*******************************************************************************
!
!! POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
!
!
!  Formula:
!
!    Denoting the centroid coordinates by (CX,CY), then
!
!      CX = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
!      CY = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
!
!    Green's theorem states that
!
!      Integral ( Polygon boundary ) ( M dx + N dy ) =
!      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
!
!    Using M = 0 and N = x**2/2, we get:
!
!      CX = 0.5 * Integral ( Polygon boundary ) x**2 dy,
!
!    which becomes
!
!      CX = 1/6 SUM ( I = 1 to N )
!        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
!
!    where, when I = N, the index "I+1" is replaced by 1.
!
!    A similar calculation gives us a formula for CY.
!
!  Reference:
!
!    Gerard Bashein and Paul Detmer,
!    Centroid of a Polygon,
!    Graphics Gems IV, edited by Paul Heckbert,
!    AP Professional, 1994.
!
!  Modified:
!
!    22 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of sides of the polygonal shape.
!
!    Input, real X(N), Y(N), the coordinates of the vertices of the shape.
!
!    Output, real CX, CY, the coordinates of the centroid of the shape.
!
  implicit none
!
  integer n
!
  real area
  real cx
  real cy
  integer i
  integer ip1
  real temp
  real x(n)
  real y(n)
!
  area = 0.0E+00
  cx = 0.0E+00
  cy = 0.0E+00

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    temp = ( x(i) * y(ip1) - x(ip1) * y(i) )

    area = area + temp
    cx = cx + ( x(ip1) + x(i) ) * temp
    cy = cy + ( y(ip1) + y(i) ) * temp

  end do

  area = area / 2.0E+00

  cx = cx / ( 6.0E+00 * area )
  cy = cy / ( 6.0E+00 * area )

  return
end
subroutine polygon_centroid_3d ( n, x, y, z, cx, cy, cz )
!
!*******************************************************************************
!
!! POLYGON_CENTROID_3D computes the centroid of a polygon in 3D.
!
!
!  Method:
!
!    The centroid is the area-weighted sum of the centroids of
!    disjoint triangles that make up the polygon.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of vertices of the polygon.
!
!    Input, real X(N), Y(N), Z(N), the coordinates of the vertices.
!
!    Output, real CX, CY, CZ, the coordinates of the centroid.
!
  implicit none
!
  integer n
!
  real area
  real areat
  real cx
  real cy
  real cz
  integer i
  real x(n)
  real y(n)
  real z(n)
!
  area = 0.0E+00
  cx = 0.0E+00
  cy = 0.0E+00
  cz = 0.0E+00

  do i = 1, n - 2

    call triangle_area_3d ( x(i), y(i), z(i), x(i+1), &
      y(i+1), z(i+1), x(n), y(n), z(n), areat )

    area = area + areat
    cx = cx + areat * ( x(i) + x(i+1) + x(n) ) / 3.0E+00
    cy = cy + areat * ( y(i) + y(i+1) + y(n) ) / 3.0E+00
    cz = cz + areat * ( z(i) + z(i+1) + z(n) ) / 3.0E+00

  end do

  cx = cx / area
  cy = cy / area
  cz = cz / area

  return
end
subroutine polygon_contains_point_2_2d ( n, xn, xval, yn, yval, inside )
!
!*******************************************************************************
!
!! POLYGON_CONTAINS_POINT_2_2D finds if a point is inside a convex polygon in 2D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes or vertices in the polygon.
!    N must be at least 3.
!
!    Input, real XN(N), the X coordinates of the vertices.
!
!    Input, real XVAL, the X coordinate of the point to be tested.
!
!    Input, real YN(N), the Y coordinates of the vertices.
!
!    Input, real YVAL, the Y coordinate of the point to be tested.
!
!    Output, logical INSIDE, is .TRUE. if ( X,Y) is inside
!    the polygon or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  integer n
!
  integer i
  logical inside
  real x1
  real x2
  real x3
  real xn(n)
  real xval
  real y1
  real y2
  real y3
  real yn(n)
  real yval
!
  inside = .false.
!
!  A point is inside a convex polygon if and only if it is inside
!  one of the triangles formed by X(1),Y(1) and any two consecutive
!  points on the polygon's circumference.
!
  x1 = xn(1)
  y1 = yn(1)

  do i = 2, n-1

    x2 = xn(i)
    y2 = yn(i)
    x3 = xn(i+1)
    y3 = yn(i+1)

    call triangle_contains_point_1_2d ( x1, y1, x2, y2, x3, y3, xval, yval, &
      inside )

    if ( inside ) then
      return
    end if

  end do

  return
end
subroutine polygon_contains_point_2d ( n, xn, xval, yn, yval, inside )
!
!*******************************************************************************
!
!! POLYGON_CONTAINS_POINT_2D finds if a point is inside a simple polygon in 2D.
!
!
!  Discussion:
!
!    A simple polygon is one whose boundary never crosses itself.
!
!  Reference:
!
!    ACM Algorithm 112.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes or vertices in the polygon.
!    N must be at least 3.
!
!    Input, real XN(N), the X coordinates of the vertices.
!
!    Input, real XVAL, the X coordinate of the point to be tested.
!
!    Input, real YN(N), the Y coordinates of the vertices.
!
!    Input, real YVAL, the Y coordinate of the point to be tested.
!
!    Output, logical INSIDE, is .TRUE. if ( X,Y) is inside
!    the polygon and .FALSE. otherwise.
!
  implicit none
!
  integer n
!
  integer i
  logical inside
  real x1
  real x2
  real xn(n)
  real xval
  real y1
  real y2
  real yn(n)
  real yval
!
  inside = .false.

  do i = 1, n

    x1 = xn(i)
    y1 = yn(i)

    if ( i < n ) then
      x2 = xn(i+1)
      y2 = yn(i+1)
    else
      x2 = xn(1)
      y2 = yn(1)
    end if

    if ( ( y1 < yval .eqv. yval <= y2 ) .and. &
      ( xval - x1 ) * abs ( y2 - y1 ) < ( x2 - x1 ) * ( yval - y1 ) ) then

      inside = .not. inside

    end if

  end do

  return
end
function polygon_convex ( n, x, y )
!
!*******************************************************************************
!
!! POLYGON_CONVEX determines whether a polygon is convex in 2D.
!
!
!  Discussion:
!
!    If the polygon has less than 3 distinct vertices, it is
!    classified as convex degenerate.
!
!    If the polygon "goes around" more than once, it is classified
!    as NOT convex.
!
!  Reference:
!
!    Peter Schorn and Frederick Fisher,
!    Testing the Convexity of a Polygon,
!    Graphics Gems, 1994.
!
!  Modified:
!
!    02 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer N, the number of vertices.
!
!    Input/output, X(N), Y(N), the coordinates of the vertices of the
!    polygon.  On output, duplicate consecutive points have been deleted,
!    and the vertices have been reordered so that the lexicographically
!    least point comes first.
!
!    Output, integer POLYGON_CONVEX:
!    -1, the polygon is not convex;
!     0, the polygon has less than 3 vertices; it is "degenerately" convex;
!     1, the polygon is convex and counterclockwise;
!     2, the polygon is convex and clockwise.
!
  implicit none
!
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real, parameter :: RAD_TO_DEG = 180.0E+00 / pi
!
  integer n
!
  real angle
  integer, parameter :: CONVEX_CCW = 1
  integer, parameter :: CONVEX_CW = 2
  real cross
  integer, parameter :: DEGENERATE_CONVEX = 0
  real dot
  real exterior_total
  integer i
  integer ip1
  integer ip2
  integer, parameter :: NOT_CONVEX = -1
  integer polygon_convex
  real sense
  real, parameter :: tol = 1.0E+00
  real x(n)
  real y(n)
!
  exterior_total = 0.0E+00
!
!  If there are not at least 3 distinct vertices, we are done.
!
  if ( n < 3 ) then
    polygon_convex = DEGENERATE_CONVEX
    return
  end if

  sense = 0.0E+00
!
!  Consider each polygonal vertex I.
!
  do i = 1, n

    ip1 = i + 1
    if ( ip1 > n ) then
      ip1 = ip1 - n
    end if

    ip2 = i + 2
    if ( ip2 > n ) then
      ip2 = ip2 - n
    end if

    dot =   ( x(ip2) - x(ip1) ) * ( x(i) - x(ip1) ) &
          + ( y(ip2) - y(ip1) ) * ( y(i) - y(ip1) )

    cross =   ( x(ip2) - x(ip1) ) * ( y(i) - y(ip1) ) &
            - ( x(i) - x(ip1) ) * ( y(ip2) - y(ip1) )

    angle = atan2 ( cross, dot )
!
!  See if the turn defined by this vertex is our first indication of
!  the "sense" of the polygon, or if it disagrees with the previously
!  defined sense.
!
    if ( sense == 0.0E+00 ) then

      if ( angle < 0.0E+00 ) then
        sense = -1.0E+00
      else if ( angle > 0.0E+00 ) then
        sense = +1.0E+00
      end if

   else if ( sense == 1.0E+00 ) then

      if ( angle < 0.0E+00 ) then
        polygon_convex = NOT_CONVEX
        return
      end if

    else if ( sense == -1.0E+00 ) then

      if ( angle > 0.0E+00 ) then
        polygon_convex = NOT_CONVEX
        return
      end if

    end if
!
!  If the exterior total is greater than 360, then the polygon is
!  going around again.
!
    angle = atan2 ( -cross, -dot )

    exterior_total = exterior_total + angle

    if ( abs ( exterior_total ) * RAD_TO_DEG > 360.0E+00 + TOL ) then
      polygon_convex = NOT_CONVEX
      return
    end if

  end do

  if ( sense == +1.0E+00 ) then
    polygon_convex = CONVEX_CCW
  else if ( sense == -1.0E+00 ) then
    polygon_convex = CONVEX_CW
  end if

  return
end
subroutine polygon_inrad_data_2d ( area, n, radin, radout, side )
!
!*******************************************************************************
!
!! POLYGON_INRAD_DATA_2D determines polygonal data from its inner radius in 2D.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real AREA, the area of the regular polygon.
!
!    Input, integer N, the number of sides of the polygon.
!    N must be at least 3.
!
!    Input, real RADIN, the inner radius of the polygon, that is,
!    the radius of the largest circle that can be inscribed within
!    the polygon.
!
!    Output, real RADOUT, the outer radius of the polygon, that is,
!    the radius of the smallest circle that can be described about
!    the polygon.
!
!    Output, real SIDE, the length of one side of the polygon.
!
  implicit none
!
  real angle
  real area
  integer n
  real r_pi
  real radin
  real radout
  real side
!
  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_INRAD_DATA_2D - Fatal error!'
    write ( *, '(a)' ) '  Input value of N must be at least 3'
    write ( *, '(a,i6)' ) '  but your input value was N = ', n
    stop
  end if

  angle = r_pi ( ) / real ( n )
  area = n * radin * radin * tan ( angle )
  side = 2.0E+00 * radin * tan ( angle )
  radout = 0.5E+00 * side / sin ( angle )

  return
end
subroutine polygon_lattice_area_2d ( i, b, area )
!
!*******************************************************************************
!
!! POLYGON_LATTICE_AREA_2D computes the area of a lattice polygon in 2D.
!
!
!  Discussion:
!
!    We define a lattice to be the 2D plane, in which the points
!    whose (X,Y) coordinates are both integers are given a special
!    status as "lattice points".
!
!    A lattice polygon is a polygon whose vertices are lattice points.
!
!    The area of a lattice polygon can be computed by Pick's Theorem:
!
!      Area = I + B / 2 - 1
!
!    where
!
!      I = the number of lattice points contained strictly inside the polygon;
!
!      B = the number of lattice points that lie exactly on the boundary.
!    
!  Reference:
!
!    Branko Gruenbaum and G C Shephard,
!    Pick's Theorem,
!    The American Mathematical Monthly,
!    Volume 100, 1993, pages 150-161.
!
!  Modified:
!
!    05 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number of interior lattice points.
!
!    Input, integer B, the number of boundary lattice points.
!
!    Output, real AREA, the area of the lattice polygon.
!
  implicit none
!
  integer n
!
  real area
  integer b
  integer i
!
  area = real ( i ) + real ( b ) / 2.0E+00 - 1.0E+00

  return
end
subroutine polygon_outrad_data_2d ( area, n, radin, radout, side )
!
!*******************************************************************************
!
!! POLYGON_OUTRAD_DATA_2D determines polygonal data from its outer radius in 2D.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real AREA, the area of the regular polygon.
!
!    Input, integer N, the number of sides of the polygon.
!    N must be at least 3.
!
!    Output, real RADIN, the inner radius of the polygon, that is,
!    the radius of the largest circle that can be inscribed
!    within the polygon.
!
!    Input, real RADOUT, the outer radius of the polygon, that is,
!    the radius of the smallest circle that can be described
!    around the polygon.
!
!    Output, real SIDE, the length of one side of the polygon.
!
  implicit none
!
  real angle
  real area
  integer n
  real r_pi
  real radin
  real radout
  real side
!
  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_OUTRAD_DATA_2D - Fatal error!'
    write ( *, '(a)' ) '  Input value of N must be at least 3'
    write ( *, '(a,i6)' ) '  but your input value was N = ', n
    stop
  end if

  angle = r_pi ( ) / real ( n )
  area = 0.5E+00 * n * radout * radout * sin ( 2.0E+00 * angle )
  side = 2.0E+00 * radout * sin ( angle )
  radin = 0.5E+00 * side / tan ( angle )

  return
end
subroutine polygon_side_data_2d ( area, n, radin, radout, side )
!
!*******************************************************************************
!
!! POLYGON_SIDE_DATA_2D determines polygonal data from its side length in 2D.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real AREA, the area of the regular polygon.
!
!    Input, integer N, the number of sides of the polygon.
!    N must be at least 3.
!
!    Output, real RADIN, the inner radius of the polygon, that is,
!    the radius of the largest circle that can be inscribed within
!    the polygon.
!
!    Output, real RADOUT, the outer radius of the polygon, that is,
!    the radius of the smallest circle that can be described about
!    the polygon.
!
!    Input, real SIDE, the length of one side of the polygon.
!
  implicit none
!
  real angle
  real area
  integer n
  real r_pi
  real radin
  real radout
  real side
!
  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYGON_SIDE_DATA_2D - Fatal error!'
    write ( *, '(a)' ) '  Input value of N must be at least 3'
    write ( *, '(a,i6)' ) '  but your input value was N = ', n
    stop
  end if

  angle = r_pi ( ) / real ( n )
  area = 0.5E+00 * n * side * side / tan ( angle )
  radin = 0.5E+00 * side / tan ( angle )
  radout = 0.5E+00 * side / sin ( angle )

  return
end
subroutine polyhedron_surface_3d ( coord, maxorder, nface, node, &
  n, order, area )
!
!*******************************************************************************
!
!! POLYHEDRON_SURFACE_3D computes the surface area of a polyhedron in 3D.
!
!
!  Restriction:
!
!    The computation is not valid unless the faces of the polyhedron
!    are planar polygons.
!
!  Reference:
!
!    Allen Van Gelder,
!    Efficient Computation of Polygon Area and Polyhedron Volume,
!    Graphics Gems V, edited by Alan Paeth,
!    AP Professional, 1995.
!
!  Modified:
!
!    18 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COORD(3,N), the 3D coordinates of the vertices.
!    The vertices may be listed in any order.
!
!    Input, integer MAXORDER, the maximum number of vertices that make
!    up a face of the polyhedron.
!
!    Input, integer NFACE, the number of faces of the polyhedron.
!
!    Input, integer NODE(NFACE,MAXORDER).  Face I is defined by
!    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
!    are listed in neighboring order.
!
!    Input, integer N, the number of points stored in COORD.
!
!    Input, integer ORDER(NFACE), the number of vertices making up each face.
!
!    Output, real AREA, the total surface area of the polyhedron.
!
  implicit none
!
  integer maxorder
  integer nface
  integer n
!
  real ainc
  real area
  real coord(3,n)
  integer i
  integer j
  integer k
  integer node(nface,maxorder)
  integer order(nface)
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
  real z1
  real z2
  real z3
  real z4
!
  area = 0.0E+00

  do i = 1, nface

    x4 = 0.0E+00
    y4 = 0.0E+00
    z4 = 0.0E+00
!
!  Compute the area vector for this face.
!
    do j = 1, order(i)

      k = node(i,j)
      x1 = coord(1,k)
      y1 = coord(2,k)
      z1 = coord(3,k)

      if ( j < order(i) ) then
        k = node(i,j+1)
      else
        k = node(i,1)
      end if

      x2 = coord(1,k)
      y2 = coord(2,k)
      z2 = coord(3,k)

      call cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )

      x4 = x4 + x3
      y4 = y4 + y3
      z4 = z4 + z3

    end do
!
!  Add the magnitude of the area vector to the sum.
!
    ainc = sqrt ( x4 * x4 + y4 * y4 + z4 * z4 )
    area = area + ainc

  end do

  area = 0.5E+00 * area

  return
end
subroutine polyhedron_volume_3d ( coord, maxorder, nface, node, &
  n, order, volume )
!
!*******************************************************************************
!
!! POLYHEDRON_VOLUME_3D computes the volume of a polyhedron in 3D.
!
!
!  Restriction:
!
!    The computation is not valid unless the faces of the polyhedron
!    are planar polygons.
!
!  Reference:
!
!    Allen Van Gelder,
!    Efficient Computation of Polygon Area and Polyhedron Volume,
!    Graphics Gems V, edited by Alan Paeth,
!    AP Professional, 1995.
!
!  Modified:
!
!    18 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real COORD(3,N), the 3D coordinates of the vertices.
!    The vertices may be listed in any order.
!
!    Input, integer MAXORDER, the maximum number of vertices that make
!    up a face of the polyhedron.
!
!    Input, integer NFACE, the number of faces of the polyhedron.
!
!    Input, integer NODE(NFACE,MAXORDER).  Face I is defined by
!    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
!    are listed in neighboring order.
!
!    Input, integer N, the number of points stored in COORD.
!
!    Input, integer ORDER(NFACE), the number of vertices making up
!    each face.
!
!    Output, real VOLUME, the volume of the polyhedron.
!
  implicit none
!
  integer maxorder
  integer nface
  integer n
!
  real coord(3,n)
  integer i
  integer j
  integer k
  integer node(nface,maxorder)
  integer order(nface)
  real volume
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
  real z1
  real z2
  real z3
  real z4
!
  volume = 0.0E+00

  do i = 1, nface

    x4 = 0.0E+00
    y4 = 0.0E+00
    z4 = 0.0E+00
!
!  Compute the area vector for this face.
!
    do j = 1, order(i)

      k = node(i,j)
      x1 = coord(1,k)
      y1 = coord(2,k)
      z1 = coord(3,k)

      if ( j < order(i) ) then
        k = node(i,j+1)
      else
        k = node(i,1)
      end if

      x2 = coord(1,k)
      y2 = coord(2,k)
      z2 = coord(3,k)

      call cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )

      x4 = x4 + x3
      y4 = y4 + y3
      z4 = z4 + z3

    end do
!
!  Area vector dot any vertex.
!
    k = node(i,1)
    volume = volume + x4 * coord(1,k) + y4 * coord(2,k) + z4 * coord(3,k)

  end do

  volume = volume / 6.0E+00

  return
end
subroutine polyline_index_point_nd ( maxpts, ndim, npts, t, xpts, x )
!
!*******************************************************************************
!
!! POLYLINE_INDEX_POINT_ND evaluates a polyline at a given arclength in ND.
!
!
!  Discussion:
!
!    The polyline is defined as the set of M-1 line segments lying
!    between a sequence of M points.  The arclength of a point lying
!    on the polyline is simply the length of the broken line from the
!    initial point.  Any point on the polyline can be found by
!    specifying its arclength.
!
!    If the given arclength coordinate is less than 0, or greater
!    than the arclength coordinate of the last given point, then
!    extrapolation is used, that is, the first and last line segments
!    are extended as necessary.
!
!    The arclength coordinate system measures the distance between
!    any two points on the polyline as the length of the segment of the
!    line that joins them.
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXPTS, the first dimension of the array
!    used to hold XPTS.  MAXPTS must be at least equal to NPTS.
!
!    Input, integer NDIM, the dimension of the space in which
!    the points lie.  The second dimension of XPTS.
!
!    Input, integer NPTS, the number of points.
!
!    Input, real T, the desired arclength coordinate.
!
!    Input, real XPTS(MAXPTS,NDIM), a set of NPTS coordinates
!    in NDIM space, describing a set of points that define
!    a polyline.
!
!    Output, real X(NDIM), a point lying on the polyline defined
!    by XPTS, and having arclength coordinate T.
!
  implicit none
!
  integer maxpts
  integer ndim
!
  integer i
  integer j
  integer npts
  real s
  real t
  real tleft
  real trite
  real xpts(maxpts,ndim)
  real x(ndim)
!
  if ( maxpts <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYLINE_INDEX_POINT_ND - Fatal error!'
    write ( *, '(a)' ) '  The input quantity MAXPTS is nonpositive.'
    write ( *, '(a,i6)' ) '  MAXPTS = ', maxpts
    stop
  end if

  if ( npts <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYLINE_INDEX_POINT_ND - Fatal error!'
    write ( *, '(a)' ) '  The input quantity NPTS is nonpositive.'
    write ( *, '(a,i6)' ) '  NPTS = ', npts
    stop
  end if

  if ( maxpts < npts ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLYLINE_INDEX_POINT_ND - Fatal error!'
    write ( *, '(a)' ) '  MAXPTS is less than NPTS.'
    write ( *, '(a,i6)' ) '  MAXPTS = ', maxpts
    write ( *, '(a,i6)' ) '  NPTS = ', npts
    stop
  end if

  if ( npts == 1 ) then

    x(1:ndim) = xpts(1,1:ndim)

  else

    trite = 0.0E+00
    do i = 1, npts - 1
!
!  Find the distance between points I and I+1.
!
      tleft = trite
      trite = trite &
        + sqrt ( sum ( ( xpts(i+1,1:ndim) - xpts(i,1:ndim) )**2 ) )
!
!  Interpolate or extrapolate in an interval.
!
      if ( t <= trite .or. i == npts - 1 ) then

        s = ( t - tleft ) / ( trite - tleft )
        x(1:ndim) = ( 1.0E+00 - s ) * xpts(i,1:ndim) + s * xpts(i+1,1:ndim)

        return
      end if
    end do
  end if

  return
end
subroutine polyline_length_nd ( maxpts, ndim, npoint, spoint, length, xpoint )
!
!*******************************************************************************
!
!! POLYLINE_LENGTH_ND computes the length of a polyline in ND.
!
!
!  Definition:
!
!    A polyline of order M is the geometric structure consisting of
!    the M-1 line segments that lie between successive elements of a list
!    of M points.
!
!    An ordinary line segment is a polyline of order 2.
!    The letter "V" is a polyline of order 3.
!    The letter "N" is a polyline of order 4, and so on.
!
!  Formula:
!
!    DIST(I+1,I) = sqrt ( SUM(j = 1,NDIM) ( X(I+1) - X(I) )**2 )
!
!    LENGTH = SUM (I = 1 to NPOINT-1) DIST(I+1,I)
!
!  Modified:
!
!    18 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXPTS, the declared first dimension of XPOINT.
!
!    Input, integer NDIM, the number of dimensions of the points.
!
!    Input, integer NPOINT, the number of points.
!
!    Output, real SPOINT(MAXPTS), the arclength coordinates
!    of each point.  The first point has SPOINT(1) = 0.
!
!    Output, real LENGTH, the length of the polyline,
!    which is also equal to SPOINT(NPOINT).
!
!    Input, real XPOINT(MAXPTS,NDIM), the coordinates of the points.
!
  implicit none
!
  integer maxpts
  integer ndim
!
  integer i
  real length
  integer npoint
  real spoint(maxpts)
  real temp
  real xpoint(maxpts,ndim)
!
  spoint(1) = 0.0E+00

  do i = 2, npoint

    spoint(i) = spoint(i-1) &
      + sqrt ( sum ( ( xpoint(i,1:ndim) - xpoint(i-1,1:ndim) )**2 ) )

  end do

  length = spoint(npoint)

  return
end
subroutine proplane2 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, npnt, xp, yp, zp, &
  alpha, beta )
!
!*******************************************************************************
!
!! PROPLANE2 produces 2D coordinates of points that lie in a plane, in 3D.
!
!
!  Discussion:
!
!    The plane is specified by three non-colinear points, which we will
!    call P1, P2 and P3.
!
!    The first thing to do is to compute two orthonormal vectors V1 and
!    V2, so that any point P that lies in the plane may be written as
!
!      P = P1 + alpha * V1 + beta * V2
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
!    coordinates of three points on the plane.  These three points
!    should not lie in a straight line, but this condition is not
!    checked.
!
!    Input, integer NPNT, the number of points to project.
!
!    Input, real XP(NPNT), YP(NPNT), ZP(NPNT), are the Cartesian
!    coordinates of points which lie on the plane spanned by the
!    three points.  These points are not checked to ensure that
!    they lie on the plane.
!
!    Output, real ALPHA(NPNT), BETA(NPNT), the "in-plane" coordinates of
!    the points.
!
  implicit none
!
  integer npnt
!
  real alpha(npnt)
  real beta(npnt)
  real dot
  integer i
  real v1(3)
  real v2(3)
  real x1
  real x2
  real x3
  real xp(npnt)
  real y1
  real y2
  real y3
  real yp(npnt)
  real z1
  real z2
  real z3
  real zp(npnt)
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1) = x2 - x1
  v1(2) = y2 - y1
  v1(3) = z2 - z1

  call vector_unit_nd ( 3, v1 )

  v2(1) = x3 - x1
  v2(2) = y3 - y1
  v2(3) = z3 - z1

  dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

  v2(1:3) = v2(1:3) - dot * v1(1:3)

  call vector_unit_nd ( 3, v2 )
!
!  Now decompose each (X,Y,Z).
!
  do i = 1, npnt

    alpha(i) = ( xp(i) - x1 ) * v1(1) +  ( yp(i) - y1 ) * v1(2) + &
               ( zp(i) - z1 ) * v1(3)

    beta(i) =  ( xp(i) - x1 ) * v2(1) + ( yp(i) - y1 ) * v2(2) + &
               ( zp(i) - z1 ) * v2(3)

  end do

  return
end
subroutine proplane3 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  npnt, xo, yo, zo, xp, yp, zp )
!
!*******************************************************************************
!
!! PROPLANE3 projects points orthographically onto a plane, in 3D.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
!    coordinates of three points on the plane.
!
!    Input, integer NPNT, the number of points to project.
!
!    Input, real XO(NPNT), YO(NPNT), ZO(NPNT), are the coordinates of
!    the object points.
!
!    Output, real XP(NPNT), YP(NPNT), ZP(NPNT), are the coordinates of
!    the projections of the object points through the focus point onto
!    the plane.
!
!    XP, YP, and ZP may share the same memory as XO, YO, and ZO, in
!    which case the projections will overwrite the original data.
!
  implicit none
!
  integer npnt
!
  real a
  real b
  real c
  real d
  integer i
  real x1
  real x2
  real x3
  real xo(npnt)
  real xp(npnt)
  real y1
  real y2
  real y3
  real yo(npnt)
  real yp(npnt)
  real z1
  real z2
  real z3
  real zo(npnt)
  real zp(npnt)
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!  For each point, its image in the plane is the nearest point
!  in the plane.
!
  do i = 1, npnt

    call plane_imp_point_near_3d ( a, b, c, d, xo(i), yo(i), &
      zo(i), xp(i), yp(i), zp(i) )

  end do

  return
end
subroutine provec ( base, m, n, vecm, vecn, vecnm )
!
!*******************************************************************************
!
!! PROVEC projects a vector from M space into N space.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real BASE(M,N).  The columns of BASE contain
!    N vectors, each of length M, which form the basis for
!    a space of dimension N.
!
!    Input, integer M, the dimension of the higher order space.
!
!    Input, integer N, the dimension of the lower order space.
!
!    Input, real VECM(M), is an M dimensional vector.
!
!    Output, real VECN(N), the projection of VECM into the
!    lower dimensional space.  These values represent
!    coordinates in the lower order space.
!
!    Output, real VECNM(M), the projection of VECM into the
!    lower dimensional space, but using coordinates in
!    the higher dimensional space.
!
  implicit none
!
  integer m
  integer n
!
  real base(m,n)
  integer i
  integer j
  integer k
  real temp
  real vecm(m)
  real vecn(n)
  real vecnm(m)
!
!  For each vector, remove all projections onto previous vectors,
!  and then normalize.  This should result in a matrix BASE
!  whose columns are orthonormal.
!
  do j = 1, n

    do i = 1, j-1

      temp = dot_product ( base(1:m,i), base(1:m,j) )

      base(1:m,j) = base(1:m,j) - temp * base(1:m,i)

    end do

    temp = sqrt ( sum ( base(1:m,j)**2 ) )

    if ( temp > 0.0E+00 ) then
      base(1:m,j) = base(1:m,j) / temp
    end if

  end do
!
!  Compute the coordinates of the projection of the vector
!  simply by taking dot products.
!
  do j = 1, n
    vecn(j) = dot_product ( vecm(1:m), base(1:m,j) )
  end do
!
!  Compute the coordinates of the projection in terms of
!  the original space.
!
  do i = 1, m
    vecnm(i) = dot_product ( base(i,1:n), vecn(1:n) )
  end do

  return
end
subroutine pyramid_volume_3d ( h, s, volume )
!
!*******************************************************************************
!
!! PYRAMID_VOLUME_3D computes the volume of a pyramid with square base in 3D.
!
!
!  Modified:
!
!    10 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real H, R, the height of the pyramid, and the length of one
!    side of the square base.
!
!    Output, real VOLUME, the volume of the pyramid.
!
  implicit none
!
  real h
  real s
  real volume
!
  volume = s * s * h / 3.0E+00

  return
end
subroutine quad_area_2d ( x1, y1, x2, y2, x3, y3, x4, y4, area )
!
!*******************************************************************************
!
!! QUAD_AREA_2D computes the area of a quadrilateral in 2D.
!
!
!  Discussion:
!
!    This algorithm should be able to handle nonconvex quadrilaterals.
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, X4, Y4, the X and Y coordinates
!    of the corners of the quadrilateral.  The corners should be
!    specified in clockwise or counterclockwise order.
!
!    Output, real AREA, the absolute area of the quadrilateral.
!
  implicit none
!
  real area
  real area1
  real area2
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
 !
  call triangle_area_2d ( x1, y1, x2, y2, x3, y3, area1 )

  call triangle_area_2d ( x3, y3, x4, y4, x1, y1, area2 )

  area = area1 + area2

  return
end
subroutine quad_contains_point_2d ( x1, y1, x2, y2, x3, y3, x4, y4, x, y, &
  inside )
!
!*******************************************************************************
!
!! QUAD_CONTAINS_POINT_2D finds if a point is inside a convex quadrilateral in 2D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, X4, Y4, the X and Y coordinates
!    of the quadrilateral.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, logical INSIDE, is .TRUE. if ( X,Y) is inside
!    the quadrilateral or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  real angle_123
  real angle_12x
  real angle_234
  real angle_23x
  real angle_341
  real angle_34x
  real angle_412
  real angle_41x
  real anglei_rad_2d
  logical inside
  real x
  real x1
  real x2
  real x3
  real x4
  real y
  real y1
  real y2
  real y3
  real y4
!
!  This will only handle convex quadrilaterals.
!
  inside = .false.

  angle_12x = anglei_rad_2d ( x1, y1, x2, y2, x, y )
  angle_123 = anglei_rad_2d ( x1, y1, x2, y2, x3, y3 )

  if ( angle_12x > angle_123 ) then
    return
  end if

  angle_23x = anglei_rad_2d ( x2, y2, x3, y3, x, y )
  angle_234 = anglei_rad_2d ( x2, y2, x3, y3, x4, y4 )

  if ( angle_23x > angle_234 ) then
    return
  end if

  angle_34x = anglei_rad_2d ( x3, y3, x4, y4, x, y )
  angle_341 = anglei_rad_2d ( x3, y3, x4, y4, x1, y1 )

  if ( angle_34x > angle_341 ) then
    return
  end if

  angle_41x = anglei_rad_2d ( x4, y4, x1, y1, x, y )
  angle_412 = anglei_rad_2d ( x4, y4, x1, y1, x2, y2 )

  if ( angle_41x > angle_412 ) then
    return
  end if

  inside = .true.

  return
end
subroutine quad_point_dist_2d ( x1, y1, x2, y2, x3, y3, x4, y4, x, y, dist )
!
!*******************************************************************************
!
!! QUAD_POINT_DIST_2D finds the distance from a point to a quadrilateral in 2D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, X4, Y4, the corners of the
!    quadrilateral.
!
!    Input, real X, Y, the point which is to be checked.
!
!    Output, real DIST, the distance from the point to the quadrilateral.
!    DIST is zero if the point lies exactly on the quadrilateral.
!
  implicit none
!
  real dist
  real dist2
  real enorm0_2d
  real x
  real x1
  real x2
  real x3
  real x4
  real y
  real y1
  real y2
  real y3
  real y4
!
  dist =             enorm0_2d ( x, y, x1, y1 )
  dist = min ( dist, enorm0_2d ( x, y, x2, y2 ) )
  dist = min ( dist, enorm0_2d ( x, y, x3, y3 ) )
  dist = min ( dist, enorm0_2d ( x, y, x4, y4 ) )

  call line_seg_point_dist_2d ( x1, y1, x2, y2, x, y, dist2 )
  dist = min ( dist, dist2 )

  call line_seg_point_dist_2d ( x2, y2, x3, y3, x, y, dist2 )
  dist = min ( dist, dist2 )

  call line_seg_point_dist_2d ( x3, y3, x4, y4, x, y, dist2 )
  dist = min ( dist, dist2 )

  call line_seg_point_dist_2d ( x4, y4, x1, y1, x, y, dist2 )
  dist = min ( dist, dist2 )

  return
end
subroutine quad_point_dist_signed_2d ( x1, y1, x2, y2, x3, y3, &
  x4, y4, x, y, dist_signed )
!
!*******************************************************************************
!
!! QUAD_POINT_DIST_SIGNED_2D: signed distanct ( quadrilateral, point ) in 2D.
!
!
!  Discussion:
!
!    The quadrilateral must be convex.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, X4, Y4, the corners of the
!    quadrilateral, given in order.  It is possible to list the points
!    in such an order that the quadrilateral is not convex.  This
!    routine does not check for convexity, and will produce erroneous
!    results in such a case.
!
!    Input, real X, Y, the point which is to be checked.
!
!    Output, real DIST_SIGNED, the signed distance from the point to the
!    convex quadrilateral.  DIST_SIGNED is actually the maximum of the
!    signed distances from the point to each of the four lines that make
!    up the quadrilateral.
!
!    (Essentially, if the point is outside the convex quadrilateral,
!    only one of the signed distances can be positive, or two can
!    be positive and equal.)
!
!    If DIST_SIGNED is
!    0, the point is on the boundary;
!    negative, the point is in the interior;
!    positive, the point is in the exterior.
!
  implicit none
!
  real dis
  real dis12
  real dis23
  real dis34
  real dis41
  real dist_signed
  real x
  real x1
  real x2
  real x3
  real x4
  real xm
  real y
  real y1
  real y2
  real y3
  real y4
  real ym
!
!  Compare the signed distance from each line segment to the point,
!  with the signed distance to the midpoint of the opposite line.
!
!  The signed distances should all be negative if the point is inside
!  the triangle.
!
  call line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dis12 )

  xm = 0.5E+00 * ( x3 + x4 )
  ym = 0.5E+00 * ( y3 + y4 )

  call line_exp_point_dist_signed_2d ( x1, y1, x2, y2, xm, ym, dis )

  if ( dis > 0.0E+00 ) then
    dis = - dis
    dis12 = - dis12
  end if

  call line_exp_point_dist_signed_2d ( x2, y2, x3, y3, x, y, dis23 )

  xm = 0.5E+00 * ( x4 + x1 )
  ym = 0.5E+00 * ( y4 + y1 )

  call line_exp_point_dist_signed_2d ( x2, y2, x3, y3, xm, ym, dis )

  if ( dis > 0.0E+00 ) then
    dis = - dis
    dis23 = - dis23
  end if

  call line_exp_point_dist_signed_2d ( x3, y3, x4, y4, x, y, dis34 )

  xm = 0.5E+00 * ( x1 + x2 )
  ym = 0.5E+00 * ( y1 + y2 )

  call line_exp_point_dist_signed_2d ( x3, y3, x4, y4, xm, ym, dis )

  if ( dis > 0.0E+00 ) then
    dis = - dis
    dis34 = - dis34
  end if

  call line_exp_point_dist_signed_2d ( x4, y4, x1, y1, x, y, dis41 )

  xm = 0.5E+00 * ( x2 + x3 )
  ym = 0.5E+00 * ( y2 + y3 )

  call line_exp_point_dist_signed_2d ( x4, y4, x1, y1, xm, ym, dis )

  if ( dis > 0.0E+00 ) then
    dis = - dis
    dis41 = - dis41
  end if

  dist_signed = max ( dis12, dis23, dis34, dis41 )

  return
end
subroutine quat_conj ( q )
!
!*******************************************************************************
!
!! QUAT_CONJ conjugates a quaternion.
!
!
!  Discussion:
!
!    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
!    may be written as
!
!      Q = A + Bi + Cj + Dk.
!
!    The conjugate of Q is
!
!      conj ( Q ) = A - Bi - Cj - Dk.
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real Q(4).  On input, the quaternion to be conjugated.
!    On output, the conjugated quaternion.
!
  implicit none
!
  real q(4)
!
  q(2) = - q(2)
  q(3) = - q(3)
  q(4) = - q(4)

  return
end
subroutine quat_inv ( q )
!
!*******************************************************************************
!
!! QUAT_INV inverts a quaternion.
!
!
!  Discussion:
!
!    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
!    may be written as
!
!      Q = A + Bi + Cj + Dk.
!
!    The inverse of Q is
!
!      inverse ( Q ) = conjugate ( Q ) / ( norm ( Q ) )**2.
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real Q(4).  On input, the quaternion to be inverted.
!    On output, the inverse of the input quaternion.
!
  implicit none
!
  real q(4)
!
  q(1:4) = q(1:4) / sum ( q(1:4)**2 ) 
  q(2:4) = - q(2:4)

  return
end
subroutine quat_mul ( q1, q2, q3 )
!
!*******************************************************************************
!
!! QUAT_MUL multiplies two quaternions.
!
!
!  Discussion:
!
!    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
!    may be written as
!
!      Q = A + Bi + Cj + Dk.
!
!    To multiply two quaternions, use the relationships:
!
!      ij = -ji = k
!      jk = -kj = i
!      ki = -ik = j
!      ii = jj = kk = -1
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Q1(4), Q2(4), the two quaternions to be multiplied.
!
!    Output, real Q3(4), the product of the two quaternions.
!
  implicit none
!
  real q1(4)
  real q2(4)
  real q3(4)
!
  q3(1) = q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3) - q1(4) * q2(4)
  q3(2) = q1(1) * q2(2) + q1(2) * q2(1) + q1(3) * q2(4) - q1(4) * q2(3)
  q3(3) = q1(1) * q2(3) - q1(2) * q2(4) + q1(3) * q2(1) + q1(4) * q2(2)
  q3(4) = q1(1) * q2(4) + q1(2) * q2(3) - q1(3) * q2(2) + q1(4) * q2(1)

  return
end
function quat_norm ( q )
!
!*******************************************************************************
!
!! QUAT_NORM computes the norm of a quaternion.
!
!
!  Discussion:
!
!    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
!    may be written as
!
!      Q = A + Bi + Cj + Dk.
!
!    The norm of Q is
!
!      norm(Q) = sqrt ( A**2 + B**2 + C**2 + D**2 ).
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Q(4), the quaternion.
!
!    Output, real QUAT_NORM, the norm of the quaternion.
!
  implicit none
!
  real q(4)
  real quat_norm
!
  quat_norm = sqrt ( sum ( q(1:4)**2 ) )

  return
end
function r_modp ( x, y )
!
!*******************************************************************************
!
!! R_MODP returns the nonnegative remainder of real division.
!
!
!  Formula:
!
!    If
!      REM = R_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R_MODP(A,360.0) is between 0 and 360, always.
!
!  Examples:
!
!        I         J     MOD  R_MODP   R_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    24 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number to be divided.
!
!    Input, real Y, the number that divides X.
!
!    Output, real R_MODP, the nonnegative remainder when X is divided by Y.
!
  implicit none
!
  real r_modp
  real x
  real y
!
  if ( y == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r_modp = mod ( x, y )

  if ( r_modp < 0.0E+00 ) then
    r_modp = r_modp + abs ( y )
  end if

  return
end
function r_pi ( )
!
!*******************************************************************************
!
!! R_PI returns the value of pi.
!
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real R_PI, the value of pi.
!
  implicit none
!
  real r_pi
!
  r_pi = 3.14159265358979323846264338327950288419716939937510E+00

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Output, real R, the randomly chosen value.
!
  implicit none
!
  logical, save :: seed = .false.
  real r
  real rhi
  real rlo
  real t
!
  if ( .not. seed ) then
    call random_seed ( )
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine r_to_bin_even2 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R_TO_BIN_EVEN2 determines the appropriate "bin" for C in [A,B].
!
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A = 5, B = 15
!
!    <-1-+-2-+-3-+-4-+-5->
!    5   7   9  11  13  15
!
!
!    C   BIN
!
!    1    1
!    3    1
!    4.9  1
!    5    1
!    6    1
!    7.1  2
!    8    2
!    9.5  3
!   12    4
!   14    5
!   15    5
!   15.1  5
!   99    5
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, real A, B, the lower and upper limits of the bin interval.
!    While A is expected to be less than B, the code should return useful
!    results if A is actually greater than B.
!
!    Input, real C, a value to be placed in a bin.
!
!    Output, integer BIN, the index of the bin to which C is assigned.
!
  implicit none
!
  real a
  real a2
  real b
  real b2
  integer bin
  real c
  integer nbin
  logical switch
!
!  Take care of special cases.
!
  if ( nbin < 1 ) then
    bin = 0
    return
  end if

  if ( nbin == 1 ) then
    bin = 1
    return
  end if

  if ( b == a ) then
    bin = 1
    return
  end if
!
!  If the limits are descending, then we switch them now, and
!  unswitch the results at the end.
!
  if ( a < b ) then
    switch = .false.
    a2 = a
    b2 = b
  else
    switch = .true.
    a2 = b
    b2 = a
  end if
!
!  Compute the bin.
!
  if ( c <= a2 ) then
    bin = 1
  else if ( c >= b2 ) then
    bin = nbin
  else
    bin = 1 + int ( real ( nbin ) * ( c - a2 ) / ( b2 - a2 ) )
    bin = max ( bin, 1 )
    bin = min ( bin, nbin )
  end if
!
!  Reverse the switching.
!
  if ( switch ) then
    bin = nbin + 1 - bin
  end if

  return
end
subroutine r2_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R2_RANDOM returns a random R2 value in a given range.
!
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO(2), RHI(2), the minimum and maximum values.
!
!    Output, real R(2), the randomly chosen value.
!
  implicit none
!
  real r(2)
  real rhi(2)
  real rlo(2)
!
  call random_number ( harvest = r(1:2) )

  r(1:2) = ( 1.0E+00 - r(1:2) ) * rlo(1:2) &
                     + r(1:2)   * rhi(1:2)

  return
end
subroutine r2_to_bin_even2 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R2_TO_BIN_EVEN2 determines the appropriate "bin" for an R2 value.
!
!
!  Discussion:
!
!    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN
!    equal subintervals or bins.  Boundary bins take care of extreme values.
!
!  Example:
!
!    NBIN = 5, A(1) = 5,  A(2) = 0,
!              B(1) = 15, B(2) = 20.
!
!   20 +    +    +    +    +    +
!        15 | 25 | 35 | 45 | 55
!   16 +----+----+----+----+----+
!        14 | 24 | 34 | 44 | 54
!   12 +----+----+----+----+----+
!        13 | 23 | 33 | 43 | 53
!    8 +----+----+----+----+----+
!        12 | 22 | 32 | 42 | 52
!    4 +----+----+----+----+----+
!        11 | 21 | 31 | 41 | 51
!    0 +    +    +    +    +    +
!      5    7    9   11   13   15
!
!      C      BIN
!   ------  ------
!    8 -2    2  1
!    0  1    1  1
!    6  9    1  3
!   10 11    3  3
!   14 23    5  5
!   25 13    5  4
!
!  Modified:
!
!    05 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Input, real C(2), a value to be placed in a bin.
!
!    Output, integer BIN(2), the index of the bin to which C is assigned.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real c(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call r_to_bin_even2 ( nbin, a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r2_to_bin_even3 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R2_TO_BIN_EVEN3 determines the appropriate "bin" for an R2 value.
!
!
!  Discussion:
!
!    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
!    or bins.
!
!  Example:
!
!    NBIN = (/ 4, 5 /),
!
!      A(1) = 1,  A(2) = 0,
!      B(1) = 17, B(2) = 20.
!
!   20 +    +    +    +    +
!        15 | 25 | 35 | 45
!   16 +----+----+----+----+
!        14 | 24 | 34 | 44
!   12 +----+----+----+----+
!        13 | 23 | 33 | 43
!    8 +----+----+----+----+
!        12 | 22 | 32 | 42
!    4 +----+----+----+----+
!        11 | 21 | 31 | 41
!    0 +    +    +    +    +
!      1    5    9   13   17
!
!      C      BIN
!   ------  ------
!    8 -2    2  1
!    0  1    1  1
!    6  9    2  3
!   10 11    3  3
!   14 23    4  5
!   25 13    4  4
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!
!    Input, real A(2), B(2), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Input, real C(2), a value to be placed in a bin.
!
!    Output, integer BIN(2), the index of the bin to which C is assigned.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real c(ndim)
  integer i
  integer nbin(ndim)
!
  do i = 1, ndim
    call r_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r2vec_bin_even2 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BIN_EVEN2 bins an R2 array into evenly spaced bins.
!
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN 1D bins in both X and Y directions, making a total
!    of NBIN**2 2D bins.  Each set of 1D bins begins and ends at user
!    specified mininum and maximum values.
!
!    The 2D bins are indexed by the X and Y bins that construct them,
!    and ordered lexicographically by these indices:
!
!      1,4 | 2,4 | 3,4 | 4,4
!      ----+-----+-----+----
!      1,3 | 2,3 | 3,3 | 4,3
!      ----+-----+-----+----
!      1,2 | 2,2 | 3,2 | 4,2
!      ----+-----+-----+----
!      1,1 | 2,1 | 3,1 | 4,1
!
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (4,4).
!
!  Modified:
!
!    09 February 2001
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real A(2,N), the R2 data to be binned.
!
!    Input, integer NBIN, the (square root of) the number of bins.
!    NBIN must be at least 1.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the
!    index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, integer BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none
!
  integer n
  integer nbin
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer i1
  integer i2
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin,1:nbin) = -1
  bin_start(1:nbin,1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r2_to_bin_even2 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)

    if ( bin_start(i1,i2) == -1 ) then
      bin_start(i1,i2) = j
    else
      k = bin_last(i1,i2)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2) = j

  end do

  return
end
subroutine r2vec_bin_even3 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BIN_EVEN3 bins an R2 array into evenly spaced bins.
!
!
!  Discussion:
!
!    A different number of bins may be used in each dimension.
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN(1) 1D bins in the X direction, NBIN(2) for Y, making a
!    total of NBIN(1) * NBIN(2) 2D bins.  Each set of 1D bins begins and
!    ends at user specified mininum and maximum values.
!
!    The 2D bins are indexed by the X and Y bins that construct them,
!    and ordered lexicographically by these indices:
!
!      1,4 | 2,4 | 3,4 | 4,4 | 5,4
!      ----+-----+-----+-----+-----
!      1,3 | 2,3 | 3,3 | 4,3 | 5,3
!      ----+-----+-----+-----+-----
!      1,2 | 2,2 | 3,2 | 4,2 | 5,2
!      ----+-----+-----+-----+-----
!      1,1 | 2,1 | 3,1 | 4,1 | 5,1
!
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (5,4).
!
!  Modified:
!
!    26 March 2001
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real A(2,N), the R2 data to be binned.
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, integer  BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    the index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, integer BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  integer n
  integer nbin(ndim)
!
  real a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2))
  real bin_max(ndim)
  real bin_min(ndim)
  integer i1
  integer i2
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin(1),1:nbin(2)) = -1
  bin_start(1:nbin(1),1:nbin(2)) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r2_to_bin_even3 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)

    if ( bin_start(i1,i2) == -1 ) then
      bin_start(i1,i2) = j
    else
      k = bin_last(i1,i2)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2) = j

  end do

  return
end
subroutine r2vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BINNED_REORDER reorders a binned R2 data vector.
!
!
!  Discussion:
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN, the (square root of the) number of bins.
!
!    Input/output, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the
!    index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Input/output, integer BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none
!
  integer n
  integer nbin
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  real a2(ndim,n)
  integer bin_last(nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin)
  real bin_max
  real bin_min
  integer i
  integer i1
  integer i2
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin

    do i2 = 1, nbin

      j = bin_start(i1,i2)

      if ( j > 0 ) then
        bin_start(i1,i2) = k + 1
      end if

      do while ( j > 0 )
        k = k + 1
        bin_last(i1,i2) = k
        a2(1:ndim,k) = a(1:ndim,j)
        j = bin_next(j)
      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin
    do i2 = 1, nbin

      k = bin_last(i1,i2)

      if ( k > 0 ) then
        bin_next(k) = 0
      end if

    end do
  end do

  return
end
subroutine r2vec_binned_reorder2 ( n, a, nbin, bin_start, bin_last, bin_next )
!
!*******************************************************************************
!
!! R2VEC_BINNED_REORDER2 reorders a binned R2 data vector.
!
!
!  Discussion:
!
!    This routine allows there to be a different number of bins in
!    each dimension.
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN(2), the number of bins in each direction.
!
!    Input/output, integer BIN_START(NBIN(1),NBIN(2)),
!    BIN_LAST(NBIN(1),NBIN(2)), the index of the first and last element of A
!    that went into each bin, or -1 if there are no entries in the bin.
!
!    Input/output, integer BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  integer n
  integer nbin(ndim)
!
  real a(ndim,n)
  real a2(ndim,n)
  integer bin_last(nbin(1),nbin(2))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2))
  real bin_max
  real bin_min
  integer i
  integer i1
  integer i2
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      j = bin_start(i1,i2)

      if ( j > 0 ) then
        bin_start(i1,i2) = k + 1
      end if

      do while ( j > 0 )
        k = k + 1
        bin_last(i1,i2) = k
        a2(1:ndim,k) = a(1:ndim,j)
        j = bin_next(j)
      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin(1)
    do i2 = 1, nbin(2)

      k = bin_last(i1,i2)

      if ( k > 0 ) then
        bin_next(k) = 0
      end if

    end do
  end do

  return
end
subroutine r2vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )
!
!*******************************************************************************
!
!! R2VEC_BINNED_SORT_A sorts each bin of an R2 binned data vector.
!
!
!  Discussion:
!
!    Presumably, the data vector was first binned by R2VEC_BIN_EVEN,
!    then reordered by R2VEC_BINNED_REORDER.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R2 vector.
!
!    What is true is that if I < J, then either the I-th element of A occurs
!    in a lexicographically smaller bin than J, or they share a bin,
!    and the I-th element is lexicographically less than or equal to
!    the J-th element.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN, the (square root of the) number of bins.
!
!    Input, integer BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the index
!    of the first and last element of A that went into each bin, or -1
!    if there are no entries in this bin.
!
  implicit none
!
  integer n
  integer nbin
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer bin_last(nbin,nbin)
  integer bin_start(nbin,nbin)
  integer i1
  integer i2
  integer j1
  integer j2
  integer n1
!
  do i1 = 1, nbin

    do i2 = 1, nbin

      j1 = bin_start(i1,i2)

      if ( j1 > 0 ) then

        j2 = bin_last(i1,i2)

        n1 = j2 + 1 - j1

        if ( n1 > 1 ) then
          call r2vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
        end if

      end if

    end do

  end do

  return
end
subroutine r2vec_binned_sort_a2 ( n, a, nbin, bin_start, bin_last )
!
!*******************************************************************************
!
!! R2VEC_BINNED_SORT_A2 sorts each bin of an R2 binned data vector.
!
!
!  Discussion:
!
!    This routine allows a different number of bins in each dimension.
!
!    Presumably, the data vector was first binned by R2VEC_BIN_EVEN3,
!    then reordered by R2VEC_BINNED_REORDER2.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R2 vector.
!
!    What is true is that if I < J, then either the I-th element of A occurs
!    in a lexicographically smaller bin than J, or they share a bin,
!    and the I-th element is lexicographically less than or equal to
!    the J-th element.
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    the index of the first and last element of A that went into each bin, or -1
!    if there are no entries in this bin.
!
  implicit none
!
  integer, parameter :: ndim = 2
!
  integer n
  integer nbin(ndim)
!
  real a(ndim,n)
  integer bin_last(nbin(1),nbin(2))
  integer bin_start(nbin(1),nbin(2))
  integer i1
  integer i2
  integer j1
  integer j2
  integer n1
!
  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      j1 = bin_start(i1,i2)

      if ( j1 > 0 ) then

        j2 = bin_last(i1,i2)

        n1 = j2 + 1 - j1

        if ( n1 > 1 ) then
          call r2vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
        end if

      end if

    end do

  end do

  return
end
subroutine r2vec_part_quick_a ( n, a, l, r )
!
!*******************************************************************************
!
!! R2VEC_PART_QUICK_A reorders an R2 vector as part of a quick sort.
!
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1:2,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
!             -----------          ----------------------------------
!             LEFT          KEY    RIGHT
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of A.
!
!    Input/output, real A(2,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1:2,1).  Then
!    I <= L                 A(1:2,I) < KEY;
!         L < I < R         A(1:2,I) = KEY;
!                 R <= I    A(1:2,I) > KEY.
!
  implicit none
!
  integer n
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer i
  real key(ndim)
  integer l
  integer m
  integer r
  logical rvec_eq
  logical rvec_gt
  logical rvec_lt
  real temp
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R2VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:ndim) = a(1:ndim,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( rvec_gt ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      r = r - 1
      call rvec_swap ( ndim, a(1:ndim,r), a(1:ndim,l+1) )
    else if ( rvec_eq ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      m = m + 1
      call rvec_swap ( ndim, a(1:ndim,m), a(1:ndim,l+1) )
      l = l + 1
    else if ( rvec_lt ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:ndim,i) = a(1:ndim,i+m)
  end do

  l = l - m

  do i = 1, ndim
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r2vec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! R2VEC_RANDOM returns a random R2 vector in a given range.
!
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO(2), AHI(2), the minimum and maximum values allowed
!    for A(1,*) and A(2,*).
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real A(2,N), the vector of randomly chosen values.
!
  implicit none
!
  integer n
!
  real a(2,n)
  real ahi(2)
  real alo(2)
  integer i
!
  call random_number ( harvest = a(1:2,1:n) )

  do i = 1, 2
    a(i,1:n) = ( 1.0E+00 - a(i,1:n) ) * alo(i) & 
                         + a(i,1:n)   * ahi(i)
  end do

  return
end
subroutine r2vec_sort_quick_a ( n, a )
!
!*******************************************************************************
!
!! R2VEC_SORT_QUICK_A ascending sorts an R2 vector using quick sort.
!
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(2,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none
!
  integer, parameter :: maxlevel = 25
  integer n
  integer, parameter :: ndim = 2
!
  real a(ndim,n)
  integer base
  integer l_segment
  integer level
  integer n_segment
  integer rsave(maxlevel)
  integer r_segment
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R2VEC_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r2vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( l_segment > 1 ) then

      if ( level > maxlevel ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R2VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i6)' ) '  Exceeding recursion maximum of ', maxlevel
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( n_segment > 0 ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine r3_to_bin_even2 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R3_TO_BIN_EVEN2 determines the appropriate "bin" for an R3 value.
!
!
!  Discussion:
!
!    The intervals [A(I),B(I)] are each divided into NBIN
!    equal subintervals or bins.  Boundary bins take care of extreme values.
!
!  Modified:
!
!    05 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real A(3), B(3), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Input, real C(3), a value to be placed in a bin.
!
!    Output, integer BIN(3), the index of the bin to which C is assigned.
!
  implicit none
!
  integer, parameter :: ndim = 3
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real c(ndim)
  integer i
  integer nbin
!
  do i = 1, ndim
    call r_to_bin_even2 ( nbin, a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r3_to_bin_even3 ( nbin, a, b, c, bin )
!
!*******************************************************************************
!
!! R3_TO_BIN_EVEN3 determines the appropriate "bin" for an R3 value.
!
!
!  Discussion:
!
!    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
!    or bins.
!
!  Example:
!
!    NBIN = (/ 4, 5, 2 /),
!
!      A(1) = 1,  A(2) = 0,  A(3) = 8
!      B(1) = 17, B(2) = 20, B(3) = 10
!
!
!            8 < Z < 9                    9 < Z < 10
!
!   20 +     +     +     +     +     20 +     +     +     +     +
!        151 | 251 | 351 | 451            152 | 252 | 352 | 452
!   16 +-----+-----+-----+-----+     16 +-----+-----+-----+-----+
!        141 | 241 | 341 | 441            142 | 242 | 342 | 442
!   12 +-----+-----+-----+-----+     12 +-----+-----+-----+-----+
!        131 | 231 | 331 | 431            132 | 232 | 332 | 432
!    8 +-----+-----+-----+-----+      8 +-----+-----+-----+-----+
!        121 | 221 | 321 | 421            122 | 222 | 322 | 422
!    4 +-----+-----+-----+-----+      4 +-----+-----+-----+-----+
!        111 | 211 | 311 | 411            112 | 212 | 312 | 412
!    0 +     +     +     +     +      0 +     +     +     +     +
!      1     5     9    13    17        1     5     9    13    17
!
!  Modified:
!
!    05 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NBIN(3), the number of bins in each dimension.
!
!    Input, real A(3), B(3), the lower and upper limits of the bin interval.
!    While A(I) is expected to be less than B(I), the code should return useful
!    results if A(I) is actually greater than B(I).
!
!    Input, real C(3), a value to be placed in a bin.
!
!    Output, integer BIN(3), the index of the bin to which C is assigned.
!
  implicit none
!
  integer, parameter :: ndim = 3
!
  real a(ndim)
  real b(ndim)
  integer bin(ndim)
  real c(ndim)
  integer i
  integer nbin(ndim)
!
  do i = 1, ndim
    call r_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r3vec_bin_even2 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )
!
!*******************************************************************************
!
!! R3VEC_BIN_EVEN2 bins an R3 array into evenly spaced bins.
!
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN 1D bins in each coordinate, making a total
!    of NBIN**NDIM bins.  Each set of 1D bins begins and ends at user
!    specified mininum and maximum values.
!
!    The bins are indexed by the 1D bins that construct them,
!    and ordered lexicographically by these indices:
!
!  Modified:
!
!    09 February 2001
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real A(3,N), the data to be binned.
!
!    Input, integer NBIN, the (cube root of) the number of bins.
!    NBIN must be at least 1.
!
!    Input, real BIN_MIN(3), BIN_MAX(3), the bin limits.
!
!    Output, integer BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
!    the index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, integer BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none
!
  integer n
  integer nbin
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin,nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin,nbin)
  real bin_max(ndim)
  real bin_min(ndim)
  integer i1
  integer i2
  integer i3
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin,1:nbin,1:nbin) = -1
  bin_start(1:nbin,1:nbin,1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r3_to_bin_even2 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)
    i3 = bin(3)

    if ( bin_start(i1,i2,i3) == -1 ) then
      bin_start(i1,i2,i3) = j
    else
      k = bin_last(i1,i2,i3)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2,i3) = j

  end do

  return
end
subroutine r3vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )
!
!*******************************************************************************
!
!! R3VEC_BINNED_REORDER reorders a binned R3 data vector.
!
!
!  Discussion:
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(3,N), the data to be sorted.
!
!    Input, integer NBIN, the (cube root of the) number of bins.
!
!    Input/output, integer BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
!    the index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Input/output, integer BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none
!
  integer n
  integer nbin
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  real a2(ndim,n)
  integer bin_last(nbin,nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin,nbin)
  real bin_max
  real bin_min
  integer i
  integer i1
  integer i2
  integer i3
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin

    do i2 = 1, nbin

      do i3 = 1, nbin

        j = bin_start(i1,i2,i3)

        if ( j > 0 ) then
          bin_start(i1,i2,i3) = k + 1
        end if

        do while ( j > 0 )
          k = k + 1
          bin_last(i1,i2,i3) = k
          a2(1:ndim,k) = a(1:ndim,j)
          j = bin_next(j)
        end do

      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin
    do i2 = 1, nbin
      do i3 = 1, nbin
        k = bin_last(i1,i2,i3)

        if ( k > 0 ) then
          bin_next(k) = 0
        end if

      end do
    end do
  end do

  return
end
subroutine r3vec_binned_reorder2 ( n, a, nbin, bin_start, bin_last, bin_next )
!
!*******************************************************************************
!
!! R3VEC_BINNED_REORDER2 reorders a binned R3 data vector.
!
!
!  Discussion:
!
!    This routine allows there to be a different number of bins in
!    each dimension.
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
!
!  Modified:
!
!    30 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(3,N), the data to be sorted.
!
!    Input, integer NBIN(3), the number of bins in each dimension.
!
!    Input/output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the
!    index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Input/output, integer BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none
!
  integer, parameter :: ndim = 3
!
  integer n
  integer nbin(ndim)
!
  real a(ndim,n)
  real a2(ndim,n)
  integer bin_last(nbin(1),nbin(2),nbin(3))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2),nbin(3))
  real bin_max
  real bin_min
  integer i
  integer i1
  integer i2
  integer i3
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      do i3 = 1, nbin(3)

        j = bin_start(i1,i2,i3)

        if ( j > 0 ) then
          bin_start(i1,i2,i3) = k + 1
        end if

        do while ( j > 0 )
          k = k + 1
          bin_last(i1,i2,i3) = k
          a2(1:ndim,k) = a(1:ndim,j)
          j = bin_next(j)
        end do

      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin(1)
    do i2 = 1, nbin(2)
      do i3 = 1, nbin(3)

        k = bin_last(i1,i2,i3)

        if ( k > 0 ) then
          bin_next(k) = 0
        end if

      end do
    end do
  end do

  return
end
subroutine r3vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )
!
!*******************************************************************************
!
!! R3VEC_BINNED_SORT_A sorts each bin of an R3 binned data vector.
!
!
!  Discussion:
!
!    Presumably, the data vector was first binned by R3VEC_BIN_EVEN,
!    then reordered by R3VEC_BINNED_REORDER.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R3 vector.
!
!    What is true is that if I < J, then either the I-th element of A occurs
!    in a lexicographically smaller bin than J, or they share a bin,
!    and the I-th element is lexicographically less than or equal to
!    the J-th element.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real A(3,N), the data to be sorted.
!
!    Input, integer NBIN, the (cube root of the) number of bins.
!
!    Input, integer BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
!    the index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in this bin.
!
  implicit none
!
  integer n
  integer nbin
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  integer bin_last(nbin,nbin,nbin)
  integer bin_start(nbin,nbin,nbin)
  integer i1
  integer i2
  integer i3
  integer j1
  integer j2
  integer n1
!
  do i1 = 1, nbin

    do i2 = 1, nbin

      do i3 = 1, nbin

        j1 = bin_start(i1,i2,i3)

        if ( j1 > 0 ) then

          j2 = bin_last(i1,i2,i3)

          n1 = j2 + 1 - j1

          if ( n1 > 1 ) then
            call r3vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
          end if

        end if

      end do

    end do

  end do

  return
end
subroutine r3vec_part_quick_a ( n, a, l, r )
!
!*******************************************************************************
!
!! R3VEC_PART_QUICK_A reorders an R3 vector as part of a quick sort.
!
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1:3,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of A.
!
!    Input/output, real A(3,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1:3,1).  Then
!    I <= L                 A(1:3,I) < KEY;
!         L < I < R         A(1:3,I) = KEY;
!                 R <= I    A(1:3,I) > KEY.
!
  implicit none
!
  integer n
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  integer i
  real key(ndim)
  integer l
  integer m
  integer r
  logical rvec_eq
  logical rvec_gt
  logical rvec_lt
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R3VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:ndim) = a(1:ndim,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( rvec_gt ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      r = r - 1
      call rvec_swap ( ndim, a(1:ndim,r), a(1:ndim,l+1) )
    else if ( rvec_eq ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      m = m + 1
      call rvec_swap ( ndim, a(1:ndim,m), a(1:ndim,l+1) )
      l = l + 1
    else if ( rvec_lt ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:ndim,i) = a(1:ndim,i+m)
  end do

  l = l - m

  do i = 1, ndim
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r3vec_random ( alo, ahi, n, a )
!
!*******************************************************************************
!
!! R3VEC_RANDOM returns a random R3 vector in a given range.
!
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO(3), AHI(3), the minimum and maximum values.
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real A(3,N), the vector of randomly chosen values.
!
  implicit none
!
  integer n
!
  real a(3,n)
  real ahi(3)
  real alo(3)
  integer i
!
  call random_number ( harvest = a(1:3,1:n) )

  do i = 1, 3
    a(i,1:n) = ( 1.0E+00 - a(i,1:n) ) * alo(i) & 
                         + a(i,1:n)   * ahi(i)
  end do

  return
end
subroutine r3vec_sort_quick_a ( n, a )
!
!*******************************************************************************
!
!! R3VEC_SORT_QUICK_A ascending sorts an R3 vector using quick sort.
!
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(3,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none
!
  integer, parameter :: maxlevel = 25
  integer n
  integer, parameter :: ndim = 3
!
  real a(ndim,n)
  integer base
  integer l_segment
  integer level
  integer n_segment
  integer rsave(maxlevel)
  integer r_segment
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R3VEC_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r3vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( l_segment > 1 ) then

      if ( level > maxlevel ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R3VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i6)' ) '  Exceeding recursion maximum of ', maxlevel
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( n_segment > 0 ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine radec_distance_3d ( ra1, dec1, ra2, dec2, theta )
!
!*******************************************************************************
!
!! RADEC_DISTANCE_3D - angular distance, astronomical units, sphere in 3D.
!
!
!  Discussion:
!
!    Right ascension is measured in hours, between 0 and 24, and
!    essentially measures longitude.
!
!    Declination measures the angle from the equator towards the north pole,
!    and ranges from -90 (South Pole) to 90 (North Pole).
!
!    On the unit sphere, the angular separation between two points is 
!    equal to their geodesic or great circle distance.  On any other
!    sphere, multiply the angular separation by the radius of the
!    sphere to get the geodesic or great circle distance.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RA1, DEC1, RA2, DEC2, the right ascension and declination
!    of the two points.
!
!    Output, real THETA, the angular separation between the points,
!    in radians.
!
  implicit none
!
  real arc_cosine
  real cos_theta
  real dec1
  real dec2
  real degrees_to_radians
  real norm_v1
  real norm_v2
  real phi1
  real phi2
  real ra1
  real ra2
  real theta
  real theta1
  real theta2
  real v1(3)
  real v2(3)
!
  theta1 = degrees_to_radians ( 15.0E+00 * ra1 )
  phi1 = degrees_to_radians ( dec1 )

  v1(1:3) = (/ cos ( theta1 ) * cos ( phi1 ), &
               sin ( theta1 ) * cos ( phi1 ), &
                                sin ( phi1 ) /)

  norm_v1 = sqrt ( sum ( v1(1:3)**2 ) )

  theta2 = degrees_to_radians ( 15.0E+00 * ra2 )
  phi2 = degrees_to_radians ( dec2 )

  v2(1:3) = (/ cos ( theta2 ) * cos ( phi2 ), &
               sin ( theta2 ) * cos ( phi2 ), &
                                sin ( phi2 ) /)

  norm_v2 = sqrt ( sum ( v2(1:3)**2 ) )

  cos_theta = dot_product ( v1(1:3), v2(1:3) ) / ( norm_v1 * norm_v2 )

  theta = arc_cosine ( cos_theta )

  return
end
subroutine radec_to_xyz ( ra, dec, x, y, z )
!
!*******************************************************************************
!
!! RADEC_TO_XYZ converts right ascension/declination to (X,Y,Z) coordinates.
!
!
!  Discussion:
!
!    Right ascension is measured in hours, between 0 and 24, and
!    essentially measures longitude.
!
!    Declination measures the angle from the equator towards the north pole,
!    and ranges from -90 (South Pole) to 90 (North Pole).
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RA, DEC, the right ascension and declination of a point.
!
!    Output, real X, Y, Z, the corresponding coordinates of a point with
!    radius 1.
!
  implicit none
!
  real dec
  real degrees_to_radians
  real phi
  real ra
  real theta
  real x
  real y
  real z
!
  theta = degrees_to_radians ( 15.0E+00 * ra )
  phi = degrees_to_radians ( dec )

  x = cos ( theta ) * cos ( phi )
  y = sin ( theta ) * cos ( phi )
  z = sin ( phi )

  return
end
function radians_to_degrees ( angle )
!
!*******************************************************************************
!
!! RADIANS_TO_DEGREES converts an angle from radians to degrees.
!
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, an angle in radians.
!
!    Output, real RADIANS_TO_DEGREES, the equivalent angle in degrees.
!
  implicit none
!
  real angle
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real radians_to_degrees
!
  radians_to_degrees = ( angle / pi ) * 180.0E+00

  return
end
subroutine radians_to_dms ( radians, degrees, minutes, seconds )
!
!*******************************************************************************
!
!! RADIANS_TO_DMS converts an angle from radians to degrees/minutes/seconds.
!
!
!  Modified:
!
!    05 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RADIANS, the angle in radians.
!
!    Output, integer DEGREES, MINUTES, SECONDS, the equivalent angle in
!    degrees, minutes, and seconds.
!
  implicit none
!
  real angle
  integer degrees
  integer minutes
  real, parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00
  real radians
  integer seconds
!
  angle = pi * radians / 180.0E+00

  degrees = int ( angle )
  angle = ( angle - real ( degrees ) ) * 60.0E+00
  minutes = int ( angle )
  angle = ( angle - real ( degrees ) ) * 60.0E+00
  seconds = nint ( angle )

  return
end
subroutine random_initialize ( seed )
!
!*******************************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator.
!
  implicit none
!
  integer count
  integer count_max
  integer count_rate
  integer i
  integer seed
  integer, allocatable :: seed_vector(:)
  integer seed_size
  real t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) '  Initialize RANDOM_NUMBER with user SEED = ', seed

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) &
      '  Initialize RANDOM_NUMBER with arbitrary SEED = ', seed

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end
function rmat2_det ( a )
!
!*******************************************************************************
!
!! RMAT2_DET computes the determinant of a 2 by 2 matrix.
!
!
!  Discussion:
!
!    The determinant is the area spanned by the vectors making up the rows
!    or columns of the matrix.
!
!  Formula:
!
!    RMAT2_DET = A(1,1) * A(2,2) - A(1,2) * A(2,1).
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(2,2), the matrix whose determinant is desired.
!
!    Output, real RMAT2_DET, the determinant of the matrix.
!
  implicit none
!
  real a(2,2)
  real rmat2_det
!
  rmat2_det = a(1,1) * a(2,2) - a(1,2) * a(2,1)

  return
end
subroutine rmat2_inverse ( a, b, det )
!
!*******************************************************************************
!
!! RMAT2_INVERSE inverts a 2 by 2 real matrix using Cramer's rule.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(2,2), the matrix to be inverted.
!
!    Output, real B(2,2), the inverse of the matrix A.
!
!    Output, real DET, the determinant of the matrix A.
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
  implicit none
!
  real a(2,2)
  real b(2,2)
  real det
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0E+00 ) then

    b(1:2,1:2) = 0.0E+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = + a(2,2) / det
  b(1,2) = - a(1,2) / det
  b(2,1) = - a(2,1) / det
  b(2,2) = + a(1,1) / det

  return
end
function rmat3_det ( a )
!
!*******************************************************************************
!
!! RMAT3_DET computes the determinant of a 3 by 3 matrix.
!
!
!  The determinant is the volume of the shape spanned by the vectors
!  making up the rows or columns of the matrix.
!
!  Formula:
!
!    det = a11 * a22 * a33 - a11 * a23 * a32
!        + a12 * a23 * a31 - a12 * a21 * a33
!        + a13 * a21 * a32 - a13 * a22 * a31
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(3,3), the matrix whose determinant is desired.
!
!    Output, real RMAT3_DET, the determinant of the matrix.
!
  implicit none
!
  real a(3,3)
  real rmat3_det
!
  rmat3_det =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
              + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
              + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

  return
end
subroutine rmat3_inverse ( a, b, det )
!
!*******************************************************************************
!
!! RMAT3_INVERSE inverts a 3 by 3 real matrix using Cramer's rule.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(3,3), the matrix to be inverted.
!
!    Output, real B(3,3), the inverse of the matrix A.
!
!    Output, real DET, the determinant of the matrix A.
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
  implicit none
!
  real a(3,3)
  real b(3,3)
  real det
  integer i
  integer j
!
!  Compute the determinant of A
!
  det =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
        + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0E+00 ) then

    b(1:3,1:3) = 0.0E+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit
!  formula.
!
  b(1,1) = + ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
  b(1,2) = - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
  b(1,3) = + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

  b(2,1) = - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
  b(2,2) = + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
  b(2,3) = - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

  b(3,1) = + ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
  b(3,2) = - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
  b(3,3) = + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

  return
end
function rmat4_det ( a )
!
!*******************************************************************************
!
!! RMAT4_DET computes the determinant of a 4 by 4 matrix.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(4,4), the matrix whose determinant is desired.
!
!    Output, real RMAT4_DET, the determinant of the matrix.
!
  implicit none
!
  real a(4,4)
  real rmat4_det
!
  rmat4_det = &
      a(1,1) * ( &
        a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
      - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
      + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
    - a(1,2) * ( &
        a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
      - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
      + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
    + a(1,3) * ( &
        a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
      - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
      + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
    - a(1,4) * ( &
        a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
      - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
      + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
function rmat5_det ( a )
!
!*******************************************************************************
!
!! RMAT5_DET computes the determinant of a 5 by 5 matrix.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(5,5), the matrix whose determinant is desired.
!
!    Output, real RMAT5_DET, the determinant of the matrix.
!
  implicit none
!
  real a(5,5)
  real b(4,4)
  real rmat4_det
  real rmat5_det
  integer i
  integer inc
  integer j
  integer k
!
!  Expand the determinant into the sum of the determinants of the
!  five 4 by 4 matrices created by dropping row 1, and column k.
!
  rmat5_det = 0.0E+00

  do k = 1, 5

    do i = 1, 4
      do j = 1, 4

        if ( j < k ) then
          inc = 0
        else
          inc = 1
        end if

        b(i,j) = a(i+1,j+inc)

      end do
    end do

    rmat5_det = rmat5_det + (-1)**( k + 1 ) * a(1,k) * rmat4_det ( b )

  end do

  return
end
subroutine rmat_print ( lda, m, n, a, title )
!
!*******************************************************************************
!
!! RMAT_PRINT prints a real matrix.
!
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer LDA, the leading dimension of A.
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer lda
  integer n
!
  real a(lda,n)
  integer i
  integer j
  integer jhi
  integer jlo
  integer m
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 5
    jhi = min ( jlo + 4, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,5(i7,7x))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine rmat_solve ( a, n, nrhs, info )
!
!*******************************************************************************
!
!! RMAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real A(N,N+NRHS), contains in rows and columns 1
!    to N the coefficient matrix, and in columns N+1 through
!    N+NRHS, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Input, integer NRHS, the number of right hand sides.  NRHS
!    must be at least 0.
!
!    Output, integer INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none
!
  integer n
  integer nrhs
!
  real a(n,n+nrhs)
  real apivot
  real factor
  integer i
  integer info
  integer ipivot
  integer j
  integer k
  real temp
!
  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j+1, n
      if ( abs ( a(i,j) ) > abs ( apivot ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0E+00 ) then
      info = j
      return
    end if
!
!  Interchange.
!
    do i = 1, n + nrhs
      call r_swap ( a(ipivot,i), a(j,i) )
    end do
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0E+00
    a(j,j+1:n+nrhs) = a(j,j+1:n+nrhs) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then

        factor = a(i,j)
        a(i,j) = 0.0E+00
        a(i,j+1:n+nrhs) = a(i,j+1:n+nrhs) - factor * a(j,j+1:n+nrhs)

      end if

    end do

  end do

  return
end
subroutine rotation_axis2mat_3d ( axis, angle, a )
!
!*******************************************************************************
!
!! ROTATION_AXIS2MAT_3D converts a rotation from axis to matrix format in 3D.
!
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Modified:
!
!    27 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real AXIS(3), the axis vector which remains unchanged by
!    the rotation.
!
!    Input, real ANGLE, the angular measurement of the rotation about
!    the axis, in radians.
!
!    Output, real A(3,3), the rotation matrix.
!
  implicit none
!
  real a(3,3)
  real angle
  real axis(3)
  real ca
  real norm
  real sa
  real v1
  real v2
  real v3
!
  v1 = axis(1)
  v2 = axis(2)
  v3 = axis(3)

  norm = sqrt ( v1 * v1 + v2 * v2 + v3 * v3 )

  if ( norm == 0.0E+00 ) then
    return
  end if

  v1 = v1 / norm
  v2 = v2 / norm
  v3 = v3 / norm

  ca = cos ( angle )
  sa = sin ( angle )

  a(1,1) =                v1 * v1 + ca * ( 1.0E+00 - v1 * v1 )
  a(1,2) = ( 1.0E+00 - ca ) * v1 * v2 - sa * v3
  a(1,3) = ( 1.0E+00 - ca ) * v1 * v3 + sa * v2

  a(2,1) = ( 1.0E+00 - ca ) * v2 * v1 + sa * v3
  a(2,2) =                v2 * v2 + ca * ( 1.0E+00 - v2 * v2 )
  a(2,3) = ( 1.0E+00 - ca ) * v2 * v3 - sa * v1

  a(3,1) = ( 1.0E+00 - ca ) * v3 * v1 - sa * v2
  a(3,2) = ( 1.0E+00 - ca ) * v3 * v2 + sa * v1
  a(3,3) =                v3 * v3 + ca * ( 1.0E+00 - v3 * v3 )

  return
end
subroutine rotation_axis2quat_3d ( axis, angle, q )
!
!*******************************************************************************
!
!! ROTATION_AXIS2QUAT_3D converts a rotation from axis to quaternion format in 3D.
!
!
!  Definition:
!
!    A rotation quaternion Q has the form:
!
!      Q = A + Bi + Cj + Dk
!
!    where A, B, C and D are real numbers, and i, j, and k are to be regarded
!    as symbolic constant basis vectors, similar to the role of the "i"
!    in the representation of imaginary numbers.
!
!    A is the cosine of half of the angle of rotation.  (B,C,D) is a
!    unit vector pointing in the direction of the axis of rotation.
!    Rotation multiplication and inversion can be carried out using
!    this format and the usual rules for quaternion multiplication
!    and inversion.
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real AXIS(3), the axis vector which remains unchanged by
!    the rotation.
!
!    Input, real ANGLE, the angular measurement of the rotation about
!    the axis, in radians.
!
!    Output, real Q(4), the quaternion representing the rotation.
!
  implicit none
!
  real axis(3)
  real angle
  real norm
  real q(4)
!
  norm = sqrt ( axis(1) * axis(1) + axis(2) * axis(2) + axis(3) * axis(3) )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROTATION_AXIS2QUAT_3D - Fatal error!'
    write ( *, '(a)' ) '  The axis vector is null.'
  end if

  q(1) = cos ( 0.5E+00 * angle )

  q(2) = axis(1) * sin ( 0.5E+00 * angle ) / norm
  q(3) = axis(2) * sin ( 0.5E+00 * angle ) / norm
  q(4) = axis(3) * sin ( 0.5E+00 * angle ) / norm

  return
end
subroutine rotation_axis_vector_3d ( axis, angle, v, w )
!
!*******************************************************************************
!
!! ROTATION_AXIS_VECTOR_3D rotates a vector around an axis vector in 3D.
!
!
!  Modified:
!
!    31 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real AXIS(3), the axis vector for the rotation.
!
!    Input, real ANGLE, the angle, in radians, of the rotation.
!
!    Input, real V(3), the vector to be rotated.
!
!    Output, real W(3), the rotated vector.
!
  implicit none
!
  real angle
  real axis(3)
  real dot
  real norm
  real norm2
  real rot(3)
  real v(3)
  real w(3)
  real xa
  real xn
  real xn2
  real xp
  real ya
  real yn
  real yn2
  real yp
  real za
  real zn
  real zn2
  real zp
!
!  Compute the length of the rotation axis.
!
  xa = axis(1)
  ya = axis(2)
  za = axis(3)

  norm = sqrt ( xa * xa + ya * ya + za * za )

  if ( norm == 0.0E+00 ) then
    return
  end if

  xa = xa / norm
  ya = ya / norm
  za = za / norm
!
!  Compute the dot product of the vector and the rotation axis.
!
  dot = v(1) * xa + v(2) * ya + v(3) * za
!
!  Compute the parallel component of the vector.
!
  xp = dot * xa
  yp = dot * ya
  zp = dot * za
!
!  Compute the normal component of the vector.
!
  xn = v(1) - xp
  yn = v(2) - yp
  zn = v(3) - zp

  norm2 = sqrt ( xn * xn + yn * yn + zn * zn )

  if ( norm2 == 0.0E+00 ) then
    w(1:3) = (/ xp, yp, zp /)
    return
  end if

  xn = xn / norm2
  yn = yn / norm2
  zn = zn / norm2
!
!  Compute a second vector, lying in the plane, perpendicular
!  to V, and forming a right-handed system.
!
  call cross_3d ( xa, ya, za, xn, yn, zn, xn2, yn2, zn2 )

  norm = sqrt ( xn2 * xn2 + yn2 * yn2 + zn2 * zn2 )

  xn2 = xn2 / norm
  yn2 = yn2 / norm
  zn2 = zn2 / norm
!
!  Rotate the normal component by the angle.
!
  rot(1) = norm2 * ( cos ( angle ) * xn + sin ( angle ) * xn2 )
  rot(2) = norm2 * ( cos ( angle ) * yn + sin ( angle ) * yn2 )
  rot(3) = norm2 * ( cos ( angle ) * zn + sin ( angle ) * zn2 )
!
!  The rotated vector is the parallel component plus the rotated component.
!
  w(1) = xp + rot(1)
  w(2) = yp + rot(2)
  w(3) = zp + rot(3)

  return
end
subroutine rotation_mat2axis_3d ( a, axis, angle )
!
!*******************************************************************************
!
!! ROTATION_MAT2AXIS_3D converts a rotation from matrix to axis format in 3D.
!
!
!  Discussion:
!
!    The computation is based on the fact that a rotation matrix must
!    have an eigenvector corresponding to the eigenvalue of 1, hence:
!
!      ( A - I ) * v = 0.
!
!    The eigenvector V is the axis of rotation.
!
!  Reference:
!
!    Jack Kuipers
!    Quaternions and Rotation Sequences,
!    Princeton, 1998.
!
!  Modified:
!
!    27 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(3,3), the rotation matrix.
!
!    Output, real AXIS(3), the axis vector which remains unchanged by
!    the rotation.
!
!    Output, real ANGLE, the angular measurement of the rotation about
!    the axis, in radians.
!
  implicit none
!
  real a(3,3)
  real axis(3)
  real angle
  real arc_cosine
  real norm
!
  norm = sqrt ( ( a(3,2) - a(2,3) )**2 + ( a(1,3) - a(3,1) )**2 &
              + ( a(2,1) - a(1,2) )**2 )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROTATION_MAT2AXIS_3D - Fatal error!'
    write ( *, '(a)' ) '  A is not a rotation matrix,'
    write ( *, '(a)' ) '  or there are multiple axes of rotation.'
  end if

  axis(1) = ( a(3,2) - a(2,3) ) / norm
  axis(2) = ( a(1,3) - a(3,1) ) / norm
  axis(3) = ( a(2,1) - a(1,2) ) / norm
!
!  Find the angle.
!
  angle = arc_cosine ( 0.5E+00 * ( a(1,1) + a(2,2) + a(3,3) - 1.0E+00 ) )

  return
end
subroutine rotation_mat2quat_3d ( a, q )
!
!*******************************************************************************
!
!! ROTATION_MAT2QUAT_3D converts a rotation from matrix to quaternion format in 3D.
!
!
!  Discussion:
!
!    The computation is based on the fact that a rotation matrix must
!    have an eigenvector corresponding to the eigenvalue of 1, hence:
!
!      ( A - I ) * v = 0.
!
!    The eigenvector V is the axis of rotation.
!
!  Reference:
!
!    Jack Kuipers
!    Quaternions and Rotation Sequences,
!    Princeton, 1998.
!
!  Modified:
!
!    27 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(3,3), the rotation matrix.
!
!    Output, real Q(4), the quaternion representing the rotation.
!
  implicit none
!
  real a(3,3)
  real angle
  real arc_cosine
  real cos_phi
  real norm
  real q(4)
  real sin_phi
!
  norm = sqrt ( ( a(3,2) - a(2,3) )**2 + ( a(1,3) - a(3,1) )**2 &
              + ( a(2,1) - a(1,2) )**2 )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROTATION_MAT2AXIS_3D - Fatal error!'
    write ( *, '(a)' ) '  A is not a rotation matrix,'
    write ( *, '(a)' ) '  or there are multiple axes of rotation.'
  end if

  angle = arc_cosine ( 0.5E+00 * ( a(1,1) + a(2,2) + a(3,3) - 1.0E+00 ) )

  cos_phi = cos ( 0.5E+00 * angle )

  sin_phi = sqrt ( 1.0E+00 - cos_phi**2 )

  q(1) = cos_phi
  q(2) = sin_phi * ( a(3,2) - a(2,3) ) / norm
  q(3) = sin_phi * ( a(1,3) - a(3,1) ) / norm
  q(4) = sin_phi * ( a(2,1) - a(1,2) ) / norm

  return
end
subroutine rotation_mat_vector_3d ( a, v, w )
!
!*******************************************************************************
!
!! ROTATION_MAT_VECTOR applies a marix rotation to a vector in 3d.
!
!
!  Modified:
!
!    30 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(3,3), the matrix defining the rotation.
!
!    Input, real V(3), the vector to be rotated.
!
!    Output, real W(3), the rotated vector.
!
  implicit none
!
  real a(3,3)
  real v(3)
  real w(3)
!
  w(1:3) = matmul ( a(1:3,1:3), v(1:3) )

  return
end
subroutine rotation_quat2axis_3d ( q, axis, angle )
!
!*******************************************************************************
!
!! ROTATION_QUAT2AXIS_3D converts a rotation from quaternion to axis format in 3D.
!
!
!  Definition:
!
!    A rotation quaternion Q has the form:
!
!      Q = A + Bi + Cj + Dk
!
!    where A, B, C and D are real numbers, and i, j, and k are to be regarded
!    as symbolic constant basis vectors, similar to the role of the "i"
!    in the representation of imaginary numbers.
!
!    A is the cosine of half of the angle of rotation.  (B,C,D) is a
!    vector pointing in the direction of the axis of rotation.
!    Rotation multiplication and inversion can be carried out using
!    this format and the usual rules for quaternion multiplication
!    and inversion.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Q(4), the quaternion representing the rotation.
!
!    Output, real AXIS(3), the axis vector which remains unchanged by
!    the rotation.
!
!    Output, real ANGLE, the angular measurement of the rotation about
!    the axis, in radians.
!
  implicit none
!
  real axis(3)
  real angle
  real cos_phi
  real q(4)
  real sin_phi
!
  sin_phi = sqrt ( sum ( q(2:4)**2 ) )

  cos_phi = q(1)

  angle = 2.0E+00 * atan2 ( sin_phi, cos_phi )

  if ( sin_phi == 0.0E+00 ) then
    axis(1:3) = (/ 1.0E+00, 0.0E+00, 0.0E+00 /)
  else
    axis(1:3) = (/ q(2), q(3), q(4) /) / sin_phi
  end if

  return
end
subroutine rotation_quat2mat_3d ( q, a )
!
!*******************************************************************************
!
!! ROTATION_QUAT2MAT_3D converts a rotation from quaternion to matrix format in 3D.
!
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Modified:
!
!    27 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real AXIS(3), the axis vector which remains unchanged by
!    the rotation.
!
!    Input, real ANGLE, the angular measurement of the rotation about
!    the axis, in radians.
!
!    Output, real A(3,3), the rotation matrix.
!
  implicit none
!
  real a(3,3)
  real angle
  real ca
  real cos_phi
  real q(4)
  real sa
  real sin_phi
  real v1
  real v2
  real v3
!
  sin_phi = sqrt ( sum ( q(2:4)**2 ) )

  cos_phi = q(1)

  angle = 2.0E+00 * atan2 ( sin_phi, cos_phi )

  if ( sin_phi == 0.0E+00 ) then
    v1 = 1.0E+00
    v2 = 0.0E+00
    v3 = 0.0E+00
  else
    v1 = q(2) / sin_phi
    v2 = q(3) / sin_phi
    v3 = q(4) / sin_phi
  end if

  ca = cos ( angle )
  sa = sin ( angle )

  a(1,1) =                    v1 * v1 + ca * ( 1.0E+00 - v1 * v1 )
  a(1,2) = ( 1.0E+00 - ca ) * v1 * v2 - sa * v3
  a(1,3) = ( 1.0E+00 - ca ) * v1 * v3 + sa * v2

  a(2,1) = ( 1.0E+00 - ca ) * v2 * v1 + sa * v3
  a(2,2) =                    v2 * v2 + ca * ( 1.0E+00 - v2 * v2 )
  a(2,3) = ( 1.0E+00 - ca ) * v2 * v3 - sa * v1

  a(3,1) = ( 1.0E+00 - ca ) * v3 * v1 - sa * v2
  a(3,2) = ( 1.0E+00 - ca ) * v3 * v2 + sa * v1
  a(3,3) =                    v3 * v3 + ca * ( 1.0E+00 - v3 * v3 )

  return
end
subroutine rotation_quat_vector_3d ( q, v, w )
!
!*******************************************************************************
!
!! ROTATION_QUAT_VECTOR applies a quaternion rotation to a vector in 3d.
!
!
!  Discussion:
!
!    If Q is a unit quaternion that encodes a rotation of ANGLE
!    radians about the vector AXIS, then for an arbitrary real
!    vector V, the result W of the rotation on V can be written as:
!
!      W = Q * V * Conj(Q)
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Q(4), the quaternion defining the rotation.
!
!    Input, real V(3), the vector to be rotated.
!
!    Output, real V(3), the rotated vector.
!
  implicit none
!
  real q(4)
  real v(3)
  real w(3)
!
  w(1) = &
         ( 2.0E+00 * ( q(1) * q(1) + q(2) * q(2) ) - 1.0E+00 ) * v(1) &
       +   2.0E+00 * ( q(2) * q(3) - q(1) * q(4) )             * v(2) &
       +   2.0E+00 * ( q(2) * q(4) + q(1) * q(3) )             * v(3)

  w(2) = &
           2.0E+00 * ( q(2) * q(3) + q(1) * q(4) )             * v(1) &
       + ( 2.0E+00 * ( q(1) * q(1) + q(3) * q(3) ) - 1.0E+00 ) * v(2) &
       +   2.0E+00 * ( q(3) * q(4) - q(1) * q(2) )             * v(3)

  w(3) = &
           2.0E+00 * ( q(2) * q(4) - q(1) * q(3) )             * v(1) &
       +   2.0E+00 * ( q(3) * q(4) + q(1) * q(2) )             * v(2) &
       + ( 2.0E+00 * ( q(1) * q(1) + q(4) * q(4) ) - 1.0E+00 ) * v(3)

  return
end
subroutine rtris2 ( point_num, point_xy, tri_num, tri_vert, tri_nabe )
!
!*******************************************************************************
!
!! RTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of vertices.
!
!    Input/output, real POINT_XY(2,POINT_NUM), the coordinates of the vertices.
!    On output, the vertices have been sorted into dictionary order.
!
!    Output, integer TRI_NUM, the number of triangles in the triangulation;
!    TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the number 
!    of boundary vertices.
!
!    Output, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
!    The elements are indices of POINT_XY.  The vertices of the triangles are
!    in counter clockwise order.
!
!    Output, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list.
!    Positive elements are indices of TIL; negative elements are used for links
!    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TRI_NABE(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3).
!
  implicit none
!
  integer point_num
!
  real cmax
  integer e
  integer i
  integer ierr
  integer j
  integer k
  integer l
  integer ledg
  integer lr
  integer lrline
  integer ltri
  integer m
  integer m1
  integer m2
  integer n
  real point_xy(2,point_num)
  integer redg
  integer rtri
  integer stack(point_num)
  integer t
  real temp
  real tol
  integer top
  integer tri_nabe(3,point_num*2)
  integer tri_num
  integer tri_vert(3,point_num*2)
!
  tol = 100.0E+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r2vec_sort_quick_a ( point_num, point_xy )
!
!  Make sure that the data points are "reasonably" distinct.
!
  m1 = 1

  do i = 2, point_num

    m = m1
    m1 = i

    k = 0

    do j = 1, 2

      cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

      if ( abs ( point_xy(j,m) - point_xy(j,m1) ) > tol * ( cmax + 1.0E+00 ) ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RTRIS2 - Fatal error!'
      write ( *, '(a,i6)' ) '  Fails for point number I = ', i
      write ( *, '(a,i6)' ) '  M = ', m
      write ( *, '(a,i6)' ) '  M1 = ', m1
      write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
      write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
      ierr = 224
      return
    end if

  end do
!
!  Starting from points M1 and M2, search for a third point M that
!  makes a "healthy" triangle (M1,M2,M)
!
  m1 = 1
  m2 = 2
  j = 3

  do

    if ( j > point_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RTRIS2 - Fatal error!'
      ierr = 225
      return
    end if

    m = j

    lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0E+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the triangle information for (M1,M2,M), and for any other
!  triangles you created because points were collinear with M1, M2.
!
  tri_num = j - 2

  if ( lr == -1 ) then

    tri_vert(1,1) = m1
    tri_vert(2,1) = m2
    tri_vert(3,1) = m
    tri_nabe(3,1) = -3

    do i = 2, tri_num

      m1 = m2
      m2 = i+1
      tri_vert(1,i) = m1
      tri_vert(2,i) = m2
      tri_vert(3,i) = m
      tri_nabe(1,i-1) = -3 * i
      tri_nabe(2,i-1) = i
      tri_nabe(3,i) = i - 1

    end do

    tri_nabe(1,tri_num) = -3 * tri_num - 1
    tri_nabe(2,tri_num) = -5
    ledg = 2
    ltri = tri_num

  else

    tri_vert(1,1) = m2
    tri_vert(2,1) = m1
    tri_vert(3,1) = m
    tri_nabe(1,1) = -4

    do i = 2, tri_num
      m1 = m2
      m2 = i+1
      tri_vert(1,i) = m2
      tri_vert(2,i) = m1
      tri_vert(3,i) = m
      tri_nabe(3,i-1) = i
      tri_nabe(1,i) = -3 * i - 3
      tri_nabe(2,i) = i - 1
    end do

    tri_nabe(3,tri_num) = -3 * tri_num
    tri_nabe(2,1) = -3 * tri_num - 2
    ledg = 2
    ltri = 1

  end if
!
!  Insert the vertices one at a time from outside the convex hull, 
!  determine visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, point_num

    m = i
    m1 = tri_vert(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = tri_vert(ledg+1,ltri)
    else
      m2 = tri_vert(1,ltri)
    end if

    lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0E+00 )

    if ( lr > 0 ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tri_nabe(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg ( point_xy(1,m), point_xy(2,m), point_num, point_xy, tri_num, &
      tri_vert, tri_nabe, ltri, ledg, rtri, redg )

    n = tri_num + 1
    l = -tri_nabe(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -tri_nabe(e,t)
      m2 = tri_vert(e,t)

      if ( e <= 2 ) then
        m1 = tri_vert(e+1,t)
      else
        m1 = tri_vert(1,t)
      end if

      tri_num = tri_num + 1
      tri_nabe(e,t) = tri_num
      tri_vert(1,tri_num) = m1
      tri_vert(2,tri_num) = m2
      tri_vert(3,tri_num) = m
      tri_nabe(1,tri_num) = t
      tri_nabe(2,tri_num) = tri_num - 1
      tri_nabe(3,tri_num) = tri_num + 1
      top = top + 1

      if ( top > point_num ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RTRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
        return
      end if

      stack(top) = tri_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    tri_nabe(ledg,ltri) = -3 * n - 1
    tri_nabe(2,n) = -3 * tri_num - 2
    tri_nabe(3,tri_num) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, point_num, point_xy, tri_num, &
      tri_vert, tri_nabe, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RTRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
      return
    end if

  end do

  return
end
subroutine rvec_bracket ( n, x, xval, left, right )
!
!*******************************************************************************
!
!! RVEC_BRACKET searches a sorted array for successive brackets of a value.
!
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, length of input array.
!
!    Input, real X(N), an array that has been sorted into ascending order.
!
!    Input, real XVAL, a value to be bracketed.
!
!    Output, integer LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      XVAL > X(N), when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none
!
  integer n
!
  integer i
  integer left
  integer right
  real x(n)
  real xval
!
  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
subroutine rvec2_print ( n, a1, a2, title )
!
!*******************************************************************************
!
!! RVEC2_PRINT prints a pair of real vectors.
!
!
!  Modified:
!
!    14 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  real a1(n)
  real a2(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,2g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine rvec3_print ( n, a1, a2, a3, title )
!
!*******************************************************************************
!
!! RVEC3_PRINT prints a trio of real vectors.
!
!
!  Modified:
!
!    16 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A1(N), A2(N), A3(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,3g14.6)' ) i, a1(i), a2(i), a3(i)
  end do

  return
end
function rvec_eq ( n, a1, a2 )
!
!*******************************************************************************
!
!! RVEC_EQ is true if every pair of entries in two vectors is equal.
!
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real A1(N), A2(N), two vectors to compare.
!
!    Output, logical RVEC_EQ.
!    RVEC_EQ is .TRUE. if every pair of elements A1(I) and A2(I) are equal,
!    and .FALSE. otherwise.
!
  implicit none
!
  integer n
!
  real a1(n)
  real a2(n)
  logical rvec_eq
!
  rvec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
function rvec_gt ( n, a1, a2 )
!
!*******************************************************************************
!
!! RVEC_GT == ( A1 > A2 ) for real vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 > A2  <=>                              A1(1) > A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real A1(N), A2(N), the vectors to be compared.
!
!    Output, logical RVEC_GT, is TRUE if and only if A1 > A2.
!
  implicit none
!
  integer n
!
  real a1(n)
  real a2(n)
  integer i
  logical rvec_gt
!
  rvec_gt = .false.

  do i = 1, n

    if ( a1(i) > a2(i) ) then
      rvec_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      rvec_gt = .false.
      exit
    end if

  end do

  return
end
function rvec_lt ( n, a1, a2 )
!
!*******************************************************************************
!
!! RVEC_LT == ( A1 < A2 ) for real vectors.
!
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 < A2  <=>                              A1(1) < A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real A1(N), A2(N), the vectors to be compared.
!
!    Output, logical RVEC_LT, is TRUE if and only if A1 < A2.
!
  implicit none
!
  integer n
!
  real a1(n)
  real a2(n)
  integer i
  logical rvec_lt
!
  rvec_lt = .false.

  do i = 1, n

    if ( a1(i) < a2(i) ) then
      rvec_lt = .true.
      exit
    else if ( a1(i) > a2(i) ) then
      rvec_lt = .false.
      exit
    end if

  end do

  return
end
subroutine rvec_print_2d ( x, y, title )
!
!*******************************************************************************
!
!! RVEC_PRINT_2D prints a 2D vector.
!
!
!  Discussion:
!
!    A format is used which suggests a coordinate pair:
!
!  Example:
!
!    Center : ( 1.23, 7.45 )
!
!  Modified:
!
!    03 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, a title to be printed, or ' '.
!
  implicit none
!
  character ( len = * ) title
  real x
  real y
!
  if ( len_trim ( title ) == 0 ) then
    write ( *, '( a1, g14.6, a1, g14.6, a1 )' ) '(', x, ',', y, ')'
  else
    write ( *, '( a, a4, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ), ' : (', x, ',', y, ')'
  end if

  return
end
subroutine rvec_print_3d ( x, y, z, title )
!
!*******************************************************************************
!
!! RVEC_PRINT_3D prints a 3D vector.
!
!
!  Discussion:
!
!    A format is used which suggests a coordinate triple:
!
!  Example:
!
!    Center : ( 1.23, 7.45, -1.45 )
!
!  Modified:
!
!    03 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, Z, the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, a title to be printed, or ' '.
!
  implicit none
!
  character ( len = * ) title
  real x
  real y
  real z
!
  if ( len_trim ( title ) == 0 ) then
    write ( *, '( a1, g14.6, a1, g14.6, a1, g14.6, a1 )' ) &
      '(', x, ',', y, ',', z, ')'
  else
    write ( *, '( a, a4, g14.6, a1, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ) , ' : (', x, ',', y, ',', z, ')'
  end if

  return
end
subroutine rvec_swap ( n, a1, a2 )
!
!*******************************************************************************
!
!! RVEC_SWAP swaps the entries of two real vectors.
!
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the arrays.
!
!    Input/output, real A1(N), A2(N), the vectors to swap.
!
  implicit none
!
  integer n
!
  real a1(n)
  real a2(n)
  real a3(n)
!
  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine shape_point_dist_2d ( xc, yc, x1, y1, nside, x, y, dist )
!
!*******************************************************************************
!
!! SHAPE_POINT_DIST_2D: distance ( regular shape, point ) in 2D.
!
!
!  Discussion:
!
!    The "regular shape" is assumed to be an equilateral and equiangular
!    polygon, such as the standard square, pentagon, hexagon, and so on.
!
!  Modified:
!
!    14 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XC, YC, the center of the shape.
!
!    Input, real X1, Y1, the first vertex of the shape.
!
!    Input, integer NSIDE, the number of sides in the shape.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, real DIST, the distance from the point to the shape.
!
  implicit none
!
  real angle
  real angle_deg_2d
  real angle2
  real degrees_to_radians
  real dist
  real dist2
  real enorm0_2d
  integer nside
  real radius
  real sector_angle
  integer sector_index
  real x
  real x1
  real xa
  real xb
  real xc
  real y
  real y1
  real ya
  real yb
  real yc
!
!  Determine the angle subtended by a single side.
!
  sector_angle = 360.0E+00 / real ( nside )
!
!  How long is the half-diagonal?
!
  radius = enorm0_2d ( x1, y1, xc, yc )
!
!  If the radius is zero, then the shape is a point and the computation is easy.
!
  if ( radius == 0.0E+00 ) then
    dist = enorm0_2d ( x, y, xc, yc )
    return
  end if
!
!  If the test point is at the center, the computation is easy.
!
  dist2 = enorm0_2d ( x, y, xc, yc )

  if ( dist2 == 0.0E+00 ) then
    dist = radius
    return
  end if
!
!  Determine the angle between the ray to the first corner,
!  and the ray to the test point.
!
  angle = angle_deg_2d ( x1, y1, xc, yc, x, y )
!
!  Determine the sector of the point.
!
  sector_index = int ( angle / sector_angle ) + 1
!
!  Generate the two corner points that terminate the SECTOR-th side.
!
  angle2 = real ( sector_index - 1 ) * sector_angle
  angle2 = degrees_to_radians ( angle2 )

  call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xa, ya )

  angle2 = real ( sector_index ) * sector_angle
  angle2 = degrees_to_radians ( angle2 )

  call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xb, yb )
!
!  Determine the distance from the test point to the line segment that
!  is the SECTOR-th side.
!
  call line_seg_point_dist_2d ( xa, ya, xb, yb, x, y, dist )

  return
end
subroutine shape_point_near_2d ( xc, yc, x1, y1, nside, x, y, xn, yn, dist )
!
!*******************************************************************************
!
!! SHAPE_POINT_NEAR_2D: nearest point ( regular shape, point ) in 2D.
!
!
!  Discussion:
!
!    The "regular shape" is assumed to be an equilateral and equiangular
!    polygon, such as the standard square, pentagon, hexagon, and so on.
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XC, YC, the center of the shape.
!
!    Input, real X1, Y1, the first vertex of the shape.
!
!    Input, integer NSIDE, the number of sides in the shape.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, real XN, YN, the point on the shape that is nearest
!    to the given point.
!
!    Output, real DIST, the distance between the points.
!
  implicit none
!
  real angle
  real angle_deg_2d
  real angle2
  real degrees_to_radians
  real dist
  real dist2
  real enorm0_2d
  integer nside
  real radius
  real sector_angle
  integer sector_index
  real t
  real x
  real x1
  real xa
  real xb
  real xc
  real xn
  real y
  real y1
  real ya
  real yb
  real yc
  real yn
!
!  Determine the angle subtended by a single side.
!
  sector_angle = 360.0E+00 / real ( nside )
!
!  How long is the half-diagonal?
!
  radius = enorm0_2d ( x1, y1, xc, yc )
!
!  If the radius is zero, then the shape is a point and the computation is easy.
!
  if ( radius == 0.0E+00 ) then
    dist = enorm0_2d ( x, y, xc, yc )
    return
  end if
!
!  If the test point is at the center, the computation is easy.
!
  dist2 = enorm0_2d ( x, y, xc, yc )

  if ( dist2 == 0.0E+00 ) then
    dist = radius
    return
  end if
!
!  Determine the angle between the ray to the first corner,
!  and the ray to the test point.
!
  angle = angle_deg_2d ( x1, y1, xc, yc, x, y )
!
!  Determine the sector of the point.
!
  sector_index = int ( angle / sector_angle ) + 1
!
!  Generate the two corner points that terminate the SECTOR-th side.
!
  angle2 = real ( sector_index - 1 ) * sector_angle
  angle2 = degrees_to_radians ( angle2 )

  call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xa, ya )

  angle2 = real ( sector_index ) * sector_angle
  angle2 = degrees_to_radians ( angle2 )

  call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xb, yb )
!
!  Determine the point on the SECTOR-th side of the shape which is
!  nearest.
!
  call line_seg_point_near_2d ( xa, ya, xb, yb, x, y, xn, yn, dist, t )

  return
end
subroutine shape_print_3d ( max_num, max_order, point_num, &
   face_num, face_order, point_coord, face_point )
!
!*******************************************************************************
!
!! SHAPE_PRINT_3D prints information about a polyhedron in 3D.
!
!
!  Modified:
!
!    11 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAX_NUM, the maximum number of faces, and of
!    points.  This quantity is used to dimension arrays.
!
!    Input, integer MAX_ORDER, the maximum order of any face.
!
!    Input, integer POINT_NUM, the number of points in the shape.
!
!    Input, integer FACE_NUM, the number of faces in the shape.
!
!    Input, integer FACE_ORDER, the number of vertices per face.
!
!    Input, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
!    the X, Y and Z coordinates of point J.
!
!    Input, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
!    contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter-clockwise direction defined
!    by the outward normal at the face.
!
  implicit none
!
  integer max_num
  integer max_order
!
  integer face_num
  integer face_order
  integer face_point(max_order,max_num)
  integer i
  integer j
  integer point_num
  real point_coord(3,max_num)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_PRINT_3D'
  write ( *, '(a)' ) '  Information about a regular polyhedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of vertices is ', point_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertex coordinates are:'
  write ( *, '(a)' ) ' '
  do i = 1, point_num
    write ( *, '(i4,2x,3f10.4)' ) i, point_coord(1:3,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of faces is ', face_num
  write ( *, '(a,i6)' ) '  The order of each face is ', face_order
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vertices in each face are:'
  write ( *, '(a)' ) ' '

  do i = 1, face_num
    write ( *, '(i5,3x,10i5)' ) i, face_point(1:face_order,i)
  end do

  return
end
subroutine shape_ray_int_2d ( xc, yc, x1, y1, nside, xa, ya, xb, yb, xi, yi )
!
!*******************************************************************************
!
!! SHAPE_RAY_INT_2D: intersection ( regular shape, ray ) in 2D.
!
!
!  Discussion:
!
!    The "regular shape" is assumed to be an equilateral and equiangular
!    polygon, such as the standard square, pentagon, hexagon, and so on.
!
!    The origin of the ray is assumed to be inside the shape.  This
!    guarantees that the ray will intersect the shape in exactly one point.
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real XC, YC, the center of the shape.
!
!    Input, real X1, Y1, the first vertex of the shape.
!
!    Input, integer NSIDE, the number of sides in the shape.
!
!    Input, real XA, YA, the origin of the ray.
!
!    Input, real XB, YB, a second point on the ray.
!
!    Output, real XI, YI, the point on the shape intersected by the ray.
!
  implicit none
!
  real angle2
  real degrees_to_radians
  real enorm0_2d
  logical inside
  integer ival
  integer nside
  real radius
  real sector_angle
  integer sector_index
  real x1
  real xa
  real xb
  real xc
  real xi
  real xv1
  real xv2
  real y1
  real ya
  real yb
  real yc
  real yi
  real yv1
  real yv2
!
!  Warning!
!  No check is made to ensure that the ray origin is inside the shape.
!  These calculations are not valid if that is not true!
!
!  Determine the angle subtended by a single side.
!
  sector_angle = 360.0E+00 / real ( nside )
!
!  How long is the half-diagonal?
!
  radius = enorm0_2d ( x1, y1, xc, yc )
!
!  If the radius is zero, refuse to continue.
!
  if ( radius == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE_RAY_INT_2D - Fatal error!'
    write ( *, '(a)' ) '  The shape has radius zero.'
    stop
  end if
!
!  Determine which sector side intersects the ray.
!
  xv2 = 0.0E+00
  yv2 = 0.0E+00

  do sector_index = 1, nside
!
!  Determine the two vertices that define this sector.
!
    if ( sector_index == 1 ) then

      angle2 = real ( sector_index - 1 ) * sector_angle
      angle2 = degrees_to_radians ( angle2 )

      call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xv1, yv1 )

    else

      xv1 = xv2
      yv1 = yv2

    end if

    angle2 = real ( sector_index ) * sector_angle
    angle2 = degrees_to_radians ( angle2 )

    call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xv2, yv2 )
!
!  Draw the angle from one vertex to the ray origin to the next vertex,
!  and see if that angle contains the ray.  If so, then the ray
!  must intersect the shape side of that sector.
!
    call angle_contains_ray_2d ( xv1, yv1, xa, ya, xv2, yv2, xb, yb, inside )

    if ( inside ) then
!
!  Determine the intersection of the lines defined by the ray and the
!  sector side.  (We're already convinced that the ray and sector line
!  segment intersect, so we can use the simpler code that treats them
!  as full lines).
!
      call lines_exp_int_2d ( xa, ya, xb, yb, xv1, yv1, xv2, yv2, ival, xi, yi )

      return

    end if

  end do
!
!  If the calculation fell through the loop, then something's wrong.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_RAY_INT_2D - Fatal error!'
  write ( *, '(a)' ) '  Cannot find intersection of ray and shape.'
  stop
end
subroutine soccer_shape_3d ( point_num, face_num, face_order, point_coord, &
  face_point )
!
!*******************************************************************************
!
!! SOCCER_SHAPE_3D describes a truncated icosahedron in 3D.
!
!
!  Discussion:
!
!    The shape is a truncated icosahedron, which is the design used
!    on a soccer ball.  There are 12 pentagons and 20 hexagons.
!
!    I think I understand the RESHAPE command, but at least with
!    2-dimensional arrays, I am getting bizarre results if I
!    have a situation like this:
!
!      integer a(10,20)
!
!      a(1:8,1:12) = reshape ( (/ ... /), (/ 8,12 /) )
!
!    What happens is that the array is returned with every column
!    being the same as the first column.  I gave up trying to choose
!    between blaming me, FORTRAN 90, or the SGI, and just dumbed
!    down the code.
!
!  Reference:
!
!    http://polyhedra.wolfram.com/uniform/u25.html
!
!  Modified:
!
!    06 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer POINT_NUM, the number of points in the shape (60).
!
!    Output, integer FACE_NUM, the number of faces in the shape (32).
!
!    Output, real POINT_COORD(3,60); POINT_COORD(*,J) contains
!    the X, Y and Z coordinates of point J.
!
!    Output, integer FACE_ORDER(32), the number of vertices per face.
!
!    Output, integer FACE_POINT(6,32); FACE_POINT(I,J)
!    contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter-clockwise direction defined
!    by the outward normal at the face.
!
  implicit none
!
  integer face_num
  integer face_order(32)
  integer face_point(6,32)
  integer point_num
  real point_coord(3,60)
!
  point_num = 60
  face_num = 32
  face_order = 6
!
!  Set point coordinates.
!
  point_coord(1:3,1:60) = reshape ( &
    (/ -1.00714,    0.153552,   0.067258, &
       -0.960284,   0.0848813, -0.33629,  &
       -0.95172,   -0.153552,   0.33629,  &
       -0.860021,   0.529326,   0.150394, &
       -0.858,     -0.290893,  -0.470806, &
       -0.849436,  -0.529326,   0.201774, &
       -0.802576,  -0.597996,  -0.201774, &
       -0.7842,     0.418215,  -0.502561, &
       -0.749174,  -0.0848813,  0.688458, &
       -0.722234,   0.692896,  -0.201774, &
       -0.657475,   0.597996,   0.502561, &
       -0.602051,   0.290893,   0.771593, &
       -0.583675,  -0.692896,   0.470806, &
       -0.579632,  -0.333333,  -0.771593, &
       -0.52171,   -0.418215,   0.771593, &
       -0.505832,   0.375774,  -0.803348, &
       -0.489955,  -0.830237,  -0.33629,  &
       -0.403548,   0.,        -0.937864, &
       -0.381901,   0.925138,  -0.201774, &
       -0.352168,  -0.666667,  -0.688458, &
       -0.317142,   0.830237,   0.502561, &
       -0.271054,  -0.925138,   0.33629,  &
       -0.227464,   0.333333,   0.937864, &
       -0.224193,  -0.993808,  -0.067258, &
       -0.179355,   0.993808,   0.150394, &
       -0.165499,   0.608015,  -0.803348, &
       -0.147123,  -0.375774,   0.937864, &
       -0.103533,   0.882697,  -0.502561, &
       -0.0513806,  0.666667,   0.771593, &
        0.0000000,  0.,         1.021,    &
        0.0000000,  0.,        -1.021,    &
        0.0513806, -0.666667,  -0.771593, &
        0.103533,  -0.882697,   0.502561, &
        0.147123,   0.375774,  -0.937864, &
        0.165499,  -0.608015,   0.803348, &
        0.179355,  -0.993808,  -0.150394, &
        0.224193,   0.993808,   0.067258, &
        0.227464,  -0.333333,  -0.937864, &
        0.271054,   0.925138,  -0.33629,  &
        0.317142,  -0.830237,  -0.502561, &
        0.352168,   0.666667,   0.688458, &
        0.381901,  -0.925138,   0.201774, &
        0.403548,   0.,         0.937864, &
        0.489955,   0.830237,   0.33629,  &
        0.505832,  -0.375774,   0.803348, &
        0.521710,   0.418215,  -0.771593, &
        0.579632,   0.333333,   0.771593, &
        0.583675,   0.692896,  -0.470806, &
        0.602051,  -0.290893,  -0.771593, &
        0.657475,  -0.597996,  -0.502561, &
        0.722234,  -0.692896,   0.201774, &
        0.749174,   0.0848813, -0.688458, &
        0.784200,  -0.418215,   0.502561, &
        0.802576,   0.597996,   0.201774, &
        0.849436,   0.529326,  -0.201774, &
        0.858000,   0.290893,   0.470806, &
        0.860021,  -0.529326,  -0.150394, &
        0.951720,   0.153552,  -0.33629,  &
        0.960284,  -0.0848813,  0.33629,  &
        1.007140,  -0.153552,  -0.067258 /), &
    (/ 3, 60 /) )
!
!  Set face orders.
!
  face_order(1:32) = (/ &
    6, 6, 5, 6, 5, 6, 5, 6, 6, 6, &
    5, 6, 5, 6, 5, 6, 6, 6, 5, 6, &
    5, 5, 6, 6, 6, 5, 6, 5, 6, 6, &
    5, 6 /)
!
!  Set faces.
!
  face_point(1:6,1:32) = reshape ( &
    (/ 30, 43, 47, 41, 29, 23, &
       30, 23, 12,  9, 15, 27, &
       30, 27, 35, 45, 43,  0, &
       43, 45, 53, 59, 56, 47, &
       23, 29, 21, 11, 12,  0, &
       27, 15, 13, 22, 33, 35, &
       47, 56, 54, 44, 41,  0, &
       45, 35, 33, 42, 51, 53, &
       12, 11,  4,  1,  3,  9, &
       29, 41, 44, 37, 25, 21, &
       15,  9,  3,  6, 13,  0, &
       56, 59, 60, 58, 55, 54, &
       53, 51, 57, 60, 59,  0, &
       11, 21, 25, 19, 10,  4, &
       33, 22, 24, 36, 42,  0, &
       13,  6,  7, 17, 24, 22, &
       54, 55, 48, 39, 37, 44, &
       51, 42, 36, 40, 50, 57, &
        4, 10,  8,  2,  1,  0, &
        3,  1,  2,  5,  7,  6, &
       25, 37, 39, 28, 19,  0, &
       55, 58, 52, 46, 48,  0, &
       60, 57, 50, 49, 52, 58, &
       10, 19, 28, 26, 16,  8, &
       36, 24, 17, 20, 32, 40, &
        7,  5, 14, 20, 17,  0, &
       48, 46, 34, 26, 28, 39, &
       50, 40, 32, 38, 49,  0, &
        8, 16, 18, 14,  5,  2, &
       46, 52, 49, 38, 31, 34, &
       16, 26, 34, 31, 18,  0, &
       32, 20, 14, 18, 31, 38 /), &
    (/ 6, 32 /) )

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )
!
!*******************************************************************************
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Modified:
!
!    19 May 1999
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if I > J;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    ISGN => 0 means I is greater than or equal to J.
!
  implicit none
!
  integer i
  integer indx
  integer isgn
  integer j
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( isgn > 0 ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  INDX > 0, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

  return
end
subroutine sphere_dia2imp_3d ( x1, y1, z1, x2, y2, z2, r, xc, yc, zc )
!
!*******************************************************************************
!
!! SPHERE_DIA2IMP_3D converts a diameter to an implicit sphere in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, are the (X,Y,Z) coordinates
!    of two points which form a diameter of the sphere.
!
!    Output, real R, the computed radius of the sphere.
!
!    Output, real XC, YC, ZC, the computed center of the sphere.
!
  implicit none
!
  real enorm0_3d
  real r
  real x1
  real x2
  real xc
  real y1
  real y2
  real yc
  real z1
  real z2
  real zc
!
  r = 0.5E+00 * enorm0_3d ( x1, y1, z1, x2, y2, z2 )

  xc = 0.5E+00 * ( x1 + x2 )
  yc = 0.5E+00 * ( y1 + y2 )
  zc = 0.5E+00 * ( z1 + z2 )

  return
end
subroutine sphere_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, &
  z3, x4, y4, z4, r, xc, yc, zc )
!
!*******************************************************************************
!
!! SPHERE_EXP2IMP_3D converts a sphere from explicit to implicit form in 3D.
!
!
!  Formula:
!
!    An explicit sphere in 3D is determined by four points,
!    which should be distinct, and not coplanar.
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4,
!    the coordinates of four distinct noncoplanar points on the sphere.
!
!    Output, real R, XC, YC, ZC, the radius and coordinates of the center
!    of the sphere.  If the linear system is
!    singular, then R = -1, XC = YC = ZC = 0.
!
  implicit none
!
  real r
  real x1
  real x2
  real x3
  real x4
  real xc
  real y1
  real y2
  real y3
  real y4
  real yc
  real z1
  real z2
  real z3
  real z4
  real zc
!
  call tetra_circumsphere_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
    x4, y4, z4, r, xc, yc, zc )

  return
end
subroutine sphere_exp_contains_point_3d ( x1, y1, z1, x2, y2, z2, &
  x3, y3, z3, x4, y4, z4, x, y, z, inside )
!
!*******************************************************************************
!
!! SPHERE_EXP_CONTAINS_POINT_3D determines if an explicit sphere contains a point in 3D.
!
!
!  Formula:
!
!    An explicit sphere in 3D is determined by four points,
!    which should be distinct, and not coplanar.
!
!  Method:
!
!    The computation checks the determinant of:
!
!      x1  y1  z1  x1**2+y1**2+z1**2  1
!      x2  y2  z2  x2**2+y2**2+z2**2  1
!      x3  y3  z3  x3**2+y3**2+z3**2  1
!      x4  y4  z4  x4**2+y4**2+z4**2  1
!      x   y   z   x**2 +y**2 +z**2   1
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
!    (X,Y,Z) coordinates of four points that lie on a circle.
!
!    Input, real X, Y, Z, the (X,Y,Z) coordinates of a point, whose
!    position relative to the sphere is desired.
!
!    Output, logical INSIDE, is TRUE if the point is in the sphere,
!    FALSE otherwise.
!
  implicit none
!
  real a(5,5)
  real det
  real rmat5_det
  logical inside
  real x
  real x1
  real x2
  real x3
  real x4
  real y
  real y1
  real y2
  real y3
  real y4
  real z
  real z1
  real z2
  real z3
  real z4
!
!  Compute the determinant.
!
  a(1,1) = x1
  a(1,2) = y1
  a(1,3) = z1
  a(1,4) = x1 * x1 + y1 * y1 + z1 * z1
  a(1,5) = 1.0E+00

  a(2,1) = x2
  a(2,2) = y2
  a(2,3) = z2
  a(2,4) = x2 * x2 + y2 * y2 + z2 * z2
  a(2,5) = 1.0E+00

  a(3,1) = x3
  a(3,2) = y3
  a(3,3) = z3
  a(3,4) = x3 * x3 + y3 * y3 + z3 * z3
  a(3,5) = 1.0E+00

  a(4,1) = x4
  a(4,2) = y4
  a(4,3) = z4
  a(4,4) = x4 * x4 + y4 * y4 + z4 * z4
  a(4,5) = 1.0E+00

  a(5,1) = x
  a(5,2) = y
  a(5,3) = z
  a(5,4) = x * x + y * y + z * z
  a(5,5) = 1.0E+00

  det = rmat5_det ( a )

  if ( det < 0.0E+00 ) then
    inside = .false.
  else if ( det >= 0.0E+00 ) then
    inside = .true.
  end if

  return
end
subroutine sphere_exp_near_point_3d ( x1, y1, z1, x2, y2, z2, &
  x3, y3, z3, x4, y4, z4, x, y, z, xn, yn, zn )
!
!*******************************************************************************
!
!! SPHERE_EXP_NEAR_POINT_3D finds the nearest point on an explicit sphere to a point in 3D.
!
!
!  Formula:
!
!    An explicit sphere in 3D is determined by four points,
!    which should be distinct, and not coplanar.
!
!  Method:
!
!    If the center of the sphere is (Xc,Yc,Zc), and the point is (X,Y,Z), then
!    the desired point lies at a positive distance R along the vector (X-Xc, Y-Yc, Z-Zc)
!    unless (X,Y,Z) = (Xc,Yc,Zc), in which case any point on the sphere is "nearest".
!
!  Modified:
!
!    14 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
!    (X,Y,Z) coordinates of four points that lie on a circle.
!
!    Input, real X, Y, Z, the (X,Y,Z) coordinates of a point, whose
!    nearest point on the sphere is desired.
!
!    Output, real XN, YN, ZN, the nearest point on the sphere.
!
  implicit none
!
  real enorm0_3d
  real norm
  real r
  real x
  real x1
  real x2
  real x3
  real x4
  real xc
  real xn
  real y
  real y1
  real y2
  real y3
  real y4
  real yc
  real yn
  real z
  real z1
  real z2
  real z3
  real z4
  real zc
  real zn
!
!  Find the center.
!
  call sphere_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, &
    z3, x4, y4, z4, r, xc, yc, zc )
!
!  If (X,Y,Z) = (XC,YC,ZC), bail out now.
!
  norm = enorm0_3d ( x, y, z, xc, yc, zc )

  if ( norm == 0.0E+00 ) then
    xn = xc + r
    yn = yc
    zn = zc
    return
  end if
!
!  Compute the nearest point.
!
  xn = xc + r * ( x - xc ) / norm
  yn = yc + r * ( y - yc ) / norm
  zn = zc + r * ( z - zc ) / norm

  return
end
subroutine sphere_imp_area_3d ( r, area )
!
!*******************************************************************************
!
!! SPHERE_IMP_AREA_3D computes the surface area of an implicit sphere in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Output, real AREA, the area of the sphere.
!
  implicit none
!
  real area
  real r_pi
  real r
!
  area = 4.0E+00 * r_pi ( ) * r * r

  return
end
subroutine sphere_imp_area_nd ( n, r, area )
!
!*******************************************************************************
!
!! SPHERE_IMP_AREA_ND computes the surface area of an implicit sphere in ND.
!
!
!  Discussion:
!
!    N   Area
!
!    2   2       * PI    * R
!    3   4       * PI    * R**2
!    4   2       * PI**2 * R**3
!    5   (8/3)   * PI**2 * R**4
!    6             PI**3 * R**5
!    7   (16/15) * PI**3 * R**6
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Input, real R, the radius of the sphere.
!
!    Output, real AREA, the area of the sphere.
!
  implicit none
!
  real area
  integer n
  real r
!
  call sphere_unit_area_nd ( n, area )

  area = area * r**(n-1)

  return
end
subroutine sphere_imp_cap_area_3d ( r, h, area )
!
!*******************************************************************************
!
!! SPHERE_IMP_CAP_AREA_3D computes the surface area of a spherical cap in 3D.
!
!
!  Definition.
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Consider the point on the radius line which is
!    H units from P.  Draw the circle that lies in the plane perpendicular to
!    the radius, and which intersects the sphere.  The circle divides the sphere
!    into two pieces, and the corresponding disk divides the solid sphere into
!    two pieces.  The spherical cap is the part of the solid sphere that
!    includes the point P.
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real H, the "height" of the spherical cap.  H must be between
!    0 and 2 * R.
!
!    Output, real AREA, the area of the spherical cap.
!
  implicit none
!
  real area
  real h
  real r_pi
  real r
!
  if ( h <= 0.0E+00 ) then
    area = 0.0E+00
  else if ( h >= 2.0E+00 * r ) then
    area = 4.0E+00 * r_pi ( ) * r**2
  else
    area = 2.0E+00 * r_pi ( ) * r * h
  end if

  return
end
subroutine sphere_imp_cap_volume_3d ( r, h, volume )
!
!*******************************************************************************
!
!! SPHERE_IMP_CAP_VOLUME_3D computes the volume of a spherical cap in 3D.
!
!
!  Definition.
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Consider the point on the radius line which is
!    H units from P.  Draw the circle that lies in the plane perpendicular to
!    the radius, and which intersects the sphere.  The circle divides the sphere
!    into two pieces, and the corresponding disk divides the solid sphere into
!    two pieces.  The spherical cap is the part of the solid sphere that includes
!    the point P.
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real H, the "height" of the spherical cap.  H must be between
!    0 and 2 * R.
!
!    Output, real VOLUME, the volume of the spherical cap.
!
  implicit none
!
  real h
  real r_pi
  real r
  real volume
!
  if ( h <= 0.0E+00 ) then
    volume = 0.0E+00
  else if ( h >= 2.0E+00 * r ) then
    volume = ( 4.0E+00 / 3.0E+00 ) * r_pi ( ) * r**3
  else
    volume = ( 1.0E+00 / 3.0E+00 ) * r_pi ( ) * h**2 * ( 3.0E+00 * r - h )
  end if

  return
end
subroutine sphere_imp_contains_point_3d ( r, xc, yc, zc, x, y, z, inside )
!
!*******************************************************************************
!
!! SPHERE_IMP_CONTAINS_POINT_3D determines if an implicit sphere contains a point in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
!
!    Input, real X, Y, Z, the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is inside or on the sphere,
!    FALSE otherwise.
!
  implicit none
!
  real enormsq0_3d
  logical inside
  real r
  real x
  real xc
  real y
  real yc
  real z
  real zc
!
  if ( enormsq0_3d ( x, y, z, xc, yc, zc ) <= r * r ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
subroutine sphere_imp_gridfaces_3d ( maxtri, nlat, nlong, ntri, tri )
!
!*******************************************************************************
!
!! SPHERE_IMP_GRIDFACES_3D produces a grid of triangles on an implicit sphere in 3D.
!
!
!  Discussion:
!
!    The point numbering system is the same used in SPHERE_IMP_GRIDPOINTS_3D,
!    and that routine may be used to compute the coordinates of the points.
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Modified:
!
!    25 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXTRI, the maximum number of triangles.
!
!    Input, integer NLAT, NLONG, the number of latitude and longitude
!    lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so NLAT = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Output, integer NTRI, the number of triangles.
!
!    Output, integer TRI(3,MAXTRI), contains triples of point indices for
!    the triangles that make up the grid.
!
  implicit none
!
  integer maxtri
!
  integer i
  integer j
  integer n
  integer n_max
  integer n_min
  integer ne
  integer nlat
  integer nlong
  integer ntri
  integer nw
  integer s
  integer s_max
  integer s_min
  integer se
  integer sw
  integer tri(3,maxtri)
!
  ntri = 0
!
!  The first row.
!
  n = 1

  sw = 2
  se = sw + 1

  s_min = 2
  s_max = nlong + 1

  do j = 0, nlong - 1

    if ( ntri < maxtri ) then
      ntri = ntri + 1
      tri(1,ntri) = sw
      tri(2,ntri) = se
      tri(3,ntri) = n
    end if

    sw = se

    if ( se == s_max ) then
      se = s_min
    else
      se = se + 1
    end if

  end do
!
!  The intermediate rows.
!
  do i = 1, nlat

    n_max = s_max
    n_min = s_min

    s_max = s_max + nlong
    s_min = s_min + nlong

    nw = n_min
    ne = nw + 1
    sw = s_min
    se = sw + 1

    do j = 0, nlong - 1

      if ( ntri < maxtri ) then
        ntri = ntri + 1
        tri(1,ntri) = sw
        tri(2,ntri) = se
        tri(3,ntri) = nw
      end if

      if ( ntri < maxtri ) then
        ntri = ntri + 1
        tri(1,ntri) = ne
        tri(2,ntri) = nw
        tri(3,ntri) = se
      end if

      sw = se
      nw = ne

      if ( se == s_max ) then
        se = s_min
      else
        se = se + 1
      end if

      if ( ne == n_max ) then
        ne = n_min
      else
        ne = ne + 1
      end if

    end do

  end do
!
!  The last row.
!
  n_max = s_max
  n_min = s_min

  s = n_max + 1

  nw = n_min
  ne = nw + 1

  do j = 0, nlong - 1

    if ( ntri < maxtri ) then
      ntri = ntri + 1
      tri(1,ntri) = ne
      tri(2,ntri) = nw
      tri(3,ntri) = s
    end if

    nw = ne

    if ( ne == n_max ) then
      ne = n_min
    else
      ne = ne + 1
    end if

  end do
  return
end
subroutine sphere_imp_gridlines_3d ( maxline, nlat, nlong, nline, line )
!
!*******************************************************************************
!
!! SPHERE_IMP_GRIDLINES_3D produces "grid lines" on an implicit sphere in 3D.
!
!
!  Discussion:
!
!    The point numbering system is the same used in SPHERE_IMP_GRIDPOINTS_3D,
!    and that routine may be used to compute the coordinates of the points.
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXLINE, the maximum number of gridlines.
!
!    Input, integer NLAT, NLONG, the number of latitude and longitude
!    lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so NLAT = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Output, integer NLINE, the number of grid lines.
!
!    Output, integer LINE(2,MAXLINE), contains pairs of point indices for
!    line segments that make up the grid.
!
  implicit none
!
  integer maxline
!
  integer i
  integer j
  integer line(2,maxline)
  integer new
  integer newcol
  integer nlat
  integer nline
  integer nlong
  integer old
!
  nline = 0
!
!  "Vertical" lines.
!
  do j = 0, nlong - 1

    old = 1
    new = j + 2

    if ( nline < maxline ) then
      nline = nline + 1
      line(1,nline) = old
      line(2,nline) = new
    end if

    do i = 1, nlat - 1

      old = new
      new = old + nlong

      if ( nline < maxline ) then
        nline = nline + 1
        line(1,nline) = old
        line(2,nline) = new
      end if

    end do

    old = new

    if ( nline < maxline ) then
      nline = nline + 1
      line(1,nline) = old
      line(2,nline) = 1 + nlat * nlong + 1
    end if

  end do
!
!  "Horizontal" lines.
!
  do i = 1, nlat

    new = 1 + ( i - 1 ) * nlong + 1

    do j = 0, nlong - 2
      old = new
      new = old + 1
      if ( nline < maxline ) then
        nline = nline + 1
        line(1,nline) = old
        line(2,nline) = new
      end if
    end do

    old = new
    new = 1 + ( i - 1 ) * nlong + 1
    if ( nline < maxline ) then
      line(1,nline) = old
      line(2,nline) = new
    end if

  end do
!
!  "Diagonal" lines.
!
  do j = 0, nlong - 1

    old = 1
    new = j + 2
    newcol = j

    do i = 1, nlat - 1

      old = new
      new = old + nlong + 1

      newcol = newcol + 1
      if ( newcol > nlong - 1 ) then
        newcol = 0
        new = new - nlong
      end if

      if ( nline < maxline ) then
        nline = nline + 1
        line(1,nline) = old
        line(2,nline) = new
      end if

    end do

  end do

  return
end
subroutine sphere_imp_gridpoints_3d ( r, xc, yc, zc, maxpoint, nlat, nlong, &
  npoint, x, y, z )
!
!*******************************************************************************
!
!! SPHERE_IMP_GRIDPOINTS_3D produces "grid points" on an implicit sphere in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Modified:
!
!    22 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
!
!    Input, integer MAXPOINT, the maximum number of grid points, which
!    should be at least 2 + NLAT * NLONG.
!
!    Input, integer NLAT, NLONG, the number of latitude and longitude
!    lines to draw.  The latitudes do not include the North and South
!    poles, which will be included automatically, so NLAT = 5, for instance,
!    will result in points along 7 lines of latitude.
!
!    Output, integer NPOINT, the number of grid points.  The number of
!    grid points depends on N as follows:
!
!      NPOINT = 2 + NLAT * NLONG.
!
!    Output, real X(MAXPOINT), Y(MAXPOINT), Z(MAXPOINT), the coordinates of the
!    grid points.
!
  implicit none
!
  integer maxpoint
!
  integer i
  integer j
  integer nlat
  integer nlong
  integer npoint
  real r_pi
  real pi_const
  real phi
  real r
  real theta
  real x(maxpoint)
  real xc
  real y(maxpoint)
  real yc
  real z(maxpoint)
  real zc
!
  npoint = 0
  pi_const = r_pi ( )
!
!  The north pole.
!
  theta = 0.0E+00
  phi = 0.0E+00
  npoint = npoint + 1
  if ( npoint <= maxpoint ) then
    call sphere_imp_local2xyz_3d ( r, xc, yc, zc, theta, phi, &
      x(npoint), y(npoint), z(npoint) )
  end if
!
!  Do each intermediate ring of latitude.
!
  do i = 1, nlat

    phi = pi_const * real ( i ) / real ( nlat + 1 )
!
!  Along that ring of latitude, compute points at various longitudes.
!
    do j = 0, nlong-1

      theta = 2.0E+00 * pi_const * real ( j ) / real ( nlong )

      npoint = npoint + 1
      if ( npoint <= maxpoint ) then
        call sphere_imp_local2xyz_3d ( r, xc, yc, zc, theta, phi, &
          x(npoint), y(npoint), z(npoint) )
      end if

    end do
  end do
!
!  The south pole.
!
  theta = 0.0E+00
  phi = pi_const
  npoint = npoint + 1
  if ( npoint <= maxpoint ) then
    call sphere_imp_local2xyz_3d ( r, xc, yc, zc, theta, phi, &
      x(npoint), y(npoint), z(npoint) )
  end if

  return
end
subroutine sphere_imp_line_project_3d ( r, xc, yc, zc, npnt, &
  x, y, z, maxpnt2, npnt2, xp, yp, zp, thetamin, thetamax )
!
!*******************************************************************************
!
!! SPHERE_IMP_LINE_PROJECT_3D projects a line onto an implicit sphere in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Discussion:
!
!    The line to be projected is specified as a sequence of points.
!    If two successive points subtend a small angle, then the second
!    point is essentially dropped.  If two successive points subtend
!    a large angle, then intermediate points are inserted, so that
!    the projected line stays closer to the sphere.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.  If R is
!    zero, (XP,YP,ZP) will be returned as (XC,YC,ZC), and if R is
!    negative, points will end up diametrically opposite from where
!    you would expect them for a positive R.
!
!    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
!
!    Input, integer NPNT, the number of points on the line that is
!    to be projected.
!
!    Input, real X(NPNT), Y(NPNT), Z(NPNT), the coordinates of the points
!    on the line that is to be projected.
!    Note: if any (X,Y,Z) coincides with the center of the sphere, then
!    its projection is mathematically undefined.  (XP,YP,ZP) will
!    be returned as (XC,YC,ZC).
!
!    Input, integer MAXPNT2, the maximum number of points on the projected
!    line.  Even if the routine thinks that more points are needed,
!    no more than MAXPNT2 will be generated.
!
!    Output, integer NPNT2, the number of points on the projected line.
!    NPNT2 can be zero, if the line has an angular projection of less
!    than THETAMIN radians.
!
!    Output, real XP(NPNT2), YP(NPNT2), ZP(NPNT2), the coordinates of the
!    points representing the projected line.  These points lie on the
!    sphere.  Successive points are separated by at least THETAMIN
!    radians, and by no more than THETAMAX radians.
!
!    Input, real THETAMIN, THETAMAX, the minimum and maximum angular
!    projections allowed between successive projected points.
!    If two successive points on the original line have projections
!    separated by more than THETAMAX radians, then intermediate points
!    will be inserted, in an attempt to keep the line closer to the
!    sphere.  If two successive points are separated by less than
!    THETAMIN radians, then the second point is dropped, and the
!    line from the first point to the next point is considered.
!
  implicit none
!
  integer maxpnt2
  integer npnt
!
  real alpha
  real ang3d
  real arc_cosine
  real dot
  real enorm0_3d
  integer i
  integer j
  integer nfill
  integer npnt2
  real r
  real thetamax
  real thetamin
  real tnorm
  real x(npnt)
  real x1
  real x2
  real xc
  real xi
  real xp(maxpnt2)
  real y(npnt)
  real y1
  real y2
  real yc
  real yi
  real yp(maxpnt2)
  real z(npnt)
  real z1
  real z2
  real zc
  real zi
  real zp(maxpnt2)
!
!  Check the input.
!
  if ( r == 0.0E+00 ) then
    npnt2 = 0
    return
  end if

  x1 = xc
  y1 = yc
  z1 = zc

  x2 = xc
  y2 = yc
  z2 = zc

  npnt2 = 0

  do i = 1, npnt

    if ( x(i) == xc .and. y(i) == yc .and. z(i) == zc ) then

    else

      x1 = x2
      y1 = y2
      z1 = z2

      alpha = enorm0_3d ( xc, yc, zc, x(i), y(i), z(i) )

      x2 = xc + r * ( x(i) - xc ) / alpha
      y2 = yc + r * ( y(i) - yc ) / alpha
      z2 = zc + r * ( z(i) - zc ) / alpha
!
!  If we haven't gotten any points yet, take this point as our start.
!
      if ( npnt2 == 0 ) then

        npnt2 = npnt2 + 1
        xp(npnt2) = x2
        yp(npnt2) = y2
        zp(npnt2) = z2
!
!  Compute the angular projection of (X1,Y1,Z1) to (X2,Y2,Z2).
!
    else if ( npnt2 >= 1 ) then

      dot = ( x1 - xc ) * ( x2 - xc ) + ( y1 - yc ) * ( y2 - yc ) + &
            ( z1 - zc ) * ( z2 - zc )
      ang3d = arc_cosine (  dot / r**2 )
!
!  If the angle is at least THETAMIN, (or it's the last point),
!  then we will draw a line segment.
!
      if ( abs ( ang3d ) > thetamin .or. i == npnt ) then
!
!  Now we check to see if the line segment is too long.
!
        if ( abs ( ang3d ) > thetamax ) then

          nfill = int ( abs ( ang3d ) / thetamax )

          do j = 1, nfill-1

            xi = ( ( nfill - j ) * ( x1 - xc ) + j * ( x2 - xc ) )
            yi = ( ( nfill - j ) * ( y1 - yc ) + j * ( y2 - yc ) )
            zi = ( ( nfill - j ) * ( z1 - zc ) + j * ( z2 - zc ) )
            tnorm = sqrt ( xi * xi + yi * yi + zi * zi )

            if ( tnorm /= 0.0E+00 ) then
              xi = xc + r * xi / tnorm
              yi = yc + r * yi / tnorm
              zi = zc + r * zi / tnorm
              npnt2 = npnt2 + 1
              xp(npnt2) = xi
              yp(npnt2) = yi
              zp(npnt2) = zi
            end if

          end do

        end if
!
!  Now tack on the projection of point 2.
!
        npnt2 = npnt2 + 1
        xp(npnt2) = x2
        yp(npnt2) = y2
        zp(npnt2) = z2

      end if

    end if

    end if

  end do

  return
end
subroutine sphere_imp_local2xyz_3d ( r, xc, yc, zc, theta, phi, x, y, z )
!
!*******************************************************************************
!
!! SPHERE_IMP_LOCAL2XYZ_3D converts local to XYZ coordinates on an implicit sphere in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!    The "local" spherical coordinates of a point are two angles, THETA and PHI.
!    PHI measures the angle that the vector from the origin to the point
!    makes with the positive Z axis.  THETA measures the angle that the
!    projection of the vector onto the XY plane makes with the positive X axis.
!
!  Modified:
!
!    19 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
!
!    Input, real THETA, PHI, the local (THETA,PHI) spherical coordinates
!    of a point on the sphere.  THETA and PHI are angles, measure in
!    radians.  Usually, 0 <= THETA < 2 * PI, and 0 <= PHI <= PI.
!
!    Output, real X, Y, Z, the XYZ coordinates of the point.
!
  implicit none
!
  real phi
  real r
  real theta
  real x
  real xc
  real y
  real yc
  real z
  real zc
!
  x = xc + r * sin ( phi ) * cos ( theta )
  y = yc + r * sin ( phi ) * sin ( theta )
  z = zc + r * cos ( phi )

  return
end
subroutine sphere_imp_near_point_3d ( r, xc, yc, zc, x, y, z, xn, yn, zn )
!
!*******************************************************************************
!
!! SPHERE_IMP_NEAR_POINT_3D finds the nearest point on an implicit sphere to a point in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Method:
!
!    If the center of the sphere is (Xc,Yc,Zc), and the point is (X,Y,Z), then
!    the desired point lies at a positive distance R along the vector 
!    (X-Xc, Y-Yc, Z-Zc) unless (X,Y,Z) = (Xc,Yc,Zc), in which case any point 
!    on the sphere is "nearest".
!
!  Modified:
!
!    14 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
!
!    Input, real X, Y, Z, the (X,Y,Z) coordinates of a point, whose
!    nearest point on the sphere is desired.
!
!    Output, real XN, YN, ZN, the nearest point on the sphere.
!
  implicit none
!
  real enorm0_3d
  real norm
  real r
  real x
  real xc
  real xn
  real y
  real yc
  real yn
  real z
  real zc
  real zn
!
!  If (X,Y,Z) = (XC,YC,ZC), bail out now.
!
  norm = enorm0_3d ( x, y, z, xc, yc, zc )

  if ( norm == 0.0E+00 ) then
    xn = xc + r
    yn = yc
    zn = zc
    return
  end if
!
!  Compute the nearest point.
!
  xn = xc + r * ( x - xc ) / norm
  yn = yc + r * ( y - yc ) / norm
  zn = zc + r * ( z - zc ) / norm

  return
end
subroutine sphere_imp_point_project_3d ( r, xc, yc, zc, x, y, z, xp, yp, zp )
!
!*******************************************************************************
!
!! SPHERE_IMP_POINT_PROJECT_3D projects a point onto an implicit sphere, in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
!
!    Input, real X, Y, Z, the coordinates of a point.
!
!    Output, real X, Y, Z, the coordinates of the point as projected
!    onto the sphere from the center.
!
  implicit none
!
  real enorm0_3d
  real norm
  real r
  real x
  real xc
  real xp
  real y
  real yc
  real yp
  real z
  real zc
  real zp
!
  if ( r == 0.0E+00 ) then

    xp = xc
    yp = yc
    zp = zc

  else if ( x == xc .and. y == yc .and. z == zc ) then

    xp = xc
    yp = yc
    zp = zc + r

  else

    norm = enorm0_3d ( x, y, z, xc, yc, zc )

    xp = xc + r * ( x - xc ) / norm
    yp = yc + r * ( y - yc ) / norm
    zp = zc + r * ( z - zc ) / norm

  end if

  return
end
subroutine sphere_imp_spiralpoints_3d ( r, xc, yc, zc, n, x, y, z )
!
!*******************************************************************************
!
!! SPHERE_IMP_SPIRALPOINTS_3D produces spiral points on an implicit sphere in 3D.
!
!
!  Discussion:
!
!    The points should be arranged on the sphere in a pleasing design.
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Reference:
!
!    E B Saff and A B J Kuijlaars,
!    Distributing Many Points on a Sphere,
!    The Mathematical Intelligencer,
!    Volume 19, Number 1, 1997, pages 5-11.
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
!
!    Input, integer N, the number of points to create.
!
!    Output, real X(N), Y(N), Z(N), the coordinates of the grid points.
!
  implicit none
!
  integer n
!
  real cosphi
  integer i
  real r
  real sinphi
  real theta
  real, parameter :: two_pi = &
    2.0E+00 * 3.14159265358979323846264338327950288419716939937510E+00
  real x(n)
  real xc
  real y(n)
  real yc
  real z(n)
  real zc
!
  do i = 1, n

    cosphi = ( real ( n - i ) * ( -1.0E+00 ) &
             + real ( i - 1 ) * ( +1.0E+00 ) ) / real ( n - 1 )

    sinphi = sqrt ( 1.0E+00 - cosphi**2 )

    if ( i == 1 .or. i == n ) then
      theta = 0.0E+00
    else
      theta = theta + 3.6E+00 / ( sinphi * sqrt ( real ( n ) ) )
      theta = mod ( theta, two_pi )
    end if

    x(i) = xc + r * sinphi * cos ( theta )
    y(i) = yc + r * sinphi * sin ( theta )
    z(i) = zc + r * cosphi

  end do

  return
end
subroutine sphere_imp_volume_3d ( r, volume )
!
!*******************************************************************************
!
!! SPHERE_IMP_VOLUME_3D computes the volume of an implicit sphere in 3D.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Modified:
!
!    26 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Output, real VOLUME, the volume of the sphere.
!
  implicit none
!
  real r_pi
  real r
  real volume
!
  volume = ( 4.0E+00 / 3.0E+00 ) * r_pi ( ) * r**3

  return
end
subroutine sphere_imp_volume_nd ( n, r, volume )
!
!*******************************************************************************
!
!! SPHERE_IMP_VOLUME_ND computes the volume of an implicit sphere in ND.
!
!
!  Formula:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
!
!  Discussion:
!
!    N  Volume
!
!    2             PI    * R**2
!    3  (4/3)    * PI    * R**3
!    4  (1/2)    * PI**2 * R**4
!    5  (8/15)   * PI**2 * R**5
!    6  (1/6)    * PI**3 * R**6
!    7  (16/105) * PI**3 * R**7
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Input, real R, the radius of the sphere.
!
!    Output, real VOLUME, the volume of the sphere.
!
  implicit none
!
  integer n
  real r
  real volume
!
  call sphere_unit_volume_nd ( n, volume )

  volume = volume * r**n

  return
end
subroutine sphere_imp_zone_area_3d ( r, h1, h2, area  )
!
!*******************************************************************************
!
!! SPHERE_IMP_ZONE_AREA_3D computes the surface area of a spherical zone in 3D.
!
!
!  Definition.
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Now choose two points on the radius line, a
!    distance H1 and H2 from the point P.  Consider all the points on or within
!    the sphere whose projection onto the radius lies between these two points.
!    These points constitute the spherical zone, which can also be considered
!    the difference of two spherical caps.
!
!  Modified:
!
!    02 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real H1, H2, the distances that define the thickness of the zone.
!    H1 and H2 must be between 0 and 2 * R.
!
!    Output, real AREA, the area of the spherical zone.
!
  implicit none
!
  real area
  real h
  real h1
  real h2
  real r_pi
  real r
!
  h = abs ( h1 - h2 )

  if ( h <= 0.0E+00 ) then
    area = 0.0E+00
  else if ( h >= 2.0E+00 * r ) then
    area = 4.0E+00 * r_pi ( ) * r**2
  else
    area = 2.0E+00 * r_pi ( ) * r * h
  end if

  return
end
subroutine sphere_imp_zone_volume_3d ( r, h1, h2, volume )
!
!*******************************************************************************
!
!! SPHERE_IMP_ZONE_VOLUME_3D computes the volume of a spherical zone in 3D.
!
!
!  Definition.
!
!    Draw any radius of the sphere and note the point P where the radius
!    intersects the sphere.  Now choose two points on the radius line, a
!    distance H1 and H2 from the point P.  Consider all the points on or within
!    the sphere whose projection onto the radius lies between these two points.
!    These points constitute the spherical zone, which can also be considered
!    the difference of two spherical caps.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real H1, H2, the distances that define the thickness of the zone.
!    H1 and H2 must be between 0 and 2 * R.
!
!    Output, real VOLUME, the volume of the spherical zone
!
  implicit none
!
  real h1
  real h11
  real h2
  real h22
  real r_pi
  real r
  real volume
!
  h11 = min ( h1, h2 )
  h11 = max ( h11, 0.0E+00 )

  if ( h11 >= 2.0E+00 * r ) then
    volume = 0.0E+00
    return
  end if

  h22 = max ( h1, h2 )
  h22 = min ( h22, 2.0E+00 * r )

  if ( h22 <= 0.0E+00 ) then
    volume = 0.0E+00
    return
  end if

  volume = ( 1.0E+00 / 3.0E+00 ) * r_pi ( ) * ( &
      h22**2 * ( 3.0E+00 * r - h22 ) &
    - h11**2 * ( 3.0E+00 * r - h11 ) )

  return
end
subroutine sphere_unit_area_nd ( n, area )
!
!*******************************************************************************
!
!! SPHERE_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
!
!
!  Discussion:
!
!    N   Area
!
!    2   2       * PI
!    3   4       * PI
!    4   2       * PI**2
!    5   (8/3)   * PI**2
!    6             PI**3
!    7   (16/15) * PI**3
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Output, real AREA, the area of the sphere.
!
  implicit none
!
  real area
  integer i
  integer m
  integer n
  real r_pi
!
  if ( mod ( n, 2 ) == 0 ) then
    m = n / 2
    area = 2.0E+00 * ( r_pi ( ) )**m
    do i = 1, m-1
      area = area / real ( i )
    end do
  else
    m = ( n - 1 ) / 2
    area = 2.0E+00**n * ( r_pi ( ) )**m
    do i = m+1, 2*m
      area = area / real ( i )
    end do
  end if

  return
end
subroutine sphere_unit_sample_2d ( x )
!
!*******************************************************************************
!
!! SPHERE_UNIT_SAMPLE_2D picks a random point on the unit sphere (circle) in 2D.
!
!
!  Modified:
!
!    15 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real X(2), the random point on the unit circle.
!
  implicit none
!
  real r_pi
  real u
  real x(2)
!
  call random_number ( harvest = u )

  x(1) = cos ( 2.0E+00 * r_pi ( ) * u )
  x(2) = sin ( 2.0E+00 * r_pi ( ) * u )

  return
end
subroutine sphere_unit_sample_3d ( x )
!
!*******************************************************************************
!
!! SPHERE_UNIT_SAMPLE_3D picks a random point on the unit sphere in 3D.
!
!
!  Modified:
!
!    03 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real X(3), the sample point.
!
  implicit none
!
  real arc_cosine
  integer i
  real phi
  real r_pi
  real theta
  real vdot
  real x(3)
!
!  Pick a uniformly random VDOT, which must be between -1 and 1.
!  This represents the dot product of the random vector with the Z unit vector.
!
!  Note: this works because the surface area of the sphere between
!  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
!  a patch of area uniformly.
!
  call random_number ( harvest = vdot )
  vdot = 2.0E+00 * vdot - 1.0E+00

  phi = arc_cosine ( vdot )
!
!  Pick a uniformly random rotation between 0 and 2 Pi around the
!  axis of the Z vector.
!
  call random_number ( harvest = theta )
  theta = 2.0E+00 * r_pi ( ) * theta

  x(1) = cos ( theta ) * sin ( phi )
  x(2) = sin ( theta ) * sin ( phi )
  x(3) = cos ( phi )

  return
end
subroutine sphere_unit_sample2_3d ( x )
!
!*******************************************************************************
!
!! SPHERE_UNIT_SAMPLE2_3D is a BAD method for sampling the unit sphere in 3D.
!
!
!  Discussion:
!
!    Points on the unit sphere have coordinates ( PHI, THETA ) where
!    PHI varies from 0 to PI, and THETA from 0 to 2 PI, so that:
!
!    x = cos ( theta ) * sin ( phi )
!    y = sin ( theta ) * sin ( phi )
!    z =                 cos ( phi )
!
!    This routine implements a sampling of the sphere that simply
!    picks PHI and THETA uniformly at random from their ranges.
!    This is a uniform sampling on the cylinder, but it is NOT
!    a uniform sampling on the sphere.  I implement it here just
!    so I can run some tests against the code in SPHERE_UNIT_SAMPLE_3D.
!
!  Modified:
!
!    17 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real X(3), the sample point.
!
  implicit none
!
  real phi
  real r_pi
  real theta
  real x(3)
!
  call random_number ( harvest = phi )
  phi = r_pi ( ) * phi

  call random_number ( harvest = theta )
  theta = 2.0E+00 * r_pi ( ) * theta

  x(1) = cos ( theta ) * sin ( phi )
  x(2) = sin ( theta ) * sin ( phi )
  x(3) = cos ( phi )

  return
end
subroutine sphere_unit_sample_nd ( n, x )
!
!*******************************************************************************
!
!! SPHERE_UNIT_SAMPLE_ND picks a random point on the unit sphere in ND.
!
!
!  Discussion:
!
!    N-1 random Givens rotations are applied to the point ( 1, 0, 0, ..., 0 ).
!
!    The I-th Givens rotation is in the plane of coordinate axes I and I+1,
!    and has the form:
!
!     [ cos ( theta )  - sin ( theta ) ] * x(i)      = x'(i)
!     [ sin ( theta )    cos ( theta ) ]   x(i+1)      x'(i+1)
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Output, real X(N), the random point.
!
  implicit none
!
  integer n
!
  integer i
  real random_cosine(1:n-1)
  real random_sign(1:n-1)
  real random_sine(1:n-1)
  real x(n)
  real xi
!
  x(1) = 1.0E+00
  x(2:n) = 0.0E+00

  call random_number ( harvest = random_cosine(1:n-1) )
  random_cosine(1:n-1) = 2.0E+00 * random_cosine(1:n-1) - 1.0E+00

  call random_number ( harvest = random_sign(1:n-1) )
  random_sign(1:n-1) = real ( 2 * int ( 2.0E+00 * random_sign(1:n-1) ) - 1 )

  random_sine(1:n-1) = &
    random_sign(1:n-1) * sqrt ( 1.0E+00 - random_cosine(1:n-1)**2 )

  do i = 1, n-1
    xi = x(i)
    x(i  ) = random_cosine(i) * xi
    x(i+1) = random_sine(i)   * xi
  end do

  return
end
subroutine sphere_unit_sample2_nd ( n, x )
!
!*******************************************************************************
!
!! SPHERE_UNIT_SAMPLE2_ND picks a random point on the unit sphere in ND.
!
!
!  Discussion:
!
!    N independent normally distributed random numbers are generated,
!    and then scaled to have unit norm.
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Output, real X(N), the random point.
!
  implicit none
!
  integer n
!
  real norm
  real x(n)
!
  call normal_01_vector ( n, x )

  norm = sqrt ( sum ( x(1:n)**2 ) )

  x(1:n) = x(1:n) / norm

  return
end
subroutine sphere_unit_sample3_nd ( n, x )
!
!*******************************************************************************
!
!! SPHERE_UNIT_SAMPLE3_ND picks a random point on the unit sphere in ND.
!
!
!  Discussion:
!
!    Points in the [-1,1] cube are generated.  Points lying outside
!    the sphere are rejected.  Points inside the unit sphere are normalized
!    to lie on the sphere.
!
!    Because the volume of the unit sphere
!    relative to the unit cube decreases drastically in higher dimensions,
!    this routine becomes increasingly inefficient at higher N.  
!    Above N = 5, this problem will become significant.
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Output, real X(N), the random point.
!
  implicit none
!
  integer n
!
  real norm
  real x(n)
!
  do

    call random_number ( harvest = x(1:n) )
    x(1:n) = 2.0E+00 * x(1:n) - 1.0E+00

    norm = sqrt ( sum ( x(1:n)**2 ) )

    if ( norm <= 1.0E00 ) then
      x(1:n) = x(1:n) / norm
      exit
    end if

  end do

  return
end
subroutine sphere_unit_volume_nd ( n, volume )
!
!*******************************************************************************
!
!! SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
!
!
!  Discussion:
!
!    N  Volume
!
!    2             PI
!    3  (4/3)    * PI
!    4  (1/2)    * PI**2
!    5  (8/15)   * PI**2
!    6  (1/6)    * PI**3
!    7  (16/105) * PI**3
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the space.
!
!    Output, real VOLUME, the volume of the sphere.
!
  implicit none
!
  integer i
  integer m
  integer n
  real r_pi
  real volume
!
  if ( mod ( n, 2 ) == 0 ) then
    m = n / 2
    volume = ( r_pi ( ) )**m
    do i = 1, m
      volume = volume / real ( i )
    end do
  else
    m = ( n - 1 ) / 2
    volume = ( r_pi ( ) )**m * 2.0E+00**n
    do i = m+1, 2*m+1
      volume = volume / real ( i )
    end do
  end if

  return
end
subroutine stri_angles_to_area_3d ( r, a, b, c, area )
!
!*******************************************************************************
!
!! STRI_ANGLES_TO_AREA_3D computes the area of a spherical triangle.
!
!
!  Formula:
!
!    A sphere in 3D satisfies the equation:
!
!      X**2 + Y**2 + Z**2 = R**2
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI ) * R**2
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real A, B, C, the angles of the triangle.
!
!    Output, real AREA, the area of the sphere.
!
  implicit none
!
  real area
  real a
  real b
  real c
  real r
  real r_pi
!
!  Apply Girard's formula.
!
  area = r**2 * ( a + b + c - r_pi ( ) )

  return
end
subroutine stri_sides_to_angles_3d ( r, as, bs, cs, a, b, c )
!
!*******************************************************************************
!
!! STRI_SIDES_TO_ANGLES_3D computes spherical triangle angles in 3D.
!
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real AS, BS, CS, the (geodesic) length of the sides of the
!    triangle.
!
!    Output, real A, B, C, the spherical angles of the triangle.
!    Angle A is opposite the side of length AS, and so on.
!
  implicit none
!
  real a
  real as
  real asu
  real b
  real bs
  real bsu
  real c
  real cs
  real csu
  real r
  real ssu
  real tan_a2
  real tan_b2
  real tan_c2
!
  asu = as / r
  bsu = bs / r
  csu = cs / r
  ssu = ( asu + bsu + csu ) / 2.0E+00

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - asu )     ) )

  a = 2.0E+00 * atan ( tan_a2 )

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) )

  b = 2.0E+00 * atan ( tan_b2 )

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / & 
                  ( sin ( ssu ) * sin ( ssu - csu )     ) )

  c = 2.0E+00 * atan ( tan_c2 )

  return
end
subroutine stri_vertices_to_area_3d ( r, v1, v2, v3, area )
!
!*******************************************************************************
!
!! STRI_VERTICES_TO_AREA_3D computes the area of a spherical triangle in 3D.
!
!
!  Formula:
!
!    A sphere in 3D satisfies the equation:
!
!      X**2 + Y**2 + Z**2 = R**2
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI ) * R**2
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real AREA, the area of the sphere.
!
  implicit none
!
  real area
  real a
  real as
  real b
  real bs
  real c
  real cs
  real r
  real v1(3)
  real v2(3)
  real v3(3)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call stri_vertices_to_sides ( r, v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call stri_sides_to_angles_3d ( r, as, bs, cs, a, b, c )
!
!  Get the area
!
  call stri_angles_to_area_3d ( r, a, b, c, area )

  return
end
subroutine stri_vertices_to_centroid_3d ( r, v1, v2, v3, vs )
!
!*******************************************************************************
!
!! STRI_VERTICES_TO_CENTROID_3D gets a spherical triangle centroid in 3D.
!
!
!  Formula:
!
!    A sphere in 3D satisfies the equation:
!
!      X**2 + Y**2 + Z**2 = R**2
!
!    A spherical triangle is specified by three points on the sphere.
!
!    The (true) centroid of a spherical triangle is the point
!
!      VT = (XT,YT,ZT) = Integral ( X, Y, Z ) dArea / Integral 1 dArea
!
!    Note that the true centroid does NOT, in general, lie on the sphere.  
!
!    The "flat" centroid VF is the centroid of the planar triangle defined by
!    the vertices of the spherical triangle.
!
!    The "spherical" centroid VS of a spherical triangle is computed by
!    the intersection of the geodesic bisectors of the triangle angles.
!    The spherical centroid lies on the sphere.
!
!    VF, VT and VS lie on a line through the center of the sphere.  We can
!    easily calculate VF by averaging the vertices, and from this determine
!    VS by normalizing.
!
!    (Of course, we still will not have actually computed VT, which lies
!    somewhere between VF and VS!)
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real VS(3), the coordinates of the "spherical centroid"
!    of the spherical triangle.
!
  implicit none
!
  real norm
  real r
  real v1(3)
  real v2(3)
  real v3(3)
  real vs(3)
!
  vs(1:3) = ( v1(1:3) + v2(1:3) + v3(1:3) ) / 3.0E+00

  norm = sqrt ( sum ( vs(1:3)**2 ) )

  vs(1:3) = r * vs(1:3) / norm

  return
end
subroutine stri_vertices_to_sides ( r, v1, v2, v3, as, bs, cs )
!
!*******************************************************************************
!
!! STRI_VERTICES_TO_SIDES_3D computes spherical triangle sides in 3D.
!
!
!  Discussion:
!
!    We can use the ACOS system call here, but the ARC_COSINE routine
!    will automatically take care of cases where the input argument is
!    (usually slightly) out of bounds.
!
!  Modified:
!
!    21 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the radius of the sphere.
!
!    Input, real V1(3), V2(3), V3(3), the vertices of the spherical
!    triangle.
!
!    Output, real AS, BS, CS, the (geodesic) length of the sides of the
!    triangle.
!
  implicit none
!
  real arc_cosine
  real as
  real bs
  real cs
  real r
  real v1(3)
  real v2(3)
  real v3(3)
!
  as = r * arc_cosine ( dot_product ( v2(1:3), v3(1:3) ) / r**2 )
  bs = r * arc_cosine ( dot_product ( v3(1:3), v1(1:3) ) / r**2 )
  cs = r * arc_cosine ( dot_product ( v1(1:3), v2(1:3) ) / r**2 )

  return
end
subroutine string_2d ( iorder, istrng, nstrng, nvec, x1vec, x2vec, y1vec, y2vec)
!
!*******************************************************************************
!
!! STRING_2D groups line segments into connected lines in 2D.
!
!
!  Discussion:
!
!    The routine receives an unordered set of line segments, described by
!    pairs of coordinates (X1,Y1) and (X2,Y2), and tries to group them
!    into ordered lists that constitute connected jagged lines.
!
!    This routine will not match two endpoints unless they are exactly equal.
!
!  Modified:
!
!    28 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IORDER(NVEC).
!
!    If two segments I and J have the same value of ISTRNG, then
!    they belong to the same string.  Then IORDER(I) and IORDER(J)
!    give the relative order of the two segments in the string.
!    Thus if IORDER(I) = -3 and IORDER(J) = 2, then when the string
!    is traversed, segment I is traversed first, then four other
!    segments are traversed, and then segment J is traversed.
!
!    The segment with IORDER(I) = 0 is the initial segment from which
!    the entire string was "grown".
!
!    Output, integer ISTRNG(NVEC), ISTRNG(I) is the number of the string
!    to which line segment I belongs, between 1 and NVEC.
!
!    Output, integer NSTRNG, the number of strings created.
!
!    Input, integer NVEC, the number of line segments to be analyzed.
!
!    Input/output, real X1VEC(NVEC), X2VEC(NVEC), Y1VEC(NVEC), Y2VEC(NVEC).
!
!    On input, line segment I has endpoints ( X1VEC(I), Y1VEC(I) )
!    and ( X2VEC(I), Y2VEC(I) ).
!
!    On output, the order of the components may have been
!    switched.  That is, for some I, ( X1VEC(I), Y1VEC(I) ) and
!    ( X2VEC(I), Y2VEC(I) ) may have been swapped.
!
!    More importantly, all the entries X1VEC(I), Y1VEC(I),
!    Y1VEC(I) and Y2VEC(I) may have been swapped with another index J.
!
!    The resulting coordinates will have been sorted in order
!    of the string to which they belong, and then by the order
!    of their traversal within that string.
!
  implicit none
!
  integer nvec
!
  integer i
  integer indx
  integer iorder(nvec)
  integer iseed
  integer isgn
  integer istrng(nvec)
  integer itemp
  integer j
  integer jval
  integer kval
  integer match
  integer nstrng
  real temp
  real x1val
  real x1vec(nvec)
  real x2val
  real x2vec(nvec)
  real y1val
  real y1vec(nvec)
  real y2val
  real y2vec(nvec)
!
!  Mark ISTRNG so that each segment is alone.
!
  iorder(1:nvec) = 0
  istrng(1:nvec) = nvec + i
!
!  Starting with the lowest numbered group of line segments,
!  see if any higher numbered groups belong.
!
  iseed = 1
  nstrng = 1
  istrng(iseed) = nstrng

  do

    x1val = x1vec(iseed)
    x2val = x2vec(iseed)
    y1val = y1vec(iseed)
    y2val = y2vec(iseed)

    jval = iorder(iseed)
    kval = iorder(iseed)

    do

      match = 0

      do j = 1, nvec

        if ( istrng(j) > nstrng ) then

          if ( x1val == x1vec(j) .and. y1val == y1vec(j) ) then

            jval = jval - 1
            iorder(j) = jval
            istrng(j) = nstrng
            x1val = x2vec(j)
            y1val = y2vec(j)
            match = match + 1

            call r_swap ( x1vec(j), x2vec(j) )
            call r_swap ( y1vec(j), y2vec(j) )

          else if ( x1val == x2vec(j) .and. y1val == y2vec(j) ) then

            jval = jval - 1
            iorder(j) = jval
            istrng(j) = nstrng
            x1val = x1vec(j)
            y1val = y1vec(j)
            match = match + 1

          else if ( x2val == x1vec(j) .and. y2val == y1vec(j) ) then

            kval = kval + 1
            iorder(j) = kval
            istrng(j) = nstrng
            x2val = x2vec(j)
            y2val = y2vec(j)
            match = match + 1

          else if ( x2val == x2vec(j) .and. y2val == y2vec(j) ) then

            kval = kval + 1
            iorder(j) = kval
            istrng(j) = nstrng
            x2val = x1vec(j)
            y2val = y1vec(j)
            match = match + 1

            call r_swap ( x1vec(j), x2vec(j) )
            call r_swap ( y1vec(j), y2vec(j) )

          end if

        end if

      end do
!
!  If the string has closed on itself, then we don't want to
!  look for any more matches for this string.
!
      if ( x1val == x2val .and. y1val == y2val ) then
        exit
      end if
!
!  If we made no matches this pass, we're done.
!
      if ( match <= 0 ) then
        exit
      end if

    end do
!
!  This string is "exhausted".  Are there any line segments we
!  haven't looked at yet?
!
    iseed = 0

    do i = 1, nvec
      if ( istrng(i) > nstrng ) then
        iseed = i
        nstrng = nstrng + 1
        istrng(i) = nstrng
        exit
      end if
    end do

    if ( iseed == 0 ) then
      exit
    end if

  end do
!
!  There are no more line segments to look at.  Renumber the
!  isolated segments.
!
!  Question: Can this ever happen?
!
  do i = 1, nvec
    if ( istrng(i) > nvec ) then
      nstrng = nstrng + 1
      istrng(i) = nstrng
    end if
  end do
!
!  Now sort the line segments by string and by order of traversal.
!
  i = 0
  isgn = 0
  j = 0

  indx = 0

  do

    call sort_heap_external ( nvec, indx, i, j, isgn )

    if ( indx > 0 ) then

      call i_swap ( iorder(i), iorder(j) )
      call i_swap ( istrng(i), istrng(j) )
      call r_swap ( x1vec(i), x1vec(j) )
      call r_swap ( y1vec(i), y1vec(j) )
      call r_swap ( x2vec(i), x2vec(j) )
      call r_swap ( y2vec(i), y2vec(j) )

    else if ( indx < 0 ) then

      if ( ( istrng(i) < istrng(j) ) .or. &
           ( istrng(i) == istrng(j) .and. iorder(i) < iorder(j) ) ) then

        isgn = - 1

      else

        isgn = + 1

      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine super_ellipse_points_2d ( x0, y0, r1, r2, expo, psi, n, x, y )
!
!*******************************************************************************
!
!! SUPER_ELLIPSE_POINTS_2D returns N points on a tilted superellipse in 2D.
!
!
!  Discussion:
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter.
!
!    The parametric formula of the (untilted) superellipse is:
!
!      X = R1 * cos**EXPO ( THETA )
!      Y = R2 * sin**EXPO ( THETA )
!
!    An implicit form of the (untilted) superellipse is:
!
!      (X/R1)**(2/EXPO) + (Y/R2)**(2/EXPO) = 1
!   
!  Reference:
!
!    Martin Gardner,
!    The Mathematical Carnival,
!    Knopf, 1975, pages 240-254.
! 
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X0, Y0, the coordinates of the center of the superellipse.
!
!    Input, real R1, R2, the "radius" of the superellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    Input, real EXPO, the exponent of the superellipse. 
!    0 = a rectangle;
!    between 0 and 1, a "rounded" rectangle;
!    1.0 = an ellipse;
!    2.0 = a diamond;
!    > 2.0 a pinched shape.
!
!    Input, real PSI, the angle that the major axis of the superellipse
!    makes with the X axis.  A value of 0.0 means that the major and
!    minor axes of the superellipse will be the X and Y coordinate axes.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real X(N), Y(N), the coordinates of points on the superellipse.
!
  implicit none
!
  integer n
!
  real act
  real ast
  integer i
  real expo
  real r_pi
  real psi
  real r1
  real r2
  real sct
  real sst
  real theta
  real x0
  real x(n)
  real y0
  real y(n)
!
  do i = 1, n

    theta = ( 2.0E+00 * r_pi ( ) * real ( i - 1 ) ) / real ( n )

    act = abs ( cos ( theta ) )
    sct = sign ( 1.0E+00, cos ( theta ) )
    ast = abs ( sin ( theta ) )
    sst = sign ( 1.0E+00, sin ( theta ) )

    x(i) = x0 + r1 * cos ( psi ) * sct * ( act )**expo &
              - r2 * sin ( psi ) * sst * ( ast )**expo

    y(i) = y0 + r1 * sin ( psi ) * sct * ( act )**expo &
              + r2 * cos ( psi ) * sst * ( ast )**expo

  end do

  return
end
subroutine swapec ( i, top, btri, bedg, point_num, point_xy, tri_num, &
  tri_vert, tri_nabe, stack, ierr )
!
!*******************************************************************************
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to the triangulation.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer I, the index of the new vertex.
!
!    Input/output, integer TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input/output, integer BTRI, BEDG; on input, if positive, are the
!    triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, intger POINT_NUM, the number of points.
!
!    Input, real POINT_XY(2,POINT_NUM), the coordinates of the points.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input/output, integer TRI_VERT(3,TRI_NUM), the triangle incidence list.  
!    May be updated on output because of swaps.
!
!    Input/output, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list; 
!    negative values are used for links of the counter-clockwise linked 
!    list of boundary edges;  May be updated on output because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer IERR is set to 8 for abnormal return.
!
  implicit none
!
  integer point_num
  integer tri_num
!
  integer a
  integer b
  integer bedg
  integer btri
  integer c
  integer diaedg
  integer e
  integer ee
  integer em1
  integer ep1
  integer f
  integer fm1
  integer fp1
  integer i
  integer ierr
  integer i_wrap
  integer l
  integer r
  integer s
  integer stack(point_num)
  integer swap
  integer t
  integer top
  integer tri_nabe(3,tri_num)
  integer tri_vert(3,tri_num)
  integer tt
  integer u
  real point_xy(2,point_num)
  real x
  real y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  x = point_xy(1,i)
  y = point_xy(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( tri_vert(1,t) == i ) then
      e = 2
      b = tri_vert(3,t)
    else if ( tri_vert(2,t) == i ) then
      e = 3
      b = tri_vert(1,t)
    else
      e = 1
      b = tri_vert(2,t)
    end if

    a = tri_vert(e,t)
    u = tri_nabe(e,t)

    if ( tri_nabe(1,u) == t ) then
      f = 1
      c = tri_vert(3,u)
    else if ( tri_nabe(2,u) == t ) then
      f = 2
      c = tri_vert(1,u)
    else
      f = 3
      c = tri_vert(2,u)
    end if

    swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
      point_xy(2,c), point_xy(1,b), point_xy(2,b) )

    if ( swap == 1 ) then

      em1 = i_wrap ( e - 1, 1, 3 )
      ep1 = i_wrap ( e + 1, 1, 3 )
      fm1 = i_wrap ( f - 1, 1, 3 )
      fp1 = i_wrap ( f + 1, 1, 3 )

      tri_vert(ep1,t) = c
      tri_vert(fp1,u) = i
      r = tri_nabe(ep1,t)
      s = tri_nabe(fp1,u)
      tri_nabe(ep1,t) = u
      tri_nabe(fp1,u) = t
      tri_nabe(e,t) = s
      tri_nabe(f,u) = r

      if ( tri_nabe(fm1,u) > 0 ) then
        top = top + 1
        stack(top) = u
      end if

      if ( s > 0 ) then

        if ( tri_nabe(1,s) == u ) then
          tri_nabe(1,s) = t
        else if ( tri_nabe(2,s) == u ) then
          tri_nabe(2,s) = t
        else
          tri_nabe(3,s) = t
        end if

        top = top + 1

        if ( top > point_num ) then
          ierr = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( tri_nabe(ee,tt) > 0 )

          tt = tri_nabe(ee,tt)

          if ( tri_vert(1,tt) == a ) then
            ee = 3
          else if ( tri_vert(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tri_nabe(ee,tt) = l

      end if

      if ( r > 0 ) then

        if ( tri_nabe(1,r) == t ) then
          tri_nabe(1,r) = u
        else if ( tri_nabe(2,r) == t ) then
          tri_nabe(2,r) = u
        else
          tri_nabe(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( tri_nabe(ee,tt) > 0 )

          tt = tri_nabe(ee,tt)

          if ( tri_vert(1,tt) == b ) then
            ee = 3
          else if ( tri_vert(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tri_nabe(ee,tt) = l

      end if

    end if

  end do

  return
end
function tan_degrees ( angle )
!
!*******************************************************************************
!
!! TAN_DEGREES returns the tangent of an angle given in degrees.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ANGLE, the angle, in degrees.
!
!    Output, real TAN_DEGREES, the tangent of the angle.
!
  implicit none
!
  real angle
  real angle_rad
  real degrees_to_radians
  real tan_degrees
!
  angle_rad = degrees_to_radians ( angle )
  tan_degrees  = sin ( angle_rad ) / cos ( angle_rad )

  return
end
subroutine tetra_barycentric_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, x, y, z, c )
!
!*******************************************************************************
!
!! TETRA_BARYCENTRIC_3D returns the barycentric coordinates of a point in 3D.
!
!
!  Discussion:
!
!    The barycentric coordinates of a point (X,Y,Z) with respect to
!    a tetrahedron are a set of four values C(1:4), each associated
!    with a vertex of the tetrahedron.  The values must sum to 1.
!    If all the values are between 0 and 1, the point is contained
!    within the tetrahedron.
!
!    The barycentric coordinate of point X related to vertex A can be
!    interpreted as the ratio of the volume of the tetrahedron with 
!    vertex A replaced by vertex X to the volume of the original 
!    tetrahedron.
!
!
!  Modified:
!
!    20 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, Z4, Y4, Z4, the
!    coordinates of the tetrahedron vertices.
!
!    Input, real X, Y, Z, the point to be checked.
!
!    Output, real C(4), the barycentric coordinates of (X,Y,Z) with
!    respect to the tetrahedron.
!
  implicit none
!
  integer, parameter :: n = 3
  integer, parameter :: nrhs = 1
!
  real a(n,n+nrhs)
  real c(4)
  integer info
  real x
  real x1
  real x2
  real x3
  real x4
  real y
  real y1
  real y2
  real y3
  real y4
  real z
  real z1
  real z2
  real z3
  real z4
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
!    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
!    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
!
!  which is satisfied by the barycentric coordinates of (X,Y,Z).
!
  a(1,1) = x2 - x1
  a(1,2) = x3 - x1
  a(1,3) = x4 - x1
  a(1,4) = x - x1

  a(2,1) = y2 - y1
  a(2,2) = y3 - y1
  a(2,3) = y4 - y1
  a(2,4) = y - y1

  a(3,1) = z2 - z1
  a(3,2) = z3 - z1
  a(3,3) = z4 - z1
  a(3,4) = z - z1
!
!  Solve the linear system.
!
  call rmat_solve ( a, n, nrhs, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRA_BARYCENTRIC_3D - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper tetrahedron.'
    stop
  end if

  c(1) = a(1,4)
  c(2) = a(2,4)
  c(3) = a(3,4)
  c(4) = 1.0E+00 - c(1) - c(2) - c(3)

  return
end
subroutine tetra_centroid_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, xc, yc, zc )
!
!*******************************************************************************
!
!! TETRA_CENTROID_3D computes the centroid of a tetrahedron in 3D.
!
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
!    coordinates of the vertices.
!
!    Output, real XC, YC, ZC, the coordinates of the centroid.
!
  implicit none
!
  real x1
  real x2
  real x3
  real x4
  real xc
  real y1
  real y2
  real y3
  real y4
  real yc
  real z1
  real z2
  real z3
  real z4
  real zc
!
  xc = 0.25E+00 * ( x1 + x2 + x3 + x4 )
  yc = 0.25E+00 * ( y1 + y2 + y3 + y4 )
  zc = 0.25E+00 * ( z1 + z2 + z3 + z4 )

  return
end
subroutine tetra_circumsphere_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, r, xc, yc, zc )
!
!*******************************************************************************
!
!! TETRA_CIRCUMSPHERE_3D computes the circumscribed sphere of a tetrahedron in 3D.
!
!
!  Discussion:
!
!    The circumscribed sphere of a tetrahedron is the sphere that
!    passes through the four vertices.  The circumsphere is not necessarily
!    the smallest sphere that contains the tetrahedron.
!
!    Surprisingly, the diameter of the sphere can be found by solving
!    a 3 by 3 linear system.  This is because the vectors P2 - P1,
!    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
!    right triangle with the diameter.  Hence, the dot product of
!    P2 - P1 with the diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
!    the diameter vector originating at P1.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4,
!    the coordinates of the vertices.
!
!    Output, real R, XC, YC, ZC, the coordinates of the center of the
!    circumscribed sphere, and its radius.  If the linear system is
!    singular, then R = -1, XC = YC = ZC = 0.
!
  implicit none
!
  integer, parameter :: n = 3
  integer, parameter :: nrhs = 1
!
  real a(n,n+nrhs)
  real enormsq0_3d
  integer info
  real r
  real x1
  real x2
  real x3
  real x4
  real xc
  real y1
  real y2
  real y3
  real y4
  real yc
  real z1
  real z2
  real z3
  real z4
  real zc
!
!  Set up the linear system.
!
  a(1,1) = x2 - x1
  a(1,2) = y2 - y1
  a(1,3) = z2 - z1
  a(1,4) = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )

  a(2,1) = x3 - x1
  a(2,2) = y3 - y1
  a(2,3) = z3 - z1
  a(2,4) = enormsq0_3d ( x1, y1, z1, x3, y3, z3 )

  a(3,1) = x4 - x1
  a(3,2) = y4 - y1
  a(3,3) = z4 - z1
  a(3,4) = enormsq0_3d ( x1, y1, z1, x4, y4, z4 )
!
!  Solve the linear system.
!
  call rmat_solve ( a, n, nrhs, info )
!
!  Compute R, X, Y.
!
  if ( info /= 0 ) then
    r = -1.0E+00
    xc = 0.0E+00
    yc = 0.0E+00
    zc = 0.0E+00
  else
    r = 0.5E+00 * sqrt ( a(1,N+1)**2 + a(2,N+1)**2 + a(3,N+1)**2 )
    xc = x1 + 0.5E+00 * a(1,N+1)
    yc = y1 + 0.5E+00 * a(2,N+1)
    zc = z1 + 0.5E+00 * a(3,N+1)
  end if

  return
end
subroutine tetra_contains_point_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, x, y, z, inside )
!
!*******************************************************************************
!
!! TETRA_CONTAINS_POINT_3D finds if a point is inside a tetrahedron in 3D.
!
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, Z4, Y4, Z4, the
!    coordinates of the tetrahedron vertices.
!
!    Input, real X, Y, Z, the point to be checked.
!
!    Output, logical INSIDE, is .TRUE. if (X,Y,Z) is inside
!    the tetrahedron or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  real c(4)
  logical inside
  real x
  real x1
  real x2
  real x3
  real x4
  real y
  real y1
  real y2
  real y3
  real y4
  real z
  real z1
  real z2
  real z3
  real z4
!
  call tetra_barycentric_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
    x4, y4, z4, x, y, z, c )
!
!  If the point is in the tetrahedron, its barycentric coordinates
!  must be nonnegative.
!
  if ( any ( c(1:4) < 0.0E+00 ) ) then
    inside = .false.
  else
    inside = .true.
  end if

  return
end
subroutine tetra_shape_3d ( max_num, max_order, point_num, face_num, &
   face_order, point_coord, face_point )
!
!*******************************************************************************
!
!! TETRA_SHAPE_3D describes a tetrahedron in 3D.
!
!
!  Discussion:
!
!    The vertices of the tetrahedron lie on the unit sphere.
!    The last point is the north pole of the unit sphere.
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAX_NUM, the maximum number of faces, and of
!    points.  This quantity is used to dimension arrays.
!    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
!    not be set.
!
!    Input, integer MAX_ORDER, the maximum order of any face.
!    This quantity is used to dimension arrays.  If FACE_ORDER
!    exceeds MAX_ORDER, the arrays will not be set.
!
!    Output, integer POINT_NUM, the number of points in the shape.
!
!    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
!    the X, Y and Z coordinates of point J.
!
!    Output, integer FACE_NUM, the number of faces in the shape.
!
!    Output, integer FACE_ORDER, the number of vertices per face.
!
!    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
!    contains the index of the I-th point in the J-th face.  The
!    points are listed in the counter-clockwise direction defined
!    by the outward normal at the face.
!
  implicit none
!
  integer max_num
  integer max_order
!
  integer face_num
  integer face_order
  integer face_point(max_order,max_num)
  integer point_num
  real point_coord(3,max_num)
!
  point_num = 4
  face_num = 4
  face_order = 3
!
!  Check.
!
  if ( point_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRA_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of vertices exceeds MAX_NUM.'
    return
  end if

  if ( face_num > max_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRA_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Number of faces exceeds MAX_NUM.'
    return
  end if

  if ( face_order > max_order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRA_SHAPE_3D - Fatal error!'
    write ( *, '(a)' ) '  Face order exceeds MAX_ORDER.'
    return
  end if
!
!  Set point coordinates.
!
  point_coord(1,1) =   0.0E+00
  point_coord(2,1) =   2.0E+00 * sqrt ( 2.0E+00 ) / 3.0E+00
  point_coord(3,1) = - 1.0E+00 / 3.0E+00

  point_coord(1,2) = - sqrt ( 2.0E+00 / 3.0E+00 )
  point_coord(2,2) = - sqrt ( 2.0E+00 ) / 3.0E+00
  point_coord(3,2) = - 1.0E+00 / 3.0E+00

  point_coord(1,3) =   sqrt ( 2.0E+00 / 3.0E+00 )
  point_coord(2,3) = - sqrt ( 2.0E+00 ) / 3.0E+00
  point_coord(3,3) = - 1.0E+00 / 3.0E+00

  point_coord(1,4) = 0.0E+00
  point_coord(2,4) = 0.0E+00
  point_coord(3,4) = 1.0E+00
!
!  Set faces.
!
  face_point(1,1) = 1
  face_point(2,1) = 3
  face_point(3,1) = 2

  face_point(1,2) = 1
  face_point(2,2) = 2
  face_point(3,2) = 4

  face_point(1,3) = 1
  face_point(2,3) = 4
  face_point(3,3) = 3

  face_point(1,4) = 2
  face_point(2,4) = 3
  face_point(3,4) = 4

  return
end
subroutine tetra_volume_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x4, y4, z4, volume )
!
!*******************************************************************************
!
!! TETRA_VOLUME_3D computes the volume of a tetrahedron in 3D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
!    coordinates of the corners of the tetrahedron.
!
!    Output, real VOLUME, the volume of the tetrahedron.
!
  implicit none
!
  real a(4,4)
  real rmat4_det
  real volume
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
  real z1
  real z2
  real z3
  real z4
!
  a(1:4,1) = (/ x1, x2, x3, x4 /)

  a(1:4,2) = (/ y1, y2, y3, y4 /)

  a(1:4,3) = (/ z1, z2, z3, z4 /)

  a(1:4,4) = (/ 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00 /)

  volume = abs ( rmat4_det ( a ) ) / 6.0E+00

  return
end
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone
!
  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tmat_init ( a )
!
!*******************************************************************************
!
!! TMAT_INIT initializes the geometric transformation matrix.
!
!
!  Definition:
!
!    The geometric transformation matrix can be thought of as a 4 by 4
!    matrix "A" having components:
!
!      r11 r12 r13 t1
!      r21 r22 r23 t2
!      r31 r32 r33 t3
!        0   0   0  1
!
!    This matrix encodes the rotations, scalings and translations that
!    are applied to graphical objects.
!
!    A point P = (x,y,z) is rewritten in "homogeneous coordinates" as
!    PH = (x,y,z,1).  Then to apply the transformations encoded in A to
!    the point P, we simply compute A * PH.
!
!    Individual transformations, such as a scaling, can be represented
!    by simple versions of the transformation matrix.  If the matrix
!    A represents the current set of transformations, and we wish to
!    apply a new transformation B, then the original points are
!    transformed twice:  B * ( A * PH ).  The new transformation B can
!    be combined with the original one A, to give a single matrix C that
!    encodes both transformations: C = B * A.
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the geometric transformation matrix.
!
  implicit none
!
  real a(4,4)
  integer i
  integer j
!
  do i = 1, 4
    do j = 1, 4
      if ( i == j ) then
        a(i,j) = 1.0E+00
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine tmat_mxm ( a, b, c )
!
!*******************************************************************************
!
!! TMAT_MXM multiplies two geometric transformation matrices.
!
!
!  Note:
!
!    The product is accumulated in a temporary array, and then assigned
!    to the result.  Therefore, it is legal for any two, or all three,
!    of the arguments to share memory.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the first geometric transformation matrix.
!
!    Input, real B(4,4), the second geometric transformation matrix.
!
!    Output, real C(4,4), the product A * B.
!
  implicit none
!
  real a(4,4)
  real b(4,4)
  real c(4,4)
!
  c(1:4,1:4) = matmul ( a(1:4,1:4), b(1:4,1:4) )

  return
end
subroutine tmat_mxp ( a, x, y )
!
!*******************************************************************************
!
!! TMAT_MXP multiplies a geometric transformation matrix times a point.
!
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the geometric transformation matrix.
!
!    Input, real X(3), the point to be multiplied.  The fourth component
!    of X is implicitly assigned the value of 1.
!
!    Output, real Y(3), the result of A*X.  The product is accumulated in
!    a temporary vector, and then assigned to the result.  Therefore, it
!    is legal for X and Y to share memory.
!
  implicit none
!
  real a(4,4)
  real x(3)
  real y(3)
!
  y(1:3) = a(1:3,4) + matmul ( a(1:3,1:3), x(1:3) )

  return
end
subroutine tmat_mxp2 ( a, x, y, n )
!
!*******************************************************************************
!
!! TMAT_MXP2 multiplies a geometric transformation matrix times N points.
!
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the geometric transformation matrix.
!
!    Input, real X(3,N), the points to be multiplied.
!
!    Output, real Y(3,N), the transformed points.  Each product is
!    accumulated in a temporary vector, and then assigned to the
!    result.  Therefore, it is legal for X and Y to share memory.
!
  implicit none
!
  integer n
!
  real a(4,4)
  integer i
  integer j
  integer k
  real x(3,n)
  real y(3,n)
  real z(3)
!
  do k = 1, n

    do i = 1, 3
      z(i) = a(i,4)
      do j = 1, 3
        z(i) = z(i) + a(i,j) * x(j,k)
      end do
    end do

    y(1:3,k) = z(1:3)

  end do

  return
end
subroutine tmat_mxv ( a, x, y )
!
!*******************************************************************************
!
!! TMAT_MXV multiplies a geometric transformation matrix times a vector.
!
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the geometric transformation matrix.
!
!    Input, real X(3), the vector to be multiplied.  The fourth component
!    of X is implicitly assigned the value of 1.
!
!    Output, real Y(3), the result of A*X.  The product is accumulated in
!    a temporary vector, and then assigned to the result.  Therefore, it
!    is legal for X and Y to share memory.
!
  implicit none
!
  real a(4,4)
  integer i
  integer j
  real x(3)
  real y(3)
  real z(3)
!
  do i = 1, 3
    z(i) = 0.0E+00
    do j = 1, 3
      z(i) = z(i) + a(i,j) * x(j)
    end do
    z(i) = z(i) + a(i,4)
  end do

  y(1:3) = z(1:3)

  return
end
subroutine tmat_rot_axis ( a, b, angle, axis )
!
!*******************************************************************************
!
!! TMAT_ROT_AXIS applies a coordinate axis rotation to the geometric transformation matrix.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the current geometric transformation matrix.
!
!    Output, real B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, real ANGLE, the angle, in degrees, of the rotation.
!
!    Input, character AXIS, is 'X', 'Y' or 'Z', specifying the coordinate
!    axis about which the rotation occurs.
!
  implicit none
!
  real a(4,4)
  real angle
  real angle_rad
  character axis
  real b(4,4)
  real c(4,4)
  real degrees_to_radians
  integer i
  integer j
!
  angle_rad = degrees_to_radians ( angle )

  call tmat_init ( c )

  if ( axis == 'X' .or. axis == 'x' ) then
    c(2,2) =   cos ( angle_rad )
    c(2,3) = - sin ( angle_rad )
    c(3,2) =   sin ( angle_rad )
    c(3,3) =   cos ( angle_rad )
  else if ( axis == 'Y' .or. axis == 'y' ) then
    c(1,1) =   cos ( angle_rad )
    c(1,3) =   sin ( angle_rad )
    c(3,1) = - sin ( angle_rad )
    c(3,3) =   cos ( angle_rad )
  else if ( axis == 'Z' .or. axis == 'z' ) then
    c(1,1) =   cos ( angle_rad )
    c(1,2) = - sin ( angle_rad )
    c(2,1) =   sin ( angle_rad )
    c(2,2) =   cos ( angle_rad )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TMAT_ROT_AXIS - Fatal error!'
    write ( *, '(a)' ) '  Illegal rotation axis: ' // axis
    write ( *, '(a)' ) '  Legal choices are ''X'', ''Y'', or ''Z''.'
    return
  end if

  b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )

  return
end
subroutine tmat_rot_vector ( a, b, angle, axis )
!
!*******************************************************************************
!
!! TMAT_ROT_VECTOR applies an arbitrary axis rotation to the geometric transformation matrix.
!
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the current geometric transformation matrix.
!
!    Output, real B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, real ANGLE, the angle, in degrees, of the rotation.
!
!    Input, real AXIS(3), the axis vector about which rotation occurs.
!    AXIS may not be the zero vector.
!
  implicit none
!
  real a(4,4)
  real angle
  real angle_rad
  real axis(3)
  real b(4,4)
  real c(4,4)
  real ca
  real degrees_to_radians
  integer i
  integer j
  real norm
  real sa
  real v1
  real v2
  real v3
!
  v1 = axis(1)
  v2 = axis(2)
  v3 = axis(3)

  norm = sqrt ( v1 * v1 + v2 * v2 + v3 * v3 )

  if ( norm == 0.0E+00 ) then
    return
  end if

  v1 = v1 / norm
  v2 = v2 / norm
  v3 = v3 / norm

  angle_rad = degrees_to_radians ( angle )
  ca = cos ( angle_rad )
  sa = sin ( angle_rad )

  call tmat_init ( c )

  c(1,1) =                    v1 * v1 + ca * ( 1.0E+00 - v1 * v1 )
  c(1,2) = ( 1.0E+00 - ca ) * v1 * v2 - sa * v3
  c(1,3) = ( 1.0E+00 - ca ) * v1 * v3 + sa * v2

  c(2,1) = ( 1.0E+00 - ca ) * v2 * v1 + sa * v3
  c(2,2) =                    v2 * v2 + ca * ( 1.0E+00 - v2 * v2 )
  c(2,3) = ( 1.0E+00 - ca ) * v2 * v3 - sa * v1

  c(3,1) = ( 1.0E+00 - ca ) * v3 * v1 - sa * v2
  c(3,2) = ( 1.0E+00 - ca ) * v3 * v2 + sa * v1
  c(3,3) =                    v3 * v3 + ca * ( 1.0E+00 - v3 * v3 )

  b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )

  return
end
subroutine tmat_scale ( a, b, sx, sy, sz )
!
!*******************************************************************************
!
!! TMAT_SCALE applies a scaling to the geometric transformation matrix.
!
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the current geometric transformation matrix.
!
!    Output, real B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, real SX, SY, SZ, the scalings to be applied to the X, Y and
!    Z coordinates.
!
  implicit none
!
  real a(4,4)
  real b(4,4)
  real c(4,4)
  integer i
  integer j
  real sx
  real sy
  real sz
!
  call tmat_init ( c )

  c(1,1) = sx
  c(2,2) = sy
  c(3,3) = sz

  b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )

  return
end
subroutine tmat_shear ( a, b, axis, s )
!
!*******************************************************************************
!
!! TMAT_SHEAR applies a shear to the geometric transformation matrix.
!
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the current geometric transformation matrix.
!
!    Output, real B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, character ( len = 2 ) AXIS, is 'XY', 'XZ', 'YX', 'YZ', 'ZX' or 'ZY',
!    specifying the shear equation:
!
!      XY:  x' = x + s * y;
!      XZ:  x' = x + s * z;
!      YX:  y' = y + s * x;
!      YZ:  y' = y + s * z;
!      ZX:  z' = z + s * x;
!      ZY:  z' = z + s * y.
!
!    Input, real S, the shear coefficient.
!
  implicit none
!
  real a(4,4)
  character ( len = 2 ) axis
  real b(4,4)
  real c(4,4)
  integer i
  integer j
  real s
!
  call tmat_init ( c )

  if ( axis == 'XY' .or. axis == 'xy' ) then
    c(1,2) = s
  else if ( axis == 'XZ' .or. axis == 'xz' ) then
    c(1,3) = s
  else if ( axis == 'YX' .or. axis == 'yx' ) then
    c(2,1) = s
  else if ( axis == 'YZ' .or. axis == 'yz' ) then
    c(2,3) = s
  else if ( axis == 'ZX' .or. axis == 'zx' ) then
    c(3,1) = s
  else if ( axis == 'ZY' .or. axis == 'zy' ) then
    c(3,2) = s
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TMAT_SHEAR - Fatal error!'
    write ( *, '(a)' ) '  Illegal shear axis: ' // axis
    write ( *, '(a)' ) '  Legal choices are XY, XZ, YX, YZ, ZX, or ZY.'
    return
  end if

  b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )

  return
end
subroutine tmat_trans ( a, b, x, y, z )
!
!*******************************************************************************
!
!! TMAT_TRANS applies a translation to the geometric transformation matrix.
!
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real A(4,4), the current geometric transformation matrix.
!
!    Output, real B(4,4), the modified transformation matrix.
!    A and B may share the same memory.
!
!    Input, real X, Y, Z, the translation.  This may be thought of as the
!    point that the origin moves to under the translation.
!
  implicit none
!
  real a(4,4)
  real b(4,4)
  real x
  real y
  real z
!
  b(1:4,1:4) = a(1:4,1:4)

  b(1:3,4) = b(1:3,4) + (/ x, y, z /)

  return
end
function torus_area_3d ( r1, r2 )
!
!*******************************************************************************
!
!! TORUS_AREA_3D returns the area of a torus in 3D.
!
!
!  Integration region:
!
!    Points (X,Y,Z) such that:
!
!      ( SQRT ( X**2 + Y**2 ) - R1 )**2 + Z**2 = R2**2.
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R1, R2, the two radii that define the torus.
!
!    Output, real TORUS_AREA_3D, the area of the torus.
!
  implicit none
!
  real r_pi
  real r1
  real r2
  real torus_area_3d
!
  torus_area_3d = 4.0E+00 * r_pi ( )**2 * r1 * r2

  return
end
subroutine torus_volume_3d ( r1, r2, volume )
!
!*******************************************************************************
!
!! TORUS_VOLUME_3D computes the volume of a torus in 3D.
!
!
!  Definition:
!
!    A torus with radii R1 and R2 is the set of points (X,Y,Z)
!    satisfying:
!
!    ( sqrt ( X**2 + Y**2 ) - R1 )**2 + Z**2 <= R2**2
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R1, R2, the "inner" and "outer" radii of the torus.
!
!    Output, real VOLUME, the volume of the torus.
!
  implicit none
!
  real r_pi
  real r1
  real r2
  real volume
!
  volume = 2.0E+00 * ( r_pi ( ) )**2 * r1 * r2**2

  return
end
subroutine triangle_angles_2d ( x1, y1, x2, y2, x3, y3, angle1, angle2, angle3 )
!
!*******************************************************************************
!
!! TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
!
!
!  Discussion:
!
!    The law of cosines is used:
!
!      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
!
!    where GAMMA is the angle opposite side C.
!
!  Modified:
!
!    15 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real ANGLE1, ANGLE2, ANGLE3, the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none
!
  real a
  real angle1
  real angle2
  real angle3
  real arc_cosine
  real b
  real c
  real enorm0_2d
  real r_pi
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
!
  a = enorm0_2d ( x2, y2, x1, y1 )
  b = enorm0_2d ( x3, y3, x2, y2 )
  c = enorm0_2d ( x1, y1, x3, y3 )
!
!  Take care of a ridiculous special case.
!
  if ( a == 0.0E+00 .and. b == 0.0E+00 .and. c == 0.0E+00 ) then
    angle1 = 2.0E+00 * r_pi ( ) / 3.0E+00
    angle2 = 2.0E+00 * r_pi ( ) / 3.0E+00
    angle3 = 2.0E+00 * r_pi ( ) / 3.0E+00
    return
  end if

  if ( c == 0.0E+00 .or. a == 0.0E+00 ) then
    angle1= r_pi ( )
  else
    angle1 = arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0E+00 * c * a ) )
  end if

  if ( a == 0.0E+00 .or. b == 0.0E+00 ) then
    angle2 = r_pi ( )
  else
    angle2 = arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0E+00 * a * b ) )
  end if

  if ( b == 0.0E+00 .or. c == 0.0E+00 ) then
    angle3 = r_pi ( )
  else
    angle3 = arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0E+00 * b * c ) )
  end if

  return
end
subroutine triangle_angles_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, angle1, &
  angle2, angle3 )
!
!*******************************************************************************
!
!! TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
!
!
!  Discussion:
!
!    The law of cosines is used:
!
!      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
!
!    where GAMMA is the angle opposite side C.
!
!  Modified:
!
!    15 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the triangle vertices.
!
!    Output, real ANGLE1, ANGLE2, ANGLE3, the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none
!
  real a
  real angle1
  real angle2
  real angle3
  real arc_cosine
  real b
  real c
  real enorm0_3d
  real r_pi
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
  a = enorm0_3d ( x2, y2, z2, x1, y1, z1 )
  b = enorm0_3d ( x3, y3, z3, x2, y2, z2 )
  c = enorm0_3d ( x1, y1, z1, x3, y3, z3 )
!
!  Take care of a ridiculous special case.
!
  if ( a == 0.0E+00 .and. b == 0.0E+00 .and. c == 0.0E+00 ) then
    angle1 = 2.0E+00 * r_pi ( ) / 3.0E+00
    angle2 = 2.0E+00 * r_pi ( ) / 3.0E+00
    angle3 = 2.0E+00 * r_pi ( ) / 3.0E+00
    return
  end if

  if ( c == 0.0E+00 .or. a == 0.0E+00 ) then
    angle1= r_pi ( )
  else
    angle1 = arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0E+00 * c * a ) )
  end if

  if ( a == 0.0E+00 .or. b == 0.0E+00 ) then
    angle2 = r_pi ( )
  else
    angle2 = arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0E+00 * a * b ) )
  end if

  if ( b == 0.0E+00 .or. c == 0.0E+00 ) then
    angle3 = r_pi ( )
  else
    angle3 = arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0E+00 * b * c ) )
  end if

  return
end
subroutine triangle_area_2_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, area )
!
!*******************************************************************************
!
!! TRIANGLE_AREA_2_3D computes the area of a triangle in 3D.
!
!
!  Discussion:
!
!    This routine computes the area "the hard way".
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the triangle vertices.
!
!    Output, real AREA, the area of the triangle.
!
  implicit none
!
  real alpha
  real area
  real base
  real dot
  real enorm0_3d
  real height
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
  real z1
  real z2
  real z3
!
!  Find the projection of (P3-P1) onto (P2-P1).
!
  dot = ( x2 - x1 ) * ( x3 - x1 ) + ( y2 - y1 ) * ( y3 - y1 ) + &
        ( z2 - z1 ) * ( z3 - z1 )

  base = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
!
!  The height of the triangle is the length of (P3-P1) after its
!  projection onto (P2-P1) has been subtracted.
!
  if ( base == 0.0E+00 ) then

    height = 0.0E+00

  else

    alpha = dot / ( base * base )

    height = enorm0_3d ( x1 + alpha * ( x2 - x1 ), y1 + alpha * ( y2 - y1 ), &
                         z1 + alpha * ( z2 - z1 ), x3, y3, z3 )

  end if

  area = 0.5E+00 * base * height

  return
end
subroutine triangle_area_2d ( x1, y1, x2, y2, x3, y3, area )
!
!*******************************************************************************
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!
!  Modified:
!
!    09 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real AREA, the absolute area of the triangle.
!
  implicit none
!
  real area
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
!
  area = 0.5E+00 * abs ( &
    ( x1 * ( y2 - y3 ) + x2 * ( y3 - y1 ) + x3 * ( y1 - y2 ) ) )

  return
end
subroutine triangle_area_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, area )
!
!*******************************************************************************
!
!! TRIANGLE_AREA_3D computes the area of a triangle in 3D.
!
!
!  Discussion:
!
!    This routine uses the fact that the norm of the cross product vector
!    is the area of the parallelogram they form.  The triangle they
!    form has half that area.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    06 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the triangle vertices.
!
!    Output, real AREA, the area of the triangle.
!
  implicit none
!
  real area
  real enorm_3d
  real norm
  real x1
  real x2
  real x3
  real x4
  real y1
  real y2
  real y3
  real y4
  real z1
  real z2
  real z3
  real z4
!
  call cross_3d ( x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1, z3 - z1, &
                  x4,      y4,      z4 )

  norm = enorm_3d ( x4, y4, z4 )

  area = 0.5E+00 * norm

  return
end
subroutine triangle_area_heron ( s1, s2, s3, area )
!
!*******************************************************************************
!
!! TRIANGLE_AREA_HERON computes the area of a triangle using Heron's formula.
!
!
!  Discussion:
!
!    The formula is valid for any spatial dimension, depending only
!    on the lengths of the sides, and not the coordinates of the vertices.
!
!  Modified:
!
!    13 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real S1, S2, S3, the lengths of the three sides.
!
!    Output, real AREA, the area of the triangle, or -1.0 if the
!    sides cannot constitute a triangle.
!
  implicit none
!
  real area
  real s1
  real s2
  real s3
!
  area = (   s1 + s2 + s3 ) &
       * ( - s1 + s2 + s3 ) &
       * (   s1 - s2 + s3 ) &
       * (   s1 + s2 - s3 )

  if ( area < 0.0E+00 ) then
    area = -1.0E+00
    return
  end if

  area = 0.25E+00 * sqrt ( area )

  return
end
subroutine triangle_area_signed_2d ( x1, y1, x2, y2, x3, y3, area_signed )
!
!*******************************************************************************
!
!! TRIANGLE_AREA_SIGNED_2D computes the signed area of a triangle in 2D.
!
!
!  Discussion:
!
!    The area is POSITIVE if the vertices are listed in counterclockwise
!    order.
!
!  Modified:
!
!    09 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real AREA_SIGNED, the signed area of the triangle.  The size
!    of AREA_SIGNED is the usual area, but the sign can be negative
!    or positive, depending on the order in which the vertices were listed.
!
  implicit none
!
  real area_signed
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
!
  area_signed = 0.5E+00 * &
    ( x1 * ( y2 - y3 ) + x2 * ( y3 - y1 ) + x3 * ( y1 - y2 ) )

  return
end
subroutine triangle_barycentric_2d ( x1, y1, x2, y2, x3, y3, x, y, c )
!
!*******************************************************************************
!
!! TRIANGLE_BARYCENTRIC_2D finds the barycentric coordinates of a point in 2D.
!
!
!  Discussion:
!
!    The barycentric coordinate of point X related to vertex A can be
!    interpreted as the ratio of the area of the triangle with 
!    vertex A replaced by vertex X to the area of the original 
!    triangle.
!
!  Modified:
!
!    20 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, real C(3), the barycentric coordinates of (X,Y) with respect
!    to the triangle.
!
  implicit none
!
  integer, parameter :: n = 2
  integer, parameter :: nrhs = 1
!
  real a(n,n+nrhs)
  real c(3)
  integer info
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1 ) C1  = X-X1
!    ( Y2-Y1  Y3-Y1 ) C2    Y-Y1
!
!  which is satisfied by the barycentric coordinates of (X,Y).
!
  a(1,1) = x2 - x1
  a(1,2) = x3 - x1
  a(1,3) = x - x1

  a(2,1) = y2 - y1
  a(2,2) = y3 - y1
  a(2,3) = y - y1
!
!  Solve the linear system.
!
  call rmat_solve ( a, n, nrhs, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_BARYCENTRIC_2D - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper triangle.'
    stop
  end if

  c(1) = a(1,3)
  c(2) = a(2,3)
  c(3) = 1.0E+00 - c(1) - c(2)

  return
end
subroutine triangle_centroid_2d ( x1, y1, x2, y2, x3, y3, x, y )
!
!*******************************************************************************
!
!! TRIANGLE_CENTROID_2D computes the centroid of a triangle in 2D.
!
!
!  Discussion:
!
!    The centroid of a triangle can also be considered the center
!    of gravity, assuming that the triangle is made of a thin uniform
!    sheet of massy material.
!
!    The centroid of a triangle is the intersection of the medians.
!
!    A median of a triangle is a line connecting a vertex to the
!    midpoint of the opposite side.
!
!    In barycentric coordinates, in which the vertices of the triangle
!    have the coordinates (1,0,0), (0,1,0) and (0,0,1), the centroid
!    has coordinates (1/3,1/3,1/3).
!
!    In geometry, the centroid of a triangle is often symbolized by "G".
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    29 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real X, Y, the coordinates of the centroid of the triangle.
!
  implicit none
!
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  x = ( x1 + x2 + x3 ) / 3.0E+00
  y = ( y1 + y2 + y3 ) / 3.0E+00

  return
end
subroutine triangle_centroid_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z )
!
!*******************************************************************************
!
!! TRIANGLE_CENTROID_3D computes the centroid of a triangle in 3D.
!
!
!  Discussion:
!
!    The centroid of a triangle can also be considered the center
!    of gravity, assuming that the triangle is made of a thin uniform
!    sheet of massy material.
!
!    The centroid of a triangle is the intersection of the medians.
!    A median of a triangle is a line connecting any vertex to the
!    midpoint of the opposite side.
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the triangle vertices.
!
!    Output, real X, Y, Z, the coordinates of the centroid.
!
  implicit none
!
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
  real z
  real z1
  real z2
  real z3
!
  x = ( x1 + x2 + x3 ) / 3.0E+00
  y = ( y1 + y2 + y3 ) / 3.0E+00
  z = ( z1 + z2 + z3 ) / 3.0E+00

  return
end
subroutine triangle_circumcenter_2d ( x1, y1, x2, y2, x3, y3, x, y )
!
!*******************************************************************************
!
!! TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
!
!
!  Discussion:
!
!    The circumcenter of a triangle is the center of the circumcircle, the
!    circle that passes through the three vertices of the triangle.
!
!    The circumcircle contains the triangle, but it is not necessarily the
!    smallest triangle to do so.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the circle.
!
!    The circumcenter is the intersection of the perpendicular bisectors
!    of the sides of the triangle.
!
!    Surprisingly, the diameter of the circle can be found by solving
!    a 2 by 2 linear system.  If we label the vertices of the triangle
!    P1, P2 and P3, then the vectors P2 - P1 and P3 - P1 are secants of
!    the circle.
!
!    In geometry, the circumcenter of a triangle is often symbolized by "O".
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    02 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real X, Y, the coordinates of the circumcenter of the triangle.
!
  implicit none
!
  integer, parameter :: n = 2
  integer, parameter :: nrhs = 1
!
  real a(n,n+nrhs)
  real enormsq0_2d
  integer info
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
!  Set up the linear system.
!
  a(1,1) = x2 - x1
  a(1,2) = y2 - y1
  a(1,3) = enormsq0_2d ( x1, y1, x2, y2 )

  a(2,1) = x3 - x1
  a(2,2) = y3 - y1
  a(2,3) = enormsq0_2d ( x1, y1, x3, y3 )
!
!  Solve the linear system.
!
  call rmat_solve ( a, n, nrhs, info )
!
!  Compute X, Y.
!
  if ( info /= 0 ) then
    x = 0.0E+00
    y = 0.0E+00
  else
    x = x1 + 0.5E+00 * a(1,n+1)
    y = y1 + 0.5E+00 * a(2,n+1)
  end if

  return
end
subroutine triangle_circumcircle_2d ( x1, y1, x2, y2, x3, y3, r, x, y )
!
!*******************************************************************************
!
!! TRIANGLE_CIRCUMCIRCLE_2D computes the circumcircle of a triangle in 2D.
!
!
!  Discussion:
!
!    The circumscribed circle of a triangle is the circle that passes through
!    the three vertices of the triangle.  The circumscribed circle contains
!    the triangle, but it is not necessarily the smallest triangle to do so.
!
!    Surprisingly, the diameter of the circle can be found by solving
!    a 2 by 2 linear system.  This is because the vectors P2 - P1
!    and P3 - P1 are secants of the circle, and each forms a right
!    triangle with the diameter.  Hence, the dot product of
!    P2 - P1 with the diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1.  This determines the
!    diameter vector originating at P1.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the circle.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    13 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real R, X, Y, the radius and coordinates of the center of the
!    circumscribed circle.  If the linear system is
!    singular, then R = -1, X = Y = 0.
!
  implicit none
!
  integer, parameter :: n = 2
  integer, parameter :: nrhs = 1
!
  real a(n,n+nrhs)
  real enormsq0_2d
  integer info
  real r
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
!  Set up the linear system.
!
  a(1,1) = x2 - x1
  a(1,2) = y2 - y1
  a(1,3) = enormsq0_2d ( x1, y1, x2, y2 )

  a(2,1) = x3 - x1
  a(2,2) = y3 - y1
  a(2,3) = enormsq0_2d ( x1, y1, x3, y3 )
!
!  Solve the linear system.
!
  call rmat_solve ( a, n, nrhs, info )
!
!  Compute R, X, Y.
!
  if ( info /= 0 ) then
    r = -1.0E+00
    x = 0.0E+00
    y = 0.0E+00
  else
    r = 0.5E+00 * sqrt ( a(1,n+1) * a(1,n+1) + a(2,n+1) * a(2,n+1) )
    x = x1 + 0.5E+00 * a(1,n+1)
    y = y1 + 0.5E+00 * a(2,n+1)
  end if

  return
end
subroutine triangle_contains_point_1_2d ( x1, y1, x2, y2, x3, y3, x, y, &
  inside )
!
!*******************************************************************************
!
!! TRIANGLE_CONTAINS_POINT_1_2D finds if a point is inside a triangle in 2D.
!
!
!  Modified:
!
!    16 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, logical INSIDE, is .TRUE. if (X,Y) is inside
!    the triangle or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  real c(3)
  logical inside
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  call triangle_barycentric_2d ( x1, y1, x2, y2, x3, y3, x, y, c )

  if ( any ( c(1:3) < 0.0E+00 ) ) then
    inside = .false.
  else
    inside = .true.
  end if

  return
end
subroutine triangle_contains_point_2_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )
!
!*******************************************************************************
!
!! TRIANGLE_CONTAINS_POINT_2_2D finds if a point is inside a triangle in 2D.
!
!
!  Discussion:
!
!    The routine assumes that the vertices are given in counter-clockwise
!    order.
!
!    The routine determines if (X,Y) is "to the right of" each of the lines
!    that bound the triangle.  It does this by computing the cross product
!    of vectors from a vertex to its next vertex, and to (X,Y).
!
!  Modified:
!
!    09 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real X, Y, the point to be checked.
!
!    Output, logical INSIDE, is .TRUE. if (X,Y) is inside
!    the triangle or on its boundary, and .FALSE. otherwise.
!
  implicit none
!
  logical inside
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  inside = .false.

  if ( ( x - x1 ) * ( y2 - y1 ) - ( y - y1 ) * ( x2 - x1 ) > 0.0E+00 ) then
    return
  end if

  if ( ( x - x2 ) * ( y3 - y2 ) - ( y - y2 ) * ( x3 - x2 ) > 0.0E+00 ) then
    return
  end if

  if ( ( x - x3 ) * ( y1 - y3 ) - ( y - y3 ) * ( x1 - x3 ) > 0.0E+00 ) then
    return
  end if

  inside = .true.

  return
end
subroutine triangle_diameter_2d ( x1, y1, x2, y2, x3, y3, diam )
!
!*******************************************************************************
!
!! TRIANGLE_DIAMETER_2D computes the diameter of a triangle in 2D.
!
!
!  Discussion:
!
!    The diameter of a triangle is the diameter of the smallest circle
!    that can be drawn around the triangle.  At least two of the vertices
!    of the triangle will intersect the circle, but not necessarily
!    all three!
!
!  Modified:
!
!    13 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real DIAM, the diameter of the triangle.
!
  implicit none
!
  integer, parameter :: n = 2
  integer, parameter :: nrhs = 1
!
  real a
  real b
  real c
  real diam
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
!
!  Compute the (squares of) the lengths of the sides.
!
  a = ( x1 - x2 )**2 + ( y1 - y2 )**2
  b = ( x2 - x3 )**2 + ( y2 - y3 )**2
  c = ( x3 - x1 )**2 + ( y3 - y1 )**2
!
!  Take care of a zero side.
!
  if ( a == 0.0E+00 ) then
    diam = sqrt ( b )
    return
  else if ( b == 0.0E+00 ) then
    diam = sqrt ( c )
    return
  else if ( c == 0.0E+00 ) then
    diam = sqrt ( a )
    return
  end if
!
!  Make A the largest.
!
  if ( a < b ) then
    call r_swap ( a, b )
  end if

  if ( a < c ) then
    call r_swap ( a, c )
  end if
!
!  If A is very large...
!
  if ( a > b + c ) then
    diam = sqrt ( a )
  else
    a = sqrt ( a )
    b = sqrt ( b )
    c = sqrt ( c )
    diam = 2.0E+00 * a * b * c / sqrt ( ( a + b + c ) * ( - a + b + c ) &
      * ( a - b + c ) * ( a + b - c ) )
  end if

  return
end
subroutine triangle_gridpoints_2d ( x1, y1, x2, y2, x3, y3, nsub, maxgrid, &
  ngrid, x, y )
!
!*******************************************************************************
!
!! TRIANGLE_GRIDPOINTS_2D computes gridpoints within a triangle in 2D.
!
!
!  Discussion:
!
!    The gridpoints are computed by repeated halving of the triangle.
!    The 0-th set of grid points is the vertices themselves.
!    The first set of grid points is the midpoints of the sides.
!    These points can be used to draw 4 triangles that make up the original
!    triangle.  The second set of grid points is the side midpoints and centers
!    of these four triangles.
!
!     NSUB                        NGRID
!    -----                        -----
!        0      1                  =  1  (centroid)
!        1      1 + 2              =  3  (vertices)
!        2      1 + 2 + 3          =  6
!        3      1 + 2 + 3 + 4      = 10
!        4      1 + 2 + 3 + 4 + 5  = 15
!
!    NGRID is the sum of the integers from 1 to NSUB+1 or
!
!      NGRID = (NSUB+1) * (NSUB+2) / 2
!
!  Modified:
!
!    18 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Input, integer NSUB, the number of subdivisions.
!
!    Input, integer MAXGRID, the maximum number of grid points.
!
!    Output, integer NGRID, the number of grid points returned.
!
!    Output, real X(MAXGRID), Y(MAXGRID), coordinates of the grid points.
!
  implicit none
!
  integer maxgrid
!
  real a
  real b
  real c
  integer i
  integer j
  integer k
  integer ngrid
  integer nsub
  real x(maxgrid)
  real x1
  real x2
  real x3
  real y(maxgrid)
  real y1
  real y2
  real y3
!
  ngrid = 0
!
!  Special case, NSUB = 0.
!
  if ( nsub == 0 ) then
    if ( maxgrid >= 1 ) then
      ngrid = 1
      x(1) = ( x1 + x2 + x3 ) / 3.0E+00
      y(1) = ( y1 + y2 + y3 ) / 3.0E+00
    end if
    return
  end if

  do i = 0, nsub

    a = real ( i ) / real ( nsub )

    do j = 0, nsub - i

      b = real ( j ) / real ( nsub )
      k = nsub - i - j
      c = real ( k ) / real ( nsub )

      if ( ngrid < maxgrid ) then

        ngrid = ngrid + 1
        x(ngrid) = ( a * x1 + b * x2 + c * x3 )
        y(ngrid) = ( a * y1 + b * y2 + c * y3 )

      end if

    end do
  end do

  return
end
subroutine triangle_incenter_2d ( x1, y1, x2, y2, x3, y3, xc, yc )
!
!*******************************************************************************
!
!! TRIANGLE_INCENTER_2D computes the incenter of a triangle in 2D.
!
!
!  Discussion:
!
!    The incenter of a triangle is the center of the inscribed circle.
!
!    The inscribed circle of a triangle is the largest circle that can
!    be drawn inside the triangle.
!
!    The inscribed circle is tangent to all three sides of the triangle.
!
!    The angle bisectors of the triangle intersect at the center of the
!    inscribed circle.
!
!    In geometry, the incenter is often represented by "I".
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real R, XC, YC, the radius and coordinates of the center of the
!    inscribed circle.
!
  implicit none
!
  real enorm0_2d
  real perim
  real s12
  real s23
  real s31
  real x1
  real x2
  real x3
  real xc
  real y1
  real y2
  real y3
  real yc
!
  s12 = enorm0_2d ( x1, y1, x2, y2 )
  s23 = enorm0_2d ( x2, y2, x3, y3 )
  s31 = enorm0_2d ( x3, y3, x1, y1 )

  perim = s12 + s23 + s31

  if ( perim == 0.0E+00 ) then
    xc = x1
    yc = y1
  else
    xc = ( s23 * x1 + s31 * x2 + s12 * x3 ) / perim
    yc = ( s23 * y1 + s31 * y2 + s12 * y3 ) / perim
  end if

  return
end
subroutine triangle_incircle_2d ( x1, y1, x2, y2, x3, y3, r, xc, yc )
!
!*******************************************************************************
!
!! TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
!
!
!  Discussion:
!
!    The inscribed circle of a triangle is the largest circle that can
!    be drawn inside the triangle.  It is tangent to all three sides,
!    and the lines from its center to the vertices bisect the angles
!    made by each vertex.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    29 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real R, XC, YC, the radius and coordinates of the center of the
!    inscribed circle.
!
  implicit none
!
  real enorm0_2d
  real perim
  real r
  real s12
  real s23
  real s31
  real x1
  real x2
  real x3
  real xc
  real y1
  real y2
  real y3
  real yc
!
  s12 = enorm0_2d ( x1, y1, x2, y2 )
  s23 = enorm0_2d ( x2, y2, x3, y3 )
  s31 = enorm0_2d ( x3, y3, x1, y1 )
  perim = s12 + s23 + s31

  if ( perim == 0.0E+00 ) then
    xc = x1
    yc = y1
    r = 0.0E+00
    return
  end if

  xc = ( s23 * x1 + s31 * x2 + s12 * x3 ) / perim
  yc = ( s23 * y1 + s31 * y2 + s12 * y3 ) / perim

  r = 0.5E+00 * sqrt ( ( - s12 + s23 + s31 ) * ( + s12 - s23 + s31 ) &
                 * ( + s12 + s23 - s31 ) / perim )

  return
end
subroutine triangle_line_imp_int_2d ( a, b, c, trix, triy, nin, xint, yint )
!
!*******************************************************************************
!
!! TRIANGLE_LINE_IMP_INT_2D finds where an implicit line intersects a triangle in 2D.
!
!
!  Definition:
!
!    An implicit line is the set of points ( X, Y ) satisfying
!
!      A * X + B * Y + C = 0
!
!    where at least one of A and B is not zero.
!
!  Discussion:
!
!    If the line happens to be identical with one of the sides
!    of the triangle, NIN will be returned as 3.
!
!    If the intersection point is one of the nodes of the triangle,
!    it will be counted twice, which will result in a value of NIN
!    that is 2 or 3, depending if the line also intersects
!    the side opposite the node.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, C, determine the equation of the line:
!    A*X + B*Y + C = 0.
!
!    Input, real TRIX(3), TRIY(3), the triangle vertices.
!
!    Output, integer NIN, the number of points of intersection
!    of the line with the triangle.  NIN may be 0, 1, 2 or 3.
!
!    Output, real XINT(3), YINT(3).  XINT(I), YINT(I) is the
!    I-th intersection point, for I = 1 to NIN.
!
  implicit none
!
  real a
  real a1
  real b
  real b1
  real c
  real c1
  logical inside
  integer ival
  integer nin
  real trix(3)
  real triy(3)
  real x
  real xint(3)
  real y
  real yint(3)
!
  nin = 0
!
!  Get the implicit form of the line through vertices 1 and 2.
!
  call line_exp2imp_2d ( trix(1), triy(1), trix(2), triy(2), a1, b1, c1 )
!
!  Seek an intersection with the original line.
!
  call lines_imp_int_2d ( a, b, c, a1, b1, c1, ival, x, y )
!
!  Determine if the intersection is inside the triangle.
!
  if ( ival == 1 ) then

    call triangle_contains_point_1_2d ( trix(1), triy(1), &
      trix(2), triy(2), trix(3), triy(3), x, y, inside )

    if ( inside ) then
      nin = nin + 1
      xint(nin) = x
      yint(nin) = y
    end if

  end if
!
!  Get the implicit form of the line through vertices 2 and 3.
!
  call line_exp2imp_2d ( trix(2), triy(2), trix(3), triy(3), a1, b1, c1 )
!
!  Seek an intersection with the original line.
!
  call lines_imp_int_2d ( a, b, c, a1, b1, c1, ival, x, y )
!
!  Determine if the intersection is inside the triangle.
!
  if ( ival == 1 ) then

    call triangle_contains_point_1_2d ( trix(1), triy(1), &
      trix(2), triy(2), trix(3), triy(3), x, y, inside )

    if ( inside ) then
      nin = nin + 1
      xint(nin) = x
      yint(nin) = y
    end if

  end if
!
!  Get the implicit form of the line through vertices 3 and 1.
!
  call line_exp2imp_2d ( trix(3), triy(3), trix(1), triy(1), a1, b1, c1 )
!
!  Seek an intersection with the original line.
!
  call lines_imp_int_2d ( a, b, c, a1, b1, c1, ival, x, y )
!
!  Determine if the intersection is inside the triangle.
!
  if ( ival == 1 ) then

    call triangle_contains_point_1_2d ( trix(1), triy(1), &
      trix(2), triy(2), trix(3), triy(3), x, y, inside )

    if ( inside ) then
      nin = nin + 1
      xint(nin) = x
      yint(nin) = y
    end if

  end if

  return
end
function triangle_orientation_2d ( x1, y1, x2, y2, x3, y3 )
!
!*******************************************************************************
!
!! TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
!
!
!  Discussion:
!
!    Three distinct non-colinear points in the plane define a circle.
!    If the points are visited in the order (x1,y1), (x2,y2), and then
!    (x3,y3), this motion defines a clockwise or counterclockwise
!    rotation along the circle.
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, integer TRIANGLE_ORIENTATION_2D, reports if the three points lie
!    clockwise on the circle that passes through them.  The possible
!    return values are:
!
!    0, the points are distinct, noncolinear, and lie counterclockwise
!    on their circle.
!
!    1, the points are distinct, noncolinear, and lie clockwise
!    on their circle.
!
!    2, the points are distinct and colinear.
!
!    3, at least two of the points are identical.
!
  implicit none
!
  real det
  integer triangle_orientation_2d
  real x1
  real x2
  real x3
  real y1
  real y2
  real y3
!
  if ( ( x1 == x2 .and. y1 == y2 ) .or. ( x2 == x3 .and. y2 == y3 ) .or. &
       ( x3 == x1 .and. y3 == y1 ) ) then
    triangle_orientation_2d = 3
    return
  end if

  det = ( x1 - x3 ) * (y2 - y3 ) - ( x2 - x3 ) * ( y1 - y3 )

  if ( det == 0.0E+00 ) then
    triangle_orientation_2d = 2
  else if ( det < 0.0E+00 ) then
    triangle_orientation_2d = 1
  else if ( det > 0.0E+00 ) then
    triangle_orientation_2d = 0
  end if

  return
end
subroutine triangle_orthocenter_2d ( x1, y1, x2, y2, x3, y3, x, y )
!
!*******************************************************************************
!
!! TRIANGLE_ORTHOCENTER_2D computes the orthocenter of a triangle in 2D.
!
!
!  Discussion:
!
!    The orthocenter is defined as the intersection of the three altitudes
!    of a triangle.
!
!    An altitude of a triangle is the line through a vertex of the triangle
!    and perpendicular to the opposite side.
!
!    In geometry, the orthocenter of a triangle is often symbolized by "H".
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Modified:
!
!    30 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real X, Y, the coordinates of the orthocenter of the triangle.
!
  implicit none
!
  integer ival
  real x
  real x1
  real x2
  real x23
  real x3
  real x31
  real y
  real y1
  real y2
  real y23
  real y3
  real y31
!
!  Determine a point (X23,Y23) common to the line ((X2,Y2),(X3,Y3)) and
!  its perpendicular through (X1,Y1).
!
  call line_exp_perp_2d ( x2, y2, x3, y3, x1, y1, x23, y23 )
!
!  Determine a point (X31,Y31) common to the line ((X3,Y3),(X1,Y1)) and
!  its perpendicular through (X2,Y2).
!
  call line_exp_perp_2d ( x3, y3, x1, y1, x2, y2, x31, y31 )
!
!  Determine (X,Y), the intersection of the lines ((X1,Y1),(X23,Y23)) and
!  ((X2,Y2),(X31,Y32)).
!
  call lines_exp_int_2d ( x1, y1, x23, y23, x2, y2, x31, y31, ival, x, y )

  return
end
subroutine triangle_point_dist_2d ( x1, y1, x2, y2, x3, y3, x, y, dist )
!
!*******************************************************************************
!
!! TRIANGLE_POINT_DIST_2D: distance ( triangle, point ) in 2D.
!
!
!  Modified:
!
!    15 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Input, real X, Y, the point which is to be checked.
!
!    Output, real DIST, the distance from the point to the triangle.
!    DIST is zero if the point lies exactly on the triangle.
!
  implicit none
!
  real dist
  real dist2
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
  call line_seg_point_dist_2d ( x1, y1, x2, y2, x, y, dist2 )
  dist = dist2

  call line_seg_point_dist_2d ( x2, y2, x3, y3, x, y, dist2 )
  dist = min ( dist, dist2 )

  call line_seg_point_dist_2d ( x3, y3, x1, y1, x, y, dist2 )
  dist = min ( dist, dist2 )

  return
end
subroutine triangle_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  x, y, z, dist )
!
!*******************************************************************************
!
!! TRIANGLE_POINT_DIST_3D: distance ( triangle, point ) in 3D.
!
!
!  Modified:
!
!    15 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the triangle vertices.
!
!    Input, real X, Y, Z, the point which is to be checked.
!
!    Output, real DIST, the distance from the point to the triangle.
!    DIST is zero if the point lies exactly on the triangle.
!
  implicit none
!
  real dist
  real dist2
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
  real z
  real z1
  real z2
  real z3
!
!  Compute the distances from the point to each of the sides.
!
  call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, x, y, z, dist2 )
  dist = dist2

  call line_seg_point_dist_3d ( x2, y2, z2, x3, y3, z3, x, y, z, dist2 )
  dist = min ( dist, dist2 )

  call line_seg_point_dist_3d ( x3, y3, z3, x1, y1, z1, x, y, z, dist2 )
  dist = min ( dist, dist2 )

  return
end
subroutine triangle_point_dist_signed_2d ( x1, y1, x2, y2, x3, y3, x, y, &
  dist_signed )
!
!*******************************************************************************
!
!! TRIANGLE_POINT_DIST_SIGNED_2D: signed distance ( triangle, point ) in 2D.
!
!
!  Modified:
!
!    09 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!    These should be given in counter clockwise order.
!
!    Input, real X, Y, the point which is to be checked.
!
!    Output, real DIST_SIGNED, the signed distance from the point to the
!    triangle.  DIST_SIGNED is the maximum of the signed distances from the
!    point to each of the three lines that make up the triangle.
!
!    If DIST_SIGNED is:
!    0, the point is on the boundary of the triangle;
!    negative, the point is in the triangle;
!    positive, the point is outside the triangle.
!
  implicit none
!
  real dis12
  real dis23
  real dis31
  real dist_signed
  real x
  real x1
  real x2
  real x3
  real y
  real y1
  real y2
  real y3
!
!  Compute the signed line-distances to the point.
!
  call line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dis12 )
  call line_exp_point_dist_signed_2d ( x2, y2, x3, y3, x, y, dis23 )
  call line_exp_point_dist_signed_2d ( x3, y3, x1, y1, x, y, dis31 )
!
!  If the point is inside the triangle, all the line-distances are negative.
!  The largest (negative) line-distance has the smallest magnitude,
!  and is the signed triangle-distance.
!
  if ( dis12 <= 0.0E+00 .and. dis23 <= 0.0E+00 .and. dis31 <= 0.0E+00 ) then
    dist_signed = max ( dis12, dis23, dis31 )
!
!  If the point is outside the triangle, then we have to compute
!  the (positive) line-segment distances and take the minimum.
!
  else

    call line_seg_point_dist_2d ( x1, y1, x2, y2, x, y, dis12 )
    call line_seg_point_dist_2d ( x2, y2, x3, y3, x, y, dis23 )
    call line_seg_point_dist_2d ( x3, y3, x1, y1, x, y, dis31 )

    dist_signed = min ( dis12, dis23, dis31 )

  end if

  return
end
subroutine triangle_point_near_2d ( x1, y1, x2, y2, x3, y3, x, y, xn, yn, &
  dist )
!
!*******************************************************************************
!
!! TRIANGLE_POINT_NEAR_2D computes the nearest triangle point to a point in 2D.
!
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Input, real X, Y, the point whose nearest neighbor
!    on the line is to be determined.
!
!    Output, real XN, YN, the nearest point to (X,Y).
!
!    Output, real DIST, the distance from the point to the triangle.
!
  implicit none
!
  real dist
  real dist12
  real dist23
  real dist31
  real t
  real x
  real x1
  real x2
  real x3
  real xn
  real xn12
  real xn23
  real xn31
  real y
  real y1
  real y2
  real y3
  real yn
  real yn12
  real yn23
  real yn31
!
!  Find the distance to each of the line segments that make up the edges
!  of the triangle.
!
  call line_seg_point_near_2d ( x1, y1, x2, y2, x, y, xn12, yn12, dist12, t )

  call line_seg_point_near_2d ( x2, y2, x3, y3, x, y, xn23, yn23, dist23, t )

  call line_seg_point_near_2d ( x3, y3, x1, y1, x, y, xn31, yn31, dist31, t )

  if ( dist12 <= dist23 .and. dist12 <= dist31 ) then
    dist = dist12
    xn = xn12
    yn = yn12
  else if ( dist23 <= dist12 .and. dist23 <= dist31 ) then
    dist = dist23
    xn = xn23
    yn = yn23
  else
    dist = dist31
    xn = xn31
    yn = yn31
  end if

  return
end
subroutine triangle_sample_2d ( x1, y1, x2, y2, x3, y3, x, y )
!
!*******************************************************************************
!
!! TRIANGLE_SAMPLE_2D returns a random point in a triangle.
!
!
!  Modified:
!
!    13 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Output, real X, Y, a random point in the triangle.
!
  implicit none
!
  real alpha
  real beta
  real r
  real x
  real x1
  real x12
  real x13
  real x2
  real x3
  real y
  real y1
  real y12
  real y13
  real y2
  real y3
!
  call random_number ( harvest = r )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
  alpha = sqrt ( r )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
  x12 = alpha * x1 + ( 1.0E+00 - alpha ) * x2
  y12 = alpha * y1 + ( 1.0E+00 - alpha ) * y2

  x13 = alpha * x1 + ( 1.0E+00 - alpha ) * x3
  y13 = alpha * y1 + ( 1.0E+00 - alpha ) * y3
!
!  Now choose, uniformly at random, a point on the line L.
!
  call random_number ( harvest = beta )

  x = beta * x12 + ( 1.0E+00 - beta ) * x13
  y = beta * y12 + ( 1.0E+00 - beta ) * y13

  return
end
subroutine triangle_xsi_to_xy_2d ( x1, y1, x2, y2, x3, y3, xsi1, xsi2, xsi3, &
  x, y )
!
!*******************************************************************************
!
!! TRIANGLE_XSI_TO_XY_2D converts from barycentric to XY coordinates in 2D.
!
!
!  Modified:
!
!    13 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Input, real XSI1, XSI2, XSI3, the barycentric coordinates of a point.
!    XSI1 + XSI2 + XSI3 should equal 1, but this is not checked.
!
!    Output, real X, Y, the XY coordinates of the point.
!
  implicit none
!
  real x
  real x1
  real x2
  real x3
  real xsi1
  real xsi2
  real xsi3
  real y
  real y1
  real y2
  real y3
!
  x = xsi1 * x1 + xsi2 * x2 + xsi3 * x3
  y = xsi1 * y1 + xsi2 * y2 + xsi3 * y3

  return
end
subroutine triangle_xy_to_xsi_2d ( x1, y1, x2, y2, x3, y3, x, y, xsi1, xsi2, &
  xsi3 )
!
!*******************************************************************************
!
!! TRIANGLE_XY_TO_XSI_2D converts from XY to barycentric in 2D.
!
!
!  Modified:
!
!    13 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle vertices.
!
!    Input, real X, Y, the XY coordinates of a point.
!
!    Output, real XSI1, XSI2, XSI3, the barycentric coordinates of the point.
!    XSI1 + XSI2 + XSI3 should equal 1.
!
  implicit none
!
  real det
  real x
  real x1
  real x2
  real x3
  real xsi1
  real xsi2
  real xsi3
  real y
  real y1
  real y2
  real y3
!
  det = ( x1 - x3 ) * ( y2 - y3 ) - ( x2 - x3 ) * ( y1 - y3 )

  xsi1 = (   ( y2 - y3 ) * ( x - x3 ) - ( x2 - x3 ) * ( y - y3 ) ) / det
  xsi2 = ( - ( y1 - y3 ) * ( x - x3 ) + ( x1 - x3 ) * ( y - y3 ) ) / det

  xsi3 = 1.0E+00 - xsi1 - xsi2

  return
end
subroutine triangulation_boundary_count ( point_num, tri_num, bound_num )
!
!*******************************************************************************
!
!! TRIANGULATION_BOUNDARY_COUNT returns the number of boundary edges.
!
!
!  Discussion:
!
!    We assume we are given information about a legal, maximal triangulation
!    of a set of points in the plane.  In a maximal triangulation, no more
!    edges can be added without crossing an existing edge.
!
!    Given the number of points and triangles, we are going to apply
!    Euler's formula to determine the number of edges that lie on the
!    convex hull of the set of points.
!
!    The number of faces, including the infinite face, is TRI_NUM + 1.
!
!    Let BOUND_NUM denote the number of edges on the boundary of the convex
!    hull.  Each of the TRI_NUM triangles uses three edges.  Every edge
!    occurs in two different faces, so the number of edges must be
!    ( 3 * TRI_NUM + BOUND_NUM ) / 2.
!
!    The number of points used in the triangulation is POINT_NUM.
!
!    Euler's formula asserts that, for a simple connected figure in the
!    plane with no edge crossings, POINT_NUM points, EDGE_NUM edges and
!    FACE_NUM faces:
!
!      POINT_NUM - EDGE_NUM + FACE_NUM = 2
!
!    In our context, this becomes
!
!      POINT_NUM - ( 3 * TRI_NUM + BOUND_NUM ) / 2 + TRI_NUM + 1 = 2
!
!    or
!
!      BOUND_NUM = 2 * POINT_NUM - TRI_NUM - 2
!
!  Modified:
!
!    25 July 2001
!
!  Reference:
!
!    de Berg, Krevald, Overmars and Schwarzkopf,
!    Computational Geometry, Section 9.1,
!    Springer, 2000.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Output, integer BOUND_NUM, the number of edges that lie on the convex
!    hull of the triangulation.
!
  implicit none
!
  integer bound_num
  integer point_num
  integer tri_num
!
  bound_num = 2 * point_num - tri_num - 2

  return
end
subroutine triangulation_check ( point_num, tri_num, tri_vert, ierror )
!
!*******************************************************************************
!
!! TRIANGULATION_CHECK makes some simple checks on a triangulation.
!
!
!  Discussion:
!
!    Because this routine does not receive the physical coordinates of
!    the nodes, it cannot ensure that the triangulation is maximal,
!    that is, that no more triangles can be created.
!
!  Modified:
!
!    31 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points or nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up the
!    triangles.  These should be listed in counterclockwise order.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred, the triangulation is not valid.
!
  implicit none
!
  integer tri_num
  integer point_num
!
  integer bound_num
  integer euler
  integer i
  integer ierror
  integer tri_vert(3,tri_num)
  integer used(point_num)
!
  ierror = 0
!
!  Checks 1 and 2:
!  POINT_NUM must be at least 3.
!  TRI_NUM must be at least 1.
!
  if ( point_num < 3 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The number of nodes is less than 3!'
    return
  end if

  if ( tri_num < 1 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The number of triangles is less than 1!'
    return
  end if
!
!  Checks 3 and 4:
!  Verify that all node values are greater than or equal to 1
!  and less than or equal to POINT_NUM.
!
  if ( any ( tri_vert(1:3,1:tri_num) < 1 ) ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some vertices are less than 1!'
    return
  end if

  if ( any ( tri_vert(1:3,1:tri_num) > point_num ) ) then
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some vertices are greater than POINT_NUM!'
    return
  end if
!
!  Check 5:
!  Verify that every node is used at least once.
!
  used(1:point_num) = 0

  used(tri_vert(1,1:tri_num)) = used(tri_vert(1,1:tri_num)) + 1
  used(tri_vert(2,1:tri_num)) = used(tri_vert(2,1:tri_num)) + 1
  used(tri_vert(3,1:tri_num)) = used(tri_vert(3,1:tri_num)) + 1

  if ( any ( used(1:point_num) == 0 ) ) then
    ierror = 5
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some nodes are never used as triangle vertices!'
    return
  end if
!
!  Check 6:
!  Verify that no node is repeated in a triangle.
!
  do i = 1, tri_num
    if ( tri_vert(1,i) == tri_vert(2,i) .or. &
         tri_vert(2,i) == tri_vert(3,i) .or. &
         tri_vert(3,i) == tri_vert(1,i) ) then
      ierror = 6
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_CHECK - Fatal error!'
      write ( *, '(a)' ) '  A triangle contains a null edge!'
      return
    end if
  end do
!
!  Check 7:
!  Verify that no edge is repeated, and that repeated edges occur in
!  negated pairs.
!
  call triangulation_edge_check ( tri_num, tri_vert, bound_num, ierror )

  if ( ierror /= 0 ) then
    ierror = 7
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_CHECK - Fatal error!'
    write ( *, '(a)' ) '  Some edges are repeated or given in the wrong direction!'
    return
  end if
!
!  Check 8:
!  Does the triangulation satisfy Euler's criterion?
!  If not, then the triangulation is not proper.  (For instance, there
!  might be a hole in the interior.)
!
  euler = bound_num + tri_num + 2 - 2 * point_num

  if ( euler /= 0 ) then
    ierror = 8
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The triangulation does not satisfy Euler''s criterion!'
    return
  end if

  return
end
subroutine triangulation_edge_check ( tri_num, tri_vert, bound_num, ierror )
!
!*******************************************************************************
!
!! TRIANGULATION_EDGE_CHECK checks the edges of a triangulation.
!
!
!  Modified:
!
!    25 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
!
!    Output, integer BOUND_NUM, the number of edges that lie on the boundary.
!
!    Output, integer IERROR, an error flag.
!    0, no errors were detected.
!    nonzero, an error occurred.
!
  implicit none
!
  integer tri_num
!
  integer bound_num
  integer i
  integer ierror
  integer j
  integer k
  integer row(3*tri_num,3)
  integer tri
  integer tri_vert(3,tri_num)
!
  ierror = 0
!
!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,+1) or (J,I,-1),
!    (J,K,+1) or (K,J,-1),
!    (K,I,+1) or (I,K,-1)
!
!  where we choose (I,J,+1) if I < J, or else (J,I,-1) and so on.
!
  do tri = 1, tri_num

    i = tri_vert(1,tri)
    j = tri_vert(2,tri)
    k = tri_vert(3,tri)

    if ( i < j ) then
      row(3*(tri-1)+1,1:3) = (/ i, j, +1 /)
    else
      row(3*(tri-1)+1,1:3) = (/ j, i, -1 /)
    end if

    if ( j < k ) then
      row(3*(tri-1)+2,1:3) = (/ j, k, +1 /)
    else
      row(3*(tri-1)+2,1:3) = (/ k, j, -1 /)
    end if

    if ( k < i ) then
      row(3*(tri-1)+3,1:3) = (/ k, i, +1 /)
    else
      row(3*(tri-1)+3,1:3) = (/ i, k, -1 /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!
  call irow_sort_a ( 3*tri_num, 3*tri_num, 3, row )
!
!  Step 3.
!
!  If any record occurs twice, we have an error.
!  Unpaired records lie on the convex hull.
!
  i = 0
  bound_num = 0

  do while ( i < 3 * tri_num )

    i = i + 1

    if ( i == 3 * tri_num ) then

      bound_num = bound_num + 1

    else

      if ( row(i,1) == row(i+1,1) .and. row(i,2) == row(i+1,2) ) then

        if ( row(i,3) == row(i+1,3) ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TRIANGULATION_EDGE_CHECK - Fatal error!'
          write ( *, '(a)' ) '  An edge occurs twice.'
          return
        else
          i = i + 1
        end if

      else

        bound_num = bound_num + 1

      end if

    end if

  end do

  return
end
subroutine triangulation_example_del_2d ( point_num, xc, tri_num, tri_vert, &
  tri_nabe )
!
!*******************************************************************************
!
!! TRIANGULATION_EXAMPLE_DEL_2D sets up a sample Delaunay triangulation.
!
!
!  Discussion:
!
!    This subroutine is designed to be called twice.
!
!    Before the first call, set POINT_NUM = TRI_NUM = 0.  Call the routine.
!    On return, POINT_NUM and TRI_NUM will have their correct values.
!
!    Now you may allocate space for XC, TRI_VERT and TRI_NABE.
!
!    Then call again, using the correct values for POINT_NUM and TRI_NUM.
!    On return, XC, TRI_VERT and TRI_NABE will have their correct values.
!
!  Modified:
!
!    15 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer POINT_NUM, the number of points.  If POINT_NUM is
!    not the correct value on input, it will be set to the correct value,
!    and the routine will return immediately.
!
!    Output, real XC(2,POINT_NUM), the coordinates of the points.
!
!    Output, integer TRI_NUM, the number of triangles.  If TRI_NUM is
!    not the correct value on input, it will be set to the correct value,
!    and the routine will return immediately.
!
!    Output, integer TRI_VERT(3,TRI_NUM), the nodes that make up the triangles.
!
!    Output, integer TRI_NABE(3,TRI_NUM), the triangle neighbors on each side.
!    Negative values indicate edges that lie on the exterior.
!
  implicit none
!
  integer, parameter :: point_num_internal = 13
  integer, parameter :: tri_num_internal = 16
!
  integer point_num
  integer tri_nabe(3,tri_num_internal)
  integer tri_num
  integer tri_vert(3,tri_num_internal)
  real xc(2,point_num_internal)
!
  if ( point_num /= point_num_internal ) then
    point_num = point_num_internal
    tri_num = tri_num_internal
    return
  end if

  if ( tri_num /= tri_num_internal ) then
    tri_num = tri_num_internal
    return
  end if

  xc = reshape ( (/ &
       0.0E+00, 0.0E+00, &
       2.0E+00, 2.0E+00, &
      -1.0E+00, 3.0E+00, &
      -2.0E+00, 2.0E+00, &
       8.0E+00, 2.0E+00, &
       9.0E+00, 5.0E+00, &
       7.0E+00, 4.0E+00, &
       5.0E+00, 6.0E+00, &
       6.0E+00, 7.0E+00, &
       8.0E+00, 8.0E+00, &
      11.0E+00, 7.0E+00, &
      10.0E+00, 4.0E+00, &
       6.0E+00, 4.0E+00 /), (/ 2, point_num /) )

  tri_vert(1:3,1:tri_num ) = reshape ( (/ &
     3,   4,   1, &
     3,   1,   2, &
     3,   2,   8, &
     2,   1,   5, &
     8,   2,  13, &
     8,  13,   9, &
     3,   8,   9, &
    13,   2,   5, &
     9,  13,   7, &
     7,  13,   5, &
     6,   7,   5, &
     9,   7,   6, &
    10,   9,   6, &
     6,   5,  12, &
    11,   6,  12, &
    10,   6,  11 /), (/ 3, tri_num /) )

  tri_nabe(1:3,1:tri_num) = reshape ( (/ &
       -4,  -13,    2, &
        1,    4,    3, &
        2,    5,    7, &
        2,  -43,    8, &
        3,    8,    6, &
        5,    9,    7, &
        3,    6,   -3, &
        5,    4,   10, &
        6,   10,   12, &
        9,    8,   11, &
       12,   10,   14, &
        9,   11,   13, &
      -23,   12,   16, &
       11,  -47,   15, &
       16,   14,  -50, &
       13,   15,  -39 /), (/ 3, tri_num /) )

  return
end
subroutine triangulation_nabe_nodes ( point_num, tri_num, tri_vert, &
  nabes_first, nabes_num, nabes_max, nabes_dim, nabes )
!
!*******************************************************************************
!
!! TRIANGULATION_NABE_NODES determines the neighbors of triangulation nodes.
!
!
!  Example:
!
!    On input, the triangle data structure is:
!
!    Triangle  Nodes
!    --------  ----------
!     1        3,   4,   1
!     2        3,   1,   2
!     3        3,   2,   6
!     4        2,   1,   5
!     5        6,   2,   5
!
!  On output, the auxilliary neighbor arrays are:
!
!    Node  Num  First
!    ----  ---  -----
!     1     4     1
!     2     4     5
!     3     4     9
!     4     2    13
!     5     3    15
!     6     3    18
!
!  and the neighbor array is:
!
!    Position  Node
!    --------  ----
!
!     1        2
!     2        3
!     3        4
!     4        5
!    -----------
!     5        1
!     6        3
!     7        5
!     8        6
!    -----------
!     9        1
!    10        2
!    11        4
!    12        6
!    -----------
!    13        1
!    14        3
!    -----------
!    15        1
!    16        2
!    17        6
!    -----------
!    18        2
!    19        3
!    20        5
!
!  Modified:
!
!    18 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
!
!    Output, integer NABES_FIRST(POINT_NUM), the index in NABES of the first
!    neighbor in the list for each node.
!
!    Output, integer NABES_NUM(POINT_NUM), the number of neighbors of each node.
!
!    Input, integer NABES_MAX, the maximum dimension of NABES.
!
!    Output, integer NABES_DIM, the dimension of NABES.
!
!    Output, integer NABES(NABES_DIM), a list of the neighbors of all the nodes.
!    Neighbors of node 1 are listed first, and so on.
!
  implicit none
!
  integer nabes_max
  integer point_num
  integer tri_num
!
  integer i
  integer i_current
  integer j
  integer k
  integer nabe
  integer nabes(nabes_max)
  integer nabes1(nabes_max)
  integer nabes_dim
  integer nabes_first(point_num)
  integer nabes_num(point_num)
  integer nuniq
  integer tri
  integer tri_vert(3,tri_num)
!
!  Step 1.  From the triangle list (I,J,K)
!  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
!
  nabes_dim = 0
  do tri = 1, tri_num
    i = tri_vert(1,tri)
    j = tri_vert(2,tri)
    k = tri_vert(3,tri)
    nabes1(nabes_dim+1:nabes_dim+6) = (/ i, i, j, j, k, k /)
    nabes(nabes_dim+1:nabes_dim+6) = (/ j, k, i, k, i, j /)
    nabes_dim = nabes_dim + 6
  end do
!
!  Step 2. Dictionary sort the neighbor relations.
!
  call ivec2_sort_a ( nabes_dim, nabes1, nabes )
!
!  Step 3. Remove duplicate entries.
!
  call ivec2_unique ( nabes_dim, nabes1, nabes, nuniq )

  nabes_dim = nuniq
!
!  Step 4. Construct the NABES_NUM and NABES_FIRST data.
!
  nabes_num(1:point_num) = 0
  nabes_first(1:point_num) = 0
  i_current = 0
  do nabe = 1, nabes_dim
    i = nabes1(nabe)
    if ( i == i_current ) then
      nabes_num(i) = nabes_num(i) + 1
    else
      i_current = i
      nabes_first(i) = nabe
      nabes_num(i) = 1
    end if
  end do

  return
end
subroutine triangulation_nabe_nodes_print ( point_num, nabes_first, &
  nabes_num, nabes_dim, nabes )
!
!*******************************************************************************
!
!! TRIANGULATION_NABE_NODES_PRINT prints a triangulation node neighbor array.
!
!
!  Modified:
!
!    11 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, integer NABES_FIRST(POINT_NUM), the index in NABES of the first
!    neighbor in the list for each node.
!
!    Input, integer NABES_NUM(POINT_NUM), the number of neighbors of each node.
!
!    Input, integer NABES_DIM, the dimension of NABES.
!
!    Input, integer NABES(NABES_DIM), a list of the neighbors of all the nodes.
!    Neighbors of node 1 are listed first, and so on.
!
  implicit none
!
  integer nabes_dim
  integer point_num
!
  integer i
  integer nabes(nabes_dim)
  integer nabes_first(point_num)
  integer nabes_num(point_num)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Node based arrays:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Node  Neighbors  Index #1'
  write ( *, '(a)' ) ' '
  do i = 1, point_num
    write ( *, '(3i6)' ) i, nabes_num(i), nabes_first(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The raw neighbor array.'
  write ( *, '(a)' ) ' '
  do i = 1, nabes_dim
    write ( *, '(2i6)' ) i, nabes(i)
  end do

  return
end
subroutine triangulation_nabe_triangles ( tri_num, tri_vert, tri_nabe )
!
!*******************************************************************************
!
!! TRIANGULATION_NABE_TRIANGLES determines triangle neighbors.
!
!
!  Discussion:
!
!    A triangulation of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each triangle.  However, in some cases, it is necessary to know
!    triangle adjacency information, that is, which triangle, if any,
!    is adjacent to a given triangle on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 3 * TRI_NUM
!    data items.
!
!  Example:
!
!    The input information from TRI_VERT:
!
!    Triangle   Nodes
!    --------   ---------------
!     1         3      4      1
!     2         3      1      2
!     3         3      2      8
!     4         2      1      5
!     5         8      2     13
!     6         8     13      9
!     7         3      8      9
!     8        13      2      5
!     9         9     13      7
!    10         7     13      5
!    11         6      7      5
!    12         9      7      6
!    13        10      9      6
!    14         6      5     12
!    15        11      6     12
!    16        10      6     11
!
!    The output information in TRI_NABE:
!
!    Triangle  Neighboring Triangles
!    --------  ---------------------
!
!     1        -1     -1      2
!     2         1      4      3
!     3         2      5      7
!     4         2     -1      8
!     5         3      8      6
!     6         5      9      7
!     7         3      6     -1
!     8         5      4     10
!     9         6     10     12
!    10         9      8     11
!    11        12     10     14
!    12         9     11     13
!    13        -1     12     16
!    14        11     -1     15
!    15        16     14     -1
!    16        13     15     -1
!
!  Modified:
!
!    17 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
!
!    Output, integer TRI_NABE(3,TRI_NUM), the three triangles that are direct
!    neighbors of a given triangle.  TRI_NABE(1,I) is the index of the triangle
!    which touches side 1, defined by nodes 2 and 3, and so on.  TRI_NABE(1,I)
!    is negative if there is no neighbor on that side.  In this case, that
!    side of the triangle lies on the boundary of the triangulation.
!
  implicit none
!
  integer tri_num
!
  integer i
  integer irow
  integer j
  integer k
  integer row(3*tri_num,4)
  integer side1
  integer side2
  integer tri_nabe(3,tri_num)
  integer tri
  integer tri_vert(3,tri_num)
  integer tri1
  integer tri2
!
!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,1,T) or (J,I,1,T),
!    (J,K,2,T) or (K,J,2,T),
!    (K,I,3,T) or (I,K,3,T)
!
!  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
!
  do tri = 1, tri_num

    i = tri_vert(1,tri)
    j = tri_vert(2,tri)
    k = tri_vert(3,tri)

    if ( i < j ) then
      row(3*(tri-1)+1,1:4) = (/ i, j, 1, tri /)
    else
      row(3*(tri-1)+1,1:4) = (/ j, i, 1, tri /)
    end if

    if ( j < k ) then
      row(3*(tri-1)+2,1:4) = (/ j, k, 2, tri /)
    else
      row(3*(tri-1)+2,1:4) = (/ k, j, 2, tri /)
    end if

    if ( k < i ) then
      row(3*(tri-1)+3,1:4) = (/ k, i, 3, tri /)
    else
      row(3*(tri-1)+3,1:4) = (/ i, k, 3, tri /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on columns 1 and 2; the routine we call here
!  sorts on columns 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two triangles share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two rows of ROW that start out ( I, J, ?, ? ).  By sorting ROW,
!  we make sure that these two rows occur consecutively.  That will
!  make it easy to notice that the triangles are neighbors.
!
  call irow_sort_a ( 3*tri_num, 3*tri_num, 4, row )
!
!  Step 3. Neighboring triangles show up as consecutive rows with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in TRI_NABE.
!
  tri_nabe(1:3,1:tri_num) = -1

  irow = 1

  do

    if ( irow >= 3*tri_num ) then
      exit
    end if

    if ( row(irow,1) /= row(irow+1,1) .or. row(irow,2) /= row(irow+1,2) ) then
      irow = irow + 1
      cycle
    end if

    side1 = row(irow,3)
    tri1 = row(irow,4)
    side2 = row(irow+1,3)
    tri2 = row(irow+1,4)

    tri_nabe(side1,tri1) = tri2
    tri_nabe(side2,tri2) = tri1

    irow = irow + 2

  end do

  return
end
subroutine triangulation_neighbor ( tri_num, tri_vert, t1, s1, t2, s2 )
!
!*******************************************************************************
!
!! TRIANGULATION_NEIGHBOR determines a neighbor of a given triangle.
!
!
!  Discussion:
!
!    A set of nodes is given.  A triangulation of the nodes has been
!    defined and recorded in TRI_VERT.  The TRI_VERT data structure records
!    triangles as sets of three nodes, N1, N2, N3, that implicitly define three
!    sides, being the line segments N1-N2, N2-N3, and N3-N1.
!
!    The nodes of the triangle are listed in counterclockwise order.
!    This means that if two triangles share a side, then the nodes
!    defining that side occur in the order (N1,N2) for one triangle,
!    and (N2,N1) for the other.
!
!    The routine is given a triangle and a side, and asked to find
!    another triangle (if any) that shares that side.  The routine
!    simply searches the TRI_VERT structure for an occurrence of the
!    nodes in the opposite order.
!
!  Modified:
!
!    23 August 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input/output, integer TRI_VERT(3,TRI_NUM), the nodes that define
!    each triangle.
!
!    Input, integer T1, the index of the triangle.
!
!    Input, integer S1, the index of the triangle side.
!
!    Output, integer T2, the index of the triangle which is the neighbor
!    to T1 on side S1, or 0 if there is no such neighbor.
!
!    Output, integer S2, the index of the side of triangle T2 which
!    is shared with triangle T1, or 0 if there is no such neighbor.
!
  implicit none
!
  integer tri_num
!
  integer i_wrap
  integer n1
  integer n2
  integer s
  integer s1
  integer s2
  integer ss
  integer t
  integer t1
  integer t2
  integer tri_vert(3,tri_num)
!
  t2 = 0
  s2 = 0

  n1 = tri_vert(s1,t1)
  ss = i_wrap ( s1+1, 1, 3 )
  n2 = tri_vert(ss,t1)

  do t = 1, tri_num
    do s = 1, 3
      if ( tri_vert(s,t) == n1 ) then
        ss = i_wrap ( s-1, 1, 3 )
        if ( tri_vert(ss,t) == n2 ) then
          t2 = t
          s2 = ss
          return
        end if
      end if
    end do
  end do

  return
end
subroutine triangulation_print ( point_num, xc, tri_num, tri_vert, tri_nabe )
!
!*******************************************************************************
!
!! TRIANGULATION_PRINT prints out information defining a Delaunay triangulation.
!
!
!  Discussion:
!
!    Triangulations created by RTRIS include extra information encoded
!    in the negative values of TRI_NABE.
!
!    Because some of the nodes counted in POINT_NUM may not actually be
!    used in the triangulation, I needed to compute the true number
!    of vertices.  I added this calculation on 13 October 2001.
!
!  Modified:
!
!    19 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, real XC(2,POINT_NUM), the point coordinates.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up the triangles.
!
!    Input, integer TRI_NABE(3,TRI_NUM), the triangle neighbors on each side.
!    If there is no triangle neighbor on a particular side, the value of
!    TRI_NABE should be negative.  If the triangulation data was created by
!    RTRIS2, then there is more information encoded in the negative values.
!
  implicit none
!
  integer point_num
  integer tri_num
!
  integer boundary_num
  integer i
  integer i_wrap
  integer j
  integer k
  integer n1
  integer n2
  integer s
  logical skip
  integer t
  integer tri_nabe(3,tri_num)
  integer tri_vert(3,tri_num)
  integer, allocatable, dimension ( : ) :: vertex_list
  integer vertex_num
  real xc(2,point_num)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_PRINT'
  write ( *, '(a)' ) '  Information defining a triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points is ', point_num

  call rmat_print ( point_num, point_num, 2, transpose ( xc ), &
    '  Point coordinates (transpose of internal array)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of triangles is ', tri_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sets of three points are used as vertices of'
  write ( *, '(a)' ) '  the triangles.  For each triangle, the points'
  write ( *, '(a)' ) '  are listed in counterclockwise order.'

  call imat_print ( tri_num, tri_num, 3, transpose ( tri_vert ), &
    '  Nodes that make up triangles (transpose of internal array)' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  On each side of a given triangle, there is either'
  write ( *, '(a)' ) '  another triangle, or a piece of the convex hull.'
  write ( *, '(a)' ) '  For each triangle, we list the indices of the three'
  write ( *, '(a)' ) '  neighbors, or (if negative) the codes of the'
  write ( *, '(a)' ) '  segments of the convex hull.'

  call imat_print ( tri_num, tri_num, 3, transpose ( tri_nabe ), &
    '  Indices of neighboring triangles (transpose of internal array)' )
!
!  Determine the number of vertices.  This is not the same as the
!  number of points!
!
  allocate ( vertex_list(1:3*tri_num) )

  vertex_list(1:3*tri_num) = reshape ( tri_vert(1:3,1:tri_num), &
    (/ 3*tri_num /) )

  call ivec_sort_heap_a ( 3*tri_num, vertex_list )

  call ivec_unique ( 3*tri_num, vertex_list, vertex_num )

  deallocate ( vertex_list )
!
!  Determine the number of boundary points.
!
  boundary_num = 2 * vertex_num - tri_num - 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of boundary points is ', boundary_num
   
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The segments that make up the convex hull can be'
  write ( *, '(a)' ) '  determined from the negative entries of the triangle'
  write ( *, '(a)' ) '  neighbor list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  # Tri Side  N1  N2'
  write ( *, '(a)' ) ' '

  skip = .false.

  k = 0

  do i = 1, tri_num

    do j = 1, 3

      if ( tri_nabe(j,i) < 0 ) then
        s = - tri_nabe(j,i)
        t = s / 3

        if ( t < 1 .or. t > tri_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Sorry, this data does not use the RTRIS2 convention'
          write ( *, '(a)' ) '  for convex hull segments.'
          skip = .true.
          exit
        end if

        s = mod ( s, 3 ) + 1
        k = k + 1
        n1 = tri_vert(s,t)
        n2 = tri_vert(i_wrap(s+1,1,3),t)
        write ( *, '(5i4)' ) k, t, s, n1, n2
      end if

    end do

    if ( skip ) then
      exit
    end if

  end do

  return
end
subroutine triangulation_sample_2d ( point_num, xc, tri_num, tri_vert, &
  num_ran, xd, td )
!
!*******************************************************************************
!
!! TRIANGULATION_SAMPLE_2D returns random points in a triangulation.
!
!
!  Discussion:
!
!    It is assumed that the triangulation consists of a set of non-overlapping
!    triangles.
!
!    The point is chosen uniformly in the area covered by the triangulation.
!
!  Modified:
!
!    13 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of points used in the triangulation.
!
!    Input, real XC(2,POINT_NUM), the point coordinates.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up the triangles.
!
!    Input, integer NUM_RAN, the number of points to sample.
!
!    Output, real XD(2,NUM_RAN), the sample points.
!
!    Output, integer TD(NUM_RAN), the triangle to which each sample point
!    belongs.
!
  implicit none
!
  integer point_num
  integer num_ran
  integer tri_num
!
  real area
  real area_cum(0:tri_num)
  real area_total
  integer i
  integer i1
  integer i2
  integer i3
  integer left
  real r
  integer right
  integer td(num_ran)
  integer tri_vert(3,tri_num)
  real x
  real x1
  real x2
  real x3
  real xc(2,point_num)
  real xd(2,num_ran)
  real y
  real y1
  real y2
  real y3
!
!  Compute the areas of the triangles.
!  Build a cumulative area vector.
!  Convert it to a relative cumulative area vector.
!
  area_cum(0) = 0.0E+00

  do i = 1, tri_num

    i1 = tri_vert(1,i)
    i2 = tri_vert(2,i)
    i3 = tri_vert(3,i)

    call triangle_area_2d ( xc(1,i1), xc(2,i1), xc(1,i2), xc(2,i2), xc(1,i3), &
      xc(2,i3), area )

    area_cum(i) = area_cum(i-1) + area

  end do

  area_total = area_cum(tri_num)

  area_cum(0:tri_num) = area_cum(0:tri_num) / area_total
!
!  Pick random values.  A random value R indicates the corresponding triangle
!  whose cumulative relative area contains R.
!
!  Bracket the random value in the cumulative relative areas,
!  indicating a triangle.
!
!  Pick a random point in the triangle.
!
  do i = 1, num_ran

    call random_number ( harvest = r )

    call rvec_bracket ( tri_num+1, area_cum, r, left, right )

    td(i) = right - 1

    i1 = tri_vert(1,td(i))
    x1 = xc(1,i1)
    y1 = xc(2,i1)

    i2 = tri_vert(2,td(i))
    x2 = xc(1,i2)
    y2 = xc(2,i2)

    i3 = tri_vert(3,td(i))
    x3 = xc(1,i3)
    y3 = xc(2,i3)

    call triangle_sample_2d ( x1, y1, x2, y2, x3, y3, x, y )

    xd(1,i) = x
    xd(2,i) = y

  end do

  return
end
subroutine triangulation_search_2d ( point_num, xc, tri_num, tri_vert, &
  tri_nabe, x, y, triangle, edge )
!
!*******************************************************************************
!
!! TRIANGULATION_SEARCH_2D searches a triangulation for a point.
!
!
!  Purpose:
!
!    Walk through neighboring triangles of a 2D Delaunay triangulation 
!    until a triangle is found containing point (X,Y), or (X,Y) is found 
!    to be outside the convex hull.  
!
!    The algorithm computes the barycentric coordinates of the point with 
!    respect to the current triangle.  If all three quantities are positive,
!    the point is contained in the triangle.  If the I-th coordinate is
!    negative, then (X,Y) lies on the far side of edge I, which is opposite
!    from vertex I.  This gives a hint as to where to search next.
!
!    For a Delaunay triangulation, the search is guaranteed to terminate.
!    For other triangulations, a cycle may occur.
!
!    Note the surprising fact that, even for a Delaunay triangulation of
!    a set of points, the nearest point to (X,Y) need not be one of the
!    vertices of the triangle containing (X,Y).  
!
!  Modified:
!
!    20 July 2001
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer POINT_NUM, the number of nodes.
!
!    Input, real XC(2,1:POINT_NUM), the coordinates of 2-D vertices.
!
!    Input, integer TRI_NUM, the number of triangles in the triangulation.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
!
!    Input, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list.
!
!    Input, real X, Y, the coordinates of a point.
!
!    Output, integer TRIANGLE, the index of the triangle where the search ended.
!    If a cycle occurred, then TRIANGLE = -1.
!
!    Output, integer EDGE, indicates the position of the point (X,Y) in
!    triangle TRIANGLE:
!    0, the interior or boundary of the triangle;
!    -1, outside the convex hull of the triangulation, past edge 1;
!    -2, outside the convex hull of the triangulation, past edge 2;
!    -3, outside the convex hull of the triangulation, past edge 3.
!
  implicit none
!
  integer point_num
  integer tri_num
!
  integer a
  real alpha
  integer b
  real beta
  integer c
  integer count
  real det
  real dx
  real dxa
  real dxb
  real dy
  real dya
  real dyb
  real gamma
  integer edge
  integer tri_vert(3,tri_num)
  integer tri_nabe(3,tri_num)
  integer triangle
  real x
  real xc(2,point_num)
  real y
!
  count = 0
  edge = 0
  call i_random ( 1, tri_num, triangle )

  do

    count = count + 1

    if ( count > tri_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_SEARCH_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm seems to be cycling.'
      triangle = -1
      edge = -1
      return
    end if
!
!  Get the vertices of triangle TRIANGLE.
!
    a = tri_vert(1,triangle)
    b = tri_vert(2,triangle)
    c = tri_vert(3,triangle)
!
!  Using vertex C as a base, compute the distances to vertices A and B,
!  and the point (X,Y).
!
    dxa = xc(1,a) - xc(1,c)
    dya = xc(2,a) - xc(2,c)

    dxb = xc(1,b) - xc(1,c)
    dyb = xc(2,b) - xc(2,c)

    dx = x - xc(1,c)
    dy = y - xc(2,c)

    det = dxa * dyb - dya * dxb
!
!  Compute the barycentric coordinates of the point (X,Y) with respect
!  to this triangle.
!
    alpha = ( dx * dyb - dy * dxb ) / det
    beta = ( dxa * dy - dya * dx ) / det
    gamma = 1.0E+00 - alpha - beta
!
!  If the barycentric coordinates are all positive, then the point
!  is inside the triangle and we're done.
!
    if ( alpha >= 0.0E+00 .and. beta >= 0.0E+00 .and. gamma >= 0.0E+00 ) then
      exit
    end if
!
!  At least one barycentric coordinate is negative.
!
!  If there is a negative barycentric coordinate for which there exists
!  an opposing triangle neighbor closer to the point, move to that triangle.
!
!  (Two coordinates could be negative, in which case we could go for the
!  most negative one, or the most negative one normalized by the actual
!  distance it represents).
!
    if ( alpha < 0.0E+00 .and. tri_nabe(2,triangle) > 0 ) then
      triangle = tri_nabe(2,triangle)
      cycle
    else if ( beta < 0.0E+00 .and. tri_nabe(3,triangle) > 0 ) then
      triangle = tri_nabe(3,triangle)
      cycle
    else if ( gamma < 0.0E+00 .and. tri_nabe(1,triangle) > 0 ) then
      triangle = tri_nabe(1,triangle)
      cycle
    end if
!
!  All negative barycentric coordinates correspond to vertices opposite
!  sides on the convex hull.
!
!  Note the edge and exit.
!
    if ( alpha < 0.0E+00 ) then
      edge = -2
      exit
    else if ( beta < 0.0E+00 ) then
      edge = -3
      exit
    else if ( gamma < 0.0E+00 ) then
      edge = -1
      exit
    end if

  end do

  return
end
subroutine tube_2d ( dist, n, x, y, x1, y1, x2, y2 )
!
!*******************************************************************************
!
!! TUBE_2D constructs a "tube" of given width around a path in 2D.
!
!
!  Discussion:
!
!    The routine is given a sequence of N points, and a distance DIST.
!
!    It returns the coordinates of the corners of the top and bottom
!    of a tube of width 2*DIST, which envelopes the line connecting
!    the points.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real DIST, the radius of the tube.
!
!    Input, integer N, the number of points defining the line.
!    N must be at least 2.
!
!    Input, real X(N), Y(N), the points which comprise the broken
!    line which is to be surrounded by the tube.  Points should
!    not be immediately repeated, that is, it should never be
!    the case that
!      X(I) = X(I+1) and Y(I) = Y(I+1).
!
!    Output, real X1(N), Y1(N), X2(N), Y2(N), the points ( X1(I), Y1(N) )
!    form one side of the tube, and ( X2(I), Y2(I) ) the other.
!
  implicit none
!
  integer n
!
  real a
  real b
  real c
  real dis1
  real dis2
  real dist
  real enorm0_2d
  integer i
  real temp
  real x(n)
  real x1(n)
  real x2(n)
  real xi
  real xim1
  real xip1
  real y(n)
  real y1(n)
  real y2(n)
  real yi
  real yim1
  real yip1
!
!  Check that N is at least 3.
!
  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUBE_2D - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 3'
    write ( *, '(a,i6)' ) '  but your input value was N = ', n
    stop
  end if
!
!  Check that consecutive points are distinct.
!
  do i = 1, n-1
    if ( x(i) == x(i+1) .and. y(i) == y(i+1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUBE_2D - Fatal error!'
      write ( *, '(a,i6)' ) '  X(I) = X(I+1) and Y(I) = Y(I+1) for I = ', i
      write ( *, '(a,2g14.6)' ) '  X(I), Y(I) = ', x(i), y(i)
      stop
    end if
  end do

  do i = 1, n

    if ( i == 1 ) then
      xim1 = x(i)
      yim1 = y(i)
    else
      xim1 = x(i-1)
      yim1 = y(i-1)
    end if

    xi = x(i)
    yi = y(i)

    if ( i < n ) then
      xip1 = x(i+1)
      yip1 = y(i+1)
    else
      xip1 = x(i)
      yip1 = y(i)
    end if

    call corpl_2d ( dist, xim1, yim1, xi, yi, xip1, yip1, x1(i), &
      y1(i), x2(i), y2(i) )
!
!  On the first and last steps, translate the corner points DIST units
!  along the line, to make an extra buffer.
!
    if ( i == 1 ) then

      temp = enorm0_2d ( x(1), y(1), x(2), y(2) )
      x1(1) = x1(1) - dist * ( x(2) - x(1) ) / temp
      y1(1) = y1(1) - dist * ( y(2) - y(1) ) / temp
      x2(1) = x2(1) - dist * ( x(2) - x(1) ) / temp
      y2(1) = y2(1) - dist * ( y(2) - y(1) ) / temp

    else if ( i == n ) then

      temp = enorm0_2d ( x(n-1), y(n-1), x(n), y(n) )
      x1(n) = x1(n) + dist * ( x(n) - x(n-1) ) / temp
      y1(n) = y1(n) + dist * ( y(n) - y(n-1) ) / temp
      x2(n) = x2(n) + dist * ( x(n) - x(n-1) ) / temp
      y2(n) = y2(n) + dist * ( y(n) - y(n-1) ) / temp

    end if
!
!  The new points ( X1(I), Y1(I) ) and ( X2(I), Y2(I) ) may need to be
!  swapped.
!
!  Compute the signed distance from the points to the line.
!
    if ( i > 1 ) then

      a = y(i-1) - y(i)
      b = x(i) - x(i-1)
      c = x(i-1) * y(i) - x(i) * y(i-1)

      dis1 = ( a * x1(i-1) + b * y1(i-1) + c ) / sqrt ( a * a + b * b )

      dis2 = ( a * x1(i) + b * y1(i) + c ) / sqrt ( a * a + b * b )

      if ( sign ( 1.0E+00, dis1 ) /= sign ( 1.0E+00, dis2 ) ) then

        call r_swap ( x1(i), x2(i) )
        call r_swap ( y1(i), y2(i) )

      end if

    end if

  end do

  return
end
subroutine tuple_next2 ( n, xmin, xmax, x, rank )
!
!*******************************************************************************
!
!! TUPLE_NEXT2 computes the next element of an integer tuple space.
!
!
!  Discussion:
!
!    The elements X are N vectors.
!
!    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
!
!    The elements are produced one at a time.
!
!    The first element is
!      (XMIN(1), XMIN(2), ..., XMIN(N)),
!    the second is (probably)
!      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
!    and the last element is
!      (XMAX(1), XMAX(2), ..., XMAX(N))
!
!    Intermediate elements are produced in a lexicographic order, with
!    the first index more important than the last, and the ordering of
!    values at a fixed index implicitly defined by the sign of
!    XMAX(I) - XMIN(I).
!
!  Examples:
!
!    N = 2,
!    XMIN = (/ 1, 10 /)
!    XMAX = (/ 3,  8 /)
!
!    RANK    X
!    ----  -----
!      1   1 10
!      2   1  9
!      3   1  8
!      4   2 10
!      5   2  9
!      6   2  8
!      7   3 10
!      8   3  9
!      9   3  8
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, integer XMIN(N), XMAX(N), the "minimum" and "maximum" entry values.
!    These values are minimum and maximum only in the sense of the lexicographic
!    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
!    than XMAX(I).
!
!    Input/output, integer X(N), on input the previous tuple.
!    On output, the next tuple.
!
!    Input/output, integer RANK, the rank of the item.  On first call,
!    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
!    there are no more items in the sequence.
!
  implicit none
!
  integer n
!
  integer i
  integer rank
  integer x(n)
  integer xmin(n)
  integer xmax(n)
!
  if ( rank < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank > product ( 1 + abs ( xmax(1:n) - xmin(1:n) ) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank == 0 ) then
    x(1:n) = xmin(1:n)
    rank = 1
    return
  end if

  rank = rank + 1
  i = n

  do

    if ( x(i) /= xmax(i) ) then
      x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
      exit
    end if

    x(i) = xmin(i)

    if ( i == 1 ) then
      rank = 0
      exit
    end if

    i = i - 1

  end do

  return
end
subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
  ltri, ledg, rtri, redg )
!
!*******************************************************************************
!
!! VBEDG determines which boundary edges are visible to a point.
!
!
!  Discussion:
!
!    The point (X,Y) is assumed to be outside the convex hull of the
!    region covered by the 2D triangulation.
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Modified:
!
!    25 August 2001
!
!  Parameters:
!
!    Input, real X, Y, the coordinates of a point outside the convex hull
!    of the current triangulation.
!
!    Input, integer POINT_NUM, the number of points.
!
!    Input, real POINT_XY(2,POINT_NUM), the coordinates of the vertices.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the triangle incidence list.
!
!    Input, integer TRI_NABE(3,TRI_NUM), the triangle neighbor list; negative
!    values are used for links of a counter clockwise linked list of boundary
!    edges;
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer LTRI, LEDG.  If LTRI /= 0 then these values are
!    assumed to be already computed and are not changed, else they are updated.
!    On output, LTRI is the index of boundary triangle to the left of the
!    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
!    edge of triangle LTRI to the left of the leftmost boundary edge visible
!    from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer RTRI.  On input, the index of the boundary triangle
!    to begin the search at.  On output, the index of the rightmost boundary
!    triangle visible from (X,Y).
!
!    Input/output, integer REDG, the edge of triangle RTRI that is visible
!    from (X,Y).  1 <= REDG <= 3.
!
  implicit none
!
  integer point_num
  integer tri_num
!
  integer a
  integer b
  integer e
  integer i_wrap
  integer l
  logical ldone
  integer ledg
  integer lr
  integer lrline
  integer ltri
  real point_xy(2,point_num)
  integer redg
  integer rtri
  integer t
  integer tri_nabe(3,tri_num)
  integer tri_vert(3,tri_num)
  real x
  real y
!
!  Find the rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -tri_nabe(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = tri_vert(e,t)

    if ( e <= 2 ) then
      b = tri_vert(e+1,t)
    else
      b = tri_vert(1,t)
    end if

    lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &     
      point_xy(2,b), 0.0E+00 )

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = tri_vert(e,t)
    e = i_wrap ( e-1, 1, 3 )

    do while ( tri_nabe(e,t) > 0 )

      t = tri_nabe(e,t)

      if ( tri_vert(1,t) == b ) then
        e = 3
      else if ( tri_vert(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = tri_vert(e,t)

    lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
       point_xy(2,b), 0.0E+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
subroutine vector_directions_2d ( x1, y1, ax, ay )
!
!*******************************************************************************
!
!! VECTOR_DIRECTIONS_2D returns the direction angles of a vector in 2D.
!
!
!  Discussion:
!
!    Let V be the vector, and let E(I) be the I-th unit coordinate axis vector.
!    The I-th direction angle is the angle between V and E(I), which is
!    the angle whose cosine is equal to the direction cosine:
!
!      Direction_Cosine(I) = V dot E(I) / |V|.
!
!    If V is the null or zero vector, then the direction cosines and
!    direction angles are undefined, and this routine simply returns
!    zeroes.
!
!  Modified:
!
!    16 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, the coordinates of the endpoint of the vector V
!    (which is presumed to begin at the origin).
!
!    Output, real AX, AY, the direction angles, in radians,
!    that the vector V makes with the X and Y axes.
!
  implicit none
!
  real ax
  real ay
  real x1
  real y1
!
  ax = atan2 ( y1, x1 )
  ay = atan2 ( x1, y1 )

  return
end
subroutine vector_directions_3d ( x1, y1, z1, ax, ay, az )
!
!*******************************************************************************
!
!! VECTOR_DIRECTIONS_3D returns the direction angles of a vector in 3D.
!
!
!  Discussion:
!
!    Let V be the vector, and let E(I) be the I-th unit coordinate axis vector.
!    The I-th direction angle is the angle between V and E(I), which is
!    the angle whose cosine is equal to the direction cosine:
!
!      Direction_Cosine(I) = V dot E(I) / |V|.
!
!    If V is the null or zero vector, then the direction cosines and
!    direction angles are undefined, and this routine simply returns
!    zeroes.
!
!  Modified:
!
!    17 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, the coordinates of the endpoint of the vector V
!    (which is presumed to begin at the origin).
!
!    Output, real AX, AY, AZ, the direction angles, in radians, that the
!    vector V makes with the X, Y and Z axes.
!
  implicit none
!
  real arc_cosine
  real ax
  real ay
  real az
  real cos_x
  real cos_y
  real cos_z
  real v1norm
  real v2norm
  real x1
  real y1
  real z1
!
!  Get the norm of the vector.
!
  v1norm = sqrt ( x1 * x1 + y1 * y1 + z1 * z1 )

  if ( v1norm == 0.0E+00 ) then
    ax = 0.0E+00
    ay = 0.0E+00
    az = 0.0E+00
    return
  end if

  v2norm = 1.0E+00
!
!  Get the direction cosines.
!
  cos_x = x1 / ( v1norm * v2norm )
  cos_y = y1 / ( v1norm * v2norm )
  cos_z = z1 / ( v1norm * v2norm )
!
!  Retrieve the direction angles.
!
  ax = arc_cosine ( cos_x )
  ay = arc_cosine ( cos_y )
  az = arc_cosine ( cos_z )

  return
end
subroutine vector_rotate_2d ( x1, y1, angle, x2, y2 )
!
!*******************************************************************************
!
!! VECTOR_ROTATE_2D rotates a vector around the origin in 2D.
!
!
!  Discussion:
!
!    To see why this formula is so, consider that the original point
!    has the form ( R cos Theta, R sin Theta ), and the rotated point
!    has the form ( R cos ( Theta + Angle ), R sin ( Theta + Angle ) ).
!    Now use the addition formulas for cosine and sine to relate
!    the new point to the old one:
!
!      ( X2 ) = ( cos Angle  - sin Angle ) * ( X1 )
!      ( Y2 )   ( sin Angle    cos Angle )   ( Y1 )
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, the components of the vector to be rotated.
!
!    Input, real ANGLE, the angle, in radians, of the rotation to be
!    carried out.  A positive angle rotates the vector in the
!    counterclockwise direction.
!
!    Output, real X2, Y2, the rotated vector.
!
  implicit none
!
  real angle
  real x1
  real x2
  real y1
  real y2
!
  x2 = cos ( angle ) * x1 - sin ( angle ) * y1
  y2 = sin ( angle ) * x1 + cos ( angle ) * y1

  return
end
subroutine vector_rotate_base_2d ( x1, y1, xb, yb, angle, x2, y2 )
!
!*******************************************************************************
!
!! VECTOR_ROTATE_BASE_2D rotates a vector around a base point in 2D.
!
!
!  Discussion:
!
!    The original vector is assumed to be ( X1-XB, Y1-YB ), and the
!    rotated vector is ( X2-XB, Y2-YB ).
!
!  Modified:
!
!    13 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, the endpoint of the original vector.
!
!    Input, real XB, YB, the location of the base point.
!
!    Input, real ANGLE, the angle, in radians, of the rotation to be
!    carried out.  A positive angle rotates the vector in the
!    counterclockwise direction.
!
!    Output, real X2, Y2, the endpoint of the rotated vector.
!
  implicit none
!
  real angle
  real x1
  real x2
  real xb
  real y1
  real y2
  real yb
!
  x2 = xb + cos ( angle ) * ( x1 - xb ) - sin ( angle ) * ( y1 - yb )
  y2 = yb + sin ( angle ) * ( x1 - xb ) + cos ( angle ) * ( y1 - yb )

  return
end
subroutine vector_separation_3d ( v1, v2, theta )
!
!*******************************************************************************
!
!! VECTOR_SEPARATION_3D finds the angular separation between vectors in 3D.
!
!
!  Discussion:
!
!    Any two vectors lie in a plane, and are separated by a plane angle.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real V1(3), V2(3), the two vectors.
!
!    Output, real THETA, the angle between the two vectors.
!
  implicit none
!
  integer, parameter :: n = 3
!
  real arc_cosine
  real cos_theta
  real theta
  real v1(n)
  real v1_norm
  real v2(n)
  real v2_norm
!
  v1_norm = sqrt ( sum ( v1(1:n)**2 ) )

  v2_norm = sqrt ( sum ( v2(1:n)**2 ) )

  cos_theta = dot_product ( v1(1:n), v2(1:n) ) / ( v1_norm * v2_norm )

  theta = arc_cosine ( cos_theta )

  return
end
subroutine vector_separation_nd ( n, v1, v2, theta )
!
!*******************************************************************************
!
!! VECTOR_SEPARATION_ND finds the angular separation between vectors in ND.
!
!
!  Discussion:
!
!    Any two vectors lie in a plane, and are separated by a plane angle.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real V1(N), V2(N), the two vectors.
!
!    Output, real THETA, the angle between the two vectors.
!
  implicit none
!
  integer n
!
  real arc_cosine
  real cos_theta
  real theta
  real v1(n)
  real v1_norm
  real v2(n)
  real v2_norm
!
  v1_norm = sqrt ( sum ( v1(1:n)**2 ) )

  v2_norm = sqrt ( sum ( v2(1:n)**2 ) )

  cos_theta = dot_product ( v1(1:n), v2(1:n) ) / ( v1_norm * v2_norm )

  theta = arc_cosine ( cos_theta )

  return
end
subroutine vector_unit_2d ( x1, y1 )
!
!*******************************************************************************
!
!! VECTOR_UNIT_2D normalizes a vector in 2D.
!
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X1, Y1, the components of a 2D vector to be
!    normalized.  On output, V should have unit Euclidean norm.
!    However, if the input vector has zero Euclidean norm, it is
!    not altered.
!
  implicit none
!
  real temp
  real x1
  real y1
!
  temp = sqrt ( x1 * x1 + y1 * y1 )

  if ( temp /= 0.0E+00 ) then
    x1 = x1 / temp
    y1 = y1 / temp
  end if

  return
end
subroutine vector_unit_3d ( x1, y1, z1 )
!
!*******************************************************************************
!
!! VECTOR_UNIT_3D normalizes a vector in 3D.
!
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X1, Y1, Z1, the components of a 3D vector to be
!    normalized.  On output, V should have unit Euclidean norm.
!    However, if the input vector has zero Euclidean norm, it is
!    not altered.
!
  implicit none
!
  real temp
  real x1
  real y1
  real z1
!
  temp = sqrt ( x1 * x1 + y1 * y1 + z1 * z1 )

  if ( temp /= 0.0E+00 ) then
    x1 = x1 / temp
    y1 = y1 / temp
    z1 = z1 / temp
  end if

  return
end
subroutine vector_unit_nd ( n, v )
!
!*******************************************************************************
!
!! VECTOR_UNIT_ND normalizes a vector in ND.
!
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, real V(N), the vector to be normalized.  On output,
!    V should have unit Euclidean norm.  However, if the input vector
!    has zero Euclidean norm, it is not altered.
!
  implicit none
!
  integer n
!
  real enorm_nd
  integer i
  real temp
  real v(n)
!
  temp = enorm_nd ( n, v )

  if ( temp /= 0.0E+00 ) then
    v(1:n) = v(1:n) / temp
  end if

  return
end
subroutine voxel_line_3d ( x1, y1, z1, x2, y2, z2 )
!
!*******************************************************************************
!
!! VOXEL_LINE_3D computes the voxels that line on a line in 3D.
!
!
!  Discussion:
!
!    This preliminary version of the routine simply prints out the
!    intermediate voxels, rather than returning them one at a time,
!    or all together.
!
!  Modified:
!
!    16 April 1999
!
!  Reference:
!
!    Daniel Cohen,
!    Voxel Traversal along a 3D Line,
!    Graphics Gems, 1994.
!
!  Parameters:
!
!    Input, integer X1, Y1, Z1, X2, Y2, Z2, the coordinates of the
!    voxels that begin and end the line.
!
  implicit none
!
  integer ax
  integer ay
  integer az
  integer exy
  integer exz
  integer ezy
  integer i
  integer n
  integer sx
  integer sy
  integer sz
  integer x
  integer x1
  integer x2
  integer y
  integer y1
  integer y2
  integer z
  integer z1
  integer z2
!
  sx = sign ( 1, x2 - x1 )
  sy = sign ( 1, y2 - y1 )
  sz = sign ( 1, z2 - z1 )

  ax = abs ( x2 - x1 )
  ay = abs ( y2 - y1 )
  az = abs ( z2 - z1 )

  exy = ay - ax
  exz = az - ax
  ezy = ay - az

  n = ax + ay + az

  x = x1
  y = y1
  z = z1

  do i = 0, n

    write ( *, '(4i6)' ) i, x, y, z

    if ( exy < 0 ) then

      if ( exz < 0 ) then
        x = x + sx
        exy = exy + 2 * ay
        exz = exz + 2 * az
      else
        z = z + sz
        exz = exz - 2 * ax
        exy = exy + 2 * ay
      end if

    else if ( ezy < 0 ) then
       z = z + sz
       exz = exz - 2 * ax
       ezy = ezy + 2 * ay
    else
      y = y + sy
      exy = exy - 2 * ax
      ezy = ezy - 2 * az
    end if

  end do

  return
end
subroutine voxel_region_3d ( ishow, list, maxlist, nlist, nregion, nx, ny, nz )
!
!*******************************************************************************
!
!! VOXEL_REGION_3D arranges a set of voxels into contiguous regions in 3D.
!
!
!  Modified:
!
!    29 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISHOW(NX,NY,NZ).
!
!    On input, ISHOW(I,J,K) has the values:
!      0, if the voxel is OFF;
!      anything else, if the voxel is ON.
!
!    On output, ISHOW(I,J,K) has the values:
!      0, if the voxel is off,
!      N, if the voxel is ON, and is part of region N.
!
!    Output, integer LIST(MAXLIST), contains, in stack form, a list
!    of the indices of the elements in each region.
!
!    The number of elements in NREGION is NELEM = LIST(NLIST).  The
!    (I,J,K) indices of the last element in this region are in
!    LIST(NLIST-3) through LIST(NLIST-1), and the first element is
!    listed in LIST(NLIST-3*NELEM), LIST(NLIST-3*NELEM+1),
!    LIST(NLIST-3*NELEM+2).
!
!    The number of elements in NREGION-1 is listed in
!    LIST(NLIST-3*NELEM-1), and so on.
!
!    Input, integer MAXLIST, the maximum length of the array used to
!    list the elements of the regions.
!
!    Output, integer NLIST, the number of entries of LIST that were used.
!    However, if NLIST > MAXLIST, then there was not enough space in
!    LIST to store the data properly, and LIST should not be used,
!    although the data in ISHOW should be correct.
!
!    Output, integer NREGION, the number of regions discovered.
!
!    Input, integer NX, NY, NZ, the number of voxels in the X, Y and
!    Z directions.
!
  implicit none
!
  integer, parameter :: maxstack = 100
!
  integer maxlist
  integer nx
  integer ny
  integer nz
!
  integer i
  integer i2
  integer ibase
  integer ihi
  integer ilo
  integer ishow(nx,ny,nz)
  integer j
  integer j2
  integer jbase
  integer jhi
  integer jlo
  integer k
  integer k2
  integer kbase
  integer khi
  integer klo
  integer list(maxlist)
  integer nabes
  integer ncan
  integer nelements
  integer nlist
  integer nregion
  integer nstack
  integer stack(maxstack)
!
!  Reset all nonzero entries of ISHOW to -1.
!
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz

        if ( ishow(i,j,k) /= 0 ) then
          ishow(i,j,k) = - 1
        end if

      end do
    end do
  end do
!
!  Start the number of items in the region list at 0.
!
  nlist = 0
!
!  Start the number of regions at 0.
!
  nregion = 0
!
!  The stack begins empty.
!
  nstack = 0
!
!  Search for an unused "ON" voxel from which we can "grow" a new region.
!
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
!
!  We found a voxel that is "ON", and does not belong to any region.
!
        if ( ishow(i,j,k) == - 1 ) then
!
!  Increase the number of regions.
!
          nregion = nregion + 1
!
!  Add this voxel to the region.
!
          ishow(i,j,k) = nregion
!
!  Add this voxel to the stack.
!
          if ( nstack + 4 > maxstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'VOXEL_REGION - Fatal error!'
            write ( *, '(a)' ) '  The internal stack overflowed.'
            write ( *, '(a)' ) '  The algorithm has failed.'
            stop
          end if

          stack(nstack+1) = i
          stack(nstack+2) = j
          stack(nstack+3) = k

          stack(nstack+4) = 1

          nstack = nstack + 4
!
!  Add this voxel to the description of the region.
!
          nelements = 1

          if ( nlist + 3 <= maxlist ) then
            list(nlist+1) = i
            list(nlist+2) = j
            list(nlist+3) = k
          end if

          nlist = nlist + 3

          do
!
!  Find all neighbors of BASE that are "ON" but unused.
!  Mark them as belonging to this region, and stack their indices.
!
            ibase = stack(nstack-3)
            jbase = stack(nstack-2)
            kbase = stack(nstack-1)

            ilo = max ( ibase-1, 1 )
            ihi = min ( ibase+1, nx )
            jlo = max ( jbase-1, 1 )
            jhi = min ( jbase+1, ny )
            klo = max ( kbase-1, 1 )
            khi = min ( kbase+1, nz )

            nabes = 0

            do i2 = ilo, ihi
              do j2 = jlo, jhi
                do k2 = klo, khi
!
!  We found a neighbor to our current search point, which is "ON" and unused.
!
                  if ( ishow(i2,j2,k2) == - 1 ) then
!
!  Increase the number of neighbors.
!
                    nabes = nabes + 1
!
!  Mark the neighbor as belonging to the region.
!
                    ishow(i2,j2,k2) = nregion
!
!  Add the neighbor to the stack.
!
                    if ( nstack+3 > maxstack ) then
                      write ( *, '(a)' ) ' '
                      write ( *, '(a)' ) 'VOXEL_REGION - Fatal error!'
                      write ( *, '(a)' ) '  The internal stack overflowed.'
                      write ( *, '(a)' ) '  The algorithm has failed.'
                      stop
                    end if

                    stack(nstack+1) = i2
                    stack(nstack+2) = j2
                    stack(nstack+3) = k2

                    nstack = nstack+3
!
!  Add the neighbor to the description of the region.
!
                    nelements = nelements + 1

                    if ( nlist+3 <= maxlist ) then
                      list(nlist+1) = i2
                      list(nlist+2) = j2
                      list(nlist+3) = k2
                    end if

                    nlist = nlist + 3

                  end if

                end do
              end do
            end do
!
!  If any new neighbors were found, take the last one as the basis
!  for a deeper search.
!
            if ( nabes > 0 ) then

              if ( nstack+1 > maxstack ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'VOXEL_REGION - Fatal error!'
                write ( *, '(a)' ) '  The internal stack overflowed.'
                write ( *, '(a)' ) '  The algorithm has failed.'
                stop
              end if

              stack(nstack+1) = nabes
              nstack = nstack + 1
              cycle

            end if
!
!  If the current search point had no new neighbors, drop it from the stack.
!
            ncan = stack(nstack) - 1
            nstack = nstack - 3
            stack(nstack) = ncan
!
!  If there are still any unused candidates at this level, take the
!  last one as the basis for a deeper search.
!
            if ( stack(nstack) > 0 ) then
              cycle
            end if
!
!  If there are no more unused candidates at this level, then we need
!  to back up a level in the stack.  If there are any candidates at
!  that earlier level, then we can still do more searching.
!
            nstack = nstack - 1

            if ( nstack <= 0 ) then
              exit
            end if

          end do
!
!  If we have exhausted the stack, we have completed this region.
!  Tag the number of elements to the end of the region description list.
!
          nlist = nlist + 1
          if ( nlist <= maxlist ) then
            list(nlist) = nelements
          end if

        end if

      end do
    end do
  end do
!
!  Print some warnings.
!
  if ( nlist > maxlist ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VOXEL_REGION - Warning!'
    write ( *, '(a)' ) '  MAXLIST was too small to list the regions.'
    write ( *, '(a)' ) '  Do not try to use the LIST array!'
    write ( *, '(a)' ) '  The ISHOW data is OK, however.'
  end if

  return
end
subroutine voxel_step_3d ( i1, j1, k1, i2, j2, k2, inc, jnc, knc )
!
!*******************************************************************************
!
!! VOXEL_STEP_3D computes voxels along a line from a given point in 3D.
!
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I1, J1, K1, the coordinates of the base voxel from
!    which the line begins.
!
!    Input/output, integer I2, J2, K2.
!
!    On input, these are the coordinates of the current voxel on
!    the line.  For the first call, these might be I1, J1 and K1.
!
!    On output, these are the coordinates of the next voxel along
!    the line.
!
!    Input, integer INC, JNC, KNC, the increments to the voxels.
!    These values define the direction along which the line proceeds.
!    However, the voxels on the line will typically be incremented
!    by a fractional value of the vector (INC,JNC,KNC), and the
!    result is essentially rounded.
!
!    If you input INC = JNC = KNC, then no movement is possible,
!    and none is made.
!
  implicit none
!
  real alpha
  real alphai
  real alphaj
  real alphak
  integer i1
  integer i2
  integer inc
  integer j1
  integer j2
  integer jnc
  integer k1
  integer k2
  integer knc
!
!  Assuming for the moment that (I,J,K) can take on real values,
!  points on the line have the form:
!
!    I = I1 + alpha * inc
!    J = J1 + alpha * jnc
!    K = K1 + alpha * knc
!
  if ( inc == 0 .and. jnc == 0 .and. knc == 0 ) then
    return
  end if

  alpha = 0.0E+00
!
!  Compute the smallest ALPHA that will change I2, J2 or K2 by +-0.5.
!
  if ( inc > 0 ) then
    alphai = ( real ( i2 - i1 ) + 0.5E+00 ) / real ( inc )
  else if ( inc < 0 ) then
    alphai = ( real ( i2 - i1 ) - 0.5E+00 ) / real ( inc )
  else
    alphai = huge ( alphai )
  end if

  if ( jnc > 0 ) then
    alphaj = ( real ( j2 - j1 ) + 0.5E+00 ) / real ( jnc )
  else if ( jnc < 0 ) then
    alphaj = ( real ( j2 - j1 ) - 0.5E+00 ) / real ( jnc )
  else
    alphaj = huge ( alphaj )
  end if

  if ( knc > 0 ) then
    alphak = ( real ( k2 - k1 ) + 0.5E+00 ) / real ( knc )
  else if ( knc < 0 ) then
    alphak = ( real ( k2 - k1 ) - 0.5E+00 ) / real ( knc )
  else
    alphaj = huge ( alphaj )
  end if
!
!  The ALPHA of smallest positive magnitude represents the closest next voxel.
!
  alpha = huge ( alpha )

  if ( alphai > 0.0E+00 ) then
    alpha = min ( alpha, alphai )
  end if

  if ( alphaj > 0.0E+00 ) then
    alpha = min ( alpha, alphaj )
  end if

  if ( alphak > 0.0E+00 ) then
    alpha = min ( alpha, alphak )
  end if
!
!  Move to the new voxel.  Whichever index just made the half
!  step must be forced to take a whole step.
!
  if ( alpha == alphai ) then
    i2 = i2 + sign ( 1, inc )
    j2 = j1 + nint ( alpha * jnc )
    k2 = k1 + nint ( alpha * knc )
  else if ( alpha == alphaj ) then
    i2 = i1 + nint ( alpha * inc )
    j2 = j2 + sign ( 1, jnc )
    k2 = k1 + nint ( alpha * knc )
  else if ( alpha == alphak ) then
    i2 = i1 + nint ( alpha * inc )
    j2 = j1 + nint ( alpha * jnc )
    k2 = k2 + sign ( 1, knc )
  end if

  return
end
subroutine xyz_to_radec ( x, y, z, ra, dec )
!
!*******************************************************************************
!
!! XYZ_TO_RADEC converts (X,Y,Z) to right ascension/declination coordinates.
!
!
!  Discussion:
!
!    Given an XYZ point, compute its distance R from the origin, and
!    regard it as lying on a sphere of radius R, whose axis is the Z
!    axis.
!
!    The right ascension of the point is the "longitude", measured in hours,
!    between 0 and 24, with the X axis having right ascension 0, and the
!    Y axis having right ascension 6.
!
!    Declination measures the angle from the equator towards the north pole,
!    and ranges from -90 (South Pole) to 90 (North Pole).
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, Y, Z, the coordinates of a point in 3D.
!
!    Output, real RA, DEC, the corresponding right ascension and declination.
!
  implicit none
!
  real atan4
  real dec
  real norm_v
  real phi
  real ra
  real radians_to_degrees
  real theta
  real x
  real y
  real z
!
  norm_v = sqrt ( x**2 + y**2 + z**2 )

  phi = asin ( z / norm_v )

  if ( cos ( phi ) == 0.0E+00 ) then
    theta = 0.0E+00
  else
    theta = atan4 ( y, x )
  end if

  dec = radians_to_degrees ( phi )
  ra = radians_to_degrees ( theta ) / 15.0E+00

  return
end
