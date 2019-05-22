module geomlib

   use genvar

   IMPLICIT NONE

contains
!
! !ARC_COSINE computes the arc cosine function, with argument truncation.
! !ANGLE_CONTAINS_RAY_2D determines if an angle contains a ray, in 2D.
! !ANGLE_DEG_2D returns the angle swept out between two rays in 2D.
! !ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
! !ANGLE_RAD_3D returns the angle in radians between two rays in 3D.
! !ANGLE_RAD_ND returns the angle in radians between two rays in ND.
! !ANGLEI_DEG_2D returns the interior angle in degrees between two rays in 2D.
! !ANGLEI_RAD_2D returns the interior angle in radians between two rays in 2D.
! !ATAN4 computes the inverse tangent of the ratio Y / X.
! !BASIS_MAP_3D computes the matrix which maps one basis to another.
! !BOX_CONTAINS_POINT_3D determines if a point is inside a parallelepiped in 3D.
! !BOX_POINT_DIST_3D: distance ( parallelepiped, point ) in 3D.
! !CIRCLE_AREA_2D computes the area of a circle in 2D.
! !CIRCLE_LUNE_AREA_2D returns the area of a circular lune in 2D.
! !CIRCLE_LUNE_CENTROID_2D returns the centroid of a circular lune in 2D.
! !CIRCLE_SECTOR_AREA_2D computes the area of a circular sector in 2D.
! !CIRCLE_SECTOR_CENTROID_2D returns the centroid of a circular sector in 2D.
! !CIRCLE_TRIANGLE_AREA_2D returns the area of a circle triangle in 2D.
! !CIRCLE_DIA2IMP_2D converts a diameter to an implicit circle in 2D.
! !CIRCLE_EXP_CONTAINS_POINT_2D determines if an explicit circle contains a point in 2D.
! !CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
! !CIRCLE_IMP_CONTAINS_POINT_2D determines if an implicit circle contains a point in 2D.
! !CIRCLE_IMP_LINE_PAR_INT_2D: ( implicit circle, parametric line ) intersection in 2D.
! !CIRCLE_IMP_POINTS_2D returns N equally spaced points on an implicit circle in 2D.
! !CIRCLE_IMP_POINTS_ARC_2D returns N points on an arc of an implicit circle in 2D.
! !CONE_AREA_3D computes the surface area of a right circular cone in 3D.
! !CONE_SECTOR_CENTROID_2D returns the centroid of a cone in 3D.
! !CONE_VOLUME_3D computes the volume of a right circular cone in 3D.
! !CONV3D converts 3D data to a 2D projection.
! !CORPL_2D "boxes" an angle defined by three points in 2D.
! !COT returns the cotangent of an angle.
! !COTD returns the cotangent of an angle given in degrees.
! !CROSS_2D finds the cross product of a pair of vectors in 2D.
! !CROSS_3D computes the cross product of two vectors in 3D.
! !CROSS0_2D finds the cross product of (P1-P0) and (P2-P0) in 2D.
! !CROSS0_3D computes the cross product of (P1-P0) and (P2-P0) in 3D.
! !CUBE_SHAPE_3D describes a cube in 3D.
! !DEGREES_TO_RADIANS converts an angle from degrees to radians.
! !DIRECTION_PERT_3D randomly perturbs a direction vector in 3D.
! !DIRECTION_RANDOM_3D picks a random direction vector in 3D.
! !DIRECTION_RANDOM_ND generates a random direction vector in ND.
! !DODEC_SHAPE_3D describes an dodecahedron in 3D.
! !DOT_2D computes the dot product of a pair of vectors in 2D.
! !DOT_3D computes the dot product of a pair of vectors in 3D.
! !DOT_ND computes the dot product of a pair of vectors in ND.
! !DOT0_2D computes the dot product of (P1-P0) and (P2-P0) in 2D.
! !DOT0_3D computes the dot product of (P1-P0) and (P2-P0) in 3D.
! !DUAL_SHAPE_3D constructs the dual of a shape in 3D.
! !ELLIPSE_AREA_2D returns the area of an ellipse in 2D.
! !ELLIPSE_POINTS_2D returns N points on an ellipse in 2D.
! !ELLIPSE_POINTS_ARC_2D returns N points on an elliptical arc in 2D.
! !ENORM_2D computes the Euclidean norm of a vector in 2D.
! !ENORM_3D computes the Euclidean norm of a vector in 3D.
! !ENORM_ND computes the Euclidean norm of a vector in ND.
! !ENORM0_2D computes the Euclidean norm of (P1-P0) in 2D.
! !ENORM0_3D computes the Euclidean norm of (P1-P0) in 3D.
! !ENORM0_ND computes the Euclidean norm of (P1-P0) in ND.
! !ENORMSQ0_2D computes the square of the Euclidean norm of (P1-P0) in 2D.
! !ENORMSQ0_3D computes the square of the Euclidean norm of (P1-P0) in 3D.
! !ENORMSQ0_ND computes the squared Euclidean norm of (P1-P0) in ND.
! !GLOB2LOC_3D converts from a global to a local coordinate system in 3D.
! !HALFSPACE_IMP_TRIANGLE_INT_3D: intersection ( implicit halfspace, triangle ) in 3D.
! !HALFSPACE_NORM_TRIANGLE_INT_3D: intersection ( normal halfspace, triangle ) in 3D.
! !HALFSPACE_TRIANGLE_INT_3D: intersection ( halfspace, triangle ) in 3D.
! !HELIX_SHAPE_3D computes points on a helix in 3D.
! !HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
! !HEXAGON_UNIT_AREA_2D returns the area of a unit regular hexagon in 2D.
! !HEXAGON_SHAPE_2D returns points on the unit hexagon in 2D.
! !I_RANDOM returns a random integer in a given range.
! !I_SWAP switches two integer values.
! !ICOL_FIND_ITEM searches a table by columns for a given value.
! !ICOL_FIND_ITEM wrap searches a table by columns for a pair of items.
! !ICOS_SHAPE_3D describes an icosahedron in 3D.
! !LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
! !LINE_EXP_POINT_DIST_3D: distance ( explicit line, point ) in 3D.
! !LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( explicit line, point ) in 2D.
! !LINE_EXP_POINT_NEAR_2D computes the point on an explicit line nearest a point in 2D.
! !LINE_EXP_POINT_NEAR_3D: nearest point on explicit line to point in 3D.
! !LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
! !LINE_EXP2PAR_2D converts a line from explicit to parametric form in 2D.
! !LINE_EXP2PAR_3D converts a line from explicit to parametric form in 3D.
! !LINE_IMP2PAR_2D converts an implicit line to parametric form in 2D.
! !LINE_IMP_POINT_DIST_2D: distance ( implicit line, point ) in 2D.
! !LINE_IMP_POINT_DIST_3D: distance ( implicit line, point ) in 3D.
! !LINE_IMP_POINT_DIST_SIGNED_2D: signed distance ( implicit line, point ) in 2D.
! !LINE_PAR2IMP_2D converts a parametric line to implicit form in 2D.
! !LINE_PAR_POINT_DIST_2D: distance ( parametric line, point ) in 2D.
! !LINE_PAR_POINT_DIST_3D: distance ( parametric line, point ) in 3D.
! !LINE_SEG_CONTAINS_POINT_1D reports if a line  contains a point in 1D.
! !LINE_SEG_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
! !LINE_SEG_POINT_DIST_2D: distance ( line segment, point ) in 2D.
! !LINE_SEG_POINT_DIST_3D: distance ( line segment, point ) in 3D.
! !LINE_SEG_POINT_NEAR_2D finds the point on a line segment nearest a point in 2D.
! !LINE_SEG_POINT_NEAR_3D finds the point on a line segment nearest a point in 3D.
! !LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
! !LINES_EXP_ANGLE_ND returns the angle between two explicit lines in ND.
! !LINES_EXP_DIST_3D computes the distance between two explicit lines in 3D.
! !LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
! !LINES_IMP_ANGLE_2D finds the angle between two implicit lines in 2D.
! !LINES_IMP_DIST_2D determines the distance between two implicit lines in 2D.
! !LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
! !LINES_PAR_ANGLE_2D finds the angle between two parametric lines in 2D.
! !LINES_PAR_ANGLE_3D finds the angle between two parametric lines in 3D.
! !LINES_PAR_DIST_3D finds the distance between two parametric lines in 3D.
! !LINES_PAR_INT_2D determines where two parametric lines intersect in 2D.
! !LINES_SEG_DIST_2D computes the distance between two line segments in 2D.
! !LINES_SEG_DIST_3D computes the distance between two line segments in 3D.
! !LINES_SEG_INT_1D computes the intersection of two line segments in 1D.
! !LINES_SEG_INT_2D computes the intersection of two line segments in 2D.
! !LOC2GLOB_3D converts from a local to global coordinate system in 3D.
! !MINABS finds a local minimum of F(X) = A * abs ( X ) + B.
! !MINQUAD finds a local minimum of F(X) = A * X * X + B * X + C.
! !NORMAL_01_SAMPLE samples the standard Normal PDF.
! !NORP2L_2D finds two points on a line normal to a given line in 2D.
! !OCTAHEDRON_SHAPE_3D describes an octahedron in 3D.
! !PARA_CONTAINS_POINT_2D determines if a point is inside a parallelogram in 2D.
! !PARA_CONTAINS_POINT_3D determines if a point is inside a parallelogram in 3D.
! !PARA_POINT_DIST_3D: distance ( parallelogram, point ) in 3D.
! !PARABOLA_EX finds the extremal point of a parabola determined by three points.
! !PARABOLA_EX2 finds the extremal point of a parabola determined by three points.
! !PI returns the value of pi.
! !PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
! !PLANE_EXP_POINT_DIST_3D: distance ( explicit plane, point ) in 3D.
! !PLANE_EXP_PROJECT_3D projects points through a point onto a plane in 3D.
! !PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
! !PLANE_EXP2NORM_3D converts an explicit plane to normal form in 3D.
! !PLANE_GRID_3D computes points and lines making up a planar grid in 3D.
! !PLANE_IMP_LINE_PAR_INT_3D: intersection ( implicit plane, parametric line ) in 3D.
! !PLANE_IMP_LINE_SEG_NEAR_3D: nearest ( implicit plane, line segment ) in 3D.
! !PLANE_IMP_POINT_DIST_3D: distance ( implicit plane, point ) in 3D.
! !PLANE_IMP_POINT_DIST_SIGNED_3D: signed distance ( implicit plane, point) in 3D.
! !PLANE_IMP_POINT_NEAR_3D: nearest point on a implicit plane to a point in 3D.
! !PLANE_IMP_TRIANGLE_INT_3D: intersection ( implicit plane, triangle ) in 3D.
! !PLANE_IMP_TRIANGLE_NEAR_3D: nearest ( implicit plane, triangle ) in 3D.
! !PLANE_IMP2EXP_3D converts an implicit plane to explicit form in 3D.
! !PLANE_IMP2NORM_3D converts an implicit plane to normal form in 3D.
! !PLANE_NORM_BASIS_3D finds two perpendicular vectors in a plane in 3D.
! !PLANE_NORM_TRIANGLE_INT_3D: intersection ( normal plane, triangle ) in 3D.
! !PLANE_NORM2EXP_3D converts a normal plane to explicit form in 3D.
! !PLANE_NORM2IMP_3D converts a normal form plane to implicit form in 3D.
! !PLANE_TRIANGLE_INT_ADD_3D is a utility for plane/triangle intersections.
! !POINTS_AVOID_POINT_NAIVE_2D determines if a point is "far enough" from a set of points in 2D.
! !POINTS_BISECT_LINE_IMP_2D finds the implicit line bisecting the line between two points in 2D.
! !POINTS_BISECT_LINE_PAR_2D finds the parametric line bisecting the line between two points in 2D.
! !POINTS_CENTROID_2D computes the discrete centroid of a point set in 2D.
! !POINTS_COLIN_2D estimates the colinearity of 3 points in 2D.
! !POINTS_COLIN_3D estimates the colinearity of 3 points in 3D.
! !POINTS_DELAUNAY_NAIVE_2D is a naive Delaunay triangulation scheme.
! !POINTS_DIST_2D finds the distance between two points in 2D.
! !POINTS_DIST_3D finds the distance between two points in 3D.
! !POINTS_DIST_ND finds the distance between two points in ND.
! !POINTS_HULL_2D computes the convex hull of a set of points in 2D.
! !POINTS_NEAREST_POINT_NAIVE_ND finds the nearest point to a given point in ND.
! !POLYGON_1_2D integrates the function 1 over a polygon in 2D.
! !POLYGON_X_2D integrates the function X over a polygon in 2D.
! !POLYGON_Y_2D integrates the function Y over a polygon in 2D.
! !POLYGON_XX_2D integrates the function X*X over a polygon in 2D.
! !POLYGON_XY_2D integrates the function X*Y over a polygon in 2D.
! !POLYGON_YY_2D integrates the function Y*Y over a polygon in 2D.
! !POLYGON_AREA_2D computes the area of a polygon in 2D.
! !POLYGON_AREA_2_2D computes the area of a polygon in 2D.
! !POLYGON_AREA_3D computes the area of a polygon in 3D.
! !POLYGON_AREA_2_3D computes the area of a polygon in 3D.
! !POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
! !POLYGON_CENTROID_2_2D computes the centroid of a polygon in 2D.
! !POLYGON_CENTROID_3D computes the centroid of a polygon in 3D.
! !POLYGON_CONTAINS_POINT_2D finds if a point is inside a simple polygon in 2D.
! !POLYGON_CONTAINS_POINT_2_2D finds if a point is inside a convex polygon in 2D.
! !POLYGON_CONVEX determines whether a polygon is convex in 2D.
! !POLYGON_INRAD_DATA_2D determines polygonal data from its inner radius in 2D.
! !POLYGON_OUTRAD_DATA_2D determines polygonal data from its outer radius in 2D.
! !POLYGON_SIDE_DATA_2D determines polygonal data from its side length in 2D.
! !POLYHEDRON_SURFACE_3D computes the surface area of a polyhedron in 3D.
! !POLYHEDRON_VOLUME_3D computes the volume of a polyhedron in 3D.
! !POLYLINE_INDEX_POINT_ND evaluates a polyline at a given arclength in ND.
! !POLYLINE_LENGTH_ND computes the length of a polyline in ND.
! !PROPLANE2 produces 2D coordinates of points that lie in a plane, in 3D.
! !PROPLANE3 projects points orthographically onto a plane, in 3D.
! !PROVEC projects a vector from M space into N space.
! !PYRAMID_VOLUME_3D computes the volume of a pyramid with square base in 3D.
! !QUAD_AREA_2D computes the area of a quadrilateral in 2D.
! !QUAD_CONTAINS_POINT_2D finds if a point is inside a convex quadrilateral in 2D.
! !QUAD_POINT_DIST_2D finds the distance from a point to a quadrilateral in 2D.
! !QUAD_POINT_DIST_SIGNED_2D: signed distanct ( quadrilateral, point ) in 2D.
! !QUAT_CONJ conjugates a quaternion.
! !QUAT_INV inverts a quaternion.
! !QUAT_MUL multiplies two quaternions.
! !QUAT_NORM computes the norm of a quaternion.
! !R_MODP returns the nonnegative remainder of real division.
! !R_RANDOM returns a random real in a given range.
! !R_SWAP switches two real values.
! !RADEC_SEPARATION: separation between right ascension-declination points.
! !RADEC_TO_XYZ converts right ascension declination to (X,Y,Z) coordinates
! !RADIANS_TO_DEGREES converts an angle from radians to degrees.
! !RMAT2_DET computes the determinant of a 2 by 2 matrix.
! !RMAT2_INVERSE inverts a 2 by 2 real matrix using Cramer's rule.
! !RMAT3_DET computes the determinant of a 3 by 3 matrix.
! !RMAT3_INVERSE inverts a 3 by 3 real matrix using Cramer's rule.
! !RMAT4_DET computes the determinant of a 4 by 4 matrix.
! !RMAT5_DET computes the determinant of a 5 by 5 matrix.
! !RMAT_PRINT prints a real matrix.
! !RMAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
! !ROTATION_AXIS_VECTOR_3D rotates a vector around an axis vector in 3D.
! !ROTATION_AXIS2MAT_3D converts a rotation from axis to matrix format in 3D.
! !ROTATION_AXIS2QUAT_3D converts a rotation from axis to quaternion format in 3D.
! !ROTATION_MAT_VECTOR applies a marix rotation to a vector in 3d.
! !ROTATION_MAT2AXIS_3D converts a rotation from matrix to axis format in 3D.
! !ROTATION_MAT2QUAT_3D converts a rotation from matrix to quaternion format in 3D.
! !ROTATION_QUAT_VECTOR applies a quaternion rotation to a vector in 3d.
! !ROTATION_QUAT2AXIS_3D converts a rotation from quaternion to axis format in 3D.
! !ROTATION_QUAT2MAT_3D converts a rotation from quaternion to matrix format in 3D.
! !RVEC_PRINT_2D prints a 2D vector.
! !RVEC_PRINT_3D prints a 3D vector.
! !RVEC2_PRINT prints a pair of real vectors.
! !RVEC3_PRINT prints a trio of real vectors.
! !SHAPE_POINT_DIST_2D: distance ( regular shape, point ) in 2D.
! !SHAPE_POINT_NEAR_2D: nearest point ( regular shape, point ) in 2D.
! !SHAPE_PRINT_3D prints information about a polyhedron in 3D.
! !SHAPE_RAY_INT_2D: intersection ( regular shape, ray ) in 2D.
! !SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
! !SPHERE_DIA2IMP_3D converts a diameter to an implicit sphere in 3D.
! !SPHERE_EXP_CONTAINS_POINT_3D determines if an explicit sphere contains a point in 3D.
! !SPHERE_EXP_NEAR_POINT_3D finds the nearest point on an explicit sphere to a point in 3D.
! !SPHERE_EXP2IMP_3D converts a sphere from explicit to implicit form in 3D.
! !SPHERE_IMP_AREA_3D computes the surface area of an implicit sphere in 3D.
! !SPHERE_IMP_AREA_ND computes the surface area of an implicit sphere in ND.
! !SPHERE_IMP_CONTAINS_POINT_3D determines if an implicit sphere contains a point in 3D.
! !SPHERE_IMP_GRIDFACES_3D produces a grid of triangles on an implicit sphere in 3D.
! !SPHERE_IMP_GRIDLINES_3D produces "grid lines" on an implicit sphere in 3D.
! !SPHERE_IMP_GRIDPOINTS_3D produces "grid points" on an implicit sphere in 3D.
! !SPHERE_IMP_LINE_PROJECT_3D projects a line onto an implicit sphere in 3D.
! !SPHERE_IMP_LOCAL2XYZ_3D converts local to XYZ coordinates on an implicit sphere in 3D.
! !SPHERE_IMP_NEAR_POINT_3D finds the nearest point on an implicit sphere to a point in 3D.
! !SPHERE_IMP_POINT_PROJECT_3D projects a point onto an implicit sphere, in 3D.
! !SPHERE_IMP_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
! !SPHERE_IMP_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
! !SPHERE_IMP_VOLUME_3D computes the volume of an implicit sphere in 3D.
! !SPHERE_IMP_VOLUME_ND computes the volume of an implicit sphere in ND.
! !STRING_2D groups line segments into connected lines in 2D.
! !TAND returns the tangent of an angle given in degrees.
! !TETRA_CENTROID_3D computes the centroid of a tetrahedron in 3D.
! !TETRA_CONTAINS_POINT_3D finds if a point is inside a tetrahedron in 3D.
! !TETRA_OUTSPHERE_3D computes the exscribed sphere of a tetrahedron in 3D.
! !TETRA_SHAPE_3D describes a tetrahedron in 3D.
! !TETRA_VOLUME_3D computes the volume of a tetrahedron in 3D.
! !TMAT_INIT initializes the geometric transformation matrix.
! !TMAT_MXM multiplies two geometric transformation matrices.
! !TMAT_MXP multiplies a geometric transformation matrix times a point.
! !TMAT_MXP2 multiplies a geometric transformation matrix times N points.
! !TMAT_MXV multiplies a geometric transformation matrix times a vector.
! !TMAT_ROT_AXIS applies a coordinate axis rotation to the geometric transformation matrix.
! !TMAT_ROT_VECTOR applies an arbitrary axis rotation to the geometric transformation matrix.
! !TMAT_SCALE applies a scaling to the geometric transformation matrix.
! !TMAT_SHEAR applies a shear to the geometric transformation matrix.
! !TMAT_TRANS applies a translation to the geometric transformation matrix.
! !TORUS_AREA_3D returns the area of a torus in 3D.
! !TORUS_VOLUME_3D computes the volume of a torus in 3D.
! !TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
! !TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
! !TRIANGLE_AREA_2D computes the area of a triangle in 2D.
! !TRIANGLE_AREA_SIGNED_2D computes the signed area of a triangle in 2D.
! !TRIANGLE_AREA_3D computes the area of a triangle in 3D.
! !TRIANGLE_AREA_2_3D computes the area of a triangle in 3D.
! !TRIANGLE_CENTROID_2D computes the centroid of a triangle in 2D.
! !TRIANGLE_CENTROID_3D computes the centroid of a triangle in 3D.
! !TRIANGLE_CONTAINS_POINT_2D finds if a point is inside a triangle in 2D.
! !TRIANGLE_DIAMETER_2D computes the diameter of a triangle in 2D.
! !TRIANGLE_GRIDPOINTS_2D computes gridpoints within a triangle in 2D.
! !TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
! !TRIANGLE_LINE_IMP_INT_2D finds where an implicit line intersects a triangle in 2D.
! !TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
! !TRIANGLE_OUTCIRCLE_2D computes the exscribed circle of a triangle in 2D.
! !TRIANGLE_POINT_DIST_2D: distance ( triangle, point ) in 2D.
! !TRIANGLE_POINT_DIST_3D: distance ( triangle, point ) in 3D.
! !TRIANGLE_POINT_DIST_SIGNED_2D: signed distance ( triangle, point ) in 2D.
! !TRIANGLE_POINT_NEAR_2D computes the nearest triangle point to a point in 2D.
! !TUBE_2D constructs a "tube" of given width around a path in 2D.
! !VECTOR_DIRECTIONS_2D returns the direction angles of a vector in 2D.
! !VECTOR_DIRECTIONS_3D returns the direction angles of a vector in 3D.
! !VECTOR_ROTATE_2D rotates a vector around the origin in 2D.
! !VECTOR_ROTATE_BASE_2D rotates a vector around a base point in 2D.
! !VECTOR_SEPARATION_3D finds the angular separation between vectors in 3D.
! !VECTOR_SEPARATION_ND finds the angular separation between vectors in ND.
! !VECTOR_UNIT_2D normalizes a vector in 2D.
! !VECTOR_UNIT_3D normalizes a vector in 3D.
! !VECTOR_UNIT_ND normalizes a vector in ND.
! !VOXEL_LINE_3D computes the voxels that line on a line in 3D.
! !VOXEL_REGION_3D arranges a set of voxels into contiguous regions in 3D.
! !VOXEL_STEP_3D computes voxels along a line from a given point in 3D.
! !XYZ TO RADEC converts (XYZ) to right ascension declination coordinates
!
!
!
!
! function arc_cosine ( c )
! real arc_cosine
! !
! !*******************************************************************************
! !
! !! ARC_COSINE computes the arc cosine function, with argument truncation.
! !
! !
! !  Discussion:
! !
! !    If you call your system ACOS routine with an input argument that is
! !    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
! !    This routine truncates arguments outside the range.
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real C, the argument.
! !
! !    Output, real ARC_COSINE, an angle whose cosine is C.
! !
!   real c
!   real c2
! !
!   c2 = c
!   c2 = max ( c2, -1.0E+00 )
!   c2 = min ( c2, +1.0E+00 )
!
!   arc_cosine = acos ( c2 )
!
!   return
! end function
!
! subroutine angle_contains_ray_2d ( inside, x1, y1, x2, y2, x3, y3, x, y )
! !
! !*******************************************************************************
! !
! !! ANGLE_CONTAINS_RAY_2D determines if an angle contains a ray, in 2D.
! !
! !
! !  Discussion:
! !
! !    The angle is defined by the sequence of points (X1,Y1), (X2,Y2)
! !    and (X3,Y3).
! !
! !    The ray is defined by the sequence of points (X2,Y2), (X,Y).
! !
! !  Modified:
! !
! !    17 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates of
! !    the angle.
! !
! !    Input, real X, Y, the end point of the ray to be checked.
! !    The ray is assumed to have an origin at (X2,Y2).
! !
! !    Output, logical INSIDE, is .TRUE. if the ray is inside
! !    the angle or on its boundary, and .FALSE. otherwise.
! !
!   real a1
!   real a2
!   logical inside
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!   a1 = angle_deg_2d ( x1, y1, x2, y2, x, y )
!   a2 = angle_deg_2d ( x1, y1, x2, y2, x3, y3 )
!
!   if ( a1 <= a2 ) then
!     inside = .true.
!   else
!     inside = .false.
!   end if
!
!   return
! end subroutine
!
! function angle_deg_2d ( x1, y1, x2, y2, x3, y3 )
! !
! !*******************************************************************************
! !
! !! ANGLE_DEG_2D returns the angle swept out between two rays in 2D.
! !
! !
! !  Discussion:
! !
! !    Except for the zero angle case, it should be true that
! !
! !      ANGLE_DEG_2D(X1,Y1,X2,Y2,X3,Y3)
! !    + ANGLE_DEG_2D(X3,Y3,X2,Y2,X1,Y1) = 360.0
! !
! !  Modified:
! !
! !    14 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, define the rays
! !    ( X1-X2, Y1-Y2 ) and ( X3-X2, Y3-Y2 ) which in turn define the
! !    angle, counterclockwise from ( X1-X2, Y1-Y2 ).
! !
! !    Output, real ANGLE_DEG_2D, the angle swept out by the rays, measured
! !    in degrees.  0 <= ANGLE_DEG_2D < 360.  If either ray has zero length,
! !    then ANGLE_DEG_2D is set to 0.
! !
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
!   real angle_rad_2d
!   real angle_deg_2d
! !
!   x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
!   y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )
!
!   if ( x == 0.0E+00 .and. y == 0.0E+00 ) then
!
!    angle_deg_2d = 0.0E+00
!
!   else
!
!   if( x == 0.0D+00) then
!    if(y>0.d0) then
!      angle_rad_2d =        pi/2.d0
!    else
!      angle_rad_2d =  3.d0* pi/2.d0
!    endif
!   else
!    angle_rad_2d = atan2 ( y, x )
!   endif
!
!   if ( angle_rad_2d < 0.0E+00 ) then
!     angle_rad_2d = angle_rad_2d + 2.0E+00 * dacos(-1.d0)
!   end if
!
!   angle_deg_2d = radians_to_degrees ( angle_rad_2d )
!
!   end if
!
!   return
! end function
!
!
!
! function angle_rad_2d ( x1, y1, x2, y2, x3, y3 )
! !
! !*******************************************************************************
! !
! !! ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
! !
! !
! !  Discussion:
! !
! !    Except for the zero angle case, it should be true that
! !
! !      ANGLE_RAD_2D(X1,Y1,X2,Y2,X3,Y3)
! !    + ANGLE_RAD_2D(X3,Y3,X2,Y2,X1,Y1) = 2 * PI
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, define the rays
! !    ( X1-X2, Y1-Y2 ) and ( X3-X2, Y3-Y2 ) which in turn define the
! !    angle, counterclockwise from ( X1-X2, Y1-Y2 ).
! !
! !    Output, real ANGLE_RAD_2D, the angle swept out by the rays,
! !    in radians.  0 <= ANGLE_RAD_2D < 2 PI.  If either ray has zero
! !    length, then ANGLE_RAD_2D is set to 0.
! !
!   real angle_rad_2d
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!   x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
!   y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )
!
!   if ( x == 0.0E+00 .and. y == 0.0E+00 ) then
!
!     angle_rad_2d = 0.0E+00
!
!   else
!
!   if( x == 0.0D+00) then
!    if(y>0.d0) then
!      angle_rad_2d =        pi/2.d0
!    else
!      angle_rad_2d =  3.d0* pi/2.d0
!    endif
!   else
!    angle_rad_2d = atan2 ( y, x )
!   endif
!
!     if ( angle_rad_2d < 0.0E+00 ) then
!       angle_rad_2d = angle_rad_2d + 2.0E+00 * pi
!     end if
!
!   end if
!
!   return
! end function
!
!
! function angle_rad_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
! real angle_rad_3d
! !
! !*******************************************************************************
! !
! !! ANGLE_RAD_3D returns the angle in radians between two rays in 3D.
! !
! !
! !  Discussion:
! !
! !    The routine always computes the SMALLER of the two angles between
! !    two rays.  Thus, if the rays make an (exterior) angle of
! !    1.5 radians, the (interior) angle of 0.5 radians will be reported.
! !
! !  Formula:
! !
! !    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are three points
! !    which define the rays.  The rays are:
! !    ( X1-X2, Y1-Y2, Z1-Z2 ) and ( X3-X2, Y3-Y2, Z3-Z2 ).
! !
! !    Output, real ANGLE_RAD_3D, the angle between the two rays, in radians.
! !    This value will always be between 0 and PI.  If either ray has
! !    zero length, then the angle is returned as zero.
! !
!   real theta
!   real v1norm
!   real v2norm
!   real x1,dot
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
!   real z1
!   real z2
!   real z3
! !
!   v1norm = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
!   v2norm = enorm0_3d ( x3, y3, z3, x2, y2, z2 )
!
!   if ( v1norm == 0.0E+00 .or. v2norm == 0.0E+00 ) then
!     angle_rad_3d = 0.0E+00
!   else
!     dot = dot0_3d ( x2, y2, y2, x1, y1, z1, x3, y3, z3 )
!     angle_rad_3d = arc_cosine ( dot / ( v1norm * v2norm ) )
!   end if
!
!   return
! end function
! function angle_rad_nd ( n, vec1, vec2 )
! real angle_rad_nd
! !
! !*******************************************************************************
! !
! !! ANGLE_RAD_ND returns the angle in radians between two rays in ND.
! !
! !
! !  Discussion:
! !
! !    This routine always computes the SMALLER of the two angles between
! !    two rays.  Thus, if the rays make an (exterior) angle of 1.5 PI,
! !    then the (interior) angle of 0.5 PI is reported.
! !
! !  Formula:
! !
! !    X dot Y = Norm(X) * Norm(Y) * Cos( Angle(X,Y) )
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of entries in the rays.
! !
! !    Input, real VEC1(N), VEC2(N), the two rays to be considered.
! !
! !    Output, real ANGLE_RAD_ND, the angle between the rays, in radians.
! !    This value will always be between 0 and PI.
! !
!   integer n
! !
!   real theta
!   real v1norm
!   real v2norm
!   real vec1(n)
!   real vec2(n),dot
! !
!   dot = dot_nd ( n, vec1, vec2 )
!
!   v1norm = enorm_nd ( n, vec1 )
!
!   v2norm = enorm_nd ( n, vec2 )
!
!   if ( v1norm == 0.0E+00 .or. v2norm == 0.0E+00 ) then
!     angle_rad_nd = 0.0E+00
!   else
!     angle_rad_nd = arc_cosine ( dot / ( v1norm * v2norm ) )
!   end if
!
!   return
! end function
!
! function anglei_deg_2d ( x1, y1, x2, y2, x3, y3 )
! real anglei_deg_2d
! !
! !*******************************************************************************
! !
! !! ANGLEI_DEG_2D returns the interior angle in degrees between two rays in 2D.
! !
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, define the rays
! !    (X1-X2,Y1-Y2) and (X3-X2,Y3-Y2) which in turn define the angle.
! !
! !    Output, real ANGLEI_DEG_2D, the angle swept out by the rays, measured
! !    in degrees.  This value satisfies 0 <= ANGLEI_DEG_2D < 180.0.  If either
! !    ray is of zero length, then ANGLEI_DEG_2D is returned as 0.
! !
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
!   real anglei_rad_2d
! !
!   x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
!   y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )
!
!   if ( x == 0.0E+00 .and. y == 0.0E+00 ) then
!
!     anglei_deg_2d = 0.0E+00
!
!   else
!
!     anglei_rad_2d = atan2 ( y, x )
!
!     if ( anglei_rad_2d < 0.0E+00 ) then
!       anglei_rad_2d = anglei_rad_2d + 2.0E+00 * pi
!     end if
!
!     anglei_deg_2d = radians_to_degrees ( anglei_rad_2d )
!
!     if ( anglei_deg_2d > 180.0E+00 ) then
!       anglei_deg_2d = 360.0E+00 - anglei_deg_2d
!     end if
!
!   end if
!
!   return
! end function
!
! function anglei_rad_2d ( x1, y1, x2, y2, x3, y3 )
! !
! !*******************************************************************************
! !
! !! ANGLEI_RAD_2D returns the interior angle in radians between two rays in 2D.
! !
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, define the rays
! !    (X1-X2,Y1-Y2) and (X3-X2,Y3-Y2) which in turn define the angle.
! !
! !    Output, real ANGLEI_RAD_2D, the angle swept out by the rays, measured
! !    in radians.  This value satisfies 0 <= ANGLEI_RAD_2D < PI.  If either
! !    ray is of zero length, then ANGLEI_RAD_2D is returned as 0.
! !
!   real anglei_rad_2d
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!   x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
!   y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )
!
!   if ( x == 0.0E+00 .and. y == 0.0E+00 ) then
!
!     anglei_rad_2d = 0.0E+00
!
!   else
!
!     anglei_rad_2d = atan2 ( y, x )
!
!     if ( anglei_rad_2d < 0.0E+00 ) then
!       anglei_rad_2d = anglei_rad_2d + 2.0E+00 * pi
!     end if
!
!   end if
!
!   return
! end function
!
! function atan4( y, x )
! real atan4
! !
! !*******************************************************************************
! !
! !! ATAN4 computes the inverse tangent of the ratio Y / X.
! !
! !
! !  Discussion:
! !
! !    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
! !    the built in functions ATAN and ATAN2 already do.
! !
! !    However:
! !
! !    * ATAN4 always returns a positive angle, between 0 and 2 PI,
! !      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
! !      and [-PI,+PI] respectively;
! !
! !    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
! !     function by contrast always returns an angle in the first or fourth
! !     quadrants.
! !
! !  Modified:
! !
! !    14 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real Y, X, two quantities which represent the tangent of
! !    an angle.  If Y is not zero, then the tangent is (Y/X).
! !
! !    Output, real ATAN4, an angle between 0 and 2 * PI, whose tangent is
! !    (Y/X), and which lies in the appropriate quadrant so that the signs
! !    of its cosine and sine match those of X and Y.
! !
!   real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510E+00
! !
!   real abs_x
!   real abs_y
!   real theta
!   real theta_0
!   real x
!   real y
! !
! !  Special cases:
! !
!   if ( x == 0.0E+00 ) then
!
!     if ( y > 0.0E+00 ) then
!       theta = PI / 2.0E+00
!     else if ( y < 0.0E+00 ) then
!       theta = 3.0E+00 * PI / 2.0E+00
!     else if ( y == 0.0E+00 ) then
!       theta = 0.0E+00
!     end if
!
!   else if ( y == 0.0E+00 ) then
!
!     if ( x > 0.0E+00 ) then
!       theta = 0.0E+00
!     else if ( x < 0.0E+00 ) then
!       theta = PI
!     end if
! !
! !  We assume that ATAN2 is correct when both arguments are positive.
! !
!   else
!
!     abs_y = abs ( y )
!     abs_x = abs ( x )
!
!     theta_0 = atan2 ( abs_y, abs_x )
!
!     if ( x > 0.0E+00 .and. y > 0.0E+00 ) then
!       theta = theta_0
!     else if ( x < 0.0E+00 .and. y > 0.0E+00 ) then
!       theta = PI - theta_0
!     else if ( x < 0.0E+00 .and. y < 0.0E+00 ) then
!       theta = PI + theta_0
!     else if ( x > 0.0E+00 .and. y < 0.0E+00 ) then
!       theta = 2.0E+00 * PI - theta_0
!     end if
!
!   end if
!
!   atan4 = theta
!
!   return
!   end function
!
! subroutine basis_map_3d ( u1, u2, u3, v1, v2, v3, a, ierror )
! !
! !*******************************************************************************
! !
! !! BASIS_MAP_3D computes the matrix which maps one basis to another.
! !
! !
! !  Discussion:
! !
! !    As long as the vectors U1, U2 and U3 are linearly independent,
! !    a matrix A will be computed that maps U1 to V1, U2 to V2, and
! !    U3 to V3.
! !
! !    Depending on the values of the vectors, A may represent a
! !    rotation, reflection, dilation, project, or a combination of these
! !    basic linear transformations.
! !
! !  Modified:
! !
! !    20 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real U1(3), U2(3), U3(3), the three "domain" or "preaimage"
! !    vectors, which should be linearly independent.
! !
! !    Input, real V1(3), V2(3), V3(3), the three "range" or "aimage" vectors.
! !
! !    Output, real A(3,3), a matrix with the property that A * U1 = V1,
! !    A * U2 = V2 and A * U3 = V3.
! !
! !    Output, integer IERROR, error flag.
! !    0, no error occurred.
! !    nonzero, the matrix [ U1 | U2 | U3 ] is exactly singular.
! !
!   real a(3,3)
!   real b(3,3)
!   real c(3,3)
!   real det
!   integer i
!   integer ierror
!   integer j
!   integer k
!   real u1(3)
!   real u2(3)
!   real u3(3)
!   real v1(3)
!   real v2(3)
!   real v3(3)
! !
!   ierror = 0
! !
! !  Set B = [ U1 | U2 | U3 ].
! !
!   b(1:3,1) = u1(1:3)
!   b(1:3,2) = u2(1:3)
!   b(1:3,3) = u3(1:3)
! !
! !  Compute C = the inverse of [ U1 | U2 | U3 ].
! !
!   call rmat3_inverse ( b, c, det )
!
!   if ( det == 0.0E+00 ) then
!     ierror = 1
!     return
!   end if
! !
! !  Set B = [ V1 | V2 | V3 ].
! !
!   b(1:3,1) = v1(1:3)
!   b(1:3,2) = v2(1:3)
!   b(1:3,3) = v3(1:3)
! !
! !  A = [ V1 | V2 | V3 ] * inverse [ U1 | U2 | U3 ].
! !
!   a(1:3,1:3) = matmul ( b(1:3,1:3), c(1:3,1:3) )
!
!   return
! end subroutine
!
! function point_inside_parallelipiped_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, x, y, z )
! !
! !*******************************************************************************
! !
! !! POINT_INSIDE_PARALLELIPIPED_3D determines if a point is inside a
! !parallelepiped in 3D.
! !
! !
! !  Definition:
! !
! !    A parallelepiped is a "slanted box", that is, opposite
! !    sides are parallel planes.
! !
! !  Diagram:
! !
! !                   *------------------*
! !                  / \                / \
! !                 /   \              /   \
! !                /     \            /     \
! !          (X4,Y4,Z4)--------------*       \
! !                \        .         \       \
! !                 \        .         \       \
! !                  \        .         \       \
! !                   \   (X2,Y2,Z2).....\-------\
! !                    \     /            \     /
! !                     \   /              \   /
! !                      \ /                \ /
! !                (X1,Y1,Z1)-----------(X3,Y3,Z3)
! !
! !  Modified:
! !
! !    04 February 199
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
! !    coordinates of four corners of the parallelepiped, which we will
! !    call P1, P2, P3 and P4.  It is assumed that P2, P3 and P4 are
! !    immediate neighbors of P1.
! !
! !    Input, real X, Y, Z, the point to be checked.
! !
! !    Output, logical POINT_INSIDE_PARALLELIPIPED_3D, is .TRUE. if (X,Y,Z)
! !    is inside the parallelepiped, or on its boundary, and .FALSE. otherwise.
! !
!   implicit none
! !
!   real p21dot
!   real p21normsq
!   real p31dot
!   real p31normsq
!   real p41dot
!   real p41normsq
!   logical point_inside_parallelipiped_3d
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
!   real z
!   real z1
!   real z2
!   real z3
!   real z4
! !
!   p21normsq = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )
!   p31normsq = enormsq0_3d ( x1, y1, z1, x3, y3, z3 )
!   p41normsq = enormsq0_3d ( x1, y1, z1, x4, y4, z4 )
!
!   p21dot = dot0_3d ( x1, y1, z1, x2, y2, z2, x, y, z )
!   p31dot = dot0_3d ( x1, y1, z1, x3, y3, z3, x, y, z )
!   p41dot = dot0_3d ( x1, y1, z1, x4, y4, z4, x, y, z )
!
!   if ( 0.0E+00 <= p21dot .and. p21dot <= p21normsq .and. &
!        0.0E+00 <= p31dot .and. p31dot <= p31normsq .and. &
!        0.0E+00 <= p41dot .and. p41dot <= p41normsq ) then
!     point_inside_parallelipiped_3d = .true.
!   else
!     point_inside_parallelipiped_3d = .false.
!   end if
!
!   return
! end function
!
! subroutine box_contains_point_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, x, y, z, inside )
! !
! !*******************************************************************************
! !
! !! BOX_CONTAINS_POINT_3D determines if a point is inside a parallelepiped in 3D.
! !
! !
! !  Definition:
! !
! !    A parallelepiped is a "slanted box", that is, opposite
! !    sides are parallel planes.
! !
! !  Diagram:
! !
! !                   *------------------*
! !                  / \                / \
! !                 /   \              /   \
! !                /     \            /     \
! !          (X4,Y4,Z4)--------------*       \
! !                \        .         \       \
! !                 \        .         \       \
! !                  \        .         \       \
! !                   \   (X2,Y2,Z2).....\-------\
! !                    \     /            \     /
! !                     \   /              \   /
! !                      \ /                \ /
! !                (X1,Y1,Z1)-----------(X3,Y3,Z3)
! !
! !  Modified:
! !
! !    04 February 199
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
! !    coordinates of four corners of the parallelepiped, which we will
! !    call P1, P2, P3 and P4.  It is assumed that P2, P3 and P4 are
! !    immediate neighbors of P1.
! !
! !    Input, real X, Y, Z, the point to be checked.
! !
! !    Output, logical INSIDE, is .TRUE. if (X,Y,Z) is inside the
! !    parallelepiped, or on its boundary, and .FALSE. otherwise.
! !
!   logical inside
!   real p21dot
!   real p21normsq
!   real p31dot
!   real p31normsq
!   real p41dot
!   real p41normsq
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
!   real z
!   real z1
!   real z2
!   real z3
!   real z4
! !
!   p21normsq = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )
!   p31normsq = enormsq0_3d ( x1, y1, z1, x3, y3, z3 )
!   p41normsq = enormsq0_3d ( x1, y1, z1, x4, y4, z4 )
!
!   p21dot = dot0_3d ( x1, y1, z1, x2, y2, z2, x, y, z )
!   p31dot = dot0_3d ( x1, y1, z1, x3, y3, z3, x, y, z )
!   p41dot = dot0_3d ( x1, y1, z1, x4, y4, z4, x, y, z )
!
!   if ( 0.0E+00 <= p21dot .and. p21dot <= p21normsq .and. &
!        0.0E+00 <= p31dot .and. p31dot <= p31normsq .and. &
!        0.0E+00 <= p41dot .and. p41dot <= p41normsq ) then
!     inside = .true.
!   else
!     inside = .false.
!   end if
!
!   return
! end subroutine
! subroutine box_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, x, y, z, dist )
! !
! !*******************************************************************************
! !
! !! BOX_POINT_DIST_3D: distance ( parallelepiped, point ) in 3D.
! !
! !
! !  Definition:
! !
! !    A parallelepiped is a "slanted box", that is, opposite
! !    sides are parallel planes.
! !
! !  Diagram:
! !
! !           7----8
! !          /|   /|
! !         / 3--/-5
! !        4----6 /
! !        |/   |/
! !        1----2
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, half of
! !    the corners of the box, from which the other corners can be
! !    deduced.  The corners should be chosen so that the first corner
! !    is directly connected to the other three.  The locations of
! !    corners 5, 6, 7 and 8 will be computed by the parallelogram
! !    relation.
! !
! !    Input, real X, Y, Z, the point which is to be checked.
! !
! !    Output, real DIST, the distance from the point to the box.  DIST is
! !    zero if the point lies exactly on the box.
! !
!   real dis
!   real dist
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real x5
!   real x6
!   real x7
!   real x8
!
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
!   real y5
!   real y6
!   real y7
!   real y8
!
!   real z
!   real z1
!   real z2
!   real z3
!   real z4
!   real z5
!   real z6
!   real z7
!   real z8
! !
! !  Fill in the other corners
! !
!   x5 = x2 + x3 - x1
!   y5 = y2 + y3 - y1
!   z5 = z2 + z3 - z1
!
!   x6 = x2 + x4 - x1
!   y6 = y2 + y4 - y1
!   z6 = z2 + z4 - z1
!
!   x7 = x3 + x4 - x1
!   y7 = y3 + y4 - y1
!   z7 = z3 + z4 - z1
!
!   x8 = x2 + x3 + x4 - 2.0E+00 * x1
!   y8 = y2 + y3 + y4 - 2.0E+00 * y1
!   z8 = z2 + z3 + z4 - 2.0E+00 * z1
! !
! !  Compute the distance from the point ( X, Y, Z ) to each of the six
! !  paralleogram faces.
! !
!   call para_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z, dis )
!
!   dist = dis
!
!   call para_point_dist_3d ( x1, y1, z1, x2, y2, z2, x4, y4, z4, x, y, z, dis )
!
!   dist = min ( dist, dis )
!
!   call para_point_dist_3d ( x1, y1, z1, x3, y3, z3, x4, y4, z4, x, y, z, dis )
!
!   dist = min ( dist, dis )
!
!   call para_point_dist_3d ( x8, y8, z8, x5, y5, z5, x6, y6, z6, x, y, z, dis )
!
!   dist = min ( dist, dis )
!
!   call para_point_dist_3d ( x8, y8, z8, x5, y5, z5, x7, y7, z7, x, y, z, dis )
!
!   dist = min ( dist, dis )
!
!   call para_point_dist_3d ( x8, y8, z8, x6, y6, z6, x7, y7, z7, x, y, z, dis )
!
!   dist = min ( dist, dis )
!
!   return
! end subroutine
! subroutine circle_area_2d ( r, area )
! !
! !*******************************************************************************
! !
! !! CIRCLE_AREA_2D computes the area of a circle in 2D.
! !
! !
! !  Modified:
! !
! !    12 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Output, real AREA, the area of the circle.
! !
!   real area
!   real r
! !
!   area = pi * r * r
!
!   return
! end subroutine
! subroutine circle_lune_area_2d ( r, theta, area )
! !
! !*******************************************************************************
! !
! !! CIRCLE_LUNE_AREA_2D returns the area of a circular lune in 2D.
! !
! !
! !  Discussion:
! !
! !    A lune is formed by drawing a circular arc, and joining its endpoints.
! !
! !  Modified:
! !
! !    12 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real THETA, the angle subtended by the arc.
! !
! !    Output, real AREA, the area of the lune.
! !
!   real area
!   real area_sector
!   real area_triangle
!   real r
!   real theta
! !
!   call circle_sector_area_2d ( r, theta, area_sector )
!   call circle_triangle_area_2d ( r, theta, area_triangle )
!
!   area = area_sector - area_triangle
!
!   return
! end subroutine
! subroutine circle_lune_centroid_2d ( r, xc, yc, theta1, theta2, x, y )
! !
! !*******************************************************************************
! !
! !! CIRCLE_LUNE_CENTROID_2D returns the centroid of a circular lune in 2D.
! !
! !
! !  Discussion:
! !
! !    A lune is formed by drawing a circular arc, and joining its endpoints.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    28 February 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real XC, YC, the coordinates of the center of the circle.
! !
! !    Input, real THETA1, THETA2, the angles of the first and last points
! !    on the circular arc.
! !
! !    Output, real X, Y, the coordinates of the centroid of the lune.
! !
!   real area
!   real d
!   real r
!   real theta
!   real theta1
!   real theta2
!   real x
!   real xc
!   real y
!   real yc
! !
!   theta = theta2 - theta1
!
!   if ( theta == 0.0E+00 ) then
!     d = r
!   else
!     d = 4.0E+00 * r * ( sin ( 0.5E+00 * theta ) )**3 / &
!       ( 3.0E+00 * ( theta - sin ( theta ) ) )
!   end if
!
!   x = xc + d * cos ( theta )
!   y = yc + d * sin ( theta )
!
!   return
! end subroutine
! subroutine circle_sector_area_2d ( r, theta, area )
! !
! !*******************************************************************************
! !
! !! CIRCLE_SECTOR_AREA_2D computes the area of a circular sector in 2D.
! !
! !
! !  Modified:
! !
! !    12 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real THETA, the angle defining the size of the sector, in radians.
! !
! !    Output, real AREA, the area of the circle.
! !
!   real area
!   real r
!   real theta
! !
!   area = 0.5E+00 * r * r * theta
!
!   return
! end subroutine
! subroutine circle_sector_centroid_2d ( r, xc, yc, theta1, theta2, x, y )
! !
! !*******************************************************************************
! !
! !! CIRCLE_SECTOR_CENTROID_2D returns the centroid of a circular sector in 2D.
! !
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    28 February 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real XC, YC, the coordinates of the center of the circle.
! !
! !    Input, real THETA1, THETA2, the angles of the first and last points
! !    on the circular arc.
! !
! !    Output, real X, Y, the coordinates of the centroid of the sector.
! !
!   real area
!   real d
!   real r
!   real theta
!   real theta1
!   real theta2
!   real x
!   real xc
!   real y
!   real yc
! !
!   theta = theta2 - theta1
!
!   if ( theta == 0.0E+00 ) then
!     d = 2.0E+00 * r / 3.0E+00
!   else
!     d = 4.0E+00 * r * sin ( 0.5E+00 * theta ) / &
!       ( 3.0E+00 * theta )
!   end if
!
!   x = xc + d * cos ( theta )
!   y = yc + d * sin ( theta )
!
!   return
! end subroutine
! subroutine circle_triangle_area_2d ( r, theta, area )
! !
! !*******************************************************************************
! !
! !! CIRCLE_TRIANGLE_AREA_2D returns the area of a circle triangle in 2D.
! !
! !
! !  Discussion:
! !
! !    A circle triangle is formed by drawing a circular arc, and considering
! !    the triangle formed by the endpoints of the arc plus the center of
! !    the circle.
! !
! !    Note that for angles greater than PI, the triangle will actually
! !    have NEGATIVE area.
! !
! !  Modified:
! !
! !    12 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real THETA, the angle subtended by the arc.
! !
! !    Output, real AREA, the (signed) area of the triangle.
! !
!   real area
!   real r
!   real theta
! !
!   area = 0.5E+00 * r**2 * sin ( theta )
!
!   return
! end subroutine
! subroutine circle_dia2imp_2d ( x1, y1, x2, y2, r, xc, yc )
! !
! !*******************************************************************************
! !
! !! CIRCLE_DIA2IMP_2D converts a diameter to an implicit circle in 2D.
! !
! !
! !  Discussion:
! !
! !    The diameter form of a circle is:
! !
! !      (X1,Y1) and (X2,Y2) are endpoints of a diameter.
! !
! !    The implicit form of a circle in 2D is:
! !
! !      (X-XC)**2 + (Y-YC)**2 = R**2
! !
! !  Modified:
! !
! !    27 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, are the X and Y coordinates
! !    of two points which form a diameter of the circle.
! !
! !    Output, real R, the computed radius of the circle.
! !
! !    Output, real XC, YC, the computed center of the circle.
! !
!   real r
!   real x1
!   real x2
!   real xc
!   real y1
!   real y2
!   real yc
! !
!   r = 0.5E+00 * enorm0_2d ( x1, y1, x2, y2 )
!
!   xc = 0.5E+00 * ( x1 + x2 )
!   yc = 0.5E+00 * ( y1 + y2 )
!
!   return
! end subroutine
! subroutine circle_exp_contains_point_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )
! !
! !*******************************************************************************
! !
! !! CIRCLE_EXP_CONTAINS_POINT_2D determines if an explicit circle contains a point in 2D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a circle in 2D is:
! !
! !      The circle passing through (X1,Y1), (X2,Y2), (X3,Y3).
! !
! !  Modified:
! !
! !    26 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the (X,Y) coordinates of three
! !    points that lie on a circle.
! !
! !    Input, real X, Y, the (X,Y) coordinates of a point, whose position
! !    relative to the circle is desired.
! !
! !    Output, integer INSIDE:
! !
! !   -1, the three points are distinct and noncolinear,
! !    and (x,y) lies inside the circle.
! !    0, the three points are distinct and noncolinear,
! !    and (x,y) lies on the circle.
! !    1, the three points are distinct and noncolinear,
! !    and (x,y) lies outside the circle.
! !
! !    2, the three points are distinct and colinear,
! !    and (x,y) lies on the line.
! !    3, the three points are distinct and colinear,
! !    and (x,y) does not lie on the line.
! !
! !    4, two points are distinct, and (x,y) lies on the line.
! !    5, two points are distinct, and (x,y) does not lie on the line.
! !
! !    6, all three points are equal, and (x,y) is equal to them,
! !    7, all three points are equal, and (x,y) is not equal to them.
! !
!   real a(4,4)
!   real det
!   integer inside
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
! !  P1 = P2?
! !
!   if ( x1 == x2 .and. y1 == y2 ) then
!
!     if ( x1 == x3 .and. y1 == y3 ) then
!
!       if ( x1 == x .and. y1 == y ) then
!         inside = 6
!       else
!         inside = 7
!       end if
!
!     else
!
!       det = ( x1 - x3 ) * ( y - y3 ) - ( x - x3 ) * ( y1 - y3 )
!
!       if ( det == 0.0E+00 ) then
!         inside = 4
!       else
!         inside = 5
!       end if
!     end if
!
!     return
!
!   end if
! !
! !  P1 does not equal P2.  Does P1 = P3?
! !
!   if ( x1 == x3 .and. y1 == y3 ) then
!
!     det = ( x1 - x2 ) * ( y - y2 ) - ( x - x2 ) * ( y1 - y2 )
!
!     if ( det == 0.0E+00 ) then
!       inside = 4
!     else
!       inside = 5
!     end if
!
!     return
!
!   end if
! !
! !  The points are distinct.  Are they colinear?
! !
!   det = ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 )
!
!   if ( det == 0.0E+00 ) then
!
!     det = ( x1 - x2 ) * ( y - y2 ) - ( x - x2 ) * ( y1 - y2 )
!
!     if ( det == 0.0E+00 ) then
!       inside = 2
!     else
!       inside = 3
!     end if
!
!     return
!
!   end if
! !
! !  The points are distinct and non-colinear.
! !
! !  Compute the determinant
! !
!   a(1,1) = x1
!   a(1,2) = y1
!   a(1,3) = x1 * x1 + y1 * y1
!   a(1,4) = 1.0E+00
!
!   a(2,1) = x2
!   a(2,2) = y2
!   a(2,3) = x2 * x2 + y2 * y2
!   a(2,4) = 1.0E+00
!
!   a(3,1) = x3
!   a(3,2) = y3
!   a(3,3) = x3 * x3 + y3 * y3
!   a(3,4) = 1.0E+00
!
!   a(4,1) = x
!   a(4,2) = y
!   a(4,3) = x * x + y * y
!   a(4,4) = 1.0E+00
!
!   det = rmat4_det ( a )
!
!   if ( det < 0.0E+00 ) then
!     inside = 1
!   else if ( det == 0.0E+00 ) then
!     inside = 0
!   else
!     inside = -1
!   end if
!
!   return
! end subroutine
! subroutine circle_exp2imp_2d ( x1, y1, x2, y2, x3, y3, r, xc, yc )
! !
! !*******************************************************************************
! !
! !! CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a circle in 2D is:
! !
! !      The circle passing through (X1,Y1), (X2,Y2), (X3,Y3).
! !
! !    The implicit form of a circle in 2D is:
! !
! !      (X-XC)**2 + (Y-YC)**2 = R**2
! !
! !  Discussion:
! !
! !    Any three points define a circle, as long as they don't lie on a straight
! !    line.  (If the points do lie on a straight line, we could stretch the
! !    definition of a circle to allow an infinite radius and a center at
! !    some infinite point.)
! !
! !    Instead of the formulas used here, you can use the linear system
! !    approach in the routine TRIANGLE_OUTCIRCLE_2D.
! !
! !    The diameter of the circle can be found by solving a 2 by 2 linear system.
! !    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
! !    and each forms a right triangle with the diameter.  Hence, the dot product
! !    of P2 - P1 with the diameter is equal to the square of the length
! !    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
! !    diameter vector originating at P1.
! !
! !  Reference:
! !
! !    Joseph O'Rourke,
! !    Computational Geometry,
! !    Cambridge University Press,
! !    Second Edition, 1998, page 187.
! !
! !  Modified:
! !
! !    12 June 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, are the X and Y coordinates
! !    of three points that lie on the circle.  These points should be
! !    distinct, and not collinear.
! !
! !    Output, real R, the radius of the circle.  Normally, R will be positive.
! !    R is returned as -1 in the unlikely event that the points are
! !    numerically collinear.
! !
! !    Output, real XC, YC, the center of the circle.
! !
!   real a
!   real b
!   real c
!   real d
!   real e
!   real f
!   real g
!   real r
!   real x1
!   real x2
!   real x3
!   real xc
!   real y1
!   real y2
!   real y3
!   real yc
! !
!   a = x2 - x1
!   b = y2 - y1
!   c = x3 - x1
!   d = y3 - y1
!
!   e = a * ( x1 + x2 ) + b * ( y1 + y2 )
!   f = c * ( x1 + x3 ) + d * ( y1 + y3 )
! !
! !  Our formula is:
! !
! !    G = a * ( d - b ) - b * ( c - a )
! !
! !  but we get slightly better results using the original data.
! !
!   g = a * ( y3 - y2 ) - b * ( x3 - x2 )
! !
! !  We check for collinearity.  A more useful check would compare the
! !  absolute value of G to a small quantity.
! !
!   if ( g == 0.0E+00 ) then
!     xc = 0.0E+00
!     yc = 0.0E+00
!     r = -1.0E+00
!     return
!   end if
! !
! !  The center is halfway along the diameter vector from (X1,Y1).
! !
!   xc = 0.5E+00 * ( d * e - b * f ) / g
!   yc = 0.5E+00 * ( a * f - c * e ) / g
! !
! !  Knowing the center, the radius is now easy to compute.
! !
!   r = sqrt ( ( x1 - xc )**2 + ( y1 - yc )**2 )
!
!   return
! end subroutine
! subroutine circle_imp_contains_point_2d ( r, xc, yc, x, y, inside )
! !
! !*******************************************************************************
! !
! !! CIRCLE_IMP_CONTAINS_POINT_2D determines if an implicit circle contains a point in 2D.
! !
! !
! !  Formula:
! !
! !    An implicit circle in 2D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 = R**2
! !
! !  Modified:
! !
! !    21 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real XC, YC, the coordinates of the center of the circle.
! !
! !    Input, real X, Y, the point to be checked.
! !
! !    Output, logical INSIDE, is TRUE if the point is inside or on the circle,
! !    FALSE otherwise.
! !
!   logical inside
!   real r
!   real x
!   real xc
!   real y
!   real yc
! !
!   if ( enormsq0_2d ( x, y, xc, yc ) <= r * r ) then
!     inside = .true.
!   else
!     inside = .false.
!   end if
!
!   return
! end subroutine
! subroutine circle_imp_line_par_int_2d ( r, xc, yc, x0, y0, f, g, num_int, x, y )
! !
! !*******************************************************************************
! !
! !! CIRCLE_IMP_LINE_PAR_INT_2D: ( implicit circle, parametric line ) intersection in 2D.
! !
! !
! !  Formula:
! !
! !    An implicit circle in 2D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 = R**2
! !
! !    The parametric form of a line in 2D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real XC, YC, the coordinates of the center of the circle.
! !
! !    Input, real F, G, X0, Y0, the parametric parameters of the line.
! !
! !    Output, integer NUM_INT, the number of intersecting points found.
! !    NUM_INT will be 0, 1 or 2.
! !
! !    Output, real X(2), Y(2), contains the X and Y coordinates of
! !    the intersecting points.
! !
!   real f
!   real g
!   integer num_int
!   real r
!   real root
!   real t
!   real x(2)
!   real x0
!   real xc
!   real y(2)
!   real y0
!   real yc
! !
!   root = r * r * ( f * f + g * g ) - ( f * ( yc - y0 ) - g * ( xc - x0 ) )**2
!
!   if ( root < 0.0E+00 ) then
!
!     num_int = 0
!
!   else if ( root == 0.0E+00 )then
!
!     num_int = 1
!
!     t = ( f * ( xc - x0 ) + g * ( yc - y0 ) ) / ( f * f + g * g )
!     x(1) = x0 + f * t
!     y(1) = y0 + g * t
!
!   else if ( root > 0.0E+00 ) then
!
!     num_int = 2
!
!     t = ( ( f * ( xc - x0 ) + g * ( yc - y0 ) ) - sqrt ( root ) ) &
!       / ( f * f + g * g )
!
!     x(1) = x0 + f * t
!     y(1) = y0 + g * t
!
!     t = ( ( f * ( xc - x0 ) + g * ( yc - y0 ) ) + sqrt ( root ) ) &
!       / ( f * f + g * g )
!
!     x(2) = x0 + f * t
!     y(2) = y0 + g * t
!
!   end if
!
!   return
! end subroutine
! subroutine circle_imp_points_2d ( r, xc, yc, n, x, y )
! !
! !*******************************************************************************
! !
! !! CIRCLE_IMP_POINTS_2D returns N equally spaced points on an implicit circle in 2D.
! !
! !
! !  Note:
! !
! !    The first point is always ( XC + R, YC ), and subsequent points
! !    proceed counterclockwise around the circle.
! !
! !  Definition:
! !
! !    An implicit circle in 2D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 = R**2
! !
! !  Modified:
! !
! !    24 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real XC, YC, the coordinates of the center of the circle.
! !
! !    Input, integer N, the number of points desired.  N must be at least 1.
! !
! !    Output, real X(N), Y(N), the coordinates of points on the circle.
! !
!   integer n
! !
!   integer i
!   real r
!   real theta
!   real x(n)
!   real xc
!   real y(n)
!   real yc
! !
!   do i = 1, n
!     theta = ( 2.0E+00 * pi * real ( i - 1 ) ) / real ( n )
!     x(i) = xc + r * cos ( theta )
!     y(i) = yc + r * sin ( theta )
!   end do
!
!   return
! end subroutine
! subroutine circle_imp_points_arc_2d ( r, xc, yc, theta1, theta2, n, x, y )
! !
! !*******************************************************************************
! !
! !! CIRCLE_IMP_POINTS_ARC_2D returns N points on an arc of an implicit circle in 2D.
! !
! !
! !  Discussion:
! !
! !    The first point is ( XC + R * COS ( THETA1 ), YC + R * SIN ( THETA1 ) );
! !    The last point is  ( XC + R * COS ( THETA2 ), YC + R * SIN ( THETA2 ) );
! !    and the intermediate points are evenly spaced in angle between these,
! !    and in counterclockwise order.
! !
! !  Definition:
! !
! !    An implicit circle in 2D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 = R**2
! !
! !  Modified:
! !
! !    24 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle.
! !
! !    Input, real XC, YC, the coordinates of the center of the circle.
! !
! !    Input, real THETA1, THETA2, the angular coordinates of the first
! !    and last points to be drawn, in radians.
! !
! !    Input, integer N, the number of points desired.  N must be at least 1.
! !
! !    Output, real X(N), Y(N), the coordinates of points on the circle.
! !
!   integer n
! !
!   integer i
!   real r
!   real theta
!   real theta1
!   real theta2
!   real theta3
!   real x(n)
!   real xc
!   real y(n)
!   real yc
! !
! !  THETA3 is the smallest angle, no less than THETA1, which
! !  coincides with THETA2.
! !
!   theta3 = theta1 + r_modp ( theta2 - theta1, 2.0E+00 * real(pi) )
!
!   do i = 1, n
!
!     if ( n > 1 ) then
!       theta = ( real ( n - i ) * theta1 + real ( i - 1 ) * theta3 ) &
!         / real ( n - 1 )
!     else
!       theta = 0.5E+00 * ( theta1 + theta3 )
!     end if
!
!     x(i) = xc + r * cos ( theta )
!     y(i) = yc + r * sin ( theta )
!
!   end do
!
!   return
! end subroutine
! subroutine cone_area_3d ( h, r, area )
! !
! !*******************************************************************************
! !
! !! CONE_AREA_3D computes the surface area of a right circular cone in 3D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real H, R, the height of the cone, and the radius of the
! !    circle that forms the base of the cone.
! !
! !    Output, real AREA, the surface area of the cone.
! !
!   real area
!   real h
!   real r
! !
!   area = pi * r * sqrt ( h * h + r * r )
!
!   return
! end subroutine
! subroutine cone_centroid_3d ( r, xc, yc, zc, xh, yh, zh, x, y, z )
! !
! !*******************************************************************************
! !
! !! CONE_CENTROID_2D returns the centroid of a cone in 3D.
! !
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    28 February 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the circle at the base of the cone.
! !
! !    Input, real XC, YC, ZC, the coordinates of the center of the circle.
! !
! !    Input, real XH, YH, ZH, the coordinates of the tip of the cone.
! !
! !    Output, real X, Y, Z, the coordinates of the centroid of the cone.
! !
!   real r
!   real x
!   real xc
!   real xh
!   real y
!   real yc
!   real yh
!   real z
!   real zc
!   real zh
! !
!   x = 0.75 * xc + 0.25 * xh
!   y = 0.75 * yc + 0.25 * yh
!   z = 0.75 * zc + 0.25 * zh
!
!   return
! end subroutine
! subroutine cone_volume_3d ( h, r, volume )
! !
! !*******************************************************************************
! !
! !! CONE_VOLUME_3D computes the volume of a right circular cone in 3D.
! !
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real H, R, the height of the cone, and the radius of the
! !    circle that forms the base of the cone.
! !
! !    Output, real VOLUME, the volume of the cone.
! !
!   real h
!   real r
!   real volume
! !
!   volume = pi * r * r * h / 3.0E+00
!
!   return
! end subroutine
! subroutine conv3d ( cor, cor2, cor3, maxcor2, maxcor3, ncor2, theta )
! !
! !*******************************************************************************
! !
! !! CONV3D converts 3D data to a 2D projection.
! !
! !
! !  Discussion:
! !
! !    A "presentation angle" THETA is used to project the 3D point
! !    (X3D, Y3D, Z3D) to the 2D projection (XVAL,YVAL).
! !
! !  Formula:
! !
! !    If COR = 'X':
! !
! !      X2D = Y3D - sin ( THETA ) * X3D
! !      Y2D = Z3D - sin ( THETA ) * X3D
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, character COR, the coordinate to be projected.
! !    COR should be 'X', 'Y', or 'Z'.
! !
! !    Output, real COR2(2,NCOR), the 2D projections.
! !
! !    Input, real COR3(3,MAXCOR3), the X, Y, and Z components of points.
! !
! !    Input, integer MAXCOR2, MAXCOR3, the maximum number of 2D and
! !    3D points allowed.
! !
! !    Input, integer NCOR2, the number of 2D values to be computed.
! !
! !    Input, real THETA, the presentation angle in degrees.
! !
!   integer maxcor2
!   integer maxcor3
! !
!   character cor
!   real cor2(2,maxcor2)
!   real cor3(3,maxcor3)
!   integer i
!   integer ncor2
!   real stheta
!   real theta
! !
!   stheta = sin ( degrees_to_radians ( theta ) )
!
!   if ( cor=='X' ) then
!
!     cor2(1,1:ncor2) = cor3(2,1:ncor2) - stheta * cor3(1,1:ncor2)
!     cor2(2,1:ncor2) = cor3(3,1:ncor2) - stheta * cor3(1,1:ncor2)
!
!   else if ( cor=='Y' ) then
!
!     cor2(1,1:ncor2) = cor3(1,1:ncor2) - stheta * cor3(2,1:ncor2)
!     cor2(2,1:ncor2) = cor3(3,1:ncor2) - stheta * cor3(2,1:ncor2)
!
!   else
!
!     cor2(1,1:ncor2) = cor3(1,1:ncor2) - stheta * cor3(3,1:ncor2)
!     cor2(2,1:ncor2) = cor3(2,1:ncor2) - stheta * cor3(3,1:ncor2)
!
!   end if
!
!   return
! end subroutine
! subroutine corpl_2d ( dist, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5 )
! !
! !*******************************************************************************
! !
! !! CORPL_2D "boxes" an angle defined by three points in 2D.
! !
! !
! !  Discussion:
! !
! !    CORPL2 is given
! !      (X1,Y1), (X2,Y2), and (X3,Y3)
! !    determining the two lines:
! !      (X1,Y1), (X2,Y2)
! !    and
! !      (X2,Y2), (X3,Y3)
! !    and a nonnegative distance
! !      DIST.
! !
! !    CORPL2 returns a pair of "corner" points
! !      (X4,Y4) and (X5,Y5)
! !    both of which are a distance DIST from both lines, and in fact,
! !    both of which are a distance DIST from (X2,Y2).
! !
! !                         /   3
! !                        /   /   /
! !     - - - - - - - - - 4 - / - 6 - - -
! !                      /   /   /
! !     1---------------/---2-----------------
! !                    /   /   /
! !     - - - - - - - 7 - / - 5 - - - - -
! !                  /   /   /
! !
! !    In the illustration, the numbers "1" "2" and "3" represent
! !    the points defining the lines.
! !
! !    The numbers "4" and "5" represent the desired "corner points", which
! !    are on the positive or negative sides of both lines.
! !
! !    The numbers "6" and "7" represent the undesired points, which
! !    are on the positive side of one line and the negative of the other.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real DIST, the nonnegative distance from (X1,Y1) to the computed
! !    points (X3,Y3) and (X4,Y4).
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3.
! !    (X1,Y1) and (X2,Y2) are distinct points that define a line.
! !    (X2,Y2) and (X3,Y3) are distinct points that define a line.
! !
! !    Special cases:
! !
! !    if ( X1,Y1) = (X2,Y2), this is the same as extending the line from
! !    (X3,Y3) through (X2,Y2) without a bend.
! !
! !    if ( X3,Y3) = (X2,Y2), this is the same as extending the line from
! !    (X1,Y1) through (X2,Y2) without a bend.
! !
! !    if ( X1,Y1) = (X2,Y2) = (X3,Y3) this is an error.
! !
! !    Output, real X4, Y4, X5, Y5.
! !    (X4,Y4) and (X5,Y5) are points which lie DIST units from
! !    the line between (X1,Y1) and (X2,Y2), and DIST units from
! !    the line between (X2,Y2) and (X3,Y3).
! !
!   real dist
!   real stheta
!   real temp
!   real ux
!   real ux1
!   real ux2
!   real uy
!   real uy1
!   real uy2
!   real x1
!   real x1copy
!   real x2
!   real x3
!   real x3copy
!   real x4
!   real x5
!   real y1
!   real y1copy
!   real y2
!   real y3
!   real y3copy
!   real y4
!   real y5
! !
! !  Fail if all three points are equal.
! !
!   if ( x1 == x2 .and. x2 == x3 .and. y1 == y2 .and. y2 == y3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'CORPL2 - Fatal error!'
!     write ( *, * ) '  Input points (X1,Y1) = (X2,Y2) = (X3,Y3).'
!     write ( *, * ) '  (X1,Y1)= ', x1, y1
!     stop
!   end if
! !
! !  If P1 = P2 or P2 = P3, extend the line through the doubled point.
! !
!   if ( x1 == x2 .and. y1 == y2 ) then
!     x1copy = 2.0E+00 * x2 - x3
!     y1copy = 2.0E+00 * y2 - y3
!   else
!     x1copy = x1
!     y1copy = y1
!   end if
!
!   if ( x2 == x3 .and. y2 == y3 ) then
!     x3copy = 2.0E+00 * x2 - x1
!     y3copy = 2.0E+00 * y2 - y1
!   else
!     x3copy = x3
!     y3copy = y3
!   end if
! !
! !  Now compute the unit normal vectors to each line.
! !  We choose the sign so that the unit normal to line 1 has
! !  a positive dot product with line 2.
! !
!   ux1 = y1copy - y2
!   uy1 = x2 - x1copy
!   temp = sqrt ( ux1 * ux1 + uy1 * uy1 )
!   ux1 = ux1 / temp
!   uy1 = uy1 / temp
!
!   if ( ux1 * ( x3copy - x2 ) + uy1 * ( y3copy - y2 ) < 0.0E+00 ) then
!     ux1 = - ux1
!     uy1 = - uy1
!   end if
!
!   ux2 = y3copy - y2
!   uy2 = x2 - x3copy
!   temp = sqrt ( ux2 * ux2 + uy2 * uy2)
!   ux2 = ux2 / temp
!   uy2 = uy2 / temp
!
!   if ( ux2 * ( x1copy - x2 ) + uy2 * ( y1copy - y2 ) < 0.0E+00 ) then
!     ux2 = - ux2
!     uy2 = - uy2
!   end if
! !
! !  Try to catch the case where we can't determine the
! !  sign of U1, because both U1 and -U1 are perpendicular
! !  to (P3-P2)...and similarly for U2 and (P1-P2).
! !
!   if ( ux1 * ( x3copy - x2) + uy1 * ( y3copy - y2 ) == 0.0E+00 .or. &
!        ux2 * ( x1copy - x2) + uy2 * ( y1copy - y2 ) == 0.0E+00 ) then
!
!     if ( ux1 * ux2 + uy1 * uy2 < 0.0E+00 ) then
!       ux1 = -ux1
!       uy1 = -uy1
!     end if
!
!   end if
! !
! !  Try to catch a line turning back on itself, evidenced by
! !    Cos(theta) = (P3-P2) dot (P2-P1) / ( norm(P3-P2) * norm(P2-P1) )
! !  being -1, or very close to -1.
! !
!   temp = ( ( x3copy - x2 ) * ( x2 - x1copy ) &
!          + ( y3copy - y2 ) * ( y2 - y1copy ) )
!
!   temp = temp / &
!         ( sqrt ( ( x3copy - x2 )**2 + ( y3copy - y2 )**2 ) &
!         * sqrt ( ( x2 - x1copy )**2 + ( y2 - y1copy )**2 ) )
!
!   if ( temp < -0.99 ) then
!     temp = sqrt ( ( x2 - x1copy)**2 + ( y2 - y1copy )**2 )
!     x4 = x2 + dist * ( x2 - x1copy ) / temp + dist * ux1
!     y4 = y2 + dist * ( y2 - y1copy ) / temp + dist * uy1
!     x5 = x2 + dist * ( x2 - x1copy ) / temp - dist * ux1
!     y5 = y2 + dist * ( y2 - y1copy ) / temp - dist * uy1
!     return
!   end if
! !
! !  Compute the "average" unit normal vector.
! !
! !  The average of the unit normals could be zero, but only when
! !  the second line has the same direction and opposite sense
! !  of the first, and we've already checked for that case.
! !
!   ux = 0.5E+00 * ( ux1 + ux2 )
!   uy = 0.5E+00 * ( uy1 + uy2 )
!   temp = sqrt ( ux * ux + uy * uy )
!   ux = ux / temp
!   uy = uy / temp
! !
! !  You must go DIST/STHETA units along this unit normal to
! !  result in a distance DIST from line1 (and line2).
! !
!   stheta = ux * ux1 + uy * uy1
!
!   x4 = x2 + dist * ux / stheta
!   y4 = y2 + dist * uy / stheta
!
!   x5 = x2 - dist * ux / stheta
!   y5 = y2 - dist * uy / stheta
!
!   return
! end subroutine
! function cot ( angle )
! !
! !*******************************************************************************
! !
! !! COT returns the cotangent of an angle.
! !
! !
! !  Modified:
! !
! !    12 July 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real ANGLE, the angle, in radians.
! !
! !    Output, real COT, the cotangent of the angle.
! !
!   real angle
!   real cot
! !
!   cot  = cos ( angle ) / sin ( angle )
!
!   return
! end function
! function cotd ( angle )
! !
! !*******************************************************************************
! !
! !! COTD returns the cotangent of an angle given in degrees.
! !
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real ANGLE, the angle, in degrees.
! !
! !    Output, real COTD, the cotangent of the angle.
! !
!   real angle
!   real angle_rad
!   real cotd
! !
!   angle_rad = degrees_to_radians ( angle )
!
!   cotd  = cos ( angle_rad ) / sin ( angle_rad )
!
!   return
! end function
! function cross_2d ( x1, y1, x2, y2 )
! !
! !*******************************************************************************
! !
! !! CROSS_2D finds the cross product of a pair of vectors in 2D.
! !
! !
! !  Discussion:
! !
! !    Strictly speaking, the vectors lie in the (X,Y) plane, and
! !    the cross product here is a vector in the Z direction.
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, the X and Y coordinates of the vectors.
! !
! !    Output, real CROSS_2D, the Z component of the cross product
! !    of (X1,Y1,0) and (X2,Y2,0).
! !
!   real cross_2d
!   real x1
!   real x2
!   real y1
!   real y2
! !
!   cross_2d = x1 * y2 - y1 * x2
!
!   return
! end function
! subroutine cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
! !
! !*******************************************************************************
! !
! !! CROSS_3D computes the cross product of two vectors in 3D.
! !
! !
! !  Definition:
! !
! !    The cross product in 3D can be regarded as the determinant of the
! !    symbolic matrix:
! !
! !          |  i  j  k |
! !      det | x1 y1 z1 |
! !          | x2 y2 z2 |
! !
! !      = ( y1 * z2 - z1 * y2 ) * i
! !      + ( z1 * x2 - x1 * z2 ) * j
! !      + ( x1 * y2 - y1 * x2 ) * k
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
! !
! !    Output, real X3, Y3, Z3, the cross product vector.
! !
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
!   real z1
!   real z2
!   real z3
! !
!   x3 = y1 * z2 - z1 * y2
!   y3 = z1 * x2 - x1 * z2
!   z3 = x1 * y2 - y1 * x2
!
!   return
! end subroutine
! function cross0_2d ( x0, y0, x1, y1, x2, y2 )
! !
! !*******************************************************************************
! !
! !! CROSS0_2D finds the cross product of (P1-P0) and (P2-P0) in 2D.
! !
! !
! !  Discussion:
! !
! !    Strictly speaking, the vectors lie in the (X,Y) plane, and
! !    the cross product here is a vector in the Z direction.
! !
! !    The vectors are specified with respect to a basis point P0.
! !    We are computing the normal to the triangle (P0,P1,P2).
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real X3, Y3, Z3, the cross product (P1-P0) x (P2-P0), or
! !    (X1-X0,Y1-Y0,Z1-Z0) x (X2-X0,Y2-Y0,Z2-Z0).
! !
! !    Input, real X0, Y0, X1, Y1, X2, Y2, the coordinates of the three
! !    points.  The basis point P0 is (X0,Y0).
! !
! !    Output, real CROSS0_2D, the Z component of the cross product
! !    (X1-X0,Y1-Y0,0) x (X2-X0,Y2-Y0,0).
! !
!   real cross0_2d
!   real x0
!   real x1
!   real x2
!   real y0
!   real y1
!   real y2
! !
!   cross0_2d = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )
!
!   return
! end function
! subroutine cross0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3 )
! !
! !*******************************************************************************
! !
! !! CROSS0_3D computes the cross product of (P1-P0) and (P2-P0) in 3D.
! !
! !
! !  Discussion:
! !
! !    The vectors are specified with respect to a basis point P0.
! !    We are computing the normal to the triangle (P0,P1,P2).
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, the coordinates of
! !    three points.  The basis point is (X0,Y0,Z0).
! !
! !    Output, real X3, Y3, Z3, the cross product (P1-P0) x (P2-P0), or
! !    (X1-X0,Y1-Y0,Z1-Z0) x (X2-X0,Y2-Y0,Z2-Z0).
! !
!   real x0
!   real x1
!   real x2
!   real x3
!   real y0
!   real y1
!   real y2
!   real y3
!   real z0
!   real z1
!   real z2
!   real z3
! !
!   x3 = ( y1 - y0 ) * ( z2 - z0 ) - ( z1 - z0 ) * ( y2 - y0 )
!   y3 = ( z1 - z0 ) * ( x2 - x0 ) - ( x1 - x0 ) * ( z2 - z0 )
!   z3 = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )
!
!   return
! end subroutine
! subroutine cube_shape_3d ( max_num, max_order, point_num, face_num, &
!   face_order, point_coord, face_point )
! !
! !*******************************************************************************
! !
! !! CUBE_SHAPE_3D describes a cube in 3D.
! !
! !
! !  Discussion:
! !
! !    The vertices of the cube lie on the unit sphere.
! !
! !  Modified:
! !
! !    11 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAX_NUM, the maximum number of faces, and of
! !    points.  This quantity is used to dimension arrays.
! !    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
! !    not be set.
! !
! !    Input, integer MAX_ORDER, the maximum order of any face.
! !    This quantity is used to dimension arrays.  If FACE_ORDER
! !    exceeds MAX_ORDER, the arrays will not be set.
! !
! !    Output, integer POINT_NUM, the number of points in the shape.
! !
! !    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
! !    the X, Y and Z coordinates of point J.
! !
! !    Output, integer FACE_NUM, the number of faces in the shape.
! !
! !    Output, integer FACE_ORDER, the number of vertices per face.
! !
! !    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
! !    contains the index of the I-th point in the J-th face.  The
! !    points are listed in the counter-clockwise direction defined
! !    by the outward normal at the face.
! !
!   integer max_num
!   integer max_order
! !
!   real a
!   integer face_num
!   integer face_order
!   integer face_point(max_order,max_num)
!   integer point_num
!   real point_coord(3,max_num)
! !
!   point_num = 8
!   face_num = 6
!   face_order = 4
! !
! !  Check.
! !
!   if ( point_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'CUBE_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of vertices exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'CUBE_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of faces exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_order > max_order ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'CUBE_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Face order exceeds MAX_ORDER.'
!     return
!   end if
! !
! !  Set point coordinates.
! !
!   a = sqrt ( 1.0E+00 / 3.0E+00 )
!
!   point_coord(1,1) = - a
!   point_coord(2,1) = - a
!   point_coord(3,1) = - a
!
!   point_coord(1,2) =   a
!   point_coord(2,2) = - a
!   point_coord(3,2) = - a
!
!   point_coord(1,3) =   a
!   point_coord(2,3) =   a
!   point_coord(3,3) = - a
!
!   point_coord(1,4) = - a
!   point_coord(2,4) =   a
!   point_coord(3,4) = - a
!
!   point_coord(1,5) = - a
!   point_coord(2,5) = - a
!   point_coord(3,5) =   a
!
!   point_coord(1,6) =   a
!   point_coord(2,6) = - a
!   point_coord(3,6) =   a
!
!   point_coord(1,7) =   a
!   point_coord(2,7) =   a
!   point_coord(3,7) =   a
!
!   point_coord(1,8) = - a
!   point_coord(2,8) =   a
!   point_coord(3,8) =   a
! !
! !  Set faces.
! !
!   face_point(1,1) = 1
!   face_point(2,1) = 4
!   face_point(3,1) = 3
!   face_point(4,1) = 2
!
!   face_point(1,2) = 1
!   face_point(2,2) = 2
!   face_point(3,2) = 6
!   face_point(4,2) = 5
!
!   face_point(1,3) = 2
!   face_point(2,3) = 3
!   face_point(3,3) = 7
!   face_point(4,3) = 6
!
!   face_point(1,4) = 3
!   face_point(2,4) = 4
!   face_point(3,4) = 8
!   face_point(4,4) = 7
!
!   face_point(1,5) = 1
!   face_point(2,5) = 5
!   face_point(3,5) = 8
!   face_point(4,5) = 4
!
!   face_point(1,6) = 5
!   face_point(2,6) = 6
!   face_point(3,6) = 7
!   face_point(4,6) = 8
!
!   return
! end subroutine
! function degrees_to_radians ( angle )
! !
! !*******************************************************************************
! !
! !! DEGREES_TO_RADIANS converts an angle from degrees to radians.
! !
! !
! !  Modified:
! !
! !    10 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real ANGLE, an angle in degrees.
! !
! !    Output, real DEGREES_TO_RADIANS, the equivalent angle
! !    in radians.
! !
!
!   real angle
!   real degrees_to_radians
!
!   degrees_to_radians = ( angle / 180.0E+00 ) * pi
!
!   return
! end function
! subroutine direction_pert_3d ( sigma, vbase, vran )
! !
! !*******************************************************************************
! !
! !! DIRECTION_PERT_3D randomly perturbs a direction vector in 3D.
! !
! !
! !  Modified:
! !
! !    01 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real SIGMA, determines the strength of the perturbation.
! !    SIGMA <= 0 results in a completely random direction.
! !    1 <= SIGMA results in VBASE.
! !    0 < SIGMA < 1 results in a perturbation from VBASE, which is
! !    large when SIGMA is near 0, and small when SIGMA is near 1.
! !
! !    Input, real VBASE(3), the base direction vector, which should have
! !    unit norm.
! !
! !    Output, real VRAN(3), the perturbed vector, which will have unit norm.
! !
!   real dphi
!   integer i
!   real phi
!   real psi
!   real r
!   real sigma
!   real theta
!   real v(3)
!   real vbase(3)
!   real vdot
!   real vran(3)
!   real x
! !
! !  SIGMA >= 0, just use the base vector.
! !
!   if ( sigma >= 1.0E+00 ) then
!
!     vran(1:3) = vbase(1:3)
!
!   else if ( sigma <= 0.0E+00 ) then
!
!     call r_random ( -1.0E+00, 1.0E+00, vdot )
!     phi = arc_cosine ( vdot )
!
!     call r_random ( 0.0E+00, 2.0E+00 * real(pi), theta )
!
!     vran(1) = cos ( theta ) * sin ( phi )
!     vran(2) = sin ( theta ) * sin ( phi )
!     vran(3) = cos ( phi )
!
!   else
!
!     phi = arc_cosine ( vbase(3) )
!     theta = atan2 ( vbase(2), vbase(1) )
! !
! !  Pick VDOT, which must be between -1 and 1.  This represents
! !  the dot product of the perturbed vector with the base vector.
! !
! !  RANDOM_NUMBER returns a uniformly random value between 0 and 1.
! !  The operations we perform on this quantity tend to bias it
! !  out towards 1, as SIGMA grows from 0 to 1.
! !
! !  VDOT, in turn, is a value between -1 and 1, which, for large
! !  SIGMA, we want biased towards 1.
! !
!     call r_random ( 0.0E+00, 1.0E+00, r )
!     x = exp ( ( 1.0E+00 - sigma ) * log ( r ) )
!     dphi = arc_cosine ( 2.0E+00 * x - 1.0E+00 )
! !
! !  Now we know enough to write down a vector that is rotated DPHI
! !  from the base vector.
! !
!     v(1) = cos ( theta ) * sin ( phi + dphi )
!     v(2) = sin ( theta ) * sin ( phi + dphi )
!     v(3) = cos ( phi + dphi )
! !
! !  Pick a uniformly random rotation between 0 and 2 Pi around the
! !  axis of the base vector.
! !
!     call r_random ( 0.0E+00, 2.0E+00 *real(pi), psi )
! !
! !  Carry out the rotation.
! !
!     call rotation_axis_vector_3d ( vbase, psi, v, vran )
!
!   end if
!
!   return
! end subroutine
! subroutine direction_random_3d ( vran )
! !
! !*******************************************************************************
! !
! !! DIRECTION_RANDOM_3D picks a random direction vector in 3D.
! !
! !
! !  Modified:
! !
! !    01 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real VRAN(3), the random direction vector, with unit norm.
! !
!   real phi
!   real r
!   real theta
!   real vdot
!   real vran(3)
! !
! !  Pick a uniformly random VDOT, which must be between -1 and 1.
! !  This represents the dot product of the random vector with the Z unit vector.
! !
! !  Note: this works because the surface area of the sphere between
! !  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
! !  a patch of area uniformly.
! !
!   call r_random ( -1.0E+00, 1.0E+00, vdot )
!   phi = arc_cosine ( vdot )
! !
! !  Pick a uniformly random rotation between 0 and 2 Pi around the
! !  axis of the Z vector.
! !
!   call r_random ( 0.0E+00, 2.0E+00 * real(pi), theta )
!
!   vran(1) = cos ( theta ) * sin ( phi )
!   vran(2) = sin ( theta ) * sin ( phi )
!   vran(3) = cos ( phi )
!
!   return
! end subroutine
! subroutine direction_random_nd ( n, w )
! !
! !*******************************************************************************
! !
! !! DIRECTION_RANDOM_ND generates a random direction vector in ND.
! !
! !
! !  Modified:
! !
! !    17 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Output, real W(N), a random direction vector, with unit norm.
! !
!   integer n
! !
!   integer i
!   real norm
!   real w(n)
! !
! !  Get N values from a standard normal distribution.
! !
!   do i = 1, n
!     call normal_01_sample ( w(i) )
!   end do
! !
! !  Compute the length of the vector.
! !
!   norm = sqrt ( sum ( w(1:n)**2 ) )
! !
! !  Normalize the vector.
! !
!   w(1:n) = w(1:n) / norm
!
!   return
! end subroutine
! subroutine dodec_shape_3d ( max_num, max_order, point_num, face_num, &
!   face_order, point_coord, face_point )
! !
! !*******************************************************************************
! !
! !! DODEC_SHAPE_3D describes an dodecahedron in 3D.
! !
! !
! !  Discussion:
! !
! !    The vertices of the dodecahedron lie on the unit sphere.
! !
! !  Modified:
! !
! !    12 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAX_NUM, the maximum number of faces, and of
! !    points.  This quantity is used to dimension arrays.
! !    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
! !    not be set.
! !
! !    Input, integer MAX_ORDER, the maximum order of any face.
! !    This quantity is used to dimension arrays.  If FACE_ORDER
! !    exceeds MAX_ORDER, the arrays will not be set.
! !
! !    Output, integer POINT_NUM, the number of points in the shape.
! !
! !    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
! !    the X, Y and Z coordinates of point J.
! !
! !    Output, integer FACE_NUM, the number of faces in the shape.
! !
! !    Output, integer FACE_ORDER, the number of vertices per face.
! !
! !    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
! !    contains the index of the I-th point in the J-th face.  The
! !    points are listed in the counter-clockwise direction defined
! !    by the outward normal at the face.
! !
!   integer max_num
!   integer max_order
! !
!   real a
!   real b
!   real c
!   integer face_num
!   integer face_order
!   integer face_point(max_order,max_num)
!   real phi
!   integer point_num
!   real point_coord(3,max_num)
! !
!   point_num = 20
!   face_num = 12
!   face_order = 5
! !
! !  Check.
! !
!   if ( point_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'DODEC_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of vertices exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'DODEC_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of faces exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_order > max_order ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'DODEC_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Face order exceeds MAX_ORDER.'
!     return
!   end if
! !
! !  Set point coordinates.
! !
!   phi = 0.5E+00 * ( sqrt ( 5.0E+00 ) + 1.0E+00 )
!
!   a = 1.0E+00 / sqrt ( 3.0E+00 )
!   b = phi / sqrt ( 3.0E+00 )
!   c = ( phi - 1.0E+00 ) / sqrt ( 3.0E+00 )
!
!   point_coord(1,1) =   a
!   point_coord(2,1) =   a
!   point_coord(3,1) =   a
!
!   point_coord(1,2) =   a
!   point_coord(2,2) =   a
!   point_coord(3,2) = - a
!
!   point_coord(1,3) =   a
!   point_coord(2,3) = - a
!   point_coord(3,3) =   a
!
!   point_coord(1,4) =   a
!   point_coord(2,4) = - a
!   point_coord(3,4) = - a
!
!   point_coord(1,5) = - a
!   point_coord(2,5) =   a
!   point_coord(3,5) =   a
!
!   point_coord(1,6) = - a
!   point_coord(2,6) =   a
!   point_coord(3,6) = - a
!
!   point_coord(1,7) = - a
!   point_coord(2,7) = - a
!   point_coord(3,7) =   a
!
!   point_coord(1,8) = - a
!   point_coord(2,8) = - a
!   point_coord(3,8) = - a
!
!   point_coord(1,9) =    c
!   point_coord(2,9) =    b
!   point_coord(3,9) =    0.0E+00
!
!   point_coord(1,10) = - c
!   point_coord(2,10) =   b
!   point_coord(3,10) =   0.0E+00
!
!   point_coord(1,11) =   c
!   point_coord(2,11) = - b
!   point_coord(3,11) =   0.0E+00
!
!   point_coord(1,12) = - c
!   point_coord(2,12) = - b
!   point_coord(3,12) =   0.0E+00
!
!   point_coord(1,13) =   b
!   point_coord(2,13) =   0.0E+00
!   point_coord(3,13) =   c
!
!   point_coord(1,14) =   b
!   point_coord(2,14) =   0.0E+00
!   point_coord(3,14) = - c
!
!   point_coord(1,15) = - b
!   point_coord(2,15) =   0.0E+00
!   point_coord(3,15) =   c
!
!   point_coord(1,16) = - b
!   point_coord(2,16) =   0.0E+00
!   point_coord(3,16) = - c
!
!   point_coord(1,17) =   0.0E+00
!   point_coord(2,17) =   c
!   point_coord(3,17) =   b
!
!   point_coord(1,18) =   0.0E+00
!   point_coord(2,18) = - c
!   point_coord(3,18) =   b
!
!   point_coord(1,19) =   0.0E+00
!   point_coord(2,19) =   c
!   point_coord(3,19) = - b
!
!   point_coord(1,20) =   0.0E+00
!   point_coord(2,20) = - c
!   point_coord(3,20) = - b
! !
! !  Set faces.
! !
!   face_point(1,1) = 2
!   face_point(2,1) = 9
!   face_point(3,1) = 1
!   face_point(4,1) = 13
!   face_point(5,1) = 14
!
!   face_point(1,2) = 5
!   face_point(2,2) = 10
!   face_point(3,2) = 6
!   face_point(4,2) = 16
!   face_point(5,2) = 15
!
!   face_point(1,3) = 3
!   face_point(2,3) = 11
!   face_point(3,3) = 4
!   face_point(4,3) = 14
!   face_point(5,3) = 13
!
!   face_point(1,4) = 8
!   face_point(2,4) = 12
!   face_point(3,4) = 7
!   face_point(4,4) = 15
!   face_point(5,4) = 16
!
!   face_point(1,5) = 3
!   face_point(2,5) = 13
!   face_point(3,5) = 1
!   face_point(4,5) = 17
!   face_point(5,5) = 18
!
!   face_point(1,6) = 2
!   face_point(2,6) = 14
!   face_point(3,6) = 4
!   face_point(4,6) = 20
!   face_point(5,6) = 19
!
!   face_point(1,7) = 5
!   face_point(2,7) = 15
!   face_point(3,7) = 7
!   face_point(4,7) = 18
!   face_point(5,7) = 17
!
!   face_point(1,8) = 8
!   face_point(2,8) = 16
!   face_point(3,8) = 6
!   face_point(4,8) = 19
!   face_point(5,8) = 20
!
!   face_point(1,9) = 5
!   face_point(2,9) = 17
!   face_point(3,9) = 1
!   face_point(4,9) = 9
!   face_point(5,9) = 10
!
!   face_point(1,10) = 3
!   face_point(2,10) = 18
!   face_point(3,10) = 7
!   face_point(4,10) = 12
!   face_point(5,10) = 11
!
!   face_point(1,11) = 2
!   face_point(2,11) = 19
!   face_point(3,11) = 6
!   face_point(4,11) = 10
!   face_point(5,11) = 9
!
!   face_point(1,12) = 8
!   face_point(2,12) = 20
!   face_point(3,12) = 4
!   face_point(4,12) = 11
!   face_point(5,12) = 12
!
!   return
! end subroutine
! function dot_2d ( x1, y1, x2, y2 )
! !
! !*******************************************************************************
! !
! !! DOT_2D computes the dot product of a pair of vectors in 2D.
! !
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, the coordinates of the vectors.
! !
! !    Output, real DOT_2D, the dot product.
! !
!   real dot_2d
!   real x1
!   real x2
!   real y1
!   real y2
! !
!   dot_2d = x1 * x2 + y1 * y2
!
!   return
! end function
! function dot_3d ( x1, y1, z1, x2, y2, z2 )
! !
! !*******************************************************************************
! !
! !! DOT_3D computes the dot product of a pair of vectors in 3D.
! !
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
! !
! !    Output, real DOT_3D, the dot product.
! !
!   real dot_3d
!   real x1
!   real x2
!   real y1
!   real y2
!   real z1
!   real z2
! !
!   dot_3d = x1 * x2 + y1 * y2 + z1 * z2
!
!   return
! end function
! function dot_nd ( n, v1, v2 )
! !
! !*******************************************************************************
! !
! !! DOT_ND computes the dot product of a pair of vectors in ND.
! !
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of entries in the vectors.
! !
! !    Input, real V1(N), V2(N), the two vectors.
! !
! !    Output, real DOT_ND, the dot product.
! !
!   integer n
! !
!   real dot_nd
!   integer i
!   real v1(n)
!   real v2(n)
! !
!   dot_nd = dot_product ( v1(1:n), v2(1:n) )
!
!   return
! end function
! function dot0_2d ( x0, y0, x1, y1, x2, y2 )
! !
! !*******************************************************************************
! !
! !! DOT0_2D computes the dot product of (P1-P0) and (P2-P0) in 2D.
! !
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, X1, Y1, X2, Y2, the coordinates of
! !    the points P0, P1 and P1.
! !
! !    Output, real DOT0_2D, the dot product of (P1-P0) and (P2-P0).
! !
!   real dot0_2d
!   real x0
!   real x1
!   real x2
!   real y0
!   real y1
!   real y2
! !
!   dot0_2d = ( x1 - x0 ) * ( x2 - x0 ) + ( y1 - y0 ) * ( y2 - y0 )
!
!   return
! end function
! function dot0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2 )
! !
! !*******************************************************************************
! !
! !! DOT0_3D computes the dot product of (P1-P0) and (P2-P0) in 3D.
! !
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, Z0, the coordinates of the point P0.
! !
! !    Input, real X1, Y1, Z1, the coordinates of the point P1.
! !
! !    Input, real X2, Y2, Z2, the coordinates of the point P2.
! !
! !    Output, real DOT0_3D, the dot product of (P1-P0) and (P2-P0).
! !
!   real dot0_3d
!   real x0
!   real x1
!   real x2
!   real y0
!   real y1
!   real y2
!   real z0
!   real z1
!   real z2
! !
!   dot0_3d = ( x1 - x0 ) * ( x2 - x0 ) + ( y1 - y0 ) * ( y2 - y0 ) + &
!             ( z1 - z0 ) * ( z2 - z0 )
!
!   return
! end function
! subroutine dual_shape_3d ( max_num, max_order, point_num1, face_num1, &
!   face_order1, point_coord1, face_point1, point_num2, face_num2, face_order2, &
!   point_coord2, face_point2 )
! !
! !*******************************************************************************
! !
! !! DUAL_SHAPE_3D constructs the dual of a shape in 3D.
! !
! !
! !  Modified:
! !
! !    09 August 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAX_NUM, the maximum number of faces, and of
! !    points.  This quantity is used to dimension arrays.
! !    If POINT_NUM1 or FACE_NUM1 exceeds MAX_NUM, the arrays will
! !    not be set.
! !
! !    Input, integer MAX_ORDER, the maximum order of any face.
! !    This quantity is used to dimension arrays.  If FACE_ORDER1
! !    or FACE_ORDER2 exceeds MAX_ORDER, the arrays will not be set.
! !
! !    Input, integer POINT_NUM1, the number of points in the shape.
! !
! !    Input, integer FACE_NUM1, the number of faces in the shape.
! !
! !    Input, integer FACE_ORDER1, the number of vertices per face.
! !
! !    Input, real POINT_COORD1(3,MAX_NUM); POINT_COORD1(*,J) contains
! !    the X, Y and Z coordinates of point J.
! !
! !    Input, integer FACE_POINT1(MAX_ORDER,MAX_NUM); FACE_POINT1(I,J)
! !    contains the index of the I-th point in the J-th face.  The
! !    points are listed in the counter-clockwise direction defined
! !    by the outward normal at the face.
! !
! !    Output, integer POINT_NUM2, the number of points in the dual.
! !
! !    Output, integer FACE_NUM2, the number of faces in the dual.
! !
! !    Output, integer FACE_ORDER2, the number of vertices per face
! !    in the dual.
! !
! !    Output, real POINT_COORD2(3,MAX_NUM), the point coordinates
! !    of the dual.
! !
! !    Output, integer FACE_POINT2(MAX_ORDER,MAX_NUM), the vertices
! !    of each face in the dual.
! !
!   integer max_num
!   integer max_order
! !
!   integer face_num1
!   integer face_num2
!   integer face_order1
!   integer face_order2
!   integer face_point1(max_order,max_num)
!   integer face_point2(max_order,max_num)
!   integer i
!   integer icol
!   integer iface
!   integer inext
!   integer iprev
!   integer irow
!   integer istop
!   integer j
!   integer k
!   real norm
!   integer point_num1
!   integer point_num2
!   real point_coord1(3,max_num)
!   real point_coord2(3,max_num)
!   real x
!   real y
!   real z
! !
!   point_num2 = face_num1
!   face_num2 = point_num1
! !
! !  This computation should really compute the center of gravity
! !  of the face, in the general case.
! !
! !  We'll also assume the vertices of the original and the dual
! !  are to lie on the unit sphere, so we can normalize the
! !  position vector of the vertex.
! !
!   do i = 1, face_num1
!
!     x = 0.0E+00
!     y = 0.0E+00
!     z = 0.0E+00
!     do j = 1, face_order1
!       k = face_point1(j,i)
!       x = x + point_coord1(1,k)
!       y = y + point_coord1(2,k)
!       z = z + point_coord1(3,k)
!     end do
!
!     norm = sqrt ( x * x + y * y + z * z )
!
!     point_coord2(1,i) = x / norm
!     point_coord2(2,i) = y / norm
!     point_coord2(3,i) = z / norm
!
!   end do
! !
! !  Now build the face in the dual associated with each node IFACE.
! !
!   do iface = 1, face_num2
! !
! !  Initialize the order.
! !
!     face_order2 = 0
! !
! !  Find the first occurrence of IFACE in an edge of polyhedron 1.
! !
!     call icol_find_item ( face_point1, max_order, face_order1, face_num1, &
!       iface, irow, icol )
!
!     if ( irow == 0 ) then
!       write ( *, * ) ' '
!       write ( *, * ) 'DUAL_SHAPE_3D - Fatal error!'
!       write ( *, * ) '  Could not find an edge using node ', iface
!       return
!     end if
! !
! !  Save the following node as ISTOP.
! !  When we encounter ISTOP again, this will mark the end of our search.
! !
!     i = irow + 1
!     if ( i > face_order1 ) then
!       i = 1
!     end if
!
!     istop = face_point1(i,icol)
! !
! !  Save the previous node as INEXT.
! !
!     do
!
!       i = irow - 1
!       if ( i < 1 ) then
!         i = i + face_order1
!       end if
!
!       inext = face_point1(i,icol)
!
!       face_order2 = face_order2 + 1
!
!       face_point2(face_order2,iface) = icol
! !
! !  If INEXT =/= ISTOP, continue.
! !
!       if ( inext == istop ) then
!         exit
!       end if
! !
! !  Set IPREV:= INEXT.
! !
!       iprev = inext
! !
! !  Search for the occurrence of the edge IFACE-IPREV.
! !
!       call icol_find_pair_wrap ( face_point1, max_order, face_order1, &
!         face_num1, iface, iprev, irow, icol )
!
!       if ( irow == 0 ) then
!         write ( *, * ) ' '
!         write ( *, * ) 'DUAL_SHAPE_3D - Fatal error!'
!         write ( *, * ) '  No edge from node ', iprev
!         write ( *, * ) '  to node ', iface
!         return
!       end if
!
!     end do
!
!   end do
!
!   return
! end subroutine
!
! subroutine ellipse_area_2d ( r1, r2, area )
! !
! !*******************************************************************************
! !
! !! ELLIPSE_AREA_2D returns the area of an ellipse in 2D.
! !
! !
! !  Modified:
! !
! !    28 May 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R1, R2, the "radius" of the ellipse in the major
! !    and minor axis directions.  A circle has these values equal.
! !
! !    Output, real AREA, the area of the ellipse.
! !
!   real area
!   real r1
!   real r2
! !
!   area = pi* r1 * r2
!
!   return
! end subroutine
! subroutine ellipse_points_2d ( x0, y0, r1, r2, psi, n, x, y )
! !
! !*******************************************************************************
! !
! !! ELLIPSE_POINTS_2D returns N points on an ellipse in 2D.
! !
! !
! !  Discussion:
! !
! !    The points are "equally spaced" in the angular sense.  They are
! !    not equally spaced along the perimeter of the ellipse.
! !
! !  Modified:
! !
! !    28 May 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, the coordinates of the center of the ellipse.
! !
! !    Input, real R1, R2, the "radius" of the ellipse in the major
! !    and minor axis directions.  A circle has these values equal.
! !
! !    Input, real PSI, the angle that the major axis of the ellipse
! !    makes with the X axis.  A value of 0.0 means that the major and
! !    minor axes of the ellipse will be the X and Y coordinate axes.
! !
! !    Input, integer N, the number of points desired.  N must be at least 1.
! !
! !    Output, real X(N), Y(N), the coordinates of points on the ellipse.
! !
!   integer n
! !
!   integer i
!   real psi
!   real r1
!   real r2
!   real theta
!   real x0
!   real x(n)
!   real y0
!   real y(n)
! !
!   do i = 1, n
!
!     theta = ( 2.0E+00 * pi * real ( i - 1 ) ) / real ( n )
!
!     x(i) = x0 + r1 * cos ( psi ) * cos ( theta ) &
!               - r2 * sin ( psi ) * sin ( theta )
!
!     y(i) = y0 + r1 * sin ( psi ) * cos ( theta ) &
!               + r2 * cos ( psi ) * sin ( theta )
!
!   end do
!
!   return
! end subroutine
! subroutine ellipse_points_arc_2d ( x0, y0, r1, r2, psi, theta1, theta2, n, &
!   x, y )
! !
! !*******************************************************************************
! !
! !! ELLIPSE_POINTS_ARC_2D returns N points on an elliptical arc in 2D.
! !
! !
! !  Discussion:
! !
! !    The points are "equally spaced" in the angular sense.  They are
! !    not equally spaced along the perimeter of the ellipse.
! !
! !  Modified:
! !
! !    29 May 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, the coordinates of the center of the ellipse.
! !
! !    Input, real R1, R2, the "radius" of the ellipse in the major
! !    and minor axis directions.  A circle has these values equal.
! !
! !    Input, real PSI, the angle that the major axis of the ellipse
! !    makes with the X axis.  A value of 0.0 means that the major and
! !    minor axes of the ellipse will be the X and Y coordinate axes.
! !
! !    Input, real THETA1, THETA2, the angular coordinates of the first
! !    and last points to be drawn, in radians.  This angle is measured
! !    with respect to the (possibly tilted) major axis.
! !
! !    Input, integer N, the number of points desired.  N must be at least 1.
! !
! !    Output, real X(N), Y(N), the coordinates of points on the ellipse.
! !
!   integer n
! !
!   integer i
!   real psi
!   real r1
!   real r2
!   real theta
!   real theta1
!   real theta2
!   real theta3
!   real x0
!   real x(n)
!   real y0
!   real y(n)
! !
! !  THETA3 is the smallest angle, no less than THETA1, which
! !  coincides with THETA2.
! !
!   theta3 = theta1 + r_modp ( theta2 - theta1, 2.0E+00 * real(pi) )
!
!   do i = 1, n
!
!     if ( n > 1 ) then
!       theta = ( real ( n - i ) * theta1 + real ( i - 1 ) * theta3 ) &
!         / real ( n - 1 )
!     else
!       theta = 0.5E+00 * ( theta1 + theta3 )
!     end if
!
!     x(i) = x0 + r1 * cos ( psi ) * cos ( theta ) &
!               - r2 * sin ( psi ) * sin ( theta )
!
!     y(i) = y0 + r1 * sin ( psi ) * cos ( theta ) &
!               + r2 * cos ( psi ) * sin ( theta )
!
!   end do
!
!   return
! end subroutine
! function enorm_2d ( x1, y1 )
! !
! !*******************************************************************************
! !
! !! ENORM_2D computes the Euclidean norm of a vector in 2D.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, the coordinates of the vector.
! !
! !    Output, real ENORM_2D, the Euclidean norm of the vector.
! !
!   real enorm_2d
!   real x1
!   real y1
! !
!   enorm_2d = sqrt ( x1 * x1 + y1 * y1 )
!
!   return
! end function
! function enorm_3d ( x1, y1, z1 )
! !
! !*******************************************************************************
! !
! !! ENORM_3D computes the Euclidean norm of a vector in 3D.
! !
! !
! !  Modified:
! !
! !    27 July 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, the coordinates of the vector.
! !
! !    Output, real ENORM_3D, the Euclidean norm of the vector.
! !
!   real enorm_3d
!   real x1
!   real y1
!   real z1
! !
!   enorm_3d = sqrt ( x1 * x1 + y1 * y1 + z1 * z1 )
!
!   return
! end function
! function enorm_nd ( n, x )
! real enorm_nd
! !
! !*******************************************************************************
! !
! !! ENORM_ND computes the Euclidean norm of a vector in ND.
! !
! !
! !  Modified:
! !
! !    27 July 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Input, real X(N), the coordinates of the vector.
! !
! !    Output, real ENORM_ND, the Euclidean norm of the vector.
! !
!   integer n
! !
!   integer i
!   real x(n)
! !
!   enorm_nd = sqrt ( sum ( x(1:n)**2 ) )
!
!   return
! end function
!
!
! function enorm0_2d ( x0, y0, x1, y1 )
! real enorm0_2d
! !
! !*******************************************************************************
! !
! !! ENORM0_2D computes the Euclidean norm of (P1-P0) in 2D.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, X1, Y1, the coordinates of the points P0 and P1.
! !
! !    Output, real ENORM0_2D, the Euclidean norm of (P1-P0).
! !
!   real x0
!   real x1
!   real y0
!   real y1
! !
!   enorm0_2d = sqrt ( ( x1 - x0 )**2 + ( y1 - y0 )**2 )
!
!   return
! end function
!
!
! function enorm0_3d ( x0, y0, z0, x1, y1, z1 )
! real enorm0_3d
! !
! !*******************************************************************************
! !
! !! ENORM0_3D computes the Euclidean norm of (P1-P0) in 3D.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, Z0, X1, Y1, Z1, the coordinates of the points
! !    P0 and P1.
! !
! !    Output, real ENORM0_3D, the Euclidean norm of (P1-P0).
! !
!   real x0
!   real x1
!   real y0
!   real y1
!   real z0
!   real z1
! !
!   enorm0_3d = sqrt ( ( x1 - x0 )**2 + ( y1 - y0 )**2 + ( z1 - z0 )**2 )
!
!   return
! end function
!
!
! function enorm0_nd ( n, x0, x1 )
! real enorm0_nd
! !
! !*******************************************************************************
! !
! !! ENORM0_ND computes the Euclidean norm of (P1-P0) in ND.
! !
! !
! !  Modified:
! !
! !    18 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Input, real X0(N), the base vector.
! !
! !    Input, real X1(N), the displacement vector.
! !
! !    Output, real ENORM0_ND, the Euclidean norm of the vector
! !    ( X1 - X0 ).
! !
!   integer n
! !
!   real x0(n)
!   real x1(n)
! !
!   enorm0_nd = sqrt ( sum ( ( x1(1:n) - x0(1:n) )**2 ) )
!
!   return
! end function
!
!
! function enormsq0_2d ( x0, y0, x1, y1 )
! real enormsq0_2d
! !
! !*******************************************************************************
! !
! !! ENORMSQ0_2D computes the square of the Euclidean norm of (P1-P0) in 2D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, X1, Y1, the coordinates of the points
! !    P0 and P1.
! !
! !    Output, real ENORMSQ0_2D, the square of the Euclidean norm of (P1-P0).
! !
!   real x0
!   real x1
!   real y0
!   real y1
! !
!   enormsq0_2d = ( x1 - x0 )**2 + ( y1 - y0 )**2
!
!   return
! end function
! function enormsq0_3d ( x0, y0, z0, x1, y1, z1 )
! !
! !*******************************************************************************
! !
! !! ENORMSQ0_3D computes the square of the Euclidean norm of (P1-P0) in 3D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X0, Y0, Z0, X1, Y1, Z1, the coordinates of the points
! !    P0 and P1.
! !
! !    Output, real ENORMSQ0_3D, the square of the Euclidean norm of (P1-P0).
! !
!   real enormsq0_3d
!   real x0
!   real x1
!   real y0
!   real y1
!   real z0
!   real z1
! !
!   enormsq0_3d =( x1 - x0 )**2 + ( y1 - y0 )**2 + ( z1 - z0 )**2
!
!   return
! end function
!
! function enormsq0_nd ( n, x0, x1 )
!
! !
! !*******************************************************************************
! !
! !! ENORMSQ0_ND computes the squared Euclidean norm of (P1-P0) in ND.
! !
! !
! !  Modified:
! !
! !    18 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Input, real X0(N), the coordinates of the base vector.
! !
! !    Input, real X1(N), the coordinates of the displacement vector.
! !
! !    Output, real ENORMSQ0_ND, the squared Euclidean norm of the vector
! !    ( X1 - X0 ).
! !
!   integer n
! !
!   real enormsq0_nd
!   real x0(n)
!   real x1(n)
! !
!   enormsq0_nd = sum ( ( x1(1:n) - x0(1:n) )**2 )
!
!   return
! end function
! subroutine glob2loc_3d ( cospitch, cosroll, cosyaw, locpts, globas, glopts, &
!   sinpitch, sinroll, sinyaw )
! !
! !*******************************************************************************
! !
! !! GLOB2LOC_3D converts from a global to a local coordinate system in 3D.
! !
! !
! !  Discussion:
! !
! !    A global coordinate system is given.
! !
! !    A local coordinate system has been translated to the point with
! !    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
! !    a roll.
! !
! !    A point has global coordinates GLOPTS, and it is desired to know
! !    the point's local coordinates LOCPTS.
! !
! !    The transformation may be written as
! !
! !      LOC = M_ROLL * M_PITCH * M_YAW * ( GLOB - GLOBAS )
! !
! !    where
! !
! !               (       1            0            0      )
! !    M_ROLL =   (       0        cos(Roll)    sin(Roll)  )
! !               (       0      - sin(Roll)    cos(Roll)  )
! !
! !               (   cos(Pitch)       0      - sin(Pitch) )
! !    M_PITCH =  (       0            1            0      )
! !               (   sin(Pitch)       0        cos(Pitch) )
! !
! !               (   cos(Yaw)     sin(Yaw)         0      )
! !    M_YAW    = ( - sin(Yaw)     cos(Yaw)         0      )
! !               (       0            0            1      )
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real COSPITCH, COSROLL, COSYAW, the cosines of the pitch,
! !    roll and yaw angles.
! !
! !    Input, real GLOBAS(3), the X, Y, and Z coordinates of the
! !    global base vector.
! !
! !    Input, real GLOPTS(3), the global ( X, Y and Z ) coordinates of
! !    the point whose coordinates are to be transformed.
! !
! !    Output, real LOCPTS(3), the local coordinates of the point
! !    whose global coordinates were given in GLOPTS.
! !
! !    Input, real SINPITCH, SINROLL, SINYAW, the sines of the pitch,
! !    roll and yaw angles.
! !
!   real cospitch
!   real cosroll
!   real cosyaw
!   real globas(3)
!   real glopts(3)
!   real locpts(3)
!   real sinpitch
!   real sinroll
!   real sinyaw
! !
!   locpts(1) = ( cosyaw * cospitch ) * ( glopts(1) - globas(1) ) &
!             + ( sinyaw * cospitch ) * ( glopts(2) - globas(2) ) &
!             -   sinpitch * ( glopts(3) - globas(3) )
!
!   locpts(2) = ( cosyaw * sinpitch * sinroll - sinyaw * cosroll ) &
!     * ( glopts(1) - globas(1) ) &
!     + ( sinyaw * sinpitch * sinroll + cosyaw * cosroll ) &
!     * ( glopts(2) - globas(2) ) &
!     +   cospitch * sinroll * ( glopts(3) - globas(3) )
!
!   locpts(3) = ( cosyaw * sinpitch * cosroll + sinyaw * sinroll ) &
!     * ( glopts(1) - globas(1) ) &
!     + ( sinyaw * sinpitch * cosroll - cosyaw * sinroll  ) &
!     * ( glopts(2) - globas(2) ) &
!     + ( cospitch * cosroll ) * ( glopts(3) - globas(3) )
!
!   return
! end subroutine
! subroutine halfspace_imp_triangle_int_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   a, b, c, d, num_int, x, y, z )
! !
! !*******************************************************************************
! !
! !! HALFSPACE_IMP_TRIANGLE_INT_3D: intersection ( implicit halfspace, triangle ) in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a half-space in 3D may be described as the set
! !    of points (X,Y,Z) on or "above" an implicit plane:
! !
! !      A * X + B * Y + C * Z + D >= 0
! !
! !    The triangle is specified by listing its three vertices.
! !
! !  Discussion:
! !
! !    The intersection may be described by the number of vertices of the
! !    triangle that are included in the halfspace, and by the location of
! !    points between vertices that separate a side of the triangle into
! !    an included part and an unincluded part.
! !
! !    0 vertices, 0 separators    (no intersection)
! !    1 vertex, 0 separators      (point intersection)
! !    2 vertices, 0 separators    (line intersection)
! !    3 vertices, 0 separators    (triangle intersection)
! !
! !    1 vertex, 2 separators,     (intersection is a triangle)
! !    2 vertices, 2 separators,   (intersection is a quadrilateral).
! !
! !  Modified:
! !
! !    02 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of the vertices of the triangle.
! !
! !    Input, real A, B, C, D, the parameters that define the implicit plane,
! !    which in turn define the implicit halfspace.
! !
! !    Output, integer NUM_INT, the number of intersection points returned,
! !    which will always be between 0 and 4.
! !
! !    Output, real X(4), Y(4), Z(4), the coordinates of the NUM_INT
! !    intersection points.  The points will lie in sequence on the triangle.
! !    Some points will be vertices, and some may be separators.
! !
!   real a
!   real b
!   real c
!   real d
!   real dist1
!   real dist2
!   real dist3
!   integer num_int
!   real x(4)
!   real x1
!   real x2
!   real x3
!   real y(4)
!   real y1
!   real y2
!   real y3
!   real z(4)
!   real z1
!   real z2
!   real z3
! !
! !  Compute the signed distances between the vertices and the plane.
! !
!   dist1 = a * x1 + b * y1 + c * z1 + d
!   dist2 = a * x2 + b * y2 + c * z2 + d
!   dist3 = a * x3 + b * y3 + c * z3 + d
! !
! !  Now we can find the intersections.
! !
!   call halfspace_triangle_int_3d ( x1, y1, z1, x2, y2, &
!     z2, x3, y3, z3, dist1, dist2, dist3, num_int, x, y, z )
!
!   return
! end subroutine
! subroutine halfspace_norm_triangle_int_3d ( x1, y1, z1, x2, y2, z2, &
!   x3, y3, z3, xp, yp, zp, xn, yn, zn, num_int, x, y, z )
! !
! !*******************************************************************************
! !
! !! HALFSPACE_NORM_TRIANGLE_INT_3D: intersection ( normal halfspace, triangle ) in 3D.
! !
! !
! !  Definition:
! !
! !    The normal form of a halfspace in 3D may be described as the set
! !    of points (X,Y,Z) on or "above" a plane described in normal form:
! !
! !      ( Xp, Yp, Zp ) is a point on the plane,
! !      ( Xn, Yn, Zn ) is the unit normal vector, pointing "out" of the halfspace.
! !
! !    The triangle is specified by listing its three vertices.
! !
! !  Discussion:
! !
! !    The intersection may be described by the number of vertices of the
! !    triangle that are included in the halfspace, and by the location of
! !    points between vertices that separate a side of the triangle into
! !    an included part and an unincluded part.
! !
! !    0 vertices, 0 separators    (no intersection)
! !    1 vertex, 0 separators      (point intersection)
! !    2 vertices, 0 separators    (line intersection)
! !    3 vertices, 0 separators    (triangle intersection)
! !
! !    1 vertex, 2 separators,     (intersection is a triangle)
! !    2 vertices, 2 separators,   (intersection is a quadrilateral).
! !
! !  Modified:
! !
! !    03 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of the vertices of the triangle.
! !
! !    Input, real XP, YP, ZP, a point on the bounding plane that defines
! !    the halfspace.
! !
! !    Input, real XN, YN, ZN, the components of the normal vector to the
! !    bounding plane that defines the halfspace.  By convention, the
! !    normal vector points "outwards" from the halfspace.
! !
! !    Output, integer NUM_INT, the number of intersection points returned,
! !    which will always be between 0 and 4.
! !
! !    Output, real X(4), Y(4), Z(4), the coordinates of the NUM_INT
! !    intersection points.  The points will lie in sequence on the triangle.
! !    Some points will be vertices, and some may be separators.
! !
!   real d
!   real dist1
!   real dist2
!   real dist3
!   integer num_int
!   real x(4)
!   real x1
!   real x2
!   real x3
!   real xn
!   real xp
!   real y(4)
!   real y1
!   real y2
!   real y3
!   real yn
!   real yp
!   real z(4)
!   real z1
!   real z2
!   real z3
!   real zn
!   real zp
! !
! !  Compute the signed distances between the vertices and the plane.
! !
!   d = - xn * xp - yn * yp - zn * zp
!
!   dist1 = xn * x1 + yn * y1 + zn * z1 + d
!   dist2 = xn * x2 + yn * y2 + zn * z2 + d
!   dist3 = xn * x3 + yn * y3 + zn * z3 + d
! !
! !  Now we can find the intersections.
! !
!   call halfspace_triangle_int_3d ( x1, y1, z1, x2, y2, &
!     z2, x3, y3, z3, dist1, dist2, dist3, num_int, x, y, z )
!
!   return
! end subroutine
! subroutine halfspace_triangle_int_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   dist1, dist2, dist3, num_int, x, y, z )
! !
! !*******************************************************************************
! !
! !! HALFSPACE_TRIANGLE_INT_3D: intersection ( halfspace, triangle ) in 3D.
! !
! !
! !  Definition:
! !
! !    The triangle is specified by listing its three vertices.
! !
! !  Discussion:
! !
! !    The halfspace is not described in the input data.  Rather, the
! !    distances from the triangle vertices to the halfspace are given.
! !
! !    The intersection may be described by the number of vertices of the
! !    triangle that are included in the halfspace, and by the location of
! !    points between vertices that separate a side of the triangle into
! !    an included part and an unincluded part.
! !
! !    0 vertices, 0 separators    (no intersection)
! !    1 vertex, 0 separators      (point intersection)
! !    2 vertices, 0 separators    (line intersection)
! !    3 vertices, 0 separators    (triangle intersection)
! !
! !    1 vertex, 2 separators,     (intersection is a triangle)
! !    2 vertices, 2 separators,   (intersection is a quadrilateral).
! !
! !  Modified:
! !
! !    03 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of the vertices of the triangle.
! !
! !    Input, real DIST1, DIST2, DIST3, the distances from each of the
! !    three vertices of the triangle to the halfspace.  The distance is
! !    zero if a vertex lies within the halfspace, or on the plane that
! !    defines the boundary of the halfspace.  Otherwise, it is the
! !    distance from that vertex to the bounding plane.
! !
! !    Output, integer NUM_INT, the number of intersection points returned,
! !    which will always be between 0 and 4.
! !
! !    Output, real X(4), Y(4), Z(4), the coordinates of the NUM_INT
! !    intersection points.  The points will lie in sequence on the triangle.
! !    Some points will be vertices, and some may be separators.
! !
!   real dist1
!   real dist2
!   real dist3
!   integer num_int
!   real x(4)
!   real x1
!   real x2
!   real x3
!   real y(4)
!   real y1
!   real y2
!   real y3
!   real z(4)
!   real z1
!   real z2
!   real z3
! !
! !  Walk around the triangle, looking for vertices that are included,
! !  and points of separation.
! !
!   num_int = 0
!
!   if ( dist1 <= 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x1
!     y(num_int) = y1
!     z(num_int) = z1
!
!   end if
!
!   if ( dist1 * dist2 < 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = ( dist1 * x2 - dist2 * x1 ) / ( dist1 - dist2 )
!     y(num_int) = ( dist1 * y2 - dist2 * y1 ) / ( dist1 - dist2 )
!     z(num_int) = ( dist1 * z2 - dist2 * z1 ) / ( dist1 - dist2 )
!
!   end if
!
!   if ( dist2 <= 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x2
!     y(num_int) = y2
!     z(num_int) = z2
!
!   end if
!
!   if ( dist2 * dist3 < 0.0E+00 ) then
!
!     num_int = num_int + 1
!
!     x(num_int) = ( dist2 * x3 - dist3 * x2 ) / ( dist2 - dist3 )
!     y(num_int) = ( dist2 * y3 - dist3 * y2 ) / ( dist2 - dist3 )
!     z(num_int) = ( dist2 * z3 - dist3 * z2 ) / ( dist2 - dist3 )
!
!   end if
!
!   if ( dist3 <= 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x3
!     y(num_int) = y3
!     z(num_int) = z3
!
!   end if
!
!   if ( dist3 * dist1 < 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = ( dist3 * x1 - dist1 * x3 ) / ( dist3 - dist1 )
!     y(num_int) = ( dist3 * y1 - dist1 * y3 ) / ( dist3 - dist1 )
!     z(num_int) = ( dist3 * z1 - dist1 * z3 ) / ( dist3 - dist1 )
!
!   end if
!
!   return
! end subroutine
! subroutine helix_shape_3d ( a, n, r, theta1, theta2, x, y, z )
! !
! !*******************************************************************************
! !
! !! HELIX_SHAPE_3D computes points on a helix in 3D.
! !
! !
! !  Discussion:
! !
! !    The user specifies the parameters A and R, the first and last
! !    THETA values, and the number of equally spaced THETA values
! !    at which (X,Y,Z) values are to be computed.
! !
! !  Formula:
! !
! !    X = R * COS ( THETA )
! !    Y = R * SIN ( THETA )
! !    Z = A * THETA
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, the rate at which Z advances with THETA.
! !
! !    Input, integer N, the number of points to compute on the helix.
! !
! !    Input, real R, the radius of the helix.
! !
! !    Input, real THETA1, THETA2, the first and last THETA values at
! !    which to compute points on the helix.  THETA is measured in
! !    radians.
! !
! !    Output, real X(N), Y(N), Z(N), the X, Y and Z coordinates
! !    of points on the helix.
! !
!   integer n
! !
!   real a
!   integer i
!   real r
!   real theta
!   real theta1
!   real theta2
!   real x(n)
!   real y(n)
!   real z(n)
! !
!   do i = 1, n
!
!     if ( n == 1 ) then
!       theta = 0.5E+00 * ( theta1 + theta2 )
!     else
!       theta = ( real ( n - i ) * theta1 + real ( i - 1 ) * theta2 ) &
!         / real ( n - 1 )
!     end if
!
!     x(i) = r * cos ( theta )
!     y(i) = r * sin ( theta )
!     z(i) = a * theta
!
!   end do
!
!   return
! end subroutine
! function hexagon_area_2d ( rad )
! !
! !*******************************************************************************
! !
! !! HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
! !
! !
! !  Modified:
! !
! !    16 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real RAD, the radius of the hexagon.
! !
! !    Output, real HEXAGON_AREA_2D, the area of the hexagon.
! !
!   real hexagon_area_2d
!   real hexagon_unit_area_2d
!   real rad
! !
!   hexagon_area_2d = rad**2 * hexagon_unit_area_2d ( )
!
!   return
! end function
! function hexagon_unit_area_2d ( )
! !
! !*******************************************************************************
! !
! !! HEXAGON_UNIT_AREA_2D returns the area of a unit regular hexagon in 2D.
! !
! !
! !  Modified:
! !
! !    07 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real HEXAGON_UNIT_AREA_2D, the area of the hexagon.
! !
!   real hexagon_unit_area_2d
!   real rad
! !
!   hexagon_unit_area_2d = 3.0D+00 * sqrt ( 3.0D+00 ) / 2.0D+00
!
!   return
! end function
! subroutine hexagon_shape_2d ( angle, x, y )
! !
! !*******************************************************************************
! !
! !! HEXAGON_SHAPE_2D returns points on the unit hexagon in 2D.
! !
! !
! !  Diagram:
! !
! !      120_____60
! !        /     \
! !    180/       \0
! !       \       /
! !        \_____/
! !      240     300
! !
! !  Discussion:
! !
! !    The unit hexagon has maximum radius 1, and is the hull of the points
! !
! !      (   1,              0 ),
! !      (   0.5,   sqrt (3)/2 ),
! !      ( - 0.5,   sqrt (3)/2 ),
! !      ( - 1,              0 ),
! !      ( - 0.5, - sqrt (3)/2 ),
! !      (   0.5, - sqrt (3)/2 ).
! !
! !  Modified:
! !
! !    12 July 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real ANGLE, the angle, in degrees, of the point.
! !
! !    Output, real X, Y, the coordinates of the point.
! !
!   real angle
!   real angle2
!   real x
!   real y
! !
!   angle2 = mod ( angle, 360.0E+00 )
!   if ( angle2 < 0.0E+00 ) then
!     angle2 = angle2 + 360.0E+00
!   end if
! !
! !  y = - sqrt(3) * x + sqrt(3)
! !
!   if ( 0.0E+00 <= angle2 .and. angle2 <= 60.0E+00 ) then
!
!     x = sqrt ( 3.0E+00 ) / ( tand ( angle2 ) + sqrt ( 3.0E+00 ) )
!     y = tand ( angle2 ) * x
! !
! !  y = sqrt(3) / 2
! !
!   else if ( angle2 <= 120.0E+00 ) then
!
!     y = sqrt ( 3.0E+00 ) / 2.0E+00
!     x = cotd ( angle2 ) * y
! !
! !  y = sqrt(3) * x + sqrt(3)
! !
!   else if ( angle2 <= 180.0E+00 ) then
!
!     x = sqrt ( 3.0E+00 ) / ( tand ( angle2 ) - sqrt ( 3.0E+00 ) )
!     y = tand ( angle2 ) * x
! !
! !  y = - sqrt(3) * x - sqrt(3)
! !
!   else if ( angle2 <= 240.0E+00 ) then
!
!     x = - sqrt ( 3.0E+00 ) / ( tand ( angle2 ) + sqrt ( 3.0E+00 ) )
!     y = tand ( angle2 ) * x
! !
! !  y = - sqrt(3) / 2
! !
!   else if ( angle2 <= 300.0E+00 ) then
!
!     y = - sqrt ( 3.0E+00 ) / 2.0E+00
!     x = y * cotd ( angle2 )
! !
! !  y = sqrt(3) * x - sqrt(3)
! !
!   else if ( angle2 <= 360.0E+00 ) then
!
!     x = - sqrt ( 3.0E+00 ) / ( tand ( angle2 ) - sqrt ( 3.0E+00 ) )
!     y = tand ( angle2 ) * x
!
!   end if
!
!   return
! end subroutine
! subroutine i_random ( ilo, ihi, i )
! !
! !*******************************************************************************
! !
! !! I_RANDOM returns a random integer in a given range.
! !
! !
! !  Modified:
! !
! !    23 September 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer ILO, IHI, the minimum and maximum acceptable values.
! !
! !    Output, integer I, the randomly chosen integer.
! !
!   logical, save :: seed = .false.
!   integer i
!   integer ihi
!   integer ilo
!   real r
!   real rhi
!   real rlo
! !
!   if ( .not. seed ) then
!     call random_seed
!     seed = .true.
!   end if
! !
! !  Pick a random number in (0,1).
! !
!   call random_number ( harvest = r )
! !
! !  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
! !  each with a "neighborhood" of width 1.
! !
!   rlo = real ( ilo ) - 0.5E+00
!   rhi = real ( ihi ) + 0.5E+00
! !
! !  Set I to the integer that is nearest the scaled value of R.
! !
!   i = nint ( ( 1.0E+00 - r ) * rlo + r * rhi )
! !
! !  In case of oddball events at the boundary, enforce the limits.
! !
!   i = max ( i, ilo )
!   i = min ( i, ihi )
!
!   return
! end subroutine
! subroutine i_swap ( i, j )
! !
! !*******************************************************************************
! !
! !! I_SWAP switches two integer values.
! !
! !
! !  Modified:
! !
! !    30 November 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, integer I, J.  On output, the values of I and
! !    J have been interchanged.
! !
!   integer i
!   integer j
!   integer k
! !
!   k = i
!   i = j
!   j = k
!
!   return
! end subroutine
! subroutine icol_find_item ( itable, maxrow, nrow, ncol, item, irow, icol )
! !
! !*******************************************************************************
! !
! !! ICOL_FIND_ITEM searches a table by columns for a given value.
! !
! !
! !  Modified:
! !
! !    11 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer ITABLE(MAXROW,NCOL), the table to search.
! !
! !    Input, integer MAXROW, the leading dimension of ITABLE.
! !
! !    Input, integer NROW, NCOL, the number of rows and columns
! !    in the table.
! !
! !    Input, integer ITEM, the value to search for.
! !
! !    Output, integer IROW, ICOL, the row and column indices
! !    of the first occurrence of the value ITEM.  The search
! !    is conducted by rows.  If the item is not found, then
! !    IROW = ICOL = 0.
! !
!   integer maxrow
!   integer ncol
! !
!   integer i
!   integer icol
!   integer irow
!   integer itable(maxrow,ncol)
!   integer item
!   integer nrow
!   integer j
! !
!   do j = 1, ncol
!     do i = 1, nrow
!       if ( itable(i,j) == item ) then
!         irow = i
!         icol = j
!         return
!       end if
!     end do
!   end do
!
!   irow = 0
!   icol = 0
!
!   return
! end subroutine
! subroutine icol_find_pair_wrap ( itable, maxrow, nrow, ncol, item1, item2, &
!   irow, icol )
! !
! !*******************************************************************************
! !
! !! ICOL_FIND_ITEM wrap searches a table by columns for a pair of items.
! !
! !
! !  Discussion:
! !
! !    The items must occur consecutively, with ITEM1 occurring
! !    first.  However, wrapping is allowed.  That is, if ITEM1
! !    occurs in the last row, and ITEM2 in the first, this
! !    is also regarded as a match.
! !
! !  Modified:
! !
! !    11 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer ITABLE(MAXROW,NCOL), the table to search.
! !
! !    Input, integer MAXROW, the leading dimension of ITABLE.
! !
! !    Input, integer NROW, NCOL, the number of rows and columns
! !    in the table.
! !
! !    Input, integer ITEM1, ITEM2, the values to search for.
! !
! !    Output, integer IROW, ICOL, the row and column indices
! !    of the first occurrence of the value ITEM1 followed immediately
! !    by ITEM2.  The search is conducted by columns.  If the pair of
! !    items is not found, then IROW = ICOL = 0.  If IROW = NROW,
! !    the ITEM1 occurs in row NROW and ITEM2 occurs in row 1.
! !
!   integer maxrow
!   integer ncol
! !
!   integer i
!   integer icol
!   integer ip1
!   integer irow
!   integer itable(maxrow,ncol)
!   integer item1
!   integer item2
!   integer j
!   integer nrow
! !
!   do j = 1, ncol
!     do i = 1, nrow
!
!       if ( itable(i,j) == item1 ) then
!
!         if ( i < nrow ) then
!           ip1 = i + 1
!         else
!           ip1 = 1
!         end if
!
!         if ( itable(ip1,j) == item2 ) then
!           irow = i
!           icol = j
!           return
!         end if
!
!       end if
!
!     end do
!   end do
!
!   irow = 0
!   icol = 0
!
!   return
! end subroutine
! subroutine icos_shape_3d ( max_num, max_order, point_num, face_num, &
!   face_order, point_coord, face_point )
! !
! !*******************************************************************************
! !
! !! ICOS_SHAPE_3D describes an icosahedron in 3D.
! !
! !
! !  Discussion:
! !
! !    The vertices of the icosahedron lie on the unit sphere.
! !
! !  Modified:
! !
! !    12 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAX_NUM, the maximum number of faces, and of
! !    points.  This quantity is used to dimension arrays.
! !    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
! !    not be set.
! !
! !    Input, integer MAX_ORDER, the maximum order of any face.
! !    This quantity is used to dimension arrays.  If FACE_ORDER
! !    exceeds MAX_ORDER, the arrays will not be set.
! !
! !    Output, integer POINT_NUM, the number of points in the shape.
! !
! !    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
! !    the X, Y and Z coordinates of point J.
! !
! !    Output, integer FACE_NUM, the number of faces in the shape.
! !
! !    Output, integer FACE_ORDER, the number of vertices per face.
! !
! !    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
! !    contains the index of the I-th point in the J-th face.  The
! !    points are listed in the counter-clockwise direction defined
! !    by the outward normal at the face.
! !
!   integer max_num
!   integer max_order
! !
!   real a
!   real b
!   integer face_num
!   integer face_order
!   integer face_point(max_order,max_num)
!   real phi
!   integer point_num
!   real point_coord(3,max_num)
! !
!   point_num = 12
!   face_num = 20
!   face_order = 3
! !
! !  Check.
! !
!   if ( point_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'ICOS_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of vertices exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'ICOS_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of faces exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_order > max_order ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'ICOS_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Face order exceeds MAX_ORDER.'
!     return
!   end if
! !
! !  Set point coordinates.
! !
!   phi = 0.5E+00 * ( sqrt ( 5.0E+00 ) + 1.0E+00 )
!   b = 1.0E+00 / sqrt ( 1.0E+00 + phi * phi )
!   a = phi * b
!
!   point_coord(1,1) =   a
!   point_coord(2,1) =   b
!   point_coord(3,1) = 0.0E+00
!
!   point_coord(1,2) = - a
!   point_coord(2,2) =   b
!   point_coord(3,2) = 0.0E+00
!
!   point_coord(1,3) =   a
!   point_coord(2,3) = - b
!   point_coord(3,3) = 0.0E+00
!
!   point_coord(1,4) = - a
!   point_coord(2,4) = - b
!   point_coord(3,4) = 0.0E+00
!
!   point_coord(1,5) =   b
!   point_coord(2,5) = 0.0E+00
!   point_coord(3,5) =   a
!
!   point_coord(1,6) =   b
!   point_coord(2,6) = 0.0E+00
!   point_coord(3,6) = - a
!
!   point_coord(1,7) = - b
!   point_coord(2,7) = 0.0E+00
!   point_coord(3,7) =   a
!
!   point_coord(1,8) = - b
!   point_coord(2,8) = 0.0E+00
!   point_coord(3,8) = - a
!
!   point_coord(1,9) =  0.0E+00
!   point_coord(2,9) =    a
!   point_coord(3,9) =    b
!
!   point_coord(1,10) = 0.0E+00
!   point_coord(2,10) = - a
!   point_coord(3,10) =   b
!
!   point_coord(1,11) = 0.0E+00
!   point_coord(2,11) =   a
!   point_coord(3,11) = - b
!
!   point_coord(1,12) = 0.0E+00
!   point_coord(2,12) = - a
!   point_coord(3,12) = - b
! !
! !  Set faces.
! !
!   face_point(1,1) = 1
!   face_point(2,1) = 9
!   face_point(3,1) = 5
!
!   face_point(1,2) = 1
!   face_point(2,2) = 6
!   face_point(3,2) = 11
!
!   face_point(1,3) = 3
!   face_point(2,3) = 5
!   face_point(3,3) = 10
!
!   face_point(1,4) = 3
!   face_point(2,4) = 12
!   face_point(3,4) = 6
!
!   face_point(1,5) = 2
!   face_point(2,5) = 7
!   face_point(3,5) = 9
!
!   face_point(1,6) = 2
!   face_point(2,6) = 11
!   face_point(3,6) = 8
!
!   face_point(1,7) = 4
!   face_point(2,7) = 10
!   face_point(3,7) = 7
!
!   face_point(1,8) = 4
!   face_point(2,8) = 8
!   face_point(3,8) = 12
!
!   face_point(1,9) = 1
!   face_point(2,9) = 11
!   face_point(3,9) = 9
!
!   face_point(1,10) = 2
!   face_point(2,10) = 9
!   face_point(3,10) = 11
!
!   face_point(1,11) = 3
!   face_point(2,11) = 10
!   face_point(3,11) = 12
!
!   face_point(1,12) = 4
!   face_point(2,12) = 12
!   face_point(3,12) = 10
!
!   face_point(1,13) = 5
!   face_point(2,13) = 3
!   face_point(3,13) = 1
!
!   face_point(1,14) = 6
!   face_point(2,14) = 1
!   face_point(3,14) = 3
!
!   face_point(1,15) = 7
!   face_point(2,15) = 2
!   face_point(3,15) = 4
!
!   face_point(1,16) = 8
!   face_point(2,16) = 4
!   face_point(3,16) = 2
!
!   face_point(1,17) = 9
!   face_point(2,17) = 7
!   face_point(3,17) = 5
!
!   face_point(1,18) = 10
!   face_point(2,18) = 5
!   face_point(3,18) = 7
!
!   face_point(1,19) = 11
!   face_point(2,19) = 6
!   face_point(3,19) = 8
!
!   face_point(1,20) = 12
!   face_point(2,20) = 8
!   face_point(3,20) = 6
!
!   return
! end subroutine
! subroutine line_exp_point_dist_2d ( x1, y1, x2, y2, x, y, dist )
! !
! !*******************************************************************************
! !
! !! LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a line in 2D is:
! !
! !      (X1,Y1), (X2,Y2).
! !
! !  Modified:
! !
! !    27 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are two points on
! !    the line.
! !
! !    Input, real X, Y, the point whose distance from the line is
! !    to be measured.
! !
! !    Output, real DIST, the distance from the point to the line.
! !
!   real bot
!   real dist
!   real dot
!   real t
!   real x
!   real xn
!   real x1
!   real x2
!   real y
!   real yn
!   real y1
!   real y2
! !
!   bot = enormsq0_2d ( x1, y1, x2, y2 )
!
!   if ( bot == 0.0E+00 ) then
!
!     xn = x1
!     yn = y1
! !
! !  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
! !
! !  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
! !  of the projection of (P-P1) onto (P2-P1).
! !
!   else
!
!     dot = ( x - x1 ) * ( x2 - x1 ) + ( y - y1 ) * ( y2 - y1 )
!
!     t = dot / bot
!
!     xn = x1 + t * ( x2 - x1 )
!     yn = y1 + t * ( y2 - y1 )
!
!   end if
!
!   dist = enorm0_2d ( xn, yn, x, y )
!
!   return
! end subroutine
! subroutine line_exp_point_dist_3d ( x1, y1, z1, x2, y2, z2, x, y, z, dist )
! !
! !*******************************************************************************
! !
! !! LINE_EXP_POINT_DIST_3D: distance ( explicit line, point ) in 3D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a line in 2D is:
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2).
! !
! !  Modified:
! !
! !    02 November 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, (X1,Y1,Z1) and (X2,Y2,Z2) are
! !    two points on the line.  If the points are identical, then
! !    the line will be treated as a single point.
! !
! !    Input, real X, Y, Z, the point whose distance from the line is
! !    to be measured.
! !
! !    Output, real DIST, the distance from the point to the line.
! !
!   real bot
!   real dist
!   real t
!   real x
!   real xn
!   real x1
!   real x2
!   real y
!   real yn
!   real y1
!   real y2
!   real z
!   real zn
!   real z1
!   real z2
! !
!   bot = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )
!
!   if ( bot == 0.0E+00 ) then
!
!     xn = x1
!     yn = y1
!     zn = z1
! !
! !  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
! !
! !  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
! !  of the projection of (P-P1) onto (P2-P1).
! !
!   else
!
!     t = ( ( x - x1 ) * ( x2 - x1 ) + &
!           ( y - y1 ) * ( y2 - y1 ) + &
!           ( z - z1 ) * ( z2 - z1 ) ) / bot
!
!     xn = x1 + t * ( x2 - x1 )
!     yn = y1 + t * ( y2 - y1 )
!     zn = z1 + t * ( z2 - z1 )
!
!   end if
! !
! !  Now compute the distance between the projection point and P.
! !
!   dist = enorm0_3d ( x, y, z, xn, yn, zn )
!
!   return
! end subroutine
! subroutine line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dist_signed )
! !
! !*******************************************************************************
! !
! !! LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( explicit line, point ) in 2D.
! !
! !
! !  Discussion:
! !
! !    The signed distance has two interesting properties:
! !
! !    *  The absolute value of the signed distance is the
! !        usual (Euclidean) distance.
! !
! !    *  Points with signed distance 0 lie on the line,
! !       points with a negative signed distance lie on one side
! !         of the line,
! !       points with a positive signed distance lie on the
! !         other side of the line.
! !
! !    Assuming that C is nonnegative, then if a point is a positive
! !    distance away from the line, it is on the same side of the
! !    line as the point (0,0), and if it is a negative distance
! !    from the line, it is on the opposite side from (0,0).
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, define the two points
! !    (X1,Y1) and (X2,Y2) that determine the line.
! !
! !    Input, real X, Y, the point (X,Y) whose signed distance is desired.
! !
! !    Output, real DIST_SIGNED, the signed distance from the point to the line.
! !
!   real a
!   real b
!   real c
!   real dist_signed
!   real x
!   real x1
!   real x2
!   real y
!   real y1
!   real y2
! !
! !  Convert the line to implicit form.
! !
!   call line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )
! !
! !  Compute the signed distance from the point to the line.
! !
!   dist_signed = ( a * x + b * y + c ) / sqrt ( a * a + b * b )
!
!   return
! end subroutine
! subroutine line_exp_point_near_2d ( x1, y1, x2, y2, x, y, xn, yn, dist, t )
! !
! !*******************************************************************************
! !
! !! LINE_EXP_POINT_NEAR_2D computes the point on an explicit line nearest a point in 2D.
! !
! !
! !  Modified:
! !
! !    04 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are
! !    two points on the line.  (X1,Y1) must be different from (X2,Y2).
! !
! !    Input, real X, Y, the point whose nearest neighbor on the line is
! !    to be determined.
! !
! !    Output, real XN, YN, the nearest point on the line to (X,Y).
! !
! !    Output, real DIST, the distance from the point to the line.
! !
! !    Output, real T, the relative position of the point
! !    (XN,YN) to the points (X1,Y1) and (X2,Y2).
! !
! !    (XN,YN) = (1-T)*(X1,Y1) + T*(X2,Y2).
! !
! !    Less than 0, (XN,YN) is furthest away from (X2,Y2).
! !    Between 0 and 1, (XN,YN) is between (X1,Y1) and (X2,Y2).
! !    Greater than 1, (XN,YN) is furthest away from (X1,Y1).
! !
!   real bot
!   real dist
!   real t
!   real x
!   real xn
!   real x1
!   real x2
!   real y
!   real yn
!   real y1
!   real y2
! !
!   bot = enormsq0_2d ( x1, y1, x2, y2 )
!
!   if ( bot == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINE_POINT_NEAR_2D - Fatal error!'
!     write ( *, * ) '  The points (X1,Y1) and (X2,Y2) are identical.'
!     stop
!   end if
! !
! !  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
! !
! !  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
! !  of the projection of (P-P1) onto (P2-P1).
! !
!   t = ( ( x1 - x ) * ( x1 - x2 ) + ( y1 - y ) * ( y1 - y2 ) ) / bot
!
!   xn = x1 + t * ( x2 - x1 )
!   yn = y1 + t * ( y2 - y1 )
!
!   dist = enorm0_2d ( xn, yn, x, y )
!
!   return
! end subroutine
! subroutine line_exp_point_near_3d ( x1, y1, z1, x2, y2, z2, x, y, z, &
!   xn, yn, zn, dist, t )
! !
! !*******************************************************************************
! !
! !! LINE_EXP_POINT_NEAR_3D: nearest point on explicit line to point in 3D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2.
! !    (X1,Y1,Z1) and (X2,Y2,Z2) are two points on the line.  They
! !    must be distinct.
! !
! !    Input, real X, Y, Z, the point whose nearest neighbor on the line
! !    is to be determined.
! !
! !    Output, real XN, YN, ZN, the point which is the nearest point on
! !    the line to (X,Y,Z).
! !
! !    Output, real DIST, the distance from the point to the nearest point
! !    on the line.
! !
! !    Output, real T, the relative position of the point
! !    (XN,YN,ZN) to the points (X1,Y1,Z1) and (X2,Y2,Z2).
! !
! !      (XN,YN,ZN) = (1-T)*(X1,Y1,Z1) + T*(X2,Y2,Z2).
! !
! !    Less than 0, (XN,YN,ZN) is furthest away from (X2,Y2,Z2).
! !    Between 0 and 1, (XN,YN,ZN) is between (X1,Y1,Z1) and (X2,Y2,Z2).
! !    Greater than 1, (XN,YN,ZN) is furthest away from (X1,Y1,Z1).
! !
!   real bot
!   real dist
!   real t
!   real x
!   real xn
!   real x1
!   real x2
!   real y
!   real yn
!   real y1
!   real y2
!   real z
!   real zn
!   real z1
!   real z2
! !
!   bot = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )
!
!   if ( bot == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINE_EXP_POINT_NEAR_3D - Fatal error!'
!     write ( *, * ) '  The points (X1,Y1,Z1) and (X2,Y2,Z2) '
!     write ( *, * ) '  are identical.'
!     stop
!   end if
! !
! !  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
! !
! !  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
! !  of the projection of (P-P1) onto (P2-P1).
! !
!   t = (   ( x - x1 ) * ( x2 - x1 ) &
!         + ( y - y1 ) * ( y2 - y1 ) &
!         + ( z - z1 ) * ( z2 - z1 ) ) / bot
! !
! !  Now compute the location of the projection point.
! !
!   xn = x1 + t * ( x2 - x1 )
!   yn = y1 + t * ( y2 - y1 )
!   zn = z1 + t * ( z2 - z1 )
! !
! !  Now compute the distance between the projection point and P.
! !
!   dist = enorm0_3d ( xn, yn, zn, x, y, z )
!
!   return
! end subroutine
!
! subroutine line_exp2par_2d ( x1, y1, x2, y2, f, g, x0, y0 )
! !
! !*******************************************************************************
! !
! !! LINE_EXP2PAR_2D converts a line from explicit to parametric form in 2D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a line in 2D is:
! !
! !      (X1,Y1), (X2,Y2).
! !
! !    The parametric form of a line in 2D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !
! !  Modified:
! !
! !    24 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, two points on the line.
! !
! !    Output, real F, G, X0, Y0, the parametric parameters of the line.
! !
!   real f
!   real g
!   real norm
!   real x0
!   real x1
!   real x2
!   real y0
!   real y1
!   real y2
! !
!   x0 = x1
!   y0 = y1
!
!   f = x2 - x1
!   g = y2 - y1
!
!   norm = sqrt ( f * f + g * g )
!
!   if ( norm /= 0.0E+00 ) then
!     f = f / norm
!     g = g / norm
!   end if
!
!   return
! end subroutine
! subroutine line_exp2par_3d ( x1, y1, z1, x2, y2, z2, f, g, h, x0, y0, z0 )
! !
! !*******************************************************************************
! !
! !! LINE_EXP2PAR_3D converts a line from explicit to parametric form in 3D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a line in 3D is:
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2).
! !
! !    The parametric form of a line in 3D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !      Z = Z0 + H * T
! !
! !  Modified:
! !
! !    30 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z3, two points on the line.
! !
! !    Output, real F, G, H, X0, Y0, Z0, the parametric parameters of the line.
! !
!   real f
!   real g
!   real h
!   real norm
!   real x0
!   real x1
!   real x2
!   real y0
!   real y1
!   real y2
!   real z0
!   real z1
!   real z2
! !
!   x0 = x1
!   y0 = y1
!   z0 = z1
!
!   f = x2 - x1
!   g = y2 - y1
!   h = z2 - z1
!
!   norm = sqrt ( f * f + g * g + h * h )
!
!   if ( norm /= 0.0E+00 ) then
!     f = f / norm
!     g = g / norm
!     h = h / norm
!   end if
!
!   return
! end subroutine
!
! subroutine line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )
! !
! !*******************************************************************************
! !
! !! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a line in 2D is:
! !
! !      (X1,Y1), (X2,Y2).
! !
! !    The implicit form of a line in 2D is:
! !
! !      A * X + B * Y + C = 0
! !
! !  Modified:
! !
! !    24 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are
! !    two points on the line. (X1,Y1) must be different
! !    from (X2,Y2).
! !
! !    Output, real A, B, C, three coefficients which describe
! !    the line that passes through (X1,Y1) and (X2,Y2).
! !
!   real a
!   real b
!   real c
!   real x1
!   real x2
!   real y1
!   real y2
! !
! !  Take care of degenerate cases.
! !
!   if ( x1 == x2 .and. y1 == y2 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINE_EXP2IMP_2D - Fatal error!'
!     write ( *, * ) '  (X1,Y1) = (X2,Y2)'
!     write ( *, * ) '  (X1,Y1) = ', x1, y1
!     write ( *, * ) '  (X2,Y2) = ', x2, y2
!     stop
!   end if
!
!   a = y2 - y1
!   b = x1 - x2
!   c = x2 * y1 - x1 * y2
!
!   return
! end subroutine
!
! subroutine line_imp2par_2d ( a, b, c, f, g, x0, y0 )
! !
! !*******************************************************************************
! !
! !! LINE_IMP2PAR_2D converts an implicit line to parametric form in 2D.
! !
! !
! !  Formula:
! !
! !    The implicit form of line in 2D is:
! !
! !      A * X + B * Y + C = 0
! !
! !    The parametric form of a line in 2D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    24 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, the implicit parameters of the line.
! !
! !    Output, real F, G, X0, Y0, the parametric parameters of the line.
! !
!   real a
!   real b
!   real c
!   real f
!   real g
!   real x0
!   real y0
! !
!   if ( a * a + b * b == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINE_IMP2PAR_2D - Fatal error!'
!     write ( *, * ) '  A * A + B * B = 0.'
!     stop
!   end if
!
!   x0 = - a * c / ( a * a + b * b )
!   y0 = - b * c / ( a * a + b * b )
!   f =   b / sqrt ( a * a + b * b )
!   g = - a / sqrt ( a * a + b * b )
!
!   return
! end subroutine
! subroutine line_imp_point_dist_2d ( a, b, c, x, y, dist )
! !
! !*******************************************************************************
! !
! !! LINE_IMP_POINT_DIST_2D: distance ( implicit line, point ) in 2D.
! !
! !
! !  Formula:
! !
! !    The implicit form of a line in 2D is:
! !
! !      A * X + B * Y + C = 0
! !
! !  Modified:
! !
! !    27 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, the implicit line parameters.
! !
! !    Input, real X, Y, the point whose distance from the line is
! !    to be measured.
! !
! !    Output, real DIST, the distance from the point to the line.
! !
!   real a
!   real b
!   real c
!   real dist
!   real x
!   real y
! !
!   if ( a * a + b * b == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINE_IMP_POINT_DIST_2D - Fatal error!'
!     write ( *, * ) '  A * A + B * B = 0.'
!     stop
!   end if
!
!   dist = abs ( a * x + b * y + c ) / sqrt ( a * a + b * b )
!
!   return
! end subroutine
! subroutine line_imp_point_dist_3d ( a, b, c, d, x, y, z, dist )
! !
! !*******************************************************************************
! !
! !! LINE_IMP_POINT_DIST_3D: distance ( implicit line, point ) in 3D.
! !
! !
! !  Formula:
! !
! !    The implicit form of a line in 2D is:
! !
! !      A * X + B * Y + C * Z + D = 0
! !
! !  Modified:
! !
! !    27 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, D, the implicit line parameters.
! !
! !    Input, real X, Y, Z, the point whose distance from the line is
! !    to be measured.
! !
! !    Output, real DIST, the distance from the point to the line.
! !
!   real a
!   real b
!   real c
!   real d
!   real dist
!   real x
!   real y
!   real z
! !
!   if ( a * a + b * b + c * c == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINE_IMP_POINT_DIST_3D - Fatal error!'
!     write ( *, * ) '  A * A + B * B + C * C = 0.'
!     stop
!   end if
!
!   dist = abs ( a * x + b * y + c * z + d ) / sqrt ( a * a + b * b + c * c )
!
!   return
! end subroutine
! subroutine line_imp_point_dist_signed_2d ( a, b, c, x, y, dist_signed )
! !
! !*******************************************************************************
! !
! !! LINE_IMP_POINT_DIST_SIGNED_2D: signed distance ( implicit line, point ) in 2D.
! !
! !
! !  Modified:
! !
! !    04 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, the equation of the line is A*X + B*Y + C = 0.
! !
! !    Input, real X, Y, the coordinates of the point.
! !
! !    Output, real DIST_SIGNED, the signed distance from the point to
! !    the line.
! !
!   real a
!   real b
!   real c
!   real dist_signed
!   real x
!   real y
! !
!   if ( a * a + b * b == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINE_IMP_POINT_DIST_SIGNED_2D - Fatal error!'
!     write ( *, * ) '  A * A + B * B = 0.'
!     stop
!   end if
!
!   dist_signed = - sign ( 1.0E+00, c ) * ( a * x + b * y + c ) / &
!     sqrt ( a * a + b * b )
!
!   return
! end subroutine
! subroutine line_par2imp_2d ( f, g, x0, y0, a, b, c )
! !
! !*******************************************************************************
! !
! !! LINE_PAR2IMP_2D converts a parametric line to implicit form in 2D.
! !
! !
! !  Formula:
! !
! !    The parametric form of a line in 2D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !
! !    The implicit form of a line in 2D is:
! !
! !      A * X + B * Y + C = 0
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    24 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real F, G, X0, Y0, the parametric parameters of the line.
! !
! !    Output, real A, B, C, the implicit parameters of the line.
! !
!   real a
!   real b
!   real c
!   real f
!   real g
!   real x0
!   real y0
! !
!   a = - g
!   b = f
!   c = g * x0 - f * y0
!
!   return
! end subroutine
! subroutine line_par_point_dist_2d ( f, g, x0, y0, x, y, dist )
! !
! !*******************************************************************************
! !
! !! LINE_PAR_POINT_DIST_2D: distance ( parametric line, point ) in 2D.
! !
! !
! !  Formula:
! !
! !    The parametric form of a line in 2D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    27 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real F, G, X0, Y0, the parametric line parameters.
! !
! !    Input, real X, Y, the point whose distance from the line is
! !    to be measured.
! !
! !    Output, real DIST, the distance from the point to the line.
! !
!   real dist
!   real dx
!   real dy
!   real f
!   real g
!   real x
!   real x0
!   real y
!   real y0
! !
!   dx =   g * g * ( x - x0 ) - f * g * ( y - y0 )
!   dy = - f * g * ( x - x0 ) + f * f * ( y - y0 )
!
!   dist = sqrt ( dx * dx + dy * dy ) / ( f * f + g * g )
!
!   return
! end subroutine
! subroutine line_par_point_dist_3d ( f, g, h, x0, y0, z0, x, y, z, dist )
! !
! !*******************************************************************************
! !
! !! LINE_PAR_POINT_DIST_3D: distance ( parametric line, point ) in 3D.
! !
! !
! !  Formula:
! !
! !    The parametric form of a line in 3D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !      Z = Z0 + H * T
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    27 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real F, G, H, X0, Y0, Z0, the parametric line parameters.
! !
! !    Input, real X, Y, Z, the point whose distance from the line is
! !    to be measured.
! !
! !    Output, real DIST, the distance from the point to the line.
! !
!   real dist
!   real dx
!   real dy
!   real dz
!   real f
!   real g
!   real h
!   real x
!   real x0
!   real y
!   real y0
!   real z
!   real z0
! !
!   dx =   g * ( f * ( y - y0 ) - g * ( x - x0 ) ) &
!        + h * ( f * ( z - z0 ) - h * ( x - x0 ) )
!
!   dy =   h * ( g * ( z - z0 ) - h * ( y - y0 ) )  &
!        - f * ( f * ( y - y0 ) - g * ( x - x0 ) )
!
!   dz = - f * ( f * ( z - z0 ) - h * ( x - x0 ) ) &
!        - g * ( g * ( z - z0 ) - h * ( y - y0 ) )
!
!   dist = sqrt ( dx * dx + dy * dy + dz * dz ) &
!     / ( f * f + g * g + h * h )
!
!   return
! end subroutine
! subroutine line_seg_contains_point_1d ( x1, x2, x3, u )
! !
! !*******************************************************************************
! !
! !! LINE_SEG_CONTAINS_POINT_1D reports if a line segment contains a point in 1D.
! !
! !
! !  Modified:
! !
! !    28 September 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, X2, two points defining a line segment.
! !    The line segment has origin at X1, and unit at X2.
! !
! !    Input, real X3, a point to be tested.
! !
! !    Output, real U, the coordinate of X3 in units of (X2-X1).
! !    The point X3 is contained in the line segment if 0 <= U <= 1.
! !
!   real u
!   real unit
!   real x1
!   real x2
!   real x3
! !
!   unit = x2 - x1
!
!   if ( unit == 0.0E+00 ) then
!
!     if ( x3 == x1 ) then
!       u = 0.5E+00
!     else if ( x3 < x1 ) then
!       u = - huge ( u )
!     else if ( x3 > x1 ) then
!       u = huge ( u )
!     end if
!
!   else
!
!     u = ( x3 - x1 ) / unit
!
!   end if
!
!   return
! end subroutine
! subroutine line_seg_contains_point_2d ( x1, y1, x2, y2, x3, y3, u, v )
! !
! !*******************************************************************************
! !
! !! LINE_SEG_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
! !
! !
! !  Discussion:
! !
! !    In exact arithmetic, point P3 = (X3,Y3) is on the line segment between
! !    P1=(X1,Y1) and P2=(X2,Y2) if and only if 0 <= U <= 1 and V = 0.
! !
! !  Modified:
! !
! !    28 September 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, the endpoints P1 and P2 of a line segment.
! !
! !    Input, real X3, Y3, a point P3 to be tested.
! !
! !    Output, real U, the coordinate of (X3,Y3) along the axis from
! !    with origin at P1 and unit at P2.
! !
! !    Output, real V, the magnitude of the off-axis portion of the
! !    vector P3-P1, measured in units of (P2-P1).
! !
!   real u
!   real unit
!   real v
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
! !
!   unit = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )
!
!   if ( unit == 0.0E+00 ) then
!
!     if ( x3 == x1 .and. y3 == y1 ) then
!       u = 0.5E+00
!       v = 0.0E+00
!     else
!       u = 0.5E+00
!       v = huge ( v )
!     end if
!
!   else
!
!     u = ( ( x3 - x1 ) * ( x2 - x1 ) + ( y3 - y1 ) * ( y2 - y1 ) ) / unit**2
!
!     v = sqrt ( ( ( u - 1.0E+00 ) * x1 - u * x2 + x3 )**2 &
!              + ( ( u - 1.0E+00 ) * y1 - u * y2 + y3 )**2 ) / unit
!
!   end if
!
!   return
! end subroutine
! subroutine line_seg_point_dist_2d ( x1, y1, x2, y2, x, y, dist )
! !*******************************************************************************
! ! LINE_SEG_POINT_DIST_2D: distance ( line segment, point ) in 2D.
! !  Modified:
! !    02 November 1998
! !  Author:
! !    John Burkardt
! !  Parameters:
! !    Input, real X1, Y1, X2, Y2, the endpoints of the line segment.
! !    Input, real X, Y, the point whose nearest neighbor on the line
! !    segment is to be determined.
! !    Output, real DIST, the distance from the point to the line segment.
! !
!   real bot
!   real dist
!   real t
!   real x
!   real xn
!   real x1
!   real x2
!   real y
!   real yn
!   real y1
!   real y2
! !
! !  If the line segment is actually a point, then the answer is easy.
! !
!   if ( x1 == x2 .and. y1 == y2 ) then
!
!     xn = x1
!     yn = y1
!
!   else
!
!     bot = enormsq0_2d ( x1, y1, x2, y2 )
!
!     t = ( ( x1 - x ) * ( x1 - x2 ) + ( y1 - y ) * ( y1 - y2 ) ) / bot
!
!     t = max ( t, 0.0E+00 )
!     t = min ( t, 1.0E+00 )
!
!     xn = x1 + t * ( x2 - x1 )
!     yn = y1 + t * ( y2 - y1 )
!
!   end if
!
!   dist = enorm0_2d ( xn, yn, x, y )
!
!   return
! end subroutine
!
! subroutine line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, x, y, z, dist )
! !*******************************************************************************
! ! LINE_SEG_POINT_DIST_3D: distance ( line segment, point ) in 3D.
! !  Modified:
! !    02 November 1998
! !  Author:
! !    John Burkardt
! !  Parameters:
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, the endpoints of the line
! !    segment.
! !    Input, real X, Y, Z, the point whose nearest neighbor on the line
! !    segment is to be determined.
! !    Output, real DIST, the distance from the point to the line segment.
!   real bot
!   real dist
!   real t
!   real x
!   real xn
!   real x1
!   real x2
!   real y
!   real yn
!   real y1
!   real y2
!   real z
!   real zn
!   real z1
!   real z2
! !
! !  If the line segment is actually a point, then the answer is easy.
! !
!   if ( x1 == x2 .and. y1 == y2 .and. z1 == z2 ) then
!
!     xn = x1
!     yn = y1
!     zn = z1
!
!   else
!
!     bot = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )
!
!     t = (   ( x1 - x ) * ( x1 - x2 ) &
!           + ( y1 - y ) * ( y1 - y2 ) &
!           + ( z1 - z ) * ( z1 - z2 ) ) / bot
!
!     t = max ( t, 0.0E+00 )
!     t = min ( t, 1.0E+00 )
!
!     xn = x1 + t * ( x2 - x1 )
!     yn = y1 + t * ( y2 - y1 )
!     zn = z1 + t * ( z2 - z1 )
!
!   end if
!
!   dist = enorm0_3d ( x, y, z, xn, yn, zn )
!
!   return
! end subroutine
! subroutine line_seg_point_near_2d ( x1, y1, x2, y2, x, y, xn, yn, dist, t )
! !
! !*******************************************************************************
! !
! !! LINE_SEG_POINT_NEAR_2D finds the point on a line segment nearest a point in 2D.
! !
! !
! !  Modified:
! !
! !    07 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, the two endpoints of the line segment.
! !
! !    (X1,Y1) should generally be different from (X2,Y2), but
! !    if they are equal, the program will still compute a
! !    meaningful result.
! !
! !    Input, real X, Y, the point whose nearest neighbor
! !    on the line segment is to be determined.
! !
! !    Output, real XN, YN, the point on the line segment which is
! !    nearest the point (X,Y).
! !
! !    Output, real DIST, the distance from the point to the nearest point
! !    on the line segment.
! !
! !    Output, real T, the relative position of the point (XN,YN) to the
! !    points (X1,Y1) and (X2,Y2).
! !
! !      (XN,YN) = (1-T) * (X1,Y1) + T * (X2,Y2).
! !
! !    T will always be between 0 and 1.
! !
!   real bot
!   real dist
!   real t
!   real x
!   real xn
!   real x1
!   real x2
!   real y
!   real yn
!   real y1
!   real y2
! !
!   if ( x1 == x2 .and. y1 == y2 ) then
!
!     t = 0.0E+00
!     xn = x1
!     yn = y1
!
!   else
!
!     bot = enormsq0_2d ( x1, y1, x2, y2 )
!
!     t = (   ( x1 - x ) * ( x1 - x2 ) &
!           + ( y1 - y ) * ( y1 - y2 ) ) / bot
!
!     if ( t < 0.0E+00 ) then
!       t = 0.0E+00
!       xn = x1
!       yn = y1
!     else if ( t > 1.0E+00 ) then
!       t = 1.0E+00
!       xn = x2
!       yn = y2
!     else
!       xn = x1 + t * ( x2 - x1 )
!       yn = y1 + t * ( y2 - y1 )
!     end if
!
!   end if
!
!   dist = enorm0_2d ( x, y, xn, yn )
!
!   return
! end subroutine
! subroutine line_seg_point_near_3d ( x1, y1, z1, x2, y2, z2, x, y, z, &
!   xn, yn, zn, dist, t )
! !
! !*******************************************************************************
! !
! !! LINE_SEG_POINT_NEAR_3D finds the point on a line segment nearest a point in 3D.
! !
! !
! !  Modified:
! !
! !    07 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, the two endpoints of the line segment.
! !
! !    (X1,Y1,Z1) should generally be different from (X2,Y2,Z2), but
! !    if they are equal, the program will still compute a meaningful
! !    result.
! !
! !    Input, real X, Y, Z, the point whose nearest neighbor
! !    on the line segment is to be determined.
! !
! !    Output, real XN, YN, ZN, the point on the line segment which is
! !    nearest the point (X,Y,Z).
! !
! !    Output, real DIST, the distance from the point to the nearest point
! !    on the line segment.
! !
! !    Output, real T, the relative position of the nearest point
! !    (XN,YN,ZN) to the defining points (X1,Y1,Z1) and (X2,Y2,Z2).
! !
! !      (XN,YN,ZN) = (1-T)*(X1,Y1,Z1) + T*(X2,Y2,Z2).
! !
! !    T will always be between 0 and 1.
! !
!   real bot
!   real dist
!   real t
!   real x
!   real xn
!   real x1
!   real x2
!   real y
!   real yn
!   real y1
!   real y2
!   real z
!   real zn
!   real z1
!   real z2
! !
!   if ( x1 == x2 .and. y1 == y2 .and. z1 == z2 ) then
!
!     t = 0.0E+00
!     xn = x1
!     yn = y1
!     zn = z1
!
!   else
!
!     bot = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )
!
!     t = (   ( x1 - x ) * ( x1 - x2 ) &
!           + ( y1 - y ) * ( y1 - y2 ) &
!           + ( z1 - z ) * ( z1 - z2 ) ) / bot
!
!     if ( t < 0.0E+00 ) then
!       t = 0.0E+00
!       xn = x1
!       yn = y1
!       zn = z1
!     else if ( t > 1.0E+00 ) then
!       t = 1.0E+00
!       xn = x2
!       yn = y2
!       zn = z2
!     else
!       xn = x1 + t * ( x2 - x1 )
!       yn = y1 + t * ( y2 - y1 )
!       zn = z1 + t * ( z2 - z1 )
!     end if
!
!   end if
!
!   dist = enorm0_3d ( xn, yn, zn, x, y, z )
!
!   return
! end subroutine
! subroutine lines_exp_angle_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, angle )
! !
! !*******************************************************************************
! !
! !! LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a line in 3D is:
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2).
! !
! !  Modified:
! !
! !    24 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, two distince points on the first line.
! !
! !    Input, real X3, Y3, Z3, X4, Y4, Z4, two distinct points on the second line.
! !
! !    Output, real ANGLE, the angle in radians between the two lines.
! !    ANGLE is computed using the ACOS function, and so lies between 0 and PI.
! !    But if one of the lines is degenerate, ANGLE is returned as -1.0.
! !
!   real angle
!   real ctheta
!   real pdotq
!   real pnorm
!   real qnorm
!   real x1
!   real x2
!   real x3
!   real x4
!   real y1
!   real y2
!   real y3
!   real y4
!   real z1
!   real z2
!   real z3
!   real z4
! !
!   pnorm = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
!   qnorm = enorm0_3d ( x3, y3, z3, x4, y4, z4 )
!
!   pdotq =    ( x2 - x1 ) * ( x4 - x3 ) &
!            + ( y2 - y1 ) * ( y4 - y3 ) &
!            + ( z2 - z1 ) * ( z4 - z3 )
!
!   if ( pnorm == 0.0E+00 .or. qnorm == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINES_EXP_ANGLE_3D - Fatal error!'
!     write ( *, * ) '  One of the lines is degenerate!'
!     angle = - 1.0E+00
!
!   else
!
!     ctheta = pdotq / ( pnorm * qnorm )
!     angle = arc_cosine ( ctheta )
!
!   end if
!
!   return
! end subroutine
! subroutine lines_exp_angle_nd ( p1, p2, q1, q2, n, angle )
! !
! !*******************************************************************************
! !
! !! LINES_EXP_ANGLE_ND returns the angle between two explicit lines in ND.
! !
! !
! !  Modified:
! !
! !    27 October 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real P1(N), P2(N), two points on the first line.
! !
! !    Input, real Q1(N), Q2(N), two points on the second line.
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Output, real ANGLE, the angle in radians between the two lines.
! !    ANGLE is computed using the ACOS function, and so lies between 0 and PI.
! !    But if one of the lines is degenerate, ANGLE is returned as -1.0.
! !
!   integer n
! !
!   real angle
!   real ctheta
!   integer i
!   real p1(n)
!   real p2(n)
!   real pdotq
!   real pnorm
!   real q1(n)
!   real q2(n)
!   real qnorm
! !
!   pnorm = enorm0_nd ( n, p1, p2 )
!   qnorm = enorm0_nd ( n, q1, q2 )
!
!   pdotq = 0.0E+00
!   do i = 1, n
!     pdotq = pdotq + ( p2(i) - p1(i) ) * ( q2(i) - q1(i) )
!   end do
!
!   if ( pnorm == 0.0E+00 .or. qnorm == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINES_EXP_ANGLE_ND - Fatal error!'
!     write ( *, * ) '  One of the lines is degenerate!'
!     angle = - 1.0E+00
!   else
!     ctheta = pdotq / ( pnorm * qnorm )
!     angle = arc_cosine ( ctheta )
!
!   end if
!
!   return
! end subroutine
! subroutine lines_exp_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
!   dist )
! !
! !*******************************************************************************
! !
! !! LINES_EXP_DIST_3D computes the distance between two explicit lines in 3D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a line in 3D is:
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2).
! !
! !  Modified:
! !
! !    24 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, (X1,Y1,Z1) and (X2,Y2,Z2) are
! !    two points on the first line.  They must be distinct.
! !
! !    Input, real X3, Y3, Z3, X4, Y4, Z4, (X3,Y3,Z3) and (X4,Y4,Z4) are
! !    two points on the second line.  They must be distinct.
! !
! !    Output, real DIST, the distance between the lines.
! !
!   real a11
!   real a12
!   real a13
!   real a21
!   real a22
!   real a23
!   real a31
!   real a32
!   real a33
!   real bot
!   real cr1
!   real cr2
!   real cr3
!   real dist
!   real top
!   real x1
!   real x2
!   real x3
!   real x4
!   real y1
!   real y2
!   real y3
!   real y4
!   real z1
!   real z2
!   real z3
!   real z4
! !
! !  The distance is found by computing the volume of a parallelipiped,
! !  and dividing by the area of its base.
! !
! !  But if the lines are parallel, we compute the distance by
! !  finding the distance between the first line and any point
! !  on the second line.
! !
!   a11 = x3 - x1
!   a12 = y3 - y1
!   a13 = z3 - z1
!
!   a21 = x2 - x1
!   a22 = y2 - y1
!   a23 = z2 - z1
!
!   a31 = x4 - x3
!   a32 = y4 - y3
!   a33 = z4 - z3
!
!   call cross_3d ( a21, a22, a23, a31, a32, a33, cr1, cr2, cr3 )
!
!   bot = sqrt ( cr1 * cr1 + cr2 * cr2 + cr3 * cr3 )
!
!   if ( bot == 0.0E+00 ) then
!
!     call line_exp_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, dist )
!
!   else
!
!     top = abs (   a11 * ( a22 * a33 - a23 * a32 ) &
!                 - a12 * ( a21 * a33 - a23 * a31 ) &
!                 + a13 * ( a21 * a32 - a22 * a31 ) )
!
!     dist = top / bot
!
!   end if
!
!   return
! end subroutine
! subroutine lines_exp_int_2d ( ival, x, y, x1, y1, x2, y2, x3, y3, x4, y4 )
! !
! !*******************************************************************************
! !
! !! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
! !
! !
! !  Formula:
! !
! !    The explicit form of a line in 2D is:
! !
! !      (X1,Y1), (X2,Y2).
! !
! !  Modified:
! !
! !    28 September 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, integer IVAL, reports on the intersection.
! !
! !     0, no intersection, the lines may be parallel or degenerate.
! !     1, one intersection point, returned in X, Y.
! !     2, infinitely many intersections, the lines are identical.
! !
! !    Output, real X, Y, if IVAl = 1, then X, Y contains
! !    the intersection point.  Otherwise, X = 0, Y = 0.
! !
! !    Input, real X1, Y1, X2, Y2, define the first line.
! !
! !    Input, real X3, Y3, X4, Y4, define the second line.
! !
!   real a1
!   real a2
!   real b1
!   real b2
!   real c1
!   real c2
!   integer ival
!   logical point_1
!   logical point_2
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
! !
!   ival = 0
!   x = 0.0E+00
!   y = 0.0E+00
! !
! !  Check whether either line is a point.
! !
!   if ( x1 == x2 .and. y1 == y2 ) then
!     point_1 = .true.
!   else
!     point_1 = .false.
!   end if
!
!   if ( x3 == x4 .and. y3 == y4 ) then
!     point_2 = .true.
!   else
!     point_2 = .false.
!   end if
! !
! !  Convert the lines to ABC format.
! !
!   if ( .not. point_1 ) then
!     call line_exp2imp_2d ( x1, y1, x2, y2, a1, b1, c1 )
!   end if
!
!   if ( .not. point_2 ) then
!     call line_exp2imp_2d ( x3, y3, x4, y4, a2, b2, c2 )
!   end if
! !
! !  Search for intersection of the lines.
! !
!   if ( point_1 .and. point_2 ) then
!     if ( x1 == x3 .and. y1 == y3 ) then
!       ival = 1
!       x = x1
!       y = y1
!     end if
!   else if ( point_1 ) then
!     if ( a2 * x1 + b2 * y1 == c2 ) then
!       ival = 1
!       x = x1
!       y = y1
!     end if
!   else if ( point_2 ) then
!     if ( a1 * x3 + b1 * y3 == c1 ) then
!       ival = 1
!       x = x3
!       y = y3
!     end if
!   else
!     call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )
!   end if
!
!   return
! end subroutine
!
! subroutine lines_imp_angle_2d ( a1, b1, c1, a2, b2, c2, theta )
! !
! !*******************************************************************************
! !
! !! LINES_IMP_ANGLE_2D finds the angle between two implicit lines in 2D.
! !
! !
! !  Formula:
! !
! !    The implicit form of a line in 2D is:
! !
! !      A * X + B * Y + C = 0
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    24 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A1, B1, C1, the implicit parameters of the first line.
! !
! !    Input, real A2, B2, C2, the implicit parameters of the second line.
! !
! !    Output, real THETA, the angle between the two lines.
! !
!   real a1
!   real a2
!   real b1
!   real b2
!   real c1
!   real c2
!   real pdotq
!   real pnorm
!   real qnorm
!   real theta
! !
!   pdotq = a1 * a2 + b1 * b2
!   pnorm = sqrt ( a1 * a1 + b1 * b1 )
!   qnorm = sqrt ( a2 * a2 + b2 * b2 )
!
!   theta = arc_cosine ( pdotq / ( pnorm * qnorm ) )
!
!   return
! end subroutine
! subroutine lines_imp_dist_2d ( a1, b1, c1, a2, b2, c2, dist )
! !
! !*******************************************************************************
! !
! !! LINES_IMP_DIST_2D determines the distance between two implicit lines in 2D.
! !
! !
! !  Discussion:
! !
! !    If the lines intersect, then their distance is zero.
! !    If the two lines are parallel, then they have a nonzero distance.
! !
! !  Formula:
! !
! !    The implicit form of a line in 2D is:
! !
! !      A * X + B * Y + C = 0
! !
! !  Modified:
! !
! !    12 January 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A1, B1, C1, define the first line.
! !    At least one of A1 and B1 must be nonzero.
! !
! !    Input, real A2, B2, C2, define the second line.
! !    At least one of A2 and B2 must be nonzero.
! !
! !    Output, real DIST, the distance between the two lines.
! !
!   real a1
!   real a2
!   real b1
!   real b2
!   real c1
!   real c2
!   real dist
! !
! !  Refuse to handle degenerate lines.
! !
!   if ( a1 == 0.0E+00 .and. b1 == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINES_IMP_DIST_2D - Fatal error!'
!     write ( *, * ) '  Line 1 is degenerate.'
!     stop
!   else if ( a2 == 0.0E+00 .and. b2 == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'LINES_IMP_DIST_2D - Fatal error!'
!     write ( *, * ) '  Line 2 is degenerate.'
!     stop
!   end if
! !
! !  Determine if the lines intersect.
! !
!   if ( a1 * b2 /= a2 * b1 ) then
!     dist = 0.0E+00
!     return
!   end if
! !
! !  Determine the distance between the parallel lines.
! !
!   dist = abs ( c2 / sqrt ( a2**2 + b2**2 ) - c1 / sqrt ( a1**2 + b1**2 ) )
!
!   return
! end subroutine
! subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )
! !
! !*******************************************************************************
! !
! !! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
! !
! !
! !  Formula:
! !
! !    The implicit form of a line in 2D is:
! !
! !      A * X + B * Y + C = 0
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A1, B1, C1, define the first line.
! !    At least one of A1 and B1 must be nonzero.
! !
! !    Input, real A2, B2, C2, define the second line.
! !    At least one of A2 and B2 must be nonzero.
! !
! !    Output, integer IVAL, reports on the intersection.
! !
! !    -1, both A1 and B1 were zero.
! !    -2, both A2 and B2 were zero.
! !     0, no intersection, the lines are parallel.
! !     1, one intersection point, returned in X, Y.
! !     2, infinitely many intersections, the lines are identical.
! !
! !    Output, real X, Y, if IVAL = 1, then X, Y contains
! !    the intersection point.  Otherwise, X = 0, Y = 0.
! !
!   real a(2,2)
!   real a1
!   real a2
!   real b(2,2)
!   real b1
!   real b2
!   real c1
!   real c2
!   real det
!   integer ival
!   real x
!   real y
! !
!   x = 0.0E+00
!   y = 0.0E+00
! !
! !  Refuse to handle degenerate lines.
! !
!   if ( a1 == 0.0E+00 .and. b1 == 0.0E+00 ) then
!     ival = -1
!     return
!   else if ( a2 == 0.0E+00 .and. b2 == 0.0E+00 ) then
!     ival = -2
!     return
!   end if
! !
! !  Set up a linear system, and compute its inverse.
! !
!   a(1,1) = a1
!   a(1,2) = b1
!   a(2,1) = a2
!   a(2,2) = b2
!
!   call rmat2_inverse ( a, b, det )
! !
! !  If the inverse exists, then the lines intersect.
! !  Multiply the inverse times -C to get the intersection point.
! !
!   if ( det /= 0.0E+00 ) then
!
!     ival = 1
!     x = - b(1,1) * c1 - b(1,2) * c2
!     y = - b(2,1) * c1 - b(2,2) * c2
! !
! !  If the inverse does not exist, then the lines are parallel
! !  or coincident.  Check for parallelism by seeing if the
! !  C entries are in the same ratio as the A or B entries.
! !
!   else
!
!     ival = 0
!
!     if ( a1 == 0.0E+00 ) then
!       if ( b2 * c1 == c2 * b1 ) then
!         ival = 2
!       end if
!     else
!       if ( a2 * c1 == c2 * a1 ) then
!         ival = 2
!       end if
!     end if
!
!   end if
!
!   return
! end subroutine
!
! subroutine lines_par_angle_2d ( f1, g1, x01, y01, f2, g2, x02, y02, theta )
! !
! !*******************************************************************************
! !
! !! LINES_PAR_ANGLE_2D finds the angle between two parametric lines in 2D.
! !
! !
! !  Formula:
! !
! !    The parametric form of a line in 2D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    24 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real F1, G1, X01, Y01, the parametric parameters of the
! !    first line.
! !
! !    Input, real F2, G2, X02, Y02, the parametric parameters of the
! !    second line.
! !
! !    Output, real THETA, the angle between the two lines.
! !
!   real f1
!   real f2
!   real g1
!   real g2
!   real pdotq
!   real pnorm
!   real qnorm
!   real theta
!   real x01
!   real x02
!   real y01
!   real y02
! !
!   pdotq = f1 * f2 + g1 * g2
!   pnorm = sqrt ( f1 * f1 + g1 * g1 )
!   qnorm = sqrt ( f2 * f2 + g2 * g2 )
!
!   theta = arc_cosine ( pdotq / ( pnorm * qnorm ) )
!
!   return
! end subroutine
! subroutine lines_par_angle_3d ( f1, g1, h1, x01, y01, z01, f2, g2, h2, &
!   x02, y02, z02, theta )
! !
! !*******************************************************************************
! !
! !! LINES_PAR_ANGLE_3D finds the angle between two parametric lines in 3D.
! !
! !
! !  Formula:
! !
! !    The parametric form of a line in 3D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !      Z = Z0 + H * T
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real F1, G1, H1, X01, Y01, Z01, the parametric parameters
! !    of the first line.
! !
! !    Input, real F2, G2, H2, X02, Y02, Z02, the parametric parameters
! !    of the second line.
! !
! !    Output, real THETA, the angle between the two lines.
! !
!   real f1
!   real f2
!   real g1
!   real g2
!   real h1
!   real h2
!   real pdotq
!   real pnorm
!   real qnorm
!   real theta
!   real x01
!   real x02
!   real y01
!   real y02
!   real z01
!   real z02
! !
!   pdotq = f1 * f2 + g1 * g2 + h1 * h2
!   pnorm = sqrt ( f1 * f1 + g1 * g1 + h1 * h1 )
!   qnorm = sqrt ( f2 * f2 + g2 * g2 + h2 * h2 )
!
!   theta = arc_cosine ( pdotq / ( pnorm * qnorm ) )
!
!   return
! end subroutine
! subroutine lines_par_dist_3d ( f1, g1, h1, x01, y01, z01, f2, g2, h2, &
!   x02, y02, z02, dist )
! !
! !*******************************************************************************
! !
! !! LINES_PAR_DIST_3D finds the distance between two parametric lines in 3D.
! !
! !
! !  Formula:
! !
! !    The parametric form of a line in 3D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !      Z = Z0 + H * T
! !
! !  Warning:
! !
! !    This code does not work for parallel or near parallel lines.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real F1, G1, H1, X01, Y01, Z01, the parametric parameters
! !    of the first line.
! !
! !    Input, real F2, G2, H2, X02, Y02, Z02, the parametric parameters
! !    of the second line.
! !
! !    Output, real DIST, the distance between the two lines.
! !
!   real dist
!   real f1
!   real f2
!   real g1
!   real g2
!   real h1
!   real h2
!   real x01
!   real x02
!   real y01
!   real y02
!   real z01
!   real z02
! !
!   dist = abs ( ( x02 - x01 ) * ( g1 * h2 - g2 * h1 ) &
!              + ( y02 - y01 ) * ( h1 * f2 - h2 * f1 ) &
!              + ( z02 - z01 ) * ( f1 * g2 - f2 * g1 ) )  / &
!              ( ( f1 * g2 - f2 * g1 )**2 &
!              + ( g1 * h2 - g2 * h1 )**2 &
!              + ( h1 * f2 - h2 * f1 )**2 )
!
!   return
! end subroutine
! subroutine lines_par_int_2d ( f1, g1, x1, y1, f2, g2, x2, y2, t1, t2, &
!   xint, yint )
! !
! !*******************************************************************************
! !
! !! LINES_PAR_INT_2D determines where two parametric lines intersect in 2D.
! !
! !
! !  Formula:
! !
! !    The parametric form of a line in 2D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    26 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real F1, G1, X1, Y1, define the first parametric line.
! !
! !    Input, real F2, G2, X2, Y2, define the second parametric line.
! !
! !    Output, real DET, the determinant of the system.  If DET = 0,
! !    then the lines are parallel and no intersection was computed.
! !
! !    Output, real T1, T2, the T parameters on the first and second
! !    lines of the intersection point.
! !
! !    Output, real XINT, YINT, the intersection point.
! !
!   real det
!   real f1
!   real f2
!   real g1
!   real g2
!   real t1
!   real t2
!   real x1
!   real x2
!   real xint
!   real y1
!   real y2
!   real yint
! !
!   det = f2 * g1 - f1 * g2
!
!   if ( det == 0.0E+00 ) then
!     t1 = 0.0E+00
!     t2 = 0.0E+00
!     xint = 0.0E+00
!     yint = 0.0E+00
!   else
!     t1 = ( f2 * ( y2 - y1 ) - g2 * ( x2 - x1 ) ) / det
!     t2 = ( f1 * ( y2 - y1 ) - g1 * ( x2 - x1 ) ) / det
!     xint = x1 + f1 * t1
!     yint = y1 + g1 * t1
!   end if
!
!   return
! end subroutine
! subroutine lines_seg_dist_2d ( x1, y1, x2, y2, x3, y3, x4, y4, dist )
! !
! !*******************************************************************************
! !
! !! LINES_SEG_DIST_2D computes the distance between two line segments in 2D.
! !
! !
! !  Modified:
! !
! !    29 October 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, the endpoints of the first segment.
! !
! !    Input, real X3, Y3, X4, Y4, the endpoints of the second segment.
! !
! !    Output, real DIST, the distance between the line segments.
! !
!   real dist
!   real dl
!   real dm
!   real dr
!   real tl
!   real tm
!   real tmin
!   real tr
!   real x1
!   real x2
!   real x3
!   real x4
!   real xm
!   real xmin
!   real y1
!   real y2
!   real y3
!   real y4
!   real ym
!   real ymin
! !
! !  Label the left, middle and right of line 1 as L1, M1, R1.
! !
! !  Compute the distance from each point to the opposite line segment.
! !
!   call line_seg_point_dist_2d ( x3, y3, x4, y4, x1, y1, dl )
!
!   xm = 0.5E+00 * ( x1 + x2 )
!   ym = 0.5E+00 * ( y1 + y2 )
!
!   call line_seg_point_dist_2d ( x3, y3, x4, y4, xm, ym, dm )
!
!   call line_seg_point_dist_2d ( x3, y3, x4, y4, x2, y2, dr )
! !
! !  Now find the "theoretical" minimum of the distance function.
! !
!   tl = 0.0E+00
!   tm = 0.5E+00
!   tr = 1.0E+00
!   call minabs ( tl, dl, tm, dm, tr, dr, tmin, dist )
! !
! !  Evaluate the distance at the minimum, to account for
! !  a flattening of the absolute value function caused by parallel
! !  or coincident portions of the lines.
! !
!   xmin = x1 + tmin * ( x2 - x1 )
!   ymin = y1 + tmin * ( y2 - y1 )
!   call line_seg_point_dist_2d ( x3, y3, x4, y4, xmin, ymin, dist )
!
!   return
! end subroutine
! subroutine lines_seg_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, dist )
! !
! !*******************************************************************************
! !
! !! LINES_SEG_DIST_3D computes the distance between two line segments in 3D.
! !
! !
! !  Modified:
! !
! !    03 November 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, the endpoints of the first segment.
! !
! !    Input, real X3, Y3, Z3, X4, Y4, Z4, the endpoints of the second segment.
! !
! !    Output, real DIST, the distance between the line segments.
! !
!   real d1
!   real d2
!   real dist
!   real dl
!   real dm
!   real dr
!   real t1
!   real t2
!   real tl
!   real tm
!   real tmin
!   real tr
!   real x1
!   real x2
!   real x3
!   real x4
!   real xn1
!   real xn2
!   real y1
!   real y2
!   real y3
!   real y4
!   real yn1
!   real yn2
!   real z1
!   real z2
!   real z3
!   real z4
!   real zn1
!   real zn2
! !
! !  Find the nearest points on line 2 to the endpoints of line 1.
! !
!   call line_seg_point_near_3d ( x3, y3, z3, x4, y4, z4, x1, y1, z1, &
!     xn1, yn1, zn1, d1, t1 )
!
!   call line_seg_point_near_3d ( x3, y3, z3, x4, y4, z4, x2, y2, z2, &
!      xn2, yn2, zn2, d2, t2 )
!
!   if ( t1 == t2 ) then
!     call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, dist )
!     return
!   end if
! !
! !  On line 2, over the interval between the points nearest to line 1,
! !  the square of the distance of any point to line 1 is a quadratic function.
! !  Evaluate it at three points, and seek its local minimum.
! !
!   call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, xn1, yn1, zn1, dl )
!
!   call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, &
!     0.5E+00 * ( xn1 + xn2 ), 0.5E+00 * ( yn1 + yn2 ), &
!     0.5E+00 * ( zn1 + zn2 ), dm )
!
!   call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, xn2, yn2, zn2, dr )
!
!   tl = 0.0E+00
!   tm = 0.5E+00
!   tr = 1.0E+00
!
!   dl = dl * dl
!   dm = dm * dm
!   dr = dr * dr
!
!   call minquad ( tl, dl, tm, dm, tr, dr, tmin, dist )
!
!   dist = sqrt ( dist )
!
!   return
! end subroutine
! subroutine lines_seg_int_1d ( x1, x2, x3, x4, flag, x5, x6 )
! !
! !*******************************************************************************
! !
! !! LINES_SEG_INT_1D computes the intersection of two line segments in 1D.
! !
! !
! !  Modified:
! !
! !    07 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, X2, the endpoints of the first segment.
! !
! !    Input, real X3, X4, the endpoints of the second segment.
! !
! !    Output, integer FLAG, records the results.
! !    0, the line segments do not intersect.
! !    1, the line segments intersect.
! !
! !    Output, real X5, X6, the endpoints of the intersection segment.
! !    If FLAG = 0, X5 = X6 = 0.
! !
!   integer flag
!   real x1
!   real x2
!   real x3
!   real x4
!   real x5
!   real x6
!   real y1
!   real y2
!   real y3
!   real y4
! !
!   y1 = min ( x1, x2 )
!   y2 = max ( x1, x2 )
!   y3 = min ( x3, x4 )
!   y4 = max ( x3, x4 )
!
!   flag = 0
!   x5 = 0.0E+00
!   x6 = 0.0E+00
!
!   if ( y4 < y1 ) then
!     return
!   else if ( y3 > y2 ) then
!     return
!   end if
!
!   flag = 1
!   x5 = max ( y1, y3 )
!   x6 = min ( y2, y4 )
!
!   return
! end subroutine
! subroutine lines_seg_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5 )
! !
! !*******************************************************************************
! !
! !! LINES_SEG_INT_2D computes the intersection of two line segments in 2D.
! !
! !
! !  Modified:
! !
! !    28 September 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, the endpoints of the first segment.
! !
! !    Input, real X3, Y3, X4, Y4, the endpoints of the second segment.
! !
! !    Output, integer FLAG, records the results.
! !    0, the line segments do not intersect.
! !    1, the line segments intersect.
! !
! !    Output, real X5, Y5.
! !    If FLAG = 0, X5 = Y5 = 0.
! !    If FLAG = 1, then (X5,Y5) is a point of intersection.
! !
!   integer flag
!   integer ival
!   real u
!   real v
!   real x1
!   real x2
!   real x3
!   real x4
!   real x5
!   real y1
!   real y2
!   real y3
!   real y4
!   real y5
! !
!   x5 = 0.0E+00
!   y5 = 0.0E+00
! !
! !  Find the intersection of the two lines.
! !
!   call lines_exp_int_2d ( ival, x5, y5, x1, y1, x2, y2, x3, y3, x4, y4 )
!
!   if ( ival == 0 ) then
!     flag = 0
!     return
!   end if
! !
! !  Is the point on the first segment?
! !
!   call line_seg_contains_point_2d ( x1, y1, x2, y2, x5, y5, u, v )
!
!   if ( u < 0.0E+00 .or. 1.0E+00 < u .or. v > 0.001E+00 ) then
!     flag = 0
!     return
!   end if
! !
! !  Is the point on the second segment?
! !
!   call line_seg_contains_point_2d ( x3, y3, x4, y4, x5, y5, u, v )
!
!   if ( u < 0.0E+00 .or. 1.0E+00 < u .or. v > 0.001E+00 ) then
!     flag = 0
!     return
!   end if
!
!   flag = 1
!
!   return
! end subroutine
! subroutine loc2glob_3d ( cospitch, cosroll, cosyaw, locpts, globas, glopts, &
!   sinpitch, sinroll, sinyaw )
! !
! !*******************************************************************************
! !
! !! LOC2GLOB_3D converts from a local to global coordinate system in 3D.
! !
! !
! !  Discussion:
! !
! !    A global coordinate system is given.
! !
! !    A local coordinate system has been translated to the point with
! !    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
! !    a roll.
! !
! !    A point has local coordinates LOCPTS, and it is desired to know
! !    the point's global coordinates GLOPTS.
! !
! !    The transformation may be written as
! !
! !      GLOB = GLOBAS + N_YAW * N_PITCH * N_ROLL * LOC
! !
! !    where
! !
! !               (  cos(Yaw)   -sin(Yaw)        0      )
! !    N_YAW    = (  sin(Yaw)    cos(Yaw)        0      )
! !               (      0           0           1      )
! !
! !               (  cos(Pitch)      0       sin(Pitch) )
! !    N_PITCH =  (      0           1           0      )
! !               ( -sin(Pitch)      0       cos(Pitch) )
! !
! !               (      1           0           0      )
! !    N_ROLL =   (      0       cos(Roll)  -sin(Roll)  )
! !               (      0       sin(Roll)   cos(Roll)  )
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real COSPITCH, COSROLL, COSYAW, the cosines of the pitch,
! !    roll and yaw angles.
! !
! !    Input, real GLOBAS(3), the global coordinates of the base vector.
! !
! !    Output, real GLOPTS(3), the global coordinates of the point.
! !
! !    Input, real LOCPTS(3), the local coordinates of the point.
! !
! !    Input, real SINPITCH, SINROLL, SINYAW, the sines of the pitch,
! !    roll and yaw angles.
! !
!   real cospitch
!   real cosroll
!   real cosyaw
!   real globas(3)
!   real glopts(3)
!   real locpts(3)
!   real sinpitch
!   real sinroll
!   real sinyaw
! !
!   glopts(1) = globas(1) + (  cosyaw * cospitch ) * locpts(1) &
!     + (  cosyaw * sinpitch * sinroll - sinyaw * cosroll ) * locpts(2) &
!     + (  cosyaw * sinpitch * cosroll + sinyaw * sinroll ) * locpts(3)
!
!   glopts(2) = globas(2) + (  sinyaw * cospitch ) * locpts(1) &
!     + (  sinyaw * sinpitch * sinroll + cosyaw * cosroll ) * locpts(2) &
!     + (  sinyaw * sinpitch * cosroll - cosyaw * sinroll ) * locpts(3)
!
!   glopts(3) = globas(3) + ( -sinpitch ) * locpts(1) &
!     + (  cospitch * sinroll ) * locpts(2) &
!     + (  cospitch * cosroll ) * locpts(3)
!
!   return
! end subroutine
! subroutine minabs ( x1, y1, x2, y2, x3, y3, xmin, ymin )
! !
! !*******************************************************************************
! !
! !! MINABS finds a local minimum of F(X) = A * abs ( X ) + B.
! !
! !
! !  Modified:
! !
! !    28 October 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, real X1, Y1, X2, Y2, X3, Y3, are three sets of data
! !    of the form ( X, F(X) ).  The three X values must be distinct.
! !    On output, the data has been sorted so that X1 < X2 < X3,
! !    and the Y values have been rearranged accordingly.
! !
! !    Output, real XMIN, YMIN.  XMIN is a point within the interval
! !    spanned by X1, X2 and X3, at which F takes its local minimum
! !    value YMIN.
! !
!   real slope
!   real slope12
!   real slope13
!   real slope23
!   real temp
!   real x1
!   real x2
!   real x3
!   real xmin
!   real y1
!   real y2
!   real y3
!   real ymin
! !
! !  Refuse to deal with coincident data.
! !
!   if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'MINABS - Fatal error!'
!     write ( *, * ) '  X values are equal.'
!     return
!   end if
! !
! !  Sort the data.
! !
!   if ( x2 < x1 ) then
!     call r_swap ( x1, x2 )
!     call r_swap ( y1, y2 )
!   end if
!
!   if ( x3 < x1 ) then
!     call r_swap ( x1, x3 )
!     call r_swap ( y1, y3 )
!   end if
!
!   if ( x3 < x2 ) then
!     call r_swap ( x2, x3 )
!     call r_swap ( y2, y3 )
!   end if
! !
! !  Now determine the slopes.
! !
!   slope12 = ( y2 - y1 ) / ( x2 - x1 )
!   slope23 = ( y3 - y2 ) / ( x3 - x2 )
!   slope13 = ( y3 - y1 ) / ( x3 - x1 )
! !
! !  Case 1: Minimum must be at an endpoint.
! !
!   if ( slope12 >= slope13 .or. slope12 >= 0.0E+00 ) then
!
!     if ( y1 < y3 ) then
!       xmin = x1
!       ymin = y1
!     else
!       xmin = x3
!       ymin = y3
!     end if
! !
! !  Case 2: The curve decreases, and decreases faster than the line
! !  joining the endpoints.
! !
! !  Whichever of SLOPE12 and SLOPE23 is the greater in magnitude
! !  represents the actual slope of the underlying function.
! !  Find where two lines of that slope, passing through the
! !  endpoint data, intersect.
! !
!   else
!
!     slope = max ( abs ( slope12 ), slope23 )
!
!     xmin = 0.5E+00 * ( x1 + x3 + ( y1 - y3 ) / slope )
!     ymin = y1 - slope * ( xmin - x1 )
!
!   end if
!
!   return
! end subroutine
! subroutine minquad ( x1, y1, x2, y2, x3, y3, xmin, ymin )
! !
! !*******************************************************************************
! !
! !! MINQUAD finds a local minimum of F(X) = A * X * X + B * X + C.
! !
! !
! !  Note:
! !
! !    MINQUAD is primarily intended as a utility routine.
! !    The square of the distance function between a point
! !    and a line segment has the form of F(X).  Hence, we can seek
! !    the line on the second segment which minimizes the square of
! !    the distance to the other line segment.
! !
! !  Modified:
! !
! !    02 November 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, real X1, Y1, X2, Y2, X3, Y3, are three sets of data
! !    of the form ( X, F(X) ).  The three X values must be distinct.
! !    On output, the data has been sorted so that X1 < X2 < X3,
! !    and the Y values have been rearranged accordingly.
! !
! !    Output, real XMIN, YMIN.  XMIN is a point within the interval
! !    spanned by X1, X2 and X3, at which F takes its local minimum
! !    value YMIN.
! !
!   integer ierror
!   real x
!   real x1
!   real x2
!   real x3
!   real xleft
!   real xmin
!   real xrite
!   real y
!   real y1
!   real y2
!   real y3
!   real ymin
! !
! !  Refuse to deal with coincident data.
! !
!   if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'MINQUAD - Fatal error!'
!     write ( *, * ) '  X values are equal.'
!     return
!   end if
! !
! !  Find the interval endpoints.
! !
!   xleft = min ( x1, x2, x3 )
!   xrite = max ( x1, x2, x3 )
! !
! !  Find the minimizer and its function value, over the three input points.
! !
!   if ( y1 <= y2 .and. y1 <= y3 ) then
!     xmin = x1
!     ymin = y1
!   else if ( y2 <= y1 .and. y2 <= y3 ) then
!     xmin = x2
!     ymin = y2
!   else
!     xmin = x3
!     ymin = y3
!   end if
! !
! !  Find the minimizer and its function value over the real line.
! !
!   call parabola_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )
! !
! !  If F is linear, then take the already computed min.
! !
!   if ( ierror == 2 ) then
! !
! !  If F has a maximum, then take the already computed min.
! !
!   else if ( y > ymin ) then
! !
! !  If the minimizer is to the left, take the already computed min.
! !
!   else if ( x < xleft ) then
! !
! !  If the minimizer is to the right, take the already computed min.
! !
!   else if ( x > xrite ) then
!
!   else
!
!     xmin = x
!     ymin = y
!
!   end if
!
!   return
! end subroutine
! subroutine normal_01_sample ( x )
! !
! !*******************************************************************************
! !
! !! NORMAL_01_SAMPLE samples the standard Normal PDF.
! !
! !
! !  Discussion:
! !
! !    The standard normal distribution has mean 0 and standard
! !    deviation 1.
! !
! !  Method:
! !
! !    The Box-Muller method is used.
! !
! !  Modified:
! !
! !    01 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real X, a sample of the PDF.
! !
!   real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510E+00
! !
!   integer, save :: iset = -1
!   real v1
!   real v2
!   real x
!   real, save :: xsave = 0.0E+00
! !
!   if ( iset == -1 ) then
!     call random_seed
!     iset = 0
!   end if
!
!   if ( iset == 0 ) then
!
!     call random_number ( harvest = v1 )
!
!     if ( v1 <= 0.0E+00 ) then
!       write ( *, * ) ' '
!       write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
!       write ( *, * ) '  V1 <= 0.'
!       write ( *, * ) '  V1 = ', v1
!       stop
!     end if
!
!     call random_number ( harvest = v2 )
!
!     if ( v2 <= 0.0E+00 ) then
!       write ( *, * ) ' '
!       write ( *, * ) 'NORMAL_01_SAMPLE - Fatal error!'
!       write ( *, * ) '  V2 <= 0.'
!       write ( *, * ) '  V2 = ', v2
!       stop
!     end if
!
!     x = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * PI * v2 )
!
!     xsave = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * PI * v2 )
!
!     iset = 1
!
!   else
!
!     x = xsave
!     iset = 0
!
!   end if
!
!   return
! end subroutine
! subroutine norp2l_2d ( dist, x1, y1, x2, y2, x3, y3, x4, y4 )
! !
! !*******************************************************************************
! !
! !! NORP2L_2D finds two points on a line normal to a given line in 2D.
! !
! !
! !  Diagram:
! !
! !    P3
! !    |
! !    P1------P2
! !    |
! !    P4
! !
! !  Discussion:
! !
! !    The routine is given points P1 = (X1,Y1) and P2 = (X2,Y2)
! !    determining a line, and a distance DIST.  It returns the pair
! !    of points P3 = (X3,Y3) and P4 = (X4,Y4) which are a distance
! !    DIST from (X1,Y1), and lie on a line perpendicular to P2-P1.
! !
! !  Modified:
! !
! !    23 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real DIST, the nonnegative distance from (X1,Y1) to the computed
! !    points (X3,Y3) and (X4,Y4).
! !
! !    Input, real X1, Y1, X2, Y2.
! !    (X1,Y1) and (X2,Y2) are distinct points that define a line.
! !
! !    Output, real X3, Y3, X4, Y4.
! !    (X3,Y3) and (X4,Y4) are points which lie DIST units from
! !    (X1,Y1), in a direction normal to the line from (X1,Y1) to
! !    (X2,Y2).  The points are distinct unless DIST = 0.
! !
!   real a
!   real b
!   real dist
!   real length
!   real x1
!   real x2
!   real x3
!   real x4
!   real y1
!   real y2
!   real y3
!   real y4
! !
! !  Compute coefficients A and B of the line A*X+B*Y = C through
! !  the points (X1,Y1) and (X2,Y2).
! !
!   a = y2 - y1
!   b = x1 - x2
!   length = sqrt ( a * a + b * b )
! !
! !  Make sure the points are distinct.
! !
!   if ( length == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'NORP2L_2D - Fatal error!'
!     write ( *, * ) '  The line segment has zero length.'
!     stop
!   end if
! !
! !  The tangent vector to the line is (-b,a).
! !  The normal vector to the line is (a,b).
! !
! !  Compute the points (X1,Y1) +- DIST * (A,B) / Length(A,B).
! !
!   x3 = x1 + a * dist / length
!   y3 = y1 + b * dist / length
!
!   x4 = x1 - a * dist / length
!   y4 = y1 - b * dist / length
!
!   return
! end subroutine
! subroutine octahedron_shape_3d ( max_num, max_order, point_num, face_num, &
!   face_order, point_coord, face_point )
! !
! !*******************************************************************************
! !
! !! OCTAHEDRON_SHAPE_3D describes an octahedron in 3D.
! !
! !
! !  Discussion:
! !
! !    The vertices of the octahedron lie on the unit sphere.
! !
! !  Modified:
! !
! !    11 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAX_NUM, the maximum number of faces, and of
! !    points.  This quantity is used to dimension arrays.
! !    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
! !    not be set.
! !
! !    Input, integer MAX_ORDER, the maximum order of any face.
! !    This quantity is used to dimension arrays.  If FACE_ORDER
! !    exceeds MAX_ORDER, the arrays will not be set.
! !
! !    Output, integer POINT_NUM, the number of points in the shape.
! !
! !    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
! !    the X, Y and Z coordinates of point J.
! !
! !    Output, integer FACE_NUM, the number of faces in the shape.
! !
! !    Output, integer FACE_ORDER, the number of vertices per face.
! !
! !    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
! !    contains the index of the I-th point in the J-th face.  The
! !    points are listed in the counter-clockwise direction defined
! !    by the outward normal at the face.
! !
!   integer max_num
!   integer max_order
! !
!   integer face_num
!   integer face_order
!   integer face_point(max_order,max_num)
!   integer point_num
!   real point_coord(3,max_num)
! !
!   point_num = 6
!   face_num = 8
!   face_order = 3
! !
! !  Check.
! !
!   if ( point_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'OCTAHEDRON_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of vertices exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'OCTAHEDRON_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of faces exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_order > max_order ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'OCTAHEDRON_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Face order exceeds MAX_ORDER.'
!     return
!   end if
! !
! !  Set point coordinates.
! !
!   point_coord(1,1) =   0.0E+00
!   point_coord(2,1) =   0.0E+00
!   point_coord(3,1) = - 1.0E+00
!
!   point_coord(1,2) =   0.0E+00
!   point_coord(2,2) = - 1.0E+00
!   point_coord(3,2) =   0.0E+00
!
!   point_coord(1,3) =   1.0E+00
!   point_coord(2,3) =   0.0E+00
!   point_coord(3,3) =   0.0E+00
!
!   point_coord(1,4) =   0.0E+00
!   point_coord(2,4) =   1.0E+00
!   point_coord(3,4) =   0.0E+00
!
!   point_coord(1,5) =   0.0E+00
!   point_coord(2,5) = - 1.0E+00
!   point_coord(3,5) =   0.0E+00
!
!   point_coord(1,6) =   0.0E+00
!   point_coord(2,6) =   0.0E+00
!   point_coord(3,6) =   1.0E+00
! !
! !  Set faces.
! !
!   face_point(1,1) = 1
!   face_point(2,1) = 3
!   face_point(3,1) = 2
!
!   face_point(1,2) = 1
!   face_point(2,2) = 4
!   face_point(3,2) = 3
!
!   face_point(1,3) = 1
!   face_point(2,3) = 5
!   face_point(3,3) = 4
!
!   face_point(1,4) = 1
!   face_point(2,4) = 2
!   face_point(3,4) = 5
!
!   face_point(1,5) = 2
!   face_point(2,5) = 3
!   face_point(3,5) = 6
!
!   face_point(1,6) = 3
!   face_point(2,6) = 4
!   face_point(3,6) = 6
!
!   face_point(1,7) = 4
!   face_point(2,7) = 5
!   face_point(3,7) = 6
!
!   face_point(1,8) = 5
!   face_point(2,8) = 2
!   face_point(3,8) = 6
!
!   return
! end subroutine
! subroutine para_contains_point_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )
! !
! !*******************************************************************************
! !
! !! PARA_CONTAINS_POINT_2D determines if a point is inside a parallelogram in 2D.
! !
! !
! !  Diagram:
! !
! !         (X3,Y3).............
! !            /              .
! !           /              .
! !          /              .
! !    (X1,Y1)--------->(X2,Y2)
! !
! !  Modified:
! !
! !    04 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of three
! !    corners of the parallelogram.  (X1,Y1) should be directly connected
! !    to (X2,Y2) and (X3,Y3).
! !
! !    Input, real X, Y, the point to be checked.
! !
! !    Output, logical INSIDE, is .TRUE. if (X,Y) is inside the
! !    parallelogram, or on its boundary, and .FALSE. otherwise.
! !
!   logical inside
!   real p21dot
!   real p21normsq
!   real p31dot
!   real p31normsq
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!   p21normsq = enormsq0_2d ( x1, y1, x2, y2 )
!   p31normsq = enormsq0_2d ( x1, y1, x3, y3 )
!
!   p21dot = dot0_2d ( x1, y1, x2, y2, x, y )
!   p31dot = dot0_2d ( x1, y1, x3, y3, x, y )
!
!   if ( 0.0E+00 <= p21dot .and. p21dot <= p21normsq .and. &
!        0.0E+00 <= p31dot .and. p31dot <= p31normsq ) then
!     inside = .true.
!   else
!     inside = .false.
!   end if
!
!   return
! end subroutine
! subroutine para_contains_point_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x, y, z, inside )
! !
! !*******************************************************************************
! !
! !! PARA_CONTAINS_POINT_3D determines if a point is inside a parallelogram in 3D.
! !
! !
! !  Diagram:
! !
! !         (X2,Y2,Z2).............
! !            /                 .
! !           /                 .
! !          /                 .
! !    (X1,Y1,Z1)--------->(X3,Y3,Z3)
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
! !    the corners of the parallelogram.
! !
! !    Input, real X, Y, Z, the point to be checked.
! !
! !    Output, logical INSIDE, is .TRUE. if (X,Y,Z) is inside the
! !    parallelogram, or on its boundary, and .FALSE. otherwise.
! !    A slight amount of leeway is allowed for error, since a three
! !    dimensional point may lie exactly in the plane of the parallelogram,
! !    and yet be computationally slightly outside it.
! !
!   real, parameter :: tol = 0.00001E+00
! !
!   real dot
!   real dotb
!   real dott
!   logical inside
!   real v
!   real x
!   real x1
!   real x2
!   real x3
!   real xn12
!   real xn23
!   real xn31
!   real y
!   real y1
!   real y2
!   real y3
!   real yn12
!   real yn23
!   real yn31
!   real z
!   real z1
!   real z2
!   real z3
!   real zn12
!   real zn23
!   real zn31
! !
! !  Compute V3, the vector normal to V1 = P2-P1 and V2 = P3-P1.
! !
!   call cross0_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xn12, yn12, zn12 )
! !
! !  If the component of V = P-P1 in the V3 direction is too large,
! !  then it does not lie in the parallelogram.
! !
!   dot = ( x - x1 ) * xn12 + ( y - y1 ) * yn12 +  ( z - z1 ) * zn12
!
!   v = enorm0_3d ( x, y, z, x2, y2, z2 )
!
!   if ( abs ( dot ) > tol * ( 1.0E+00 + v ) ) then
!     inside = .false.
!     return
!   end if
! !
! !  Compute V23, the vector normal to V2 and V3.
! !
!   call cross_3d ( x3-x1, y3-y1, z3-z1, xn12, yn12, zn12, xn23, yn23, zn23 )
! !
! !  Compute ALPHA = ( V dot V23 ) / ( V1 dot V23 )
! !
!   dott = ( x - x1 ) * xn23 + ( y - y1 ) * yn23 + ( z - z1 ) * zn23
!
!   dotb = ( x2 - x1 ) * xn23 + ( y2 - y1 ) * yn23 + ( z2 - z1 ) * zn23
!
!   if ( dotb < 0.0E+00 ) then
!     dott = - dott
!     dotb = - dotb
!   end if
!
!   if ( dott < 0.0E+00 .or. dott > dotb ) then
!     inside = .false.
!     return
!   end if
! !
! !  Compute V31, the vector normal to V3 and V1.
! !
!   call cross_3d ( xn12, yn12, zn12, x2-x1, y2-y1, z2-z1, xn31, yn31, zn31 )
! !
! !  Compute BETA = ( V dot V31 ) / ( V2 dot V31 )
! !
!   dott = ( x - x1 ) * xn31 + ( y - y1 ) * yn31 + ( z - z1 ) * zn31
!
!   dotb = ( x3 - x1 ) * xn31 + ( y3 - y1 ) * yn31 + ( z3 - z1 ) * zn31
!
!   if ( dotb < 0.0E+00 ) then
!     dott = - dott
!     dotb = - dotb
!   end if
!
!   if ( dott < 0.0E+00 .or. dott > dotb ) then
!     inside = .false.
!     return
!   end if
! !
! !  V = ALPHA * V1 + BETA * V2, where both ALPHA and BETA are between
! !  0 and 1.
! !
!   inside = .true.
!
!   return
! end subroutine
! subroutine para_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z, &
!   dist )
! !
! !*******************************************************************************
! !
! !! PARA_POINT_DIST_3D: distance ( parallelogram, point ) in 3D.
! !
! !
! !  Diagram:
! !
! !         (X2,Y2,Z2).............
! !            /                 .
! !           /                 .
! !          /                 .
! !    (X1,Y1,Z1)--------->(X3,Y3,Z3)
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, determine the
! !    parallelogram, generated by the vectors from (X1,Y1) to (X2,Y2)
! !    and from (X1,Y1) to (X3,Y3).
! !
! !    Input, real X, Y, Z, the point which is to be checked.
! !
! !    Output, real DIST, the distance from the point to the
! !    parallelogram.  DIST is zero if the point lies exactly on the
! !    parallelogram.
! !
!   real dis13
!   real dis21
!   real dis34
!   real dis42
!   real dist
!   logical inside
!   real t
!   real temp
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real xn
!   real xp
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
!   real yn
!   real yp
!   real z
!   real z1
!   real z2
!   real z3
!   real z4
!   real zn
!   real zp
! !
! !  Compute P, the unit normal to X2-X1 and X3-X1:
! !
!   call cross0_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xp, yp, zp )
!
!   temp = sqrt ( xp * xp + yp * yp + zp * zp )
!
!   if ( temp == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PARA_POINT_DIST_3D - Fatal error!'
!     write ( *, * ) '  The normal vector is zero.'
!     stop
!   end if
!
!   xp = xp / temp
!   yp = yp / temp
!   zp = zp / temp
! !
! !  Find ( XN, YN, ZN), the nearest point to ( X, Y, Z ) in the plane.
! !
!   t = xp * ( x - x1 ) + yp * ( y - y1 ) + zp * ( z - z1 )
!
!   xn = x - xp * t
!   yn = y - yp * t
!   zn = z - zp * t
! !
! !  If ( XN, YN, ZN ) lies WITHIN the parallelogram, we're done.
! !
!   call para_contains_point_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!     x, y, z, inside )
!
!   if ( inside ) then
!     dist = enorm0_3d ( x, y, z, xn, yn, zn )
!     return
!   end if
! !
! !  Otherwise, find the distance between ( X, Y, Z ) and each of the
! !  four line segments that make up the boundary of the parallelogram.
! !
!   x4 = x2 + x3 - x1
!   y4 = y2 + y3 - y1
!   z4 = z2 + z3 - z1
!
!   call line_seg_point_dist_3d ( x1, y1, z1, x3, y3, z3, x, y, z, dis13 )
!
!   call line_seg_point_dist_3d ( x3, y3, z3, x4, y4, z4, x, y, z, dis34 )
!
!   call line_seg_point_dist_3d ( x4, y4, z4, x2, y2, z2, x, y, z, dis42 )
!
!   call line_seg_point_dist_3d ( x2, y2, z2, x1, y1, z1, x, y, z, dis21 )
!
!   dist = min ( dis13, dis34, dis42, dis21 )
!
!   return
! end subroutine
! subroutine parabola_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )
! !
! !*******************************************************************************
! !
! !! PARABOLA_EX finds the extremal point of a parabola determined by three points.
! !
! !
! !  Modified:
! !
! !    02 November 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
! !    on the parabola.  X1, X2 and X3 must be distinct.
! !
! !    Output, real X, Y, the X coordinate of the extremal point of the
! !    parabola, and the value of the parabola at that point.
! !
! !    Output, integer IERROR, error flag.
! !    0, no error.
! !    1, two of the X values are equal.
! !    2, the data lies on a straight line; there is no finite extremal point.
! !
!   real bot
!   integer ierror
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!   ierror = 0
!
!   if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
!     ierror = 1
!     return
!   end if
!
!   if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
!     x = x1
!     y = y1
!     return
!   end if
!
!   bot = ( x2 - x3 ) * y1 - ( x1 - x3 ) * y2 + ( x1 - x2 ) * y3
!
!   if ( bot == 0.0E+00 ) then
!     ierror = 2
!     return
!   end if
!
!   x = 0.5E+00 * (    x1 * x1 * ( y3 - y2 ) &
!                + x2 * x2 * ( y1 - y3 ) &
!                + x3 * x3 * ( y2 - y1 ) ) / bot
!
!   y = (    ( x - x2 ) * ( x - x3 ) * ( x2 - x3 ) * y1 &
!          - ( x - x1 ) * ( x - x3 ) * ( x1 - x3 ) * y2 &
!           + ( x - x1 ) * ( x - x2 ) * ( x1 - x2 ) * y3 ) / &
!           ( ( x1 - x2 ) * ( x2 - x3 ) * ( x1 - x3 ) )
!
!   return
! end subroutine
! subroutine parabola_ex2 ( x1, y1, x2, y2, x3, y3, x, y, a, b, c, ierror )
! !
! !*******************************************************************************
! !
! !! PARABOLA_EX2 finds the extremal point of a parabola determined by three points.
! !
! !
! !  Modified:
! !
! !    29 October 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
! !    on the parabola.  X1, X2 and X3 must be distinct.
! !
! !    Output, real X, Y, the X coordinate of the extremal point of the
! !    parabola, and the value of the parabola at that point.
! !
! !    Output, real A, B, C, the coefficients that define the parabola:
! !    P(X) = A * X * X + B * X + C.
! !
! !    Output, integer IERROR, error flag.
! !    0, no error.
! !    1, two of the X values are equal.
! !    2, the data lies on a straight line; there is no finite extremal
! !    point.
! !
!   real a
!   real b
!   real c
!   real det
!   integer ierror
!   real v(3,3)
!   real w(3,3)
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!   ierror = 0
!
!   if ( x1 == x2 .or. x2 == x3 .or. x3 == x1 ) then
!     ierror = 1
!     return
!   end if
!
!   if ( y1 == y2 .and. y2 == y3 .and. y3 == y1 ) then
!     x = x1
!     y = y1
!     return
!   end if
! !
! !  Set up the Vandermonde matrix.
! !
!   v(1,1) = 1.0E+00
!   v(1,2) = x1
!   v(1,3) = x1 * x1
!
!   v(2,1) = 1.0E+00
!   v(2,2) = x2
!   v(2,3) = x2 * x2
!
!   v(3,1) = 1.0E+00
!   v(3,2) = x3
!   v(3,3) = x3 * x3
! !
! !  Get the inverse.
! !
!   call rmat3_inverse ( v, w, det )
! !
! !  Compute the parabolic coefficients.
! !
!   c = w(1,1) * y1 + w(1,2) * y2 + w(1,3) * y3
!   b = w(2,1) * y1 + w(2,2) * y2 + w(2,3) * y3
!   a = w(3,1) * y1 + w(3,2) * y2 + w(3,3) * y3
! !
! !  Determine the extremal point.
! !
!   if ( a == 0.0E+00 ) then
!     ierror = 2
!     return
!   end if
!
!   x = - b / ( 2.0E+00 * a )
!   y = a * x * x + b * x + c
!
!   return
! end subroutine
!
! subroutine plane_exp_normal_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xn, yn, &
!   zn )
! !
! !*******************************************************************************
! !
! !! PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
! !
! !
! !  Definition:
! !
! !    The explicit form of a plane in 3D is
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of three points that constitute a line.  These points should not
! !    be identical, nor colinear.
! !
! !    Output, real XN, YN, ZN, the coordinates of the unit normal
! !    vector to the plane containing the three points.
! !
!   real temp
!   real x1
!   real x2
!   real x3
!   real xn
!   real y1
!   real y2
!   real y3
!   real yn
!   real z1
!   real z2
!   real z3
!   real zn
! !
! !  The cross product (P2-P1) x (P3-P1) is a vector normal to
! !  (P2-P1) and (P3-P1).
! !
!   call cross0_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xn, yn, zn )
!
!   temp = sqrt ( xn * xn + yn * yn + zn * zn )
!
!   if ( temp == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_EXP_NORMAL_3D - Fatal error!'
!     write ( *, * ) '  The plane is poorly defined.'
!     stop
!   else
!     xn = xn / temp
!     yn = yn / temp
!     zn = zn / temp
!   end if
!
!   return
! end subroutine
   subroutine plane_exp_point_dist_3d(x1, y1, z1, x2, y2, z2, x3, y3, z3, &
                                      x, y, z, dist)
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
      call plane_exp2imp_3d(x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d)

      call plane_imp_point_dist_3d(a, b, c, d, x, y, z, dist)

      return
   end subroutine
! subroutine plane_exp_project_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   xf, yf, zf, npnt, xo, yo, zo, xp, yp, zp, ivis )
! !
! !*******************************************************************************
! !
! !! PLANE_EXP_PROJECT_3D projects points through a point onto a plane in 3D.
! !
! !
! !  Definition:
! !
! !    The explicit form of a plane in 3D is
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
! !    coordinates of three points on the plane.
! !
! !    Input, real XF, YF, ZF, are the coordinates of the focus point.
! !
! !    Input, integer NPNT, the number of points to project.
! !
! !    Input, real XO(NPNT), YO(NPNT), ZO(NPNT), are the coordinates of
! !    the object points.
! !
! !    Output, real XP(NPNT), YP(NPNT), ZP(NPNT), are the coordinates of the
! !    projections of the object points through the focus point onto
! !    the plane.  XP, YP, and ZP may share the same memory as XO, YO,
! !    and ZO, in which case the projections will overwrite the original data.
! !
! !    Output, integer IVIS(NPNT), visibility indicator:
! !    3, the object was behind the plane;
! !    2, the object was already on the plane;
! !    1, the object was between the focus and the plane;
! !    0, the line from the object to the focus is parallel to the plane,
! !    so the object is "invisible".
! !    -1, the focus is between the object and the plane.  The object
! !    might be considered invisible.
! !
!   integer npnt
! !
!   real a
!   real alpha
!   real b
!   real beta
!   real c
!   real d
!   real disfo
!   real disfn
!   integer i
!   integer ivis(npnt)
!   real x1
!   real x2
!   real x3
!   real xf
!   real xn
!   real xo(npnt)
!   real xp(npnt)
!   real y1
!   real y2
!   real y3
!   real yf
!   real yn
!   real yo(npnt)
!   real yp(npnt)
!   real z1
!   real z2
!   real z3
!   real zf
!   real zn
!   real zo(npnt)
!   real zp(npnt)
! !
! !  Put the plane into ABCD form.
! !
!   call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
! !
! !  Get the nearest point on the plane to the focus.
! !
!   call plane_imp_point_near_3d ( a, b, c, d, xf, yf, zf, xn, yn, zn )
! !
! !  Get the distance from the focus to the plane.
! !
!   call points_dist_3d ( xf, yf, zf, xn, yn, zn, disfn )
! !
! !  If the focus lies in the plane, this is bad.  We could still
! !  project points that actually lie in the plane, but we'll
! !  just bail out.
! !
!   if ( disfn == 0.0E+00 ) then
!     ivis(1:npnt) = 0
!     xp(1:npnt) = xf
!     yp(1:npnt) = yf
!     zp(1:npnt) = zf
!     return
!   end if
! !
! !  Process the points.
! !
!   do i = 1, npnt
! !
! !  Get the distance from the focus to the object.
! !
!     call points_dist_3d ( xf, yf, zf, xo(i), yo(i), zo(i), disfo )
!
!     if ( disfo == 0.0E+00 ) then
!
!       ivis(i) = 0
!       xp(i) = xn
!       yp(i) = yn
!       zp(i) = zn
!
!     else
! !
! !  Compute ALPHA, the angle between (OBJECT-FOCUS) and (NEAREST-FOCUS).
! !
!       alpha = angle_rad_3d ( xo(i), yo(i), zo(i), xf, yf, zf, xn, yn, zn )
!
!       if ( cos ( alpha ) == 0.0E+00 ) then
!
!         ivis(i) = 0
!         xp(i) = xn
!         yp(i) = yn
!         zp(i) = zn
!
!       else
! !
! !  BETA is Dist(NEAREST-FOCUS) / ( Cos(ALPHA)*Dist(OBJECT-FOCUS) )
! !
!         beta = disfn / ( cos ( alpha ) * disfo )
!
!         if ( beta > 1.0E+00 ) then
!           ivis(i) = 1
!         else if ( beta == 1.0E+00 ) then
!           ivis(i) = 2
!         else if ( beta > 0.0E+00 ) then
!           ivis(i) = 3
!         else
!           ivis(i) = -1
!         end if
! !
! !  Set the projected point.
! !
!         xp(i) = xf + beta * ( xo(i) - xf )
!         yp(i) = yf + beta * ( yo(i) - yf )
!         zp(i) = zf + beta * ( zo(i) - zf )
!
!       end if
!
!     end if
!
!   end do
!
!   return
! end subroutine
   subroutine plane_exp2imp_3d(x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d)
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

      a = (y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1)
      b = (z2 - z1)*(x3 - x1) - (x2 - x1)*(z3 - z1)
      c = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1)
      d = -x2*a - y2*b - z2*c

      return
   end subroutine
! subroutine plane_exp2norm_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xp, yp, zp, &
!   xn, yn, zn )
! !
! !*******************************************************************************
! !
! !! PLANE_EXP2NORM_3D converts an explicit plane to normal form in 3D.
! !
! !
! !  Definition:
! !
! !    The explicit form of a plane in 3D is
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
! !
! !    The normal form of a plane in 3D is
! !
! !      (Xp,Yp,Zp), a point on the plane, and
! !      (Xn,Yn,Zn), the unit normal to the plane.
! !
! !  Modified:
! !
! !    02 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, X2, X3, Y3, Z3, are three points
! !    on the plane, which must be distinct, and not collinear.
! !
! !    Output, real XP, YP, ZP, a point on the plane.
! !
! !    Output, real XN, YN, ZN, the unit normal vector to the plane.
! !
!   real norm
!   real x1
!   real x2
!   real x3
!   real xn
!   real xp
!   real y1
!   real y2
!   real y3
!   real yn
!   real yp
!   real z1
!   real z2
!   real z3
!   real zn
!   real zp
!
!   xp = x1
!   yp = y1
!   zp = z1
!
!   xn = ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 )
!   yn = ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 )
!   zn = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )
!
!   norm = sqrt ( xn * xn + yn * yn + zn * zn )
!
!   if ( norm == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_EXP2NORM_3D - Fatal error!'
!     write ( *, * ) '  The normal vector is null.'
!     write ( *, * ) '  Two points coincide, or nearly so.'
!     stop
!   end if
!
!   xn = xn / norm
!   yn = yn / norm
!   zn = zn / norm
!
!   return
! end subroutine
! subroutine plane_grid_3d ( cor3, ierror, lines, maxcor3, maxline, &
!   ncor3, nline, x1, y1, z1, x2, y2, z2, x3, y3, z3 )
! !
! !*******************************************************************************
! !
! !! PLANE_GRID_3D computes points and lines making up a planar grid in 3D.
! !
! !
! !  Note:
! !
! !    The data format used is that of SGI Inventor.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real COR3(3,MAXCOR3), the X, Y, Z coordinates of points
! !    used in the grid.
! !
! !    Output, integer IERROR, error indicator.
! !    0, no error.
! !    1, more space for point coordinates is needed.
! !    2, more space for line data is needed.
! !
! !    Output, integer LINES(MAXLINE), the indices of points used in
! !    the lines of the grid.  Successive entries of LINES are joined
! !    by a line, unless an entry equals -1.  Note that indices begin
! !    with 0.
! !
! !    Input, integer MAXCOR3, the maximum number of points that can be
! !    stored.
! !
! !    Input, integer MAXLINE, the maximum number of line items that can
! !    be stored.
! !
! !    Input/output, integer NCOR3, the number of points stored in COR3.
! !
! !    On input, if NCOR3 is zero, then the data computed by this routine
! !    will be stored normally in COR3.  If NCOR3 is nonzero, then it
! !    is assumed that COR3 already contains some useful data.  The
! !    new data is appended to COR3.
! !
! !    On output, NCOR3 is increased by the number of points computed
! !    by this routine.
! !
! !    Input/output, integer NLINE, the number of line data items.
! !
! !    On input, if NLINE is zero, then the data computed by this routine
! !    will be stored normally in LINES.  If NLINE is nonzero, then it
! !    is assumed that LINES already contains some useful data.  The
! !    new data is appended to LINES.
! !
! !    On output, NLINE is increased by the number of points computed
! !    by this routine.
! !
! !    Input, real X1, Y1, Z1, X2, Y2, X2, X3, Y3, Z3, are three points
! !    on the plane, which must be distinct, and not collinear.
! !
!   integer, parameter :: nx = 5
!   integer, parameter :: ny = 5
! !
!   integer maxcor3
!   integer maxline
! !
!   real a
!   real amax
!   real amin
!   real b
!   real bmax
!   real bmin
!   real cor3(3,maxcor3)
!   real dot
!   integer i
!   integer ierror
!   integer j
!   integer k
!   integer lines(maxline)
!   integer nbase
!   integer ncor3
!   integer nline
!   real v1(3)
!   real v2(3)
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
!   real z1
!   real z2
!   real z3
! !
!   ierror = 0
!
!   nbase = ncor3
! !
! !  Compute the two basis vectors for the affine plane.
! !
!   v1(1) = x2 - x1
!   v1(2) = y2 - y1
!   v1(3) = z2 - z1
!
!   call vector_unit_nd ( 3, v1 )
!
!   v2(1) = x3 - x1
!   v2(2) = y3 - y1
!   v2(3) = z3 - z1
!
!   dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
!
!   v2(1:3) = v2(1:3) - dot * v1(1:3)
!
!   call vector_unit_nd ( 3, v2 )
! !
! !  Compute the (V1,V2) coordinate range of the data.
! !
!   amax = 0.0E+00
!   amin = 0.0E+00
!   bmax = 0.0E+00
!   bmin = 0.0E+00
!
!   do i = 1, ncor3
!
!     a = 0.0E+00
!     b = 0.0E+00
!     do j = 1, 3
!       a = a + v1(j) * cor3(j,i)
!       b = b + v2(j) * cor3(j,i)
!     end do
!
!     if ( i == 1 ) then
!       amax = a
!       amin = a
!       bmax = b
!       bmin = b
!     else
!       amax = max ( amax, a )
!       amin = min ( amin, a )
!       bmax = max ( bmax, b )
!       bmin = min ( bmin, b )
!     end if
!
!   end do
! !
! !  Generate the points we will use.
! !
!   if ( ncor3 + nx * ny > maxcor3 ) then
!     ierror = 1
!     return
!   end if
!
!   do j = 1, ny
!
!     b = ( real ( ny - j ) * bmin + real ( j - 1 ) * bmax ) / real ( ny - 1 )
!
!     do i = 1, nx
!
!       a = ( real ( nx - i ) * amin + real ( i - 1 ) * amax ) / real ( nx - 1 )
!
!       ncor3 = ncor3 + 1
!       cor3(1:3,ncor3) = a * v1(1:3) + b * v2(1:3)
!
!     end do
!
!   end do
! !
! !  Do the "horizontals".
! !
!   do i = 1, nx
!
!     do j = 1, ny
!
!       if ( nline >= maxline ) then
!         ierror = 2
!         return
!       end if
!
!       nline = nline + 1
!       lines(nline) = nbase + ( j - 1 ) * nx + i
!
!     end do
!
!     if ( nline >= maxline ) then
!       ierror = 2
!       return
!     end if
!
!     nline = nline + 1
!     lines(nline) = 0
!
!   end do
! !
! !  Do the "verticals".
! !
!   do j = 1, ny
!
!     do i = 1, nx
!
!       if ( nline >= maxline ) then
!         ierror = 2
!         return
!       end if
!
!       nline = nline + 1
!       lines(nline) = nbase + ( j - 1 ) * nx + i
!
!     end do
!
!     if ( nline >= maxline ) then
!       ierror = 2
!       return
!     end if
!
!     nline = nline + 1
!     lines(nline) = 0
!
!   end do
!
!   return
! end subroutine
! subroutine plane_imp_line_par_int_3d ( a, b, c, d, x0, y0, z0, f, g, h, &
!   intersect, x, y, z )
! !
! !*******************************************************************************
! !
! !! PLANE_IMP_LINE_PAR_INT_3D: intersection ( implicit plane, parametric line ) in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a plane in 3D is:
! !
! !      A * X + B * Y + C * Z + D = 0
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983, page 111.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, D, parameters that define the implicit plane.
! !
! !    Input, real X0, Y0, Z0, F, G, H, parameters that define the
! !    parametric line.
! !
! !    Output, logical INTERSECT, is TRUE if the line and the plane
! !    intersect, and false otherwise.
! !
! !    Output, real X, Y, Z, is a point of intersection of the line
! !    and the plane, if INTERSECT is TRUE.
! !
!   real, parameter :: TOL = 0.00001E+00
! !
!   real a
!   real b
!   real c
!   real d
!   real denom
!   real f
!   real g
!   real h
!   logical intersect
!   real norm1
!   real norm2
!   real t
!   real x
!   real x0
!   real y
!   real y0
!   real z
!   real z0
! !
! !  Check.
! !
!   norm1 = sqrt ( a * a + b * b + c * c )
!
!   if ( norm1 == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
!     write ( *, * ) '  The plane normal vector is null.'
!     stop
!   end if
!
!   norm2 = sqrt ( f * f + g * g + h * h )
!
!   if ( norm2 == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
!     write ( *, * ) '  The line direction vector is null.'
!     stop
!   end if
!
!   denom = a * f + b * g + c * h
! !
! !  The line and the plane may be parallel.
! !
!   if ( abs ( denom ) < TOL * norm1 * norm2 ) then
!
!     if ( a * x0 + b * y0 + c * z0 + d == 0.0E+00 ) then
!       intersect = .TRUE.
!       x = x0
!       y = y0
!       z = z0
!     else
!       intersect = .FALSE.
!       x = 0.0E+00
!       y = 0.0E+00
!       z = 0.0E+00
!     end if
! !
! !  If they are not parallel, they must intersect.
! !
!   else
!
!     intersect = .TRUE.
!     t = - ( a * x0 + b * y0 + c * z0 + d ) / denom
!     x = x0 + t * f
!     y = y0 + t * g
!     z = z0 + t * h
!
!   end if
!
!   return
! end subroutine
! subroutine plane_imp_line_seg_near_3d ( x1, y1, z1, x2, y2, z2, &
!   a, b, c, d, dist, xp, yp, zp, xls, yls, zls )
! !
! !*******************************************************************************
! !
! !! PLANE_IMP_LINE_SEG_NEAR_3D: nearest ( implicit plane, line segment ) in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a plane in 3D is:
! !
! !      A * X + B * Y + C * Z + D = 0
! !
! !  Modified:
! !
! !    17 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, the endpoints of the line
! !    segment.
! !
! !    Input, real A, B, C, D, the parameters that define the implicit
! !    plane.
! !
! !    Output, real DIST, the distance between the line segment and
! !    the plane.
! !
! !    Output, real XP, YP, ZP, the nearest point on the plane.
! !
! !    Output, real XLS, YLS, ZLS, the nearest point on the line segment
! !    to the plane.  If DIST is zero, the (XLS,YLS,ZLS) is a point of
! !    intersection of the plane and the line segment.
! !
!   real a
!   real alpha
!   real an
!   real b
!   real bn
!   real c
!   real cn
!   real d
!   real dist
!   real dn
!   real norm
!   real p1
!   real p2
!   real x1
!   real x2
!   real xls
!   real xp
!   real y1
!   real y2
!   real yls
!   real yp
!   real z1
!   real z2
!   real zls
!   real zp
! !
!   xls = 0.0E+00
!   yls = 0.0E+00
!   zls = 0.0E+00
!   xp = 0.0E+00
!   yp = 0.0E+00
!   zp = 0.0E+00
!
!   norm = sqrt ( a * a + b * b + c * c )
!
!   if ( norm == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_IMP_LINE_SEG_NEAR_3D - Fatal error!'
!     write ( *, * ) '  Plane normal vector is null.'
!     stop
!   end if
! !
! !  The normalized coefficients allow us to compute the (signed) distance.
! !
!   an = a / norm
!   bn = b / norm
!   cn = c / norm
!   dn = d / norm
! !
! !  If the line segment is actually a point, then the answer is easy.
! !
!   if ( x1 == x2 .and. y1 == y2 .and. z1 == z2 ) then
!
!     p1 = an * x1 + bn * y1 + cn * z1 + dn
!     dist = abs ( p1 )
!     xls = x1
!     yls = y1
!     zls = z1
!     xp = xls - an * p1
!     yp = yls - bn * p1
!     zp = zls - cn * p1
!     return
!
!   end if
! !
! !  Compute the projections of the two points onto the normal vector.
! !
!   p1 = an * x1 + bn * y1 + cn * z1 + dn
!   p2 = an * x2 + bn * y2 + cn * z2 + dn
! !
! !  If these have the same sign, then the line segment does not
! !  cross the plane, and one endpoint is the nearest point.
! !
!   if ( ( p1 > 0.0E+00 .and. p2 > 0.0E+00 ) .or. &
!        ( p1 < 0.0E+00 .and. p2 < 0.0E+00 ) ) then
!
!     p1 = abs ( p1 )
!     p2 = abs ( p2 )
!
!     if ( p1 < p2 ) then
!       xls = x1
!       yls = y1
!       zls = z1
!       xp = xls - an * p1
!       yp = yls - bn * p1
!       zp = zls - cn * p1
!       dist = p1
!     else
!       xls = x2
!       yls = y2
!       zls = z2
!       dist = p2
!       xp = xls - an * p2
!       yp = yls - bn * p2
!       zp = zls - cn * p2
!     end if
! !
! !  If the projections differ in sign, the line segment crosses the plane.
! !
!   else
!
!     if ( p1 == 0.0E+00 ) then
!       alpha = 0.0E+00
!     else if ( p2 == 0.0E+00 ) then
!       alpha = 1.0E+00
!     else
!       alpha = p2 / ( p2 - p1 )
!     end if
!
!     xls = alpha * x1 + ( 1.0E+00 - alpha ) * x2
!     yls = alpha * y1 + ( 1.0E+00 - alpha ) * y2
!     zls = alpha * z1 + ( 1.0E+00 - alpha ) * z2
!     xp = xls
!     yp = yls
!     zp = zls
!
!     dist = 0.0E+00
!
!   end if
!
!   return
! end subroutine
   subroutine plane_imp_point_dist_3d(a, b, c, d, x, y, z, dist)
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
      real a
      real b
      real c
      real d
      real dist
      real x
      real y
      real z
!
      dist = abs((a*x + b*y + c*z + d))/sqrt(a*a + b*b + c*c)

      return
   end subroutine
! subroutine plane_imp_point_dist_signed_3d ( a, b, c, d, x, y, z, dist_signed )
! !
! !*******************************************************************************
! !
! !! PLANE_IMP_POINT_DIST_SIGNED_3D: signed distance ( implicit plane, point) in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a plane in 3D is:
! !
! !      A * X + B * Y + C * Z + D = 0
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, D, determine the equation of the
! !    plane, which is:
! !      A*X + B*Y + C*Z + D = 0.
! !
! !    Input, real X, Y, Z, the coordinates of the point.
! !
! !    Output, real DIST_SIGNED, the signed distance from the point to
! !    the plane.
! !
!   real a
!   real b
!   real c
!   real d
!   real dist_signed
!   real x
!   real y
!   real z
! !
!   dist_signed = ( a * x + b * y + c * z + d ) / &
!     ( - sign ( 1.0E+00, d ) * sqrt ( a * a + b * b + c * c ) )
!
!   return
! end subroutine
! subroutine plane_imp_point_near_3d ( a, b, c, d, x, y, z, xn, yn, zn )
! !
! !*******************************************************************************
! !
! !! PLANE_IMP_POINT_NEAR_3D: nearest point on a implicit plane to a point in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a plane in 3D is:
! !
! !      A * X + B * Y + C * Z + D = 0
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, D, coefficients that define the plane as
! !    the set of points for which A*X+B*Y+C*Z+D = 0.
! !
! !    Input, real X, Y, Z, the coordinates of the point.
! !
! !    Output, real XN, YN, ZN, the coordinates of the nearest point on
! !    the plane.
! !
!   real a
!   real b
!   real c
!   real d
!   real t
!   real x
!   real xn
!   real y
!   real yn
!   real z
!   real zn
! !
!   if ( a == 0.0E+00 .and. b == 0.0E+00 .and. c == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_IMP_POINT_NEAR_3D - Fatal error!'
!     write ( *, * ) '  A = B = C = 0.'
!     stop
!   end if
! !
! !  The normal N to the plane is (A,B,C).
! !
! !  The line defined by (XN-X)/A = (YN-Y)/B = (ZN-Z)/C = T
! !  goes through (X,Y,Z) and is parallel to N.
! !
! !  Solving for the point (XN,YN,ZN) we get
! !
! !    XN = A*T+X
! !    YN = B*T+Y
! !    ZN = C*T+Z
! !
! !  Now place these values in the equation for the plane:
! !
! !    A*(A*T+X) + B*(B*T+Y) + C*(C*T+Z) + D = 0
! !
! !  and solve for T:
! !
! !    T = (-A*X-B*Y-C*Z-D) / (A * A + B * B + C * C )
! !
!   t = - ( a * x + b * y + c * z + d ) / ( a * a + b * b + c * c )
!
!   xn = x + a * t
!   yn = y + b * t
!   zn = z + c * t
!
!   return
! end subroutine
! subroutine plane_imp_triangle_int_3d ( x1, y1, z1, x2, y2, z2, &
!   x3, y3, z3, a, b, c, d, num_int, x, y, z )
! !
! !*******************************************************************************
! !
! !! PLANE_IMP_TRIANGLE_INT_3D: intersection ( implicit plane, triangle ) in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a plane in 3D is:
! !
! !      A * X + B * Y + C * Z + D = 0
! !
! !  Discussion:
! !
! !    There may be 0, 1, 2 or 3 points of intersection returned.
! !
! !    If two intersection points are returned, then the entire line
! !    between them comprises points of intersection.
! !
! !    If three intersection points are returned, then all points of
! !    the triangle intersect the plane.
! !
! !  Modified:
! !
! !    02 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of the vertices of the triangle.
! !
! !    Input, real A, B, C, D, the parameters that define the implicit plane.
! !
! !    Output, integer NUM_INT, the number of intersection points returned.
! !
! !    Output, real X(3), Y(3), Z(3), the coordinates of the intersection points.
! !
!   real a
!   real b
!   real c
!   real d
!   real dist1
!   real dist2
!   real dist3
!   integer num_int
!   real x(3)
!   real x1
!   real x2
!   real x3
!   real y(3)
!   real y1
!   real y2
!   real y3
!   real z(3)
!   real z1
!   real z2
!   real z3
! !
!   num_int = 0
! !
! !  Compute the signed distances between the vertices and the plane.
! !
!   dist1 = a * x1 + b * y1 + c * z1 + d
!   dist2 = a * x2 + b * y2 + c * z2 + d
!   dist3 = a * x3 + b * y3 + c * z3 + d
! !
! !  Consider any zero distances.
! !
!   if ( dist1 == 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x1
!     y(num_int) = y1
!     z(num_int) = z1
!
!   end if
!
!   if ( dist2 == 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x2
!     y(num_int) = y2
!     z(num_int) = z2
!
!   end if
!
!   if ( dist3 == 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x3
!     y(num_int) = y3
!     z(num_int) = z3
!
!   end if
! !
! !  If 2 or 3 of the nodes intersect, we're already done.
! !
!   if ( num_int >= 2 ) then
!     return
!   end if
! !
! !  If one node intersects, then we're done unless the other two
! !  are of opposite signs.
! !
!   if ( num_int == 1 ) then
!
!     if ( dist1 == 0.0E+00 ) then
!
!       call plane_triangle_int_add_3d ( x2, y2, z2, x3, y3, z3, &
!         dist2, dist3, num_int, x, y, z )
!
!     else if ( dist2 == 0.0E+00 ) then
!
!       call plane_triangle_int_add_3d ( x1, y1, z1, x3, y3, z3, &
!         dist1, dist3, num_int, x, y, z )
!
!     else if ( dist3 == 0.0E+00 ) then
!
!       call plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, &
!         dist1, dist2, num_int, x, y, z )
!
!     end if
!
!     return
!
!   end if
! !
! !  All nodal distances are nonzero, and there is at least one
! !  positive and one negative.
! !
!   if ( dist1 * dist2 < 0.0E+00 .and. dist1 * dist3 < 0.0E+00 ) then
!
!     call plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, &
!       dist1, dist2, num_int, x, y, z )
!
!     call plane_triangle_int_add_3d ( x1, y1, z1, x3, y3, z3, &
!       dist1, dist3, num_int, x, y, z )
!
!   else if ( dist2 * dist1 < 0.0E+00 .and.  dist2 * dist3 < 0.0E+00 ) then
!
!     call plane_triangle_int_add_3d ( x2, y2, z2, x1, y1, z1, &
!       dist2, dist1, num_int, x, y, z )
!
!     call plane_triangle_int_add_3d ( x2, y2, z2, x3, y3, z3, &
!       dist2, dist3, num_int, x, y, z )
!
!   else if ( dist3 * dist1 < 0.0E+00 .and. dist3 * dist2 < 0.0E+00 ) then
!
!     call plane_triangle_int_add_3d ( x3, y3, z3, x1, y1, z1, &
!       dist3, dist1, num_int, x, y, z )
!
!     call plane_triangle_int_add_3d ( x3, y3, z3, x2, y2, z2, &
!       dist3, dist2, num_int, x, y, z )
!
!   end if
!
!   return
! end subroutine
! subroutine plane_imp_triangle_near_3d ( x1, y1, z1, x2, y2, z2, &
!   x3, y3, z3, a, b, c, d, dist, num_near, x, y, z )
! !
! !*******************************************************************************
! !
! !! PLANE_IMP_TRIANGLE_NEAR_3D: nearest ( implicit plane, triangle ) in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a plane in 3D is:
! !
! !      A * X + B * Y + C * Z + D = 0
! !
! !  Comments:
! !
! !    Please see to it that the underlying distance routine always returns
! !    one of the endpoints if the entire line segment is at zero distance.
! !
! !  Modified:
! !
! !    18 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of the vertices of the triangle.
! !
! !    Input, real A, B, C, D, the parameters that define the implicit plane.
! !
! !    Output, real DIST, the distance between the triangle and the plane.
! !
! !    Output, integer NUM_NEAR, the number of nearest points returned.
! !
! !    Output, real X(6), Y(6), Z(6), a collection of nearest points.
! !
! !    If DIST = 0, then each point is a point of intersection, and there
! !    will be at most 3 such points returned.
! !
! !    If DIST > 0, then the points are listed in pairs, with the first
! !    being on the triangle, and the second on the plane.  Two points will
! !    be listed in the most common case, but possibly 4 or 6.
! !
!   real a
!   real b
!   real c
!   real d
!   real dist
!   real dist12
!   real dist23
!   real dist31
!   integer num_near
!   real x(6)
!   real x1
!   real x2
!   real x3
!   real xp
!   real xt
!   real y(6)
!   real y1
!   real y2
!   real y3
!   real yp
!   real yt
!   real z(6)
!   real z1
!   real z2
!   real z3
!   real zp
!   real zt
! !
!   num_near = 0
! !
! !  Consider the line segment P1 - P2.
! !
!   call plane_imp_line_seg_near_3d ( x1, y1, z1, x2, y2, z2, &
!     a, b, c, d, dist12, xp, yp, zp, xt, yt, zt )
!
!   dist = dist12
!
!   num_near = num_near + 1
!   x(num_near) = xt
!   y(num_near) = yt
!   z(num_near) = zt
!
!   if ( dist12 > 0.0E+00 ) then
!     num_near = num_near + 1
!     x(num_near) = xp
!     y(num_near) = yp
!     z(num_near) = zp
!   end if
! !
! !  Consider the line segment P2 - P3.
! !
!   call plane_imp_line_seg_near_3d ( x2, y2, z2, x3, y3, z3, &
!     a, b, c, d, dist23, xp, yp, zp, xt, yt, zt )
!
!   if ( dist23 < dist ) then
!
!     num_near = 0
!     dist = dist23
!
!     num_near = num_near + 1
!     x(num_near) = xt
!     y(num_near) = yt
!     z(num_near) = zt
!
!     if ( dist23 > 0.0E+00 ) then
!       num_near = num_near + 1
!       x(num_near) = xp
!       y(num_near) = yp
!       z(num_near) = zp
!     end if
!
!   else if ( dist23 == dist ) then
!
!     num_near = num_near + 1
!     x(num_near) = xt
!     y(num_near) = yt
!     z(num_near) = zt
!
!     if ( dist23 > 0.0E+00 ) then
!       num_near = num_near + 1
!       x(num_near) = xp
!       y(num_near) = yp
!       z(num_near) = zp
!     end if
!
!   end if
! !
! !  Consider the line segment P3 - P1.
! !
!   call plane_imp_line_seg_near_3d ( x3, y3, z3, x1, y1, z1, &
!     a, b, c, d, dist31, xp, yp, zp, xt, yt, zt )
!
!   if ( dist31 < dist ) then
!
!     num_near = 0
!     dist = dist31
!
!     num_near = num_near + 1
!     x(num_near) = xt
!     y(num_near) = yt
!     z(num_near) = zt
!
!     if ( dist31 > 0.0E+00 ) then
!       num_near = num_near + 1
!       x(num_near) = xp
!       y(num_near) = yp
!       z(num_near) = zp
!     end if
!
!   else if ( dist31 == dist ) then
!
!     num_near = num_near + 1
!     x(num_near) = xt
!     y(num_near) = yt
!     z(num_near) = zt
!
!     if ( dist31 > 0.0E+00 ) then
!       num_near = num_near + 1
!       x(num_near) = xp
!       y(num_near) = yp
!       z(num_near) = zp
!     end if
!
!   end if
!
!   return
! end subroutine
! subroutine plane_imp2exp_3d ( a, b, c, d, x1, y1, z1, x2, y2, z2, x3, y3, z3 )
! !
! !*******************************************************************************
! !
! !! PLANE_IMP2EXP_3D converts an implicit plane to explicit form in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a plane in 3D is
! !
! !      A * X + B * Y + C * Z + D = 0.
! !
! !    The explicit form of a plane in 3D is
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
! !
! !  Modified:
! !
! !    27 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, D, parameters that define the implicit plane.
! !
! !    Output, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
! !    three points that lie on the plane.
! !
!   real a
!   real b
!   real c
!   real d
!   real x1
!   real x2
!   real x3
!   real xn
!   real xp
!   real y1
!   real y2
!   real y3
!   real yn
!   real yp
!   real z1
!   real z2
!   real z3
!   real zn
!   real zp
! !
!   call plane_imp2norm_3d ( a, b, c, d, xp, yp, zp, xn, yn, zn )
!
!   call plane_norm2exp_3d ( xp, yp, zp, xn, yn, zn, x1, y1, z1, &
!     x2, y2, z2, x3, y3, z3 )
!
!   return
! end subroutine
! subroutine plane_imp2norm_3d ( a, b, c, d, xp, yp, zp, xn, yn, zn )
! !
! !*******************************************************************************
! !
! !! PLANE_IMP2NORM_3D converts an implicit plane to normal form in 3D.
! !
! !
! !  Definition:
! !
! !    The implicit form of a plane in 3D is
! !
! !      A * X + B * Y + C * Z + D = 0.
! !
! !    The normal form of a plane in 3D is
! !
! !      (Xp,Yp,Zp), a point on the plane, and
! !      (Xn,Yn,Zn), the unit normal to the plane.
! !
! !  Modified:
! !
! !    02 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, D, parameters that define the implicit plane.
! !
! !    Output, real XP, YP, ZP, a point on the plane.
! !
! !    Output, real XN, YN, ZN, the unit normal vector to the plane.
! !
!   real a
!   real b
!   real c
!   real d
!   real norm
!   real xn
!   real xp
!   real yn
!   real yp
!   real zn
!   real zp
! !
!   norm = sqrt ( a * a + b * b + c * c )
!
!   if ( norm == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_IMP2NORM_3D - Fatal error!'
!     write ( *, * ) '  The (A,B,C) has zero norm.'
!     stop
!   end if
!
!   xn = a / norm
!   yn = b / norm
!   zn = c / norm
!
!   if ( a /= 0.0E+00 ) then
!     xp = - d / a
!     yp = 0.0E+00
!     zp = 0.0E+00
!   else if ( b /= 0.0E+00 ) then
!     xp = 0.0E+00
!     yp = - d / b
!     zp = 0.0E+00
!   else if ( c /= 0.0E+00 ) then
!     xp = 0.0E+00
!     yp = 0.0E+00
!     zp = - d / c
!   else
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_IMP2NORM_3D - Fatal error!'
!     write ( *, * ) '  The (A,B,C) vector is null.'
!     stop
!   end if
!
!   return
! end subroutine
! subroutine plane_norm_basis_3d ( xp, yp, zp, xn, yn, zn, xq, yq, zq, &
!   xr, yr, zr )
! !
! !*******************************************************************************
! !
! !! PLANE_NORM_BASIS_3D finds two perpendicular vectors in a plane in 3D.
! !
! !
! !  Discussion:
! !
! !    Given a plane in point, normal form P = (XP,YP,ZP) and N = (XN,YN,ZN),
! !    any point in that plane can be described in terms of the point P
! !    plus a weighted sum of two vectors Q = (XQ,YQ,ZQ) and R = (XR,YR,ZR):
! !
! !      (X,Y,Z) = (XP,YP,ZP) + a * (XQ,YQ,ZQ) + b * (XR,YR,ZR).
! !
! !    The vector Q has unit length, and is perpendicular to P and R;
! !    the vector R has unit length, and is perpendicular to P and Q.
! !
! !  Modified:
! !
! !    24 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real XP, YP, ZP, a point on the plane.  (Actually, we never
! !    need to know these values to do the calculation!)
! !
! !    Input, real XN, YN, ZN, a normal vector N to the plane.  The
! !    vector must not have zero length, but it is not necessary for N
! !    to have unit length.
! !
! !    Output, real XQ, YQ, ZQ, a vector of unit length, perpendicular
! !    to the vector N and the vector R.
! !
! !    Output, real XR, YR, ZR, a vector of unit length, perpendicular
! !    to the vector N and the vector Q.
! !
!   real dot
!   real min_com
!   real norm_n
!   real norm_q
!   real xn
!   real xp
!   real xq
!   real xr
!   real yn
!   real yp
!   real yq
!   real yr
!   real zn
!   real zp
!   real zq
!   real zr
! !
! !  Compute the length of N = (XN,YN,ZN).
! !
!   norm_n = sqrt ( xn * xn + yn * yn + zn * zn )
!
!   if ( norm_n == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'PLANE_NORM_BASIS_3D - Fatal error!'
!     write ( *, * ) '  The normal vector is 0.'
!     stop
!   end if
! !
! !  To find a vector distinct from N, find the minimum component
! !  of N, and set the corresponding component of Q to 1, the
! !  rest to zero.
! !
!   xq = 0.0E+00
!   yq = 0.0E+00
!   zq = 0.0E+00
!
!   min_com = abs ( xn )
!
!   if ( abs ( yn ) < min_com ) then
!     min_com = abs ( yn )
!   end if
!
!   if ( abs ( zn ) < min_com ) then
!     min_com = abs ( zn )
!   end if
!
!   if ( min_com == abs ( xn ) ) then
!     xq = 1.0E+00
!
!   else if ( min_com == abs ( yn ) ) then
!     yq = 1.0E+00
!   else if ( min_com == abs ( zn ) ) then
!     zq = 1.0E+00
!   end if
! !
! !  Now subtract off the component of Q in the direction of N,
! !  computing Q = Q - Q dot N / || N ||,
! !  and then normalize.
! !
!   dot = ( xq * xn + yq * yn + zq * zn ) / norm_n
!
!   xq = xq - dot * xn / norm_n
!   yq = yq - dot * yn / norm_n
!   zq = zq - dot * zn / norm_n
!
!   norm_q = sqrt ( xq * xq + yq * yq + zq * zq )
!
!   xq = xq / norm_q
!   yq = yq / norm_q
!   zq = zq / norm_q
! !
! !  Now just take the cross product N x Q to get the R vector.
! !  Plus, if we did things right, R will already have unit length.
! !
!   xr = ( yn * zq - zn * yq ) / norm_n
!   yr = ( zn * xq - xn * zq ) / norm_n
!   zr = ( xn * yq - yn * xq ) / norm_n
!
!   return
! end subroutine
! subroutine plane_norm_triangle_int_3d ( x1, y1, z1, x2, y2, z2, &
!   x3, y3, z3, xp, yp, zp, xn, yn, zn, num_int, x, y, z )
! !
! !*******************************************************************************
! !
! !! PLANE_NORM_TRIANGLE_INT_3D: intersection ( normal plane, triangle ) in 3D.
! !
! !
! !  Definition:
! !
! !    The normal form of a plane in 3D is:
! !
! !      (Xp,Yp,Zp) is a point on the plane,
! !      (Xn,Yn,Zn) is a normal vector to the plane.
! !
! !  Discussion:
! !
! !    There may be 0, 1, 2 or 3 points of intersection returned.
! !
! !    If two intersection points are returned, then the entire line
! !    between them comprises points of intersection.
! !
! !    If three intersection points are returned, then all points of
! !    the triangle intersect the plane.
! !
! !  Modified:
! !
! !    03 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of the vertices of the triangle.
! !
! !    Input, real XP, YP, ZP, a point on the plane.
! !
! !    Input, real XN, YN, ZN, a normal vector to the plane.
! !
! !    Output, integer NUM_INT, the number of intersection points returned.
! !
! !    Output, real X(3), Y(3), Z(3), the coordinates of the intersection points.
! !
!   real d
!   real dist1
!   real dist2
!   real dist3
!   integer num_int
!   real x(3)
!   real x1
!   real x2
!   real x3
!   real xn
!   real xp
!   real y(3)
!   real y1
!   real y2
!   real y3
!   real yn
!   real yp
!   real z(3)
!   real z1
!   real z2
!   real z3
!   real zn
!   real zp
! !
!   num_int = 0
! !
! !  Compute the signed distances between the vertices and the plane.
! !
!   d = - xn * xp - yn * yp - zn * zp
!   dist1 = xn * x1 + yn * y1 + zn * z1 + d
!   dist2 = xn * x2 + yn * y2 + zn * z2 + d
!   dist3 = xn * x3 + yn * y3 + zn * z3 + d
! !
! !  Consider any zero distances.
! !
!   if ( dist1 == 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x1
!     y(num_int) = y1
!     z(num_int) = z1
!
!   end if
!
!   if ( dist2 == 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x2
!     y(num_int) = y2
!     z(num_int) = z2
!
!   end if
!
!   if ( dist3 == 0.0E+00 ) then
!
!     num_int = num_int + 1
!     x(num_int) = x3
!     y(num_int) = y3
!     z(num_int) = z3
!
!   end if
! !
! !  If 2 or 3 of the nodes intersect, we're already done.
! !
!   if ( num_int >= 2 ) then
!     return
!   end if
! !
! !  If one node intersects, then we're done unless the other two
! !  are of opposite signs.
! !
!   if ( num_int == 1 ) then
!
!     if ( dist1 == 0.0E+00 ) then
!
!       call plane_triangle_int_add_3d ( x2, y2, z2, x3, y3, z3, &
!         dist2, dist3, num_int, x, y, z )
!
!     else if ( dist2 == 0.0E+00 ) then
!
!       call plane_triangle_int_add_3d ( x1, y1, z1, x3, y3, z3, &
!         dist1, dist3, num_int, x, y, z )
!
!     else if ( dist3 == 0.0E+00 ) then
!
!       call plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, &
!         dist1, dist2, num_int, x, y, z )
!
!     end if
!
!     return
!
!   end if
! !
! !  All nodal distances are nonzero, and there is at least one
! !  positive and one negative.
! !
!   if ( dist1 * dist2 < 0.0E+00 .and. dist1 * dist3 < 0.0E+00 ) then
!
!     call plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, &
!       dist1, dist2, num_int, x, y, z )
!
!     call plane_triangle_int_add_3d ( x1, y1, z1, x3, y3, z3, &
!       dist1, dist3, num_int, x, y, z )
!
!   else if ( dist2 * dist1 < 0.0E+00 .and. dist2 * dist3 < 0.0E+00 ) then
!
!     call plane_triangle_int_add_3d ( x2, y2, z2, x1, y1, z1, &
!       dist2, dist1, num_int, x, y, z )
!
!     call plane_triangle_int_add_3d ( x2, y2, z2, x3, y3, z3, &
!       dist2, dist3, num_int, x, y, z )
!
!   else if ( dist3 * dist1 < 0.0E+00 .and. dist3 * dist2 < 0.0E+00 ) then
!
!     call plane_triangle_int_add_3d ( x3, y3, z3, x1, y1, z1, &
!       dist3, dist1, num_int, x, y, z )
!
!     call plane_triangle_int_add_3d ( x3, y3, z3, x2, y2, z2, &
!       dist3, dist2, num_int, x, y, z )
!
!   end if
!
!   return
! end subroutine
! subroutine plane_norm2exp_3d ( xp, yp, zp, xn, yn, zn, x1, y1, z1, &
!   x2, y2, z2, x3, y3, z3 )
! !
! !*******************************************************************************
! !
! !! PLANE_NORM2EXP_3D converts a normal plane to explicit form in 3D.
! !
! !
! !  Definition:
! !
! !    The normal form of a plane in 3D is
! !
! !      (Xp,Yp,Zp), a point on the plane, and
! !      (Xn,Yn,Zn), the unit normal to the plane.
! !
! !    The explicit form of a plane in 3D is
! !
! !      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
! !
! !  Modified:
! !
! !    27 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real XP, YP, ZP, a point on the plane.  (Actually, we never
! !    need to know these values to do the calculation!)
! !
! !    Input, real XN, YN, ZN, a normal vector N to the plane.  The
! !    vector must not have zero length, but it is not necessary for N
! !    to have unit length.
! !
! !    Output, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
! !    three points that lie on the plane.
! !
!   real x1
!   real x2
!   real x3
!   real xn
!   real xp
!   real xq
!   real xr
!   real y1
!   real y2
!   real y3
!   real yn
!   real yp
!   real yq
!   real yr
!   real z1
!   real z2
!   real z3
!   real zn
!   real zp
!   real zq
!   real zr
!
!   call plane_norm_basis_3d ( xp, yp, zp, xn, yn, zn, xq, yq, zq, xr, yr, zr )
!
!   x1 = xp
!   y1 = yp
!   z1 = zp
!
!   x2 = xp + xq
!   y2 = yp + yq
!   z2 = zp + zq
!
!   x3 = xp + xr
!   y3 = yp + yr
!   z3 = zp + zr
!
!   return
! end subroutine
! subroutine plane_norm2imp_3d ( xp, yp, zp, xn, yn, zn, a, b, c, d )
! !
! !*******************************************************************************
! !
! !! PLANE_NORM2IMP_3D converts a normal form plane to implicit form in 3D.
! !
! !
! !  Definition:
! !
! !    The normal form of a plane in 3D is
! !
! !      (Xp,Yp,Zp), a point on the plane, and
! !      (Xn,Yn,Zn), the unit normal to the plane.
! !
! !    The implicit form of a plane in 3D is
! !
! !      A * X + B * Y + C * Z + D = 0.
! !
! !  Modified:
! !
! !    02 June 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real XP, YP, ZP, a point on the plane.
! !
! !    Input, real XN, YN, ZN, the unit normal vector to the plane.
! !
! !    Output, real A, B, C, D, parameters that define the implicit plane.
! !
!   real a
!   real b
!   real c
!   real d
!   real xn
!   real xp
!   real yn
!   real yp
!   real zn
!   real zp
!
!   a = xn
!   b = yn
!   c = zn
!   d = - a * xp - b * yp - c * zp
!
!   return
! end subroutine
! subroutine plane_triangle_int_add_3d ( x1, y1, z1, x2, y2, z2, dist1, dist2, &
!   num_int, x, y, z )
! !
! !*******************************************************************************
! !
! !! PLANE_TRIANGLE_INT_ADD_3D is a utility for plane/triangle intersections.
! !
! !
! !  Discussion:
! !
! !    This routine is called to consider the value of the signed distance
! !    from a plane of two nodes of a triangle.  If the two values
! !    have opposite signs, then there is a point of intersection between
! !    them.  The routine computes this point and adds it to the list.
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of two vertices of a triangle.
! !
! !    Input, real DIST1, DIST2, the signed distances of the two vertices
! !    from a plane.
! !
! !    Input/output, integer NUM_INT, the number of intersection points.
! !
! !    Input/output, real X(NUM_INT), Y(NUM_INT), Z(NUM_INT), the coordinates
! !    of the intersection points.
! !
!   real alpha
!   real dist1
!   real dist2
!   integer num_int
!   real x(*)
!   real x1
!   real x2
!   real y(*)
!   real y1
!   real y2
!   real z(*)
!   real z1
!   real z2
! !
!   if ( dist1 == 0.0E+00 ) then
!     num_int = num_int + 1
!     x(num_int) = x1
!     y(num_int) = y1
!     z(num_int) = z1
!   else if ( dist2 == 0.0E+00 ) then
!     num_int = num_int + 1
!     x(num_int) = x2
!     y(num_int) = y2
!     z(num_int) = z2
!   else if ( dist1 * dist2 < 0.0E+00 ) then
!     alpha = dist2 / ( dist2 - dist1 )
!     num_int = num_int + 1
!     x(num_int) = alpha * x1 + ( 1.0E+00 - alpha ) * x2
!     y(num_int) = alpha * y1 + ( 1.0E+00 - alpha ) * y2
!     z(num_int) = alpha * z1 + ( 1.0E+00 - alpha ) * z2
!   end if
!
!   return
! end subroutine
! function points_avoid_point_naive_2d ( n, xy_set, xy_test )
! !
! !*******************************************************************************
! !
! !! POINTS_AVOID_POINT_NAIVE_2D determines if a point is "far enough" from a set of points in 2D.
! !
! !
! !  Discussion:
! !
! !    The routine discards points that are too close to other points.
! !    The method used to check this is quadratic in the number of points,
! !    and may take an inordinate amount of time if there are a large
! !    number of points.  But in that case, what do you want?  If you want
! !    lots of points, you don't want to delete any because it won't matter.
! !
! !    The test point is "far enough" from an accepted point if
! !    the Euclidean distance is at least 100 times EPSILON.
! !
! !  Modified:
! !
! !    24 February 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of accepted points.
! !
! !    Input, real XY_SET(2,N), the accepted points.
! !
! !    Input, real XY_TEST(2), a point to be tested.
! !
! !    Output, logical POINTS_AVOID_POINT_NAIVE_2D, is TRUE if XY_TEST is
! !    "far enough" from all the accepted points.
! !
!   integer n
!   integer, parameter :: ndim = 2
! !
!   integer j
!   logical points_avoid_point_naive_2d
!   real tol
!   real xy_set(ndim,n)
!   real xy_test(ndim)
! !
!   tol = 100.0E+00 * epsilon ( tol )
!
!   points_avoid_point_naive_2d = .true.
!
!   do j = 1, n
!
!     if ( sqrt ( sum ( ( xy_set(1:ndim,j) - xy_test(1:ndim) )**2 ) ) < tol ) then
!       points_avoid_point_naive_2d = .false.
!       return
!     end if
!
!   end do
!
!   return
! end function
! subroutine points_bisect_line_imp_2d ( x1, y1, x2, y2, a, b, c )
! !
! !*******************************************************************************
! !
! !! POINTS_BISECT_LINE_IMP_2D finds the implicit line bisecting the line between two points in 2D.
! !
! !
! !  Formula:
! !
! !    The implicit form of a line in 2D is:
! !
! !      A * X + B * Y + C = 0
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    26 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, the coordinates of two points.
! !
! !    Output, real A, B, C, the parameters of the implicit line
! !    equidistant from both points.
! !
!   real a
!   real b
!   real c
!   real x1
!   real x2
!   real y1
!   real y2
! !
!   a = x1 - x2
!   b = y1 - y2
!   c = - 0.5E+00 * ( ( x1 * x1 + y1 * y1 ) - ( x2 * x2 + y2 * y2 ) )
!
!   return
! end subroutine
! subroutine points_bisect_line_par_2d ( x1, y1, x2, y2, f, g, x, y )
! !
! !*******************************************************************************
! !
! !! POINTS_BISECT_LINE_PAR_2D finds the parametric line bisecting the line between two points in 2D.
! !
! !
! !  Formula:
! !
! !    The parametric form of a line in 2D is:
! !
! !      X = X0 + F * T
! !      Y = Y0 + G * T
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    26 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, the coordinates of two points.
! !
! !    Output, real F, G, X, Y, the parameters of the parametric line
! !    equidistant from both points.
! !
!   real f
!   real g
!   real x
!   real x1
!   real x2
!   real y
!   real y1
!   real y2
! !
!   f = 0.5E+00 * ( x1 + x2 )
!   g = 0.5E+00 * ( y1 + y2 )
!   x = - ( y2 - y1 )
!   y = + ( x2 - x1 )
!
!   return
! end subroutine
! subroutine points_centroid_2d ( n, x, y, cent )
! !
! !*******************************************************************************
! !
! !! POINTS_CENTROID_2D computes the discrete centroid of a point set in 2D.
! !
! !
! !  Formula:
! !
! !    Given a discrete set of points S, the discrete centroid z is defined by
! !
! !                           Sum ( x in S ) ( x - z )**2
! !        = min ( y in S ) { Sum ( x in S ) ( x - y )**2
! !
! !    In other words, the discrete centroid is a point in the set whose distance
! !    to the other points is minimized.  The discrete centroid of a point set
! !    need not be unique.  Consider a point set that comprises the
! !    vertices of an equilateral triangle.
! !
! !    This discrete centroid may also be referred to as the K-means cluster.
! !
! !  Modified:
! !
! !    16 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of points.
! !
! !    Input, real X(N), Y(N), the coordinates of the points.
! !
! !    Output, integer CENT, the index of a discrete centroid of the set.
! !
!   integer n
! !
!   integer cent
!   real dist
!   real dist_min
!   integer i
!   integer j
!   real x(n)
!   real y(n)
! !
!   dist_min = 0.0E+00
!   cent = 0
!
!   do i = 1, n
!
!     dist = 0.0E+00
!     do j = 1, n
!       dist = dist + ( x(i) - x(j) )**2 + ( y(i) - y(j) )**2
!     end do
!
!     if ( i == 1 ) then
!       dist_min = dist
!       cent = i
!     else if ( dist < dist_min ) then
!       dist_min = dist
!       cent = i
!     end if
!
!   end do
!
!   return
! end subroutine
! subroutine points_colin_2d ( x1, y1, x2, y2, x3, y3, colin )
! !
! !*******************************************************************************
! !
! !! POINTS_COLIN_2D estimates the colinearity of 3 points in 2D.
! !
! !
! !  Modified:
! !
! !    13 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of the points.
! !
! !    Output, real COLIN, an estimate of colinearity, namely, the ratio
! !    of the area of the triangle spanned by the points to the area
! !    of the equilateral triangle with the same perimeter.
! !    COLIN is 1 if the points are maximally noncolinear, 0 if the
! !    points are exactly colinear, and otherwise is closer to 1 or 0 depending
! !    on whether the points are far or close to colinearity.
! !
!   real area
!   real area2
!   real colin
!   real perim
!   real side
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
! !
!   call triangle_area_2d ( x1, y1, x2, y2, x3, y3, area )
!
!   if ( area == 0.0E+00 ) then
!
!     colin = 0.0E+00
!
!   else
!
!     perim =   enorm0_2d ( x1, y1, x2, y2 ) &
!             + enorm0_2d ( x2, y2, x3, y3 ) &
!             + enorm0_2d ( x3, y3, x1, y1 )
!
!     side = perim / 3.0E+00
!
!     area2 = 0.25E+00 * sqrt ( 3.0E+00 ) * side * side
!
!     colin = area / area2
!
!   end if
!
!   return
! end subroutine
! subroutine points_colin_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, colin )
! !
! !*******************************************************************************
! !
! !! POINTS_COLIN_3D estimates the colinearity of 3 points in 3D.
! !
! !
! !  Modified:
! !
! !    13 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates of
! !    the points.
! !
! !    Output, real COLIN, an estimate of colinearity, namely, the ratio
! !    of the area of the triangle spanned by the points to the area
! !    of the equilateral triangle with the same perimeter.
! !    COLIN is 1 if the points are maximally noncolinear, 0 if the
! !    points are exactly colinear, and otherwise is closer to 1 or 0 depending
! !    on whether the points are far or close to colinearity.
! !
!   real area
!   real area2
!   real colin
!   real perim
!   real side
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
!   real z1
!   real z2
!   real z3
! !
!   call triangle_area_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, area )
!
!   area = abs ( area )
!
!   if ( area == 0.0E+00 ) then
!
!     colin = 0.0E+00
!
!   else
!
!     perim = enorm0_3d ( x1, y1, z1, x2, y2, z2 ) &
!           + enorm0_3d ( x2, y2, z2, x3, y3, z3 ) &
!           + enorm0_3d ( x3, y3, z3, x1, y1, z1 )
!
!     side = perim / 3.0E+00
!
!     area2 = 0.25E+00 * sqrt ( 3.0E+00 ) * side * side
!
!     colin = area / area2
!
!   end if
!
!   return
! end subroutine
! subroutine points_delaunay_naive_2d ( n, x, y, maxtri, ntri, tri )
! !
! !*******************************************************************************
! !
! !! POINTS_DELAUNAY_NAIVE_2D is a naive Delaunay triangulation scheme.
! !
! !
! !  Discussion:
! !
! !    This routine is only suitable as a demonstration code for small
! !    problems.  Its running time is of order N**4.  Much faster
! !    algorithms are available.
! !
! !    Given a set of nodes in the plane, a triangulation is set of
! !    triples of distinct nodes, forming triangles, so that every
! !    point within the convex hull of the set of nodes is either
! !    one of the nodes, or lies on an edge of one or more triangles,
! !    or lies within exactly one triangle.
! !
! !    A Delaunay triangulation is a triangulation with additional
! !    properties.
! !
! !  Reference:
! !
! !    Joseph O'Rourke,
! !    Computational Geometry,
! !    Cambridge University Press,
! !    Second Edition, 1998, page 187.
! !
! !  Modified:
! !
! !    08 November 2000
! !
! !  Parameters:
! !
! !    Input, integer N, the number of nodes.  N must be at least 3.
! !
! !    Input, real X(N), Y(N), the coordinates of the nodes.
! !
! !    Input, integer MAXTRI, the maximum number of triangles.
! !
! !    Output, integer NTRI, the number of triangles in the triangulation.
! !
! !    Output, integer TRI(3,MAXTRI), contains in TRI(1,I), TRI(2,I)
! !    and TRI(3,I) the indices of the nodes making the I-th triangle.
! !
!   integer maxtri
!   integer n
! !
!   logical flag
!   integer i
!   integer j
!   integer k
!   integer m
!   integer ntri
!   integer tri(3,maxtri)
!   real x(n)
!   real xn
!   real y(n)
!   real yn
!   real z(n)
!   real zn
! !
!   ntri = 0
!
!   if ( n < 3 ) then
!     return
!   end if
! !
! !  Compute Z = X*X + Y*Y.
! !
!   z(1:n) = x(1:n)**2 + y(1:n)**2
! !
! !  For each triple (I,J,K):
! !
!   do i = 1, n-2
!     do j = i+1, n
!       do k = i+1, n
!
!         if ( j /= k ) then
!
!           xn = ( y(j) - y(i) ) * ( z(k) - z(i) ) &
!              - ( y(k) - y(i) ) * ( z(j) - z(i) )
!
!           yn = ( x(k) - x(i) ) * ( z(j) - z(i) ) &
!              - ( x(j) - x(i) ) * ( z(k) - z(i) )
!
!           zn = ( x(j) - x(i) ) * ( y(k) - y(i) ) &
!              - ( x(k) - x(i) ) * ( y(j) - y(i) )
!
!           flag = ( zn < 0.0E+00 )
!
!           if ( flag ) then
!             do m = 1, n
!               flag = flag .and. &
!                 ( ( x(m) - x(i) ) * xn + ( y(m) - y(i) ) * yn &
!                   + ( z(m) - z(i) ) * zn <= 0.0E+00 )
!             end do
!           end if
!
!           if ( flag ) then
!             if ( ntri < maxtri ) then
!               ntri = ntri + 1
!               tri(1,ntri) = i
!               tri(2,ntri) = j
!               tri(3,ntri) = k
!             end if
!           end if
!
!         end if
!
!       end do
!     end do
!   end do
!
!   return
! end subroutine
! subroutine points_dist_2d ( x1, y1, x2, y2, dist )
! !
! !*******************************************************************************
! !
! !! POINTS_DIST_2D finds the distance between two points in 2D.
! !
! !
! !  Modified:
! !
! !    17 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, determines the pair of
! !    points (X1,Y1) and (X2,Y2) whose distance apart is be determined.
! !
! !    Output, real DIST, the distance between the points.
! !
!   real dist
!   real x1
!   real x2
!   real y1
!   real y2
! !
!   dist = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 )
!
!   return
! end subroutine
! subroutine points_dist_3d ( x1, y1, z1, x2, y2, z2, dist )
! !
! !*******************************************************************************
! !
! !! POINTS_DIST_3D finds the distance between two points in 3D.
! !
! !
! !  Modified:
! !
! !    27 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, determines the pair of points
! !    (X1,Y1,Z1) and (X2,Y2,Z2) whose distance apart is be determined.
! !
! !    Output, real DIST, the distance between the points.
! !
!   real dist
!   real x1
!   real x2
!   real y1
!   real y2
!   real z1
!   real z2
! !
!   dist = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 + ( z1 - z2 )**2 )
!
!   return
! end subroutine
! subroutine points_dist_nd ( x1, x2, n, dist )
! !
! !*******************************************************************************
! !
! !! POINTS_DIST_ND finds the distance between two points in ND.
! !
! !
! !  Modified:
! !
! !    31 January 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1(N), X2(N), the coordinates of two points.
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Output, real DIST, the distance between the points.
! !
!   integer n
! !
!   real dist
!   real x1(n)
!   real x2(n)
! !
!   dist = sqrt ( sum ( ( x1(1:n) - x2(1:n) )**2 ) )
!
!   return
! end subroutine
! subroutine points_hull_2d ( ival, n, nval, x, y )
! !
! !*******************************************************************************
! !
! !! POINTS_HULL_2D computes the convex hull of a set of points in 2D.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, integer IVAL(N).  Entries 1 through NVAL of IVAL contain
! !    the indices of the nodes that form the convex hull, in order.
! !
! !    Input, integer N, the number of nodes.
! !
! !    Output, integer NVAL, the number of nodes that lie on the convex hull.
! !
! !    Input, real X(N), Y(N), the X and Y coordinates of the nodes.
! !
!   integer n
! !
!   real ai
!   real angmax
!   real di
!   real dr
!   integer i
!   integer iq
!   integer ir
!   integer istart
!   integer ival(n)
!   integer nval
!   real x(n)
!   real xp
!   real xq
!   real xr
!   real y(n)
!   real yp
!   real yq
!   real yr
! !
!   if ( n < 1 ) then
!     nval = 0
!     return
!   end if
! !
! !  If N = 1, the hull is the point.
! !
!   if ( n == 1 ) then
!     nval = 1
!     ival(1) = 1
!     return
!   end if
! !
! !  If N = 2, then the convex hull is either the two distinct points,
! !  or possibly a single (repeated) point.
! !
!   if ( n == 2 ) then
!
!     if ( x(1) /= x(2) .or. y(1) /= y(2) ) then
!       nval = 2
!       ival(1) = 1
!       ival(2) = 2
!     else
!       nval = 1
!       ival(1) = 1
!     end if
!
!     return
!
!   end if
! !
! !  Find the leftmost point, and take the bottom-most in a tie.
! !  Call it "Q".
! !
!   iq = 1
!   do i = 2, n
!     if ( x(i) < x(iq) .or. ( x(i) == x(iq) .and. y(i) < y(iq) ) ) then
!       iq = i
!     end if
!   end do
!
!   xq = x(iq)
!   yq = y(iq)
! !
! !  Remember the starting point.
! !
!   istart = iq
!   nval = 1
!   ival(1) = iq
! !
! !  For the first point, make a dummy previous point, 1 unit south,
! !  and call it "P".
! !
!   xp = xq
!   yp = yq - 1.0E+00
! !
! !  Now, having old point P, and current point Q, find the new point R
! !  so the angle PQR is maximal.
! !
! !  Watch out for the possibility that the two nodes are identical.
! !
!   do
!
!     ir = 0
!     angmax = 0.0E+00
!
!     do i = 1, n
!
!       if ( i /= iq .and. ( x(i) /= xq .or. y(i) /= yq ) ) then
!
!         ai = angle_deg_2d ( xp, yp, xq, yq, x(i), y(i) )
!
!         if ( ir == 0 .or. ai > angmax ) then
!
!           ir = i
!           xr = x(ir)
!           yr = y(ir)
!           angmax = ai
! !
! !  In case of ties, choose the nearer point.
! !
!         else if ( ir /= 0 .and. ai == angmax ) then
!
!           di = enorm0_2d ( xq, yq, x(i), y(i) )
!           dr = enorm0_2d ( xq, yq, xr, yr )
!
!           if ( di < dr ) then
!             ir = i
!             xr = x(ir)
!             yr = y(ir)
!             angmax = ai
!           end if
!
!         end if
!
!       end if
!
!     end do
! !
! !  If we've returned to our starting point, exit.
! !
!     if ( ir == istart ) then
!       exit
!     end if
!
!     nval = nval + 1
!
!     if ( nval > n ) then
!       write ( *, * ) ' '
!       write ( *, * ) 'POINTS_HULL_2D - Fatal error!'
!       write ( *, * ) '  The algorithm failed.'
!       stop
!     end if
! !
! !  Set Q := P, P := R, and repeat.
! !
!     ival(nval) = ir
!     xp = xq
!     yp = yq
!     iq = ir
!     xq = xr
!     yq = yr
!
!   end do
!
!   return
! end subroutine
! subroutine points_nearest_point_naive_nd ( ndim, n, pset, p, i_min, d_min )
! !
! !*******************************************************************************
! !
! !! POINTS_NEAREST_POINT_NAIVE_ND finds the nearest point to a given point in ND.
! !
! !
! !  Discussion:
! !
! !    A naive algorithm is used.  No attempt is made to optimize the
! !    calculation, so there will be N distance calculations done.
! !
! !    For a large dataset, it would be better to group the points into
! !    clusters, so that far fewer distance calculations need to be done.
! !
! !  Modified:
! !
! !    31 January 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer NDIM, the dimension of the points.
! !
! !    Input, integer N, the number of points in the set.
! !
! !    Input, real PSET(NDIM,N), the coordinates of the points in the set.
! !
! !    Input, real P(NDIM), the point to be tested.
! !
! !    Output, integer I_MIN, the index of the nearest point in PSET to P.
! !
! !    Output, real D_MIN, the distance between P and PSET(*,I_MIN).
! !
!   integer n
!   integer ndim
! !
!   real d
!   real d_min
!   integer i
!   integer i_min
!   real p(ndim)
!   real pset(ndim,n)
! !
!   d_min = huge ( d_min )
!   i_min = 0
!
!   do i = 1, n
!
!     d = sum ( ( p(1:ndim) - pset(1:ndim,i) )**2 )
!
!     if ( d < d_min ) then
!       d_min = d
!       i_min = i
!     end if
!
!   end do
! !
! !  We save a little work by waiting til the end to take the square root.
! !
!   d_min = sqrt ( d_min )
!
!   return
! end subroutine
! subroutine polygon_1_2d ( result, n, x, y )
! !
! !*******************************************************************************
! !
! !! POLYGON_1_2D integrates the function 1 over a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    INTEGRAL = 0.5 * SUM ( I = 1 to N )
! !      ( X(I) + X(I-1) ) * ( Y(I) - Y(I-1) )
! !
! !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
! !
! !  Note:
! !
! !    The integral of 1 over a polygon is the area of the polygon.
! !
! !  Reference:
! !
! !    S F Bockman,
! !    Generalizing the Formula for Areas of Polygons to Moments,
! !    American Mathematical Society Monthly,
! !    1989, pages 131-132.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real RESULT, the value of the integral.
! !
! !    Input, integer N, the number of vertices of the polygon.
! !    N should be at least 3 for a nonzero result.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices of the polygon.
! !    These vertices should be given in counter-clockwise order.
! !
!   integer n
! !
!   integer i
!   integer im1
!   real result
!   real x(n)
!   real y(n)
! !
!   result = 0.0E+00
!
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_1_2D - Warning!'
!     write ( *, * ) '  The number of vertices must be at least 3.'
!     write ( *, * ) '  The input value of N = ', n
!     return
!   end if
!
!   do i = 1, n
!
!     if ( i == 1 ) then
!       im1 = n
!     else
!       im1 = i - 1
!     end if
!
!     result = result + 0.5E+00 * ( x(i) + x(im1) ) * ( y(i) - y(im1) )
!
!   end do
!
!   return
! end subroutine
! subroutine polygon_x_2d ( result, n, x, y )
! !
! !*******************************************************************************
! !
! !! POLYGON_X_2D integrates the function X over a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    INTEGRAL = (1/6) * SUM ( I = 1 to N )
! !      ( X(I)**2 + X(I) * X(I-1) + X(I-1)**2 ) * ( Y(I) - Y(I-1) )
! !
! !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    S F Bockman,
! !    Generalizing the Formula for Areas of Polygons to Moments,
! !    American Mathematical Society Monthly,
! !    1989, pages 131-132.
! !
! !  Parameters:
! !
! !    Output, real RESULT, the value of the integral.
! !
! !    Input, integer N, the number of vertices of the polygon.
! !    N should be at least 3 for a nonzero result.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices of the polygon.
! !    These vertices should be given in counter-clockwise order.
! !
!   integer n
! !
!   integer i
!   integer im1
!   real result
!   real x(n)
!   real y(n)
! !
!   result = 0.0E+00
!
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_X_2D - Warning!'
!     write ( *, * ) '  The number of vertices must be at least 3.'
!     write ( *, * ) '  The input value of N = ', n
!     return
!   end if
!
!   do i = 1, n
!
!     if ( i == 1 ) then
!       im1 = n
!     else
!       im1 = i - 1
!     end if
!
!     result = result + ( x(i)**2 + x(i) * x(im1) + x(im1)**2 ) * ( y(i) - y(im1) )
!
!   end do
!
!   result = result / 6.0E+00
!
!   return
! end subroutine
! subroutine polygon_y_2d ( result, n, x, y )
! !
! !*******************************************************************************
! !
! !! POLYGON_Y_2D integrates the function Y over a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    INTEGRAL = (1/6) * SUM ( I = 1 to N )
! !      - ( Y(I)**2 + Y(I) * Y(I-1) + Y(I-1)**2 ) * ( X(I) - X(I-1) )
! !
! !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
! !
! !  Reference:
! !
! !    S F Bockman,
! !    Generalizing the Formula for Areas of Polygons to Moments,
! !    American Mathematical Society Monthly,
! !    1989, pages 131-132.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real RESULT, the value of the integral.
! !
! !    Input, integer N, the number of vertices of the polygon.
! !    N should be at least 3 for a nonzero result.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices
! !    of the polygon.  These vertices should be given in
! !    counter-clockwise order.
! !
!   integer n
! !
!   integer i
!   integer im1
!   real result
!   real x(n)
!   real y(n)
! !
!   result = 0.0E+00
!
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_Y_2D - Warning!'
!     write ( *, * ) '  The number of vertices must be at least 3.'
!     write ( *, * ) '  The input value of N = ', n
!     return
!   end if
!
!   do i = 1, n
!
!     if ( i == 1 ) then
!       im1 = n
!     else
!       im1 = i - 1
!     end if
!
!     result = result - ( y(i)**2 + y(i) * y(im1) + y(im1)**2 ) * ( x(i) - x(im1) )
!
!   end do
!
!   result = result / 6.0E+00
!
!   return
! end subroutine
! subroutine polygon_xx_2d ( result, n, x, y )
! !
! !*******************************************************************************
! !
! !! POLYGON_XX_2D integrates the function X*X over a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    INTEGRAL = (1/12) * SUM ( I = 1 to N )
! !      ( X(I)**3 + X(I)**2 * X(I-1) + X(I) * X(I-1)**2 + X(I-1)**3 )
! !      * ( Y(I) - Y(I-1) )
! !
! !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
! !
! !  Reference:
! !
! !    S F Bockman,
! !    Generalizing the Formula for Areas of Polygons to Moments,
! !    American Mathematical Society Monthly,
! !    1989, pages 131-132.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real RESULT, the value of the integral.
! !
! !    Input, integer N, the number of vertices of the polygon.
! !    N should be at least 3 for a nonzero result.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices
! !    of the polygon.  These vertices should be given in
! !    counter-clockwise order.
! !
!   integer n
! !
!   integer i
!   integer im1
!   real result
!   real x(n)
!   real y(n)
! !
!   result = 0.0E+00
!
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_XX_2D - Warning!'
!     write ( *, * ) '  The number of vertices must be at least 3.'
!     write ( *, * ) '  The input value of N = ', n
!     return
!   end if
!
!   do i = 1, n
!
!     if ( i == 1 ) then
!       im1 = n
!     else
!       im1 = i - 1
!     end if
!
!     result = result + ( x(i)**3 + x(i)**2 * x(im1) &
!       + x(i) * x(im1)**2 + x(im1)**3 ) * ( y(i) - y(im1) )
!
!   end do
!
!   result = result / 12.0E+00
!
!   return
! end subroutine
! subroutine polygon_xy_2d ( result, n, x, y )
! !
! !*******************************************************************************
! !
! !! POLYGON_XY_2D integrates the function X*Y over a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    INTEGRAL = (1/24) * SUM ( I = 1 to N )
! !      ( Y(I)   * ( 3 * X(I)**2 + 2 * X(I) * X(I-1) +     X(I-1)**2 )
! !      + Y(I-1) * (     X(I)**2 + 2 * X(I) * X(I-1) + 3 * X(I-1)**2 ) )
! !      * ( Y(I) - Y(I-1) )
! !
! !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
! !
! !  Reference:
! !
! !    S F Bockman,
! !    Generalizing the Formula for Areas of Polygons to Moments,
! !    American Mathematical Society Monthly,
! !    1989, pages 131-132.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real RESULT, the value of the integral.
! !
! !    Input, integer N, the number of vertices of the polygon.
! !    N should be at least 3 for a nonzero result.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices
! !    of the polygon.  These vertices should be given in
! !    counter-clockwise order.
! !
!   integer n
! !
!   integer i
!   integer im1
!   real result
!   real x(n)
!   real y(n)
! !
!   result = 0.0E+00
!
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_XY_2D - Warning!'
!     write ( *, * ) '  The number of vertices must be at least 3.'
!     write ( *, * ) '  The input value of N = ', n
!     return
!   end if
!
!   do i = 1, n
!
!     if ( i == 1 ) then
!       im1 = n
!     else
!       im1 = i - 1
!     end if
!
!     result = result + ( &
!       y(i) * ( 3.0E+00 * x(i)**2 + 2.0E+00 * x(i) * x(im1) + x(im1)**2 ) &
!       + y(im1) * ( x(i)**2 + 2.0E+00 * x(i) * x(im1) + 3.0E+00 * x(im1)**2 ) &
!       ) * ( y(i) - y(im1) )
!
!   end do
!
!   result = result / 24.0E+00
!
!   return
! end subroutine
! subroutine polygon_yy_2d ( result, n, x, y )
! !
! !*******************************************************************************
! !
! !! POLYGON_YY_2D integrates the function Y*Y over a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    INTEGRAL = (1/12) * SUM ( I = 1 to N )
! !      - ( Y(I)**3 + Y(I)**2 * Y(I-1) + Y(I) * Y(I-1)**2 + Y(I-1)**3 )
! !      * ( X(I) - X(I-1) )
! !
! !    where X(0) and Y(0) should be replaced by X(N) and Y(N).
! !
! !  Reference:
! !
! !    S F Bockman,
! !    Generalizing the Formula for Areas of Polygons to Moments,
! !    American Mathematical Society Monthly,
! !    1989, pages 131-132.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real RESULT, the value of the integral.
! !
! !    Input, integer N, the number of vertices of the polygon.
! !    N should be at least 3 for a nonzero result.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices
! !    of the polygon.  These vertices should be given in
! !    counter-clockwise order.
! !
!   integer n
! !
!   integer i
!   integer im1
!   real result
!   real x(n)
!   real y(n)
! !
!   result = 0.0E+00
!
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_YY_2D - Warning!'
!     write ( *, * ) '  The number of polygonal vertices must be '
!     write ( *, * ) '  at least 3, but the input polygon has N = ', n
!     return
!   end if
!
!   do i = 1, n
!
!     if ( i == 1 ) then
!       im1 = n
!     else
!       im1 = i - 1
!     end if
!
!     result = result - ( y(i)**3 + y(i)**2 * y(im1) + y(i) * y(im1)**2 &
!       + y(im1)**3 ) * ( x(i) - x(im1) )
!
!   end do
!
!   result = result / 12.0E+00
!
!   return
! end subroutine
! subroutine polygon_area_2d ( n, x, y, area )
! !
! !*******************************************************************************
! !
! !! POLYGON_AREA_2D computes the area of a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    AREA = ABS ( 0.5 * SUM ( I = 1 to N ) X(I) * ( Y(I+1) - Y(I-1) ) )
! !    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of vertices of the polygon.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices.
! !
! !    Output, real AREA, the absolute area of the polygon.
! !
!   integer n
! !
!   real area
!   integer i
!   integer im1
!   integer ip1
!   real x(n)
!   real y(n)
! !
!   area = 0.0E+00
!
!   do i = 1, n
!
!     if ( i > 1 ) then
!       im1 = i - 1
!     else
!       im1 = n
!     end if
!
!     if ( i < n ) then
!       ip1 = i + 1
!     else
!       ip1 = 1
!     end if
!
!     area = area + x(i) * ( y(ip1) - y(im1) )
!
!   end do
!
!   area = 0.5E+00 * abs ( area )
!
!   return
! end subroutine
! subroutine polygon_area_2_2d ( n, x, y, area )
! !
! !*******************************************************************************
! !
! !! POLYGON_AREA_2_2D computes the area of a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    The area is the sum of the areas of the triangles formed by
! !    node N with consecutive pairs of nodes.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of vertices of the polygon.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices.
! !
! !    Output, real AREA, the absolute area of the polygon.
! !
!   integer n
! !
!   real area
!   real areat
!   integer i
!   real x(n)
!   real y(n)
! !
!   area = 0.0E+00
!
!   do i = 1, n - 2
!
!     call triangle_area_2d ( x(i), y(i), x(i+1), y(i+1), x(n), y(n), areat )
!
!     area = area + areat
!
!   end do
!
!   return
! end subroutine
! subroutine polygon_area_3d ( n, x, y, z, area, normal )
! !
! !*******************************************************************************
! !
! !! POLYGON_AREA_3D computes the area of a polygon in 3D.
! !
! !
! !  Restriction:
! !
! !    The computation is not valid unless the vertices really do lie
! !    in a plane, so that the polygon that is defined is "flat".
! !    The polygon does not have to be "regular", that is, neither its
! !    sides nor its angles need to be equal.
! !
! !  Reference:
! !
! !    Allen Van Gelder,
! !    Efficient Computation of Polygon Area and Polyhedron Volume,
! !    Graphics Gems V, edited by Alan Paeth,
! !    AP Professional, 1995.
! !
! !  Modified:
! !
! !    30 September 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of vertices.
! !
! !    Input, real X(N), Y(N), Z(N), the coordinates of the vertices.
! !    The vertices should be listed in neighboring order.
! !
! !    Output, real AREA, the area of the polygon.
! !
! !    Output, real NORMAL(3), the unit normal vector to the polygon.
! !
!   integer n
! !
!   real area
!   integer i
!   integer ip1
!   real normal(3)
!   real x(n)
!   real x1
!   real x2
!   real x3
!   real y(n)
!   real y1
!   real y2
!   real y3
!   real z(n)
!   real z1
!   real z2
!   real z3
! !
!   normal(1) = 0.0E+00
!   normal(2) = 0.0E+00
!   normal(3) = 0.0E+00
!
!   do i = 1, n
!
!     x1 = x(i)
!     y1 = y(i)
!     z1 = z(i)
!
!     if ( i < n ) then
!       ip1 = i + 1
!     else
!       ip1 = 1
!     end if
!
!     x2 = x(ip1)
!     y2 = y(ip1)
!     z2 = z(ip1)
!
!     call cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!     normal(1) = normal(1) + x3
!     normal(2) = normal(2) + y3
!     normal(3) = normal(3) + z3
!
!   end do
!
!   area = sqrt ( normal(1)**2 + normal(2)**2 + normal(3)**2 )
!
!   if ( area /= 0.0E+00 ) then
!     normal(1) = normal(1) / area
!     normal(2) = normal(2) / area
!     normal(3) = normal(3) / area
!   else
!     normal(1) = 1.0E+00
!     normal(2) = 0.0E+00
!     normal(3) = 0.0E+00
!   end if
!
!   area = 0.5E+00 * area
!
!   return
! end subroutine
! subroutine polygon_area_2_3d ( n, x, y, z, area )
! !
! !*******************************************************************************
! !
! !! POLYGON_AREA_2_3D computes the area of a polygon in 3D.
! !
! !
! !  Formula:
! !
! !    The area is the sum of the areas of the triangles formed by
! !    node N with consecutive pairs of nodes.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of vertices of the polygon.
! !
! !    Input, real X(N), Y(N), Z(N), the coordinates of the vertices.
! !
! !    Output, real AREA, the absolute area of the polygon.
! !
!   integer n
! !
!   real area
!   real areat
!   integer i
!   real x(n)
!   real y(n)
!   real z(n)
! !
!   area = 0.0E+00
!
!   do i = 1, n - 2
!
!     call triangle_area_3d ( x(i), y(i), z(i), x(i+1), y(i+1), z(i+1), &
!       x(n), y(n), z(n), areat )
!
!     area = area + areat
!
!   end do
!
!   return
! end subroutine
! subroutine polygon_centroid_2d ( n, x, y, cx, cy )
! !
! !*******************************************************************************
! !
! !! POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
! !
! !
! !  Formula:
! !
! !    Denoting the centroid coordinates by (CX,CY), then
! !
! !      CX = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
! !      CY = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
! !
! !    Green's theorem states that
! !
! !      Integral ( Polygon boundary ) ( M dx + N dy ) =
! !      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
! !
! !    Using M = 0 and N = x**2/2, we get:
! !
! !      CX = 0.5 * Integral ( Polygon boundary ) x**2 dy,
! !
! !    which becomes
! !
! !      CX = 1/6 SUM ( I = 1 to N )
! !        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
! !
! !    where, when I = N, the index "I+1" is replaced by 1.
! !
! !    A similar calculation gives us a formula for CY.
! !
! !  Reference:
! !
! !    Gerard Bashein and Paul Detmer,
! !    Centroid of a Polygon,
! !    Graphics Gems IV, edited by Paul Heckbert,
! !    AP Professional, 1994.
! !
! !  Modified:
! !
! !    22 September 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of sides of the polygonal shape.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices of the shape.
! !
! !    Output, real CX, CY, the coordinates of the centroid of the shape.
! !
!   integer n
! !
!   real area
!   real cx
!   real cy
!   integer i
!   integer ip1
!   real temp
!   real x(n)
!   real y(n)
! !
!   area = 0.0E+00
!   cx = 0.0E+00
!   cy = 0.0E+00
!
!   do i = 1, n
!
!     if ( i < n ) then
!       ip1 = i + 1
!     else
!       ip1 = 1
!     end if
!
!     temp = ( x(i) * y(ip1) - x(ip1) * y(i) )
!
!     area = area + temp
!     cx = cx + ( x(ip1) + x(i) ) * temp
!     cy = cy + ( y(ip1) + y(i) ) * temp
!
!   end do
!
!   area = area / 2.0E+00
!
!   cx = cx / ( 6.0E+00 * area )
!   cy = cy / ( 6.0E+00 * area )
!
!   return
! end subroutine
! subroutine polygon_centroid_2_2d ( n, x, y, cx, cy )
! !
! !*******************************************************************************
! !
! !! POLYGON_CENTROID_2_2D computes the centroid of a polygon in 2D.
! !
! !
! !  Method:
! !
! !    The centroid is the area-weighted sum of the centroids of
! !    disjoint triangles that make up the polygon.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of vertices of the polygon.
! !
! !    Input, real X(N), Y(N), the coordinates of the vertices.
! !
! !    Output, real CX, CY, the coordinates of the centroid.
! !
!   integer n
! !
!   real area
!   real areat
!   real cx
!   real cy
!   integer i
!   real x(n)
!   real y(n)
! !
!   area = 0.0E+00
!   cx = 0.0E+00
!   cy = 0.0E+00
!
!   do i = 1, n - 2
!
!     call triangle_area_2d ( x(i), y(i), x(i+1), y(i+1), x(n), y(n), areat )
!
!     area = area + areat
!     cx = cx + areat * ( x(i) + x(i+1) + x(n) ) / 3.0E+00
!     cy = cy + areat * ( y(i) + y(i+1) + y(n) ) / 3.0E+00
!
!   end do
!
!   cx = cx / area
!   cy = cy / area
!
!   return
! end subroutine
! subroutine polygon_centroid_3d ( n, x, y, z, cx, cy, cz )
! !
! !*******************************************************************************
! !
! !! POLYGON_CENTROID_3D computes the centroid of a polygon in 3D.
! !
! !
! !  Method:
! !
! !    The centroid is the area-weighted sum of the centroids of
! !    disjoint triangles that make up the polygon.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of vertices of the polygon.
! !
! !    Input, real X(N), Y(N), Z(N), the coordinates of the vertices.
! !
! !    Output, real CX, CY, CZ, the coordinates of the centroid.
! !
!   integer n
! !
!   real area
!   real areat
!   real cx
!   real cy
!   real cz
!   integer i
!   real x(n)
!   real y(n)
!   real z(n)
! !
!   area = 0.0E+00
!   cx = 0.0E+00
!   cy = 0.0E+00
!   cz = 0.0E+00
!
!   do i = 1, n - 2
!
!     call triangle_area_3d ( x(i), y(i), z(i), x(i+1), &
!       y(i+1), z(i+1), x(n), y(n), z(n), areat )
!
!     area = area + areat
!     cx = cx + areat * ( x(i) + x(i+1) + x(n) ) / 3.0E+00
!     cy = cy + areat * ( y(i) + y(i+1) + y(n) ) / 3.0E+00
!     cz = cz + areat * ( z(i) + z(i+1) + z(n) ) / 3.0E+00
!
!   end do
!
!   cx = cx / area
!   cy = cy / area
!   cz = cz / area
!
!   return
! end subroutine
!
! subroutine polygon_contains_point_2d ( n, xn, xval, yn, yval, inside )
! !
! !*******************************************************************************
! !
! !! POLYGON_CONTAINS_POINT_2D finds if a point is inside a simple polygon in 2D.
! !
! !
! !  Discussion:
! !
! !    A simple polygon is one whose boundary never crosses itself.
! !
! !  Reference:
! !
! !    ACM Algorithm 112.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of nodes or vertices in the polygon.
! !    N must be at least 3.
! !
! !    Input, real XN(N), the X coordinates of the vertices.
! !
! !    Input, real XVAL, the X coordinate of the point to be tested.
! !
! !    Input, real YN(N), the Y coordinates of the vertices.
! !
! !    Input, real YVAL, the Y coordinate of the point to be tested.
! !
! !    Output, logical INSIDE, is .TRUE. if ( X,Y) is inside
! !    the polygon and .FALSE. otherwise.
! !
!   integer n
! !
!   integer i
!   logical inside
!   real x1
!   real x2
!   real xn(n)
!   real xval
!   real y1
!   real y2
!   real yn(n)
!   real yval
! !
!   inside = .false.
!
!   do i = 1, n
!
!     x1 = xn(i)
!     y1 = yn(i)
!
!     if ( i < n ) then
!       x2 = xn(i+1)
!       y2 = yn(i+1)
!     else
!       x2 = xn(1)
!       y2 = yn(1)
!     end if
!
!     if ( ( y1 < yval .eqv. yval <= y2 ) .and. &
!       ( xval - x1 ) * abs ( y2 - y1 ) < ( x2 - x1 ) * ( yval - y1 ) ) then
!
!       inside = .not. inside
!
!     end if
!
!   end do
!
!   return
! end subroutine
!
! subroutine polygon_contains_point_2_2d ( n, xn, xval, yn, yval, inside )
! !
! !*******************************************************************************
! !
! !! POLYGON_CONTAINS_POINT_2_2D finds if a point is inside a convex polygon in 2D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of nodes or vertices in the polygon.
! !    N must be at least 3.
! !
! !    Input, real XN(N), the X coordinates of the vertices.
! !
! !    Input, real XVAL, the X coordinate of the point to be tested.
! !
! !    Input, real YN(N), the Y coordinates of the vertices.
! !
! !    Input, real YVAL, the Y coordinate of the point to be tested.
! !
! !    Output, logical INSIDE, is .TRUE. if ( X,Y) is inside
! !    the polygon or on its boundary, and .FALSE. otherwise.
! !
!   integer n
! !
!   integer i
!   logical inside
!   real x1
!   real x2
!   real x3
!   real xn(n)
!   real xval
!   real y1
!   real y2
!   real y3
!   real yn(n)
!   real yval
! !
!   inside = .false.
! !
! !  A point is inside a convex polygon if and only if it is inside
! !  one of the triangles formed by X(1),Y(1) and any two consecutive
! !  points on the polygon's circumference.
! !
!   x1 = xn(1)
!   y1 = yn(1)
!
!   do i = 2, n-1
!
!     x2 = xn(i)
!     y2 = yn(i)
!     x3 = xn(i+1)
!     y3 = yn(i+1)
!
!     call triangle_contains_point_2d ( x1, y1, x2, y2, x3, y3, xval, yval, &
!       inside )
!
!     if ( inside ) then
!       return
!     end if
!
!   end do
!
!   return
! end subroutine
! function polygon_convex ( n, x, y )
! !
! !*******************************************************************************
! !
! !! POLYGON_CONVEX determines whether a polygon is convex in 2D.
! !
! !
! !  Discussion:
! !
! !    If the polygon has less than 3 distinct vertices, it is
! !    classified as convex degenerate.
! !
! !    If the polygon "goes around" more than once, it is classified
! !    as NOT convex.
! !
! !  Reference:
! !
! !    Peter Schorn and Frederick Fisher,
! !    Testing the Convexity of a Polygon,
! !    Graphics Gems, 1994.
! !
! !  Modified:
! !
! !    02 May 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters
! !
! !    Input, integer N, the number of vertices.
! !
! !    Input/output, X(N), Y(N), the coordinates of the vertices of the
! !    polygon.  On output, duplicate consecutive points have been deleted,
! !    and the vertices have been reordered so that the lexicographically
! !    least point comes first.
! !
! !    Output, integer POLYGON_CONVEX:
! !    -1, the polygon is not convex;
! !     0, the polygon has less than 3 vertices; it is "degenerately" convex;
! !     1, the polygon is convex and counterclockwise;
! !     2, the polygon is convex and clockwise.
! !
!   integer, parameter :: NOT_CONVEX = -1
!   integer, parameter :: DEGENERATE_CONVEX = 0
!   integer, parameter :: CONVEX_CCW = 1
!   integer, parameter :: CONVEX_CW = 2
!   real, parameter :: PI = 3.14159265358979323846264338327950288419716939937510E+00
!   real, parameter :: RAD_TO_DEG = 180.0E+00 / PI
!   real, parameter :: TOL = 1.0E+00
! !
!   integer n
! !
!   real angle
!   real cross
!   real dot
!   real exterior_total
!   integer i
!   integer ip1
!   integer ip2
!   integer polygon_convex
!   real sense
!   real x(n)
!   real y(n)
! !
!   exterior_total = 0.0E+00
! !
! !  If there are not at least 3 distinct vertices, we are done.
! !
!   if ( n < 3 ) then
!     polygon_convex = DEGENERATE_CONVEX
!     return
!   end if
!
!   sense = 0.0E+00
! !
! !  Consider each polygonal vertex I.
! !
!   do i = 1, n
!
!     ip1 = i + 1
!     if ( ip1 > n ) then
!       ip1 = ip1 - n
!     end if
!
!     ip2 = i + 2
!     if ( ip2 > n ) then
!       ip2 = ip2 - n
!     end if
!
!     dot =   ( x(ip2) - x(ip1) ) * ( x(i) - x(ip1) ) &
!           + ( y(ip2) - y(ip1) ) * ( y(i) - y(ip1) )
!
!     cross =   ( x(ip2) - x(ip1) ) * ( y(i) - y(ip1) ) &
!             - ( x(i) - x(ip1) ) * ( y(ip2) - y(ip1) )
!
!     angle = atan2 ( cross, dot )
! !
! !  See if the turn defined by this vertex is our first indication of
! !  the "sense" of the polygon, or if it disagrees with the previously
! !  defined sense.
! !
!     if ( sense == 0.0E+00 ) then
!
!       if ( angle < 0.0E+00 ) then
!         sense = -1.0E+00
!       else if ( angle > 0.0E+00 ) then
!         sense = +1.0E+00
!       end if
!
!    else if ( sense == 1.0E+00 ) then
!
!       if ( angle < 0.0E+00 ) then
!         polygon_convex = NOT_CONVEX
!         return
!       end if
!
!     else if ( sense == -1.0E+00 ) then
!
!       if ( angle > 0.0E+00 ) then
!         polygon_convex = NOT_CONVEX
!         return
!       end if
!
!     end if
! !
! !  If the exterior total is greater than 360, then the polygon is
! !  going around again.
! !
!     angle = atan2 ( -cross, -dot )
!
!     exterior_total = exterior_total + angle
!
!     if ( abs ( exterior_total ) * RAD_TO_DEG > 360.0E+00 + TOL ) then
!       polygon_convex = NOT_CONVEX
!       return
!     end if
!
!   end do
!
!   if ( sense == +1.0E+00 ) then
!     polygon_convex = CONVEX_CCW
!   else if ( sense == -1.0E+00 ) then
!     polygon_convex = CONVEX_CW
!   end if
!
!   return
! end function
! subroutine polygon_inrad_data_2d ( area, n, radin, radout, side )
! !
! !*******************************************************************************
! !
! !! POLYGON_INRAD_DATA_2D determines polygonal data from its inner radius in 2D.
! !
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real AREA, the area of the regular polygon.
! !
! !    Input, integer N, the number of sides of the polygon.
! !    N must be at least 3.
! !
! !    Input, real RADIN, the inner radius of the polygon, that is,
! !    the radius of the largest circle that can be inscribed within
! !    the polygon.
! !
! !    Output, real RADOUT, the outer radius of the polygon, that is,
! !    the radius of the smallest circle that can be described about
! !    the polygon.
! !
! !    Output, real SIDE, the length of one side of the polygon.
! !
!   real angle
!   real area
!   integer n
!   real radin
!   real radout
!   real side
! !
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_INRAD_DATA_2D - Fatal error!'
!     write ( *, * ) '  Input value of N must be at least 3'
!     write ( *, * ) '  but your input value was N = ', n
!     stop
!   end if
!
!   angle = pi / real ( n )
!   area = n * radin * radin * tan ( angle )
!   side = 2.0E+00 * radin * tan ( angle )
!   radout = 0.5E+00 * side / sin ( angle )
!
!   return
! end subroutine
! subroutine polygon_outrad_data_2d ( area, n, radin, radout, side )
! !
! !*******************************************************************************
! !
! !! POLYGON_OUTRAD_DATA_2D determines polygonal data from its outer radius in 2D.
! !
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real AREA, the area of the regular polygon.
! !
! !    Input, integer N, the number of sides of the polygon.
! !    N must be at least 3.
! !
! !    Output, real RADIN, the inner radius of the polygon, that is,
! !    the radius of the largest circle that can be inscribed
! !    within the polygon.
! !
! !    Input, real RADOUT, the outer radius of the polygon, that is,
! !    the radius of the smallest circle that can be described
! !    around the polygon.
! !
! !    Output, real SIDE, the length of one side of the polygon.
! !
!   real angle
!   real area
!   integer n
!   real radin
!   real radout
!   real side
! !
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_OUTRAD_DATA_2D - Fatal error!'
!     write ( *, * ) '  Input value of N must be at least 3'
!     write ( *, * ) '  but your input value was N = ', n
!     stop
!   end if
!
!   angle = pi / real ( n )
!   area = 0.5E+00 * n * radout * radout * sin ( 2.0E+00 * angle )
!   side = 2.0E+00 * radout * sin ( angle )
!   radin = 0.5E+00 * side / tan ( angle )
!
!   return
! end subroutine
! subroutine polygon_side_data_2d ( area, n, radin, radout, side )
! !
! !*******************************************************************************
! !
! !! POLYGON_SIDE_DATA_2D determines polygonal data from its side length in 2D.
! !
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, real AREA, the area of the regular polygon.
! !
! !    Input, integer N, the number of sides of the polygon.
! !    N must be at least 3.
! !
! !    Output, real RADIN, the inner radius of the polygon, that is,
! !    the radius of the largest circle that can be inscribed within
! !    the polygon.
! !
! !    Output, real RADOUT, the outer radius of the polygon, that is,
! !    the radius of the smallest circle that can be described about
! !    the polygon.
! !
! !    Input, real SIDE, the length of one side of the polygon.
! !
!   real angle
!   real area
!   integer n
!   real radin
!   real radout
!   real side
! !
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYGON_SIDE_DATA_2D - Fatal error!'
!     write ( *, * ) '  Input value of N must be at least 3'
!     write ( *, * ) '  but your input value was N = ', n
!     stop
!   end if
!
!   angle = pi / real ( n )
!   area = 0.5E+00 * n * side * side / tan ( angle )
!   radin = 0.5E+00 * side / tan ( angle )
!   radout = 0.5E+00 * side / sin ( angle )
!
!   return
! end subroutine
! subroutine polyhedron_surface_3d ( coord, maxorder, nface, node, &
!   n, order, area )
! !
! !*******************************************************************************
! !
! !! POLYHEDRON_SURFACE_3D computes the surface area of a polyhedron in 3D.
! !
! !
! !  Restriction:
! !
! !    The computation is not valid unless the faces of the polyhedron
! !    are planar polygons.
! !
! !  Reference:
! !
! !    Allen Van Gelder,
! !    Efficient Computation of Polygon Area and Polyhedron Volume,
! !    Graphics Gems V, edited by Alan Paeth,
! !    AP Professional, 1995.
! !
! !  Modified:
! !
! !    18 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real COORD(3,N), the 3D coordinates of the vertices.
! !    The vertices may be listed in any order.
! !
! !    Input, integer MAXORDER, the maximum number of vertices that make
! !    up a face of the polyhedron.
! !
! !    Input, integer NFACE, the number of faces of the polyhedron.
! !
! !    Input, integer NODE(NFACE,MAXORDER).  Face I is defined by
! !    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
! !    are listed in neighboring order.
! !
! !    Input, integer N, the number of points stored in COORD.
! !
! !    Input, integer ORDER(NFACE), the number of vertices making up each face.
! !
! !    Output, real AREA, the total surface area of the polyhedron.
! !
!   integer maxorder
!   integer nface
!   integer n
! !
!   real ainc
!   real area
!   real coord(3,n)
!   integer i
!   integer j
!   integer k
!   integer node(nface,maxorder)
!   integer order(nface)
!   real x1
!   real x2
!   real x3
!   real x4
!   real y1
!   real y2
!   real y3
!   real y4
!   real z1
!   real z2
!   real z3
!   real z4
! !
!   area = 0.0E+00
!
!   do i = 1, nface
!
!     x4 = 0.0E+00
!     y4 = 0.0E+00
!     z4 = 0.0E+00
! !
! !  Compute the area vector for this face.
! !
!     do j = 1, order(i)
!
!       k = node(i,j)
!       x1 = coord(1,k)
!       y1 = coord(2,k)
!       z1 = coord(3,k)
!
!       if ( j < order(i) ) then
!         k = node(i,j+1)
!       else
!         k = node(i,1)
!       end if
!
!       x2 = coord(1,k)
!       y2 = coord(2,k)
!       z2 = coord(3,k)
!
!       call cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!       x4 = x4 + x3
!       y4 = y4 + y3
!       z4 = z4 + z3
!
!     end do
! !
! !  Add the magnitude of the area vector to the sum.
! !
!     ainc = sqrt ( x4 * x4 + y4 * y4 + z4 * z4 )
!     area = area + ainc
!
!   end do
!
!   area = 0.5E+00 * area
!
!   return
! end subroutine
! subroutine polyhedron_volume_3d ( coord, maxorder, nface, node, &
!   n, order, volume )
! !
! !*******************************************************************************
! !
! !! POLYHEDRON_VOLUME_3D computes the volume of a polyhedron in 3D.
! !
! !
! !  Restriction:
! !
! !    The computation is not valid unless the faces of the polyhedron
! !    are planar polygons.
! !
! !  Reference:
! !
! !    Allen Van Gelder,
! !    Efficient Computation of Polygon Area and Polyhedron Volume,
! !    Graphics Gems V, edited by Alan Paeth,
! !    AP Professional, 1995.
! !
! !  Modified:
! !
! !    18 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real COORD(3,N), the 3D coordinates of the vertices.
! !    The vertices may be listed in any order.
! !
! !    Input, integer MAXORDER, the maximum number of vertices that make
! !    up a face of the polyhedron.
! !
! !    Input, integer NFACE, the number of faces of the polyhedron.
! !
! !    Input, integer NODE(NFACE,MAXORDER).  Face I is defined by
! !    the vertices NODE(I,1) through NODE(I,ORDER(I)).  These vertices
! !    are listed in neighboring order.
! !
! !    Input, integer N, the number of points stored in COORD.
! !
! !    Input, integer ORDER(NFACE), the number of vertices making up
! !    each face.
! !
! !    Output, real VOLUME, the volume of the polyhedron.
! !
!   integer maxorder
!   integer nface
!   integer n
! !
!   real coord(3,n)
!   integer i
!   integer j
!   integer k
!   integer node(nface,maxorder)
!   integer order(nface)
!   real volume
!   real x1
!   real x2
!   real x3
!   real x4
!   real y1
!   real y2
!   real y3
!   real y4
!   real z1
!   real z2
!   real z3
!   real z4
! !
!   volume = 0.0E+00
!
!   do i = 1, nface
!
!     x4 = 0.0E+00
!     y4 = 0.0E+00
!     z4 = 0.0E+00
! !
! !  Compute the area vector for this face.
! !
!     do j = 1, order(i)
!
!       k = node(i,j)
!       x1 = coord(1,k)
!       y1 = coord(2,k)
!       z1 = coord(3,k)
!
!       if ( j < order(i) ) then
!         k = node(i,j+1)
!       else
!         k = node(i,1)
!       end if
!
!       x2 = coord(1,k)
!       y2 = coord(2,k)
!       z2 = coord(3,k)
!
!       call cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!       x4 = x4 + x3
!       y4 = y4 + y3
!       z4 = z4 + z3
!
!     end do
! !
! !  Area vector dot any vertex.
! !
!     k = node(i,1)
!     volume = volume + x4 * coord(1,k) + y4 * coord(2,k) + z4 * coord(3,k)
!
!   end do
!
!   volume = volume / 6.0E+00
!
!   return
! end subroutine
! subroutine polyline_index_point_nd ( maxpts, ndim, npts, t, xpts, x )
! !
! !*******************************************************************************
! !
! !! POLYLINE_INDEX_POINT_ND evaluates a polyline at a given arclength in ND.
! !
! !
! !  Discussion:
! !
! !    The polyline is defined as the set of M-1 line segments lying
! !    between a sequence of M points.  The arclength of a point lying
! !    on the polyline is simply the length of the broken line from the
! !    initial point.  Any point on the polyline can be found by
! !    specifying its arclength.
! !
! !    If the given arclength coordinate is less than 0, or greater
! !    than the arclength coordinate of the last given point, then
! !    extrapolation is used, that is, the first and last line segments
! !    are extended as necessary.
! !
! !    The arclength coordinate system measures the distance between
! !    any two points on the polyline as the length of the segment of the
! !    line that joins them.
! !
! !  Modified:
! !
! !    18 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAXPTS, the first dimension of the array
! !    used to hold XPTS.  MAXPTS must be at least equal to NPTS.
! !
! !    Input, integer NDIM, the dimension of the space in which
! !    the points lie.  The second dimension of XPTS.
! !
! !    Input, integer NPTS, the number of points.
! !
! !    Input, real T, the desired arclength coordinate.
! !
! !    Input, real XPTS(MAXPTS,NDIM), a set of NPTS coordinates
! !    in NDIM space, describing a set of points that define
! !    a polyline.
! !
! !    Output, real X(NDIM), a point lying on the polyline defined
! !    by XPTS, and having arclength coordinate T.
! !
!   integer maxpts
!   integer ndim
! !
!   integer i
!   integer j
!   integer npts
!   real s
!   real sum
!   real t
!   real tleft
!   real trite
!   real xpts(maxpts,ndim)
!   real x(ndim)
! !
!   if ( maxpts <= 0 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYLINE_INDEX_POINT_ND - Fatal error!'
!     write ( *, * ) '  The input quantity MAXPTS is nonpositive.'
!     write ( *, * ) '  MAXPTS = ', maxpts
!     stop
!   end if
!
!   if ( npts <= 0 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYLINE_INDEX_POINT_ND - Fatal error!'
!     write ( *, * ) '  Input quantity NPTS nonpositive.'
!     write ( *, * ) '  NPTS = ', npts
!     stop
!   end if
!
!   if ( maxpts < npts ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'POLYLINE_INDEX_POINT_ND - Fatal error!'
!     write ( *, * ) '  MAXPTS less than NPTS.'
!     write ( *, * ) '  MAXPTS = ', maxpts, ' NPTS = ', npts
!     stop
!   end if
!
!   if ( npts == 1 ) then
!
!     x(1:ndim) = xpts(1,1:ndim)
!
!   else
!
!     trite = 0.0E+00
!     do i = 1, npts - 1
! !
! !  Find the distance between points I and I+1.
! !
!       sum = 0.0E+00
!       do j = 1, ndim
!         sum = sum + ( xpts(i+1,j) - xpts(i,j) )**2
!       end do
!
!       tleft = trite
!       trite = trite + sqrt ( sum )
! !
! !  Interpolate or extrapolate in an interval.
! !
!       if ( t <= trite .or. i == npts - 1 ) then
!
!         s = ( t - tleft ) / ( trite - tleft )
!         x(1:ndim) = ( 1.0E+00 - s ) * xpts(i,1:ndim) + s * xpts(i+1,1:ndim)
!
!         return
!       end if
!     end do
!   end if
!
!   return
! end subroutine
! subroutine polyline_length_nd ( maxpts, ndim, npoint, spoint, length, xpoint )
! !
! !*******************************************************************************
! !
! !! POLYLINE_LENGTH_ND computes the length of a polyline in ND.
! !
! !
! !  Definition:
! !
! !    A polyline of order M is the geometric structure consisting of
! !    the M-1 line segments that lie between successive elements of a list
! !    of M points.
! !
! !    An ordinary line segment is a polyline of order 2.
! !    The letter "V" is a polyline of order 3.
! !    The letter "N" is a polyline of order 4, and so on.
! !
! !  Formula:
! !
! !    DIST(I+1,I) = sqrt ( SUM(j = 1,NDIM) ( X(I+1) - X(I) )**2 )
! !
! !    LENGTH = SUM (I = 1 to NPOINT-1) DIST(I+1,I)
! !
! !  Modified:
! !
! !    18 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAXPTS, the declared first dimension of XPOINT.
! !
! !    Input, integer NDIM, the number of dimensions of the points.
! !
! !    Input, integer NPOINT, the number of points.
! !
! !    Output, real SPOINT(MAXPTS), the arclength coordinates
! !    of each point.  The first point has SPOINT(1) = 0.
! !
! !    Output, real LENGTH, the length of the polyline,
! !    which is also equal to SPOINT(NPOINT).
! !
! !    Input, real XPOINT(MAXPTS,NDIM), the coordinates of the points.
! !
!   integer maxpts
!   integer ndim
! !
!   integer i
!   integer j
!   real length
!   integer npoint
!   real spoint(maxpts)
!   real temp
!   real xpoint(maxpts,ndim)
! !
!   length = 0.0E+00
!   spoint(1) = 0.0E+00
!
!   do i = 2, npoint
!
!     temp = 0.0E+00
!
!     do j = 1, ndim
!       temp = temp + ( xpoint(i,j) - xpoint(i-1,j) )**2
!     end do
!
!     temp = sqrt ( temp )
!
!     length = length + temp
!     spoint(i) = length
!
!   end do
!
!   return
! end subroutine
! subroutine proplane2 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, npnt, xp, yp, zp, &
!   alpha, beta )
! !
! !*******************************************************************************
! !
! !! PROPLANE2 produces 2D coordinates of points that lie in a plane, in 3D.
! !
! !
! !  Discussion:
! !
! !    The plane is specified by three non-colinear points, which we will
! !    call P1, P2 and P3.
! !
! !    The first thing to do is to compute two orthonormal vectors V1 and
! !    V2, so that any point P that lies in the plane may be written as
! !
! !      P = P1 + alpha * V1 + beta * V2
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
! !    coordinates of three points on the plane.  These three points
! !    should not lie in a straight line, but this condition is not
! !    checked.
! !
! !    Input, integer NPNT, the number of points to project.
! !
! !    Input, real XP(NPNT), YP(NPNT), ZP(NPNT), are the Cartesian
! !    coordinates of points which lie on the plane spanned by the
! !    three points.  These points are not checked to ensure that
! !    they lie on the plane.
! !
! !    Output, real ALPHA(NPNT), BETA(NPNT), the "in-plane" coordinates of
! !    the points.
! !
!   integer npnt
! !
!   real alpha(npnt)
!   real beta(npnt)
!   real dot
!   integer i
!   real v1(3)
!   real v2(3)
!   real x1
!   real x2
!   real x3
!   real xp(npnt)
!   real y1
!   real y2
!   real y3
!   real yp(npnt)
!   real z1
!   real z2
!   real z3
!   real zp(npnt)
! !
! !  Compute the two basis vectors for the affine plane.
! !
!   v1(1) = x2 - x1
!   v1(2) = y2 - y1
!   v1(3) = z2 - z1
!
!   call vector_unit_nd ( 3, v1 )
!
!   v2(1) = x3 - x1
!   v2(2) = y3 - y1
!   v2(3) = z3 - z1
!
!   dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
!
!   v2(1:3) = v2(1:3) - dot * v1(1:3)
!
!   call vector_unit_nd ( 3, v2 )
! !
! !  Now decompose each (X,Y,Z).
! !
!   do i = 1, npnt
!
!     alpha(i) = ( xp(i) - x1 ) * v1(1) +  ( yp(i) - y1 ) * v1(2) + &
!                ( zp(i) - z1 ) * v1(3)
!
!     beta(i) =  ( xp(i) - x1 ) * v2(1) + ( yp(i) - y1 ) * v2(2) + &
!                ( zp(i) - z1 ) * v2(3)
!
!   end do
!
!   return
! end subroutine
! subroutine proplane3 ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   npnt, xo, yo, zo, xp, yp, zp )
! !
! !*******************************************************************************
! !
! !! PROPLANE3 projects points orthographically onto a plane, in 3D.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the
! !    coordinates of three points on the plane.
! !
! !    Input, integer NPNT, the number of points to project.
! !
! !    Input, real XO(NPNT), YO(NPNT), ZO(NPNT), are the coordinates of
! !    the object points.
! !
! !    Output, real XP(NPNT), YP(NPNT), ZP(NPNT), are the coordinates of
! !    the projections of the object points through the focus point onto
! !    the plane.
! !
! !    XP, YP, and ZP may share the same memory as XO, YO, and ZO, in
! !    which case the projections will overwrite the original data.
! !
!   integer npnt
! !
!   real a
!   real b
!   real c
!   real d
!   integer i
!   real x1
!   real x2
!   real x3
!   real xo(npnt)
!   real xp(npnt)
!   real y1
!   real y2
!   real y3
!   real yo(npnt)
!   real yp(npnt)
!   real z1
!   real z2
!   real z3
!   real zo(npnt)
!   real zp(npnt)
! !
! !  Put the plane into ABCD form.
! !
!   call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
! !
! !  For each point, its aimage in the plane is the nearest point
! !  in the plane.
! !
!   do i = 1, npnt
!
!     call plane_imp_point_near_3d ( a, b, c, d, xo(i), yo(i), &
!       zo(i), xp(i), yp(i), zp(i) )
!
!   end do
!
!   return
! end subroutine
! subroutine provec ( base, m, n, vecm, vecn, vecnm )
! !
! !*******************************************************************************
! !
! !! PROVEC projects a vector from M space into N space.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real BASE(M,N).  The columns of BASE contain
! !    N vectors, each of length M, which form the basis for
! !    a space of dimension N.
! !
! !    Input, integer M, the dimension of the higher order space.
! !
! !    Input, integer N, the dimension of the lower order space.
! !
! !    Input, real VECM(M), is an M dimensional vector.
! !
! !    Output, real VECN(N), the projection of VECM into the
! !    lower dimensional space.  These values represent
! !    coordinates in the lower order space.
! !
! !    Output, real VECNM(M), the projection of VECM into the
! !    lower dimensional space, but using coordinates in
! !    the higher dimensional space.
! !
!   integer m
!   integer n
! !
!   real base(m,n)
!   integer i
!   integer j
!   integer k
!   real temp
!   real vecm(m)
!   real vecn(n)
!   real vecnm(m)
! !
! !  For each vector, remove all projections onto previous vectors,
! !  and then normalize.  This should result in a matrix BASE
! !  whose columns are orthonormal.
! !
!   do j = 1, n
!
!     do i = 1, j-1
!
!       temp = 0.0E+00
!       do k = 1, m
!         temp = temp + base(k,i) * base(k,j)
!       end do
!
!       base(1:m,j) = base(1:m,j) - temp * base(1:m,i)
!
!     end do
!
!     temp = 0.0E+00
!     do i = 1, m
!       temp = temp + base(i,j)**2
!     end do
!     temp = sqrt ( temp )
!
!     if ( temp > 0.0E+00 ) then
!       base(1:m,j) = base(1:m,j) / temp
!     end if
!
!   end do
! !
! !  Compute the coordinates of the projection of the vector
! !  simply by taking dot products.
! !
!   do j = 1, n
!     temp = 0.0E+00
!     do k = 1, m
!       temp = temp + base(k,j) * vecm(k)
!     end do
!     vecn(j) = temp
!   end do
! !
! !  Compute the coordinates of the projection in terms of
! !  the original space.
! !
!   do i = 1, m
!     vecnm(i) = 0.0E+00
!     do j = 1, n
!       vecnm(i) = vecnm(i) + base(i,j) * vecn(j)
!     end do
!   end do
!
!   return
! end subroutine
! subroutine pyramid_volume_3d ( h, s, volume )
! !
! !*******************************************************************************
! !
! !! PYRAMID_VOLUME_3D computes the volume of a pyramid with square base in 3D.
! !
! !
! !  Modified:
! !
! !    10 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real H, R, the height of the pyramid, and the length of one
! !    side of the square base.
! !
! !    Output, real VOLUME, the volume of the pyramid.
! !
!   real h
!   real s
!   real volume
! !
!   volume = s * s * h / 3.0E+00
!
!   return
! end subroutine
! subroutine quad_area_2d ( x1, y1, x2, y2, x3, y3, x4, y4, area )
! !
! !*******************************************************************************
! !
! !! QUAD_AREA_2D computes the area of a quadrilateral in 2D.
! !
! !
! !  Discussion:
! !
! !    This algorithm should be able to handle nonconvex quadrilaterals.
! !
! !  Modified:
! !
! !    23 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, X4, Y4, the X and Y coordinates
! !    of the corners of the quadrilateral.  The corners should be
! !    specified in clockwise or counterclockwise order.
! !
! !    Output, real AREA, the absolute area of the quadrilateral.
! !
!   real area
!   real area1
!   real area2
!   real x1
!   real x2
!   real x3
!   real x4
!   real y1
!   real y2
!   real y3
!   real y4
!  !
!   call triangle_area_2d ( x1, y1, x2, y2, x3, y3, area1 )
!
!   call triangle_area_2d ( x3, y3, x4, y4, x1, y1, area2 )
!
!   area = area1 + area2
!
!   return
! end subroutine
! subroutine quad_contains_point_2d ( x1, y1, x2, y2, x3, y3, x4, y4, x, y, &
!   inside )
! !
! !*******************************************************************************
! !
! !! QUAD_CONTAINS_POINT_2D finds if a point is inside a convex quadrilateral in 2D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, X4, Y4, the X and Y coordinates
! !    of the quadrilateral.
! !
! !    Input, real X, Y, the point to be checked.
! !
! !    Output, logical INSIDE, is .TRUE. if ( X,Y) is inside
! !    the quadrilateral or on its boundary, and .FALSE. otherwise.
! !
!   real angle_123
!   real angle_12x
!   real angle_234
!   real angle_23x
!   real angle_341
!   real angle_34x
!   real angle_412
!   real angle_41x
!   logical inside
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
! !
! !  This will only handle convex quadrilaterals.
! !
!   inside = .false.
!
!   angle_12x = anglei_rad_2d ( x1, y1, x2, y2, x, y )
!   angle_123 = anglei_rad_2d ( x1, y1, x2, y2, x3, y3 )
!
!   if ( angle_12x > angle_123 ) then
!     return
!   end if
!
!   angle_23x = anglei_rad_2d ( x2, y2, x3, y3, x, y )
!   angle_234 = anglei_rad_2d ( x2, y2, x3, y3, x4, y4 )
!
!   if ( angle_23x > angle_234 ) then
!     return
!   end if
!
!   angle_34x = anglei_rad_2d ( x3, y3, x4, y4, x, y )
!   angle_341 = anglei_rad_2d ( x3, y3, x4, y4, x1, y1 )
!
!   if ( angle_34x > angle_341 ) then
!     return
!   end if
!
!   angle_41x = anglei_rad_2d ( x4, y4, x1, y1, x, y )
!   angle_412 = anglei_rad_2d ( x4, y4, x1, y1, x2, y2 )
!
!   if ( angle_41x > angle_412 ) then
!     return
!   end if
!
!   inside = .true.
!
!   return
! end subroutine
! subroutine quad_point_dist_2d ( x1, y1, x2, y2, x3, y3, x4, y4, x, y, dist )
! !
! !*******************************************************************************
! !
! !! QUAD_POINT_DIST_2D finds the distance from a point to a quadrilateral in 2D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, X4, Y4, the corners of the
! !    quadrilateral.
! !
! !    Input, real X, Y, the point which is to be checked.
! !
! !    Output, real DIST, the distance from the point to the quadrilateral.
! !    DIST is zero if the point lies exactly on the quadrilateral.
! !
!   real dist
!   real dist2
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
! !
!   dist =             enorm0_2d ( x, y, x1, y1 )
!   dist = min ( dist, enorm0_2d ( x, y, x2, y2 ) )
!   dist = min ( dist, enorm0_2d ( x, y, x3, y3 ) )
!   dist = min ( dist, enorm0_2d ( x, y, x4, y4 ) )
!
!   call line_seg_point_dist_2d ( x1, y1, x2, y2, x, y, dist2 )
!   dist = min ( dist, dist2 )
!
!   call line_seg_point_dist_2d ( x2, y2, x3, y3, x, y, dist2 )
!   dist = min ( dist, dist2 )
!
!   call line_seg_point_dist_2d ( x3, y3, x4, y4, x, y, dist2 )
!   dist = min ( dist, dist2 )
!
!   call line_seg_point_dist_2d ( x4, y4, x1, y1, x, y, dist2 )
!   dist = min ( dist, dist2 )
!
!   return
! end subroutine
! subroutine quad_point_dist_signed_2d ( x1, y1, x2, y2, x3, y3, &
!   x4, y4, x, y, dist_signed )
! !
! !*******************************************************************************
! !
! !! QUAD_POINT_DIST_SIGNED_2D: signed distanct ( quadrilateral, point ) in 2D.
! !
! !
! !  Discussion:
! !
! !    The quadrilateral must be convex.
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, X4, Y4, the corners of the
! !    quadrilateral, given in order.  It is possible to list the points
! !    in such an order that the quadrilateral is not convex.  This
! !    routine does not check for convexity, and will produce erroneous
! !    results in such a case.
! !
! !    Input, real X, Y, the point which is to be checked.
! !
! !    Output, real DIST_SIGNED, the signed distance from the point to the
! !    convex quadrilateral.  DIST_SIGNED is actually the maximum of the
! !    signed distances from the point to each of the four lines that make
! !    up the quadrilateral.
! !
! !    (Essentially, if the point is outside the convex quadrilateral,
! !    only one of the signed distances can be positive, or two can
! !    be positive and equal.)
! !
! !    If DIST_SIGNED is
! !    0, the point is on the boundary;
! !    negative, the point is in the interior;
! !    positive, the point is in the exterior.
! !
!   real dis
!   real dis12
!   real dis23
!   real dis34
!   real dis41
!   real dist_signed
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real xm
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
!   real ym
! !
! !  Compare the signed distance from each line segment to the point,
! !  with the signed distance to the midpoint of the opposite line.
! !
! !  The signed distances should all be negative if the point is inside
! !  the triangle.
! !
!   call line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dis12 )
!
!   xm = 0.5E+00 * ( x3 + x4 )
!   ym = 0.5E+00 * ( y3 + y4 )
!
!   call line_exp_point_dist_signed_2d ( x1, y1, x2, y2, xm, ym, dis )
!
!   if ( dis > 0.0E+00 ) then
!     dis = - dis
!     dis12 = - dis12
!   end if
!
!   call line_exp_point_dist_signed_2d ( x2, y2, x3, y3, x, y, dis23 )
!
!   xm = 0.5E+00 * ( x4 + x1 )
!   ym = 0.5E+00 * ( y4 + y1 )
!
!   call line_exp_point_dist_signed_2d ( x2, y2, x3, y3, xm, ym, dis )
!
!   if ( dis > 0.0E+00 ) then
!     dis = - dis
!     dis23 = - dis23
!   end if
!
!   call line_exp_point_dist_signed_2d ( x3, y3, x4, y4, x, y, dis34 )
!
!   xm = 0.5E+00 * ( x1 + x2 )
!   ym = 0.5E+00 * ( y1 + y2 )
!
!   call line_exp_point_dist_signed_2d ( x3, y3, x4, y4, xm, ym, dis )
!
!   if ( dis > 0.0E+00 ) then
!     dis = - dis
!     dis34 = - dis34
!   end if
!
!   call line_exp_point_dist_signed_2d ( x4, y4, x1, y1, x, y, dis41 )
!
!   xm = 0.5E+00 * ( x2 + x3 )
!   ym = 0.5E+00 * ( y2 + y3 )
!
!   call line_exp_point_dist_signed_2d ( x4, y4, x1, y1, xm, ym, dis )
!
!   if ( dis > 0.0E+00 ) then
!     dis = - dis
!     dis41 = - dis41
!   end if
!
!   dist_signed = max ( dis12, dis23, dis34, dis41 )
!
!   return
! end subroutine
! subroutine quat_conj ( q )
! !
! !*******************************************************************************
! !
! !! QUAT_CONJ conjugates a quaternion.
! !
! !
! !  Discussion:
! !
! !    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
! !    may be written as
! !
! !      Q = A + Bi + Cj + Dk.
! !
! !    The conjugate of Q is
! !
! !      conj ( Q ) = A - Bi -Cj -Dk.
! !
! !  Modified:
! !
! !    29 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, real Q(4).  On input, the quaternion to be conjugated.
! !    On output, the conjugated quaternion.
! !
!   real q(4)
! !
!   q(2) = - q(2)
!   q(3) = - q(3)
!   q(4) = - q(4)
!
!   return
! end subroutine
! subroutine quat_inv ( q )
! !
! !*******************************************************************************
! !
! !! QUAT_INV inverts a quaternion.
! !
! !
! !  Discussion:
! !
! !    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
! !    may be written as
! !
! !      Q = A + Bi + Cj + Dk.
! !
! !    The inverse of Q is
! !
! !      inverse ( Q ) = conjugate ( Q ) / ( norm ( Q ) )**2.
! !
! !  Modified:
! !
! !    29 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, real Q(4).  On input, the quaternion to be inverted.
! !    On output, the inverse of the input quaternion.
! !
!   real normsq
!   real q(4)
!   real quat_norm
! !
!   normsq = ( quat_norm ( q ) )**2
!
!   q(1:4) = q(1:4) / normsq
!   q(1:3) = - q(1:3)
!
!   return
! end subroutine
! subroutine quat_mul ( q1, q2, q3 )
! !
! !*******************************************************************************
! !
! !! QUAT_MUL multiplies two quaternions.
! !
! !
! !  Discussion:
! !
! !    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
! !    may be written as
! !
! !      Q = A + Bi + Cj + Dk.
! !
! !    To multiply two quaternions, use the relationships:
! !
! !      ij = -ji = k
! !      jk = -kj = i
! !      ki = -ik = j
! !      ii = jj = kk = -1
! !
! !  Modified:
! !
! !    29 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real Q1(4), Q2(4), the two quaternions to be multiplied.
! !
! !    Output, real Q3(4), the product of the two quaternions.
! !
!   real q1(4)
!   real q2(4)
!   real q3(4)
! !
!   q3(1) = q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3) - q1(4) * q2(4)
!   q3(2) = q1(1) * q2(2) + q1(2) * q2(1) + q1(3) * q2(4) - q1(4) * q2(3)
!   q3(3) = q1(1) * q2(3) - q1(2) * q2(4) + q1(3) * q2(1) + q1(4) * q2(2)
!   q3(4) = q1(1) * q2(4) + q1(2) * q2(3) - q1(3) * q2(2) + q1(4) * q2(1)
!
!   return
! end subroutine
! function quat_norm ( q )
! !
! !*******************************************************************************
! !
! !! QUAT_NORM computes the norm of a quaternion.
! !
! !
! !  Discussion:
! !
! !    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
! !    may be written as
! !
! !      Q = A + Bi + Cj + Dk.
! !
! !    The norm of Q is
! !
! !      norm(Q) = sqrt ( A**2 + B**2 + C**2 + D**2 ).
! !
! !  Modified:
! !
! !    29 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real Q(4), the quaternion.
! !
! !    Output, real QUAT_NORM, the norm of the quaternion.
! !
!   real q(4)
!   real quat_norm
! !
!   quat_norm = sqrt ( q(1) * q(1) + q(2) * q(2) + q(3) * q(3) + q(4) * q(4) )
!
!   return
! end function
!
! function r_modp ( x, y )
! !
! !*******************************************************************************
! !
! !! R_MODP returns the nonnegative remainder of real division.
! !
! !
! !  Formula:
! !
! !    If
! !      REM = R_MODP ( X, Y )
! !      RMULT = ( X - REM ) / Y
! !    then
! !      X = Y * RMULT + REM
! !    where REM is always nonnegative.
! !
! !  Comments:
! !
! !    The MOD function computes a result with the same sign as the
! !    quantity being divided.  Thus, suppose you had an angle A,
! !    and you wanted to ensure that it was between 0 and 360.
! !    Then mod(A,360.0) would do, if A was positive, but if A
! !    was negative, your result would be between -360 and 0.
! !
! !    On the other hand, R_MODP(A,360.0) is between 0 and 360, always.
! !
! !  Examples:
! !
! !        I         J     MOD  R_MODP   R_MODP Factorization
! !
! !      107        50       7       7    107 =  2 *  50 + 7
! !      107       -50       7       7    107 = -2 * -50 + 7
! !     -107        50      -7      43   -107 = -3 *  50 + 43
! !     -107       -50      -7      43   -107 =  3 * -50 + 43
! !
! !  Modified:
! !
! !    24 July 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X, the number to be divided.
! !
! !    Input, real Y, the number that divides X.
! !
! !    Output, real R_MODP, the nonnegative remainder when X is divided by Y.
! !
!   real r_modp
!   real x
!   real y
! !
!   if ( y == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'R_MODP - Fatal error!'
!     write ( *, * ) '  R_MODP ( X, Y ) called with Y = ', y
!     stop
!   end if
!
!   r_modp = mod ( x, y )
!
!   if ( r_modp < 0.0E+00 ) then
!     r_modp = r_modp + abs ( y )
!   end if
!
!   return
! end function
! subroutine r_random ( rlo, rhi, r )
! !
! !*******************************************************************************
! !
! !! R_RANDOM returns a random real in a given range.
! !
! !
! !  Modified:
! !
! !    01 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real RLO, RHI, the minimum and maximum values.
! !
! !    Output, real R, the randomly chosen value.
! !
!   logical, save :: seed = .false.
!   real r
!   real rhi
!   real rlo
!   real t
! !
!   if ( .not. seed ) then
!     call random_seed
!     seed = .true.
!   end if
! !
! !  Pick a random number in (0,1).
! !
!   call random_number ( harvest = t )
! !
! !  Set R.
! !
!   r = ( 1.0E+00 - t ) * rlo + t * rhi
!
!   return
! end subroutine
! subroutine r_swap ( x, y )
! !
! !*******************************************************************************
! !
! !! R_SWAP switches two real values.
! !
! !
! !  Modified:
! !
! !    01 May 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, real X, Y.  On output, the values of X and
! !    Y have been interchanged.
! !
!   real x
!   real y
!   real z
! !
!   z = x
!   x = y
!   y = z
!
!   return
! end subroutine
!
!
! subroutine radec_separation ( ra1, dec1, ra2, dec2, theta )
! !
! !*******************************************************************************
! !
! !! RADEC_SEPARATION: separation between right ascension/declination points.
! !
! !
! !  Discussion:
! !
! !    Right ascension is measured in hours, between 0 and 24, and
! !    essentially measures longitude.
! !
! !    Declination measures the angle from the equator towards the north pole,
! !    and ranges from -90 (South Pole) to 90 (North Pole).
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real RA1, DEC1, RA2, DEC2, the right ascension and declination
! !    of the two points.
! !
! !    Output, real THETA, the angular separation between the points,
! !    in radians.
! !
!   real dec1
!   real cos_theta
!   real dec2
!   real norm_v1
!   real norm_v2
!   real phi1
!   real phi2
!   real ra1
!   real ra2
!   real theta
!   real theta1
!   real theta2
!   real v1(3)
!   real v2(3)
! !
!   theta1 = degrees_to_radians ( 15.0E+00 * ra1 )
!   phi1 = degrees_to_radians ( dec1 )
!
!   v1(1:3) = (/ cos ( theta1 ) * cos ( phi1 ), &
!                sin ( theta1 ) * cos ( phi1 ), &
!                                 sin ( phi1 ) /)
!
!   norm_v1 = sqrt ( sum ( v1(1:3)**2 ) )
!
!   theta2 = degrees_to_radians ( 15.0E+00 * ra2 )
!   phi2 = degrees_to_radians ( dec2 )
!
!   v2(1:3) = (/ cos ( theta2 ) * cos ( phi2 ), &
!                sin ( theta2 ) * cos ( phi2 ), &
!                                 sin ( phi2 ) /)
!
!   norm_v2 = sqrt ( sum ( v2(1:3)**2 ) )
!
!   cos_theta = dot_product ( v1(1:3), v2(1:3) ) / ( norm_v1 * norm_v2 )
!
!   theta = arc_cosine ( cos_theta )
!
!   return
! end subroutine
!
!
!
! subroutine radec_to_xyz ( ra, dec, x, y, z )
! !
! !*******************************************************************************
! !
! !! RADEC_TO_XYZ converts right ascension/declination to (X,Y,Z) coordinates.
! !
! !
! !  Discussion:
! !
! !    Right ascension is measured in hours, between 0 and 24, and
! !    essentially measures longitude.
! !
! !    Declination measures the angle from the equator towards the north pole,
! !    and ranges from -90 (South Pole) to 90 (North Pole).
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real RA, DEC, the right ascension and declination of a point.
! !
! !    Output, real X, Y, Z, the corresponding coordinates of a point with
! !    radius 1.
! !
!   real dec
!   real phi
!   real ra
!   real theta
!   real x
!   real y
!   real z
! !
!   theta = degrees_to_radians ( 15.0E+00 * ra )
!   phi = degrees_to_radians ( dec )
!
!   x = cos ( theta ) * cos ( phi )
!   y = sin ( theta ) * cos ( phi )
!   z = sin ( phi )
!
!   return
! end subroutine
!
! function radians_to_degrees ( angle )
! !
! !*******************************************************************************
! !
! !! RADIANS_TO_DEGREES converts an angle from radians to degrees.
! !
! !
! !  Modified:
! !
! !    10 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real ANGLE, an angle in radians.
! !
! !    Output, real RADIANS_TO_DEGREES, the equivalent angle in degrees.
! !
!   real radians_to_degrees,angle
!   radians_to_degrees = ( angle / pi ) * 180.0E+00
!
!   return
! end function
!
! subroutine rmat2_inverse ( a, b, det )
! !
! !*******************************************************************************
! !
! !! RMAT2_INVERSE inverts a 2 by 2 real matrix using Cramer's rule.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(2,2), the matrix to be inverted.
! !
! !    Output, real B(2,2), the inverse of the matrix A.
! !
! !    Output, real DET, the determinant of the matrix A.
! !
! !    If DET is zero, then A is singular, and does not have an
! !    inverse.  In that case, B is simply set to zero, and a
! !    message is printed.
! !
! !    If DET is nonzero, then its value is roughly an estimate
! !    of how nonsingular the matrix A is.
! !
!   real a(2,2)
!   real b(2,2)
!   real det
! !
! !  Compute the determinant.
! !
!   det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
! !
! !  If the determinant is zero, bail out.
! !
!   if ( det == 0.0E+00 ) then
!
!     b(1:2,1:2) = 0.0E+00
!
!     return
!   end if
! !
! !  Compute the entries of the inverse matrix using an explicit formula.
! !
!   b(1,1) = + a(2,2) / det
!   b(1,2) = - a(1,2) / det
!   b(2,1) = - a(2,1) / det
!   b(2,2) = + a(1,1) / det
!
!   return
! end subroutine
! function rmat2_det ( a )
! !
! !*******************************************************************************
! !
! !! RMAT2_DET computes the determinant of a 2 by 2 matrix.
! !
! !
! !  Discussion:
! !
! !    The determinant is the area spanned by the vectors making up the rows
! !    or columns of the matrix.
! !
! !  Formula:
! !
! !    RMAT2_DET = A(1,1) * A(2,2) - A(1,2) * A(2,1).
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(2,2), the matrix whose determinant is desired.
! !
! !    Output, real RMAT2_DET, the determinant of the matrix.
! !
!   real a(2,2)
!   real rmat2_det
! !
!   rmat2_det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!   return
! end function
!
! function rmat3_det ( a )
! !
! !*******************************************************************************
! !
! !! RMAT3_DET computes the determinant of a 3 by 3 matrix.
! !
! !
! !  The determinant is the volume of the shape spanned by the vectors
! !  making up the rows or columns of the matrix.
! !
! !  Formula:
! !
! !    det = a11 * a22 * a33 - a11 * a23 * a32
! !        + a12 * a23 * a31 - a12 * a21 * a33
! !        + a13 * a21 * a32 - a13 * a22 * a31
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(3,3), the matrix whose determinant is desired.
! !
! !    Output, real RMAT3_DET, the determinant of the matrix.
! !
!   real a(3,3)
!   real rmat3_det
! !
!   rmat3_det =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
!               + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
!               + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!
!   return
! end function
! subroutine rmat3_inverse ( a, b, det )
! !
! !*******************************************************************************
! !
! !! RMAT3_INVERSE inverts a 3 by 3 real matrix using Cramer's rule.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(3,3), the matrix to be inverted.
! !
! !    Output, real B(3,3), the inverse of the matrix A.
! !
! !    Output, real DET, the determinant of the matrix A.
! !
! !    If DET is zero, then A is singular, and does not have an
! !    inverse.  In that case, B is simply set to zero, and a
! !    message is printed.
! !
! !    If DET is nonzero, then its value is roughly an estimate
! !    of how nonsingular the matrix A is.
! !
!   real a(3,3)
!   real b(3,3)
!   real det
!   integer i
!   integer j
! !
! !  Compute the determinant of A
! !
!   det =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
!         + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
!         + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
! !
! !  If the determinant is zero, bail out.
! !
!   if ( det == 0.0E+00 ) then
!
!     b(1:3,1:3) = 0.0E+00
!
!     return
!   end if
! !
! !  Compute the entries of the inverse matrix using an explicit
! !  formula.
! !
!   b(1,1) = + ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
!   b(1,2) = - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
!   b(1,3) = + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det
!
!   b(2,1) = - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
!   b(2,2) = + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
!   b(2,3) = - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det
!
!   b(3,1) = + ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
!   b(3,2) = - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
!   b(3,3) = + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det
!
!   return
! end subroutine
! function rmat4_det ( a )
! !
! !*******************************************************************************
! !
! !! RMAT4_DET computes the determinant of a 4 by 4 matrix.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the matrix whose determinant is desired.
! !
! !    Output, real RMAT4_DET, the determinant of the matrix.
! !
!   real a(4,4)
!   real rmat4_det
! !
!   rmat4_det = &
!       a(1,1) * ( &
!         a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
!       - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
!       + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
!     - a(1,2) * ( &
!         a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
!       - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
!       + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
!     + a(1,3) * ( &
!         a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
!       - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
!       + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
!     - a(1,4) * ( &
!         a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
!       - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
!       + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )
!
!   return
! end function
! function rmat5_det ( a )
! !
! !*******************************************************************************
! !
! !! RMAT5_DET computes the determinant of a 5 by 5 matrix.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(5,5), the matrix whose determinant is desired.
! !
! !    Output, real RMAT5_DET, the determinant of the matrix.
! !
!   real a(5,5)
!   real b(4,4)
!   real rmat5_det
!   integer i
!   integer inc
!   integer j
!   integer k
! !
! !  Expand the determinant into the sum of the determinants of the
! !  five 4 by 4 matrices created by dropping row 1, and column k.
! !
!   rmat5_det = 0.0E+00
!
!   do k = 1, 5
!
!     do i = 1, 4
!       do j = 1, 4
!
!         if ( j < k ) then
!           inc = 0
!         else
!           inc = 1
!         end if
!
!         b(i,j) = a(i+1,j+inc)
!
!       end do
!     end do
!
!     rmat5_det = rmat5_det + (-1)**( k + 1 ) * a(1,k) * rmat4_det ( b )
!
!   end do
!
!   return
! end function
! subroutine rmat_print ( lda, m, n, a, title )
! !
! !*******************************************************************************
! !
! !! RMAT_PRINT prints a real matrix.
! !
! !
! !  Modified:
! !
! !    24 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer LDA, the leading dimension of A.
! !
! !    Input, integer M, the number of rows in A.
! !
! !    Input, integer N, the number of columns in A.
! !
! !    Input, real A(LDA,N), the matrix to be printed.
! !
! !    Input, character ( len = * ) TITLE, a title to be printed first.
! !    TITLE may be blank.
! !
!   integer lda
!   integer n
! !
!   real a(lda,n)
!   integer i
!   integer j
!   integer jhi
!   integer jlo
!   integer m
!   character ( len = * ) title
! !
!   if ( title /= ' ' ) then
!     write ( *, * ) ' '
!     write ( *, '(a)' ) trim ( title )
!   end if
!
!   do jlo = 1, n, 5
!     jhi = min ( jlo + 4, n )
!     write ( *, * ) ' '
!     write ( *, '(6x,5(i7,7x))' ) ( j, j = jlo, jhi )
!     write ( *, * ) ' '
!     do i = 1, m
!       write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi)
!     end do
!   end do
!
!   return
! end subroutine
! subroutine rmat_solve ( a, n, nrhs, info )
! !
! !*******************************************************************************
! !
! !! RMAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
! !
! !
! !  Modified:
! !
! !    08 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, real A(N,N+NRHS), contains in rows and columns 1
! !    to N the coefficient matrix, and in columns N+1 through
! !    N+NRHS, the right hand sides.  On output, the coefficient matrix
! !    area has been destroyed, while the right hand sides have
! !    been overwritten with the corresponding solutions.
! !
! !    Input, integer NRHS, the number of right hand sides.  NRHS
! !    must be at least 0.
! !
! !    Output, integer INFO, singularity flag.
! !    0, the matrix was not singular, the solutions were computed;
! !    J, factorization failed on step J, and the solutions could not
! !    be computed.
! !
!   integer n
!   integer nrhs
! !
!   real a(n,n+nrhs)
!   real apivot
!   real factor
!   integer i
!   integer info
!   integer ipivot
!   integer j
!   integer k
!   real temp
! !
!   info = 0
!
!   do j = 1, n
! !
! !  Choose a pivot row.
! !
!     ipivot = j
!     apivot = a(j,j)
!
!     do i = j+1, n
!       if ( abs ( a(i,j) ) > abs ( apivot ) ) then
!         apivot = a(i,j)
!         ipivot = i
!       end if
!     end do
!
!     if ( apivot == 0.0E+00 ) then
!       info = j
!       return
!     end if
! !
! !  Interchange.
! !
!     do i = 1, n + nrhs
!       call r_swap ( a(ipivot,i), a(j,i) )
!     end do
! !
! !  A(J,J) becomes 1.
! !
!     a(j,j) = 1.0E+00
!     a(j,j+1:n+nrhs) = a(j,j+1:n+nrhs) / apivot
! !
! !  A(I,J) becomes 0.
! !
!     do i = 1, n
!
!       if ( i /= j ) then
!
!         factor = a(i,j)
!         a(i,j) = 0.0E+00
!         a(i,j+1:n+nrhs) = a(i,j+1:n+nrhs) - factor * a(j,j+1:n+nrhs)
!
!       end if
!
!     end do
!
!   end do
!
!   return
! end subroutine
! subroutine rotation_axis_vector_3d ( axis, angle, v, w )
! !
! !*******************************************************************************
! !
! !! ROTATION_AXIS_VECTOR_3D rotates a vector around an axis vector in 3D.
! !
! !
! !  Modified:
! !
! !    31 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real AXIS(3), the axis vector for the rotation.
! !
! !    Input, real ANGLE, the angle, in radians, of the rotation.
! !
! !    Input, real V(3), the vector to be rotated.
! !
! !    Output, real W(3), the rotated vector.
! !
!   real angle
!   real axis(3)
!   real dot
!   real norm
!   real norm2
!   real rot(3)
!   real v(3)
!   real w(3)
!   real xa
!   real xn
!   real xn2
!   real xp
!   real ya
!   real yn
!   real yn2
!   real yp
!   real za
!   real zn
!   real zn2
!   real zp
! !
! !  Compute the length of the rotation axis.
! !
!   xa = axis(1)
!   ya = axis(2)
!   za = axis(3)
!
!   norm = sqrt ( xa * xa + ya * ya + za * za )
!
!   if ( norm == 0.0E+00 ) then
!     return
!   end if
!
!   xa = xa / norm
!   ya = ya / norm
!   za = za / norm
! !
! !  Compute the dot product of the vector and the rotation axis.
! !
!   dot = v(1) * xa + v(2) * ya + v(3) * za
! !
! !  Compute the parallel component of the vector.
! !
!   xp = dot * xa
!   yp = dot * ya
!   zp = dot * za
! !
! !  Compute the normal component of the vector.
! !
!   xn = v(1) - xp
!   yn = v(2) - yp
!   zn = v(3) - zp
!
!   norm2 = sqrt ( xn * xn + yn * yn + zn * zn )
!
!   if ( norm2 == 0.0E+00 ) then
!     w(1:3) = (/ xp, yp, zp /)
!     return
!   end if
!
!   xn = xn / norm2
!   yn = yn / norm2
!   zn = zn / norm2
! !
! !  Compute a second vector, lying in the plane, perpendicular
! !  to V, and forming a right-handed system.
! !
!   call cross_3d ( xa, ya, za, xn, yn, zn, xn2, yn2, zn2 )
!
!   norm = sqrt ( xn2 * xn2 + yn2 * yn2 + zn2 * zn2 )
!
!   xn2 = xn2 / norm
!   yn2 = yn2 / norm
!   zn2 = zn2 / norm
! !
! !  Rotate the normal component by the angle.
! !
!   rot(1) = norm2 * ( cos ( angle ) * xn + sin ( angle ) * xn2 )
!   rot(2) = norm2 * ( cos ( angle ) * yn + sin ( angle ) * yn2 )
!   rot(3) = norm2 * ( cos ( angle ) * zn + sin ( angle ) * zn2 )
! !
! !  The rotated vector is the parallel component plus the rotated component.
! !
!   w(1) = xp + rot(1)
!   w(2) = yp + rot(2)
!   w(3) = zp + rot(3)
!
!   return
! end subroutine
! subroutine rotation_axis2mat_3d ( axis, angle, a )
! !
! !*******************************************************************************
! !
! !! ROTATION_AXIS2MAT_3D converts a rotation from axis to matrix format in 3D.
! !
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Modified:
! !
! !    27 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real AXIS(3), the axis vector which remains unchanged by
! !    the rotation.
! !
! !    Input, real ANGLE, the angular measurement of the rotation about
! !    the axis, in radians.
! !
! !    Output, real A(3,3), the rotation matrix.
! !
!   real a(3,3)
!   real angle
!   real axis(3)
!   real ca
!   real norm
!   real sa
!   real v1
!   real v2
!   real v3
! !
!   v1 = axis(1)
!   v2 = axis(2)
!   v3 = axis(3)
!
!   norm = sqrt ( v1 * v1 + v2 * v2 + v3 * v3 )
!
!   if ( norm == 0.0E+00 ) then
!     return
!   end if
!
!   v1 = v1 / norm
!   v2 = v2 / norm
!   v3 = v3 / norm
!
!   ca = cos ( angle )
!   sa = sin ( angle )
!
!   a(1,1) =                v1 * v1 + ca * ( 1.0E+00 - v1 * v1 )
!   a(1,2) = ( 1.0E+00 - ca ) * v1 * v2 - sa * v3
!   a(1,3) = ( 1.0E+00 - ca ) * v1 * v3 + sa * v2
!
!   a(2,1) = ( 1.0E+00 - ca ) * v2 * v1 + sa * v3
!   a(2,2) =                v2 * v2 + ca * ( 1.0E+00 - v2 * v2 )
!   a(2,3) = ( 1.0E+00 - ca ) * v2 * v3 - sa * v1
!
!   a(3,1) = ( 1.0E+00 - ca ) * v3 * v1 - sa * v2
!   a(3,2) = ( 1.0E+00 - ca ) * v3 * v2 + sa * v1
!   a(3,3) =                v3 * v3 + ca * ( 1.0E+00 - v3 * v3 )
!
!   return
! end subroutine
! subroutine rotation_axis2quat_3d ( axis, angle, q )
! !
! !*******************************************************************************
! !
! !! ROTATION_AXIS2QUAT_3D converts a rotation from axis to quaternion format in 3D.
! !
! !
! !  Definition:
! !
! !    A rotation quaternion Q has the form:
! !
! !      Q = A + Bi + Cj + Dk
! !
! !    where A, B, C and D are real numbers, and i, j, and k are to be regarded
! !    as symbolic constant basis vectors, similar to the role of the "i"
! !    in the representation of aimaginary numbers.
! !
! !    A is the cosine of half of the angle of rotation.  (B,C,D) is a
! !    unit vector pointing in the direction of the axis of rotation.
! !    Rotation multiplication and inversion can be carried out using
! !    this format and the usual rules for quaternion multiplication
! !    and inversion.
! !
! !  Modified:
! !
! !    24 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real AXIS(3), the axis vector which remains unchanged by
! !    the rotation.
! !
! !    Input, real ANGLE, the angular measurement of the rotation about
! !    the axis, in radians.
! !
! !    Output, real Q(4), the quaternion representing the rotation.
! !
!   real axis(3)
!   real angle
!   real norm
!   real q(4)
! !
!   norm = sqrt ( axis(1) * axis(1) + axis(2) * axis(2) + axis(3) * axis(3) )
!
!   if ( norm == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'ROTATION_AXIS2QUAT_3D - Fatal error!'
!     write ( *, * ) '  The axis vector is null.'
!   end if
!
!   q(1) = cos ( 0.5E+00 * angle )
!
!   q(2) = axis(1) * sin ( 0.5E+00 * angle ) / norm
!   q(3) = axis(2) * sin ( 0.5E+00 * angle ) / norm
!   q(4) = axis(3) * sin ( 0.5E+00 * angle ) / norm
!
!   return
! end subroutine
! subroutine rotation_mat_vector_3d ( a, v, w )
! !
! !*******************************************************************************
! !
! !! ROTATION_MAT_VECTOR applies a marix rotation to a vector in 3d.
! !
! !
! !  Modified:
! !
! !    30 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(3,3), the matrix defining the rotation.
! !
! !    Input, real V(3), the vector to be rotated.
! !
! !    Output, real W(3), the rotated vector.
! !
!   real a(3,3)
!   real v(3)
!   real w(3)
! !
!   w(1:3) = matmul ( a(1:3,1:3), v(1:3) )
!
!   return
! end subroutine
! subroutine rotation_mat2axis_3d ( a, axis, angle )
! !
! !*******************************************************************************
! !
! !! ROTATION_MAT2AXIS_3D converts a rotation from matrix to axis format in 3D.
! !
! !
! !  Discussion:
! !
! !    The computation is based on the fact that a rotation matrix must
! !    have an eigenvector corresponding to the eigenvalue of 1, hence:
! !
! !      ( A - I ) * v = 0.
! !
! !    The eigenvector V is the axis of rotation.
! !
! !  Reference:
! !
! !    Jack Kuipers
! !    Quaternions and Rotation Sequences,
! !    Princeton, 1998.
! !
! !  Modified:
! !
! !    27 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(3,3), the rotation matrix.
! !
! !    Output, real AXIS(3), the axis vector which remains unchanged by
! !    the rotation.
! !
! !    Output, real ANGLE, the angular measurement of the rotation about
! !    the axis, in radians.
! !
!   real a(3,3)
!   real axis(3)
!   real angle
!   real norm
! !
!   norm = sqrt ( ( a(3,2) - a(2,3) )**2 + ( a(1,3) - a(3,1) )**2 &
!               + ( a(2,1) - a(1,2) )**2 )
!
!   if ( norm == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'ROTATION_MAT2AXIS_3D - Fatal error!'
!     write ( *, * ) '  A is not a rotation matrix,'
!     write ( *, * ) '  or there are multiple axes of rotation.'
!   end if
!
!   axis(1) = ( a(3,2) - a(2,3) ) / norm
!   axis(2) = ( a(1,3) - a(3,1) ) / norm
!   axis(3) = ( a(2,1) - a(1,2) ) / norm
! !
! !  Find the angle.
! !
!   angle = arc_cosine ( 0.5E+00 * ( a(1,1) + a(2,2) + a(3,3) - 1.0E+00 ) )
!
!   return
! end subroutine
! subroutine rotation_mat2quat_3d ( a, q )
! !
! !*******************************************************************************
! !
! !! ROTATION_MAT2QUAT_3D converts a rotation from matrix to quaternion format in 3D.
! !
! !
! !  Discussion:
! !
! !    The computation is based on the fact that a rotation matrix must
! !    have an eigenvector corresponding to the eigenvalue of 1, hence:
! !
! !      ( A - I ) * v = 0.
! !
! !    The eigenvector V is the axis of rotation.
! !
! !  Reference:
! !
! !    Jack Kuipers
! !    Quaternions and Rotation Sequences,
! !    Princeton, 1998.
! !
! !  Modified:
! !
! !    27 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A(3,3), the rotation matrix.
! !
! !    Output, real Q(4), the quaternion representing the rotation.
! !
!   real a(3,3)
!   real angle
!   real cos_phi
!   real norm
!   real q(4)
!   real sin_phi
! !
!   norm = sqrt ( ( a(3,2) - a(2,3) )**2 + ( a(1,3) - a(3,1) )**2 &
!               + ( a(2,1) - a(1,2) )**2 )
!
!   if ( norm == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'ROTATION_MAT2AXIS_3D - Fatal error!'
!     write ( *, * ) '  A is not a rotation matrix,'
!     write ( *, * ) '  or there are multiple axes of rotation.'
!   end if
!
!   angle = arc_cosine ( 0.5E+00 * ( a(1,1) + a(2,2) + a(3,3) - 1.0E+00 ) )
!
!   cos_phi = cos ( 0.5E+00 * angle )
!
!   sin_phi = sqrt ( 1.0E+00 - cos_phi * cos_phi )
!
!   q(1) = cos_phi
!   q(2) = sin_phi * ( a(3,2) - a(2,3) ) / norm
!   q(3) = sin_phi * ( a(1,3) - a(3,1) ) / norm
!   q(4) = sin_phi * ( a(2,1) - a(1,2) ) / norm
!
!   return
! end subroutine
! subroutine rotation_quat_vector_3d ( q, v, w )
! !
! !*******************************************************************************
! !
! !! ROTATION_QUAT_VECTOR applies a quaternion rotation to a vector in 3d.
! !
! !
! !  Discussion:
! !
! !    If Q is a unit quaternion that encodes a rotation of ANGLE
! !    radians about the vector AXIS, then for an arbitrary real
! !    vector V, the result W of the rotation on V can be written as:
! !
! !      W = Q * V * Conj(Q)
! !
! !  Modified:
! !
! !    29 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real Q(4), the quaternion defining the rotation.
! !
! !    Input, real V(3), the vector to be rotated.
! !
! !    Output, real V(3), the rotated vector.
! !
!   real q(4)
!   real v(3)
!   real w(3)
! !
!   w(1) = &
!          ( 2.0E+00 * ( q(1) * q(1) + q(2) * q(2) ) - 1.0E+00 ) * v(1) &
!        +   2.0E+00 * ( q(2) * q(3) - q(1) * q(4) )             * v(2) &
!        +   2.0E+00 * ( q(2) * q(4) + q(1) * q(3) )             * v(3)
!
!   w(2) = &
!            2.0E+00 * ( q(2) * q(3) + q(1) * q(4) )             * v(1) &
!        + ( 2.0E+00 * ( q(1) * q(1) + q(3) * q(3) ) - 1.0E+00 ) * v(2) &
!        +   2.0E+00 * ( q(3) * q(4) - q(1) * q(2) )             * v(3)
!
!   w(3) = &
!            2.0E+00 * ( q(2) * q(4) - q(1) * q(3) )             * v(1) &
!        +   2.0E+00 * ( q(3) * q(4) + q(1) * q(2) )             * v(2) &
!        + ( 2.0E+00 * ( q(1) * q(1) + q(4) * q(4) ) - 1.0E+00 ) * v(3)
!
!   return
! end subroutine
! subroutine rotation_quat2axis_3d ( q, axis, angle )
! !
! !*******************************************************************************
! !
! !! ROTATION_QUAT2AXIS_3D converts a rotation from quaternion to axis format in 3D.
! !
! !
! !  Definition:
! !
! !    A rotation quaternion Q has the form:
! !
! !      Q = A + Bi + Cj + Dk
! !
! !    where A, B, C and D are real numbers, and i, j, and k are to be regarded
! !    as symbolic constant basis vectors, similar to the role of the "i"
! !    in the representation of aimaginary numbers.
! !
! !    A is the cosine of half of the angle of rotation.  (B,C,D) is a
! !    vector pointing in the direction of the axis of rotation.
! !    Rotation multiplication and inversion can be carried out using
! !    this format and the usual rules for quaternion multiplication
! !    and inversion.
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real Q(4), the quaternion representing the rotation.
! !
! !    Output, real AXIS(3), the axis vector which remains unchanged by
! !    the rotation.
! !
! !    Output, real ANGLE, the angular measurement of the rotation about
! !    the axis, in radians.
! !
!   real axis(3)
!   real angle
!   real cos_phi
!   real q(4)
!   real sin_phi
! !
!   sin_phi = sqrt ( q(2) * q(2) + q(3) * q(3) + q(4) * q(4) )
!
!   cos_phi = q(1)
!
!   angle = 2.0E+00 * atan2 ( sin_phi, cos_phi )
!
!   if ( sin_phi == 0.0E+00 ) then
!     axis(1:3) = (/ 1.0E+00, 0.0E+00, 0.0E+00 /)
!   else
!     axis(1:3) = (/ q(2), q(3), q(4) /) / sin_phi
!   end if
!
!   return
! end subroutine
!
!
! subroutine rotation_quat2mat_3d ( q, a )
! !
! !*******************************************************************************
! !
! !! ROTATION_QUAT2MAT_3D converts a rotation from quaternion to matrix format in 3D.
! !
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Modified:
! !
! !    27 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real AXIS(3), the axis vector which remains unchanged by
! !    the rotation.
! !
! !    Input, real ANGLE, the angular measurement of the rotation about
! !    the axis, in radians.
! !
! !    Output, real A(3,3), the rotation matrix.
! !
!   real a(3,3)
!   real angle
!   real ca
!   real cos_phi
!   real q(4)
!   real sa
!   real sin_phi
!   real v1
!   real v2
!   real v3
! !
!   sin_phi = sqrt ( q(2) * q(2) + q(3) * q(3) + q(4) * q(4) )
!
!   cos_phi = q(1)
!
!   angle = 2.0E+00 * atan2 ( sin_phi, cos_phi )
!
!   if ( sin_phi == 0.0E+00 ) then
!     v1 = 1.0E+00
!     v2 = 0.0E+00
!     v3 = 0.0E+00
!   else
!     v1 = q(2) / sin_phi
!     v2 = q(3) / sin_phi
!     v3 = q(4) / sin_phi
!   end if
!
!   ca = cos ( angle )
!   sa = sin ( angle )
!
!   a(1,1) =                    v1 * v1 + ca * ( 1.0E+00 - v1 * v1 )
!   a(1,2) = ( 1.0E+00 - ca ) * v1 * v2 - sa * v3
!   a(1,3) = ( 1.0E+00 - ca ) * v1 * v3 + sa * v2
!
!   a(2,1) = ( 1.0E+00 - ca ) * v2 * v1 + sa * v3
!   a(2,2) =                    v2 * v2 + ca * ( 1.0E+00 - v2 * v2 )
!   a(2,3) = ( 1.0E+00 - ca ) * v2 * v3 - sa * v1
!
!   a(3,1) = ( 1.0E+00 - ca ) * v3 * v1 - sa * v2
!   a(3,2) = ( 1.0E+00 - ca ) * v3 * v2 + sa * v1
!   a(3,3) =                    v3 * v3 + ca * ( 1.0E+00 - v3 * v3 )
!
!   return
! end subroutine
! subroutine rvec_print_2d ( x, y )
! !
! !*******************************************************************************
! !
! !! RVEC_PRINT_2D prints a 2D vector.
! !
! !
! !  Discussion:
! !
! !    A format is used which suggests a coordinate pair:
! !
! !  Example:
! !
! !    ( 1.23, 7.45 )
! !
! !  Modified:
! !
! !    17 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X, Y, the coordinates of the vector.
! !
!   real x
!   real y
!
!   write ( *, '( a1, g14.6, a1, g14.6, a1 )' ) '(', x, ',', y, ')'
!
!   return
! end subroutine
! subroutine rvec_print_3d ( x, y, z )
! !
! !*******************************************************************************
! !
! !! RVEC_PRINT_3D prints a 3D vector.
! !
! !
! !  Discussion:
! !
! !    A format is used which suggests a coordinate triple:
! !
! !  Example:
! !
! !    ( 1.23, 7.45, -1.45 )
! !
! !  Modified:
! !
! !    17 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X, Y, Z, the coordinates of the vector.
! !
!   real x
!   real y
!   real z
!
!   write ( *, '( a1, g14.6, a1, g14.6, a1, g14.6, a1 )' ) &
!     '(', x, ',', y, ',', z, ')'
!
!   return
! end subroutine
! subroutine rvec2_print ( n, a1, a2, title )
! !
! !*******************************************************************************
! !
! !! RVEC2_PRINT prints a pair of real vectors.
! !
! !
! !  Modified:
! !
! !    14 June 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of components of the vector.
! !
! !    Input, real A1(N), A2(N), the vectors to be printed.
! !
! !    Input, character ( len = * ) TITLE, a title to be printed first.
! !    TITLE may be blank.
! !
!   integer n
! !
!   real a1(n)
!   real a2(n)
!   integer i
!   character ( len = * ) title
! !
!   if ( title /= ' ' ) then
!     write ( *, * ) ' '
!     write ( *, '(a)' ) trim ( title )
!   end if
!
!   write ( *, * ) ' '
!   do i = 1, n
!     write ( *, '(i6,2g14.6)' ) i, a1(i), a2(i)
!   end do
!
!   return
! end subroutine
! subroutine rvec3_print ( n, a1, a2, a3, title )
! !
! !*******************************************************************************
! !
! !! RVEC3_PRINT prints a trio of real vectors.
! !
! !
! !  Modified:
! !
! !    16 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of components of the vector.
! !
! !    Input, real A1(N), A2(N), A3(N), the vectors to be printed.
! !
! !    Input, character ( len = * ) TITLE, a title to be printed first.
! !    TITLE may be blank.
! !
!   integer n
! !
!   real a1(n)
!   real a2(n)
!   real a3(n)
!   integer i
!   character ( len = * ) title
! !
!   if ( title /= ' ' ) then
!     write ( *, * ) ' '
!     write ( *, '(a)' ) trim ( title )
!   end if
!
!   write ( *, * ) ' '
!   do i = 1, n
!     write ( *, '(i6,3g14.6)' ) i, a1(i), a2(i), a3(i)
!   end do
!
!   return
! end subroutine
! subroutine shape_point_dist_2d ( xc, yc, x1, y1, nside, x, y, dist )
! !
! !*******************************************************************************
! !
! !! SHAPE_POINT_DIST_2D: distance ( regular shape, point ) in 2D.
! !
! !
! !  Discussion:
! !
! !    The "regular shape" is assumed to be an equilateral and equiangular
! !    polygon, such as the standard square, pentagon, hexagon, and so on.
! !
! !  Modified:
! !
! !    14 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real XC, YC, the center of the shape.
! !
! !    Input, real X1, Y1, the first vertex of the shape.
! !
! !    Input, integer NSIDE, the number of sides in the shape.
! !
! !    Input, real X, Y, the point to be checked.
! !
! !    Output, real DIST, the distance from the point to the shape.
! !
!   real angle
!   real angle2
!   real dist
!   real dist2
!   integer nside
!   real radius
!   real sector_angle
!   integer sector_index
!   real x
!   real x1
!   real xa
!   real xb
!   real xc
!   real y
!   real y1
!   real ya
!   real yb
!   real yc
! !
! !  Determine the angle subtended by a single side.
! !
!   sector_angle = 360.0E+00 / real ( nside )
! !
! !  How long is the half-diagonal?
! !
!   radius = enorm0_2d ( x1, y1, xc, yc )
! !
! !  If the radius is zero, then the shape is a point and the computation is easy.
! !
!   if ( radius == 0.0E+00 ) then
!     dist = enorm0_2d ( x, y, xc, yc )
!     return
!   end if
! !
! !  If the test point is at the center, the computation is easy.
! !
!   dist2 = enorm0_2d ( x, y, xc, yc )
!
!   if ( dist2 == 0.0E+00 ) then
!     dist = radius
!     return
!   end if
! !
! !  Determine the angle between the ray to the first corner,
! !  and the ray to the test point.
! !
!   angle = angle_deg_2d ( x1, y1, xc, yc, x, y )
! !
! !  Determine the sector of the point.
! !
!   sector_index = int ( angle / sector_angle ) + 1
! !
! !  Generate the two corner points that terminate the SECTOR-th side.
! !
!   angle2 = real ( sector_index - 1 ) * sector_angle
!   angle2 = degrees_to_radians ( angle2 )
!
!   call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xa, ya )
!
!   angle2 = real ( sector_index ) * sector_angle
!   angle2 = degrees_to_radians ( angle2 )
!
!   call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xb, yb )
! !
! !  Determine the distance from the test point to the line segment that
! !  is the SECTOR-th side.
! !
!   call line_seg_point_dist_2d ( xa, ya, xb, yb, x, y, dist )
!
!   return
! end subroutine
! subroutine shape_point_near_2d ( xc, yc, x1, y1, nside, x, y, xn, yn, dist )
! !
! !*******************************************************************************
! !
! !! SHAPE_POINT_NEAR_2D: nearest point ( regular shape, point ) in 2D.
! !
! !
! !  Discussion:
! !
! !    The "regular shape" is assumed to be an equilateral and equiangular
! !    polygon, such as the standard square, pentagon, hexagon, and so on.
! !
! !  Modified:
! !
! !    17 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real XC, YC, the center of the shape.
! !
! !    Input, real X1, Y1, the first vertex of the shape.
! !
! !    Input, integer NSIDE, the number of sides in the shape.
! !
! !    Input, real X, Y, the point to be checked.
! !
! !    Output, real XN, YN, the point on the shape that is nearest
! !    to the given point.
! !
! !    Output, real DIST, the distance between the points.
! !
!   real angle
!   real angle2
!   real dist
!   real dist2
!   integer nside
!   real radius
!   real sector_angle
!   integer sector_index
!   real t
!   real x
!   real x1
!   real xa
!   real xb
!   real xc
!   real xn
!   real y
!   real y1
!   real ya
!   real yb
!   real yc
!   real yn
! !
! !  Determine the angle subtended by a single side.
! !
!   sector_angle = 360.0E+00 / real ( nside )
! !
! !  How long is the half-diagonal?
! !
!   radius = enorm0_2d ( x1, y1, xc, yc )
! !
! !  If the radius is zero, then the shape is a point and the computation is easy.
! !
!   if ( radius == 0.0E+00 ) then
!     dist = enorm0_2d ( x, y, xc, yc )
!     return
!   end if
! !
! !  If the test point is at the center, the computation is easy.
! !
!   dist2 = enorm0_2d ( x, y, xc, yc )
!
!   if ( dist2 == 0.0E+00 ) then
!     dist = radius
!     return
!   end if
! !
! !  Determine the angle between the ray to the first corner,
! !  and the ray to the test point.
! !
!   angle = angle_deg_2d ( x1, y1, xc, yc, x, y )
! !
! !  Determine the sector of the point.
! !
!   sector_index = int ( angle / sector_angle ) + 1
! !
! !  Generate the two corner points that terminate the SECTOR-th side.
! !
!   angle2 = real ( sector_index - 1 ) * sector_angle
!   angle2 = degrees_to_radians ( angle2 )
!
!   call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xa, ya )
!
!   angle2 = real ( sector_index ) * sector_angle
!   angle2 = degrees_to_radians ( angle2 )
!
!   call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xb, yb )
! !
! !  Determine the point on the SECTOR-th side of the shape which is
! !  nearest.
! !
!   call line_seg_point_near_2d ( xa, ya, xb, yb, x, y, xn, yn, dist, t )
!
!   return
! end subroutine
! subroutine shape_print_3d ( max_num, max_order, point_num, &
!    face_num, face_order, point_coord, face_point )
! !
! !*******************************************************************************
! !
! !! SHAPE_PRINT_3D prints information about a polyhedron in 3D.
! !
! !
! !  Modified:
! !
! !    11 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAX_NUM, the maximum number of faces, and of
! !    points.  This quantity is used to dimension arrays.
! !
! !    Input, integer MAX_ORDER, the maximum order of any face.
! !
! !    Input, integer POINT_NUM, the number of points in the shape.
! !
! !    Input, integer FACE_NUM, the number of faces in the shape.
! !
! !    Input, integer FACE_ORDER, the number of vertices per face.
! !
! !    Input, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
! !    the X, Y and Z coordinates of point J.
! !
! !    Input, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
! !    contains the index of the I-th point in the J-th face.  The
! !    points are listed in the counter-clockwise direction defined
! !    by the outward normal at the face.
! !
!   integer max_num
!   integer max_order
! !
!   integer face_num
!   integer face_order
!   integer face_point(max_order,max_num)
!   integer i
!   integer j
!   integer point_num
!   real point_coord(3,max_num)
! !
!   write ( *, * ) ' '
!   write ( *, * ) 'SHAPE_PRINT_3D'
!   write ( *, * ) '  Information about a regular polyhedron.'
!   write ( *, * ) ' '
!   write ( *, * ) '  The number of vertices is ', point_num
!   write ( *, * ) ' '
!   write ( *, * ) '  Vertex coordinates are:'
!   write ( *, * ) ' '
!   do i = 1, point_num
!     write ( *, '(i4,2x,3f10.4)' ) i, point_coord(1:3,i)
!   end do
!
!   write ( *, * ) ' '
!   write ( *, * ) '  The number of faces is ', face_num
!   write ( *, * ) '  The order of each face is ', face_order
!   write ( *, * ) ' '
!   write ( *, * ) '  The vertices in each face are:'
!   write ( *, * ) ' '
!
!   do i = 1, face_num
!     write ( *, '(i5,3x,10i5)' ) i, face_point(1:face_order,i)
!   end do
!
!   return
! end subroutine
! subroutine shape_ray_int_2d ( xc, yc, x1, y1, nside, xa, ya, xb, yb, xi, yi )
! !
! !*******************************************************************************
! !
! !! SHAPE_RAY_INT_2D: intersection ( regular shape, ray ) in 2D.
! !
! !
! !  Discussion:
! !
! !    The "regular shape" is assumed to be an equilateral and equiangular
! !    polygon, such as the standard square, pentagon, hexagon, and so on.
! !
! !    The origin of the ray is assumed to be inside the shape.  This
! !    guarantees that the ray will intersect the shape in exactly one point.
! !
! !  Modified:
! !
! !    17 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real XC, YC, the center of the shape.
! !
! !    Input, real X1, Y1, the first vertex of the shape.
! !
! !    Input, integer NSIDE, the number of sides in the shape.
! !
! !    Input, real XA, YA, the origin of the ray.
! !
! !    Input, real XB, YB, a second point on the ray.
! !
! !    Output, real XI, YI, the point on the shape intersected by the ray.
! !
!   real angle2
!   logical inside
!   integer ival
!   integer nside
!   real radius
!   real sector_angle
!   integer sector_index
!   real x1
!   real xa
!   real xb
!   real xc
!   real xi
!   real xv1
!   real xv2
!   real y1
!   real ya
!   real yb
!   real yc
!   real yi
!   real yv1
!   real yv2
! !
! !  Warning!
! !  No check is made to ensure that the ray origin is inside the shape.
! !  These calculations are not valid if that is not true!
! !
! !  Determine the angle subtended by a single side.
! !
!   sector_angle = 360.0E+00 / real ( nside )
! !
! !  How long is the half-diagonal?
! !
!   radius = enorm0_2d ( x1, y1, xc, yc )
! !
! !  If the radius is zero, refuse to continue.
! !
!   if ( radius == 0.0E+00 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'SHAPE_RAY_INT_2D - Fatal error!'
!     write ( *, * ) '  The shape has radius zero.'
!     stop
!   end if
! !
! !  Determine which sector side intersects the ray.
! !
!   xv2 = 0.0E+00
!   yv2 = 0.0E+00
!
!   do sector_index = 1, nside
! !
! !  Determine the two vertices that define this sector.
! !
!     if ( sector_index == 1 ) then
!
!       angle2 = real ( sector_index - 1 ) * sector_angle
!       angle2 = degrees_to_radians ( angle2 )
!
!       call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xv1, yv1 )
!
!     else
!
!       xv1 = xv2
!       yv1 = yv2
!
!     end if
!
!     angle2 = real ( sector_index ) * sector_angle
!     angle2 = degrees_to_radians ( angle2 )
!
!     call vector_rotate_base_2d ( x1, y1, xc, yc, angle2, xv2, yv2 )
! !
! !  Draw the angle from one vertex to the ray origin to the next vertex,
! !  and see if that angle contains the ray.  If so, then the ray
! !  must intersect the shape side of that sector.
! !
!     call angle_contains_ray_2d ( inside, xv1, yv1, xa, ya, xv2, yv2, xb, yb )
!
!     if ( inside ) then
! !
! !  Determine the intersection of the lines defined by the ray and the
! !  sector side.  (We're already convinced that the ray and sector line
! !  segment intersect, so we can use the simpler code that treats them
! !  as full lines).
! !
!       call lines_exp_int_2d ( ival, xi, yi, xa, ya, xb, yb, xv1, yv1, xv2, yv2 )
!
!       return
!
!     end if
!
!   end do
! !
! !  If the calculation fell through the loop, then something's wrong.
! !
!   write ( *, * ) ' '
!   write ( *, * ) 'SHAPE_RAY_INT_2D - Fatal error!'
!   write ( *, * ) '  Cannot find intersection of ray and shape.'
!   stop
! end subroutine
! subroutine sort_heap_external ( n, indx, i, j, isgn )
! !
! !*******************************************************************************
! !
! !! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
! !
! !
! !  Discussion:
! !
! !    The actual list of data is not passed to the routine.  Hence this
! !    routine may be used to sort integers, reals, numbers, names,
! !    dates, shoe sizes, and so on.  After each call, the routine asks
! !    the user to compare or interchange two items, until a special
! !    return value signals that the sorting is completed.
! !
! !  Modified:
! !
! !    19 May 1999
! !
! !  Reference:
! !
! !    A Nijenhuis and H Wilf,
! !    Combinatorial Algorithms,
! !    Academic Press, 1978, second edition,
! !    ISBN 0-12-519260-6.
! !
! !  Parameters:
! !
! !    Input, integer N, the number of items to be sorted.
! !
! !    Input/output, integer INDX, the main communication signal.
! !
! !    The user must set INDX to 0 before the first call.
! !    Thereafter, the user should not change the value of INDX until
! !    the sorting is done.
! !
! !    On return, if INDX is
! !
! !      greater than 0,
! !      * interchange items I and J;
! !      * call again.
! !
! !      less than 0,
! !      * compare items I and J;
! !      * set ISGN = -1 if I < J, ISGN = +1 if I > J;
! !      * call again.
! !
! !      equal to 0, the sorting is done.
! !
! !    Output, integer I, J, the indices of two items.
! !    On return with INDX positive, elements I and J should be interchanged.
! !    On return with INDX negative, elements I and J should be compared, and
! !    the result reported in ISGN on the next call.
! !
! !    Input, integer ISGN, results of comparison of elements I and J.
! !    (Used only when the previous call returned INDX less than 0).
! !    ISGN <= 0 means I is less than or equal to J;
! !    ISGN => 0 means I is greater than or equal to J.
! !
!   integer i
!   integer indx
!   integer isgn
!   integer j
!   integer, save :: k = 0
!   integer, save :: k1 = 0
!   integer n
!   integer, save :: n1 = 0
! !
! !  INDX = 0: This is the first call.
! !
!   if ( indx == 0 ) then
!
!     n1 = n
!     k = n / 2
!     k1 = k
! !
! !  INDX < 0: The user is returning the results of a comparison.
! !
!   else if ( indx < 0 ) then
!
!     if ( indx == -2 ) then
!
!       if ( isgn < 0 ) then
!         i = i + 1
!       end if
!
!       j = k1
!       k1 = i
!       indx = - 1
!       return
!
!     end if
!
!     if ( isgn > 0 ) then
!       indx = 2
!       return
!     end if
!
!     if ( k <= 1 ) then
!
!       if ( n1 == 1 ) then
!         indx = 0
!       else
!         i = n1
!         n1 = n1 - 1
!         j = 1
!         indx = 1
!       end if
!
!       return
!
!     end if
!
!     k = k - 1
!     k1 = k
! !
! !  INDX > 0, the user was asked to make an interchange.
! !
!   else if ( indx == 1 ) then
!
!     k1 = k
!
!   end if
!
!   do
!
!     i = 2 * k1
!
!     if ( i == n1 ) then
!       j = k1
!       k1 = i
!       indx = - 1
!       return
!     else if ( i <= n1 ) then
!       j = i + 1
!       indx = - 2
!       return
!     end if
!
!     if ( k <= 1 ) then
!       exit
!     end if
!
!     k = k - 1
!     k1 = k
!
!   end do
!
!   if ( n1 == 1 ) then
!     indx = 0
!   else
!     i = n1
!     n1 = n1 - 1
!     j = 1
!     indx = 1
!   end if
!
!   return
! end subroutine
! subroutine sphere_dia2imp_3d ( x1, y1, z1, x2, y2, z2, r, xc, yc, zc )
! !
! !*******************************************************************************
! !
! !! SPHERE_DIA2IMP_3D converts a diameter to an implicit sphere in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, are the (X,Y,Z) coordinates
! !    of two points which form a diameter of the sphere.
! !
! !    Output, real R, the computed radius of the sphere.
! !
! !    Output, real XC, YC, ZC, the computed center of the sphere.
! !
!   real r
!   real x1
!   real x2
!   real xc
!   real y1
!   real y2
!   real yc
!   real z1
!   real z2
!   real zc
! !
!   r = 0.5E+00 * enorm0_3d ( x1, y1, z1, x2, y2, z2 )
!
!   xc = 0.5E+00 * ( x1 + x2 )
!   yc = 0.5E+00 * ( y1 + y2 )
!   zc = 0.5E+00 * ( z1 + z2 )
!
!   return
! end subroutine
! subroutine sphere_exp_contains_point_3d ( x1, y1, z1, x2, y2, z2, &
!   x3, y3, z3, x4, y4, z4, x, y, z, inside )
! !
! !*******************************************************************************
! !
! !! SPHERE_EXP_CONTAINS_POINT_3D determines if an explicit sphere contains a point in 3D.
! !
! !
! !  Formula:
! !
! !    An explicit sphere in 3D is determined by four points,
! !    which should be distinct, and not coplanar.
! !
! !  Method:
! !
! !    The computation checks the determinant of:
! !
! !      x1  y1  z1  x1**2+y1**2+z1**2  1
! !      x2  y2  z2  x2**2+y2**2+z2**2  1
! !      x3  y3  z3  x3**2+y3**2+z3**2  1
! !      x4  y4  z4  x4**2+y4**2+z4**2  1
! !      x   y   z   x**2 +y**2 +z**2   1
! !
! !  Modified:
! !
! !    26 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
! !    (X,Y,Z) coordinates of four points that lie on a circle.
! !
! !    Input, real X, Y, Z, the (X,Y,Z) coordinates of a point, whose
! !    position relative to the sphere is desired.
! !
! !    Output, logical INSIDE, is TRUE if the point is in the sphere,
! !    FALSE otherwise.
! !
!   real a(5,5)
!   real det
!   logical inside
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
!   real z
!   real z1
!   real z2
!   real z3
!   real z4
! !
! !  Compute the determinant.
! !
!   a(1,1) = x1
!   a(1,2) = y1
!   a(1,3) = z1
!   a(1,4) = x1 * x1 + y1 * y1 + z1 * z1
!   a(1,5) = 1.0E+00
!
!   a(2,1) = x2
!   a(2,2) = y2
!   a(2,3) = z2
!   a(2,4) = x2 * x2 + y2 * y2 + z2 * z2
!   a(2,5) = 1.0E+00
!
!   a(3,1) = x3
!   a(3,2) = y3
!   a(3,3) = z3
!   a(3,4) = x3 * x3 + y3 * y3 + z3 * z3
!   a(3,5) = 1.0E+00
!
!   a(4,1) = x4
!   a(4,2) = y4
!   a(4,3) = z4
!   a(4,4) = x4 * x4 + y4 * y4 + z4 * z4
!   a(4,5) = 1.0E+00
!
!   a(5,1) = x
!   a(5,2) = y
!   a(5,3) = z
!   a(5,4) = x * x + y * y + z * z
!   a(5,5) = 1.0E+00
!
!   det = rmat5_det ( a )
!
!   if ( det < 0.0E+00 ) then
!     inside = .false.
!   else if ( det >= 0.0E+00 ) then
!     inside = .true.
!   end if
!
!   return
! end subroutine
! subroutine sphere_exp_near_point_3d ( x1, y1, z1, x2, y2, z2, &
!   x3, y3, z3, x4, y4, z4, x, y, z, xn, yn, zn )
! !
! !*******************************************************************************
! !
! !! SPHERE_EXP_NEAR_POINT_3D finds the nearest point on an explicit sphere to a point in 3D.
! !
! !
! !  Formula:
! !
! !    An explicit sphere in 3D is determined by four points,
! !    which should be distinct, and not coplanar.
! !
! !  Method:
! !
! !    If the center of the sphere is (Xc,Yc,Zc), and the point is (X,Y,Z), then
! !    the desired point lies at a positive distance R along the vector (X-Xc, Y-Yc, Z-Zc)
! !    unless (X,Y,Z) = (Xc,Yc,Zc), in which case any point on the sphere is "nearest".
! !
! !  Modified:
! !
! !    14 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
! !    (X,Y,Z) coordinates of four points that lie on a circle.
! !
! !    Input, real X, Y, Z, the (X,Y,Z) coordinates of a point, whose
! !    nearest point on the sphere is desired.
! !
! !    Output, real XN, YN, ZN, the nearest point on the sphere.
! !
!   real norm
!   real r
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real xc
!   real xn
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
!   real yc
!   real yn
!   real z
!   real z1
!   real z2
!   real z3
!   real z4
!   real zc
!   real zn
! !
! !  Find the center.
! !
!   call sphere_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, &
!     z3, x4, y4, z4, r, xc, yc, zc )
! !
! !  If (X,Y,Z) = (XC,YC,ZC), bail out now.
! !
!   norm = enorm0_3d ( x, y, z, xc, yc, zc )
!
!   if ( norm == 0.0E+00 ) then
!     xn = xc + r
!     yn = yc
!     zn = zc
!     return
!   end if
! !
! !  Compute the nearest point.
! !
!   xn = xc + r * ( x - xc ) / norm
!   yn = yc + r * ( y - yc ) / norm
!   zn = zc + r * ( z - zc ) / norm
!
!   return
! end subroutine
! subroutine sphere_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, &
!   z3, x4, y4, z4, r, xc, yc, zc )
! !
! !*******************************************************************************
! !
! !! SPHERE_EXP2IMP_3D converts a sphere from explicit to implicit form in 3D.
! !
! !
! !  Formula:
! !
! !    An explicit sphere in 3D is determined by four points,
! !    which should be distinct, and not coplanar.
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    05 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4,
! !    the coordinates of four distinct noncoplanar points on the sphere.
! !
! !    Output, real R, XC, YC, ZC, the radius and coordinates of the center
! !    of the sphere.  If the linear system is
! !    singular, then R = -1, XC = YC = ZC = 0.
! !
!   real r
!   real x1
!   real x2
!   real x3
!   real x4
!   real xc
!   real y1
!   real y2
!   real y3
!   real y4
!   real yc
!   real z1
!   real z2
!   real z3
!   real z4
!   real zc
! !
!   call tetra_outsphere_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
!     r, xc, yc, zc )
!
!   return
! end subroutine
! subroutine sphere_imp_area_3d ( r, area )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_AREA_3D computes the surface area of an implicit sphere in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Modified:
! !
! !    26 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Output, real AREA, the area of the sphere.
! !
!   real area
!   real r
! !
!   area = 4.0E+00 * pi * r * r
!
!   return
! end subroutine
! subroutine sphere_imp_area_nd ( n, r, area )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_AREA_ND computes the surface area of an implicit sphere in ND.
! !
! !
! !  Discussion:
! !
! !    N   Area
! !
! !    2   2       * PI    * R
! !    3   4       * PI    * R**2
! !    4   2       * PI**2 * R**3
! !    5   (8/3)   * PI**2 * R**4
! !    6             PI**3 * R**5
! !    7   (16/15) * PI**3 * R**6
! !
! !  Modified:
! !
! !    26 October 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Output, real AREA, the area of the sphere.
! !
!   real area
!   integer n
!   real r
! !
!   call sphere_imp_unit_area_nd ( n, area )
!
!   area = area * r**(n-1)
!
!   return
! end subroutine
! subroutine sphere_imp_contains_point_3d ( r, xc, yc, zc, x, y, z, inside )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_CONTAINS_POINT_3D determines if an implicit sphere contains a point in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Modified:
! !
! !    05 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
! !
! !    Input, real X, Y, Z, the point to be checked.
! !
! !    Output, logical INSIDE, is TRUE if the point is inside or on the sphere,
! !    FALSE otherwise.
! !
!   logical inside
!   real r
!   real x
!   real xc
!   real y
!   real yc
!   real z
!   real zc
! !
!   if ( enormsq0_3d ( x, y, z, xc, yc, zc ) <= r * r ) then
!     inside = .true.
!   else
!     inside = .false.
!   end if
!
!   return
! end subroutine
! subroutine sphere_imp_gridfaces_3d ( maxtri, nlat, nlong, ntri, tri )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_GRIDFACES_3D produces a grid of triangles on an implicit sphere in 3D.
! !
! !
! !  Discussion:
! !
! !    The point numbering system is the same used in SPHERE_IMP_GRIDPOINTS_3D,
! !    and that routine may be used to compute the coordinates of the points.
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Modified:
! !
! !    25 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAXTRI, the maximum number of triangles.
! !
! !    Input, integer NLAT, NLONG, the number of latitude and longitude
! !    lines to draw.  The latitudes do not include the North and South
! !    poles, which will be included automatically, so NLAT = 5, for instance,
! !    will result in points along 7 lines of latitude.
! !
! !    Output, integer NTRI, the number of triangles.
! !
! !    Output, integer TRI(3,MAXTRI), contains triples of point indices for
! !    the triangles that make up the grid.
! !
!   integer maxtri
! !
!   integer i
!   integer j
!   integer n
!   integer n_max
!   integer n_min
!   integer ne
!   integer nlat
!   integer nlong
!   integer ntri
!   integer nw
!   integer s
!   integer s_max
!   integer s_min
!   integer se
!   integer sw
!   integer tri(3,maxtri)
! !
!   ntri = 0
! !
! !  The first row.
! !
!   n = 1
!
!   sw = 2
!   se = sw + 1
!
!   s_min = 2
!   s_max = nlong + 1
!
!   do j = 0, nlong - 1
!
!     if ( ntri < maxtri ) then
!       ntri = ntri + 1
!       tri(1,ntri) = sw
!       tri(2,ntri) = se
!       tri(3,ntri) = n
!     end if
!
!     sw = se
!
!     if ( se == s_max ) then
!       se = s_min
!     else
!       se = se + 1
!     end if
!
!   end do
! !
! !  The intermediate rows.
! !
!   do i = 1, nlat
!
!     n_max = s_max
!     n_min = s_min
!
!     s_max = s_max + nlong
!     s_min = s_min + nlong
!
!     nw = n_min
!     ne = nw + 1
!     sw = s_min
!     se = sw + 1
!
!     do j = 0, nlong - 1
!
!       if ( ntri < maxtri ) then
!         ntri = ntri + 1
!         tri(1,ntri) = sw
!         tri(2,ntri) = se
!         tri(3,ntri) = nw
!       end if
!
!       if ( ntri < maxtri ) then
!         ntri = ntri + 1
!         tri(1,ntri) = ne
!         tri(2,ntri) = nw
!         tri(3,ntri) = se
!       end if
!
!       sw = se
!       nw = ne
!
!       if ( se == s_max ) then
!         se = s_min
!       else
!         se = se + 1
!       end if
!
!       if ( ne == n_max ) then
!         ne = n_min
!       else
!         ne = ne + 1
!       end if
!
!     end do
!
!   end do
! !
! !  The last row.
! !
!   n_max = s_max
!   n_min = s_min
!
!   s = n_max + 1
!
!   nw = n_min
!   ne = nw + 1
!
!   do j = 0, nlong - 1
!
!     if ( ntri < maxtri ) then
!       ntri = ntri + 1
!       tri(1,ntri) = ne
!       tri(2,ntri) = nw
!       tri(3,ntri) = s
!     end if
!
!     nw = ne
!
!     if ( ne == n_max ) then
!       ne = n_min
!     else
!       ne = ne + 1
!     end if
!
!   end do
!   return
! end subroutine
! subroutine sphere_imp_gridlines_3d ( maxline, nlat, nlong, nline, line )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_GRIDLINES_3D produces "grid lines" on an implicit sphere in 3D.
! !
! !
! !  Discussion:
! !
! !    The point numbering system is the same used in SPHERE_IMP_GRIDPOINTS_3D,
! !    and that routine may be used to compute the coordinates of the points.
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Modified:
! !
! !    24 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAXLINE, the maximum number of gridlines.
! !
! !    Input, integer NLAT, NLONG, the number of latitude and longitude
! !    lines to draw.  The latitudes do not include the North and South
! !    poles, which will be included automatically, so NLAT = 5, for instance,
! !    will result in points along 7 lines of latitude.
! !
! !    Output, integer NLINE, the number of grid lines.
! !
! !    Output, integer LINE(2,MAXLINE), contains pairs of point indices for
! !    line segments that make up the grid.
! !
!   integer maxline
! !
!   integer i
!   integer j
!   integer line(2,maxline)
!   integer new
!   integer newcol
!   integer nlat
!   integer nline
!   integer nlong
!   integer old
! !
!   nline = 0
! !
! !  "Vertical" lines.
! !
!   do j = 0, nlong - 1
!
!     old = 1
!     new = j + 2
!
!     if ( nline < maxline ) then
!       nline = nline + 1
!       line(1,nline) = old
!       line(2,nline) = new
!     end if
!
!     do i = 1, nlat - 1
!
!       old = new
!       new = old + nlong
!
!       if ( nline < maxline ) then
!         nline = nline + 1
!         line(1,nline) = old
!         line(2,nline) = new
!       end if
!
!     end do
!
!     old = new
!
!     if ( nline < maxline ) then
!       nline = nline + 1
!       line(1,nline) = old
!       line(2,nline) = 1 + nlat * nlong + 1
!     end if
!
!   end do
! !
! !  "Horizontal" lines.
! !
!   do i = 1, nlat
!
!     new = 1 + ( i - 1 ) * nlong + 1
!
!     do j = 0, nlong - 2
!       old = new
!       new = old + 1
!       if ( nline < maxline ) then
!         nline = nline + 1
!         line(1,nline) = old
!         line(2,nline) = new
!       end if
!     end do
!
!     old = new
!     new = 1 + ( i - 1 ) * nlong + 1
!     if ( nline < maxline ) then
!       line(1,nline) = old
!       line(2,nline) = new
!     end if
!
!   end do
! !
! !  "Diagonal" lines.
! !
!   do j = 0, nlong - 1
!
!     old = 1
!     new = j + 2
!     newcol = j
!
!     do i = 1, nlat - 1
!
!       old = new
!       new = old + nlong + 1
!
!       newcol = newcol + 1
!       if ( newcol > nlong - 1 ) then
!         newcol = 0
!         new = new - nlong
!       end if
!
!       if ( nline < maxline ) then
!         nline = nline + 1
!         line(1,nline) = old
!         line(2,nline) = new
!       end if
!
!     end do
!
!   end do
!
!   return
! end subroutine
! subroutine sphere_imp_gridpoints_3d ( r, xc, yc, zc, maxpoint, nlat, nlong, &
!   npoint, x, y, z )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_GRIDPOINTS_3D produces "grid points" on an implicit sphere in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Modified:
! !
! !    22 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
! !
! !    Input, integer MAXPOINT, the maximum number of grid points, which
! !    should be at least 2 + NLAT * NLONG.
! !
! !    Input, integer NLAT, NLONG, the number of latitude and longitude
! !    lines to draw.  The latitudes do not include the North and South
! !    poles, which will be included automatically, so NLAT = 5, for instance,
! !    will result in points along 7 lines of latitude.
! !
! !    Output, integer NPOINT, the number of grid points.  The number of
! !    grid points depends on N as follows:
! !
! !      NPOINT = 2 + NLAT * NLONG.
! !
! !    Output, real X(MAXPOINT), Y(MAXPOINT), Z(MAXPOINT), the coordinates of the
! !    grid points.
! !
!   integer maxpoint
! !
!   integer i
!   integer j
!   integer nlat
!   integer nlong
!   integer npoint
!   real pi_const
!   real phi
!   real r
!   real theta
!   real x(maxpoint)
!   real xc
!   real y(maxpoint)
!   real yc
!   real z(maxpoint)
!   real zc
! !
!   npoint = 0
!   pi_const = pi
! !
! !  The north pole.
! !
!   theta = 0.0E+00
!   phi = 0.0E+00
!   npoint = npoint + 1
!   if ( npoint <= maxpoint ) then
!     call sphere_imp_local2xyz_3d ( r, xc, yc, zc, theta, phi, &
!       x(npoint), y(npoint), z(npoint) )
!   end if
! !
! !  Do each intermediate ring of latitude.
! !
!   do i = 1, nlat
!
!     phi = pi_const * real ( i ) / real ( nlat + 1 )
! !
! !  Along that ring of latitude, compute points at various longitudes.
! !
!     do j = 0, nlong-1
!
!       theta = 2.0E+00 * pi_const * real ( j ) / real ( nlong )
!
!       npoint = npoint + 1
!       if ( npoint <= maxpoint ) then
!         call sphere_imp_local2xyz_3d ( r, xc, yc, zc, theta, phi, &
!           x(npoint), y(npoint), z(npoint) )
!       end if
!
!     end do
!   end do
! !
! !  The south pole.
! !
!   theta = 0.0E+00
!   phi = pi_const
!   npoint = npoint + 1
!   if ( npoint <= maxpoint ) then
!     call sphere_imp_local2xyz_3d ( r, xc, yc, zc, theta, phi, &
!       x(npoint), y(npoint), z(npoint) )
!   end if
!
!   return
! end subroutine
! subroutine sphere_imp_line_project_3d ( r, xc, yc, zc, npnt, &
!   x, y, z, maxpnt2, npnt2, xp, yp, zp, thetamin, thetamax )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_LINE_PROJECT_3D projects a line onto an implicit sphere in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Discussion:
! !
! !    The line to be projected is specified as a sequence of points.
! !    If two successive points subtend a small angle, then the second
! !    point is essentially dropped.  If two successive points subtend
! !    a large angle, then intermediate points are inserted, so that
! !    the projected line stays closer to the sphere.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the sphere.  If R is
! !    zero, (XP,YP,ZP) will be returned as (XC,YC,ZC), and if R is
! !    negative, points will end up diametrically opposite from where
! !    you would expect them for a positive R.
! !
! !    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
! !
! !    Input, integer NPNT, the number of points on the line that is
! !    to be projected.
! !
! !    Input, real X(NPNT), Y(NPNT), Z(NPNT), the coordinates of the points
! !    on the line that is to be projected.
! !    Note: if any (X,Y,Z) coincides with the center of the sphere, then
! !    its projection is mathematically undefined.  (XP,YP,ZP) will
! !    be returned as (XC,YC,ZC).
! !
! !    Input, integer MAXPNT2, the maximum number of points on the projected
! !    line.  Even if the routine thinks that more points are needed,
! !    no more than MAXPNT2 will be generated.
! !
! !    Output, integer NPNT2, the number of points on the projected line.
! !    NPNT2 can be zero, if the line has an angular projection of less
! !    than THETAMIN radians.
! !
! !    Output, real XP(NPNT2), YP(NPNT2), ZP(NPNT2), the coordinates of the
! !    points representing the projected line.  These points lie on the
! !    sphere.  Successive points are separated by at least THETAMIN
! !    radians, and by no more than THETAMAX radians.
! !
! !    Input, real THETAMIN, THETAMAX, the minimum and maximum angular
! !    projections allowed between successive projected points.
! !    If two successive points on the original line have projections
! !    separated by more than THETAMAX radians, then intermediate points
! !    will be inserted, in an attempt to keep the line closer to the
! !    sphere.  If two successive points are separated by less than
! !    THETAMIN radians, then the second point is dropped, and the
! !    line from the first point to the next point is considered.
! !
!   integer maxpnt2
!   integer npnt
! !
!   real alpha
!   real ang3d
!   integer i
!   integer j
!   integer nfill
!   integer npnt2
!   real r
!   real thetamax
!   real thetamin
!   real tnorm
!   real x(npnt)
!   real x1
!   real x2
!   real xc
!   real xi
!   real xp(maxpnt2)
!   real y(npnt)
!   real y1
!   real y2
!   real yc
!   real yi
!   real yp(maxpnt2)
!   real z(npnt)
!   real z1
!   real z2
!   real zc
!   real zi
!   real zp(maxpnt2)
!   real dot
!
! !
! !  Check the input.
! !
!   if ( r == 0.0E+00 ) then
!     npnt2 = 0
!     return
!   end if
!
!   x1 = xc
!   y1 = yc
!   z1 = zc
!
!   x2 = xc
!   y2 = yc
!   z2 = zc
!
!   npnt2 = 0
!
!   do i = 1, npnt
!
!     if ( x(i) == xc .and. y(i) == yc .and. z(i) == zc ) then
!
!     else
!
!       x1 = x2
!       y1 = y2
!       z1 = z2
!
!       alpha = enorm0_3d ( xc, yc, zc, x(i), y(i), z(i) )
!
!       x2 = xc + r * ( x(i) - xc ) / alpha
!       y2 = yc + r * ( y(i) - yc ) / alpha
!       z2 = zc + r * ( z(i) - zc ) / alpha
! !
! !  If we haven't gotten any points yet, take this point as our start.
! !
!       if ( npnt2 == 0 ) then
!
!         npnt2 = npnt2 + 1
!         xp(npnt2) = x2
!         yp(npnt2) = y2
!         zp(npnt2) = z2
! !
! !  Compute the angular projection of (X1,Y1,Z1) to (X2,Y2,Z2).
! !
!     else if ( npnt2 >= 1 ) then
!
!       dot = ( x1 - xc ) * ( x2 - xc ) + ( y1 - yc ) * ( y2 - yc ) + &
!             ( z1 - zc ) * ( z2 - zc )
!       ang3d = arc_cosine (  dot / r**2 )
! !
! !  If the angle is at least THETAMIN, (or it's the last point),
! !  then we will draw a line segment.
! !
!       if ( abs ( ang3d ) > thetamin .or. i == npnt ) then
! !
! !  Now we check to see if the line segment is too long.
! !
!         if ( abs ( ang3d ) > thetamax ) then
!
!           nfill = int ( abs ( ang3d ) / thetamax )
!
!           do j = 1, nfill-1
!
!             xi = ( ( nfill - j ) * ( x1 - xc ) + j * ( x2 - xc ) )
!             yi = ( ( nfill - j ) * ( y1 - yc ) + j * ( y2 - yc ) )
!             zi = ( ( nfill - j ) * ( z1 - zc ) + j * ( z2 - zc ) )
!             tnorm = sqrt ( xi * xi + yi * yi + zi * zi )
!
!             if ( tnorm /= 0.0E+00 ) then
!               xi = xc + r * xi / tnorm
!               yi = yc + r * yi / tnorm
!               zi = zc + r * zi / tnorm
!               npnt2 = npnt2 + 1
!               xp(npnt2) = xi
!               yp(npnt2) = yi
!               zp(npnt2) = zi
!             end if
!
!           end do
!
!         end if
! !
! !  Now tack on the projection of point 2.
! !
!         npnt2 = npnt2 + 1
!         xp(npnt2) = x2
!         yp(npnt2) = y2
!         zp(npnt2) = z2
!
!       end if
!
!     end if
!
!     end if
!
!   end do
!
!   return
! end subroutine
! subroutine sphere_imp_local2xyz_3d ( r, xc, yc, zc, theta, phi, x, y, z )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_LOCAL2XYZ_3D converts local to XYZ coordinates on an implicit sphere in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !    The "local" spherical coordinates of a point are two angles, THETA and PHI.
! !    PHI measures the angle that the vector from the origin to the point
! !    makes with the positive Z axis.  THETA measures the angle that the
! !    projection of the vector onto the XY plane makes with the positive X axis.
! !
! !  Modified:
! !
! !    19 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
! !
! !    Input, real THETA, PHI, the local (THETA,PHI) spherical coordinates
! !    of a point on the sphere.  THETA and PHI are angles, measure in
! !    radians.  Usually, 0 <= THETA < 2 * PI, and 0 <= PHI <= PI.
! !
! !    Output, real X, Y, Z, the XYZ coordinates of the point.
! !
!   real phi
!   real r
!   real theta
!   real x
!   real xc
!   real y
!   real yc
!   real z
!   real zc
! !
!   x = xc + r * sin ( phi ) * cos ( theta )
!   y = yc + r * sin ( phi ) * sin ( theta )
!   z = zc + r * cos ( phi )
!
!   return
! end subroutine
! subroutine sphere_imp_near_point_3d ( r, xc, yc, zc, x, y, z, xn, yn, zn )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_NEAR_POINT_3D finds the nearest point on an implicit sphere to a point in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Method:
! !
! !    If the center of the sphere is (Xc,Yc,Zc), and the point is (X,Y,Z), then
! !    the desired point lies at a positive distance R along the vector (X-Xc, Y-Yc, Z-Zc)
! !    unless (X,Y,Z) = (Xc,Yc,Zc), in which case any point on the sphere is "nearest".
! !
! !  Modified:
! !
! !    14 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
! !
! !    Input, real X, Y, Z, the (X,Y,Z) coordinates of a point, whose
! !    nearest point on the sphere is desired.
! !
! !    Output, real XN, YN, ZN, the nearest point on the sphere.
! !
!   real norm
!   real r
!   real x
!   real xc
!   real xn
!   real y
!   real yc
!   real yn
!   real z
!   real zc
!   real zn
! !
! !  If (X,Y,Z) = (XC,YC,ZC), bail out now.
! !
!   norm = enorm0_3d ( x, y, z, xc, yc, zc )
!
!   if ( norm == 0.0E+00 ) then
!     xn = xc + r
!     yn = yc
!     zn = zc
!     return
!   end if
! !
! !  Compute the nearest point.
! !
!   xn = xc + r * ( x - xc ) / norm
!   yn = yc + r * ( y - yc ) / norm
!   zn = zc + r * ( z - zc ) / norm
!
!   return
! end subroutine
! subroutine sphere_imp_point_project_3d ( r, xc, yc, zc, x, y, z, xp, yp, zp )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_POINT_PROJECT_3D projects a point onto an implicit sphere, in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Modified:
! !
! !    24 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Input, real XC, YC, ZC, the coordinates of the center of the sphere.
! !
! !    Input, real X, Y, Z, the coordinates of a point.
! !
! !    Output, real X, Y, Z, the coordinates of the point as projected
! !    onto the sphere from the center.
! !
!   real norm
!   real r
!   real x
!   real xc
!   real xp
!   real y
!   real yc
!   real yp
!   real z
!   real zc
!   real zp
! !
!   if ( r == 0.0E+00 ) then
!
!     xp = xc
!     yp = yc
!     zp = zc
!
!   else if ( x == xc .and. y == yc .and. z == zc ) then
!
!     xp = xc
!     yp = yc
!     zp = zc + r
!
!   else
!
!     norm = enorm0_3d ( x, y, z, xc, yc, zc )
!
!     xp = xc + r * ( x - xc ) / norm
!     yp = yc + r * ( y - yc ) / norm
!     zp = zc + r * ( z - zc ) / norm
!
!   end if
!
!   return
! end subroutine
! subroutine sphere_imp_unit_area_nd ( n, area )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
! !
! !
! !  Discussion:
! !
! !    N   Area
! !
! !    2   2       * PI
! !    3   4       * PI
! !    4   2       * PI**2
! !    5   (8/3)   * PI**2
! !    6             PI**3
! !    7   (16/15) * PI**3
! !
! !  Modified:
! !
! !    26 October 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Output, real AREA, the area of the sphere.
! !
!   real area
!   integer i
!   integer m
!   integer n
! !
!   if ( mod ( n, 2 ) == 0 ) then
!     m = n / 2
!     area = 2.0E+00 * ( pi )**m
!     do i = 1, m-1
!       area = area / real ( i )
!     end do
!   else
!     m = ( n - 1 ) / 2
!     area = 2.0E+00**n * ( pi )**m
!     do i = m+1, 2*m
!       area = area / real ( i )
!     end do
!   end if
!
!   return
! end subroutine
! subroutine sphere_imp_unit_volume_nd ( n, volume )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
! !
! !
! !  Discussion:
! !
! !    N  Volume
! !
! !    2             PI
! !    3  (4/3)    * PI
! !    4  (1/2)    * PI**2
! !    5  (8/15)   * PI**2
! !    6  (1/6)    * PI**3
! !    7  (16/105) * PI**3
! !
! !  Modified:
! !
! !    26 October 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Output, real VOLUME, the volume of the sphere.
! !
!   integer i
!   integer m
!   integer n
!   real volume
! !
!   if ( mod ( n, 2 ) == 0 ) then
!     m = n / 2
!     volume = ( pi )**m
!     do i = 1, m
!       volume = volume / real ( i )
!     end do
!   else
!     m = ( n - 1 ) / 2
!     volume = ( pi )**m * 2.0E+00**n
!     do i = m+1, 2*m+1
!       volume = volume / real ( i )
!     end do
!   end if
!
!   return
! end subroutine
! subroutine sphere_imp_volume_3d ( r, volume )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_VOLUME_3D computes the volume of an implicit sphere in 3D.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Modified:
! !
! !    26 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Output, real VOLUME, the volume of the sphere.
! !
!   real r
!   real volume
! !
!   volume = ( 4.0E+00 / 3.0E+00 ) * pi * r**3
!
!   return
! end subroutine
! subroutine sphere_imp_volume_nd ( n, r, volume )
! !
! !*******************************************************************************
! !
! !! SPHERE_IMP_VOLUME_ND computes the volume of an implicit sphere in ND.
! !
! !
! !  Formula:
! !
! !    An implicit sphere in 3D satisfies the equation:
! !
! !      ( X - XC )**2 + ( Y - YC )**2 + ( Z - ZC )**2 = R**2
! !
! !  Discussion:
! !
! !    N  Volume
! !
! !    2             PI    * R**2
! !    3  (4/3)    * PI    * R**3
! !    4  (1/2)    * PI**2 * R**4
! !    5  (8/15)   * PI**2 * R**5
! !    6  (1/6)    * PI**3 * R**6
! !    7  (16/105) * PI**3 * R**7
! !
! !  Modified:
! !
! !    26 October 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the space.
! !
! !    Input, real R, the radius of the sphere.
! !
! !    Output, real VOLUME, the volume of the sphere.
! !
!   integer n
!   real r
!   real volume
! !
!   call sphere_imp_unit_volume_nd ( n, volume )
!
!   volume = volume * r**n
!
!   return
! end subroutine
! subroutine string_2d ( iorder, istrng, nstrng, nvec, x1vec, x2vec, y1vec, y2vec)
! !
! !*******************************************************************************
! !
! !! STRING_2D groups line segments into connected lines in 2D.
! !
! !
! !  Discussion:
! !
! !    The routine receives an unordered set of line segments, described by
! !    pairs of coordinates (X1,Y1) and (X2,Y2), and tries to group them
! !    into ordered lists that constitute connected jagged lines.
! !
! !    This routine will not match two endpoints unless they are exactly equal.
! !
! !  Modified:
! !
! !    28 February 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Output, integer IORDER(NVEC).
! !
! !    If two segments I and J have the same value of ISTRNG, then
! !    they belong to the same string.  Then IORDER(I) and IORDER(J)
! !    give the relative order of the two segments in the string.
! !    Thus if IORDER(I) = -3 and IORDER(J) = 2, then when the string
! !    is traversed, segment I is traversed first, then four other
! !    segments are traversed, and then segment J is traversed.
! !
! !    The segment with IORDER(I) = 0 is the initial segment from which
! !    the entire string was "grown".
! !
! !    Output, integer ISTRNG(NVEC), ISTRNG(I) is the number of the string
! !    to which line segment I belongs, between 1 and NVEC.
! !
! !    Output, integer NSTRNG, the number of strings created.
! !
! !    Input, integer NVEC, the number of line segments to be analyzed.
! !
! !    Input/output, real X1VEC(NVEC), X2VEC(NVEC), Y1VEC(NVEC), Y2VEC(NVEC).
! !
! !    On input, line segment I has endpoints ( X1VEC(I), Y1VEC(I) )
! !    and ( X2VEC(I), Y2VEC(I) ).
! !
! !    On output, the order of the components may have been
! !    switched.  That is, for some I, ( X1VEC(I), Y1VEC(I) ) and
! !    ( X2VEC(I), Y2VEC(I) ) may have been swapped.
! !
! !    More importantly, all the entries X1VEC(I), Y1VEC(I),
! !    Y1VEC(I) and Y2VEC(I) may have been swapped with another index J.
! !
! !    The resulting coordinates will have been sorted in order
! !    of the string to which they belong, and then by the order
! !    of their traversal within that string.
! !
!   integer nvec
! !
!   integer i
!   integer indx
!   integer iorder(nvec)
!   integer iseed
!   integer isgn
!   integer istrng(nvec)
!   integer itemp
!   integer j
!   integer jval
!   integer kval
!   integer match
!   integer nstrng
!   real temp
!   real x1val
!   real x1vec(nvec)
!   real x2val
!   real x2vec(nvec)
!   real y1val
!   real y1vec(nvec)
!   real y2val
!   real y2vec(nvec)
! !
! !  Mark ISTRNG so that each segment is alone.
! !
!   iorder(1:nvec) = 0
!   istrng(1:nvec) = nvec + i
! !
! !  Starting with the lowest numbered group of line segments,
! !  see if any higher numbered groups belong.
! !
!   iseed = 1
!   nstrng = 1
!   istrng(iseed) = nstrng
!
!   do
!
!     x1val = x1vec(iseed)
!     x2val = x2vec(iseed)
!     y1val = y1vec(iseed)
!     y2val = y2vec(iseed)
!
!     jval = iorder(iseed)
!     kval = iorder(iseed)
!
!     do
!
!       match = 0
!
!       do j = 1, nvec
!
!         if ( istrng(j) > nstrng ) then
!
!           if ( x1val == x1vec(j) .and. y1val == y1vec(j) ) then
!
!             jval = jval - 1
!             iorder(j) = jval
!             istrng(j) = nstrng
!             x1val = x2vec(j)
!             y1val = y2vec(j)
!             match = match + 1
!
!             call r_swap ( x1vec(j), x2vec(j) )
!             call r_swap ( y1vec(j), y2vec(j) )
!
!           else if ( x1val == x2vec(j) .and. y1val == y2vec(j) ) then
!
!             jval = jval - 1
!             iorder(j) = jval
!             istrng(j) = nstrng
!             x1val = x1vec(j)
!             y1val = y1vec(j)
!             match = match + 1
!
!           else if ( x2val == x1vec(j) .and. y2val == y1vec(j) ) then
!
!             kval = kval + 1
!             iorder(j) = kval
!             istrng(j) = nstrng
!             x2val = x2vec(j)
!             y2val = y2vec(j)
!             match = match + 1
!
!           else if ( x2val == x2vec(j) .and. y2val == y2vec(j) ) then
!
!             kval = kval + 1
!             iorder(j) = kval
!             istrng(j) = nstrng
!             x2val = x1vec(j)
!             y2val = y1vec(j)
!             match = match + 1
!
!             call r_swap ( x1vec(j), x2vec(j) )
!             call r_swap ( y1vec(j), y2vec(j) )
!
!           end if
!
!         end if
!
!       end do
! !
! !  If the string has closed on itself, then we don't want to
! !  look for any more matches for this string.
! !
!       if ( x1val == x2val .and. y1val == y2val ) then
!         exit
!       end if
! !
! !  If we made no matches this pass, we're done.
! !
!       if ( match <= 0 ) then
!         exit
!       end if
!
!     end do
! !
! !  This string is "exhausted".  Are there any line segments we
! !  haven't looked at yet?
! !
!     iseed = 0
!
!     do i = 1, nvec
!       if ( istrng(i) > nstrng ) then
!         iseed = i
!         nstrng = nstrng + 1
!         istrng(i) = nstrng
!         exit
!       end if
!     end do
!
!     if ( iseed == 0 ) then
!       exit
!     end if
!
!   end do
! !
! !  There are no more line segments to look at.  Renumber the
! !  isolated segments.
! !
! !  Question: Can this ever happen?
! !
!   do i = 1, nvec
!     if ( istrng(i) > nvec ) then
!       nstrng = nstrng + 1
!       istrng(i) = nstrng
!     end if
!   end do
! !
! !  Now sort the line segments by string and by order of traversal.
! !
!   i = 0
!   isgn = 0
!   j = 0
!
!   indx = 0
!
!   do
!
!     call sort_heap_external ( nvec, indx, i, j, isgn )
!
!     if ( indx > 0 ) then
!
!       call i_swap ( iorder(i), iorder(j) )
!       call i_swap ( istrng(i), istrng(j) )
!       call r_swap ( x1vec(i), x1vec(j) )
!       call r_swap ( y1vec(i), y1vec(j) )
!       call r_swap ( x2vec(i), x2vec(j) )
!       call r_swap ( y2vec(i), y2vec(j) )
!
!     else if ( indx < 0 ) then
!
!       if ( ( istrng(i) < istrng(j) ) .or. &
!            ( istrng(i) == istrng(j) .and. iorder(i) < iorder(j) ) ) then
!
!         isgn = - 1
!
!       else
!
!         isgn = + 1
!
!       end if
!
!     else if ( indx == 0 ) then
!
!       exit
!
!     end if
!
!   end do
!
!   return
! end subroutine
!
!
! function tand ( angle )
! !
! !*******************************************************************************
! !
! !! TAND returns the tangent of an angle given in degrees.
! !
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real ANGLE, the angle, in degrees.
! !
! !    Output, real TAND, the tangent of the angle.
! !
!   real angle
!   real angle_rad
!   real tand
! !
!   angle_rad = degrees_to_radians ( angle )
!   tand  = sin ( angle_rad ) / cos ( angle_rad )
!
!   return
! end function
!
!
! subroutine tetra_centroid_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, xc, yc, zc )
! !
! !*******************************************************************************
! !
! !! TETRA_CENTROID_3D computes the centroid of a tetrahedron in 3D.
! !
! !
! !  Modified:
! !
! !    05 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
! !    coordinates of the vertices.
! !
! !    Output, real XC, YC, ZC, the coordinates of the centroid.
! !
!   real x1
!   real x2
!   real x3
!   real x4
!   real xc
!   real y1
!   real y2
!   real y3
!   real y4
!   real yc
!   real z1
!   real z2
!   real z3
!   real z4
!   real zc
! !
!   xc = 0.25E+00 * ( x1 + x2 + x3 + x4 )
!   yc = 0.25E+00 * ( y1 + y2 + y3 + y4 )
!   zc = 0.25E+00 * ( z1 + z2 + z3 + z4 )
!
!   return
! end subroutine
!
! subroutine tetra_contains_point_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, x, y, z, inside )
! !
! !*******************************************************************************
! !
! !! TETRA_CONTAINS_POINT_3D finds if a point is inside a tetrahedron in 3D.
! !
! !
! !  Modified:
! !
! !    07 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, Z4, Y4, Z4, the
! !    coordinates of the tetrahedron vertices.
! !
! !    Input, real X, Y, Z, the point to be checked.
! !
! !    Output, logical INSIDE, is .TRUE. if (X,Y,Z) is inside
! !    the tetrahedron or on its boundary, and .FALSE. otherwise.
! !
!   integer, parameter :: N = 3
!   integer, parameter :: NRHS = 1
! !
!   real a(N,N+NRHS)
!   real c1
!   real c2
!   real c3
!   integer info
!   logical inside
!   real x
!   real x1
!   real x2
!   real x3
!   real x4
!   real y
!   real y1
!   real y2
!   real y3
!   real y4
!   real z
!   real z1
!   real z2
!   real z3
!   real z4
! !
! !  Set up the linear system
! !
! !    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
! !    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
! !    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
! !
! !  which is satisfied by the barycentric coordinates of (X,Y,Z).
! !
!   a(1,1) = x2 - x1
!   a(1,2) = x3 - x1
!   a(1,3) = x4 - x1
!   a(1,4) = x - x1
!
!   a(2,1) = y2 - y1
!   a(2,2) = y3 - y1
!   a(2,3) = y4 - y1
!   a(2,4) = y - y1
!
!   a(3,1) = z2 - z1
!   a(3,2) = z3 - z1
!   a(3,3) = z4 - z1
!   a(3,4) = z - z1
! !
! !  Solve the linear system.
! !
!   call rmat_solve ( a, N, NRHS, info )
!
!   if ( info /= 0 ) then
!     inside = .false.
!     return
!   end if
!
!   c1 = a(1,4)
!   c2 = a(2,4)
!   c3 = a(3,4)
! !
! !  If the point is in the tetrahedron, its barycentric coordinates
! !  must be nonnegative, and sum to no more than 1.
! !
!   if ( c1 < 0.0E+00 .or. c2 < 0.0E+00 .or. c3 < 0.0E+00 ) then
!     inside = .false.
!   else if ( c1 + c2 + c3 > 1.0E+00 ) then
!     inside = .false.
!   else
!     inside = .true.
!   end if
!
!   return
! end subroutine
! subroutine tetra_outsphere_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, r, xc, yc, zc )
! !
! !*******************************************************************************
! !
! !! TETRA_OUTSPHERE_3D computes the exscribed sphere of a tetrahedron in 3D.
! !
! !
! !  Discussion:
! !
! !    The exscribed sphere of a tetrahedron is the smallest sphere that
! !    can be drawn around the tetrahedron.  It passes through the four
! !    vertices.
! !
! !    Surprisingly, the diameter of the sphere can be found by solving
! !    a 3 by 3 linear system.  This is because the vectors P2 - P1,
! !    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
! !    right triangle with the diameter.  Hence, the dot product of
! !    P2 - P1 with the diameter is equal to the square of the length
! !    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
! !    the diameter vector originating at P1.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    05 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4,
! !    the coordinates of the vertices.
! !
! !    Output, real R, XC, YC, ZC, the coordinates of the center of the
! !    exscribed sphere, and its radius.  If the linear system is
! !    singular, then R = -1, XC = YC = ZC = 0.
! !
!   integer, parameter :: N = 3
!   integer, parameter :: NRHS = 1
! !
!   real a(N,N+NRHS)
!   integer info
!   real r
!   real x1
!   real x2
!   real x3
!   real x4
!   real xc
!   real y1
!   real y2
!   real y3
!   real y4
!   real yc
!   real z1
!   real z2
!   real z3
!   real z4
!   real zc
! !
! !  Set up the linear system.
! !
!   a(1,1) = x2 - x1
!   a(1,2) = y2 - y1
!   a(1,3) = z2 - z1
!   a(1,4) = enormsq0_3d ( x1, y1, z1, x2, y2, z2 )
!
!   a(2,1) = x3 - x1
!   a(2,2) = y3 - y1
!   a(2,3) = z3 - z1
!   a(2,4) = enormsq0_3d ( x1, y1, z1, x3, y3, z3 )
!
!   a(3,1) = x4 - x1
!   a(3,2) = y4 - y1
!   a(3,3) = z4 - z1
!   a(3,4) = enormsq0_3d ( x1, y1, z1, x4, y4, z4 )
! !
! !  Solve the linear system.
! !
!   call rmat_solve ( a, N, NRHS, info )
! !
! !  Compute R, X, Y.
! !
!   if ( info /= 0 ) then
!     r = -1.0E+00
!     xc = 0.0E+00
!     yc = 0.0E+00
!     zc = 0.0E+00
!   else
!     r = 0.5E+00 * sqrt ( a(1,N+1)**2 + a(2,N+1)**2 + a(3,N+1)**2 )
!     xc = x1 + 0.5E+00 * a(1,N+1)
!     yc = y1 + 0.5E+00 * a(2,N+1)
!     zc = z1 + 0.5E+00 * a(3,N+1)
!   end if
!
!   return
! end subroutine
! subroutine tetra_shape_3d ( max_num, max_order, point_num, face_num, &
!    face_order, point_coord, face_point )
! !
! !*******************************************************************************
! !
! !! TETRA_SHAPE_3D describes a tetrahedron in 3D.
! !
! !
! !  Discussion:
! !
! !    The vertices of the tetrahedron lie on the unit sphere.
! !    The last point is the north pole of the unit sphere.
! !
! !  Modified:
! !
! !    12 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer MAX_NUM, the maximum number of faces, and of
! !    points.  This quantity is used to dimension arrays.
! !    If POINT_NUM or FACE_NUM exceeds MAX_NUM, the arrays will
! !    not be set.
! !
! !    Input, integer MAX_ORDER, the maximum order of any face.
! !    This quantity is used to dimension arrays.  If FACE_ORDER
! !    exceeds MAX_ORDER, the arrays will not be set.
! !
! !    Output, integer POINT_NUM, the number of points in the shape.
! !
! !    Output, real POINT_COORD(3,MAX_NUM); POINT_COORD(*,J) contains
! !    the X, Y and Z coordinates of point J.
! !
! !    Output, integer FACE_NUM, the number of faces in the shape.
! !
! !    Output, integer FACE_ORDER, the number of vertices per face.
! !
! !    Output, integer FACE_POINT(MAX_ORDER,MAX_NUM); FACE_POINT(I,J)
! !    contains the index of the I-th point in the J-th face.  The
! !    points are listed in the counter-clockwise direction defined
! !    by the outward normal at the face.
! !
!   integer max_num
!   integer max_order
! !
!   integer face_num
!   integer face_order
!   integer face_point(max_order,max_num)
!   integer point_num
!   real point_coord(3,max_num)
! !
!   point_num = 4
!   face_num = 4
!   face_order = 3
! !
! !  Check.
! !
!   if ( point_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'TETRA_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of vertices exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_num > max_num ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'TETRA_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Number of faces exceeds MAX_NUM.'
!     return
!   end if
!
!   if ( face_order > max_order ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'TETRA_SHAPE_3D - Fatal error!'
!     write ( *, * ) '  Face order exceeds MAX_ORDER.'
!     return
!   end if
! !
! !  Set point coordinates.
! !
!   point_coord(1,1) =   0.0E+00
!   point_coord(2,1) =   2.0E+00 * sqrt ( 2.0E+00 ) / 3.0E+00
!   point_coord(3,1) = - 1.0E+00 / 3.0E+00
!
!   point_coord(1,2) = - sqrt ( 2.0E+00 / 3.0E+00 )
!   point_coord(2,2) = - sqrt ( 2.0E+00 ) / 3.0E+00
!   point_coord(3,2) = - 1.0E+00 / 3.0E+00
!
!   point_coord(1,3) =   sqrt ( 2.0E+00 / 3.0E+00 )
!   point_coord(2,3) = - sqrt ( 2.0E+00 ) / 3.0E+00
!   point_coord(3,3) = - 1.0E+00 / 3.0E+00
!
!   point_coord(1,4) = 0.0E+00
!   point_coord(2,4) = 0.0E+00
!   point_coord(3,4) = 1.0E+00
! !
! !  Set faces.
! !
!   face_point(1,1) = 1
!   face_point(2,1) = 3
!   face_point(3,1) = 2
!
!   face_point(1,2) = 1
!   face_point(2,2) = 2
!   face_point(3,2) = 4
!
!   face_point(1,3) = 1
!   face_point(2,3) = 4
!   face_point(3,3) = 3
!
!   face_point(1,4) = 2
!   face_point(2,4) = 3
!   face_point(3,4) = 4
!
!   return
! end subroutine
! subroutine tetra_volume_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x4, y4, z4, volume )
! !
! !*******************************************************************************
! !
! !! TETRA_VOLUME_3D computes the volume of a tetrahedron in 3D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, the
! !    coordinates of the corners of the tetrahedron.
! !
! !    Output, real VOLUME, the volume of the tetrahedron.
! !
!   real a(4,4)
!   real volume
!   real x1
!   real x2
!   real x3
!   real x4
!   real y1
!   real y2
!   real y3
!   real y4
!   real z1
!   real z2
!   real z3
!   real z4
! !
!   a(1:4,1) = (/ x1, x2, x3, x4 /)
!
!   a(1:4,2) = (/ y1, y2, y3, y4 /)
!
!   a(1:4,3) = (/ z1, z2, z3, z4 /)
!
!   a(1:4,4) = (/ 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00 /)
!
!   volume = abs ( rmat4_det ( a ) ) / 6.0E+00
!
!   return
! end subroutine
! subroutine tmat_init ( a )
! !
! !*******************************************************************************
! !
! !! TMAT_INIT initializes the geometric transformation matrix.
! !
! !
! !  Definition:
! !
! !    The geometric transformation matrix can be thought of as a 4 by 4
! !    matrix "A" having components:
! !
! !      r11 r12 r13 t1
! !      r21 r22 r23 t2
! !      r31 r32 r33 t3
! !        0   0   0  1
! !
! !    This matrix encodes the rotations, scalings and translations that
! !    are applied to graphical objects.
! !
! !    A point P = (x,y,z) is rewritten in "homogeneous coordinates" as
! !    PH = (x,y,z,1).  Then to apply the transformations encoded in A to
! !    the point P, we simply compute A * PH.
! !
! !    Individual transformations, such as a scaling, can be represented
! !    by simple versions of the transformation matrix.  If the matrix
! !    A represents the current set of transformations, and we wish to
! !    apply a new transformation B, then the original points are
! !    transformed twice:  B * ( A * PH ).  The new transformation B can
! !    be combined with the original one A, to give a single matrix C that
! !    encodes both transformations: C = B * A.
! !
! !  Modified:
! !
! !    19 October 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the geometric transformation matrix.
! !
!   real a(4,4)
!   integer i
!   integer j
! !
!   do i = 1, 4
!     do j = 1, 4
!       if ( i == j ) then
!         a(i,j) = 1.0E+00
!       else
!         a(i,j) = 0.0E+00
!       end if
!     end do
!   end do
!
!   return
! end subroutine
! subroutine tmat_mxm ( a, b, c )
! !
! !*******************************************************************************
! !
! !! TMAT_MXM multiplies two geometric transformation matrices.
! !
! !
! !  Note:
! !
! !    The product is accumulated in a temporary array, and then assigned
! !    to the result.  Therefore, it is legal for any two, or all three,
! !    of the arguments to share memory.
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the first geometric transformation matrix.
! !
! !    Input, real B(4,4), the second geometric transformation matrix.
! !
! !    Output, real C(4,4), the product A * B.
! !
!   real a(4,4)
!   real b(4,4)
!   real c(4,4)
! !
!   c(1:4,1:4) = matmul ( a(1:4,1:4), b(1:4,1:4) )
!
!   return
! end subroutine
! subroutine tmat_mxp ( a, x, y )
! !
! !*******************************************************************************
! !
! !! TMAT_MXP multiplies a geometric transformation matrix times a point.
! !
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the geometric transformation matrix.
! !
! !    Input, real X(3), the point to be multiplied.  The fourth component
! !    of X is implicitly assigned the value of 1.
! !
! !    Output, real Y(3), the result of A*X.  The product is accumulated in
! !    a temporary vector, and then assigned to the result.  Therefore, it
! !    is legal for X and Y to share memory.
! !
!   real a(4,4)
!   real x(3)
!   real y(3)
!   real z(4)
! !
!   z(1:4) = (/ x(1:3), 1.0E+00 /)
!
!   y(1:3) = matmul ( a(1:3,1:4), z(1:4) )
!
!   return
! end subroutine
! subroutine tmat_mxp2 ( a, x, y, n )
! !
! !*******************************************************************************
! !
! !! TMAT_MXP2 multiplies a geometric transformation matrix times N points.
! !
! !
! !  Modified:
! !
! !    20 October 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the geometric transformation matrix.
! !
! !    Input, real X(3,N), the points to be multiplied.
! !
! !    Output, real Y(3,N), the transformed points.  Each product is
! !    accumulated in a temporary vector, and then assigned to the
! !    result.  Therefore, it is legal for X and Y to share memory.
! !
!   integer n
! !
!   real a(4,4)
!   integer i
!   integer j
!   integer k
!   real x(3,n)
!   real y(3,n)
!   real z(3)
! !
!   do k = 1, n
!
!     do i = 1, 3
!       z(i) = a(i,4)
!       do j = 1, 3
!         z(i) = z(i) + a(i,j) * x(j,k)
!       end do
!     end do
!
!     y(1:3,k) = z(1:3)
!
!   end do
!
!   return
! end subroutine
! subroutine tmat_mxv ( a, x, y )
! !
! !*******************************************************************************
! !
! !! TMAT_MXV multiplies a geometric transformation matrix times a vector.
! !
! !
! !  Modified:
! !
! !    12 August 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the geometric transformation matrix.
! !
! !    Input, real X(3), the vector to be multiplied.  The fourth component
! !    of X is implicitly assigned the value of 1.
! !
! !    Output, real Y(3), the result of A*X.  The product is accumulated in
! !    a temporary vector, and then assigned to the result.  Therefore, it
! !    is legal for X and Y to share memory.
! !
!   real a(4,4)
!   integer i
!   integer j
!   real x(3)
!   real y(3)
!   real z(3)
! !
!   do i = 1, 3
!     z(i) = 0.0E+00
!     do j = 1, 3
!       z(i) = z(i) + a(i,j) * x(j)
!     end do
!     z(i) = z(i) + a(i,4)
!   end do
!
!   y(1:3) = z(1:3)
!
!   return
! end subroutine
!
!
! subroutine tmat_rot_axis ( a, b, angle, axis )
! !
! !*******************************************************************************
! !
! !! TMAT_ROT_AXIS applies a coordinate axis rotation to the geometric transformation matrix.
! !
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the current geometric transformation matrix.
! !
! !    Output, real B(4,4), the modified geometric transformation matrix.
! !    A and B may share the same memory.
! !
! !    Input, real ANGLE, the angle, in degrees, of the rotation.
! !
! !    Input, character AXIS, is 'X', 'Y' or 'Z', specifying the coordinate
! !    axis about which the rotation occurs.
! !
!   real a(4,4)
!   real angle
!   real angle_rad
!   character axis
!   real b(4,4)
!   real c(4,4)
!   integer i
!   integer j
! !
!   angle_rad = degrees_to_radians ( angle )
!
!   call tmat_init ( c )
!
!   if ( axis == 'X' .or. axis == 'x' ) then
!     c(2,2) =   cos ( angle_rad )
!     c(2,3) = - sin ( angle_rad )
!     c(3,2) =   sin ( angle_rad )
!     c(3,3) =   cos ( angle_rad )
!   else if ( axis == 'Y' .or. axis == 'y' ) then
!     c(1,1) =   cos ( angle_rad )
!     c(1,3) =   sin ( angle_rad )
!     c(3,1) = - sin ( angle_rad )
!     c(3,3) =   cos ( angle_rad )
!   else if ( axis == 'Z' .or. axis == 'z' ) then
!     c(1,1) =   cos ( angle_rad )
!     c(1,2) = - sin ( angle_rad )
!     c(2,1) =   sin ( angle_rad )
!     c(2,2) =   cos ( angle_rad )
!   else
!     write ( *, * ) ' '
!     write ( *, * ) 'TMAT_ROT_AXIS - Fatal error!'
!     write ( *, * ) '  Illegal rotation axis: ', axis
!     write ( *, * ) '  Legal choices are ''X'', ''Y'', or ''Z''.'
!     return
!   end if
!
!   b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )
!
!   return
! end subroutine
! subroutine tmat_rot_vector ( a, b, angle, axis )
! !
! !*******************************************************************************
! !
! !! TMAT_ROT_VECTOR applies an arbitrary axis rotation to the geometric transformation matrix.
! !
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the current geometric transformation matrix.
! !
! !    Output, real B(4,4), the modified geometric transformation matrix.
! !    A and B may share the same memory.
! !
! !    Input, real ANGLE, the angle, in degrees, of the rotation.
! !
! !    Input, real AXIS(3), the axis vector about which rotation occurs.
! !    AXIS may not be the zero vector.
! !
!   real a(4,4)
!   real angle
!   real angle_rad
!   real axis(3)
!   real b(4,4)
!   real c(4,4)
!   real ca
!   integer i
!   integer j
!   real norm
!   real sa
!   real v1
!   real v2
!   real v3
! !
!   v1 = axis(1)
!   v2 = axis(2)
!   v3 = axis(3)
!
!   norm = sqrt ( v1 * v1 + v2 * v2 + v3 * v3 )
!
!   if ( norm == 0.0E+00 ) then
!     return
!   end if
!
!   v1 = v1 / norm
!   v2 = v2 / norm
!   v3 = v3 / norm
!
!   angle_rad = degrees_to_radians ( angle )
!   ca = cos ( angle_rad )
!   sa = sin ( angle_rad )
!
!   call tmat_init ( c )
!
!   c(1,1) =                    v1 * v1 + ca * ( 1.0E+00 - v1 * v1 )
!   c(1,2) = ( 1.0E+00 - ca ) * v1 * v2 - sa * v3
!   c(1,3) = ( 1.0E+00 - ca ) * v1 * v3 + sa * v2
!
!   c(2,1) = ( 1.0E+00 - ca ) * v2 * v1 + sa * v3
!   c(2,2) =                    v2 * v2 + ca * ( 1.0E+00 - v2 * v2 )
!   c(2,3) = ( 1.0E+00 - ca ) * v2 * v3 - sa * v1
!
!   c(3,1) = ( 1.0E+00 - ca ) * v3 * v1 - sa * v2
!   c(3,2) = ( 1.0E+00 - ca ) * v3 * v2 + sa * v1
!   c(3,3) =                    v3 * v3 + ca * ( 1.0E+00 - v3 * v3 )
!
!   b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )
!
!   return
! end subroutine
! subroutine tmat_scale ( a, b, sx, sy, sz )
! !
! !*******************************************************************************
! !
! !! TMAT_SCALE applies a scaling to the geometric transformation matrix.
! !
! !
! !  Modified:
! !
! !    19 October 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the current geometric transformation matrix.
! !
! !    Output, real B(4,4), the modified geometric transformation matrix.
! !    A and B may share the same memory.
! !
! !    Input, real SX, SY, SZ, the scalings to be applied to the X, Y and
! !    Z coordinates.
! !
!   real a(4,4)
!   real b(4,4)
!   real c(4,4)
!   integer i
!   integer j
!   real sx
!   real sy
!   real sz
! !
!   call tmat_init ( c )
!
!   c(1,1) = sx
!   c(2,2) = sy
!   c(3,3) = sz
!
!   b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )
!
!   return
! end subroutine
! subroutine tmat_shear ( a, b, axis, s )
! !
! !*******************************************************************************
! !
! !! TMAT_SHEAR applies a shear to the geometric transformation matrix.
! !
! !
! !  Modified:
! !
! !    19 October 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the current geometric transformation matrix.
! !
! !    Output, real B(4,4), the modified geometric transformation matrix.
! !    A and B may share the same memory.
! !
! !    Input, character ( len = 2 ) AXIS, is 'XY', 'XZ', 'YX', 'YZ', 'ZX' or 'ZY',
! !    specifying the shear equation:
! !
! !      XY:  x' = x + s * y;
! !      XZ:  x' = x + s * z;
! !      YX:  y' = y + s * x;
! !      YZ:  y' = y + s * z;
! !      ZX:  z' = z + s * x;
! !      ZY:  z' = z + s * y.
! !
! !    Input, real S, the shear coefficient.
! !
!   real a(4,4)
!   character ( len = 2 ) axis
!   real b(4,4)
!   real c(4,4)
!   integer i
!   integer j
!   real s
! !
!   call tmat_init ( c )
!
!   if ( axis == 'XY' .or. axis == 'xy' ) then
!     c(1,2) = s
!   else if ( axis == 'XZ' .or. axis == 'xz' ) then
!     c(1,3) = s
!   else if ( axis == 'YX' .or. axis == 'yx' ) then
!     c(2,1) = s
!   else if ( axis == 'YZ' .or. axis == 'yz' ) then
!     c(2,3) = s
!   else if ( axis == 'ZX' .or. axis == 'zx' ) then
!     c(3,1) = s
!   else if ( axis == 'ZY' .or. axis == 'zy' ) then
!     c(3,2) = s
!   else
!     write ( *, * ) ' '
!     write ( *, * ) 'TMAT_SHEAR - Fatal error!'
!     write ( *, * ) '  Illegal shear axis: ', axis
!     write ( *, * ) '  Legal choices are XY, XZ, YX, YZ, ZX, or ZY.'
!     return
!   end if
!
!   b(1:4,1:4) = matmul ( c(1:4,1:4), a(1:4,1:4) )
!
!   return
! end subroutine
! subroutine tmat_trans ( a, b, x, y, z )
! !
! !*******************************************************************************
! !
! !! TMAT_TRANS applies a translation to the geometric transformation matrix.
! !
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Reference:
! !
! !    Foley, van Dam, Feiner, Hughes,
! !    Computer Graphics, Principles and Practice,
! !    Addison Wesley, Second Edition, 1990.
! !
! !  Parameters:
! !
! !    Input, real A(4,4), the current geometric transformation matrix.
! !
! !    Output, real B(4,4), the modified transformation matrix.
! !    A and B may share the same memory.
! !
! !    Input, real X, Y, Z, the translation.  This may be thought of as the
! !    point that the origin moves to under the translation.
! !
!   real a(4,4)
!   real b(4,4)
!   real x
!   real y
!   real z
! !
!   b(1:4,1:4) = a(1:4,1:4)
!
!   b(1:3,4) = b(1:3,4) + (/ x, y, z /)
!
!   return
! end subroutine
! function torus_area_3d ( r1, r2 )
! !
! !*******************************************************************************
! !
! !! TORUS_AREA_3D returns the area of a torus in 3D.
! !
! !
! !  Integration region:
! !
! !    Points (X,Y,Z) such that:
! !
! !      ( SQRT ( X**2 + Y**2 ) - R1 )**2 + Z**2 = R2**2.
! !
! !  Modified:
! !
! !    07 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R1, R2, the two radii that define the torus.
! !
! !    Output, real TORUS_AREA_3D, the area of the torus.
! !
!   real r1
!   real r2
!   real torus_area_3d
! !
!   torus_area_3d = 4.0E+00 * pi**2 * r1 * r2
!
!   return
! end function
! subroutine torus_volume_3d ( r1, r2, volume )
! !
! !*******************************************************************************
! !
! !! TORUS_VOLUME_3D computes the volume of a torus in 3D.
! !
! !
! !  Definition:
! !
! !    A torus with radii R1 and R2 is the set of points (X,Y,Z)
! !    satisfying:
! !
! !    ( sqrt ( X**2 + Y**2 ) - R1 )**2 + Z**2 <= R2**2
! !
! !  Modified:
! !
! !    11 December 1998
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real R1, R2, the "inner" and "outer" radii of the torus.
! !
! !    Output, real VOLUME, the volume of the torus.
! !
!   real r1
!   real r2
!   real volume
! !
!   volume = 2.0E+00 * ( pi )**2 * r1 * r2**2
!
!   return
! end subroutine
! subroutine triangle_angles_2d ( x1, y1, x2, y2, x3, y3, angle1, angle2, angle3 )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
! !
! !
! !  Discussion:
! !
! !    The law of cosines is used:
! !
! !      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
! !
! !    where GAMMA is the angle opposite side C.
! !
! !  Modified:
! !
! !    15 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates
! !    of the corners P1, P2 and P3 of the triangle.
! !
! !    Output, real ANGLE1, ANGLE2, ANGLE3, the angles opposite
! !    sides P1-P2, P2-P3 and P3-P1, in radians.
! !
!   real a
!   real angle1
!   real angle2
!   real angle3
!   real b
!   real c
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
! !
!   a = enorm0_2d ( x2, y2, x1, y1 )
!   b = enorm0_2d ( x3, y3, x2, y2 )
!   c = enorm0_2d ( x1, y1, x3, y3 )
! !
! !  Take care of a ridiculous special case.
! !
!   if ( a == 0.0E+00 .and. b == 0.0E+00 .and. c == 0.0E+00 ) then
!     angle1 = 2.0E+00 * pi / 3.0E+00
!     angle2 = 2.0E+00 * pi / 3.0E+00
!     angle3 = 2.0E+00 * pi / 3.0E+00
!     return
!   end if
!
!   if ( c == 0.0E+00 .or. a == 0.0E+00 ) then
!     angle1= pi
!   else
!     angle1 = arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0E+00 * c * a ) )
!   end if
!
!   if ( a == 0.0E+00 .or. b == 0.0E+00 ) then
!     angle2 = pi
!   else
!     angle2 = arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0E+00 * a * b ) )
!   end if
!
!   if ( b == 0.0E+00 .or. c == 0.0E+00 ) then
!     angle3 = pi
!   else
!     angle3 = arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0E+00 * b * c ) )
!   end if
!
!   return
! end subroutine
! subroutine triangle_angles_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, angle1, &
!   angle2, angle3 )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
! !
! !
! !  Discussion:
! !
! !    The law of cosines is used:
! !
! !      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
! !
! !    where GAMMA is the angle opposite side C.
! !
! !  Modified:
! !
! !    15 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of the corners P1, P2 and P3 of the triangle.
! !
! !    Output, real ANGLE1, ANGLE2, ANGLE3, the angles opposite
! !    sides P1-P2, P2-P3 and P3-P1, in radians.
! !
!   real a
!   real angle1
!   real angle2
!   real angle3
!   real b
!   real c
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
!   real z1
!   real z2
!   real z3
! !
!   a = enorm0_3d ( x2, y2, z2, x1, y1, z1 )
!   b = enorm0_3d ( x3, y3, z3, x2, y2, z2 )
!   c = enorm0_3d ( x1, y1, z1, x3, y3, z3 )
! !
! !  Take care of a ridiculous special case.
! !
!   if ( a == 0.0E+00 .and. b == 0.0E+00 .and. c == 0.0E+00 ) then
!     angle1 = 2.0E+00 * pi / 3.0E+00
!     angle2 = 2.0E+00 * pi / 3.0E+00
!     angle3 = 2.0E+00 * pi / 3.0E+00
!     return
!   end if
!
!   if ( c == 0.0E+00 .or. a == 0.0E+00 ) then
!     angle1= pi
!   else
!     angle1 = arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0E+00 * c * a ) )
!   end if
!
!   if ( a == 0.0E+00 .or. b == 0.0E+00 ) then
!     angle2 = pi
!   else
!     angle2 = arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0E+00 * a * b ) )
!   end if
!
!   if ( b == 0.0E+00 .or. c == 0.0E+00 ) then
!     angle3 = pi
!   else
!     angle3 = arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0E+00 * b * c ) )
!   end if
!
!   return
! end subroutine
! subroutine triangle_area_2d ( x1, y1, x2, y2, x3, y3, area )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates
! !    of the corners of the triangle.
! !
! !    Output, real AREA, the absolute area of the triangle.
! !
!   real area
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
! !
!   area = abs ( 0.5E+00 * ( x1 * ( y3 - y2 ) + x2 * ( y1 - y3 ) + &
!                        x3 * ( y2 - y1 ) ) )
!
!   return
! end subroutine
! subroutine triangle_area_signed_2d ( x1, y1, x2, y2, x3, y3, area_signed )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_AREA_SIGNED_2D computes the signed area of a triangle in 2D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates
! !    of the corners of the triangle.
! !
! !    Output, real AREA_SIGNED, the signed area of the triangle.  The size
! !    of AREA_SIGNED is the usual area, but the sign can be negative
! !    or positive, depending on the order in which the corners of the
! !    triangle were listed.
! !
!   real area_signed
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
! !
!   area_signed = 0.5E+00 * ( x1 * ( y3 - y2 ) + x2 * ( y1 - y3 ) + &
!                         x3 * ( y2 - y1 ) )
!
!   return
! end subroutine
! subroutine triangle_area_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, area )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_AREA_3D computes the area of a triangle in 3D.
! !
! !
! !  Discussion:
! !
! !    This routine uses the fact that the norm of the cross product vector
! !    is the area of the parallelogram they form.  The triangle they
! !    form has half that area.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    06 August 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the (X,Y,Z)
! !    coordinates of the corners of the triangle.
! !
! !    Output, real AREA, the area of the triangle.
! !
!   real area
!   real norm
!   real x1
!   real x2
!   real x3
!   real x4
!   real y1
!   real y2
!   real y3
!   real y4
!   real z1
!   real z2
!   real z3
!   real z4
! !
!   call cross_3d ( x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1, z3 - z1, &
!                   x4,      y4,      z4 )
!
!   norm = enorm_3d ( x4, y4, z4 )
!
!   area = 0.5E+00 * norm
!
!   return
! end subroutine
! subroutine triangle_area_2_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, area )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_AREA_2_3D computes the area of a triangle in 3D.
! !
! !
! !  Discussion:
! !
! !    This routine computes the area "the hard way".
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the (X,Y,Z)
! !    coordinates of the corners of the triangle.
! !
! !    Output, real AREA, the area of the triangle.
! !
!   real alpha
!   real area
!   real base
!   real dot
!   real height
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
!   real z1
!   real z2
!   real z3
! !
! !  Find the projection of (P3-P1) onto (P2-P1).
! !
!   dot = ( x2 - x1 ) * ( x3 - x1 ) + ( y2 - y1 ) * ( y3 - y1 ) + &
!         ( z2 - z1 ) * ( z3 - z1 )
!
!   base = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
! !
! !  The height of the triangle is the length of (P3-P1) after its
! !  projection onto (P2-P1) has been subtracted.
! !
!   if ( base == 0.0E+00 ) then
!
!     height = 0.0E+00
!
!   else
!
!     alpha = dot / ( base * base )
!
!     height = enorm0_3d ( x1 + alpha * ( x2 - x1 ), y1 + alpha * ( y2 - y1 ), &
!                          z1 + alpha * ( z2 - z1 ), x3, y3, z3 )
!
!   end if
!
!   area = 0.5E+00 * base * height
!
!   return
! end subroutine
! subroutine triangle_centroid_2d ( x1, y1, x2, y2, x3, y3, x, y )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_CENTROID_2D computes the centroid of a triangle in 2D.
! !
! !
! !  Discussion:
! !
! !    The centroid of a triangle can also be considered the center
! !    of gravity, assuming that the triangle is made of a thin uniform
! !    sheet of massy material.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    29 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates
! !    of the corners of the triangle.
! !
! !    Output, real X, Y, the coordinates of the centroid of the triangle.
! !
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!   x = ( x1 + x2 + x3 ) / 3.0E+00
!   y = ( y1 + y2 + y3 ) / 3.0E+00
!
!   return
! end subroutine
! subroutine triangle_centroid_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, x, y, z )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_CENTROID_3D computes the centroid of a triangle in 3D.
! !
! !
! !  Discussion:
! !
! !    The centroid of a triangle can also be considered the center
! !    of gravity, assuming that the triangle is made of a thin uniform
! !    sheet of massy material.
! !
! !  Modified:
! !
! !    04 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the coordinates
! !    of the corners of the triangle.
! !
! !    Output, real X, Y, Z, the coordinates of the centroid.
! !
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
!   real z
!   real z1
!   real z2
!   real z3
! !
!   x = ( x1 + x2 + x3 ) / 3.0E+00
!   y = ( y1 + y2 + y3 ) / 3.0E+00
!   z = ( z1 + z2 + z3 ) / 3.0E+00
!
!   return
! end subroutine
!
!
! subroutine triangle_contains_point_2d ( x1, y1, x2, y2, x3, y3, x, y, inside )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_CONTAINS_POINT_2D finds if a point is inside a triangle in 2D.
! !
! !
! !  Modified:
! !
! !    12 November 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of
! !    the corners of the triangle.
! !
! !    Input, real X, Y, the point to be checked.
! !
! !    Output, logical INSIDE, is .TRUE. if (X,Y) is inside
! !    the triangle or on its boundary, and .FALSE. otherwise.
! !
!   integer, parameter :: N = 2
!   integer, parameter :: NRHS = 1
! !
!   real a(N,N+NRHS)
!   real c1
!   real c2
!   integer info
!   logical inside
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!
! !  Set up the linear system
! !
! !    ( X2-X1  X3-X1 ) C1  = X-X1
! !    ( Y2-Y1  Y3-Y1 ) C2    Y-Y1
! !
! !  which is satisfied by the barycentric coordinates of (X,Y).
! !
!   a(1,1) = x2 - x1
!   a(1,2) = x3 - x1
!   a(1,3) = x - x1
!
!   a(2,1) = y2 - y1
!   a(2,2) = y3 - y1
!   a(2,3) = y - y1
! !
! !  Solve the linear system.
! !
!   call rmat_solve ( a, N, NRHS, info )
!
!   if ( info /= 0 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'TRIANGLE_CONTAINS_POINT - Fatal error!'
!     write ( *, * ) '  The linear system is singular.'
!     write ( *, * ) '  The input data does not form a proper triangle.'
!     stop
!   end if
!
!   c1 = a(1,3)
!   c2 = a(2,3)
! !
! !  If the point is in the triangle, its barycentric coordinates
! !  must both be nonnegative, and sum to no more than 1.
! !
!   if ( c1 < 0.0E+00 .or. c2 < 0.0E+00 ) then
!     inside = .false.
!   else if ( c1 + c2 > 1.0E+00 ) then
!     inside = .false.
!   else
!     inside = .true.
!   end if
!
!   return
! end subroutine
! subroutine triangle_diameter_2d ( x1, y1, x2, y2, x3, y3, diam )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_DIAMETER_2D computes the diameter of a triangle in 2D.
! !
! !
! !  Discussion:
! !
! !    The diameter of a triangle is the diameter of the smallest circle
! !    that can be drawn around the triangle.  At least two of the vertices
! !    of the triangle will intersect the circle, but not necessarily
! !    all three!
! !
! !  Modified:
! !
! !    13 February 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle corners.
! !
! !    Output, real DIAM, the diameter of the triangle.
! !
!   integer, parameter :: N = 2
!   integer, parameter :: NRHS = 1
! !
!   real a
!   real b
!   real c
!   real diam
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
! !
! !  Compute the (squares of) the lengths of the sides.
! !
!   a = ( x1 - x2 )**2 + ( y1 - y2 )**2
!   b = ( x2 - x3 )**2 + ( y2 - y3 )**2
!   c = ( x3 - x1 )**2 + ( y3 - y1 )**2
! !
! !  Take care of a zero side.
! !
!   if ( a == 0.0E+00 ) then
!     diam = sqrt ( b )
!     return
!   else if ( b == 0.0E+00 ) then
!     diam = sqrt ( c )
!     return
!   else if ( c == 0.0E+00 ) then
!     diam = sqrt ( a )
!     return
!   end if
! !
! !  Make A the largest.
! !
!   if ( a < b ) then
!     call r_swap ( a, b )
!   end if
!
!   if ( a < c ) then
!     call r_swap ( a, c )
!   end if
! !
! !  If A is very large...
! !
!   if ( a > b + c ) then
!     diam = sqrt ( a )
!   else
!     a = sqrt ( a )
!     b = sqrt ( b )
!     c = sqrt ( c )
!     diam = 2.0E+00 * a * b * c / sqrt ( ( a + b + c ) * ( - a + b + c ) &
!       * ( a - b + c ) * ( a + b - c ) )
!   end if
!
!   return
! end subroutine
! subroutine triangle_gridpoints_2d ( x1, y1, x2, y2, x3, y3, nsub, maxgrid, &
!   ngrid, x, y )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_GRIDPOINTS_2D computes gridpoints within a triangle in 2D.
! !
! !
! !  Discussion:
! !
! !    The gridpoints are computed by repeated halving of the triangle.
! !    The 0-th set of grid points is the corners themselves.
! !    The 1st set of grid points is the midpoints of the sides.
! !    These points can be used to draw 4 triangles that make up the original
! !    triangle.  The 2nd set of grid points is the side midpoints and centers
! !    of these four triangles.
! !
! !     NSUB                        NGRID
! !    -----                        -----
! !        0      1                  =  1  (centroid)
! !        1      1 + 2              =  3  (corners)
! !        2      1 + 2 + 3          =  6
! !        3      1 + 2 + 3 + 4      = 10
! !        4      1 + 2 + 3 + 4 + 5  = 15
! !
! !    NGRID is the sum of the integers from 1 to NSUB+1 or
! !
! !      NGRID = (NSUB+1) * (NSUB+2) / 2
! !
! !  Modified:
! !
! !    18 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates
! !    of the corners P1, P2 and P3 of the triangle.
! !
! !    Input, integer NSUB, the number of subdivisions.
! !
! !    Input, integer MAXGRID, the maximum number of grid points.
! !
! !    Output, integer NGRID, the number of grid points returned.
! !
! !    Output, real X(MAXGRID), Y(MAXGRID), coordinates of the grid points.
! !
!   integer maxgrid
! !
!   real a
!   real b
!   real c
!   integer i
!   integer j
!   integer k
!   integer ngrid
!   integer nsub
!   real x(maxgrid)
!   real x1
!   real x2
!   real x3
!   real y(maxgrid)
!   real y1
!   real y2
!   real y3
! !
!   ngrid = 0
! !
! !  Special case, NSUB = 0.
! !
!   if ( nsub == 0 ) then
!     if ( maxgrid >= 1 ) then
!       ngrid = 1
!       x(1) = ( x1 + x2 + x3 ) / 3.0E+00
!       y(1) = ( y1 + y2 + y3 ) / 3.0E+00
!     end if
!     return
!   end if
!
!   do i = 0, nsub
!
!     a = real ( i ) / real ( nsub )
!
!     do j = 0, nsub - i
!
!       b = real ( j ) / real ( nsub )
!       k = nsub - i - j
!       c = real ( k ) / real ( nsub )
!
!       if ( ngrid < maxgrid ) then
!
!         ngrid = ngrid + 1
!         x(ngrid) = ( a * x1 + b * x2 + c * x3 )
!         y(ngrid) = ( a * y1 + b * y2 + c * y3 )
!
!       end if
!
!     end do
!   end do
!
!   return
! end subroutine
! subroutine triangle_incircle_2d ( x1, y1, x2, y2, x3, y3, r, xc, yc )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
! !
! !
! !  Discussion:
! !
! !    The inscribed circle of a triangle is the largest circle that can
! !    be drawn inside the triangle.  It is tangent to all three sides,
! !    and the lines from its center to the vertices bisect the angles
! !    made by each vertex.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    29 January 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates
! !    of the corners of the triangle.
! !
! !    Output, real R, XC, YC, the radius and coordinates of the center of the
! !    inscribed circle.
! !
!   real perim
!   real r
!   real s12
!   real s23
!   real s31
!   real x1
!   real x2
!   real x3
!   real xc
!   real y1
!   real y2
!   real y3
!   real yc
! !
!   s12 = enorm0_2d ( x1, y1, x2, y2 )
!   s23 = enorm0_2d ( x2, y2, x3, y3 )
!   s31 = enorm0_2d ( x3, y3, x1, y1 )
!   perim = s12 + s23 + s31
!
!   xc = ( s23 * x1 + s31 * x2 + s12 * x3 ) / perim
!   yc = ( s23 * y1 + s31 * y2 + s12 * y3 ) / perim
!
!   r = 0.5E+00 * sqrt ( ( - s12 + s23 + s31 ) * ( + s12 - s23 + s31 ) &
!                  * ( + s12 + s23 - s31 ) / perim )
!
!   return
! end subroutine
! subroutine triangle_line_imp_int_2d ( a, b, c, trix, triy, nin, xint, yint )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_LINE_IMP_INT_2D finds where an implicit line intersects a triangle in 2D.
! !
! !
! !  Definition:
! !
! !    An implicit line is the set of points ( X, Y ) satisfying
! !
! !      A * X + B * Y + C = 0
! !
! !    where at least one of A and B is not zero.
! !
! !  Discussion:
! !
! !    If the line happens to be identical with one of the sides
! !    of the triangle, NIN will be returned as 3.
! !
! !    If the intersection point is one of the nodes of the triangle,
! !    it will be counted twice, which will result in a value of NIN
! !    that is 2 or 3, depending if the line also intersects
! !    the side opposite the node.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real A, B, C, determine the equation of the line:
! !    A*X + B*Y + C = 0.
! !
! !    Input, real TRIX(3), TRIY(3), the coordinates of the triangle corners.
! !
! !    Output, integer NIN, the number of points of intersection
! !    of the line with the triangle.  NIN may be 0, 1, 2 or 3.
! !
! !    Output, real XINT(3), YINT(3).  XINT(I), YINT(I) is the
! !    I-th intersection point, for I = 1 to NIN.
! !
!   real a
!   real a1
!   real b
!   real b1
!   real c
!   real c1
!   logical inside
!   integer ival
!   integer nin
!   real trix(3)
!   real triy(3)
!   real x
!   real xint(3)
!   real y
!   real yint(3)
! !
!   nin = 0
! !
! !  Get the implicit form of the line through vertices 1 and 2.
! !
!   call line_exp2imp_2d ( trix(1), triy(1), trix(2), triy(2), a1, b1, c1 )
! !
! !  Seek an intersection with the original line.
! !
!   call lines_imp_int_2d ( a, b, c, a1, b1, c1, ival, x, y )
! !
! !  Determine if the intersection is inside the triangle.
! !
!   if ( ival == 1 ) then
!
!     call triangle_contains_point_2d ( trix(1), triy(1), &
!       trix(2), triy(2), trix(3), triy(3), x, y, inside )
!
!     if ( inside ) then
!       nin = nin + 1
!       xint(nin) = x
!       yint(nin) = y
!     end if
!
!   end if
! !
! !  Get the implicit form of the line through vertices 2 and 3.
! !
!   call line_exp2imp_2d ( trix(2), triy(2), trix(3), triy(3), a1, b1, c1 )
! !
! !  Seek an intersection with the original line.
! !
!   call lines_imp_int_2d ( a, b, c, a1, b1, c1, ival, x, y )
! !
! !  Determine if the intersection is inside the triangle.
! !
!   if ( ival == 1 ) then
!
!     call triangle_contains_point_2d ( trix(1), triy(1), &
!       trix(2), triy(2), trix(3), triy(3), x, y, inside )
!
!     if ( inside ) then
!       nin = nin + 1
!       xint(nin) = x
!       yint(nin) = y
!     end if
!
!   end if
! !
! !  Get the implicit form of the line through vertices 3 and 1.
! !
!   call line_exp2imp_2d ( trix(3), triy(3), trix(1), triy(1), a1, b1, c1 )
! !
! !  Seek an intersection with the original line.
! !
!   call lines_imp_int_2d ( a, b, c, a1, b1, c1, ival, x, y )
! !
! !  Determine if the intersection is inside the triangle.
! !
!   if ( ival == 1 ) then
!
!     call triangle_contains_point_2d ( trix(1), triy(1), &
!       trix(2), triy(2), trix(3), triy(3), x, y, inside )
!
!     if ( inside ) then
!       nin = nin + 1
!       xint(nin) = x
!       yint(nin) = y
!     end if
!
!   end if
!
!   return
! end subroutine
!
! function triangle_orientation_2d ( x1, y1, x2, y2, x3, y3 )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
! !
! !
! !  Discussion:
! !
! !    Three distinct non-colinear points in the plane define a circle.
! !    If the points are visited in the order (x1,y1), (x2,y2), and then
! !    (x3,y3), this motion defines a clockwise or counterclockwise
! !    rotation along the circle.
! !
! !  Modified:
! !
! !    23 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates of
! !    three points, given in order.
! !
! !    Output, integer TRIANGLE_ORIENTATION_2D, reports if the three points lie
! !    clockwise on the circle that passes through them.  The possible
! !    return values are:
! !
! !    0, the points are distinct, noncolinear, and lie counterclockwise
! !    on their circle.
! !
! !    1, the points are distinct, noncolinear, and lie clockwise
! !    on their circle.
! !
! !    2, the points are distinct and colinear.
! !
! !    3, at least two of the points are identical.
! !
!   real det
!   integer triangle_orientation_2d
!   real x1
!   real x2
!   real x3
!   real y1
!   real y2
!   real y3
! !
!   if ( ( x1 == x2 .and. y1 == y2 ) .or. ( x2 == x3 .and. y2 == y3 ) .or. &
!        ( x3 == x1 .and. y3 == y1 ) ) then
!     triangle_orientation_2d = 3
!     return
!   end if
!
!   det = ( x1 - x3 ) * (y2 - y3 ) - ( x2 - x3 ) * ( y1 - y3 )
!
!   if ( det == 0.0E+00 ) then
!     triangle_orientation_2d = 2
!   else if ( det < 0.0E+00 ) then
!     triangle_orientation_2d = 1
!   else if ( det > 0.0E+00 ) then
!     triangle_orientation_2d = 0
!   end if
!
!   return
! end function
!
! subroutine triangle_outcircle_2d ( x1, y1, x2, y2, x3, y3, r, x, y )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_OUTCIRCLE_2D computes the exscribed circle of a triangle in 2D.
! !
! !
! !  Discussion:
! !
! !    The exscribed circle of a triangle is the circle that passes through
! !    the three vertices of the triangle.  The exscribed circle contains
! !    the triangle, but it is not necessarily the smallest triangle to do so.
! !
! !    Surprisingly, the diameter of the circle can be found by solving
! !    a 2 by 2 linear system.  This is because the vectors P2 - P1
! !    and P3 - P1 are secants of the circle, and each forms a right
! !    triangle with the diameter.  Hence, the dot product of
! !    P2 - P1 with the diameter is equal to the square of the length
! !    of P2 - P1, and similarly for P3 - P1.  This determines the
! !    diameter vector originating at P1.
! !
! !  Reference:
! !
! !    Adrian Bowyer and John Woodwark,
! !    A Programmer's Geometry,
! !    Butterworths, 1983.
! !
! !  Modified:
! !
! !    13 February 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle corners.
! !
! !    Output, real R, X, Y, the radius and coordinates of the center of the
! !    exscribed circle.  If the linear system is
! !    singular, then R = -1, X = Y = 0.
! !
!   integer, parameter :: N = 2
!   integer, parameter :: NRHS = 1
! !
!   real a(N,N+NRHS)
!   integer info
!   real r
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
! !  Set up the linear system.
! !
!   a(1,1) = x2 - x1
!   a(1,2) = y2 - y1
!   a(1,3) = enormsq0_2d ( x1, y1, x2, y2 )
!
!   a(2,1) = x3 - x1
!   a(2,2) = y3 - y1
!   a(2,3) = enormsq0_2d ( x1, y1, x3, y3 )
! !
! !  Solve the linear system.
! !
!   call rmat_solve ( a, N, NRHS, info )
! !
! !  Compute R, X, Y.
! !
!   if ( info /= 0 ) then
!     r = -1.0E+00
!     x = 0.0E+00
!     y = 0.0E+00
!   else
!     r = 0.5E+00 * sqrt ( a(1,N+1) * a(1,N+1) + a(2,N+1) * a(2,N+1) )
!     x = x1 + 0.5E+00 * a(1,N+1)
!     y = y1 + 0.5E+00 * a(2,N+1)
!   end if
!
!   return
! end subroutine
! subroutine triangle_point_dist_2d ( x1, y1, x2, y2, x3, y3, x, y, dist )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_POINT_DIST_2D: distance ( triangle, point ) in 2D.
! !
! !
! !  Modified:
! !
! !    15 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the corners of the triangle.
! !
! !    Input, real X, Y, the point which is to be checked.
! !
! !    Output, real DIST, the distance from the point to the triangle.
! !    DIST is zero if the point lies exactly on the triangle.
! !
!   real dist
!   real dist2
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
!   call line_seg_point_dist_2d ( x1, y1, x2, y2, x, y, dist2 )
!   dist = dist2
!
!   call line_seg_point_dist_2d ( x2, y2, x3, y3, x, y, dist2 )
!   dist = min ( dist, dist2 )
!
!   call line_seg_point_dist_2d ( x3, y3, x1, y1, x, y, dist2 )
!   dist = min ( dist, dist2 )
!
!   return
! end subroutine
! subroutine triangle_point_dist_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
!   x, y, z, dist )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_POINT_DIST_3D: distance ( triangle, point ) in 3D.
! !
! !
! !  Modified:
! !
! !    15 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, the triangle corners.
! !
! !    Input, real X, Y, Z, the point which is to be checked.
! !
! !    Output, real DIST, the distance from the point to the triangle.
! !    DIST is zero if the point lies exactly on the triangle.
! !
!   real dist
!   real dist2
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
!   real z
!   real z1
!   real z2
!   real z3
! !
! !  Compute the distances from the point to each of the sides.
! !
!   call line_seg_point_dist_3d ( x1, y1, z1, x2, y2, z2, x, y, z, dist2 )
!   dist = dist2
!
!   call line_seg_point_dist_3d ( x2, y2, z2, x3, y3, z3, x, y, z, dist2 )
!   dist = min ( dist, dist2 )
!
!   call line_seg_point_dist_3d ( x3, y3, z3, x1, y1, z1, x, y, z, dist2 )
!   dist = min ( dist, dist2 )
!
!   return
! end subroutine
! subroutine triangle_point_dist_signed_2d ( x1, y1, x2, y2, x3, y3, x, y, &
!   dist_signed )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_POINT_DIST_SIGNED_2D: signed distance ( triangle, point ) in 2D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle corners.
! !
! !    Input, real X, Y, the point which is to be checked.
! !
! !    Output, real DIST_SIGNED, the signed distance from the point to the
! !    triangle.  DIST_SIGNED is the maximum of the signeded distances from the
! !    point to each of the three lines that make up the triangle.
! !
! !    If DIST_SIGNED is:
! !    0, the point is on the boundary of the triangle;
! !    negative, the point is in the triangle;
! !    positive, the point is outside the triangle.
! !
!   real dis
!   real dis12
!   real dis23
!   real dis31
!   real dist_signed
!   real x
!   real x1
!   real x2
!   real x3
!   real y
!   real y1
!   real y2
!   real y3
! !
! !  Compare the signed distance from each line segment to the point,
! !  with the signed distance to the opposite vertex.
! !
! !  The signed distances should all be negative if the point is inside
! !  the triangle.
! !
!   call line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dis12 )
!   call line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x3, y3, dis )
!
!   if ( dis > 0.0E+00 ) then
!     dis = - dis
!     dis12 = - dis12
!   end if
!
!   call line_exp_point_dist_signed_2d ( x2, y2, x3, y3, x, y, dis23 )
!
!   call line_exp_point_dist_signed_2d ( x2, y2, x3, y3, x1, y1, dis )
!
!   if ( dis > 0.0E+00 ) then
!     dis = - dis
!     dis23 = - dis23
!   end if
!
!   call line_exp_point_dist_signed_2d ( x3, y3, x1, y1, x, y, dis31 )
!
!   call line_exp_point_dist_signed_2d ( x3, y3, x1, y1, x2, y2, dis )
!
!   if ( dis > 0.0E+00 ) then
!     dis = - dis
!     dis31 = - dis31
!   end if
!
!   dist_signed = max ( dis12, dis23, dis31 )
!
!   return
! end subroutine
! subroutine triangle_point_near_2d ( x1, y1, x2, y2, x3, y3, x, y, xn, yn, &
!   dist )
! !
! !*******************************************************************************
! !
! !! TRIANGLE_POINT_NEAR_2D computes the nearest triangle point to a point in 2D.
! !
! !
! !  Modified:
! !
! !    06 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, X2, Y2, X3, Y3, the triangle corners.
! !
! !    Input, real X, Y, the point whose nearest neighbor
! !    on the line is to be determined.
! !
! !    Output, real XN, YN, the nearest point to (X,Y).
! !
! !    Output, real DIST, the distance from the point to the triangle.
! !
!   real dist
!   real dist12
!   real dist23
!   real dist31
!   real t
!   real x
!   real x1
!   real x2
!   real x3
!   real xn
!   real xn12
!   real xn23
!   real xn31
!   real y
!   real y1
!   real y2
!   real y3
!   real yn
!   real yn12
!   real yn23
!   real yn31
! !
! !  Find the distance to each of the line segments that make up the edges
! !  of the triangle.
! !
!   call line_seg_point_near_2d ( x1, y1, x2, y2, x, y, xn12, yn12, dist12, t )
!
!   call line_seg_point_near_2d ( x2, y2, x3, y3, x, y, xn23, yn23, dist23, t )
!
!   call line_seg_point_near_2d ( x3, y3, x1, y1, x, y, xn31, yn31, dist31, t )
!
!   if ( dist12 <= dist23 .and. dist12 <= dist31 ) then
!     dist = dist12
!     xn = xn12
!     yn = yn12
!   else if ( dist23 <= dist12 .and. dist23 <= dist31 ) then
!     dist = dist23
!     xn = xn23
!     yn = yn23
!   else
!     dist = dist31
!     xn = xn31
!     yn = yn31
!   end if
!
!   return
! end subroutine
! subroutine tube_2d ( dist, n, x, y, x1, y1, x2, y2 )
! !
! !*******************************************************************************
! !
! !! TUBE_2D constructs a "tube" of given width around a path in 2D.
! !
! !
! !  Discussion:
! !
! !    The routine is given a sequence of N points, and a distance DIST.
! !
! !    It returns the coordinates of the corners of the top and bottom
! !    of a tube of width 2*DIST, which envelopes the line connecting
! !    the points.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real DIST, the radius of the tube.
! !
! !    Input, integer N, the number of points defining the line.
! !    N must be at least 2.
! !
! !    Input, real X(N), Y(N), the points which comprise the broken
! !    line which is to be surrounded by the tube.  Points should
! !    not be immediately repeated, that is, it should never be
! !    the case that
! !      X(I) = X(I+1) and Y(I) = Y(I+1).
! !
! !    Output, real X1(N), Y1(N), X2(N), Y2(N), the points ( X1(I), Y1(N) )
! !    form one side of the tube, and ( X2(I), Y2(I) ) the other.
! !
!   integer n
! !
!   real a
!   real b
!   real c
!   real dis1
!   real dis2
!   real dist
!   integer i
!   real temp
!   real x(n)
!   real x1(n)
!   real x2(n)
!   real xi
!   real xim1
!   real xip1
!   real y(n)
!   real y1(n)
!   real y2(n)
!   real yi
!   real yim1
!   real yip1
! !
! !  Check that N is at least 3.
! !
!   if ( n < 3 ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'TUBE_2D - Fatal error!'
!     write ( *, * ) '  N must be at least 3'
!     write ( *, * ) '  but your input value was N = ', n
!     stop
!   end if
! !
! !  Check that consecutive points are distinct.
! !
!   do i = 1, n-1
!     if ( x(i) == x(i+1) .and. y(i) == y(i+1) ) then
!       write ( *, * ) ' '
!       write ( *, * ) 'TUBE_2D - Fatal error!'
!       write ( *, * ) '  X(I) = X(I+1) and Y(I) = Y(I+1) for I = ', i
!       write ( *, * ) '  X(I), Y(I) = ', x(i), y(i)
!       stop
!     end if
!   end do
!
!   do i = 1, n
!
!     if ( i == 1 ) then
!       xim1 = x(i)
!       yim1 = y(i)
!     else
!       xim1 = x(i-1)
!       yim1 = y(i-1)
!     end if
!
!     xi = x(i)
!     yi = y(i)
!
!     if ( i < n ) then
!       xip1 = x(i+1)
!       yip1 = y(i+1)
!     else
!       xip1 = x(i)
!       yip1 = y(i)
!     end if
!
!     call corpl_2d ( dist, xim1, yim1, xi, yi, xip1, yip1, x1(i), &
!       y1(i), x2(i), y2(i) )
! !
! !  On the first and last steps, translate the corner points DIST units
! !  along the line, to make an extra buffer.
! !
!     if ( i == 1 ) then
!
!       temp = enorm0_2d ( x(1), y(1), x(2), y(2) )
!       x1(1) = x1(1) - dist * ( x(2) - x(1) ) / temp
!       y1(1) = y1(1) - dist * ( y(2) - y(1) ) / temp
!       x2(1) = x2(1) - dist * ( x(2) - x(1) ) / temp
!       y2(1) = y2(1) - dist * ( y(2) - y(1) ) / temp
!
!     else if ( i == n ) then
!
!       temp = enorm0_2d ( x(n-1), y(n-1), x(n), y(n) )
!       x1(n) = x1(n) + dist * ( x(n) - x(n-1) ) / temp
!       y1(n) = y1(n) + dist * ( y(n) - y(n-1) ) / temp
!       x2(n) = x2(n) + dist * ( x(n) - x(n-1) ) / temp
!       y2(n) = y2(n) + dist * ( y(n) - y(n-1) ) / temp
!
!     end if
! !
! !  The new points ( X1(I), Y1(I) ) and ( X2(I), Y2(I) ) may need to be
! !  swapped.
! !
! !  Compute the signed distance from the points to the line.
! !
!     if ( i > 1 ) then
!
!       a = y(i-1) - y(i)
!       b = x(i) - x(i-1)
!       c = x(i-1) * y(i) - x(i) * y(i-1)
!
!       dis1 = ( a * x1(i-1) + b * y1(i-1) + c ) / sqrt ( a * a + b * b )
!
!       dis2 = ( a * x1(i) + b * y1(i) + c ) / sqrt ( a * a + b * b )
!
!       if ( sign ( 1.0E+00, dis1 ) /= sign ( 1.0E+00, dis2 ) ) then
!
!         call r_swap ( x1(i), x2(i) )
!         call r_swap ( y1(i), y2(i) )
!
!       end if
!
!     end if
!
!   end do
!
!   return
! end subroutine
! subroutine vector_directions_2d ( x1, y1, ax, ay )
! !
! !*******************************************************************************
! !
! !! VECTOR_DIRECTIONS_2D returns the direction angles of a vector in 2D.
! !
! !
! !  Discussion:
! !
! !    Let V be the vector, and let E(I) be the I-th unit coordinate axis vector.
! !    The I-th direction angle is the angle between V and E(I), which is
! !    the angle whose cosine is equal to the direction cosine:
! !
! !      Direction_Cosine(I) = V dot E(I) / |V|.
! !
! !    If V is the null or zero vector, then the direction cosines and
! !    direction angles are undefined, and this routine simply returns
! !    zeroes.
! !
! !  Modified:
! !
! !    16 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, the coordinates of the endpoint of the vector V
! !    (which is presumed to begin at the origin).
! !
! !    Output, real AX, AY, the direction angles, in radians,
! !    that the vector V makes with the X and Y axes.
! !
!   real ax
!   real ay
!   real x1
!   real y1
! !
!   ax = atan2 ( y1, x1 )
!   ay = atan2 ( x1, y1 )
!
!   return
! end subroutine
!
!
! subroutine vector_directions_3d ( x1, y1, z1, ax, ay, az )
! !
! !*******************************************************************************
! !
! !! VECTOR_DIRECTIONS_3D returns the direction angles of a vector in 3D.
! !
! !
! !  Discussion:
! !
! !    Let V be the vector, and let E(I) be the I-th unit coordinate axis vector.
! !    The I-th direction angle is the angle between V and E(I), which is
! !    the angle whose cosine is equal to the direction cosine:
! !
! !      Direction_Cosine(I) = V dot E(I) / |V|.
! !
! !    If V is the null or zero vector, then the direction cosines and
! !    direction angles are undefined, and this routine simply returns
! !    zeroes.
! !
! !  Modified:
! !
! !    17 July 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, Z1, the coordinates of the endpoint of the vector V
! !    (which is presumed to begin at the origin).
! !
! !    Output, real AX, AY, AZ, the direction angles, in radians, that the
! !    vector V makes with the X, Y and Z axes.
! !
!   real ax
!   real ay
!   real az
!   real cos_x
!   real cos_y
!   real cos_z
!   real v1norm
!   real v2norm
!   real x1
!   real y1
!   real z1
! !
! !  Get the norm of the vector.
! !
!   v1norm = sqrt ( x1 * x1 + y1 * y1 + z1 * z1 )
!
!   if ( v1norm == 0.0E+00 ) then
!     ax = 0.0E+00
!     ay = 0.0E+00
!     az = 0.0E+00
!     return
!   end if
!
!   v2norm = 1.0E+00
! !
! !  Get the direction cosines.
! !
!   cos_x = x1 / ( v1norm * v2norm )
!   cos_y = y1 / ( v1norm * v2norm )
!   cos_z = z1 / ( v1norm * v2norm )
! !
! !  Retrieve the direction angles.
! !
!   ax = arc_cosine ( cos_x )
!   ay = arc_cosine ( cos_y )
!   az = arc_cosine ( cos_z )
!
!   return
! end subroutine
! subroutine vector_rotate_2d ( x1, y1, angle, x2, y2 )
! !
! !*******************************************************************************
! !
! !! VECTOR_ROTATE_2D rotates a vector around the origin in 2D.
! !
! !
! !  Discussion:
! !
! !    To see why this formula is so, consider that the original point
! !    has the form ( R cos Theta, R sin Theta ), and the rotated point
! !    has the form ( R cos ( Theta + Angle ), R sin ( Theta + Angle ) ).
! !    Now use the addition formulas for cosine and sine to relate
! !    the new point to the old one:
! !
! !      ( X2 ) = ( cos Angle  - sin Angle ) * ( X1 )
! !      ( Y2 )   ( sin Angle    cos Angle )   ( Y1 )
! !
! !  Modified:
! !
! !    19 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, the components of the vector to be rotated.
! !
! !    Input, real ANGLE, the angle, in radians, of the rotation to be
! !    carried out.  A positive angle rotates the vector in the
! !    counterclockwise direction.
! !
! !    Output, real X2, Y2, the rotated vector.
! !
!   real angle
!   real x1
!   real x2
!   real y1
!   real y2
! !
!   x2 = cos ( angle ) * x1 - sin ( angle ) * y1
!   y2 = sin ( angle ) * x1 + cos ( angle ) * y1
!
!   return
! end subroutine
! subroutine vector_rotate_base_2d ( x1, y1, xb, yb, angle, x2, y2 )
! !
! !*******************************************************************************
! !
! !! VECTOR_ROTATE_BASE_2D rotates a vector around a base point in 2D.
! !
! !
! !  Discussion:
! !
! !    The original vector is assumed to be ( X1-XB, Y1-YB ), and the
! !    rotated vector is ( X2-XB, Y2-YB ).
! !
! !  Modified:
! !
! !    13 March 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X1, Y1, the endpoint of the original vector.
! !
! !    Input, real XB, YB, the location of the base point.
! !
! !    Input, real ANGLE, the angle, in radians, of the rotation to be
! !    carried out.  A positive angle rotates the vector in the
! !    counterclockwise direction.
! !
! !    Output, real X2, Y2, the endpoint of the rotated vector.
! !
!   real angle
!   real x1
!   real x2
!   real xb
!   real y1
!   real y2
!   real yb
! !
!   x2 = xb + cos ( angle ) * ( x1 - xb ) - sin ( angle ) * ( y1 - yb )
!   y2 = yb + sin ( angle ) * ( x1 - xb ) + cos ( angle ) * ( y1 - yb )
!
!   return
! end subroutine
! subroutine vector_separation_3d ( v1, v2, theta )
! !
! !*******************************************************************************
! !
! !! VECTOR_SEPARATION_3D finds the angular separation between vectors in 3D.
! !
! !
! !  Discussion:
! !
! !    Any two vectors lie in a plane, and are separated by a plane angle.
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real V1(3), V2(3), the two vectors.
! !
! !    Output, real THETA, the angle between the two vectors.
! !
!   integer, parameter :: n = 3
! !
!   real cos_theta
!   real theta
!   real v1(n)
!   real v1_norm
!   real v2(n)
!   real v2_norm
! !
!   v1_norm = sqrt ( sum ( v1(1:n)**2 ) )
!
!   v2_norm = sqrt ( sum ( v2(1:n)**2 ) )
!
!   cos_theta = dot_product ( v1(1:n), v2(1:n) ) / ( v1_norm * v2_norm )
!
!   theta = arc_cosine ( cos_theta )
!
!   return
! end subroutine
! subroutine vector_separation_nd ( n, v1, v2, theta )
! !
! !*******************************************************************************
! !
! !! VECTOR_SEPARATION_ND finds the angular separation between vectors in ND.
! !
! !
! !  Discussion:
! !
! !    Any two vectors lie in a plane, and are separated by a plane angle.
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the dimension of the vectors.
! !
! !    Input, real V1(N), V2(N), the two vectors.
! !
! !    Output, real THETA, the angle between the two vectors.
! !
!   integer n
! !
!   real cos_theta
!   real theta
!   real v1(n)
!   real v1_norm
!   real v2(n)
!   real v2_norm
! !
!   v1_norm = sqrt ( sum ( v1(1:n)**2 ) )
!
!   v2_norm = sqrt ( sum ( v2(1:n)**2 ) )
!
!   cos_theta = dot_product ( v1(1:n), v2(1:n) ) / ( v1_norm * v2_norm )
!
!   theta = arc_cosine ( cos_theta )
!
!   return
! end subroutine
! subroutine vector_unit_2d ( x1, y1 )
! !
! !*******************************************************************************
! !
! !! VECTOR_UNIT_2D normalizes a vector in 2D.
! !
! !
! !  Modified:
! !
! !    07 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, real X1, Y1, the components of a 2D vector to be
! !    normalized.  On output, V should have unit Euclidean norm.
! !    However, if the input vector has zero Euclidean norm, it is
! !    not altered.
! !
!   real temp
!   real x1
!   real y1
! !
!   temp = sqrt ( x1 * x1 + y1 * y1 )
!
!   if ( temp /= 0.0E+00 ) then
!     x1 = x1 / temp
!     y1 = y1 / temp
!   end if
!
!   return
! end subroutine
! subroutine vector_unit_3d ( x1, y1, z1 )
! !
! !*******************************************************************************
! !
! !! VECTOR_UNIT_3D normalizes a vector in 3D.
! !
! !
! !  Modified:
! !
! !    07 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, real X1, Y1, Z1, the components of a 3D vector to be
! !    normalized.  On output, V should have unit Euclidean norm.
! !    However, if the input vector has zero Euclidean norm, it is
! !    not altered.
! !
!   real temp
!   real x1
!   real y1
!   real z1
! !
!   temp = sqrt ( x1 * x1 + y1 * y1 + z1 * z1 )
!
!   if ( temp /= 0.0E+00 ) then
!     x1 = x1 / temp
!     y1 = y1 / temp
!     z1 = z1 / temp
!   end if
!
!   return
! end subroutine
! subroutine vector_unit_nd ( n, v )
! !
! !*******************************************************************************
! !
! !! VECTOR_UNIT_ND normalizes a vector in ND.
! !
! !
! !  Modified:
! !
! !    07 February 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer N, the number of entries in the vector.
! !
! !    Input/output, real V(N), the vector to be normalized.  On output,
! !    V should have unit Euclidean norm.  However, if the input vector
! !    has zero Euclidean norm, it is not altered.
! !
!   integer n
! !
!   integer i
!   real temp
!   real v(n)
! !
!   temp = enorm_nd ( n, v )
!
!   if ( temp /= 0.0E+00 ) then
!     v(1:n) = v(1:n) / temp
!   end if
!
!   return
! end subroutine
! subroutine voxel_line_3d ( x1, y1, z1, x2, y2, z2 )
! !
! !*******************************************************************************
! !
! !! VOXEL_LINE_3D computes the voxels that line on a line in 3D.
! !
! !
! !  Discussion:
! !
! !    This preliminary version of the routine simply prints out the
! !    intermediate voxels, rather than returning them one at a time,
! !    or all together.
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Reference:
! !
! !    Daniel Cohen,
! !    Voxel Traversal along a 3D Line,
! !    Graphics Gems, 1994.
! !
! !  Parameters:
! !
! !    Input, integer X1, Y1, Z1, X2, Y2, Z2, the coordinates of the
! !    voxels that begin and end the line.
! !
!   integer ax
!   integer ay
!   integer az
!   integer exy
!   integer exz
!   integer ezy
!   integer i
!   integer n
!   integer sx
!   integer sy
!   integer sz
!   integer x
!   integer x1
!   integer x2
!   integer y
!   integer y1
!   integer y2
!   integer z
!   integer z1
!   integer z2
! !
!   sx = sign ( 1, x2 - x1 )
!   sy = sign ( 1, y2 - y1 )
!   sz = sign ( 1, z2 - z1 )
!
!   ax = abs ( x2 - x1 )
!   ay = abs ( y2 - y1 )
!   az = abs ( z2 - z1 )
!
!   exy = ay - ax
!   exz = az - ax
!   ezy = ay - az
!
!   n = ax + ay + az
!
!   x = x1
!   y = y1
!   z = z1
!
!   do i = 0, n
!
!     write ( *, '(4i6)' ) i, x, y, z
!
!     if ( exy < 0 ) then
!
!       if ( exz < 0 ) then
!         x = x + sx
!         exy = exy + 2 * ay
!         exz = exz + 2 * az
!       else
!         z = z + sz
!         exz = exz - 2 * ax
!         exy = exy + 2 * ay
!       end if
!
!     else if ( ezy < 0 ) then
!        z = z + sz
!        exz = exz - 2 * ax
!        ezy = ezy + 2 * ay
!     else
!       y = y + sy
!       exy = exy - 2 * ax
!       ezy = ezy - 2 * az
!     end if
!
!   end do
!
!   return
! end subroutine
! subroutine voxel_region_3d ( ishow, list, maxlist, nlist, nregion, nx, ny, nz )
! !
! !*******************************************************************************
! !
! !! VOXEL_REGION_3D arranges a set of voxels into contiguous regions in 3D.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input/output, integer ISHOW(NX,NY,NZ).
! !
! !    On input, ISHOW(I,J,K) has the values:
! !      0, if the voxel is OFF;
! !      anything else, if the voxel is ON.
! !
! !    On output, ISHOW(I,J,K) has the values:
! !      0, if the voxel is off,
! !      N, if the voxel is ON, and is part of region N.
! !
! !    Output, integer LIST(MAXLIST), contains, in stack form, a list
! !    of the indices of the elements in each region.
! !
! !    The number of elements in NREGION is NELEM = LIST(NLIST).  The
! !    (I,J,K) indices of the last element in this region are in
! !    LIST(NLIST-3) through LIST(NLIST-1), and the first element is
! !    listed in LIST(NLIST-3*NELEM), LIST(NLIST-3*NELEM+1),
! !    LIST(NLIST-3*NELEM+2).
! !
! !    The number of elements in NREGION-1 is listed in
! !    LIST(NLIST-3*NELEM-1), and so on.
! !
! !    Input, integer MAXLIST, the maximum length of the array used to
! !    list the elements of the regions.
! !
! !    Output, integer NLIST, the number of entries of LIST that were used.
! !    However, if NLIST > MAXLIST, then there was not enough space in
! !    LIST to store the data properly, and LIST should not be used,
! !    although the data in ISHOW should be correct.
! !
! !    Output, integer NREGION, the number of regions discovered.
! !
! !    Input, integer NX, NY, NZ, the number of voxels in the X, Y and
! !    Z directions.
! !
!   integer, parameter :: maxstack = 100
! !
!   integer maxlist
!   integer nx
!   integer ny
!   integer nz
! !
!   integer i
!   integer i2
!   integer ibase
!   integer ihi
!   integer ilo
!   integer ishow(nx,ny,nz)
!   integer j
!   integer j2
!   integer jbase
!   integer jhi
!   integer jlo
!   integer k
!   integer k2
!   integer kbase
!   integer khi
!   integer klo
!   integer list(maxlist)
!   integer nabes
!   integer ncan
!   integer nelements
!   integer nlist
!   integer nregion
!   integer nstack
!   integer stack(maxstack)
! !
! !  Reset all nonzero entries of ISHOW to -1.
! !
!   do i = 1, nx
!     do j = 1, ny
!       do k = 1, nz
!
!         if ( ishow(i,j,k) /= 0 ) then
!           ishow(i,j,k) = - 1
!         end if
!
!       end do
!     end do
!   end do
! !
! !  Start the number of items in the region list at 0.
! !
!   nlist = 0
! !
! !  Start the number of regions at 0.
! !
!   nregion = 0
! !
! !  The stack begins empty.
! !
!   nstack = 0
! !
! !  Search for an unused "ON" voxel from which we can "grow" a new region.
! !
!   do i = 1, nx
!     do j = 1, ny
!       do k = 1, nz
! !
! !  We found a voxel that is "ON", and does not belong to any region.
! !
!         if ( ishow(i,j,k) == - 1 ) then
! !
! !  Increase the number of regions.
! !
!           nregion = nregion + 1
! !
! !  Add this voxel to the region.
! !
!           ishow(i,j,k) = nregion
! !
! !  Add this voxel to the stack.
! !
!           if ( nstack + 4 > maxstack ) then
!             go to 20
!           end if
!
!           stack(nstack+1) = i
!           stack(nstack+2) = j
!           stack(nstack+3) = k
!
!           stack(nstack+4) = 1
!
!           nstack = nstack + 4
! !
! !  Add this voxel to the description of the region.
! !
!           nelements = 1
!
!           if ( nlist + 3 <= maxlist ) then
!             list(nlist+1) = i
!             list(nlist+2) = j
!             list(nlist+3) = k
!           end if
!
!           nlist = nlist + 3
!
! 10            continue
! !
! !  Find all neighbors of BASE that are "ON" but unused.
! !  Mark them as belonging to this region, and stack their indices.
! !
!           ibase = stack(nstack-3)
!           jbase = stack(nstack-2)
!           kbase = stack(nstack-1)
!
!           ilo = max ( ibase-1, 1 )
!           ihi = min ( ibase+1, nx )
!           jlo = max ( jbase-1, 1 )
!           jhi = min ( jbase+1, ny )
!           klo = max ( kbase-1, 1 )
!           khi = min ( kbase+1, nz )
!
!           nabes = 0
!
!           do i2 = ilo, ihi
!             do j2 = jlo, jhi
!               do k2 = klo, khi
! !
! !  We found a neighbor to our current search point, which is "ON" and unused.
! !
!                 if ( ishow(i2,j2,k2) == - 1 ) then
! !
! !  Increase the number of neighbors.
! !
!                   nabes = nabes + 1
! !
! !  Mark the neighbor as belonging to the region.
! !
!                   ishow(i2,j2,k2) = nregion
! !
! !  Add the neighbor to the stack.
! !
!                   if ( nstack+3 > maxstack ) then
!                     go to 20
!                   end if
!
!                   stack(nstack+1) = i2
!                   stack(nstack+2) = j2
!                   stack(nstack+3) = k2
!
!                   nstack = nstack+3
! !
! !  Add the neighbor to the description of the region.
! !
!                   nelements = nelements + 1
!
!                   if ( nlist+3 <= maxlist ) then
!                     list(nlist+1) = i2
!                     list(nlist+2) = j2
!                     list(nlist+3) = k2
!                   end if
!
!                   nlist = nlist + 3
!
!                 end if
!
!               end do
!             end do
!           end do
! !
! !  If any new neighbors were found, take the last one as the basis
! !  for a deeper search.
! !
!           if ( nabes > 0 ) then
!
!             if ( nstack+1 > maxstack ) then
!               go to 20
!             end if
!
!             stack(nstack+1) = nabes
!             nstack = nstack + 1
!             go to 10
!
!           end if
! !
! !  If the current search point had no new neighbors, drop it from the stack.
! !
!           ncan = stack(nstack) - 1
!           nstack = nstack - 3
!           stack(nstack) = ncan
! !
! !  If there are still any unused candidates at this level, take the
! !  last one as the basis for a deeper search.
! !
!           if ( stack(nstack) > 0 ) then
!             go to 10
!           end if
! !
! !  If there are no more unused candidates at this level, then we need
! !  to back up a level in the stack.  If there are any candidates at
! !  that earlier level, then we can still do more searching.
! !
!           nstack = nstack - 1
!
!           if ( nstack > 0 ) then
!             go to 10
!           end if
! !
! !  If we have exhausted the stack, we have completed this region.
! !  Tag the number of elements to the end of the region description list.
! !
!           nlist = nlist + 1
!           if ( nlist <= maxlist ) then
!             list(nlist) = nelements
!           end if
!
!         end if
!
!       end do
!     end do
!   end do
! !
! !  Print some warnings.
! !
!   if ( nlist > maxlist ) then
!     write ( *, * ) ' '
!     write ( *, * ) 'VOXEL_REGION - Warning!'
!     write ( *, * ) '  MAXLIST was too small to list the regions.'
!     write ( *, * ) '  Do not try to use the LIST array!'
!     write ( *, * ) '  The ISHOW data is OK, however.'
!   end if
!
!   return
! !
! !  Stack overflow.
! !
! 20    continue
!
!   write ( *, * ) ' '
!   write ( *, * ) 'VOXEL_REGION - Fatal error!'
!   write ( *, * ) '  The internal stack overflowed.'
!   write ( *, * ) '  The algorithm has failed.'
!
!   stop
! end subroutine
! subroutine voxel_step_3d ( i1, j1, k1, i2, j2, k2, inc, jnc, knc )
! !
! !*******************************************************************************
! !
! !! VOXEL_STEP_3D computes voxels along a line from a given point in 3D.
! !
! !
! !  Modified:
! !
! !    16 April 1999
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, integer I1, J1, K1, the coordinates of the base voxel from
! !    which the line begins.
! !
! !    Input/output, integer I2, J2, K2.
! !
! !    On input, these are the coordinates of the current voxel on
! !    the line.  For the first call, these might be I1, J1 and K1.
! !
! !    On output, these are the coordinates of the next voxel along
! !    the line.
! !
! !    Input, integer INC, JNC, KNC, the increments to the voxels.
! !    These values define the direction along which the line proceeds.
! !    However, the voxels on the line will typically be incremented
! !    by a fractional value of the vector (INC,JNC,KNC), and the
! !    result is essentially rounded.
! !
! !    If you input INC = JNC = KNC, then no movement is possible,
! !    and none is made.
! !
!   real alpha
!   real alphai
!   real alphaj
!   real alphak
!   integer i1
!   integer i2
!   integer inc
!   integer j1
!   integer j2
!   integer jnc
!   integer k1
!   integer k2
!   integer knc
! !
! !  Assuming for the moment that (I,J,K) can take on real values,
! !  points on the line have the form:
! !
! !    I = I1 + alpha * inc
! !    J = J1 + alpha * jnc
! !    K = K1 + alpha * knc
! !
!   if ( inc == 0 .and. jnc == 0 .and. knc == 0 ) then
!     return
!   end if
!
!   alpha = 0.0E+00
! !
! !  Compute the smallest ALPHA that will change I2, J2 or K2 by +-0.5.
! !
!   if ( inc > 0 ) then
!     alphai = ( real ( i2 - i1 ) + 0.5E+00 ) / real ( inc )
!   else if ( inc < 0 ) then
!     alphai = ( real ( i2 - i1 ) - 0.5E+00 ) / real ( inc )
!   else
!     alphai = huge ( alphai )
!   end if
!
!   if ( jnc > 0 ) then
!     alphaj = ( real ( j2 - j1 ) + 0.5E+00 ) / real ( jnc )
!   else if ( jnc < 0 ) then
!     alphaj = ( real ( j2 - j1 ) - 0.5E+00 ) / real ( jnc )
!   else
!     alphaj = huge ( alphaj )
!   end if
!
!   if ( knc > 0 ) then
!     alphak = ( real ( k2 - k1 ) + 0.5E+00 ) / real ( knc )
!   else if ( knc < 0 ) then
!     alphak = ( real ( k2 - k1 ) - 0.5E+00 ) / real ( knc )
!   else
!     alphaj = huge ( alphaj )
!   end if
! !
! !  The ALPHA of smallest positive magnitude represents the closest next voxel.
! !
!   alpha = huge ( alpha )
!
!   if ( alphai > 0.0E+00 ) then
!     alpha = min ( alpha, alphai )
!   end if
!
!   if ( alphaj > 0.0E+00 ) then
!     alpha = min ( alpha, alphaj )
!   end if
!
!   if ( alphak > 0.0E+00 ) then
!     alpha = min ( alpha, alphak )
!   end if
! !
! !  Move to the new voxel.  Whichever index just made the half
! !  step must be forced to take a whole step.
! !
!   if ( alpha == alphai ) then
!     i2 = i2 + sign ( 1, inc )
!     j2 = j1 + nint ( alpha * jnc )
!     k2 = k1 + nint ( alpha * knc )
!   else if ( alpha == alphaj ) then
!     i2 = i1 + nint ( alpha * inc )
!     j2 = j2 + sign ( 1, jnc )
!     k2 = k1 + nint ( alpha * knc )
!   else if ( alpha == alphak ) then
!     i2 = i1 + nint ( alpha * inc )
!     j2 = j1 + nint ( alpha * jnc )
!     k2 = k2 + sign ( 1, knc )
!   end if
!
!   return
! end subroutine
! subroutine xyz_to_radec ( x, y, z, ra, dec )
! !
! !*******************************************************************************
! !
! !! XYZ_TO_RADEC converts (X,Y,Z) to right ascension/declination coordinates.
! !
! !
! !  Discussion:
! !
! !    Given an XYZ point, compute its distance R from the origin, and
! !    regard it as lying on a sphere of radius R, whose axis is the Z
! !    axis.
! !
! !    The right ascension of the point is the "longitude", measured in hours,
! !    between 0 and 24, with the X axis having right ascension 0, and the
! !    Y axis having right ascension 6.
! !
! !    Declination measures the angle from the equator towards the north pole,
! !    and ranges from -90 (South Pole) to 90 (North Pole).
! !
! !  Modified:
! !
! !    02 December 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, real X, Y, Z, the coordinates of a point in 3D.
! !
! !    Output, real RA, DEC, the corresponding right ascension and declination.
! !
!   real dec
!   real norm_v
!   real phi
!   real ra
!   real theta
!   real x
!   real y
!   real z
! !
!   norm_v = sqrt ( x**2 + y**2 + z**2 )
!
!   phi = asin ( z / norm_v )
!
!   if ( cos ( phi ) == 0.0E+00 ) then
!     theta = 0.0E+00
!   else
!     theta = atan4 ( y, x )
!   end if
!
!   dec = radians_to_degrees ( phi )
!   ra = radians_to_degrees ( theta ) / 15.0E+00
!
!   return
! end subroutine

end module

